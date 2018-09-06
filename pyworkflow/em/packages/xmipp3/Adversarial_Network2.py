
import numpy as np
import time
from tensorflow.examples.tutorials.mnist import input_data

from keras.models import Sequential, Model, load_model
from keras.layers import Dense, Activation, Flatten, Reshape
from keras.layers import Conv2D, Conv2DTranspose, UpSampling2D, Input, \
    MaxPooling2D, Subtract, AveragePooling2D
from keras.layers import LeakyReLU, Dropout, GaussianNoise
from keras.layers import BatchNormalization
from keras.optimizers import Adam, RMSprop, SGD
from keras.datasets import mnist
import xmipp
import cv2
import tensorflow as tf
import keras
from keras.callbacks import History
from skimage.transform import rotate, AffineTransform, warp
from skimage.feature import match_template

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

path1 = '/home/javiermota/ScipionUserData/projects/CNBScipionCourse/Runs/023211_XmippProtGenerateReprojections/extra/anglesCont.xmd'
root = '/home/javiermota/ScipionUserData/projects/CNBScipionCourse/'
path2 = '/home/javiermota/rinchen3/ScipionUserData/projects/coss/Runs/037537_XmippProtGenerateReprojections/extra/anglesCont.xmd'
root2 = '/home/javiermota/rinchen3/ScipionUserData/projects/coss/'

cnbCourseparticles = '/home/javiermota/ScipionUserData/projects/CNBScipionCourse/Runs/024405_XmippProtCropResizeParticles/extra/output_images.xmd'

path1noise = '/home/javiermota/ScipionUserData/projects/CNBScipionCourse/Runs/023591_XmippProtAddNoiseParticles/extra/Noisy.xmd'

path2noise = '/home/javiermota/rinchen3/ScipionUserData/projects/coss/Runs' \
             '/037624_XmippProtAddNoiseParticles/extra/Noisy.xmd'

pathGallery = '/home/javiermota/rinchen3/ScipionUserData/projects/coss/Runs' \
              '/037568_XmippProtCreateGallery/images.xmd'

galleryCourse = '/home/javiermota/ScipionUserData/projects/CNBScipionCourse' \
                '/Runs/023535_XmippProtCreateGallery/images.xmd'

pathtest = '/home/javiermota/rinchen3/ScipionUserData/projects/coss/Runs/020417_XmippProtExtractParticles/images.xmd'
cossparticles = '/home/javiermota/rinchen3/ScipionUserData/projects/coss/Runs/009654_XmippProtExtractParticles/images.xmd'
highresParticlesCoss = '/home/javiermota/rinchen3/ScipionUserData/projects' \
                       '/coss/Runs/037851_XmippProtCropResizeParticles/extra' \
                       '/output_images.xmd'

reprojectionParticles = '/home/javiermota/rinchen3/ScipionUserData/projects/coss/Runs/038007_XmippProtCropResizeParticles/extra/output_images.xmd'

class GAN():
    def __init__(self):
        self.img_rows = 100 #40
        self.img_cols = 100 #40
        self.channels = 1
        self.shape = self.img_rows*self.img_cols
        self.img_shape = (self.img_rows, self.img_cols, self.channels)

        optimizer = Adam(0.0001, 0.5) #0.00007, 0.5

        # Build and compile the discriminator
        self.discriminator = self.build_discriminator()
        self.discriminator.compile(loss='binary_crossentropy',
                                   optimizer=optimizer,
                                   metrics=['accuracy'])

        # Build and compile the generator
        self.generator = self.build_generator()
        self.generator.compile(loss='mean_squared_error', optimizer=optimizer)

        # The generator takes noise as input and generated imgs
        z = Input(shape=self.img_shape)#Input(shape=(self.shape,))
        img = self.generator(z)

        # For the combined model we will only train the generator
        self.discriminator.trainable = False

        # The valid takes generated images as input and determines validity
        valid = self.discriminator(img)

        # The combined model  (stacked generator and discriminator) takes
        # noise as input => generates images => determines validity
        self.combined = Model(z, valid)
        self.combined.compile(loss='binary_crossentropy', optimizer=optimizer)

    def ExtractInfoMetadata(self, path, root, label, size=0, norm=-1):
        metadata = xmipp.MetaData(path)
        Image = []
        I = xmipp.Image()
        cont = 0
        for itemId in metadata:
            fn = metadata.getValue(label, itemId)
            n = fn.split('@')
            fn = n[0] + '@' + root + n[1]
            I.read(fn)
            #I.write('/home/javiermota/image')
            Data = I.getData()
            if size == 0:
                Imresize = Data
            else:
                Imresize = cv2.resize(Data, (size, size),
                                    interpolation=cv2.INTER_CUBIC)

            if norm == -1:
                Imnormalize = 2*(Imresize-np.min(Imresize))/(np.max(Imresize)-np.min(
                                  Imresize))-1
            elif norm == 0:
                Imnormalize = Imresize
            elif norm == 1:
                Imnormalize = (Imresize-np.min(Imresize))/(np.max(Imresize)-np.min(
                                  Imresize))
            else:
                Imnormalize = (Imresize-np.mean(Imresize))/np.std(Imresize)

            Image.append(Imnormalize)
            if cont > 2000:
                break
            cont += 1

        Image = np.array(Image).astype('float')
        Image = Image.reshape(len(Image), Image.shape[1], Image.shape[2], 1)

        return Image

    def normalization(self, image, type='mean'):

        NormalizedImage = []
        for im in image:
            if type == 'mean':
                Imnormalize = (im - np.mean(im)) / np.std(im)

            if type == -1:
                Imnormalize = 2 * (im - np.min(im)) / (
                            np.max(im) - np.min(im)) - 1

            if type == 1:
                Imnormalize = (im - np.min(im)) / (
                        np.max(im) - np.min(im))

            if type == 'RGB':
                Imnormalize = np.floor(im*255)

            NormalizedImage.append(Imnormalize)

        NormalizedImage = np.array(NormalizedImage).astype('float')
        NormalizedImage = NormalizedImage.reshape(len(NormalizedImage), NormalizedImage.shape[1], NormalizedImage.shape[2], 1)

        return NormalizedImage

    def addNoise(self, image):

        levelsNoise = np.arange(0.5,3.5,0.1)
        k = np.random.randint(0,len(levelsNoise))
        noise = np.random.normal(0.0, levelsNoise[k], image.shape)
        imageNoise = image + noise

        return imageNoise

    def applyTransform(self, image):

        angle = np.random.randint(-180,180)
        shiftx = np.random.randint(-10,10)
        shifty = np.random.randint(-10,10)

        imRotate = rotate(image,angle,mode='wrap')
        shift = AffineTransform(translation=[shiftx, shifty])
        imShift = warp(imRotate, shift, mode='wrap',
                       preserve_range=True)

        return imShift

    def generate_data(self, images, batch_size):

        proj = []
        noiseImage = []
        for j in range(0,batch_size):
            idx = np.random.randint(0, X_train.shape[0])
            img = images[idx]

            projection = self.applyTransform(img)
            noise = self.addNoise(projection)
            proj.append(projection)
            noiseImage.append(noise)

        projections = np.asarray(proj).astype('float32')
        imageNoise = np.asarray(noiseImage).astype('float32')

        return projections, imageNoise

    def build_generator(self):

        #noise_shape = (self.shape,)

        #model = Sequential()

        input_img = Input(shape=self.img_shape,
                          name='input')
        x = Conv2D(32,(21, 21),padding='same')(input_img)
        x = BatchNormalization(momentum=0.8)(x)
        x = LeakyReLU(alpha=0.2)(x)
        x = AveragePooling2D((2, 2), padding='same')(x)
        x = Conv2D(64, (15, 15), padding='same')(x)
        x = BatchNormalization(momentum=0.8)(x)
        x = LeakyReLU(alpha=0.2)(x)
        encoded = AveragePooling2D((2, 2), padding='same', name='encoder')(x)
        x = Conv2D(64, (3, 3), padding='same')(encoded)
        x = BatchNormalization(momentum=0.8)(x)
        x = LeakyReLU(alpha=0.2)(x)
        x = UpSampling2D(2)(x)
        x = Conv2D(32, (7, 7), padding='same')(x)
        x = BatchNormalization(momentum=0.8)(x)
        x = LeakyReLU(alpha=0.2)(x)
        x = UpSampling2D(2)(x)
        '''x = Conv2D(32, (3, 3), padding='same')(x)
        x = BatchNormalization(momentum=0.8)(x)
        x = LeakyReLU(alpha=0.2)(x)'''
        x = Conv2D(1, (9, 9), padding='same')(x)
        decoded = Activation('linear')(x)
        '''
        model.add(Dense(256, input_shape=noise_shape))
        model.add(LeakyReLU(alpha=0.2))
        model.add(BatchNormalization(momentum=0.8))
        model.add(Dense(512))
        model.add(LeakyReLU(alpha=0.2))
        model.add(BatchNormalization(momentum=0.8))
        model.add(Dense(1024))
        model.add(LeakyReLU(alpha=0.2))
        model.add(BatchNormalization(momentum=0.8))
        model.add(Dense(np.prod(self.img_shape), activation='tanh'))
        model.add(Reshape(self.img_shape))'''

        model = Model(input_img, decoded)
        model.summary()

        noise = Input(shape=self.img_shape)#Input(shape=noise_shape)
        img = model(noise)

        return Model(noise, img)


    def build_discriminator(self):

        img_shape = (self.img_rows, self.img_cols, self.channels)

        model = Sequential()

        '''model.add(Conv2D(64, (3, 3), padding='same', input_shape=img_shape))
        model.add(LeakyReLU(alpha=0.2))
        model.add(Conv2D(128, (3, 3), padding='same'))
        model.add(LeakyReLU(alpha=0.2))
        model.add(Conv2D(32, (3, 3), padding='same'))
        #model.add(BatchNormalization(momentum=0.8))
        model.add(LeakyReLU(alpha=0.2))'''
        model.add(Flatten(input_shape=img_shape))
        model.add(Dense(512))
        model.add(LeakyReLU(alpha=0.2))
        model.add(Dense(256))
        model.add(LeakyReLU(alpha=0.2))
        model.add(Dense(128))
        model.add(LeakyReLU(alpha=0.2))
        model.add(Dense(64))
        model.add(LeakyReLU(alpha=0.2))
        model.add(Dense(1, activation='sigmoid'))
        model.summary()

        img = Input(shape=img_shape)
        validity = model(img)

        return Model(img, validity)

    def train(self, X_train, epochs, batch_size=128, save_interval=50):

        half_batch = int(batch_size / 2)

        self.lossD = []
        self.lossG = []
        lossEpoch = []
        for epoch in range(epochs):

            # ---------------------
            #  Train Discriminator
            # ---------------------
            imgs, noise1 = self.generate_data(X_train, half_batch)
            # Select a random half batch of images

            # Generate a half batch of new images
            gen_imgs = self.generator.predict(noise1)
            # Train the discriminator
            d_loss_real = self.discriminator.train_on_batch(imgs, np.round(
                np.random.uniform(0.8,1.0,half_batch),1))#np.ones((half_batch, 1)))
            d_loss_fake = self.discriminator.train_on_batch(gen_imgs,
                                                            np.round(
                                                                np.random.uniform(0.0,0.2,half_batch),1)) #np.zeros((half_batch, 1)))
            d_loss = 0.5 * np.add(d_loss_real, d_loss_fake)

            # ---------------------
            #  Train Generator
            # ---------------------
            imgs2, noise2 = self.generate_data(X_train, batch_size)
            #noise2 = noise2.reshape(len(noise2), 1600)
            #noise8 = np.random.normal(0, 1, (batch_size, self.shape))

            # The generator wants the discriminator to label the generated samples
            # as valid (ones)
            valid_y = np.array([1] * batch_size)

            # Train the generator
            g_loss = self.combined.train_on_batch(noise2, valid_y)

            # Plot the progress
            print ("%d [D loss: %f, acc.: %.2f%%] [G loss: %f]" % (
            epoch, d_loss[0], 100 * d_loss[1], g_loss))

            self.lossD.append(d_loss[0])
            self.lossG.append(g_loss)
            lossEpoch.append(d_loss[0])
            # If at save interval => save generated image samples
            if epoch % save_interval == 0:
                print "MeanLoss = ", np.mean(lossEpoch)
                self.save_imgs(X_train, epoch)
                lossEpoch = []

        self.generator.save('AdversarialDenoisingCoss100CNB.h5')

    def predict(self, model, data):

        if isinstance(data, basestring):
            test = self.ExtractInfoMetadata(data, root2, xmipp.MDL_IMAGE_ORIGINAL,
                                            80,1)
        else:
            test = data
        prediction = model.predict(test)

        return prediction

    def save_imgs(self, X_train,epoch):
        filename = "denoise_%d.png"
        true, noise = self.generate_data(X_train, 10)

        gen_imgs = self.generator.predict(noise)

        # Rescale images 0 - 1
        gen_imgs = 0.5 * gen_imgs + 0.5

        fig, axs = plt.subplots(10, 3)
        cnt = 0
        for i in range(10):
            axs[i, 0].imshow(true[cnt, :, :, 0], cmap='gray')
            axs[i, 0].axis('off')
            axs[i, 1].imshow(noise[cnt, :, :, 0], cmap='gray')
            axs[i, 1].axis('off')
            axs[i, 2].imshow(gen_imgs[cnt, :, :, 0], cmap='gray')
            axs[i, 2].axis('off')
            cnt += 1
        plt.savefig(filename % epoch)
        plt.close()

        filenameD = 'loss_D.png'
        plt.plot(np.arange(0,len(self.lossD)), np.array(self.lossD))
        plt.savefig(filenameD)
        plt.close()

        filenameG = 'loss_G.png'
        plt.plot(np.arange(0, len(self.lossG)), np.array(self.lossG))
        plt.savefig(filenameG)
        plt.close()

if __name__ == '__main__':
    gan = GAN()

    '''X_train = gan.ExtractInfoMetadata(path1noise, root, xmipp.MDL_IMAGE,
                                        5, 70)'''
    X_train = gan.ExtractInfoMetadata(path1, root, xmipp.MDL_IMAGE_REF, 100,
                                       -1)
    #gan.train(X_train, epochs=10000, batch_size=32, save_interval=200)
    noise1 = gan.ExtractInfoMetadata(reprojectionParticles, root2,
                                     xmipp.MDL_IMAGE, 0, 2)
    projection = gan.ExtractInfoMetadata(reprojectionParticles, root2,
                                         xmipp.MDL_IMAGE_REF, 100, 1)
    model = load_model('AdversarialDenoisingCoss100.h5')
    predict = gan.predict(model, noise1)
    '''model2 = load_model('AdversarialDenoisingCoss100CNB.h5')
    predict2 = gan.predict(model2, noise1)'''

    NoiseEnhanced = []
    for i, im in enumerate(predict):

        NewNoise = (im + noise1[i])/2.0
        NoiseEnhanced.append(NewNoise)

    NoiseEnhanced = np.asarray(NoiseEnhanced).astype('float32')
    NoiseEnhanced = NoiseEnhanced.reshape(len(NoiseEnhanced),
                                          NoiseEnhanced.shape[1],
                                          NoiseEnhanced.shape[2],1)
    predictEnhanced = gan.predict(model, NoiseEnhanced)
    predictEnhanced = gan.normalization(predictEnhanced, 1)

    '''NoiseEnhanced2 = []
    for i, im in enumerate(predict2):
        NewNoise = (im + noise1[i]) / 2.0
        NoiseEnhanced2.append(NewNoise)

    NoiseEnhanced2 = np.asarray(NoiseEnhanced2).astype('float32')
    NoiseEnhanced2 = NoiseEnhanced2.reshape(len(NoiseEnhanced2),
                                          NoiseEnhanced2.shape[1],
                                          NoiseEnhanced2.shape[2], 1)
    predictEnhanced2 = gan.predict(model, NoiseEnhanced2)
    predictEnhanced2 = gan.normalization(predictEnhanced2, 1)'''

    noise1 = gan.normalization(noise1, 1)
    for i,im in enumerate(predictEnhanced):
        diff = im - projection[i]
        m_norm = np.sum(abs(diff))
        print m_norm
        correlation = match_template(im, projection[i])
        #correlation2 = match_template(predictEnhanced2[i], projection[i])
        correlation2 = match_template(im, noise1[i])
        #correlation3 = match_template(projection[i], noise1[i])
        print correlation[0], correlation2[0]#, correlation3[0]
        plt.subplot(2,2,1)
        plt.imshow(np.squeeze(im), cmap='gray')
        plt.subplot(2, 2, 2)
        plt.imshow(np.squeeze(noise1[i]), cmap='gray')
        '''plt.subplot(2, 2, 3)
        plt.imshow(np.squeeze(predictEnhanced2[i]), cmap='gray')'''
        plt.subplot(2, 2, 3)
        plt.imshow(np.squeeze(projection[i]), cmap='gray')
        plt.show()


