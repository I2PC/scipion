
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

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

path1 = '/home/javiermota/ScipionUserData/projects/CNBScipionCourse/Runs/023211_XmippProtGenerateReprojections/extra/anglesCont.xmd'
root = '/home/javiermota/ScipionUserData/projects/CNBScipionCourse/'
path2 = '/home/javiermota/rinchen3/ScipionUserData/projects/coss/Runs/037132_XmippProtGenerateReprojections/extra/anglesCont.xmd'
root2 = '/home/javiermota/rinchen3/ScipionUserData/projects/coss/'

cnbCourseparticles = '/home/javiermota/ScipionUserData/projects' \
                     '/CNBScipionCourse/Runs/006500_XmippProtCropResizeParticles/extra/output_images.xmd'
path3 = '/home/javiermota/rinchen3/ScipionUserData/projects/coss/Runs' \
        '/037132_XmippProtGenerateReprojections/extra/anglesCont.xmd'
path3noise = '/home/javiermota/rinchen3/ScipionUserData/projects/coss/Runs/037219_XmippProtAddNoiseParticles/extra/Noisy.xmd'

pathGallery = '/home/javiermota/rinchen3/ScipionUserData/projects/coss/Runs' \
              '/037163_XmippProtCreateGallery/images.xmd'

galleryCourse = '/home/javiermota/ScipionUserData/projects/CNBScipionCourse' \
                '/Runs/023535_XmippProtCreateGallery/images.xmd'

pathtest = '/home/javiermota/rinchen3/ScipionUserData/projects/coss/Runs/020417_XmippProtExtractParticles/images.xmd'
cossparticles = '/home/javiermota/rinchen3/ScipionUserData/projects/coss/Runs\
/009654_XmippProtExtractParticles/images.xmd'

class GAN():
    def __init__(self):
        self.img_rows = 60 #40
        self.img_cols = 60 #40
        self.channels = 1
        self.shape = self.img_rows*self.img_cols
        self.img_shape = (self.img_rows, self.img_cols, self.channels)

        optimizer = Adam(0.00007, 0.5)

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

    def ExtractInfoMetadata(self, path, root, label, crop, size, norm=-1):
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
            Imresize = cv2.resize(Data, (size, size),
                                  interpolation=cv2.INTER_CUBIC)
            Imresize = Imresize[crop:-crop, crop:-crop]

            if norm == -1:
                Imnormalize = 2*(Imresize-np.min(Imresize))/(np.max(Imresize)-np.min(
                                  Imresize))-1
            elif norm == 0:
                Imnormalize = Imresize
            else:
                Imnormalize = (Imresize-np.mean(Imresize))/np.std(Imresize)

            Image.append(Imnormalize)
            if cont > 2000:
                break
            cont += 1

        Image = np.array(Image).astype('float')
        Image = Image.reshape(len(Image), Image.shape[1], Image.shape[2], 1)

        return Image

    def addNoise(self, image):

        imageNoise = []
        projections = []
        for std in np.arange(0.5,7,0.5):

            noise = np.random.normal(0.0, std, image.shape)
            imageNoise.append(image + noise)
            projections.append(image)

        projections = np.asarray(projections).astype('float32')
        projections = projections.reshape(len(projections) * projections.shape[
            1], projections.shape[2], projections.shape[3], 1)
        imageNoise = np.asarray(imageNoise).astype('float32')
        imageNoise = imageNoise.reshape(len(imageNoise) * imageNoise.shape[
            1],imageNoise.shape[2],imageNoise.shape[3], 1)

        return projections, imageNoise


    def normalization(self, image, type='mean'):

        if type == 'mean':
            Imnormalize = (image - np.mean(image)) / np.std(image)

        if type == -1:
            Imnormalize = 2 * (image - np.min(image)) / (
                        np.max(image) - np.min(
                    image)) - 1

        return Imnormalize

    def applyTransform(self, image):

        noise = []
        proj = []
        for i,im in enumerate(image):

            angle = np.random.randint(-180,180,5).astype('float64')
            shiftx = np.random.randint(-5,5,5)
            shifty = np.random.randint(-5, 5, 5)

            for j, r in enumerate(angle):

                imRotate = rotate(im,r,mode='wrap')
                shift = AffineTransform(translation=[shiftx[j], shifty[j]])
                imShift = warp(imRotate, shift, mode='wrap',
                               preserve_range=True)
                '''imRotateNoise = rotate(imageN[i], r, mode='wrap')
                imShiftNoise = warp(imRotateNoise, shift, mode='wrap',
                               preserve_range=True)'''

                proj.append(self.normalization(imShift,-1))
                proj.append(self.normalization(imRotate, -1))
                '''noise.append(self.normalization(imShiftNoise, 'mean'))
                noise.append(self.normalization(imRotateNoise, 'mean'))'''

        proj = np.asarray(proj).astype('float64')

        return proj


    def build_generator(self):

        #noise_shape = (self.shape,)

        model = Sequential()

        input_img = Input(shape=self.img_shape,
                          name='input')
        x = Conv2D(32,(9,9),padding='same')(input_img)
        x = BatchNormalization(momentum=0.8)(x)
        x = LeakyReLU(alpha=0.2)(x)
        x = AveragePooling2D((2, 2), padding='same')(x)
        x = Conv2D(64, (5, 5), padding='same')(x)
        x = BatchNormalization(momentum=0.8)(x)
        x = LeakyReLU(alpha=0.2)(x)
        encoded = AveragePooling2D((2, 2), padding='same', name='encoder')(x)
        x = Conv2D(128, (3, 3), padding='same')(encoded)
        x = BatchNormalization(momentum=0.8)(x)
        x = LeakyReLU(alpha=0.2)(x)
        x = UpSampling2D(2)(x)
        x = Conv2D(64, (7, 7), padding='same')(x)
        x = BatchNormalization(momentum=0.8)(x)
        x = LeakyReLU(alpha=0.2)(x)
        x = UpSampling2D(2)(x)
        '''x = Conv2D(32, (3, 3), padding='same')(x)
        x = BatchNormalization(momentum=0.8)(x)
        x = LeakyReLU(alpha=0.2)(x)'''
        x = Conv2D(1, (7, 7), padding='same')(x)
        decoded = Activation('tanh')(x)
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

        model.add(Conv2D(32, (5,5), padding='same', input_shape=self.img_shape))
        model.add(LeakyReLU())
        model.add(Conv2D(64, (3, 3), padding='same'))
        model.add(BatchNormalization())
        model.add(LeakyReLU())
        model.add(Conv2D(128, (3, 3), padding='same'))
        model.add(BatchNormalization())
        model.add(LeakyReLU())
        model.add(Conv2D(256, (3, 3), padding='same'))
        model.add(BatchNormalization())
        model.add(LeakyReLU())
        model.add(Conv2D(1, (3, 3), padding='same'))
        model.add(Activation('sigmoid'))

        '''model.add(Flatten(input_shape=img_shape))
        model.add(Dense(512))
        model.add(LeakyReLU(alpha=0.2))
        #model.add(Dropout(0.5))
        model.add(Dense(256))
        model.add(LeakyReLU(alpha=0.2))
        #model.add(Dropout(0.5))
        model.add(Dense(1, activation='sigmoid'))'''
        model.summary()

        img = Input(shape=img_shape)
        validity = model(img)

        return Model(img, validity)

    def train(self, X_train, noise, epochs, batch_size=128, save_interval=50):

        # Load the dataset
        #(X_train, _), (_, _) = mnist.load_data()

        # Rescale -1 to 1
        #X_train = (X_train.astype(np.float32) - 127.5) / 127.5
        #X_train = np.expand_dims(X_train, axis=3)

        half_batch = int(batch_size / 2)

        self.lossD = []
        self.lossG = []
        for epoch in range(epochs):

            # ---------------------
            #  Train Discriminator
            # ---------------------

            # Select a random half batch of images
            idx = np.random.randint(0, X_train.shape[0], half_batch)
            imgs = X_train[idx]

            #noise7 = np.random.normal(0, 1, (half_batch, self.shape))
            noise1 = noise[idx]
            #noise1 = noise1.reshape(len(noise1), 1600)

            # Generate a half batch of new images
            gen_imgs = self.generator.predict(noise1)
            '''plt.imshow(np.squeeze(gen_imgs[0]))
            plt.gray()
            plt.show()'''
            # Train the discriminator
            d_loss_real = self.discriminator.train_on_batch(imgs, np.ones(
                (half_batch, 1)))
            d_loss_fake = self.discriminator.train_on_batch(gen_imgs, np.zeros(
                (half_batch, 1)))
            d_loss = 0.5 * np.add(d_loss_real, d_loss_fake)

            # ---------------------
            #  Train Generator
            # ---------------------
            idy = np.random.randint(0, X_train.shape[0], batch_size)
            noise2 = noise[idy]
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
            # If at save interval => save generated image samples
            if epoch % save_interval == 0:
                self.save_imgs(X_train, noise, epoch)

        self.generator.save('AdversarialDenoisingCoss.h5')

    def predict(self, model, data):

        if isinstance(data, basestring):
            test = self.ExtractInfoMetadata(data, root2, xmipp.MDL_IMAGE, 5,
                                        70,1)
        else:
            test = data
        prediction = model.predict(test)

        '''for i, im in enumerate(prediction):
            plt.subplot(1, 2, 1)
            plt.imshow(np.squeeze(im))
            plt.gray()
            plt.subplot(1, 2, 2)
            plt.imshow(np.squeeze(test[i]))
            plt.gray()
            plt.show()'''

        return prediction

    def save_imgs(self, X_train, noise, epoch):
        filename = "denoise_%d.png"
        idx = np.random.randint(0, X_train.shape[0], 10)
        noise = noise[idx]#np.random.normal(0, 1, (r * c, self.shape))
        #noise2 = noise.reshape(len(noise), 1600)
        true = X_train[idx]
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

    '''X_train1 = gan.ExtractInfoMetadata(path1, root, xmipp.MDL_IMAGE_REF,
                                        5, 70)'''
    X_train2 = gan.ExtractInfoMetadata(path3, root2, xmipp.MDL_IMAGE,
                                        5, 70)

    # self.X_train = X_train2#np.concatenate((X_train1, X_train2), axis=0)
    '''noise1 = gan.ExtractInfoMetadata(path1, root, xmipp.MDL_IMAGE,
                                      5, 70, 1)'''
    '''noise2 = gan.ExtractInfoMetadata(path3noise, root2, xmipp.MDL_IMAGE,
                                      5, 70, 1)'''
    # self.noise = noise2#np.concatenate((noise1, noise2), axis=0)

    X = gan.applyTransform(X_train2)
    X_train, noise = gan.addNoise(X)
    #X_train = np.concatenate((X_train2, X), axis=0)
    #noise = np.concatenate((noise2, y), axis=0)
    #gan.predict()
    #image = gan.ExtractInfoMetadata(path2, root2, xmipp.MDL_IMAGE_REF,5, 70)
    gan.train(X_train, noise, epochs=10000, batch_size=32, save_interval=200)
    noise1 = gan.ExtractInfoMetadata(cnbCourseparticles, root, xmipp.MDL_IMAGE,
                                     5, 70, 1)
    model = load_model('AdversarialDenoisingCoss.h5')
    predict = gan.predict(model, noise1)

    for i,im in enumerate(predict):
        plt.subplot(1,2,1)
        plt.imshow(np.squeeze(im), cmap='gray')
        plt.subplot(1, 2, 2)
        plt.imshow(np.squeeze(noise1[i]), cmap='gray')
        '''plt.subplot(1, 3, 3)
        plt.imshow(np.squeeze(X_train1[i]), cmap='gray')'''
        plt.show()

    '''
    NoiseEnhanced = []
    for i, im in enumerate(predict):

        NewNoise = im + noise[i]
        #NewNoise = gan.normalization(NewNoise,'mean')
        NoiseEnhanced.append(NewNoise)

    NoiseEnhanced = np.array(NoiseEnhanced).astype('float')
    NoiseEnhanced = NoiseEnhanced.reshape(len(NoiseEnhanced),
                                          NoiseEnhanced.shape[1],
                                          NoiseEnhanced.shape[2],1)
    gan.train(X_train, noise, epochs=20000, batch_size=32, save_interval=200)'''

