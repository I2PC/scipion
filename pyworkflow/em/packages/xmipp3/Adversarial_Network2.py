
import numpy as np
import time
from tensorflow.examples.tutorials.mnist import input_data
from tensorflow.python.client import device_lib
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
from keras import backend as K

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import glob
from skimage.util import view_as_blocks
import skimage
import hashlib
from scipy.stats import pearsonr

path1 = '/home/javiermota/ScipionUserData/projects/CNBScipionCourse/Runs/023211_XmippProtGenerateReprojections/extra/anglesCont.xmd'
root = '/home/javiermota/ScipionUserData/projects/CNBScipionCourse/'
path2 = '/home/javiermota/rinchen3/ScipionUserData/projects/coss/Runs/037537_XmippProtGenerateReprojections/extra/anglesCont.xmd'
root2 = '/home/javiermota/rinchen3/ScipionUserData/projects/coss/'
root4 = '/home/javiermota/rinchen3/ScipionUserData/projects/10061_deep/'

cnbCourseparticles = '/home/javiermota/ScipionUserData/projects/CNBScipionCourse/Runs/025220_XmippProtCropResizeParticles/extra/output_images.xmd'

path1noise = '/home/javiermota/ScipionUserData/projects/CNBScipionCourse/Runs/023591_XmippProtAddNoiseParticles/extra/Noisy.xmd'

path2noise = '/home/javiermota/rinchen3/ScipionUserData/projects/coss/Runs' \
             '/037624_XmippProtAddNoiseParticles/extra/Noisy.xmd'

galleryBetagal = '/home/javiermota/rinchen3/ScipionUserData/projects' \
                '/10061_deep/Runs/020521_XmippProtCreateGallery/images.xmd'
galleryBetagalHighres = '/home/javiermota/rinchen3/ScipionUserData/projects/10061_deep/Runs/021245_XmippProtCreateGallery/images.xmd'
betagal = '/home/javiermota/rinchen3/ScipionUserData/projects/10061_deep/Runs' \
          '/020481_XmippProtCropResizeParticles/extra/output_images.xmd'
galleryCourse = '/home/javiermota/ScipionUserData/projects/CNBScipionCourse' \
                '/Runs/026647_XmippProtCreateGallery/images.xmd'
betagalRepro = '/home/javiermota/rinchen3/ScipionUserData/projects/10061_deep' \
               '/Runs/020811_XmippProtGenerateReprojections/extra/anglesCont' \
               '.xmd'
pathtest = '/home/javiermota/rinchen3/ScipionUserData/projects/coss/Runs/020417_XmippProtExtractParticles/images.xmd'
cossparticles = '/home/javiermota/rinchen3/ScipionUserData/projects/coss/Runs/009654_XmippProtExtractParticles/images.xmd'
highresParticlesCoss = '/home/javiermota/rinchen3/ScipionUserData/projects' \
                       '/coss/Runs/037851_XmippProtCropResizeParticles/extra' \
                       '/output_images.xmd'

reprojectionParticles = '/home/javiermota/rinchen3/ScipionUserData/projects/coss/Runs/038386_XmippProtCropResizeParticles/extra/output_images.xmd'

groel = '/home/javiermota/ScipionUserData/projects/Frank10000_70S/Runs' \
        '/006028_XmippProtGenerateReprojections/extra/anglesCont.xmd'
root3 = '/home/javiermota/ScipionUserData/projects/Frank10000_70S/'

galleryCNB = '/home/javiermota/ScipionUserData/projects/CNBScipionCourse/Runs' \
             '/026647_XmippProtCreateGallery/images.xmd'
micrographsCoss = '/home/javiermota/rinchen3/ScipionUserData/projects/coss/Runs/040262_XmippProtPreprocessMicrographs/extra/'
micrographsCNB = '/home/javiermota/ScipionUserData/projects/CNBScipionCourse/Runs/027537_XmippProtPreprocessMicrographs/extra/'

class GAN():
    def __init__(self):
        self.img_rows = 100 #40
        self.img_cols = 100 #40
        self.channels = 1
        self.shape = self.img_rows*self.img_cols
        self.img_shape = (self.img_rows, self.img_cols, self.channels)

        optimizer = Adam(0.0001, 0.5)

        # Build and compile the discriminator
        self.discriminator = self.build_discriminator()
        self.discriminator.compile(loss='binary_crossentropy',
                                   optimizer=optimizer,
                                   metrics=['accuracy'])

        # Build and compile the generator
        self.generator = self.build_generator()
        self.generator.compile(loss=self.dice_coef, optimizer=optimizer)

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

    def ExtractInfoMetadata(self, path, root, label, size=0, norm=-1,
                            reshape=True):
        metadata = xmipp.MetaData(path)
        Image = []
        I = xmipp.Image()
        cont = 0
        for itemId in metadata:
            fn = metadata.getValue(label, itemId)
            n = fn.split('@')
            fn = n[0] + '@' + root + n[1]
            I.read(fn)
            Data = I.getData()
            #I.write('/home/javiermota/image.mrc')
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
        if reshape:
            Image = Image.reshape(len(Image), Image.shape[1], Image.shape[2], 1)

        return Image

    def importMicrographs(self, micrograph):
        I = xmipp.Image()
        micrographs = []
        for im in glob.glob(micrograph+'*.mrc'):
            I.read(im)
            Data = I.getData()
            micrographs.append(Data)

        micrographs = np.asarray(micrographs).astype('float32')
        micrographs = micrographs.reshape(len(micrographs),
                                          micrographs.shape[1],
                                          micrographs.shape[2],1)
        return  micrographs



    def normalization(self, image, type='mean', reshape=True):

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
        if reshape:
            if  len(np.shape(NormalizedImage)) > 2:
                NormalizedImage = NormalizedImage.reshape(len(NormalizedImage), NormalizedImage.shape[1], NormalizedImage.shape[2], 1)
            else:
                NormalizedImage = NormalizedImage.reshape(1,
                                                          NormalizedImage.shape[
                                                              0], NormalizedImage.shape[1], 1)

        return NormalizedImage

    def addNoise(self, image):

        levelsNoise = np.arange(0.5,2.5,0.1)
        k = np.random.randint(0,len(levelsNoise))
        noise = np.random.normal(0.0, levelsNoise[k], image.shape)
        imageNoise = image + noise

        return imageNoise

    def applyTransform(self, image):

        angle = np.random.randint(-180,180)
        shiftx = np.random.randint(-15,15)
        shifty = np.random.randint(-15,15)

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
        x = Conv2D(64,(3, 3),padding='same')(input_img)
        x = LeakyReLU(alpha=0.2)(x)
        encoded = BatchNormalization(momentum=0.8)(x)
        '''x = Conv2D(32, (5, 5), padding='same')(x)
        x = BatchNormalization(momentum=0.8)(x)
        x = LeakyReLU(alpha=0.2)(x)'''
        '''x = AveragePooling2D((2, 2), padding='same')(x)
        x = Conv2D(64, (15, 15), padding='same')(x)
        x = BatchNormalization(momentum=0.8)(x)
        x = LeakyReLU(alpha=0.2)(x)'''
        #encoded = AveragePooling2D((2, 2), padding='same', name='encoder')(x)
        x = Conv2D(128, (3, 3), padding='same')(encoded)
        x = LeakyReLU(alpha=0.2)(x)
        x = BatchNormalization(momentum=0.8)(x)
        #x = UpSampling2D(2)(x)
        '''x = Conv2D(32, (7, 7), padding='same')(x)
        x = BatchNormalization(momentum=0.8)(x)
        x = LeakyReLU(alpha=0.2)(x)
        x = UpSampling2D(2)(x)'''
        x = Conv2D(1, (5, 5), padding='same')(x)
        decoded = Activation('linear')(x)


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
        model.add(Dense(1024))
        model.add(LeakyReLU(alpha=0.2))
        model.add(Dense(512))
        model.add(LeakyReLU(alpha=0.2))
        '''model.add(Dense(512))
        model.add(LeakyReLU(alpha=0.2))'''
        '''model.add(Dense(256))
        model.add(LeakyReLU(alpha=0.2))'''
        '''model.add(Dense(64))
        model.add(LeakyReLU(alpha=0.2))'''
        model.add(Dense(1, activation='sigmoid'))
        model.summary()

        img = Input(shape=img_shape)
        validity = model(img)

        return Model(img, validity)

    def train(self, X_train, epochs,batch_size=128,
              save_interval=50):

        half_batch = int(batch_size / 2)

        self.lossD = []
        self.lossG = []
        lossEpoch = []
        self.true, self.noise = self.generate_data(X_train, 100)
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
                np.random.uniform(0.9,1.0,half_batch),1))#np.ones((
            # half_batch, 1)))
            d_loss_fake = self.discriminator.train_on_batch(gen_imgs,
                                                            np.round(
                                                                np.random.uniform(0.0,0.1,half_batch),1)) #np.zeros((half_batch, 1)))
            d_loss = 0.5 * np.add(d_loss_real, d_loss_fake)

            # ---------------------
            #  Train Generator
            # ---------------------
            imgs2, noise2 = self.generate_data(X_train, batch_size)

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

        self.generator.save('AdversarialDenoisingGeneralNoBatch.h5')

    def predict(self, model, data):

        if isinstance(data, basestring):
            test = self.ExtractInfoMetadata(data, root2, xmipp.MDL_IMAGE_ORIGINAL,
                                            80,1)
        else:
            test = data
        prediction = model.predict(test)

        '''NoiseEnhanced = []
        for i, im in enumerate(prediction):
            NewNoise = (im + noise1[i]) / 2.0
            NoiseEnhanced.append(NewNoise)

        NoiseEnhanced = np.asarray(NoiseEnhanced).astype('float32')
        NoiseEnhanced = NoiseEnhanced.reshape(len(NoiseEnhanced),
                                              NoiseEnhanced.shape[1],
                                              NoiseEnhanced.shape[2], 1)
        predictEnhanced = model.predict(NoiseEnhanced)'''
        predictEnhanced = self.normalization(prediction, 1)

        return predictEnhanced

    def dice_coef(self, y_true, y_pred):
        intersection = K.sum(y_true * y_pred, axis=[1, 2])
        union = K.sum(y_true, axis=[1, 2]) + K.sum(y_pred, axis=[1, 2])
        return -K.mean((2. * intersection + 1) / (union + 1), axis=0)


    def save_imgs(self, X_train,epoch):
        filename = "denoise_%d.png"

        evaluate = self.generator.evaluate(self.noise,self.true)
        gen_imgs = self.generator.predict(self.noise)
        print "Validation =", evaluate

        # Rescale images 0 - 1
        gen_imgs = 0.5 * gen_imgs + 0.5

        fig, axs = plt.subplots(10, 3)
        cnt = 0
        for i in range(10):
            axs[i, 0].imshow(self.true[cnt, :, :, 0], cmap='gray')
            axs[i, 0].axis('off')
            axs[i, 1].imshow(self.noise[cnt, :, :, 0], cmap='gray')
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
    micrographs = gan.importMicrographs(micrographsCNB)
    X_train1 = gan.ExtractInfoMetadata(path1, root, xmipp.MDL_IMAGE_REF, 100,
      -1)
    X_train2 = gan.ExtractInfoMetadata(betagalRepro, root4,
                                      xmipp.MDL_IMAGE_REF, 100,
                                       -1)
    X_train3 = gan.ExtractInfoMetadata(path2, root2,
                                       xmipp.MDL_IMAGE_REF, 100,
                                       -1)

    X_train = np.concatenate((X_train1, X_train2, X_train3), axis=0)
    gan.train(X_train, epochs=100000,batch_size=32, save_interval=1000)
    '''particle = gan.ExtractInfoMetadata(reprojectionParticles, root2,
                                     xmipp.MDL_IMAGE, 0, 0, False)'''
    noise1 = gan.ExtractInfoMetadata(groel, root3,
                                     xmipp.MDL_IMAGE, 100, 2)
    model2 = load_model('AdversarialDenoisingGeneral.h5',
                        custom_objects={
                            'dice_coef': gan.dice_coef})
    projection = gan.ExtractInfoMetadata(groel, root3,
                                         xmipp.MDL_IMAGE_REF, 100, 1)
    model = load_model('AdversarialDenoisingGeneralNoBatch.h5',
                       custom_objects={
                    'dice_coef': gan.dice_coef})

    m = hashlib.md5('AdversarialDenoisingBetagal.h5').hexdigest()
    print m
    predict = gan.predict(model, noise1)
    predict2 = gan.predict(model2, noise1)

    noise1 = gan.normalization(noise1, 1)
    values = np.zeros((len(predict), 2))

    '''for k, micro in enumerate(micrographs):
        microDenoise = np.zeros((len(micrographs),micro.shape[0],
                                 micro.shape[1]))
        micro = gan.normalization(micro, 'mean', False)
        for i in range(0,micro.shape[0]-100, 100):
            for j in range(0,micro.shape[1]-100, 100):

                block = micro[i:100+i,j:100+j]
                #block = gan.normalization(block, 'mean', False)
                block = np.reshape(block, (1, 100, 100, 1))
                predictblock = gan.predict(model, block)
                predictblock = gan.normalization(np.squeeze(predictblock),1,
                                                 False)
                plt.imshow(np.squeeze(predictblock), cmap='gray')
                plt.show()
                microDenoise[k,i:100+i,j:100+j] = predictblock

        plt.subplot(1,2,1)
        plt.imshow(np.squeeze(microDenoise[k]), cmap='gray')
        plt.subplot(1,2,2)
        plt.imshow(np.squeeze(micro), cmap='gray')
        plt.show()'''

    cont = 1
    for i,im in enumerate(predict):
        #diff = im - projection[i]
        #m_norm = np.sum(abs(diff))
        #print m_norm
        correlation = pearsonr(im.ravel(), projection[i].ravel())
        correlation2 = pearsonr(predict2[i].ravel(), projection[i].ravel())
        correlation3 = pearsonr(predict2[i].ravel(), noise1[i].ravel())
        #correlation3 = match_template(projection[i], noise1[i])

        print correlation[0], correlation2[0], correlation3[0]
        plt.subplot(2,2,1)
        plt.imshow(np.squeeze(im), cmap='gray')
        plt.subplot(2,2, 2)
        plt.imshow(np.squeeze(noise1[i]), cmap='gray')
        plt.subplot(2, 2, 3)
        plt.imshow(np.squeeze(predict2[i]), cmap='gray')
        plt.subplot(2, 2, 4)
        plt.imshow(np.squeeze(projection[i]), cmap='gray')
        plt.show()

