
import numpy as np
import time
from tensorflow.examples.tutorials.mnist import input_data

from keras.models import Sequential, Model, load_model
from keras.layers import Dense, Activation, Flatten, Reshape
from keras.layers import Conv2D, Conv2DTranspose, UpSampling2D, Input, \
    MaxPooling2D, Subtract, AveragePooling2D
from keras.layers import LeakyReLU, Dropout
from keras.layers import BatchNormalization
from keras.optimizers import Adam, RMSprop, SGD
from keras.datasets import mnist
import xmipp
import cv2
import tensorflow as tf
import keras
from keras.callbacks import History

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

path1 = '/home/javiermota/ScipionUserData/projects/CNBScipionCourse/Runs/023211_XmippProtGenerateReprojections/extra/anglesCont.xmd'
root = '/home/javiermota/ScipionUserData/projects/CNBScipionCourse/'
path2 = '/home/javiermota/rinchen3/ScipionUserData/projects/coss/Runs/036738_XmippProtGenerateReprojections/extra/anglesCont.xmd'
root2 = '/home/javiermota/rinchen3/ScipionUserData/projects/coss/'

class GAN():
    def __init__(self):
        self.img_rows = 60 #40
        self.img_cols = 60 #40
        self.channels = 1
        self.shape = self.img_rows*self.img_cols
        self.img_shape = (self.img_rows, self.img_cols, self.channels)

        optimizer = Adam(0.0002, 0.5)

        # Build and compile the discriminator
        self.discriminator = self.build_discriminator()
        self.discriminator.compile(loss='binary_crossentropy',
                                   optimizer=optimizer,
                                   metrics=['accuracy'])

        # Build and compile the generator
        self.generator = self.build_generator()
        self.generator.compile(loss='binary_crossentropy', optimizer=optimizer)

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
            Data = I.getData()
            Imresize = cv2.resize(Data, (size, size),
                                  interpolation=cv2.INTER_CUBIC)

            if norm == -1:
                Imnormalize = 2*(Imresize[crop:-crop, crop:-crop]-np.min(Imresize[
                                                                      crop:-crop,
                                                                      crop:-crop]))/\
                              (np.max(Imresize[crop:-crop, crop:-crop])-np.min(
                                  Imresize[crop:-crop, crop:-crop]))-1
            else:
                Imnormalize = (Imresize[crop:-crop, crop:-crop]-np.mean(
                    Imresize[crop:-crop, crop:-crop]))/np.std(Imresize[crop:-crop, crop:-crop])

            Image.append(Imnormalize)
            if cont > 2000:
                break
            cont += 1

        Image = np.array(Image).astype('float')
        Image = Image.reshape(len(Image), Image.shape[1], Image.shape[2], 1)

        return Image

    def build_generator(self):

        #noise_shape = (self.shape,)

        model = Sequential()

        model.add(Conv2D(32,(7,7),padding='same', input_shape=self.img_shape))
        model.add(BatchNormalization(momentum=0.8))
        model.add(LeakyReLU(alpha=0.2))
        model.add(AveragePooling2D(2,2))
        model.add(
            Conv2D(64, (5, 5), padding='same'))
        model.add(BatchNormalization(momentum=0.8))
        model.add(LeakyReLU(alpha=0.2))
        model.add(AveragePooling2D(2, 2))
        model.add(Conv2D(128, (5, 5), padding='same'))
        model.add(BatchNormalization(momentum=0.8))
        model.add(LeakyReLU(alpha=0.2))
        model.add(UpSampling2D(2))
        model.add(Conv2D(64, (5, 5), padding='same'))
        model.add(BatchNormalization(momentum=0.8))
        model.add(LeakyReLU(alpha=0.2))
        model.add(UpSampling2D(2))
        model.add(Conv2D(32, (5, 5), padding='same'))
        model.add(BatchNormalization(momentum=0.8))
        model.add(LeakyReLU(alpha=0.2))
        model.add(Conv2D(1, (3, 3), padding='same'))
        model.add(Activation('tanh'))
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

        model.summary()

        noise = Input(shape=self.img_shape)#Input(shape=noise_shape)
        img = model(noise)

        return Model(noise, img)


    def build_discriminator(self):

        img_shape = (self.img_rows, self.img_cols, self.channels)

        model = Sequential()

        model.add(Flatten(input_shape=img_shape))
        model.add(Dense(512))
        model.add(LeakyReLU(alpha=0.2))
        model.add(Dense(256))
        model.add(LeakyReLU(alpha=0.2))
        model.add(Dense(1, activation='sigmoid'))
        model.summary()

        img = Input(shape=img_shape)
        validity = model(img)

        return Model(img, validity)

    def train(self, epochs, batch_size=128, save_interval=50):

        # Load the dataset
        #(X_train, _), (_, _) = mnist.load_data()

        # Rescale -1 to 1
        #X_train = (X_train.astype(np.float32) - 127.5) / 127.5
        #X_train = np.expand_dims(X_train, axis=3)


        X_train1 = self.ExtractInfoMetadata(path1, root,xmipp.MDL_IMAGE_REF,
         5, 70)
        X_train2 = self.ExtractInfoMetadata(path2, root2, xmipp.MDL_IMAGE_REF,
                                            5, 70)
        self.X_train = X_train1#np.concatenate((X_train1, X_train2), axis=0)
        noise1 = self.ExtractInfoMetadata(path1, root,xmipp.MDL_IMAGE,
         5, 70, 1)
        noise2 = self.ExtractInfoMetadata(path2, root2, xmipp.MDL_IMAGE,
                                          5, 70, 1)
        self.noise = noise1#np.concatenate((noise1, noise2), axis=0)
        half_batch = int(batch_size / 2)

        for epoch in range(epochs):

            # ---------------------
            #  Train Discriminator
            # ---------------------

            # Select a random half batch of images
            idx = np.random.randint(0, self.X_train.shape[0], half_batch)
            imgs = self.X_train[idx]

            #noise7 = np.random.normal(0, 1, (half_batch, self.shape))
            noise1 = self.noise[idx]
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
            idy = np.random.randint(0, self.X_train.shape[0], batch_size)
            noise2 = self.noise[idy]
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

            # If at save interval => save generated image samples
            if epoch % save_interval == 0:
                self.save_imgs(epoch)

        self.generator.save('AdversarialDenoising.h5')

    def predict(self):

        model = load_model('AdversarialDenoising.h5')
        pathtest = '/home/javiermota/ScipionUserData/projects/CNBScipionCourse/Runs/011294_XmippProtCropResizeParticles/extra/output_images.xmd'
        test = self.ExtractInfoMetadata(pathtest, root, xmipp.MDL_IMAGE, 10,
                                        60,1)
        prediction = model.predict(test)

        for i, im in enumerate(prediction):
            plt.subplot(1, 2, 1)
            plt.imshow(np.squeeze(im))
            plt.gray()
            plt.subplot(1, 2, 2)
            plt.imshow(np.squeeze(test[i]))
            plt.gray()
            plt.show()

    def save_imgs(self, epoch):
        filename = "mnist_%d.png"
        r, c = 5, 5
        idx = np.random.randint(0, self.X_train.shape[0], 10)
        noise = self.noise[idx]#np.random.normal(0, 1, (r * c, self.shape))
        #noise2 = noise.reshape(len(noise), 1600)
        true = self.X_train[idx]
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

if __name__ == '__main__':
    gan = GAN()
    gan.train(epochs=50000, batch_size=32, save_interval=500)