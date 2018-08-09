'''
DCGAN on MNIST using Keras
Author: Rowel Atienza
Project: https://github.com/roatienza/Deep-Learning-Experiments
Dependencies: tensorflow 1.0 and keras 2.0
Usage: python3 dcgan_mnist.py
'''

import numpy as np
import time
from tensorflow.examples.tutorials.mnist import input_data

from keras.models import Sequential, Model
from keras.layers import Dense, Activation, Flatten, Reshape
from keras.layers import Conv2D, Conv2DTranspose, UpSampling2D, Input, \
    MaxPooling2D, Subtract, AveragePooling2D
from keras.layers import LeakyReLU, Dropout
from keras.layers import BatchNormalization
from keras.optimizers import Adam, RMSprop, SGD
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
path2 = '/home/javiermota/ScipionUserData/projects/Frank10000_70S/Runs/006028_XmippProtGenerateReprojections/extra/anglesCont.xmd'
root2 = '/home/javiermota/ScipionUserData/projects/Frank10000_70S/'

class ElapsedTimer(object):
    def __init__(self):
        self.start_time = time.time()
    def elapsed(self,sec):
        if sec < 60:
            return str(sec) + " sec"
        elif sec < (60 * 60):
            return str(sec / 60) + " min"
        else:
            return str(sec / (60 * 60)) + " hr"
    def elapsed_time(self):
        print("Elapsed: %s " % self.elapsed(time.time() - self.start_time) )

class DCGAN(object):
    def __init__(self, img_rows=28, img_cols=28, channel=1):

        self.img_rows = img_rows
        self.img_cols = img_cols
        self.channel = channel
        self.D = None   # discriminator
        self.G = None   # generator
        self.AM = None  # adversarial model
        self.DM = None  # discriminator model
        self.batch_size = 32


    def discriminator(self):
        if self.D:
            return self.D
        self.D = Sequential()
        depth = 64
        dropout = 0.4
        # In: 28 x 28 x 1, depth = 1
        # Out: 14 x 14 x 1, depth=64
        input_shape = (self.img_rows, self.img_cols, self.channel)
        self.D.add(Conv2D(depth*1, 5, strides=2, input_shape=input_shape,\
            padding='same'))
        self.D.add(LeakyReLU(alpha=0.2))
        self.D.add(Dropout(dropout))

        self.D.add(Conv2D(depth*2, 5, strides=2, padding='same'))
        self.D.add(LeakyReLU(alpha=0.2))
        self.D.add(Dropout(dropout))

        self.D.add(Conv2D(depth*4, 5, strides=2, padding='same'))
        self.D.add(LeakyReLU(alpha=0.2))
        self.D.add(Dropout(dropout))

        self.D.add(Conv2D(depth*8, 5, strides=1, padding='same'))
        self.D.add(LeakyReLU(alpha=0.2))
        self.D.add(Dropout(dropout))

        # Out: 1-dim probability
        self.D.add(Flatten())
        self.D.add(Dense(1))
        self.D.add(Activation('tanh'))
        self.D.summary()
        return self.D

    def generator(self):
        if self.G:
            return self.G
        #self.G = Sequential()
        '''dropout = 0.4
        depth = 64+64+64+64
        dim = 10
        # In: 100
        # Out: dim x dim x depth
        self.G.add(Dense(dim*dim*depth, input_dim=1600))
        self.G.add(BatchNormalization(momentum=0.9))
        self.G.add(Activation('relu'))
        self.G.add(Reshape((dim, dim, depth)))
        self.G.add(Dropout(dropout))

        # In: dim x dim x depth
        # Out: 2*dim x 2*dim x depth/2
        self.G.add(UpSampling2D())
        self.G.add(Conv2DTranspose(int(depth/2), 5, padding='same'))
        self.G.add(BatchNormalization(momentum=0.9))
        self.G.add(Activation('relu'))

        self.G.add(UpSampling2D())
        self.G.add(Conv2DTranspose(int(depth/4), 5, padding='same'))
        self.G.add(BatchNormalization(momentum=0.9))
        self.G.add(Activation('relu'))

        self.G.add(Conv2DTranspose(int(depth/8), 5, padding='same'))
        self.G.add(BatchNormalization(momentum=0.9))
        self.G.add(Activation('relu'))

        # Out: 28 x 28 x 1 grayscale image [0.0,1.0] per pix
        self.G.add(Conv2DTranspose(1, 5, padding='same'))
        self.G.add(Activation('sigmoid'))'''

        input_img = Input(shape=(self.img_rows, self.img_rows, 1),
                          name='input')
        # auxiliary_input = Input(shape=(size,size,1), name='input2')
        x = Conv2D(self.batch_size, (15, 15), activation='linear',
                   kernel_initializer='random_normal',
                   padding='same')(input_img)
        x1 = Conv2D(self.batch_size, (7, 7), activation='linear',
                    kernel_initializer='random_normal',
                    padding='same')(
            x)
        x = keras.layers.subtract([x, x1])
        x = LeakyReLU()(x)#Activation('relu')(x)
        x = BatchNormalization()(x)
        x = AveragePooling2D((2, 2), padding='same')(x)
        x = Conv2D(self.batch_size*2, (5, 5), activation='linear',
                   kernel_initializer='random_normal',
                   padding='same')(x)
        x1 = Conv2D(self.batch_size*2, (3, 3), activation='linear',
                    kernel_initializer='random_normal',
                    padding='same')(x)
        x = keras.layers.subtract([x, x1])
        x = LeakyReLU()(x)#Activation('relu')(x)
        x = BatchNormalization()(x)
        encoded = AveragePooling2D((2, 2), padding='same', name='encoded')(x)

        # at this point the representation is (7, 7, 32)

        x = Conv2DTranspose(64, kernel_size=5, strides=2,
                   kernel_initializer='random_normal',
                   padding='same')(encoded)
        x = LeakyReLU(alpha=0.2)(x)#Activation('relu')(x)
        x = BatchNormalization()(x)
        x = Dropout(0.4)(x)
        #x = UpSampling2D((2, 2))(x)
        x = Conv2DTranspose(32, kernel_size=5, strides=2,padding='same')(x)
        x = LeakyReLU(alpha=0.2)(x)#Activation('relu')(x)
        x = BatchNormalization()(x)
        #x = UpSampling2D((2, 2))(x)
        decoded = Conv2DTranspose(1, kernel_size=5, strides=1,padding='same',
                         name='decoded')(x)
        decoded = Activation('tanh')(decoded)

        self.G = Model(input_img,
                            decoded)  # Model(inputs=[input_img,auxiliary_input],

        self.G.summary()
        return self.G

    def discriminator_model(self):
        if self.DM:
            return self.DM
        optimizer = Adam(0.00002)#RMSprop(lr=0.000002,
        # decay=6e-8)
        self.DM = Sequential()
        self.DM.add(self.discriminator())
        self.DM.compile(loss='binary_crossentropy', optimizer=optimizer,\
            metrics=['accuracy'])
        return self.DM

    def adversarial_model(self):
        if self.AM:
            return self.AM
        optimizer = Adam(0.00001, 0.5) #RMSprop(lr=0.000001,
        # decay=3e-8)
        self.AM = Sequential()
        self.AM.add(self.generator())
        self.AM.add(self.discriminator())
        self.AM.compile(loss='binary_crossentropy', optimizer=optimizer,\
            metrics=['accuracy'])
        return self.AM

class MNIST_DCGAN(object):
    def __init__(self):

        '''self.x_train = input_data.read_data_sets("mnist",\
        	one_hot=True).train.images
        self.x_train = self.x_train.reshape(-1, self.img_rows,\
        	self.img_cols, 1).astype(np.float32)'''
        self.x_train = self.ExtractInfoMetadata(path1, root,
                                                xmipp.MDL_IMAGE_REF, 10, 80)

        self.img_rows = self.x_train.shape[1]
        self.img_cols = self.x_train.shape[2]
        self.channel = 1

        self.DCGAN = DCGAN(self.img_rows, self.img_cols)
        self.discriminator =  self.DCGAN.discriminator_model()
        self.adversarial = self.DCGAN.adversarial_model()
        self.generator = self.DCGAN.generator()

    def ExtractInfoMetadata(self, path, root, label, crop, size):
        metadata = xmipp.MetaData(path)
        Image = []
        I = xmipp.Image()
        for itemId in metadata:
            fn = metadata.getValue(label, itemId)
            n = fn.split('@')
            fn = n[0] + '@' + root + n[1]
            I.read(fn)
            Data = I.getData()
            Imresize = cv2.resize(Data, (size, size),
                                  interpolation=cv2.INTER_CUBIC)
            Imnormalize = 2*(Imresize[crop:-crop, crop:-crop]-np.min(Imresize[
                                                                  crop:-crop,
                                                                  crop:-crop]))/\
                          (np.max(Imresize[crop:-crop, crop:-crop])-np.min(
                              Imresize[crop:-crop, crop:-crop]))-1
            Image.append(Imnormalize)

        Image = np.array(Image).astype('float')
        Image = Image.reshape(len(Image), Image.shape[1], Image.shape[2], 1)

        return Image

    def train(self, train_steps=2000, batch_size=256, save_interval=0):
        noise_input = None
        self.noise2 = self.ExtractInfoMetadata(path1, root, xmipp.MDL_IMAGE,
                                          10, 80)
        if save_interval>0:
            #noise_input = np.random.uniform(0.0, 1.0, size=[16, 1600])
            n = np.random.randint(0, self.x_train.shape[0], size=batch_size)
            noise_input = self.noise2[n, :, :, :]

        for i in range(train_steps):
            k = np.random.randint(0,self.x_train.shape[0], size=batch_size)
            images_train = self.x_train[k, :, :, :]
            noise = np.random.uniform(-1.0, 1.0, size=[batch_size, 100])
            noise3 = self.noise2[k,:,:,:]
            #noise3 = noise3.reshape(-1, noise2.shape[1]*noise2.shape[2])
            #images_fake = self.generator.predict(noise)
            images_fake = self.generator.predict(noise3)
            '''plt.subplot(1,2,1)
            plt.imshow(np.squeeze(noise3[0]))
            plt.gray()
            plt.subplot(1, 2, 2)
            plt.imshow(np.squeeze(images_fake[0]))
            plt.gray()
            plt.show()'''
            x = np.concatenate((images_train, images_fake))
            y = np.ones([2*batch_size, 1])
            y[batch_size:, :] = 0
            d_loss = self.discriminator.train_on_batch(x, y)

            y = np.ones([batch_size, 1])
            #noise = np.random.uniform(-1.0, 1.0, size=[batch_size, 100])
            #j = np.random.randint(0, self.x_train.shape[0], size=batch_size)
            #noise3 = self.noise2[j, :, :, :]
            a_loss = self.adversarial.train_on_batch(noise3, y)
            log_mesg = "%d: [D loss: %f, acc: %f]" % (i, d_loss[0], d_loss[1])
            log_mesg = "%s  [A loss: %f, acc: %f]" % (log_mesg, a_loss[0], a_loss[1])
            print(log_mesg)
            if save_interval>0:
                if (i+1)%save_interval==0:
                    self.plot_images(save2file=True, samples=noise_input.shape[0],\
                        noise=noise_input, step=(i+1))
        self.generator.save('AdversarialDenoising.h5')

    def write_log(callback, names, logs, batch_no):
        for name, value in zip(names, logs):
            summary = tf.Summary()
            summary_value = summary.value.add()
            summary_value.simple_value = value
            summary_value.tag = name
            callback.writer.add_summary(summary, batch_no)
            callback.writer.flush()

    def plot_images(self, save2file=False, fake=True, samples=16, noise=None, step=0):
        filename = 'mnist.png'
        if fake:
            if noise is None:
                #noise = np.random.uniform(-1.0, 1.0, size=[samples, 1600])
                n = np.random.randint(0, self.x_train.shape[0], size=32)
                noise = self.noise2[n, :, :, :]
            else:
                filename = "mnist_%d.png" % step
            images = self.generator.predict(noise)
        else:
            i = np.random.randint(0, self.x_train.shape[0], samples)
            images = self.x_train[i, :, :, :]

        plt.figure(figsize=(10,10))
        for i in range(16):
            plt.subplot(4, 4, i+1)
            image = images[i, :, :, :]
            image = np.reshape(image, [self.img_rows, self.img_cols])
            plt.imshow(image, cmap='gray')
            plt.axis('off')
        plt.tight_layout()
        if save2file:
            plt.savefig(filename)
            plt.close('all')
        else:
            plt.show()

if __name__ == '__main__':
    mnist_dcgan = MNIST_DCGAN()
    timer = ElapsedTimer()
    mnist_dcgan.train(train_steps=2000, batch_size=32, save_interval=250)
    timer.elapsed_time()
    mnist_dcgan.plot_images(fake=True)
    mnist_dcgan.plot_images(fake=False, save2file=True)