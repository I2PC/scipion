import numpy as np
import keras
import tensorflow as tf
from keras import backend as K
from keras.preprocessing.image import ImageDataGenerator
from keras.models import Sequential, Model, load_model
from keras.layers import Dense, Dropout, Flatten, BatchNormalization, Activation, LeakyReLU
from keras.layers import Conv2D, MaxPooling2D, ZeroPadding2D, Input, UpSampling2D, Lambda, Reshape, AveragePooling2D
from keras.optimizers import SGD
from keras.callbacks import TensorBoard
from keras.datasets import cifar10, mnist
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import mrcfile
import glob
from keras.applications import VGG19
import cv2
import shutil
import os

shutil.rmtree('/tmp/tb')
size = 75
train_data_dir = '/home/javiermota/Downloads/images/Train/cats/*.jpg'

with mrcfile.open('/home/javiermota/Downloads/archive/10184/data/17dec27a_aldolase_00003gr_00032sq_v02_00002hln_00002esn-a-DW.mrc') as mrc:
    noisyPatch = mrc.data

noisyPatch = noisyPatch[500:575,500:575]
noisyPatch = cv2.normalize(noisyPatch.astype('float'), None, 0.0, 1.0, cv2.NORM_MINMAX)

noisyImages = []
images = []
for im in glob.glob(train_data_dir):
    try:
        image = cv2.imread(im,0)
        image = cv2.resize(image,(size,size))
        image = cv2.normalize(image.astype('float'), None, 0.0, 1.0, cv2.NORM_MINMAX)
        #plt.imshow(0.2*image+noisyPatch)
        #plt.gray()
        #plt.show()
        images.append(image)
        noisyImages.append(0.2*image+noisyPatch)
    except:
        print im

noisyImages = np.asarray(noisyImages)
images = np.asarray(images)
noisyImages = noisyImages.astype('float32')
images = images.astype('float32')
noisyImages = noisyImages.reshape((len(noisyImages), size,size,1))
images = images.reshape(len(images), size, size, 1)

input_img = Input(shape=(size,size,1),name='input')

x = Conv2D(64, (21, 21), activation='relu', padding='same')(input_img)
x = MaxPooling2D((5, 5), padding='same')(x)
x = BatchNormalization()(x)
x = Conv2D(64, (15, 15), activation='relu', padding='same')(x)
encoded = MaxPooling2D((3, 3), padding='same', name='encoded')(x)
#encoded = Conv2D(32,(3, 3), activation='relu',padding='same')(encoded)
#encoded = BatchNormalization()(encoded)

# at this point the representation is (7, 7, 32)

x = Conv2D(32, (3, 3), activation='relu', padding='same')(encoded)
x = UpSampling2D((3, 3))(x)
x = BatchNormalization()(x)
x = Conv2D(32, (5, 5), activation='relu', padding='same')(x)
x = UpSampling2D((5, 5))(x)
x = BatchNormalization()(x)
decoded = Conv2D(1, (5, 5), activation='sigmoid', padding='same', name='decoded')(x)

autoencoder = Model(input_img, decoded)

autoencoder.compile(optimizer='adam', loss='binary_crossentropy')
tf.summary.image("image",autoencoder.get_layer('decoded').output)
'''for i in range(0,31):
    tf.summary.image("image"+str(i),K.expand_dims(autoencoder.get_layer('decoded').output[:,:,:,i]))'''
#autoencoder.fit_generator(data_gen(paths),steps_per_epoch=100,nb_epoch=50)
encoded_layer = Model(inputs=autoencoder.input,outputs=autoencoder.get_layer('encoded').output)
#autoencoder.save_weights('/home/javiermota/weights_improved.h5')
autoencoder.save('/home/javiermota/autoencoder.h5')
autoencoder.fit(noisyImages, images,
                epochs=100,
                batch_size=32,
                shuffle=True,
                validation_split=0.1,
                callbacks=[TensorBoard(log_dir='/tmp/tb', histogram_freq=1, write_graph=False, write_images=True)])

