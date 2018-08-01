import numpy as np
import keras
import tensorflow as tf

from keras.optimizers import Adam
from keras import backend as K
from keras.preprocessing.image import ImageDataGenerator
from keras.models import Sequential, Model, load_model
from keras.layers import Dense, Dropout, Flatten, BatchNormalization, Activation, LeakyReLU
from keras.layers import Conv2D, MaxPooling2D, ZeroPadding2D, Input, UpSampling2D, Lambda, Reshape, AveragePooling2D
from keras.optimizers import SGD
from keras.callbacks import TensorBoard, ModelCheckpoint, EarlyStopping
from keras.datasets import cifar10, mnist
from keras.regularizers import l2
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
import xmipp
import random
from itertools import izip
from keras.constraints import maxnorm

size = 81
batch_size = 32
seed = 42

def ExtractInfoMetadata(path, root, label, crop, size):
    metadata = xmipp.MetaData(path)
    Image = []
    I = xmipp.Image()
    for itemId in metadata:
        fn = metadata.getValue(label, itemId)
        n = fn.split('@')
        fn = n[0] + '@' + root + n[1]
        I.read(fn)
        Data = I.getData()
        Imresize = cv2.resize(Data,(size,size),interpolation=cv2.INTER_CUBIC)
        Image.append(
            cv2.normalize(np.asarray(Imresize[crop:-crop, crop:-crop]), None,
                          0.0, 1.0, cv2.NORM_MINMAX))

    Image = np.array(Image).astype('float')
    Image = Image.reshape(len(Image), Image.shape[1],Image.shape[2], 1)

    return Image

def validation_data(x,y):

    valx = []
    valy = []

    for i in range(100):
        k = np.random.randint(0, len(x))
        valx.append(x[k])
        valy.append(y[k])
    return np.asarray(valx).astype('float32'), np.asarray(valy).astype(
        'float32')

def generate_model():
    input_img = Input(shape=(NoisyImage.shape[1],NoisyImage.shape[2],1),
                      name='input')
    #auxiliary_input = Input(shape=(size,size,1), name='input2')
    x = Dropout(0.2)(input_img)
    #x = keras.layers.concatenate([input_img,auxiliary_input])
    x = Conv2D(batch_size, (15, 15), activation='linear',
               kernel_initializer='glorot_normal',
               padding='same')(
        x)
    x1 = Conv2D(batch_size, (7, 7), activation='linear',
                kernel_initializer='glorot_normal',
                padding='same')(
        x)
    x = keras.layers.subtract([x,x1])
    #x = Dropout(0.5)(x)
    x = Activation('relu')(x)
    x = BatchNormalization()(x)
    x = MaxPooling2D((2, 2), padding='same')(x)
    x = Conv2D(batch_size, (5, 5), activation='linear',
               kernel_initializer='glorot_normal',
               padding='same')(x)
    x1 = Conv2D(batch_size, (3, 3), activation='linear',
                kernel_initializer='glorot_normal',
                padding='same')(x)
    x = keras.layers.subtract([x,x1])
    #x = Dropout(0.5)(x)
    x = Activation('relu')(x)
    x = BatchNormalization()(x)
    '''x = MaxPooling2D((2, 2), padding='same')(x)
    x = Conv2D(batch_size, (5, 5), activation='linear',
               kernel_initializer='random_uniform', padding='same')(x)
    x1 = Conv2D(batch_size, (3, 3), activation='linear',
                kernel_initializer='random_uniform', padding='same')(x)
    x = keras.layers.subtract([x, x1])
    x = Activation('relu')(x)
    x = BatchNormalization()(x)'''
    encoded = MaxPooling2D((2, 2), padding='same', name='encoded')(x)

    # at this point the representation is (7, 7, 32)

    x = Conv2D(batch_size, (3, 3), activation='linear',
               kernel_initializer='glorot_normal',
               padding='same')(
     encoded)
    #x = Dropout(0.5)(x)
    x = Activation('relu')(x)
    x = BatchNormalization()(x)
    x = UpSampling2D((2, 2))(x)
    x = Conv2D(batch_size, (5, 5), activation='linear',
               kernel_initializer='glorot_normal',
               padding='same')(
        x)
    #x = Dropout(0.5)(x)
    x = Activation('relu')(x)
    #x = BatchNormalization()(x)
    x = UpSampling2D((2, 2))(x)
    '''x = Conv2D(batch_size, (5, 5), activation='linear',
               kernel_initializer='glorot_normal', padding='same')(
        x)
    x = Activation('relu')(x)
    #x = BatchNormalization()(x)
    x = UpSampling2D((2, 2))(x)'''
    #x = Dropout(0.25)(x)
    decoded = Conv2D(1, (9, 9), activation='sigmoid',
                     kernel_initializer='glorot_normal', padding='same',
                     name='decoded')(x)

    autoencoder = Model(input_img, decoded)#Model(inputs=[input_img,auxiliary_input],
      #outputs= [decoded,auxiliary_input])

    return autoencoder


path1 = '/home/javiermota/ScipionUserData/projects/CNBScipionCourse/Runs/023211_XmippProtGenerateReprojections/extra/anglesCont.xmd'
root = '/home/javiermota/ScipionUserData/projects/CNBScipionCourse/'
path2 = '/home/javiermota/ScipionUserData/projects/Frank10000_70S/Runs/006028_XmippProtGenerateReprojections/extra/anglesCont.xmd'
root2 = '/home/javiermota/ScipionUserData/projects/Frank10000_70S/'

NoisyImage = ExtractInfoMetadata(path1, root, xmipp.MDL_IMAGE, 20, 120)
Projection = ExtractInfoMetadata(path1, root, xmipp.MDL_IMAGE_REF, 20, 120)

NoisyImage2 = ExtractInfoMetadata(path2, root2, xmipp.MDL_IMAGE, 20, 120)
Projection2 = ExtractInfoMetadata(path2, root2, xmipp.MDL_IMAGE_REF, 20, 120)
'''
for i,im in enumerate(NoisyImage):

    plt.subplot(1,2,1)
    plt.imshow(np.squeeze(im))
    plt.gray()
    plt.subplot(1,2,2)
    plt.imshow(np.squeeze(Projection[i]))
    plt.gray()
    plt.show()'''


NoisyImage = np.concatenate((NoisyImage, NoisyImage2), axis=0)
Projection = np.concatenate((Projection, Projection2), axis=0)



x_train_data = ImageDataGenerator(horizontal_flip=True,
                                vertical_flip=True, width_shift_range=0.1,
                                height_shift_range=0.1,fill_mode='reflect')

y_train_data = ImageDataGenerator(horizontal_flip=True,
                                vertical_flip=True, width_shift_range=0.1,
                                height_shift_range=0.1, fill_mode='reflect')

x_train_data.fit(NoisyImage, augment=True, seed=seed)
y_train_data.fit(Projection, augment=True, seed=seed)

generator = izip(x_train_data.flow(NoisyImage,seed=seed,
                                      batch_size=batch_size),
                y_train_data.flow(Projection,seed=seed, batch_size=batch_size))
'''
for X_batch, y_batch in generator :
    for i in range(0, 9):
        plt.subplot(1,2,1)
        plt.imshow(X_batch[i].reshape(X_batch.shape[1], X_batch.shape[
            2]),cmap=plt.get_cmap('gray'))
        plt.subplot(1, 2, 2)
        plt.imshow(y_batch[i].reshape(y_batch.shape[1], y_batch.shape[
            2]), cmap=plt.get_cmap('gray'))
        plt.show()


    break'''

if os.path.isdir('/tmp/tb'):
    shutil.rmtree('/tmp/tb')

pathtest = '/home/javiermota/ScipionUserData/projects/CNBScipionCourse/Runs/011294_XmippProtCropResizeParticles/extra/output_images.xmd'
test = ExtractInfoMetadata(pathtest, root, xmipp.MDL_IMAGE, 20, 120)

xval, yval = validation_data(NoisyImage, Projection)

autoencoder = generate_model()
autoencoder.compile(optimizer=Adam(lr=0.0001),
                    loss='binary_crossentropy')
tf.summary.image("input",autoencoder.get_layer('input').output)
tf.summary.image("image",autoencoder.get_layer('decoded').output)
checkpointer = ModelCheckpoint(
    filepath='DenoisingParticlesPrueba.h5',
                               monitor='val_loss',
                               save_best_only=True)

'''autoencoder.fit(NoisyImage,Projection, batch_size,
                epochs=50,
                     validation_split=0.1,
                     callbacks=[TensorBoard(log_dir='/tmp/tb', histogram_freq=1,
                                            write_graph=False,
                                            write_images=False), checkpointer])'''

autoencoder.fit_generator(generator,steps_per_epoch=len(
    NoisyImage)/batch_size, epochs = 150,
                          validation_data=(xval, yval),
                          validation_steps=20,
                           callbacks=[TensorBoard(log_dir='/tmp/tb', histogram_freq=1,
                write_graph=False, write_images=False),checkpointer])
K.clear_session()


model = load_model('DenoisingParticlesPrueba.h5')
prediction = model.predict(NoisyImage)

for i,im in enumerate(prediction):
    plt.subplot(1,2,1)
    plt.imshow(np.squeeze(im))
    plt.gray()
    plt.subplot(1,2,2)
    plt.imshow(np.squeeze(NoisyImage[i]))
    plt.gray()
    plt.show()

ImageEnhanced = prediction + NoisyImage
ImageNormalize = []
for im in ImageEnhanced:
    ImageNormalize.append(cv2.normalize(np.asarray(im), None, 0.0, 1.0,
                                       cv2.NORM_MINMAX))
ImageNormalize = np.asarray(ImageNormalize).astype('float32')
ImageNormalize = ImageNormalize.reshape(len(ImageNormalize), ImageNormalize.shape[1], ImageNormalize.shape[2], 1)

x_train_data.fit(ImageNormalize, augment=True, seed=seed)

EnhanceGenerator = izip(x_train_data.flow(ImageNormalize,seed=seed,
                                      batch_size=batch_size),
                y_train_data.flow(Projection,seed=seed, batch_size=batch_size))
K.clear_session()

xval2, yval2 = validation_data(ImageNormalize, Projection)

autoencoder2 = generate_model()
autoencoder2.compile(optimizer=Adam(lr=0.0001),
                    loss='binary_crossentropy')
tf.summary.image("input",autoencoder2.get_layer('input').output)
tf.summary.image("image",autoencoder2.get_layer('decoded').output)

checkpointer2 = ModelCheckpoint(
    filepath='DenoisingParticlesEnhancedPrueba.h5',
                               monitor='val_loss',
                               save_best_only=True)
autoencoder2.fit_generator(EnhanceGenerator,steps_per_epoch=len(
    ImageNormalize)/batch_size, epochs = 50,
                          validation_data=(xval2, yval2),
                          validation_steps=20,
                           callbacks=[TensorBoard(log_dir='/tmp/tb', histogram_freq=1,
                write_graph=False, write_images=False),checkpointer2])

'''
model2 = load_model('DenoisingParticlesEnhanced.h5')
prediction2 = model2.predict(test)

plt.figure()
for i,im in enumerate(prediction2):
    plt.subplot(1,2,1)
    plt.imshow(np.squeeze(im))
    plt.gray()
    plt.subplot(1,2,2)
    plt.imshow(np.squeeze(ImageNormalize[i]))
    plt.show()
'''

