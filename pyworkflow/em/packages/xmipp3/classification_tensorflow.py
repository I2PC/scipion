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

if os.path.isdir('/tmp/tb'):
    shutil.rmtree('/tmp/tb')
size = 81
batch_size = 32
train_data_dir = '/home/javiermota/Downloads/images/Train/cats/*.jpg'


paths = []
for image in glob.glob('/home/javiermota/Frank10000_70S/particle_images'
                       '/*.jpg'):

    paths.append(image)

def data_val(paths):
    batch_img = []
    labels = []
    for i in range(100):
        k = np.random.randint(0,len(paths))
        image = cv2.imread(paths[k],0)
        #with mrcfile.open(paths[k]) as mrc:
            #for im in mrc.data:
        image = cv2.resize(image,(size,size))
        image = cv2.normalize(image.astype('float'), None, 0.0, 1.0,
                              cv2.NORM_MINMAX)
        image = np.reshape(image,(size,size,1))
        batch_img.append(image)
        labels.append(image)

    return batch_img, batch_img

noisyPatch, other = data_val(paths)
#noisyPatch = np.asarray(noisyPatch)
#plt.imshow(noisyPatch.reshape(135,135))
#plt.gray()
#plt.show()

#with mrcfile.open('/home/javiermota/Downloads/archive/10184/data
    # /17dec27a_aldolase_00003gr_00032sq_v02_00002hln_00002esn-a-DW.mrc') as mrc:
    #noisyPatch = mrc.data

#noisyPatch = noisyPatch[500:500+size,500:500+size]
#noisyPatch = cv2.normalize(noisyPatch.astype('float'), None, 0.0, 1.0,
# cv2.NORM_MINMAX)

noisyImages = []
images = []
for im in glob.glob(train_data_dir):
    try:
        k = np.random.randint(0, 100)
        image = cv2.imread(im,0)
        image = cv2.resize(image,(size,size))
        image = cv2.normalize(image.astype('float'), None, 0.0, 1.0,
         cv2.NORM_MINMAX)
        #image = np.asarray(0.1*image+np.reshape(noisyPatch[k],(size,size)))
        #plt.imshow(image)
        #plt.gray()
        #plt.show()
        noise = cv2.normalize(np.asarray(noisyPatch[k]).astype('float'),None,
                              0.0, 1.0, cv2.NORM_MINMAX)
        imageNoise = 0.2*image+noise
        #imageNoise = cv2.medianBlur(imageNoise.astype('float32'), 5)
        imageNoise = cv2.normalize(np.asarray(imageNoise).astype('float'),None,
                              0.0, 1.0, cv2.NORM_MINMAX)
        #cv2.imwrite('/home/javiermota/entropy.jpg', imageNoise)
        images.append(image)
        noisyImages.append(imageNoise)
    except:
        pass


def validation_data(x,y):

    valx = []
    valy = []

    for i in range(100):
        k = np.random.randint(0, len(x))
        valx.append(x[k])
        valy.append(y[k])
    return np.asarray(valx), np.asarray(valy)

noisyImages = np.asarray(noisyImages)
images = np.asarray(images)
noisyImages = noisyImages.astype('float32')
'''plt.imshow(noisyImages[0])
plt.gray()
plt.show()'''
images = images.astype('float32')

val_x, val_y = validation_data(noisyImages, images)

val_x = val_x.reshape((len(val_x),size,size,1))
val_y = val_y.reshape((len(val_y),size,size,1))
noisyImages = noisyImages.reshape((len(noisyImages), size,size,1))
images = images.reshape(len(images), size, size, 1)

train_data = ImageDataGenerator()
#train_data.fit(noisyImages)
train = train_data.flow(noisyImages,images, batch_size=batch_size)

input_img = Input(shape=(size,size,1),name='input')

x = Conv2D(batch_size, (21, 21), activation='linear', padding='same')(input_img)
x1 = Conv2D(batch_size, (17, 17), activation='linear', padding='same')(
    input_img)
x = keras.layers.subtract([x,x1])
x = Activation('relu')(x)
x = BatchNormalization()(x)
x = MaxPooling2D((3, 3), padding='same')(x)
x = Conv2D(batch_size, (15, 15), activation='linear', padding='same')(x)
x1 = Conv2D(batch_size, (11, 11), activation='linear', padding='same')(x)
x = keras.layers.subtract([x,x1])
x = Activation('relu')(x)
x = BatchNormalization()(x)
encoded = MaxPooling2D((3, 3), padding='same', name='encoded')(x)
#encoded = Conv2D(32,(3, 3), activation='relu',padding='same')(encoded)
#encoded = BatchNormalization()(encoded)

# at this point the representation is (7, 7, 32)

x = Conv2D(batch_size, (9, 9), activation='linear', padding='same')(encoded)
#x1 = Conv2D(batch_size, (11, 11), activation='linear', padding='same')(x)
#x = keras.layers.add([x,x1])
x = Activation('relu')(x)
x = BatchNormalization()(x)
x = UpSampling2D((3, 3))(x)
x = Conv2D(batch_size, (15, 15), activation='linear', padding='same')(x)
#x1 = Conv2D(batch_size, (19, 19), activation='linear', padding='same')(x)
#x = keras.layers.add([x,x1])
x = Activation('relu')(x)
x = BatchNormalization()(x)
x = UpSampling2D((3, 3))(x)
#x = Dropout(0.25)(x)
decoded = Conv2D(1, (15, 15), activation='sigmoid', padding='same',
                 name='decoded')(x)

autoencoder = Model(input_img, decoded)

autoencoder.compile(optimizer=Adam(lr=0.0001), loss='binary_crossentropy')
tf.summary.image("input",autoencoder.get_layer('input').output)
tf.summary.image("image",autoencoder.get_layer('decoded').output)
#for i in range(0,31):
    #tf.summary.image("image"+str(i),K.expand_dims(autoencoder.get_layer(
    #'decoded').output[:,:,:,i]))
encoded_layer = Model(inputs=autoencoder.input,outputs=autoencoder.get_layer('encoded').output)
#autoencoder.save('/home/javiermota/autoencoder_randomNoise.h5')
checkpointer = ModelCheckpoint(
    filepath='autoencoder_500.h5',
                               monitor='val_loss',
                               save_best_only=True)
autoencoder.fit_generator(train,steps_per_epoch=len(
    noisyImages)/batch_size, epochs = 500, validation_data=(val_x,val_y),
                           callbacks=[TensorBoard(log_dir='/tmp/tb', histogram_freq=1,
                write_graph=False, write_images=False),checkpointer])
'''
autoencoder.fit(noisyImages, images,
                epochs=150,
                batch_size=batch_size,
                shuffle=True,
                validation_split=0.1,
                callbacks=[TensorBoard(log_dir='/tmp/tb', histogram_freq=1,
                write_graph=False, write_images=False),checkpointer])'''

def data_gen(paths):
    while True:
        batch_img = []
        labels = []
        for i in range(batch_size):
            k = np.random.randint(0,len(paths))
            with mrcfile.open(paths[k]) as mrc:
                for im in mrc.data:
                    image = cv2.resize(im,(size,size))
                    image = cv2.normalize(image.astype('float'), None, 0.0, 1.0,
                                          cv2.NORM_MINMAX)
                    image = np.reshape(image,(size,size,1))
                    batch_img.append(image)
                    labels.append(image)

        yield np.asarray(batch_img).astype('float32'), np.asarray(
            batch_img).astype('float32')

def predict_data(paths):
    batch_img = []
    labels = []
    #for i in range(10):
    k = np.random.randint(0,len(paths))
    #with mrcfile.open(paths[k]) as mrc:
        #for im in mrc.data:
    image = cv2.imread(paths[0],0)
    image = cv2.resize(image,(size,size))
    image = cv2.normalize(image.astype('float'), None, 0.0, 1.0,
                          cv2.NORM_MINMAX)
    image = np.reshape(image,(size,size,1))
    batch_img.append(image)
    labels.append(image)

    return np.asarray(batch_img).astype('float32')

test_data = predict_data(paths)
model = load_model('autoencoder_500.h5')
prediction = model.predict(test_data)
print prediction.shape
fig= plt.figure()
fig.add_subplot(1,2,1)
plt.imshow(np.squeeze(test_data))
plt.gray()
fig.add_subplot(1,2,2)
plt.imshow(np.squeeze(prediction))
plt.gray()
plt.show()
for layers in model.layers[:-1]:
    layers.trainable = False

model.layers.pop()
X_val, y_val = data_val(paths)
X_val = np.asarray(X_val).astype('float32')
y_val = np.asarray(y_val).astype('float32')
micrographs = Sequential()
micrographs.add(model)
micrographs.add(Conv2D(1, (15, 15), activation='sigmoid', padding='same',
                 name='decoded'))
micrographs.compile(optimizer= Adam(lr=0.001),
                    loss='binary_crossentropy')
tf.summary.image("image",micrographs.get_layer('decoded').output)
'''
Early = EarlyStopping(monitor='val_loss', patience=2)
micrographs.fit(noisyImages, images,
                epochs=50,
                batch_size=batch_size,
                shuffle=True,
                validation_split=0.1,
                callbacks=[TensorBoard(log_dir='/tmp/tb', histogram_freq=1,
                                       write_graph=False,
                                       write_images=False), Early])

micrographs.compile(optimizer='sgd', loss='binary_crossentropy')
tf.summary.image("image2",micrographs.get_layer('decoded').output)
micrographs.fit(noisyImages, images,
                epochs=10,
                batch_size=batch_size,
                shuffle=True,
                validation_split=0.1,
                callbacks=[TensorBoard(log_dir='/tmp/tb', histogram_freq=1,
                                       write_graph=False,
                                       write_images=False)])


micrographs.fit_generator(data_gen(paths),steps_per_epoch=10, epochs=10,
                          validation_data=(X_val,y_val), validation_steps=10,
                          callbacks=[TensorBoard(log_dir='/tmp/tb',
                                                 histogram_freq=1,
                                                 write_graph=False,
                                                 write_images=False)],
                          use_multiprocessing=True)'''
