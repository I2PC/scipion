
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
import mrcfile
import glob
import matplotlib.image as mpimg
from keras.applications import VGG19
import cv2
import shutil
import os

#shutil.rmtree('/tmp/tb')
# dimensions of our images.
size = 75

train_data_dir = '/home/javiermota/Downloads/images/Train/cats/*.jpg'
validation_data_dir = '/home/javiermota/images_tensorflow/Validation'
test_data_dir = '/home/javiermota/images_tensorflow/Test'
nb_train_samples = 30
nb_validation_samples = 18
epochs = 50
batch_size = 32

with mrcfile.open('/home/javiermota/Downloads/archive/10184/data/17dec27a_aldolase_00003gr_00032sq_v02_00002hln_00002esn-a-DW.mrc') as mrc:
    noisyPatch = mrc.data




(x_train, _), (x_test, _) = mnist.load_data()
noisyPatch = noisyPatch[500:575,500:575]
noisyPatch = cv2.normalize(noisyPatch.astype('float'), None, 0.0, 1.0, cv2.NORM_MINMAX)

imagesTrain_resized = []
imagesTest_resized = []
imagesNoisy = []
for im in x_train:
    im_resize = cv2.resize(im,(size,size))
    im_resize = cv2.normalize(im_resize.astype('float'), None, 0.0, 1.0, cv2.NORM_MINMAX)
    imagesTrain_resized.append(im_resize)
    imagesNoisy.append(im_resize+noisyPatch)
    plt.imshow(im_resize+noisyPatch)
    plt.gray()
    plt.show()
for im in x_test:
    im_resizeTest = cv2.resize(im,(size,size))
    imagesTest_resized.append(im_resizeTest)
x_train = np.asarray(imagesTrain_resized)
x_test = np.asarray(imagesTest_resized)
x_train = x_train.astype('float32')
x_test = x_test.astype('float32')
x_train = x_train.reshape((len(x_train), size,size,1))
x_test = x_test.reshape((len(x_test), size,size,1))
print x_train.shape
print x_test.shape

noise_factor = 2
x_train_noisy = np.asarray(imagesNoisy)#x_train + noise_factor * np.random.normal(loc=0.0, scale=1.0, size=x_train.shape)
x_test_noisy = x_test + noise_factor * np.random.normal(loc=0.0, scale=1.0, size=x_test.shape)

'''x_train_noisy = np.clip(x_train_noisy, 0., 1.)
x_test_noisy = np.clip(x_test_noisy, 0., 1.)
plt.imshow(x_train_noisy[0].reshape(size,size))
plt.gray()
plt.show()'''

paths = []
for image in glob.glob('/home/javiermota/Downloads/archive/10184/data/*.mrc'):
    paths.append(image)

def data_gen(paths):
    while True:
        batch_img = []
        labels = []
        for i in range(batch_size):
            k = np.random.randint(0,len(paths))
            with mrcfile.open(paths[k]) as mrc:
                image = cv2.resize(mrc.data,(size,size))
                image = np.reshape(image,(size,size,1))
                batch_img.append(image)
                labels.append(image)

        yield np.asarray(batch_img).astype('float32')/255., np.asarray(labels)

model = load_model('/home/javiermota/autoencoder.h5')

for layers in model.layers[:-1]:
    layers.trainable = False
    print layers, layers.trainable

micrographs = Sequential()
micrographs.add(model)
micrographs.compile(optimizer='adam',loss='binary_crossentropy')
micrographs.fit_generator(data_gen(paths),steps_per_epoch=50, epochs=10,
                          callbacks=[TensorBoard(log_dir='/tmp/tb', histogram_freq=1, write_graph=False, write_images=True)])

'''noisyImage = []
for im in glob.glob(train_data_dir):
    image = cv2.imread(im,0)
    image.astype('float32')/255.
    noise = np.random.normal(loc=0.0, scale=1, size=np.shape(image))
    noiseimage = cv2.resize(image+noise,(128,128))
    noisyImage.append(np.ravel(noiseimage))
    #plt.imshow(noiseimage,cmap = plt.get_cmap('gray'))

noisyImage = np.asarray(noisyImage)
noisyImage = np.reshape(noisyImage,(len(noisyImage),128,128,1))'''


input_img = Input(shape=(size,size,1),name='input')

x = Conv2D(32, (15, 15), activation='relu', padding='same')(input_img)
x = MaxPooling2D((5, 5), padding='same')(x)
x = BatchNormalization()(x)
x = Conv2D(32, (7, 7), activation='relu', padding='same')(x)
encoded = MaxPooling2D((3, 3), padding='same', name='encoded')(x)
encoded = BatchNormalization()(encoded)

# at this point the representation is (7, 7, 32)

x = Conv2D(32, (7, 7), activation='relu', padding='same')(encoded)
x = UpSampling2D((3, 3))(x)
x = BatchNormalization()(x)
x = Conv2D(32, (15, 15), activation='relu', padding='same')(x)
x = UpSampling2D((5, 5))(x)
x = BatchNormalization()(x)
decoded = Conv2D(1, (15, 15), activation='sigmoid', padding='same', name='decoded')(x)

'''x = Dense(32,activation='relu')(input_img)
x = Dense(16,activation='relu')(x)
x = Dense(8,activation='relu', name='encoded')(x)
x = Dense(16,activation='relu')(x)
x = Dense(32,activation='relu')(x)
decoded = Dense(size,activation='sigmoid', name='decoded')(x)'''


autoencoder = Model(input_img, decoded)

autoencoder.compile(optimizer='adam', loss='binary_crossentropy')
tf.summary.image("image",autoencoder.get_layer('decoded').output)
'''for i in range(0,31):
    tf.summary.image("image"+str(i),K.expand_dims(autoencoder.get_layer('decoded').output[:,:,:,i]))'''
#autoencoder.fit_generator(data_gen(paths),steps_per_epoch=100,nb_epoch=50)
encoded_layer = Model(inputs=autoencoder.input,outputs=autoencoder.get_layer('encoded').output)
#autoencoder.save_weights('/home/javiermota/weights_improved.h5')
autoencoder.save('/home/javiermota/autoencoder.h5')
autoencoder.fit(x_train_noisy, x_train,
                epochs=10,
                batch_size=32,
                shuffle=True,
                validation_data=(x_test_noisy, x_test),
                callbacks=[TensorBoard(log_dir='/tmp/tb', histogram_freq=1, write_graph=False, write_images=True)])

#train = ImageDataGenerator().flow_from_directory(train_data_dir, target_size=(img_width,img_height),  classes=['views'], batch_size=batch_size)
#validation = ImageDataGenerator().flow_from_directory(validation_data_dir, target_size=(img_width,img_height), classes=['views'], batch_size=batch_size/2)
#test = ImageDataGenerator().flow_from_directory(test_data_dir, target_size=(img_width,img_height), classes=['views'], batch_size=batch_size)



'''
def getImages(path):
    for image in glob.glob(path):
        with mrcfile.open(image) as mrc:
            a = mrc.data
        yield a

images = getImages('/home/javiermota/CCTData/Extract/job057/Micrographs/*.mrcs')

for i in images:
    print i


vgg_conv = VGG19(weights='imagenet',include_top=False)

for layer in vgg_conv.layers[:-4]:
    layer.trainable = False
    print layer, layer.trainable

model = Sequential()

model.add(vgg_conv)
model.add(Flatten())
model.add(Dense(1024, activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(3, activation='softmax'))

num_classes = 10
(x_train, y_train), (x_test, y_test) = cifar10.load_data()

y_train = keras.utils.to_categorical(y_train, num_classes)
y_test = keras.utils.to_categorical(y_test, num_classes)
model = Sequential()
# input: 100x100 images with 3 channels -> (100, 100, 3) tensors.
# this applies 32 convolution filters of size 3x3 each.
model.add(Conv2D(32, (5, 5),padding='same',input_shape=x_train.shape[1:]))
model.add(BatchNormalization())
#model.add(MaxPooling2D(pool_size=(3, 3)))
model.add(Dropout(0.1))
model.add(ZeroPadding2D(padding=(5,5)))
model.add(Conv2D(32, (5, 5)))
model.add(BatchNormalization())
model.add(Activation('leakyrelu'))
model.add(MaxPooling2D(pool_size=(3, 3)))
model.add(Dropout(0.25))

model.add(ZeroPadding2D(padding=(3,3)))
model.add(Conv2D(64, (3, 3)))
model.add(BatchNormalization())
model.add(Dropout(0.1))
#model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(ZeroPadding2D(padding=(3,3)))
model.add(Conv2D(64, (3, 3)))
model.add(BatchNormalization())
model.add(Activation('relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))

model.add(Flatten())
model.add(Dense(256, activation='relu'))
model.add(BatchNormalization())
model.add(Dropout(0.5))
model.add(Dense(10, activation='softmax'))

#sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])

history = model.fit(x_train, y_train, validation_split=0.2, batch_size=32, epochs=50)
print history.history.keys()
score = model.evaluate(x_test, y_test, batch_size=32)
print score
prediction = model.predict_classes(x_test, batch_size=32)
ok = 0
for i, x in enumerate(prediction):
    if x == np.where(y_test[i]==np.max(y_test[i]))[0]:
        ok += 1

accuracy = ok/float(len(y_test))
print accuracy
plt.figure()
plt.plot(history.history['loss'])
plt.figure()
plt.plot(history.history['acc'])
plt.figure()
plt.plot(history.history['val_acc'])
plt.show()'''