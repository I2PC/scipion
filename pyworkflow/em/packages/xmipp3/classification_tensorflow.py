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


#from pyworkflow.em.metadata.utils import iterRows
#import pyworkflow.em.metadata as md

(x_train, y_train), (x_test, y_test) = cifar10.load_data()
metadata = xmipp.MetaData('/home/javiermota/ScipionUserData/projects'
                      '/CNBScipionCourse/Runs/023211_XmippProtGenerateReprojections/extra/anglesCont.xmd')
Inoise=xmipp.Image()
Iproj=xmipp.Image()
NoisyImage = []
Projection = []
root = '/home/javiermota/ScipionUserData/projects/CNBScipionCourse/'
for itemId in metadata:
    fn = metadata.getValue(xmipp.MDL_IMAGE, itemId)
    n = fn.split('@')
    fnNoise = n[0]+'@'+root+n[1]
    Inoise.read(fnNoise)
    NoisyImage.append(cv2.normalize(np.asarray(Inoise.getData()), None, 0.0, \
                                                                       1.0, cv2.NORM_MINMAX))
    #plt.imshow(NoisyImage[0])
    #plt.gray()
    #plt.show()
    fp = metadata.getValue(xmipp.MDL_IMAGE_REF, itemId)
    p = fp.split('@')
    fnProj = p[0]+'@'+root+p[1]
    Iproj.read(fnProj)
    Projection.append(cv2.normalize(np.asarray(Iproj.getData()), None, 0.0, \
                                                                       1.0, cv2.NORM_MINMAX))

NoisyImage = np.asarray(NoisyImage).astype('float32')
NoisyImage = NoisyImage.reshape(len(NoisyImage), NoisyImage.shape[1], NoisyImage.shape[2], 1)
Projection = np.asarray(Projection).astype('float32')
Projection = Projection.reshape(len(Projection),Projection.shape[1], Projection.shape[2], 1)

train_data = ImageDataGenerator(shear_range=20, horizontal_flip=True,
                                vertical_flip=True, width_shift_range=0.2,
                                height_shift_range=0.2)
train_data.fit(NoisyImage)

NoisyImage = train_data.flow(NoisyImage, Projection, batch_size=32)

if os.path.isdir('/tmp/tb'):
    shutil.rmtree('/tmp/tb')
size = 81
batch_size = 32
#train_data_dir = '/home/javiermota/Downloads/images/Train/cats/*.jpg'
'''
metadata2 = xmipp.MetaData('/home/javiermota/ScipionUserData/projects'
                          '/CNBScipionCourse/Runs/011294_XmippProtCropResizeParticles/extra/output_images.xmd')

Itest = xmipp.Image()
test = []
for itemId in metadata2:
    fn = metadata2.getValue(xmipp.MDL_IMAGE, itemId)
    n = fn.split('@')
    fntest = n[0]+'@'+root+n[1]
    Itest.read(fntest)
    test.append(cv2.normalize(np.asarray(Itest.getData()), None, 0.0, \
                                                                       1.0, cv2.NORM_MINMAX))

test = np.asarray(test).astype('float32')
test = test.reshape(len(test), test.shape[1], test.shape[2], 1)
model = load_model('DenoisingParticles0.h5')
prediction = model.predict(test)

plt.figure()
for i,im in enumerate(prediction):
    plt.subplot(1,2,1)
    plt.imshow(np.squeeze(im))
    plt.gray()
    plt.subplot(1,2,2)
    plt.imshow(np.squeeze(test[i]))
    plt.show()
'''
'''
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
        image = cv2.resize(image,(size,size), interpolation=cv2.INTER_CUBIC)
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
        image = cv2.resize(image,(size,size),interpolation=cv2.INTER_CUBIC)
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
plt.imshow(noisyImages[0])
plt.gray()
plt.show()
images = images.astype('float32')

val_x, val_y = validation_data(noisyImages, images)

val_x = val_x.reshape((len(val_x),size,size,1))
val_y = val_y.reshape((len(val_y),size,size,1))
noisyImages = noisyImages.reshape((len(noisyImages), size,size,1))
images = images.reshape(len(images), size, size, 1)

train_data = ImageDataGenerator()
#train_data.fit(noisyImages)
train = train_data.flow(noisyImages,images, batch_size=batch_size)'''
def generate_model():
    input_img = Input(shape=(NoisyImage.shape[1],NoisyImage.shape[2],1),
                      name='input')
    auxiliary_input = Input(shape=(size,size,1), name='input2')

    #x = keras.layers.concatenate([input_img,auxiliary_input])
    x = Conv2D(batch_size, (15, 15), activation='linear',
               kernel_initializer='random_uniform', padding='same')(
        input_img)
    x1 = Conv2D(batch_size, (7, 7), activation='linear',
                kernel_initializer='random_uniform', padding='same')(
        input_img)
    x = keras.layers.subtract([x,x1])
    x = Activation('relu')(x)
    x = BatchNormalization()(x)
    x = MaxPooling2D((2, 2), padding='same')(x)
    x = Conv2D(batch_size, (5, 5), activation='linear',
               kernel_initializer='random_uniform', padding='same')(x)
    x1 = Conv2D(batch_size, (3, 3), activation='linear',
                kernel_initializer='random_uniform', padding='same')(x)
    x = keras.layers.subtract([x,x1])
    x = Activation('relu')(x)
    x = BatchNormalization()(x)
    x = MaxPooling2D((2, 2), padding='same')(x)
    x = Conv2D(batch_size, (5, 5), activation='linear',
               kernel_initializer='random_uniform', padding='same')(x)
    x1 = Conv2D(batch_size, (3, 3), activation='linear',
                kernel_initializer='random_uniform', padding='same')(x)
    x = keras.layers.subtract([x, x1])
    x = Activation('relu')(x)
    x = BatchNormalization()(x)
    encoded = MaxPooling2D((2, 2), padding='same', name='encoded')(x)
    #encoded = Conv2D(32,(3, 3), activation='relu',padding='same')(encoded)
    #encoded = BatchNormalization()(encoded)

    # at this point the representation is (7, 7, 32)

    x = Conv2D(batch_size, (3, 3), activation='linear',
               kernel_initializer='random_uniform', padding='same')(
     encoded)
    #x1 = Conv2D(batch_size, (11, 11), activation='linear', padding='same')(x)
    #x = keras.layers.add([x,x1])
    x = Activation('relu')(x)
    x = BatchNormalization()(x)
    x = UpSampling2D((2, 2))(x)
    x = Conv2D(batch_size, (5, 5), activation='linear',
               kernel_initializer='random_uniform', padding='same')(
        x)
    #x1 = Conv2D(batch_size, (19, 19), activation='linear', padding='same')(x)
    #x = keras.layers.add([x,x1])
    x = Activation('relu')(x)
    x = BatchNormalization()(x)
    x = UpSampling2D((2, 2))(x)
    x = Conv2D(batch_size, (5, 5), activation='linear',
               kernel_initializer='random_uniform', padding='same')(
        x)
    x = Activation('relu')(x)
    x = BatchNormalization()(x)
    x = UpSampling2D((2, 2))(x)
    #x = Dropout(0.25)(x)
    decoded = Conv2D(1, (9, 9), activation='sigmoid',
                     kernel_initializer='random_uniform', padding='same',
                     name='decoded')(x)

    autoencoder = Model(input_img, decoded)#Model(inputs=[input_img,auxiliary_input],
      #outputs= [decoded,auxiliary_input])

    return autoencoder

def predict_data(paths):
    batch_img = []
    labels = []
    for i in paths:
    #k = np.random.randint(0,len(paths))
    #with mrcfile.open(paths[k]) as mrc:
        #for im in mrc.data:
        image = cv2.imread(i,0)
        image = cv2.resize(image,(size,size),interpolation=cv2.INTER_CUBIC)
        image = cv2.normalize(image.astype('float'), None, 0.0, 1.0,
                          cv2.NORM_MINMAX)
        image = np.reshape(image,(size,size,1))
        batch_img.append(image)
        labels.append(image)

    return np.asarray(batch_img).astype('float32')


autoencoder = generate_model()
autoencoder.compile(optimizer=Adam(lr=0.0001),
                    loss='binary_crossentropy')
tf.summary.image("input",autoencoder.get_layer('input').output)
tf.summary.image("image",autoencoder.get_layer('decoded').output)
tf.summary.image("imageOri", autoencoder.get_layer('input2').output)
#for i in range(0,31):
    #tf.summary.image("image"+str(i),K.expand_dims(autoencoder.get_layer(
    #'decoded').output[:,:,:,i]))
#encoded_layer = Model(inputs=autoencoder.input,
# outputs=autoencoder.get_layer('encoded').output)
#autoencoder.save('/home/javiermota/autoencoder_randomNoise.h5')
checkpointer = ModelCheckpoint(
    filepath='DenoisingParticles0.h5',
                               monitor='val_loss',
                               save_best_only=True)
autoencoder.fit([NoisyImage, Projection],[Projection, Projection], batch_size,
                epochs=50,
                     validation_split=0.1,
                     callbacks=[TensorBoard(log_dir='/tmp/tb', histogram_freq=1,
                                            write_graph=False,
                                            write_images=False), checkpointer])
'''autoencoder.fit_generator(train,steps_per_epoch=len(
    noisyImages)/batch_size, epochs = 50, validation_data=(val_x,val_y),
                           callbacks=[TensorBoard(log_dir='/tmp/tb', histogram_freq=1,
                write_graph=False, write_images=False),checkpointer])'''
K.clear_session()

'''
for i in range(10):
    if os.path.isdir('/tmp/tb'):
        shutil.rmtree('/tmp/tb')
    if i == 0:
        predictData = noisyImages
    else:
        predictData = noisyImages
    #model = load_model('autoencoder%i.h5'%i)
    if i == 0:
        predict = model.predict(np.asarray(predictData).astype('float32'))
        train2 = predict + noisyImages
        train2 = cv2.normalize(train2.astype('float32'), None, 0.0, 1.0,
                               cv2.NORM_MINMAX)
    else:
        predict = model.predict([np.asarray(predictData).astype('float32'),
                                 np.asarray(predictData).astype('float32')])
        train2 = predict[0] + noisyImages
        train2 = cv2.normalize(train2.astype('float32'), None, 0.0, 1.0,
                               cv2.NORM_MINMAX)


    autoencoder2 = generate_model()
    autoencoder2.compile(optimizer=Adam(lr=0.0001), loss='binary_crossentropy')
    tf.summary.image("input2", autoencoder2.get_layer('input').output)
    tf.summary.image("image2", autoencoder2.get_layer('decoded').output)
    tf.summary.image("imageOri", autoencoder2.get_layer('input2').output)
    checkpointer = ModelCheckpoint(
        filepath='autoencoder%i.h5'%(i+1),
        monitor='val_loss',
        save_best_only=True)
    autoencoder2.fit([predictData, images], [images, images], batch_size, epochs=50,
                     validation_split=0.2,
                     callbacks=[TensorBoard(log_dir='/tmp/tb', histogram_freq=1,
                                            write_graph=False,
                                            write_images=False), checkpointer])
    K.clear_session()
'''
'''
trainF = []
for i in range(9):

    model = load_model('autoencoder%i.h5'%(i+1))
    train = noisyImages
    if i == 0:
        predict = model.predict(np.asarray(noisyImages).astype('float32'))
        train = predict + noisyImages
        train = cv2.normalize(train.astype('float32'), None, 0.0, 1.0,
                               cv2.NORM_MINMAX)

        trainF = predict + noisyImages
    else:
    predict = model.predict([np.asarray(train).astype('float32'),
                             np.asarray(train).astype('float32')])
    predict = predict[0]
    #train = predict + noisyImages
    #train = cv2.normalize(train.astype('float32'), None, 0.0, 1.0,
    #                       cv2.NORM_MINMAX)
    trainF.append(predict)
    #plt.figure()
    #plt.subplot(1,2,1)
    #plt.imshow(np.squeeze(train[i]))
    #plt.gray()
    #plt.subplot(1,2,2)
    #plt.imshow(np.squeeze(images[i]))
    #plt.gray()
    #plt.show()

trainF = np.sum(trainF,axis=0)/(10.0)

for im in trainF:

    plt.imshow(np.squeeze(im))
    plt.gray()
    plt.show()


'''
'''
model2 = load_model('autoencoder2.h5')
predict2 = model2.predict([np.asarray(train2).astype('float32'), np.asarray(
    train2).astype('float32')])
sum = predict + predict2[0]
#plt.imshow(np.squeeze(sum[0]))
#plt.gray()
#plt.show()
train3 = (sum + noisyImages)/3
train3 = cv2.normalize(train3.astype('float32'), None, 0.0, 1.0,
                       cv2.NORM_MINMAX)


autoencoder.fit(noisyImages, images,
                epochs=150,
                batch_size=batch_size,
                shuffle=True,
                validation_split=0.1,
                callbacks=[TensorBoard(log_dir='/tmp/tb', histogram_freq=1,
                write_graph=False, write_images=False),checkpointer])

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


test_data = predict_data(paths)
model = load_model('autoencoder_multiply.h5')
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

'''
for i,image in enumerate(NoisyImage):
    plt.imshow(image[0:10,:])
    plt.gray()
    plt.show()

    x = random.randint(0,20)
    y = random.randint(0,20)
    M = np.float32([[1, 0, x], [0, 1, y]])
    imageShift = cv2.warpAffine(image,M, (image.shape[0], image.shape[1]))
    plt.imshow(imageShift)
    plt.gray()
    plt.show()
    patchx = image[:,-x:]
    patchy = image[-y:,:]
    imageShift[:,0:x] = patchx[:,:-x]
    imageShift[0:y,:] = patchy[:-y,:]
    plt.figure()
    plt.imshow(imageShift)
    plt.gray()
    plt.show()
'''