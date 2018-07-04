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

train_data_dir = '/home/javiermota/Downloads/images/Train/cats/*.jpg'

noisyImage = []
for im in glob.glob(train_data_dir):
    image = cv2.imread(im,0)
    image = cv2.resize(image,(128,128))
    image.astype('float32')/255.
    noise = np.random.normal(loc=0.0, scale=1, size=np.shape(image))
    noisyImage.append(image+noise)

noisyImage = np.asarray(noisyImage)
noisyImage = np.reshape(noisyImage,(len(noisyImage),128,128,1))
