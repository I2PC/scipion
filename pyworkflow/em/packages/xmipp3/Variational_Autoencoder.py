'''Example of VAE on MNIST dataset using MLP
The VAE has a modular design. The encoder, decoder and VAE
are 3 models that share weights. After training the VAE model,
the encoder can be used to  generate latent vectors.
The decoder can be used to generate MNIST digits by sampling the
latent vector from a Gaussian distribution with mean=0 and std=1.
# Reference
[1] Kingma, Diederik P., and Max Welling.
"Auto-encoding variational bayes."
https://arxiv.org/abs/1312.6114
'''

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from keras.layers import Lambda, Input, Dense
from keras.models import Model
from keras.datasets import mnist
from keras.losses import mse, binary_crossentropy
from keras.utils import plot_model
from keras import backend as K
from keras.callbacks import TensorBoard

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import numpy as np
import argparse
import os
import glob
import cv2
import shutil
import tensorflow as tf
import xmipp

if os.path.isdir('/tmp/tb'):
        shutil.rmtree('/tmp/tb')
# reparameterization trick
# instead of sampling from Q(z|X), sample eps = N(0,I)
# z = z_mean + sqrt(var)*eps
def sampling(args):
    """Reparameterization trick by sampling fr an isotropic unit Gaussian.
    # Arguments:
        args (tensor): mean and log of variance of Q(z|X)
    # Returns:
        z (tensor): sampled latent vector
    """

    z_mean, z_log_var = args
    batch = K.shape(z_mean)[0]
    dim = K.int_shape(z_mean)[1]
    # by default, random_normal has mean=0 and std=1.0
    epsilon = K.random_normal(shape=(batch, dim))
    return z_mean + K.exp(0.5 * z_log_var) * epsilon


def plot_results(models,
                 data,
                 batch_size=128,
                 model_name="vae_mnist"):
    """Plots labels and MNIST digits as function of 2-dim latent vector
    # Arguments:
        models (tuple): encoder and decoder models
        data (tuple): test data and label
        batch_size (int): prediction batch size
        model_name (string): which model is using this function
    """

    encoder, decoder = models
    x_test, y_test = data
    os.makedirs(model_name, exist_ok=True)

    filename = os.path.join(model_name, "vae_mean.png")
    # display a 2D plot of the digit cfrom keras.callbacks import TensorBoardlasses in the latent space
    z_mean, _, _ = encoder.predict(x_test,
                                   batch_size=batch_size)
    plt.figure(figsize=(12, 10))
    plt.scatter(z_mean[:, 0], z_mean[:, 1], c=y_test)
    plt.colorbar()
    plt.xlabel("z[0]")
    plt.ylabel("z[1]")
    plt.savefig(filename)
    plt.show()

    filename = os.path.join(model_name, "digits_over_latent.png")
    # display a 30x30 2D manifold of digits
    n = 30
    digit_size = 28
    figure = np.zeros((digit_size * n, digit_size * n))
    # linearly spaced coordinates corresponding to the 2D plot
    # of digit classes in the latent space
    grid_x = np.linspace(-4, 4, n)
    grid_y = np.linspace(-4, 4, n)[::-1]

    for i, yi in enumerate(grid_y):
        for j, xi in enumerate(grid_x):
            z_sample = np.array([[xi, yi]])
            x_decoded = decoder.predict(z_sample)
            digit = x_decoded[0].reshape(digit_size, digit_size)
            figure[i * digit_size: (i + 1) * digit_size,
                   j * digit_size: (j + 1) * digit_size] = digit

    plt.figure(figsize=(10, 10))
    start_range = digit_size // 2
    end_range = n * digit_size + start_range + 1
    pixel_range = np.arange(start_range, end_range, digit_size)
    sample_range_x = np.round(grid_x, 1)
    sample_range_y = np.round(grid_y, 1)
    plt.xticks(pixel_range, sample_range_x)
    plt.yticks(pixel_range, sample_range_y)
    plt.xlabel("z[0]")
    plt.ylabel("z[1]")
    plt.imshow(figure, cmap='Greys_r')
    plt.savefig(filename)
    plt.show()


# MNIST dataset
#(x_train, y_train), (x_test, y_test) = mnist.load_data()
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


path2 = '/home/javiermota/ScipionUserData/projects/Frank10000_70S/Runs/006028_XmippProtGenerateReprojections/extra/anglesCont.xmd'
root2 = '/home/javiermota/ScipionUserData/projects/Frank10000_70S/'

x_train = ExtractInfoMetadata(path2, root2, xmipp.MDL_IMAGE, 20, 120)
x_test = ExtractInfoMetadata(path2, root2, xmipp.MDL_IMAGE_REF, 20, 120)

'''for im in glob.glob(train_data_dir):
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
        x_test.append(image)
        x_train.append(imageNoise)
    except:
        pass'''

x_train = np.asarray(x_train).astype('float32')
x_test = np.asarray(x_test).astype('float32')
image_size = x_train.shape[1]
original_dim = x_train.shape[1] * x_train.shape[1]
x_train = np.reshape(x_train, [-1, original_dim])
x_test = np.reshape(x_test, [-1, original_dim])
#x_train = x_train.astype('float32') / 255
#x_test = x_test.astype('float32') / 255

# network parameters
input_shape = (original_dim, )
intermediate_dim = 512
batch_size = 32
latent_dim = 5
epochs = 20

# VAE model = encoder + decoder
# build encoder model
inputs = Input(shape=input_shape, name='encoder_input')
x = Dense(intermediate_dim, activation='relu')(inputs)
z_mean = Dense(latent_dim, name='z_mean')(x)
z_log_var = Dense(latent_dim, name='z_log_var')(x)

# use reparameterization trick to push the sampling out as input
# note that "output_shape" isn't necessary with the TensorFlow backend
z = Lambda(sampling, output_shape=(latent_dim,), name='z')([z_mean, z_log_var])

# instantiate encoder model
encoder = Model(inputs, [z_mean, z_log_var, z], name='encoder')
encoder.summary()
#plot_model(encoder, to_file='vae_mlp_encoder.png', show_shapes=True)

# build decoder model
latent_inputs = Input(shape=(latent_dim,), name='z_sampling')
x = Dense(intermediate_dim, activation='relu')(latent_inputs)
outputs = Dense(original_dim, activation='sigmoid', name='decoder')(x)

# instantiate decoder model
decoder = Model(latent_inputs, outputs, name='decoded')
decoder.summary()
#plot_model(decoder, to_file='vae_mlp_decoder.png', show_shapes=True)

# instantiate VAE model
outputs = decoder(encoder(inputs)[2])
vae = Model(inputs, outputs, name='vae_mlp')

def vae_loss(inputs, outputs):
    reconstruction_loss = binary_crossentropy(inputs,
                                              outputs)
    reconstruction_loss *= original_dim
    kl_loss = 1 + z_log_var - K.square(z_mean) - K.exp(z_log_var)
    kl_loss = K.sum(kl_loss, axis=-1)
    kl_loss *= -0.5
    vae_loss = K.mean(reconstruction_loss + kl_loss)
    return vae_loss

#vae.add_loss(vae_loss)
vae.compile(optimizer='adam', loss=vae_loss)
#tf.summary.image("input",tf.reshape(vae.get_layer('encoder_input').output,
 #                                   [1, size, size, 1]))
#tf.summary.image("image",tf.reshape(decoder.get_layer('decoder').output,[81,
# 81,1]))
vae.summary()
vae.fit(x_train, x_test,
                epochs=epochs,
                batch_size=batch_size,
                validation_split=0.2, callbacks=[TensorBoard(
        log_dir='/tmp/tb', histogram_freq=1,
                                            write_graph=False,
                                            write_images=False)])

prediction = vae.predict(x_train)

for im in prediction:
    plt.imshow(im.reshape(image_size,image_size))
    plt.gray()
    plt.show()
'''plot_model(vae,
           to_file='vae_mlp.png',
           show_shapes=True)'''

'''if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    help_ = "Load h5 model trained weights"
    parser.add_argument("-w", "--weights", help=help_)
    help_ = "Use mse loss instead of binary cross entropy (default)"
    parser.add_argument("-m",
                        "--mse",
                        help=help_, action='store_true')
    args = parser.parse_args()
    models = (encoder, decoder)
    data = (x_test, y_test)

    # VAE loss = mse_loss or xent_loss + kl_loss
    if args.mse:
        reconstruction_loss = mse(inputs, outputs)
    else:
        reconstruction_loss = binary_crossentropy(inputs,
                                                  outputs)

    reconstruction_loss *= original_dim
    kl_loss = 1 + z_log_var - K.square(z_mean) - K.exp(z_log_var)
    kl_loss = K.sum(kl_loss, axis=-1)
    kl_loss *= -0.5
    vae_loss = K.mean(reconstruction_loss + kl_loss)
    vae.add_loss(vae_loss)
    vae.compile(optimizer='adam')
    vae.summary()
    plot_model(vae,
               to_file='vae_mlp.png',
               show_shapes=True)

    if args.weights:
        vae = vae.load_weights(args.weights)
    else:
        # train the autoencoder
        vae.fit(x_train,
                epochs=epochs,
                batch_size=batch_size,
                validation_data=(x_test, None))
        vae.save_weights('vae_mlp_mnist.h5')

    plot_results(models,
                 data,
                 batch_size=batch_size,
model_name="vae_mlp")'''
