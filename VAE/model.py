
from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.layers import Dense, Flatten
from tensorflow.keras import layers as kl
from tensorflow.keras import regularizers as kr
from tensorflow.keras import models as km
from tensorflow.keras import backend as K
from sklearn import preprocessing

class Sampling(kl.Layer):
    """Uses (z_mean, z_log_var) to sample z."""

    def call(self, inputs):
        z_mean, z_log_var = inputs
        batch = tf.shape(z_mean)[0]
        dim = tf.shape(z_mean)[1]
        epsilon = K.random_normal(shape=(batch, dim))
        return z_mean + tf.exp(0.5 * z_log_var) * epsilon


class Encoder(kl.Layer):
    
    def __init__(self, latent_dim=2):
        super(Encoder, self).__init__()
        self.sampling = Sampling()
        self.dense1 = kl.Dense(128, activation='relu')
        self.dense2 = kl.Dense(16, activation='relu')
        self.dense_mean = kl.Dense(latent_dim)
        self.dense_log_var = kl.Dense(latent_dim)
    
    def call(self, inputs):
        x = self.dense1(inputs)
        x = self.dense2(x)
        z_mean = self.dense_mean(x)
        z_log_var = self.dense_log_var(x)
        z = self.sampling([z_mean, z_log_var])
        return z_mean, z_log_var, z
    
    
class Decoder(kl.Layer):
    
    def __init__(self, input_dim):
        super(Decoder, self).__init__()
        self.dense1 = kl.Dense(16, activation='relu')
        self.dense2 = kl.Dense(128, activation='relu')
        self.dense_output = kl.Dense(input_dim, activation='relu')
    
    def call(self, inputs):
        x = self.dense1(inputs)
        x = self.dense2(x)
        outputs = self.dense_output(x)
        return outputs


# beta-VAE model
class VAE(keras.Model):
    def __init__(self, orig_dim, beta=1):
        super(VAE, self).__init__()
        self.orig_dim = orig_dim
        self.beta = beta
        self.encoder = Encoder()
        self.decoder = Decoder(orig_dim)

    def call(self, inputs):
        z_mean, z_log_var, z = self.encoder(inputs)
        reconstructed = self.decoder(z)
        # KL divergence regularization loss.
        kl_loss = - self.beta * 0.5 * tf.reduce_mean(1 + z_log_var - tf.square(z_mean) - tf.exp(z_log_var))
        self.add_loss(kl_loss)
        return reconstructed
    
    
def train_step(epochs, vae, X, batch_size=32):
    
    train_dataset = tf.data.Dataset.from_tensor_slices(X)
    train_dataset = train_dataset.shuffle(buffer_size=2048).batch(batch_size)

    optimizer = tf.keras.optimizers.Adam(learning_rate=1e-4)
    loss_metric = tf.keras.metrics.Mean()

    for epoch in range(epochs):
        if epoch % 10 == 0:
            print("Start of epoch " + str(epoch))

        # Iterate over the batches of the dataset.
        for step, x_batch_train in enumerate(train_dataset):
            with tf.GradientTape() as tape:
                reconstruction = vae(x_batch_train)
                # Compute reconstruction loss
                loss = tf.reduce_mean(keras.losses.mse(x_batch_train, reconstruction)) * X.shape[1]
                loss += sum(vae.losses)  # Add KLD regularization loss

            grads = tape.gradient(loss, vae.trainable_weights)
            optimizer.apply_gradients(zip(grads, vae.trainable_weights))

            loss_metric(loss)

            if epoch % 10 == 0 and step == 0:
                print("step %d: mean loss = %.4f, KL loss = %.4f" % (step, loss_metric.result(), sum(vae.losses)))