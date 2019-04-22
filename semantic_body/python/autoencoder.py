import numpy as np
from tensorflow.keras.layers import Dense, Input, Dropout
from tensorflow.keras.models import Model
RS = np.fromfile('../data/train/RS')[2:]
RS.resize(1400, 225000)
RS = np.log(RS)

D = np.fromfile('../data/train/roughExact')[2:]
D.resize(1400, 26)
from sklearn.model_selection import train_test_split
x_train, x_test, y_train, y_test = train_test_split(D, V, test_size=0.1, random_state=2)



input_dim = 225000
latent_dim = 60

inputs = Input(shape=(input_dim,))
encoded = Dense(latent_dim, activation='sigmoid')(inputs)
decoded = Dense(input_dim, activation='linear')(encoded)

autoencoder = Model(inputs=inputs, outputs=decoded)
# this model maps an input to its encoded representation
encoder = Model(inputs, encoded)
# create a placeholder for an encoded (32-dimensional) input
encoded_input = Input(shape=(latent_dim,))
# retrieve the last layer of the autoencoder model
decoder_layer = autoencoder.layers[-1]
# create the decoder model
decoder = Model(encoded_input, decoder_layer(encoded_input))

import tensorflow.keras.backend as K
def myloss(y_pred, y_true):
    return K.sum(K.square(y_pred-y_true))

# train
autoencoder.compile(optimizer='adam', loss=myloss)
history = autoencoder.fit(y_train, y_train, epochs=1000, batch_size=1260, shuffle=True, validation_data=(y_test, y_test))


encoder.save('../export/rs_encoder.h5')
decoder.save('../export/rs_decoder.h5')

import matplotlib.pyplot as plt
history = autoencoder.history
plt.plot(history.history['loss'][0:], 'r')
plt.plot(history.history['val_loss'][0:], 'b')
plt.legend(['train', 'test'])
plt.title('Wide Autoencoder Training Process_Adam_MSE')
plt.savefig('../images/chartgraph/wide_ae_training_process_adam_mse.png', bbox_inches='tight')