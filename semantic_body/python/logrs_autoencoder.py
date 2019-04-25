import numpy as np
from tensorflow.keras.layers import Dense, Input, Dropout
from tensorflow.keras.models import Model
RS = np.fromfile('../data/train/logRS')[2:]
RS.resize(1400, 300000)

D = np.fromfile('../data/train/roughExact')[2:]
D.resize(1400, 26)
from sklearn.model_selection import train_test_split
x_train, x_test, y_train, y_test = train_test_split(D, RS, test_size=0.1, random_state=2)



input_dim = 300000
latent_dim = 300

inputs = Input(shape=(input_dim,))
encoded = Dense(latent_dim, activation='linear')(inputs)
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
    return K.mean(K.square(y_pred-y_true))

# train
autoencoder.compile(optimizer='adam', loss=myloss)
history = autoencoder.fit(y_train, y_train, epochs=10, batch_size=60, shuffle=True, validation_data=(y_test, y_test))

print(((y_test- decoder.predict(encoder.predict(y_test))**2).mean()))

logRS=np.fromfile('../data/test/logRS')[2:]
logRS.resize(111, 300000)
out = decoder.predict(encoder.predict(logRS))
out = out.astype(np.float64)
out.tofile('../data/recover/newlogRS')

encoder.save('../export/logrs_encoder.h5')
decoder.save('../export/logrs_decoder.h5')

import matplotlib.pyplot as plt
history = autoencoder.history
plt.plot(history.history['loss'][0:], 'r')
plt.plot(history.history['val_loss'][0:], 'b')
plt.legend(['train', 'test'])
plt.title('logRS Autoencoder Training Process_Adam_MSE')
plt.savefig('../images/chartgraph/log_ae_training_process_adam_mse.png', bbox_inches='tight')