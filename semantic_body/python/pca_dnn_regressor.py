import numpy as np
# load data
V = np.fromfile('../data/train/dV')[2:]
D = np.fromfile('../data/train/roughExact')[2:]

V.resize(1400, 37500)
D.resize(1400, 26)

#PCA 
from sklearn.decomposition import PCA
pca = PCA(n_components=36)
newV = pca.fit_transform(V) #(1400, 36)

#Dnn Regressor
import os
import tensorflow as tf
import tensorflow.contrib.learn as skflow
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' 
## Data Loader
from sklearn.model_selection import train_test_split
x_train, x_test, y_train, y_test = train_test_split(D, newV, test_size=0.1, random_state=2)

def train_input_fn(batch_size=1):
    return tf.estimator.inputs.numpy_input_fn(
                    x={'X':x_train.astype(np.float32)}, 
                    y=y_train.astype(np.float32),
                    batch_size=batch_size,
                    shuffle=True)

def test_input_fn(batch_size=1):
    return tf.estimator.inputs.numpy_input_fn(
                    x={'X':x_test.astype(np.float32)}, 
                    y=y_test.astype(np.float32),
                    batch_size=batch_size,
                    shuffle=True)

feature_columns = [tf.feature_column.numeric_column('X', shape=(26,))]

## Dense Neural Network
hidden_layers = [32, 64, 128, 64]
dropout = 0.0

MODEL_PATH='../export/dnnregressor/'
for hl in hidden_layers:
    MODEL_PATH += '%s_' % hl
MODEL_PATH += 'D0%s' % (int(dropout*10))

validation_metrics = {"MSE": tf.contrib.metrics.streaming_mean_squared_error}


strategy = tf.contrib.distribute.OneDeviceStrategy(device='/gpu:0')
config = tf.estimator.RunConfig(train_distribute=strategy)

regressor = skflow.DNNRegressor(feature_columns=feature_columns,
                label_dimension=36,
                config=config,
                hidden_units=hidden_layers,
                model_dir=MODEL_PATH,
                dropout=dropout)

## Train

for epoch in range(500):
    regressor.fit(input_fn=train_input_fn(batch_size=1), steps=100)
    
    if(epoch%10==0):
        ev = regressor.evaluate(input_fn=test_input_fn(), metrics=validation_metrics)
        
        print('Epoch %i: %.5f MSE' % (epoch+1, ev['MSE']))

# evaluate
from sklearn import metrics
pred = np.dot(regressor.predict(x={'X', D}), pca.components_)
print("MSE:")
print(metrics.mean_squared_error(V, pred))

# save model
feature_spec = tf.feature_column.make_parse_example_spec(feature_columns)
serving_input_receiver_fn = tf.contrib.learn.build_parsing_serving_input_fn(feature_spec)
export_dir = regressor.export_savedmodel('../export/', serving_input_receiver_fn)
print('model saved!!!')