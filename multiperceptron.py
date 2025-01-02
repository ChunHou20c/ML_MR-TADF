from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
from sklearn.model_selection import GridSearchCV

import matplotlib.pyplot as plt
import math
import seaborn as sns
import pandas as pd
from sklearn.ensemble import RandomForestRegressor

import tensorflow as tf 
from tensorflow.keras.models import Sequential 
from tensorflow.keras.layers import Flatten 
from tensorflow.keras.layers import Dense 
from tensorflow.keras.layers import Activation

class model_mlp:

  def lr_exp_decay(epoch, lr):
      return 0.012 * math.exp(-0.0018*epoch)

  checkpoint_filepath = 'checkpoint/cp.mlp.weights.h5'
  model_checkpoint_callback = tf.keras.callbacks.ModelCheckpoint(
      filepath=checkpoint_filepath,
      save_weights_only=True,
      save_best_only=False)

  callbacks = []
  callbacks.append(tf.keras.callbacks.LearningRateScheduler(lr_exp_decay))
  callbacks.append(model_checkpoint_callback)

  opt = tf.keras.optimizers.Adam(learning_rate=0.0012)

  def __init__(self, shape) -> None:
    
    input_matrix = tf.keras.Input(shape=shape, name = 'inter')
    x = tf.keras.layers.Normalization(axis=None)(input_matrix)
    x = tf.keras.layers.Flatten()(x)
    x = tf.keras.layers.Dense(128, activation = 'swish')(x)
    x = tf.keras.layers.BatchNormalization()(x)
    x = tf.keras.layers.Dense(64, activation= 'swish')(x)
    x = tf.keras.layers.BatchNormalization()(x)
    x = tf.keras.layers.Dense(32, activation= 'swish')(x)
    x = tf.keras.layers.BatchNormalization()(x)
    x = tf.keras.layers.Dense(16, activation= 'swish')(x)
    x = tf.keras.layers.Dense(1, activation='linear')(x)
    x = tf.keras.layers.Lambda(abs)(x)
    self.model = tf.keras.Model(input_matrix, x, name = 'x')

  def compile(self):
    self.model.compile(optimizer=self.opt, loss='mse', metrics=['mse'])

  def summary(self):
    self.model.summary()

  def train(self, training_feature, training_label, testing_feature, testing_label, epoches):
    history = self.model.fit(training_feature,
        training_label,
        batch_size = 25,
        epochs = epoches,
        validation_data = (testing_feature, testing_label),
        callbacks=self.callbacks)
    return history

  def evaluate_model(self, training_feature, training_label, testing_feature, testing_label):
    
    self.model.load_weights(self.checkpoint_filepath)

    predicted_training_coupling = self.model.predict(training_feature)
    predicted_testing_coupling = self.model.predict(testing_feature)

    mse_training, r2_training = self.model.evaluate(training_feature, training_label)
    mse_testing, r2_testing = self.model.evaluate(testing_feature, testing_label)

    training_dataset = pd.DataFrame(list(zip(predicted_training_coupling.flatten(), training_label.flatten())), columns = ['predicted value', 'actual value'])
    testing_dataset = pd.DataFrame(list(zip(predicted_testing_coupling.flatten(), testing_label.flatten())), columns = ['predicted value', 'actual value'])

    print(training_dataset)
    print(testing_dataset)
    sns.set_style('dark')
    sns.set_theme(rc={'figure.figsize':(6,6)})
    #plot = sns.scatterplot(data=training_dataset, x = 'predicted value', y = 'actual value')
    
    fig, ax = plt.subplots(figsize=(6, 6))

    fig.suptitle('Training Data', fontsize=15)
    ax.set_title(f"$R^2$ = {r2_training:.3f} \nmse = {mse_training:.4f}", fontsize=11, loc='left')
    sns.scatterplot(
        data=training_dataset,
        x="predicted value",
        y="actual value",
        color="blue",
        ax=ax,
    )
    sns.kdeplot(
        data=training_dataset,
        x="predicted value",
        y="actual value",
        levels=5,
        fill=True,
        alpha=0.6,
        cut=2,
        ax=ax,
    )
    ax.plot([0,0.75], [0,0.75],color = 'b')
    ax.set_xlabel('predicted electronic coupling (eV)')
    ax.set_ylabel('actual electronic coupling (eV)')
    
    fig.savefig('training_set_mlp.png')
    plt.show(fig)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 6))

    fig.suptitle('Testing Data', fontsize=15)
    ax.set_title(f"$R^2$ = {r2_testing:.3f} \nmse = {mse_testing:.4f}", fontsize=11, loc='left')
    sns.scatterplot(
        data=testing_dataset,
        x="predicted value",
        y="actual value",
        color="blue",
        ax=ax,
    )
    sns.kdeplot(
        data=testing_dataset,
        x="predicted value",
        y="actual value",
        levels=5,
        fill=True,
        alpha=0.6,
        cut=2,
        ax=ax,
    )
    ax.plot([0,0.75], [0,0.75],color = 'b')
    ax.set_xlabel('predicted electronic coupling (eV)')
    ax.set_ylabel('actual electronic coupling (eV)')
    fig.savefig('testing_set_mlp.png')
    plt.show(fig)
    plt.close(fig)

if __name__ == "__main__":
    #load data

    input = pd.read_csv('./target/training_data/input.csv')
    output = pd.read_csv('./target/training_data/output.csv')

    X_train, X_test, y_train, y_test = train_test_split(input, output, test_size = 0.2, random_state = 0)

    shape = (X_train.to_numpy().shape[1],)
    model = model_mlp(shape)
    model.summary()
    model.compile()

    history = model.train(X_train.to_numpy(), y_train.to_numpy(), X_test.to_numpy(), y_test.to_numpy(), 100)

    model.evaluate_model(X_train.to_numpy(), y_train.to_numpy(), X_test.to_numpy(), y_test.to_numpy())


