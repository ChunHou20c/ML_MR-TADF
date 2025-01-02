from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
from sklearn.model_selection import GridSearchCV

import matplotlib.pyplot as plt
import math
import seaborn as sns
import pandas as pd
from sklearn.ensemble import RandomForestRegressor

class model_rf_regressor:

  def __init__(self) -> None:

    rf_regressor = RandomForestRegressor(max_depth=250,n_estimators = 400, oob_score=True, random_state = 20)
    self.model = rf_regressor

  def summary(self):
    self.model.get_params()

  def train(self, training_feature, training_label):
        self.model.fit(training_feature, training_label)

  def evaluate_model(self, training_feature, training_label, testing_feature, testing_label):
    
    predicted_training_value = self.model.predict(training_feature)
    predicted_testing_value = self.model.predict(testing_feature)

    r2_training = r2_score(training_label, predicted_training_value)
    r2_testing = r2_score(testing_label, predicted_testing_value)

    print(predicted_training_value)
    print(predicted_testing_value)
    print(r2_training)
    print(r2_testing)
    training_dataset = pd.DataFrame(list(zip(predicted_training_value, training_label.flatten())), columns = ['predicted value', 'actual value'])
    testing_dataset = pd.DataFrame(list(zip(predicted_testing_value, testing_label.flatten())), columns = ['predicted value', 'actual value'])

    print(training_dataset)
    sns.set_style('dark')
    sns.set(rc={'figure.figsize':(6,6)})
    fig, ax = plt.subplots(figsize=(6, 6))

    fig.suptitle('Training Data', fontsize=15)
    ax.set_title(f"$R^2$ = {r2_training:.3f}", fontsize=11, loc='left')
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
    ax.set_xlabel('predicted value')
    ax.set_ylabel('actual value')

    fig.savefig('training_set_rf.png')
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 6))

    fig.suptitle('Testing Data', fontsize=15)
    ax.set_title(f"$R^2$ = {r2_testing}", fontsize=11, loc='left')
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
    ax.set_xlabel('predicted value')
    ax.set_ylabel('actual value')
    fig.savefig('testing_set_rf.png')
    plt.close(fig)

if __name__ == "__main__":
    #load data

    input = pd.read_csv('./target/training_data/input.csv')
    output = pd.read_csv('./target/training_data/output.csv')

    X_train, X_test, y_train, y_test = train_test_split(input, output, test_size = 0.2, random_state = 0)

    model = model_rf_regressor()
    model.summary()

    model.train(X_train.to_numpy(), y_train.to_numpy())

    model.evaluate_model(X_train.to_numpy(), y_train.to_numpy(), X_test.to_numpy(), y_test.to_numpy())

