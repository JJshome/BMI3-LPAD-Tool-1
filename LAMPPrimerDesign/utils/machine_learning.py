from sklearn.linear_model import LinearRegression
import numpy as np

def train_model(X, y):
    """训练简单的线性回归模型"""
    model = LinearRegression()
    model.fit(X, y)
    return model

def predict_model(model, X):
    """使用训练好的模型进行预测"""
    return model.predict(X)
