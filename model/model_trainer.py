import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report
import joblib
import matplotlib.pyplot as plt
import seaborn as sns

def train_random_forest(data_path, output_dir='./model/'):
    """
    Train a Random Forest model on primer design data and save the model and evaluation metrics.
    
    :param data_path: Path to the training data CSV file
    :param output_dir: Directory to save model and evaluation files
    :return: Trained model and accuracy on test set
    """
    # Load data
    data = pd.read_csv(data_path)
    
    # Separate features and target
    X = data.drop('quality', axis=1)  # Replace 'quality' with your actual target column name
    y = data['quality']  # Replace with your actual target column name
    
    # Split data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    # Create and train the model
    model = RandomForestClassifier(n_estimators=100, random_state=42)
    model.fit(X_train, y_train)
    
    # Evaluate model
    y_pred = model.predict(X_test)
    accuracy = accuracy_score(y_test, y_pred)
    
    print(f"Model Accuracy: {accuracy:.4f}")
    print("\nClassification Report:")
    print(classification_report(y_test, y_pred))
    
    # Save model
    model_path = f"{output_dir}random_forest_model.pkl"
    joblib.dump(model, model_path)
    print(f"Model saved to {model_path}")
    
    # Feature importance
    feature_importance = pd.DataFrame({
        'feature': X.columns,
        'importance': model.feature_importances_
    }).sort_values('importance', ascending=False)
    
    # Save feature importance to CSV
    feature_importance.to_csv(f"{output_dir}feature_importances.csv", index=False)
    
    # Plot feature importance
    plt.figure(figsize=(12, 8))
    sns.barplot(x='importance', y='feature', data=feature_importance[:20])  # Top 20 features
    plt.title('Feature Importance')
    plt.tight_layout()
    plt.savefig(f"{output_dir}feature_importances.png")
    plt.savefig(f"{output_dir}feature_importances.pdf")
    
    # Plot confusion matrix
    cm = confusion_matrix(y_test, y_pred)
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues')
    plt.title('Confusion Matrix')
    plt.ylabel('True Label')
    plt.xlabel('Predicted Label')
    plt.tight_layout()
    plt.savefig(f"{output_dir}confusion_matrix.png")
    plt.savefig(f"{output_dir}confusion_matrix.pdf")
    
    return model, accuracy

if __name__ == "__main__":
    # Example usage 
    data_path = "./data/rfmodel_data/training_data.csv"  # Replace with actual path
    model, accuracy = train_random_forest(data_path)
    print(f"Final model accuracy: {accuracy:.4f}")
