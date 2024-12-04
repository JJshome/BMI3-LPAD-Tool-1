import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
import joblib

# 1. 加载数据
primers = pd.read_csv('data/primers_scoring_output.csv')
nonsense_primers = pd.read_csv('data/nonsense_primers_scoring_output.csv')

# 2. 数据预处理：合并数据并标记类别
primers['label'] = 1  # 正常引物标记为1
nonsense_primers['label'] = 0  # 无效引物标记为0

# 3. 合并数据集
data = pd.concat([primers, nonsense_primers])

# 4. 提取特征（从第二列到最后一列）和标签（第一列是lamp_id，最后一列是label）
X = data.iloc[:, 1:-1]  # 特征
y = data['label']  # 标签

# 5. 划分训练集和测试集
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# 6. 创建并训练随机森林模型
rf_model = RandomForestClassifier(n_estimators=100, random_state=42)
rf_model.fit(X_train, y_train)

# 7. 在测试集上进行预测
y_pred = rf_model.predict(X_test)

# 8. 输出模型评估结果
print("Classification Report:\n", classification_report(y_test, y_pred))

# 9. 保存训练好的模型
joblib.dump(rf_model, 'model/random_forest_model.pkl')

# 10. 输出特征的重要性（权重）
feature_importances = pd.DataFrame(rf_model.feature_importances_, index=X.columns, columns=["importance"])
print("\nFeature Importances:\n", feature_importances)

# 11. 保存特征重要性到CSV文件
feature_importances.to_csv('model/feature_importances.csv')
