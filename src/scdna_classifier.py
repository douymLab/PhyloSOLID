#!/usr/bin/env python3

###################################################################################################
######################### 构建分类器并进行实时预测（放松限制版本） #######################
###################################################################################################

import os
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
import joblib
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns

# 设置路径
# features_file_labeled = "phylosolid/models/data_labeling_for_classifier_and_ROC.txt"
features_file_labeled = "classifier/scdna/data_labeling_for_classifier_and_ROC.txt"

# 定义特征列
SELECTED_FEATURES = [
    "VAF_all", 
    "mutant_cell_frac", 
    "unmut_sc_AF_mean", 
    "sc_AF_mean", 
    "max_sc_mutant_read_count_normalized", 
    "sc_mismatches_mean_max", 
    "pseudo_bulk_mismatches_alt_mean", 
    "mutant_popAF", 
    "pseudo_bulk_indels_ratio", 
    "sc_softclippedreads_ratio_max",
]

def load_training_data():
    """
    加载训练数据
    """
    print("=== 加载训练数据 ===")
    df_all_features = pd.read_csv(features_file_labeled, sep="\t", index_col=0)
    print(f"数据形状: {df_all_features.shape}")
    print(f"样本分布: {Counter(df_all_features['sampleid'])}")
    print(f"类别分布: {Counter(df_all_features['label'])}")
    
    return df_all_features

def build_relaxed_classifier_excluding_sample(df_training, exclude_sample_id):
    """
    构建放松限制的分类器（排除特定样本）
    
    Args:
        df_training: 训练数据
        exclude_sample_id: 要排除的样本ID
        
    Returns:
        tuple: (model, pipeline)
    """
    print(f"\n=== 构建放松限制的分类器 (排除样本 {exclude_sample_id}) ===")
    
    # 排除目标样本
    df_train = df_training[df_training['sampleid'] != exclude_sample_id].copy()
    print(f"训练数据大小: {df_train.shape}")
    print(f"训练数据类别分布: {Counter(df_train['label'])}")
    
    # 准备特征和标签
    X_train = df_train[SELECTED_FEATURES].copy()
    y_train = df_train['label'].copy()
    
    # 创建预处理管道
    pipeline = make_pipeline(
        SimpleImputer(strategy='median'),
        StandardScaler()
    )
    
    # 数据预处理和训练
    print("数据预处理和模型训练...")
    X_train_processed = pipeline.fit_transform(X_train)
    
    # 使用放松的参数设置
    rf_model = RandomForestClassifier(
        n_estimators=1000,  # 减少树的数量，降低过拟合
        random_state=42,
        class_weight={  # 手动设置类别权重，偏向mosaic
            'mosaic': 1.0,     # 给mosaic更高的权重
            'germline_het': 1.0,   # 降低germline的权重
            'repeat': 1.0      # 保持repeat权重不变
        },
        max_depth=15,           # 增加树深度，捕捉更多模式
        min_samples_split=2,    # 减少最小分割样本数
        min_samples_leaf=1,     # 减少叶子节点最小样本数
        max_features='sqrt',    # 使用更少的特征进行分割
        bootstrap=True          # 使用bootstrap采样
    )
    rf_model.fit(X_train_processed, y_train)
    
    print(f"模型训练完成，类别: {rf_model.classes_}")
    
    return rf_model, pipeline

def predict_with_relaxed_threshold(model, pipeline, df_new, sample_id, output_file=None):
    """
    使用放松的阈值进行预测
    
    Args:
        model: 训练好的模型
        pipeline: 预处理管道
        df_new: 新数据框
        sample_id: 样本ID
        output_file: 输出文件路径
        
    Returns:
        pd.DataFrame: 预测结果
    """
    print(f"\n=== 为样本 {sample_id} 预测新突变 (放松阈值) ===")
    print(f"新数据形状: {df_new.shape}")
    
    # 检查必需的特征列
    missing_features = [feat for feat in SELECTED_FEATURES if feat not in df_new.columns]
    if missing_features:
        raise ValueError(f"缺少必需的特征列: {missing_features}")
    
    # 检查mutation_id列
    if 'mutation_id' not in df_new.columns:
        raise ValueError("输入数据必须包含 'mutation_id' 列")
    
    # 准备特征数据
    X_new = df_new[SELECTED_FEATURES].copy()
    
    # 使用相同的预处理管道
    X_new_processed = pipeline.transform(X_new)
    
    # 进行预测 - 使用概率而不是硬分类
    probabilities = model.predict_proba(X_new_processed)
    class_labels = model.classes_
    
    # 放松的预测逻辑：如果mosaic概率 > 0.5 就预测为mosaic
    predictions = []
    for i, prob_vector in enumerate(probabilities):
        mosaic_prob = prob_vector[list(class_labels).index('mosaic')] if 'mosaic' in class_labels else 0
        germline_prob = prob_vector[list(class_labels).index('germline')] if 'germline' in class_labels else 0
        repeat_prob = prob_vector[list(class_labels).index('repeat')] if 'repeat' in class_labels else 0
        
        # 放松的决策规则
        if mosaic_prob > 0.5:  # 降低mosaic的阈值
            predictions.append('mosaic')
        elif germline_prob > 0.6:  # 提高germline的阈值
            predictions.append('germline')
        elif repeat_prob > 0.6:   # 提高repeat的阈值
            predictions.append('repeat')
        else:
            # 如果都不满足，再放宽 mosaic 的条件
            if mosaic_prob > 0.2:  # 降低mosaic的阈值
                predictions.append('mosaic')
            else:
                # 如果 mosaic 最低都不满足就选择概率最高的
                predictions.append(class_labels[np.argmax(prob_vector)])
    
    # 创建结果DataFrame
    results_df = pd.DataFrame({
        'mutation_id': df_new['mutation_id'].values,
        'sample_id': sample_id,
        'predicted_label': predictions
    })
    
    # 添加概率分数
    for i, class_name in enumerate(class_labels):
        results_df[f'probability_{class_name}'] = probabilities[:, i]
    
    # 添加决策信息
    results_df['decision_rule'] = 'relaxed_threshold'
    
    # 打印预测统计
    prediction_counts = pd.Series(predictions).value_counts().to_dict()
    print(f"\n放松阈值预测统计:")
    for label, count in prediction_counts.items():
        percentage = count / len(predictions) * 100
        print(f"  {label}: {count} 个位点 ({percentage:.1f}%)")
    
    # 保存结果
    if output_file:
        results_df.to_csv(output_file, index=False)
        print(f"预测结果保存到: {output_file}")
    
    return results_df

def analyze_feature_importance(model, output_path):
    """
    分析特征重要性
    """
    print(f"\n=== 特征重要性分析 ===")
    
    importances = model.feature_importances_
    feature_imp_df = pd.DataFrame({
        'feature': SELECTED_FEATURES,
        'importance': importances
    }).sort_values('importance', ascending=False)
    
    print("特征重要性排序:")
    for _, row in feature_imp_df.iterrows():
        print(f"  {row['feature']}: {row['importance']:.4f}")
    
    return feature_imp_df

def generate_relaxed_prediction_report(results_df, sample_id, output_path):
    """
    生成放松预测的详细报告
    """
    report_file = os.path.join(output_path, f"relaxed_prediction_report_{sample_id}.txt")
    
    with open(report_file, 'w') as f:
        f.write(f"放松限制分类器预测报告 - 样本 {sample_id}\n")
        f.write("=" * 50 + "\n")
        f.write(f"总位点数: {len(results_df)}\n\n")
        
        # 预测统计
        pred_counts = results_df['predicted_label'].value_counts()
        f.write("预测结果统计:\n")
        for label, count in pred_counts.items():
            percentage = count / len(results_df) * 100
            f.write(f"  {label}: {count} ({percentage:.1f}%)\n")
        
        f.write("\n放松策略说明:\n")
        f.write("1. mosaic概率 > 0.5 即预测为mosaic\n")
        f.write("2. germline概率 > 0.5 才预测为germline\n") 
        f.write("3. repeat概率 > 0.5 才预测为repeat\n")
        f.write("4. 偏向于预测为mosaic，减少假阴性\n")
        
        # 分析mosaic概率分布
        if 'probability_mosaic' in results_df.columns:
            mosaic_probs = results_df['probability_mosaic']
            f.write(f"\nMosaic概率分布:\n")
            f.write(f"  最小值: {mosaic_probs.min():.3f}\n")
            f.write(f"  最大值: {mosaic_probs.max():.3f}\n")
            f.write(f"  平均值: {mosaic_probs.mean():.3f}\n")
            f.write(f"  中位数: {mosaic_probs.median():.3f}\n")
            
            # 不同概率区间的统计
            bins = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
            for i in range(len(bins)-1):
                count = ((mosaic_probs >= bins[i]) & (mosaic_probs < bins[i+1])).sum()
                f.write(f"  [{bins[i]:.1f}-{bins[i+1]:.1f}): {count} 个位点\n")
    
    print(f"放松预测报告生成: {report_file}")

def real_time_classifier_predict(df_for_classifier, sample_id, output_path):
    """
    放松限制的实时分类器预测
    
    Args:
        df_for_classifier: 需要预测的数据框
        sample_id: 样本ID
        output_path: 输出路径
        
    Returns:
        pd.DataFrame: 预测结果
    """
    print(f"\n{'='*60}")
    print(f"开始放松限制的实时分类器预测 - 样本: {sample_id}")
    print(f"{'='*60}")
    
    # 确保输出目录存在
    os.makedirs(output_path, exist_ok=True)
    
    # 1. 加载训练数据
    df_training = load_training_data()
    
    # 2. 检查样本是否在训练数据中
    training_samples = set(df_training['sampleid'].unique())
    if sample_id not in training_samples:
        print(f"警告: 样本 {sample_id} 不在训练数据中")
        exclude_sample_id = sample_id
    else:
        exclude_sample_id = sample_id
        print(f"样本 {sample_id} 在训练数据中，使用leave-one-out策略")
    
    # 3. 构建放松限制的分类器
    model, pipeline = build_relaxed_classifier_excluding_sample(df_training, exclude_sample_id)
    
    # 4. 保存模型到pkl文件
    model_file = os.path.join(output_path, f"relaxed_classifier_{sample_id}.pkl")
    joblib.dump({
        'model': model,
        'pipeline': pipeline,
        'sample_id': sample_id,
        'excluded_sample': exclude_sample_id,
        'feature_names': SELECTED_FEATURES,
        'training_date': pd.Timestamp.now()
    }, model_file)
    print(f"模型保存到: {model_file}")
    
    # 5. 使用放松的阈值进行预测
    output_file = os.path.join(output_path, f"relaxed_predictions_{sample_id}.csv")
    results = predict_with_relaxed_threshold(model, pipeline, df_for_classifier, sample_id, output_file)
    
    # 6. 分析特征重要性
    feature_imp_df = analyze_feature_importance(model, output_path)
    
    # 7. 保存特征重要性
    feature_imp_file = os.path.join(output_path, f"feature_importance_{sample_id}.csv")
    feature_imp_df.to_csv(feature_imp_file, index=False)
    print(f"特征重要性保存到: {feature_imp_file}")
    
    # 8. 生成详细报告
    generate_relaxed_prediction_report(results, sample_id, output_path)
    
    return results, model, pipeline

# 添加模型加载函数
def load_relaxed_classifier(model_path):
    """
    加载保存的放松分类器
    
    Args:
        model_path: 模型文件路径
        
    Returns:
        dict: 包含模型、管道等信息的字典
    """
    print(f"加载模型: {model_path}")
    classifier_data = joblib.load(model_path)
    
    model = classifier_data['model']
    pipeline = classifier_data['pipeline']
    sample_id = classifier_data['sample_id']
    
    print(f"加载的模型信息:")
    print(f"  样本ID: {sample_id}")
    print(f"  排除的样本: {classifier_data['excluded_sample']}")
    print(f"  特征数量: {len(classifier_data['feature_names'])}")
    print(f"  训练日期: {classifier_data['training_date']}")
    print(f"  模型类别: {model.classes_}")
    
    return classifier_data

# 添加使用已保存模型的预测函数
def predict_with_saved_classifier(df_for_classifier, model_path, sample_id, output_file=None):
    """
    使用已保存的模型进行预测
    
    Args:
        df_for_classifier: 需要预测的数据框
        model_path: 模型文件路径
        sample_id: 样本ID
        output_file: 输出文件路径
        
    Returns:
        pd.DataFrame: 预测结果
    """
    # 加载模型
    classifier_data = load_relaxed_classifier(model_path)
    model = classifier_data['model']
    pipeline = classifier_data['pipeline']
    
    # 使用放松的阈值进行预测
    results = predict_with_relaxed_threshold(model, pipeline, df_for_classifier, sample_id, output_file)
    
    return results


# # 修改使用示例
# def main():
#     """
#     主函数 - 使用示例
#     """
#     # 您的数据框
#     df_for_classifier = pd.DataFrame({
#         'mutation_id': [
#             'X_135442114_T_G', 'X_138790531_G_A', 'X_70632344_G_A',
#             '1_106340339_G_T', '1_212476040_A_G', '19_38404452_C_T',
#             '20_13387559_C_T', '20_165000_A_G', '20_44421717_C_T', '20_57865731_C_T'
#         ],
#         'VAF_all': [0.1, 0.2, 0.05, 0.15, 0.08, 0.12, 0.03, 0.18, 0.09, 0.11],
#         'mutant_cell_frac': [0.3, 0.4, 0.1, 0.35, 0.25, 0.28, 0.08, 0.42, 0.22, 0.31],
#         'unmut_sc_AF_mean': [0.01, 0.02, 0.005, 0.015, 0.008, 0.012, 0.003, 0.018, 0.009, 0.011],
#         'sc_AF_mean': [0.15, 0.25, 0.08, 0.18, 0.12, 0.16, 0.05, 0.28, 0.14, 0.17],
#         'max_sc_mutant_read_count_normalized': [5, 8, 2, 6, 4, 5, 1, 9, 3, 6],
#         'sc_mismatches_mean_max': [0.1, 0.15, 0.05, 0.12, 0.08, 0.11, 0.04, 0.16, 0.09, 0.13],
#         'pseudo_bulk_mismatches_alt_mean': [0.02, 0.03, 0.01, 0.025, 0.015, 0.022, 0.008, 0.035, 0.018, 0.024],
#         'mutant_popAF': [0.12, 0.22, 0.06, 0.16, 0.1, 0.14, 0.04, 0.25, 0.13, 0.15],
#         'pseudo_bulk_indels_ratio': [0.01, 0.02, 0.005, 0.015, 0.008, 0.012, 0.004, 0.022, 0.01, 0.014],
#         'sc_softclippedreads_ratio_max': [0.05, 0.08, 0.02, 0.06, 0.04, 0.055, 0.015, 0.09, 0.045, 0.065]
#     })
    
#     sample_id = "UMB1465"
#     output_path = "./relaxed_predictions"
    
#     # 方法1: 实时构建分类器并进行预测
#     print("=== 方法1: 实时构建分类器 ===")
#     results, model, pipeline = real_time_classifier_predict(df_for_classifier, sample_id, output_path)
    
#     print(f"\n最终预测结果预览:")
#     print(results[['mutation_id', 'predicted_label', 'probability_mosaic']].head())
    
#     # 方法2: 使用已保存的模型进行预测
#     print("\n=== 方法2: 使用已保存的模型 ===")
#     model_path = os.path.join(output_path, f"relaxed_classifier_{sample_id}.pkl")
#     if os.path.exists(model_path):
#         results_saved = predict_with_saved_classifier(
#             df_for_classifier, 
#             model_path, 
#             sample_id,
#             output_file=os.path.join(output_path, f"predictions_with_saved_model_{sample_id}.csv")
#         )
#         print(f"使用保存模型的预测结果预览:")
#         print(results_saved[['mutation_id', 'predicted_label', 'probability_mosaic']].head())
#     else:
#         print(f"模型文件不存在: {model_path}")

# if __name__ == "__main__":
#     main()

