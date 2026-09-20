#!/usr/bin/env python3

###################################################################################################
######################### Relaxed-threshold classifier for real-time prediction ###################
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
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from src.reproducibility import set_seed, deterministic_permutation


# Training labels live next to this module, independent of the process cwd.
# src/classifier/scdna/data_labeling_for_classifier_and_ROC.txt
_SCRIPT_DIR = Path(__file__).resolve().parent
features_file_labeled = _SCRIPT_DIR / "classifier" / "scdna" / "data_labeling_for_classifier_and_ROC.txt"

# Feature columns used by the classifier
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
    Load labeled training data for the scDNA classifier.
    """
    print("=== Loading training data ===")
    if not features_file_labeled.exists():
        raise FileNotFoundError(
            f"scDNA classifier training file not found: {features_file_labeled}"
        )
    print(f"Training data path: {features_file_labeled}")
    df_all_features = pd.read_csv(features_file_labeled, sep="\t", index_col=0)
    print(f"Data shape: {df_all_features.shape}")
    print(f"Sample distribution: {Counter(df_all_features['sampleid'])}")
    print(f"Class distribution: {Counter(df_all_features['label'])}")
    
    return df_all_features

def build_relaxed_classifier_excluding_sample(df_training, exclude_sample_id):
    """
    Train a relaxed-threshold classifier while excluding a target sample.
    
    Args:
        df_training: Training dataframe
        exclude_sample_id: Sample ID to leave out of training
        
    Returns:
        tuple: (model, pipeline)
    """
    print(f"\n=== Building relaxed-threshold classifier (excluding sample {exclude_sample_id}) ===")
    
    # Exclude the target sample (leave-one-out)
    df_train = df_training[df_training['sampleid'] != exclude_sample_id].copy()
    print(f"Training data size: {df_train.shape}")
    print(f"Training class distribution: {Counter(df_train['label'])}")
    
    # Prepare features and labels
    X_train = df_train[SELECTED_FEATURES].copy()
    y_train = df_train['label'].copy()
    
    # Build the preprocessing pipeline
    pipeline = make_pipeline(
        SimpleImputer(strategy='median'),
        StandardScaler()
    )
    
    # Preprocess data and train the model
    print("Preprocessing data and training the model...")
    X_train_processed = pipeline.fit_transform(X_train)
    
    # Relaxed random-forest settings
    rf_model = RandomForestClassifier(
        n_estimators=1000,  # fewer trees to reduce overfitting
        random_state=42,
        class_weight={  # class weights
            'mosaic': 1.0,
            'germline_het': 1.0,
            'repeat': 1.0
        },
        max_depth=15,           # deeper trees to capture more patterns
        min_samples_split=2,    # smaller split size
        min_samples_leaf=1,     # smaller leaf size
        max_features='sqrt',    # use a subset of features at each split
        bootstrap=True          # bootstrap sampling
    )
    rf_model.fit(X_train_processed, y_train)
    
    print(f"Model training completed. Classes: {rf_model.classes_}")
    
    return rf_model, pipeline

def predict_with_relaxed_threshold(model, pipeline, df_new, sample_id, output_file=None):
    """
    Predict labels using relaxed probability thresholds.
    
    Args:
        model: Trained classifier
        pipeline: Fitted preprocessing pipeline
        df_new: New feature dataframe
        sample_id: Sample ID
        output_file: Optional path to write predictions
        
    Returns:
        pd.DataFrame: Prediction results
    """
    print(f"\n=== Predicting mutations for sample {sample_id} (relaxed thresholds) ===")
    print(f"New data shape: {df_new.shape}")
    
    # Check required feature columns
    missing_features = [feat for feat in SELECTED_FEATURES if feat not in df_new.columns]
    if missing_features:
        raise ValueError(f"Missing required feature columns: {missing_features}")
    
    # Check mutation_id column
    if 'mutation_id' not in df_new.columns:
        raise ValueError("Input data must contain a 'mutation_id' column")
    
    # Prepare feature matrix
    X_new = df_new[SELECTED_FEATURES].copy()
    
    # Apply the same preprocessing pipeline
    X_new_processed = pipeline.transform(X_new)
    
    # Predict using class probabilities rather than hard labels
    probabilities = model.predict_proba(X_new_processed)
    class_labels = model.classes_
    
    # Relaxed prediction rule: call mosaic if mosaic probability > 0.5
    predictions = []
    for i, prob_vector in enumerate(probabilities):
        mosaic_prob = prob_vector[list(class_labels).index('mosaic')] if 'mosaic' in class_labels else 0
        germline_prob = prob_vector[list(class_labels).index('germline')] if 'germline' in class_labels else 0
        repeat_prob = prob_vector[list(class_labels).index('repeat')] if 'repeat' in class_labels else 0
        
        # Relaxed decision rule
        if mosaic_prob > 0.5:  # lower mosaic threshold
            predictions.append('mosaic')
        elif germline_prob > 0.6:  # higher germline threshold
            predictions.append('germline')
        elif repeat_prob > 0.6:   # higher repeat threshold
            predictions.append('repeat')
        else:
            # If none of the above apply, further relax the mosaic cutoff
            if mosaic_prob > 0.2:  # lower mosaic threshold again
                predictions.append('mosaic')
            else:
                # Otherwise take the class with the highest probability
                predictions.append(class_labels[np.argmax(prob_vector)])
    
    # Build the results dataframe
    results_df = pd.DataFrame({
        'mutation_id': df_new['mutation_id'].values,
        'sample_id': sample_id,
        'predicted_label': predictions
    })
    
    # Add class probability scores
    for i, class_name in enumerate(class_labels):
        results_df[f'probability_{class_name}'] = probabilities[:, i]
    
    # Record the decision rule used
    results_df['decision_rule'] = 'relaxed_threshold'
    
    # Print prediction statistics
    prediction_counts = pd.Series(predictions).value_counts().to_dict()
    print(f"\nRelaxed-threshold prediction summary:")
    for label, count in prediction_counts.items():
        percentage = count / len(predictions) * 100
        print(f"  {label}: {count} sites ({percentage:.1f}%)")
    
    # Save results
    if output_file:
        results_df.to_csv(output_file, index=False)
        print(f"Predictions saved to: {output_file}")
    
    return results_df

def analyze_feature_importance(model, output_path):
    """
    Rank and report feature importances from the trained model.
    """
    print(f"\n=== Feature importance analysis ===")
    
    importances = model.feature_importances_
    feature_imp_df = pd.DataFrame({
        'feature': SELECTED_FEATURES,
        'importance': importances
    }).sort_values('importance', ascending=False)
    
    print("Feature importance ranking:")
    for _, row in feature_imp_df.iterrows():
        print(f"  {row['feature']}: {row['importance']:.4f}")
    
    return feature_imp_df

def generate_relaxed_prediction_report(results_df, sample_id, output_path):
    """
    Write a detailed report for relaxed-threshold predictions.
    """
    report_file = os.path.join(output_path, f"relaxed_prediction_report_{sample_id}.txt")
    
    with open(report_file, 'w') as f:
        f.write(f"Relaxed-threshold classifier prediction report - sample {sample_id}\n")
        f.write("=" * 50 + "\n")
        f.write(f"Total sites: {len(results_df)}\n\n")
        
        # Prediction statistics
        pred_counts = results_df['predicted_label'].value_counts()
        f.write("Prediction summary:\n")
        for label, count in pred_counts.items():
            percentage = count / len(results_df) * 100
            f.write(f"  {label}: {count} ({percentage:.1f}%)\n")
        
        f.write("\nRelaxed decision rules:\n")
        f.write("1. Predict mosaic if mosaic probability > 0.5\n")
        f.write("2. Predict germline if germline probability > 0.6\n")
        f.write("3. Predict repeat if repeat probability > 0.6\n")
        f.write("4. Otherwise prefer mosaic (probability > 0.2) to reduce false negatives\n")
        
        # Mosaic probability distribution
        if 'probability_mosaic' in results_df.columns:
            mosaic_probs = results_df['probability_mosaic']
            f.write(f"\nMosaic probability distribution:\n")
            f.write(f"  Minimum: {mosaic_probs.min():.3f}\n")
            f.write(f"  Maximum: {mosaic_probs.max():.3f}\n")
            f.write(f"  Mean: {mosaic_probs.mean():.3f}\n")
            f.write(f"  Median: {mosaic_probs.median():.3f}\n")
            
            # Counts by probability bin
            bins = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
            for i in range(len(bins)-1):
                count = ((mosaic_probs >= bins[i]) & (mosaic_probs < bins[i+1])).sum()
                f.write(f"  [{bins[i]:.1f}-{bins[i+1]:.1f}): {count} sites\n")
    
    print(f"Relaxed-prediction report written to: {report_file}")

def real_time_classifier_predict(df_for_classifier, sample_id, output_path):
    """
    Run relaxed-threshold real-time classifier prediction.
    
    Args:
        df_for_classifier: Feature dataframe to classify
        sample_id: Sample ID
        output_path: Output directory
        
    Returns:
        tuple: (results_df, model, pipeline)
    """
    print(f"\n{'='*60}")
    print(f"Starting relaxed-threshold real-time classifier prediction - sample: {sample_id}")
    print(f"{'='*60}")
    
    # Ensure the output directory exists
    os.makedirs(output_path, exist_ok=True)
    
    # 1. Load training data
    df_training = load_training_data()
    
    # 2. Check whether the sample is present in the training set
    training_samples = set(df_training['sampleid'].unique())
    if sample_id not in training_samples:
        print(f"Warning: sample {sample_id} is not in the training data")
        exclude_sample_id = sample_id
    else:
        exclude_sample_id = sample_id
        print(f"Sample {sample_id} is in the training data; using a leave-one-out strategy")
    
    # 3. Train the relaxed-threshold classifier
    model, pipeline = build_relaxed_classifier_excluding_sample(df_training, exclude_sample_id)
    
    # 4. Save the model to a pickle file
    model_file = os.path.join(output_path, f"relaxed_classifier_{sample_id}.pkl")
    joblib.dump({
        'model': model,
        'pipeline': pipeline,
        'sample_id': sample_id,
        'excluded_sample': exclude_sample_id,
        'feature_names': SELECTED_FEATURES,
        'training_date': pd.Timestamp.now()
    }, model_file)
    print(f"Model saved to: {model_file}")
    
    # 5. Predict with relaxed thresholds
    output_file = os.path.join(output_path, f"relaxed_predictions_{sample_id}.csv")
    results = predict_with_relaxed_threshold(model, pipeline, df_for_classifier, sample_id, output_file)
    
    # 6. Analyze feature importance
    feature_imp_df = analyze_feature_importance(model, output_path)
    
    # 7. Save feature importance
    feature_imp_file = os.path.join(output_path, f"feature_importance_{sample_id}.csv")
    feature_imp_df.to_csv(feature_imp_file, index=False)
    print(f"Feature importance saved to: {feature_imp_file}")
    
    # 8. Write a detailed report
    generate_relaxed_prediction_report(results, sample_id, output_path)
    
    return results, model, pipeline

def load_relaxed_classifier(model_path):
    """
    Load a previously saved relaxed-threshold classifier.
    
    Args:
        model_path: Path to the saved model file
        
    Returns:
        dict: Dictionary containing the model, pipeline, and metadata
    """
    print(f"Loading model: {model_path}")
    classifier_data = joblib.load(model_path)
    
    model = classifier_data['model']
    pipeline = classifier_data['pipeline']
    sample_id = classifier_data['sample_id']
    
    print(f"Loaded model information:")
    print(f"  Sample ID: {sample_id}")
    print(f"  Excluded sample: {classifier_data['excluded_sample']}")
    print(f"  Number of features: {len(classifier_data['feature_names'])}")
    print(f"  Training date: {classifier_data['training_date']}")
    print(f"  Model classes: {model.classes_}")
    
    return classifier_data

def predict_with_saved_classifier(df_for_classifier, model_path, sample_id, output_file=None):
    """
    Predict using a previously saved classifier.
    
    Args:
        df_for_classifier: Feature dataframe to classify
        model_path: Path to the saved model file
        sample_id: Sample ID
        output_file: Optional path to write predictions
        
    Returns:
        pd.DataFrame: Prediction results
    """
    # Load the saved model
    classifier_data = load_relaxed_classifier(model_path)
    model = classifier_data['model']
    pipeline = classifier_data['pipeline']
    
    # Predict with relaxed thresholds
    results = predict_with_relaxed_threshold(model, pipeline, df_for_classifier, sample_id, output_file)
    
    return results


# # Usage example
# def main():
#     """
#     Main function - usage example
#     """
#     # Example feature dataframe
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
    
#     # Method 1: train a classifier in real time and predict
#     print("=== Method 1: train a classifier in real time ===")
#     results, model, pipeline = real_time_classifier_predict(df_for_classifier, sample_id, output_path)
    
#     print(f"\nPreview of final predictions:")
#     print(results[['mutation_id', 'predicted_label', 'probability_mosaic']].head())
    
#     # Method 2: predict with a saved model
#     print("\n=== Method 2: use a saved model ===")
#     model_path = os.path.join(output_path, f"relaxed_classifier_{sample_id}.pkl")
#     if os.path.exists(model_path):
#         results_saved = predict_with_saved_classifier(
#             df_for_classifier, 
#             model_path, 
#             sample_id,
#             output_file=os.path.join(output_path, f"predictions_with_saved_model_{sample_id}.csv")
#         )
#         print(f"Preview of predictions from the saved model:")
#         print(results_saved[['mutation_id', 'predicted_label', 'probability_mosaic']].head())
#     else:
#         print(f"Model file not found: {model_path}")

# if __name__ == "__main__":
#     main()
