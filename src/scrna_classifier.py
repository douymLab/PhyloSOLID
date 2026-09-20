#!/usr/bin/env python3

###################################################################################################
######################### Core functions for the relaxed scRNA classifier #########################
###################################################################################################

import os
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
import joblib
from collections import Counter
from pathlib import Path
from src.reproducibility import set_seed, deterministic_permutation


# Feature columns used by the classifier
SELECTED_FEATURES = [
    'falt', 
    'mutant_cell_fraction', 
    'vaf_in_pseudobulk',
    'vaf_mutant_avg', 
    'right_ref_querypos_num_remove_clip',
    'alt_UMI_avg_consistence_remove_single_read',
    'alt_mismatches_mean',
    'alt_consistence_soft_prop', 
    'alt_dp_mean_diff',
    'norm_alt_dp_in_pseudobulk', 
    'mismatches_p_adj',
    'sig_pvalue', 
    'alt_read_number_perUMI_median',
    'baseq_p_adj'
]

def build_relaxed_classifier_excluding_sample(df_training, exclude_sample_id):
    """
    Train a relaxed-threshold classifier while excluding a target sample.
    """
    print(f"\n=== Building relaxed-threshold classifier (excluding sample {exclude_sample_id}) ===")
    
    # Exclude the target sample (leave-one-out)
    df_train = df_training[df_training['sampleid'] != exclude_sample_id].copy()
    print(f"Training data size: {df_train.shape}")
    print(f"Training class distribution: {Counter(df_train['label2'])}")
    
    # Prepare features and labels
    X_train = df_train[SELECTED_FEATURES].copy()
    y_train = df_train['label2'].copy()
    
    # Build the preprocessing pipeline
    pipeline = make_pipeline(
        SimpleImputer(strategy='median'),
        StandardScaler()
    )
    
    # Preprocess data and train the model
    X_train_processed = pipeline.fit_transform(X_train)
    
    # Relaxed random-forest settings
    rf_model = RandomForestClassifier(
        n_estimators=1000,
        random_state=42,
        class_weight={'mosaic': 1.0, 'artifact': 0.8},
        max_depth=15,
        min_samples_split=2,
        min_samples_leaf=1,
        max_features='sqrt',
        bootstrap=True
    )
    rf_model.fit(X_train_processed, y_train)
    
    print(f"Model training completed. Classes: {rf_model.classes_}")
    return rf_model, pipeline

def predict_with_relaxed_threshold(model, pipeline, df_new, sample_id, mutation_ids=None):
    """
    Predict labels using relaxed probability thresholds.
    """
    print(f"\n=== Predicting mutations for sample {sample_id} (relaxed thresholds) ===")
    print(f"New data shape: {df_new.shape}")
    
    # Check required feature columns
    missing_features = [feat for feat in SELECTED_FEATURES if feat not in df_new.columns]
    if missing_features:
        raise ValueError(f"Missing required feature columns: {missing_features}")
    
    # Prepare feature matrix
    X_new = df_new[SELECTED_FEATURES].copy()
    
    # Apply the same preprocessing pipeline
    X_new_processed = pipeline.transform(X_new)
    
    # Predict using class probabilities rather than hard labels
    probabilities = model.predict_proba(X_new_processed)
    class_labels = model.classes_
    
    # Relaxed prediction logic
    predictions = []
    for i, prob_vector in enumerate(probabilities):
        mosaic_prob = prob_vector[list(class_labels).index('mosaic')] if 'mosaic' in class_labels else 0
        artifact_prob = prob_vector[list(class_labels).index('artifact')] if 'artifact' in class_labels else 0
        
        # Relaxed decision rule
        if mosaic_prob > 0.4:  # lower mosaic threshold
            predictions.append('mosaic')
        elif artifact_prob > 0.6:  # higher artifact threshold
            predictions.append('artifact')
        else:
            if mosaic_prob > 0.2:  # further lower the mosaic threshold
                predictions.append('mosaic')
            else:
                predictions.append(class_labels[np.argmax(prob_vector)])
    
    # Build the results dataframe
    results_df = pd.DataFrame({
        'sample_id': sample_id,
        'predicted_label': predictions
    })
    
    # Add mutation_id
    if mutation_ids is not None:
        results_df['mutation_id'] = mutation_ids
    elif 'mutation_id' in df_new.columns:
        results_df['mutation_id'] = df_new['mutation_id'].values
    elif 'identifier' in df_new.columns:
        results_df['mutation_id'] = df_new['identifier'].values
    else:
        # If no mutation_id column exists, use a positional identifier
        results_df['mutation_id'] = [f'mutation_{i+1}' for i in range(len(results_df))]
    
    # Add class probability scores
    for i, class_name in enumerate(class_labels):
        results_df[f'probability_{class_name}'] = probabilities[:, i]
    
    # Put mutation_id first
    cols = ['mutation_id', 'sample_id', 'predicted_label'] + [f'probability_{cls}' for cls in class_labels]
    results_df = results_df[cols]
    
    # Print prediction statistics
    prediction_counts = pd.Series(predictions).value_counts().to_dict()
    print(f"Relaxed-threshold prediction summary:")
    for label, count in prediction_counts.items():
        percentage = count / len(predictions) * 100
        print(f"  {label}: {count} sites ({percentage:.1f}%)")
    
    return results_df

def real_time_classifier_predict(df_for_classifier_all, sampleid, outputpath):
    """
    Core prediction function for the relaxed-threshold scRNA classifier.
    
    Args:
        df_for_classifier_all: Feature dataframe to classify
        sampleid: Sample ID
        outputpath: Output directory
        
    Returns:
        pd.DataFrame: Prediction results
    """
    print(f"\n{'='*60}")
    print(f"Starting relaxed-threshold scRNA classifier prediction - sample: {sampleid}")
    print(f"{'='*60}")
    
    # Ensure the output directory exists
    os.makedirs(outputpath, exist_ok=True)
    
    # 1. Load training data
    print("=== Loading training data ===")
    script_dir = Path(__file__).parent
    features_file_labeled = script_dir / 'classifier' / 'scrna' / 'data_labeling_sampling.ratio_155_space.csv'
    df_training = pd.read_csv(features_file_labeled, sep="\t")
    print(f"Training data shape: {df_training.shape}")
    print(f"Training class distribution: {Counter(df_training['label2'])}")
    
    # 2. Preprocess features
    print("=== Preprocessing data ===")
    df_features = df_for_classifier_all.copy()
    
    # Extract mutation_id if a relevant column exists
    mutation_ids = None
    if 'mutation_id' in df_features.columns:
        mutation_ids = df_features['mutation_id'].values
    elif 'identifier' in df_features.columns:
        mutation_ids = df_features['identifier'].values
    
    # Clean feature values
    df_features_selected = df_features[SELECTED_FEATURES].copy()
    df_features_selected.replace('no', np.nan, inplace=True)
    
    for col in df_features_selected.columns:
        df_features_selected[col] = pd.to_numeric(df_features_selected[col], errors='coerce')
        median_value = df_features_selected[col].median(skipna=True)
        df_features_selected[col] = df_features_selected[col].fillna(median_value)
        finite_max = df_features_selected[col][np.isfinite(df_features_selected[col])].max()
        df_features_selected[col] = df_features_selected[col].replace([np.inf, -np.inf], finite_max)
    
    print(f"Preprocessed data shape: {df_features_selected.shape}")
    
    # 3. Train the relaxed-threshold classifier
    model, pipeline = build_relaxed_classifier_excluding_sample(df_training, sampleid)
    
    # 4. Predict with relaxed thresholds
    results = predict_with_relaxed_threshold(model, pipeline, df_features_selected, sampleid, mutation_ids)
    
    # 5. Save prediction results
    print("\n=== Saving prediction results ===")
    
    # Save predictions for all sites
    all_sites_file = os.path.join(outputpath, f"{sampleid}.feature_and_prediction.allsites.txt")
    results.to_csv(all_sites_file, index=False, sep="\t")
    print(f"All-site predictions saved to: {all_sites_file}")
    
    # Save mosaic site list
    df_mosaic = results[results["predicted_label"] == "mosaic"]
    mosaic_list_file = os.path.join(outputpath, f"{sampleid}_mosaic_prediction.list.txt")
    df_mosaic['mutation_id'].to_csv(mosaic_list_file, index=False, header=False)
    print(f"Mosaic site list saved to: {mosaic_list_file}")
    
    # Save detailed mosaic results
    mosaic_detailed_file = os.path.join(outputpath, f"{sampleid}_mosaic_detailed_results.txt")
    df_mosaic.to_csv(mosaic_detailed_file, index=False, sep="\t")
    print(f"Detailed mosaic results saved to: {mosaic_detailed_file}")
    
    # Save summary statistics
    mosaic_count_file = os.path.join(outputpath, f"{sampleid}_prediction_summary.txt")
    with open(mosaic_count_file, 'w') as f:
        f.write(f"Prediction summary for sample {sampleid}\n")
        f.write("=" * 40 + "\n")
        f.write(f"Total sites: {len(results)}\n")
        f.write(f"Mosaic sites: {len(df_mosaic)}\n")
        f.write(f"Artifact sites: {len(results) - len(df_mosaic)}\n")
        f.write(f"Mosaic fraction: {len(df_mosaic)/len(results)*100:.1f}%\n")
        
        # Add probability statistics
        if 'probability_mosaic' in results.columns:
            f.write(f"\nMosaic probability statistics:\n")
            f.write(f"  Mean: {results['probability_mosaic'].mean():.3f}\n")
            f.write(f"  Median: {results['probability_mosaic'].median():.3f}\n")
            f.write(f"  Maximum: {results['probability_mosaic'].max():.3f}\n")
            f.write(f"  Minimum: {results['probability_mosaic'].min():.3f}\n")
    
    print(f"Prediction summary saved to: {mosaic_count_file}")
    
    # Print final statistics
    mosaic_count = len(df_mosaic)
    total_count = len(results)
    print(f"\nFinal prediction summary: {mosaic_count}/{total_count} sites predicted as mosaic ({mosaic_count/total_count*100:.1f}%)")
    
    print(f"\n{'='*60}")
    print(f"Relaxed-threshold scRNA classifier prediction completed - sample: {sampleid}")
    print(f"{'='*60}")
    
    return results

# # Usage example
# if __name__ == "__main__":
#     # Example feature dataframe (including mutation_id)
#     df_for_classifier_all = pd.DataFrame({
#         'mutation_id': [
#             'chr20_14326432_C_A', 'chr21_9001082_G_C', 'chr6_137868589_G_A',
#             'chr1_106340339_G_T', 'chr1_212476040_A_G', 'chr19_38404452_C_T'
#         ],
#         'falt': [0.1, 0.2, 0.05, 0.15, 0.08, 0.12],
#         'mutant_cell_fraction': [0.3, 0.4, 0.1, 0.35, 0.25, 0.28],
#         'vaf_in_pseudobulk': [0.01, 0.02, 0.005, 0.015, 0.008, 0.012],
#         'vaf_mutant_avg': [0.15, 0.25, 0.08, 0.18, 0.12, 0.16],
#         'right_ref_querypos_num_remove_clip': [5, 8, 2, 6, 4, 5],
#         'alt_UMI_avg_consistence_remove_single_read': [0.1, 0.15, 0.05, 0.12, 0.08, 0.11],
#         'alt_mismatches_mean': [0.02, 0.03, 0.01, 0.025, 0.015, 0.022],
#         'alt_consistence_soft_prop': [0.05, 0.08, 0.02, 0.06, 0.04, 0.055],
#         'alt_dp_mean_diff': [0.12, 0.22, 0.06, 0.16, 0.1, 0.14],
#         'norm_alt_dp_in_pseudobulk': [0.01, 0.02, 0.005, 0.015, 0.008, 0.012],
#         'mismatches_p_adj': [0.05, 0.08, 0.02, 0.06, 0.04, 0.055],
#         'sig_pvalue': [0.01, 0.02, 0.005, 0.015, 0.008, 0.012],
#         'alt_read_number_perUMI_median': [2, 3, 1, 2, 2, 2],
#         'baseq_p_adj': [0.05, 0.08, 0.02, 0.06, 0.04, 0.055]
#     })
    
#     sampleid = "10k"
#     outputpath = "./scRNA_relaxed_predictions"
    
#     # Run prediction
#     results = real_time_classifier_predict(df_for_classifier_all, sampleid, outputpath)
    
#     print(f"\nPreview of predictions:")
#     print(results.head())
