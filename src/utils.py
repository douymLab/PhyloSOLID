import os
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist, squareform
from collections import defaultdict
import random

_BULK_BARCODES = {"bulk", "pseudo_bulk", "ROOT"}


def save_celltype_table(df_celltype: pd.DataFrame, path: str) -> str:
    """Write cell-type annotation as TSV without a pandas row index."""
    df = df_celltype.copy()
    df.columns = df.columns.astype(str)
    df = df.loc[:, ~df.columns.str.match(r"^Unnamed")]
    if "barcode" in df.columns:
        df = df[~df["barcode"].astype(str).isin(_BULK_BARCODES)].copy()
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    return path

class BootstrapSupportCalculator:
    """
    Compute bootstrap support for phylogenetic-tree branches.
    """
    
    def __init__(self, mutation_matrix, mutation_names=None, n_bootstrap=1000, random_seed=42):
        """
        Parameters:
        - mutation_matrix: mutation matrix, rows=mutations, columns=samples/cells, values=0/1 (1 means the mutation is present)
        - mutation_names: mutation name list, e.g. ['M1', 'M2', ...]
        - n_bootstrap: number of bootstrap replicates
        - random_seed: random seed for reproducibility
        """
        self.original_matrix = np.array(mutation_matrix)
        self.n_mutations = self.original_matrix.shape[0]
        self.n_samples = self.original_matrix.shape[1]
        self.mutation_names = mutation_names if mutation_names else [f'M{i+1}' for i in range(self.n_mutations)]
        self.n_bootstrap = n_bootstrap
        self.random_seed = random_seed
        
        # Build the original tree and identify first-level branches
        self.original_tree = None
        self.main_branches = None  # mutation sets for first-level branches
        
    def build_tree_and_get_branches(self, distance_matrix=None, threshold=0.5):
        """
        Build a tree from a distance matrix and identify first-level branches.
        
        Parameters:
        - distance_matrix: pairwise mutation distance matrix (if None, compute from co-occurrence)
        - threshold: clustering threshold used to define first-level branches
        """
        if distance_matrix is None:
            # Compute pairwise mutation distances with Jaccard distance
            distance_matrix = self._compute_jaccard_distance()
        
        # Hierarchical clustering with UPGMA
        condensed_dist = squareform(distance_matrix)
        linkage_matrix = linkage(condensed_dist, method='average')
        
        # Cut the tree at the threshold to obtain first-level branches
        # If unspecified, use a multiple of the median of positive distances
        if threshold is None:
            threshold = np.median(distance_matrix[distance_matrix > 0]) * 0.8
        
        cluster_labels = fcluster(linkage_matrix, t=threshold, criterion='distance')
        
        # Extract mutations in each branch
        branches = defaultdict(list)
        for i, label in enumerate(cluster_labels):
            branches[label].append(i)
        
        # Keep only branches that contain at least two mutations
        self.main_branches = {f'Branch_{k}': sorted(indices) 
                             for k, indices in branches.items() if len(indices) >= 2}
        
        return self.main_branches
    
    def _compute_jaccard_distance(self):
        """Compute pairwise Jaccard distances among mutations."""
        n = self.n_mutations
        distance_matrix = np.zeros((n, n))
        
        for i in range(n):
            for j in range(i+1, n):
                # Compute co-occurrence of the two mutations
                intersection = np.sum(self.original_matrix[i] & self.original_matrix[j])
                union = np.sum(self.original_matrix[i] | self.original_matrix[j])
                
                if union == 0:
                    distance = 1.0  # neither mutation is present; maximum distance
                else:
                    jaccard_sim = intersection / union
                    distance = 1 - jaccard_sim  # Jaccard distance
                
                distance_matrix[i, j] = distance
                distance_matrix[j, i] = distance
        
        return distance_matrix
    
    def bootstrap_resample(self):
        """
        Sample mutations with replacement to generate a bootstrap replicate.
        """
        # Draw with replacement, matching the original number of mutations
        sampled_indices = random.choices(range(self.n_mutations), k=self.n_mutations)
        bootstrap_matrix = self.original_matrix[sampled_indices]
        return bootstrap_matrix, sampled_indices
    
    def compute_branch_support(self):
        """
        Compute bootstrap support for each first-level branch.
        """
        # First build the original tree and identify first-level branches
        if self.main_branches is None:
            self.build_tree_and_get_branches()
        
        print(f"Identified {len(self.main_branches)} first-level branches:")
        for branch_name, mutations in self.main_branches.items():
            mutation_labels = [self.mutation_names[i] for i in mutations]
            print(f"  {branch_name}: {mutation_labels}")
        
        # Store support counts for each branch
        support_counts = {branch_name: 0 for branch_name in self.main_branches}
        
        # Set random seeds
        random.seed(self.random_seed)
        np.random.seed(self.random_seed)
        
        print(f"\nStarting bootstrap analysis ({self.n_bootstrap} replicates)...")
        
        for bootstrap_iter in range(self.n_bootstrap):
            if (bootstrap_iter + 1) % 100 == 0:
                print(f"  Completed {bootstrap_iter + 1}/{self.n_bootstrap} replicates")
            
            # Generate a bootstrap sample
            bootstrap_matrix, sampled_indices = self.bootstrap_resample()
            
            # Check whether each main branch still appears in the bootstrap tree
            for branch_name, original_mutations in self.main_branches.items():
                if self._check_branch_exists(bootstrap_matrix, original_mutations, sampled_indices):
                    support_counts[branch_name] += 1
        
        # Convert counts to support percentages
        support_values = {}
        for branch_name, count in support_counts.items():
            support_percent = (count / self.n_bootstrap) * 100
            support_values[branch_name] = support_percent
        
        return support_values
    
    def _check_branch_exists(self, bootstrap_matrix, original_mutations, sampled_indices):
        """
        Check whether mutations from an original branch still cluster together in the bootstrap tree.
        
        Parameters:
        - bootstrap_matrix: mutation matrix after bootstrap resampling
        - original_mutations: mutation indices in the original branch
        - sampled_indices: original indices drawn during resampling
        """
        # Get rows in the bootstrap matrix that correspond to the original mutations
        # Note: rows in the bootstrap matrix correspond to indices in sampled_indices
        
        # First find which original mutations were sampled in this bootstrap replicate
        # and their positions in the bootstrap matrix
        mutation_positions = []
        for orig_idx in original_mutations:
            # Find all positions in sampled_indices that equal orig_idx
            positions = [i for i, idx in enumerate(sampled_indices) if idx == orig_idx]
            if positions:
                mutation_positions.extend(positions)
        
        # If too many branch mutations were lost from the bootstrap sample, treat the branch as absent
        if len(mutation_positions) < len(original_mutations) * 0.5:
            return False
        
        # Compute distances among these mutations in the bootstrap matrix
        if len(mutation_positions) < 2:
            return False
        
        # Extract data for these mutations
        sub_matrix = bootstrap_matrix[mutation_positions]
        
        # Compute the mean pairwise Jaccard distance among them
        distances = []
        for i in range(len(sub_matrix)):
            for j in range(i+1, len(sub_matrix)):
                intersection = np.sum(sub_matrix[i] & sub_matrix[j])
                union = np.sum(sub_matrix[i] | sub_matrix[j])
                if union == 0:
                    dist = 1.0
                else:
                    dist = 1 - (intersection / union)
                distances.append(dist)
        
        avg_distance = np.mean(distances) if distances else 1.0
        
        # Compute the mean distance between these mutations and all other mutations
        other_distances = []
        all_positions = list(range(len(bootstrap_matrix)))
        other_positions = [p for p in all_positions if p not in mutation_positions]
        
        if other_positions:
            for pos in mutation_positions:
                for other in other_positions:
                    intersection = np.sum(bootstrap_matrix[pos] & bootstrap_matrix[other])
                    union = np.sum(bootstrap_matrix[pos] | bootstrap_matrix[other])
                    if union == 0:
                        dist = 1.0
                    else:
                        dist = 1 - (intersection / union)
                    other_distances.append(dist)
            
            avg_other_distance = np.mean(other_distances) if other_distances else 1.0
            
            # Treat the branch as present if within-branch distance is smaller than between-branch distance
            return avg_distance < avg_other_distance
        else:
            # If there are no other mutations to compare against, treat the branch as present
            return True
    
    def print_results(self, support_values):
        """
        Print bootstrap support results.
        """
        print("\n" + "="*60)
        print("Bootstrap branch-support results")
        print("="*60)
        
        for branch_name, support in sorted(support_values.items(), key=lambda x: x[1], reverse=True):
            mutations = [self.mutation_names[i] for i in self.main_branches[branch_name]]
            stars = '*' * int(support / 5)  # one star per 5%
            print(f"{branch_name}: {support:.1f}%  {stars}")
            print(f"  Mutations: {mutations}")
            print()
        
        # Support summary
        high_support = sum(1 for v in support_values.values() if v >= 70)
        medium_support = sum(1 for v in support_values.values() if 50 <= v < 70)
        low_support = sum(1 for v in support_values.values() if v < 50)
        
        print("-"*60)
        print(f"Summary:")
        print(f"  High support (>=70%): {high_support} branches")
        print(f"  Medium support (50-69%): {medium_support} branches")
        print(f"  Low support (<50%): {low_support} branches")
        print("="*60)


# ============================================================
# Example usage
# ============================================================

def example_usage():
    """
    Example: how to use BootstrapSupportCalculator.
    """
    
    # Create example data
    # Assume 30 mutations (M1-M30)
    # 3 main branches: Branch1 (M1-M15), Branch2 (M16-M22), Branch3 (M23-M30)
    # 10 samples/cells
    
    np.random.seed(42)
    
    n_mutations = 30
    n_samples = 10
    
    # Create the mutation matrix
    mutation_matrix = np.zeros((n_mutations, n_samples), dtype=int)
    
    # Generate characteristic patterns for each branch
    # Branch1: M1-M15 present in samples 0-4
    for i in range(0, 15):
        mutation_matrix[i, :5] = np.random.binomial(1, 0.8, 5)
        mutation_matrix[i, 5:] = np.random.binomial(1, 0.1, 5)
    
    # Branch2: M16-M22 present in samples 3-7
    for i in range(15, 22):
        mutation_matrix[i, 3:8] = np.random.binomial(1, 0.8, 5)
        mutation_matrix[i, :3] = np.random.binomial(1, 0.1, 3)
        mutation_matrix[i, 8:] = np.random.binomial(1, 0.1, 2)
    
    # Branch3: M23-M30 present in samples 6-9
    for i in range(22, 30):
        mutation_matrix[i, 6:10] = np.random.binomial(1, 0.8, 4)
        mutation_matrix[i, :6] = np.random.binomial(1, 0.1, 6)
    
    # Add some random noise
    noise_mask = np.random.random((n_mutations, n_samples)) < 0.05
    mutation_matrix[noise_mask] = 1 - mutation_matrix[noise_mask]
    
    # Mutation names
    mutation_names = [f'M{i+1}' for i in range(n_mutations)]
    
    # Create the calculator
    calculator = BootstrapSupportCalculator(
        mutation_matrix=mutation_matrix,
        mutation_names=mutation_names,
        n_bootstrap=1000,  # can be set to 100 for a quicker test
        random_seed=42
    )
    
    # Compute bootstrap support
    support_values = calculator.compute_branch_support()
    
    # Print results
    calculator.print_results(support_values)
    
    return calculator, support_values


if __name__ == "__main__":
    # Run the example
    calculator, support_values = example_usage()
