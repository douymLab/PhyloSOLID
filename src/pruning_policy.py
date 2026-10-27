"""Shared policy for zero-discordance and empty post-QC outputs."""
import math

PRUNING_EPSILON = 1e-8

def pruning_ratio(omega_pre, omega_final):
    pre, final = float(omega_pre), float(omega_final)
    if not math.isfinite(pre) or not math.isfinite(final) or pre < 0 or final < 0:
        raise ValueError("Omega must be finite and nonnegative")
    return (pre - final) / max(final, PRUNING_EPSILON)

def has_observed_tree_content(matrix):
    """ROOT alone and all-zero mutation matrices do not constitute an output."""
    mutations = matrix.drop(columns=['ROOT'], errors='ignore')
    return not mutations.empty and bool(mutations.eq(1).any().any())
