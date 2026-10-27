"""Experimental bounded resolution backoff on a fixed initial scaffold graph.

The conflict fraction is a singlet/main-clone proxy, not a doublet diagnosis.
No cells are removed here. An explicit risk budget is required; no universal
default is claimed. NA observations never count as positive or negative calls.
"""
import math
import numpy as np
import pandas as pd


def partition_conflicts(I_matrix, mutation_group, fp_cutoff=0.5):
    """Evaluate the same grouped mutations/cells for each candidate partition."""
    if not math.isfinite(fp_cutoff) or not 0 < fp_cutoff <= 1:
        raise ValueError('fp_cutoff must be in (0, 1]')
    names = list(mutation_group)
    if not I_matrix.index.is_unique or not I_matrix.columns.is_unique:
        raise ValueError('Cell and mutation identifiers must be unique')
    missing = set(names) - set(I_matrix.columns)
    if missing:
        raise ValueError('Grouped mutations missing from I_matrix: ' + str(sorted(missing)))
    observations = I_matrix.loc[:, names]
    if ((observations.notna()) & ~observations.isin([0, 1])).to_numpy().any():
        raise ValueError('I_matrix must contain only 0, 1 and missing values')
    ones = observations.eq(1).fillna(False).to_numpy(dtype=np.int64)
    labels = list(dict.fromkeys(mutation_group.values()))
    group_counts = np.zeros((len(I_matrix), len(labels)), dtype=np.int64)
    assignments = np.array([labels.index(mutation_group[m]) for m in names])
    for column in range(len(labels)):
        group_counts[:, column] = ones[:, assignments == column].sum(axis=1)
    total = ones.sum(axis=1)
    maximum = group_counts.max(axis=1) if labels else np.zeros(len(I_matrix), dtype=np.int64)
    occupied = (group_counts > 0).sum(axis=1)
    informative = total > 0
    minimum_fp = total - maximum
    ratios = np.full(len(I_matrix), np.nan)
    np.divide(minimum_fp, total, out=ratios, where=informative)
    high_conflict = informative & (ratios >= fp_cutoff)
    n_informative = int(informative.sum())
    summary = dict(
        n_input_cells=len(I_matrix), n_grouped_mutations=len(names),
        n_groups=len(labels), n_informative_cells=n_informative,
        n_multi_group_cells=int((occupied >= 2).sum()),
        multi_group_fraction=float((occupied >= 2).sum()/n_informative) if n_informative else None,
        n_high_conflict_cells=int(high_conflict.sum()),
        high_conflict_fraction=float(high_conflict.sum()/n_informative) if n_informative else None,
        fp_cutoff=float(fp_cutoff), total_positive_calls=int(total.sum()),
        minimum_out_of_group_fp=int(minimum_fp.sum()),
    )
    per_cell = pd.DataFrame(dict(
        positive_calls=total, occupied_groups=occupied, largest_group_calls=maximum,
        minimum_out_of_group_fp=minimum_fp, fp_proxy=ratios,
        informative=informative, high_conflict=high_conflict,
    ), index=I_matrix.index)
    per_cell.index.name = 'cell_id'
    return summary, per_cell


def refine_resolution(graph, I_matrix, initial_partition, initial_resolution,
                      threshold, min_resolution=0.5, step=0.1, max_steps=20,
                      seed=42, fp_cutoff=0.5, partitioner=None, logger_obj=None):
    """Accept the first downward candidate satisfying H <= threshold.

    The graph, vertex set and denominator remain fixed. H need not be monotone.
    Bounds guarantee finite termination, not optimality. If no candidate meets
    the budget, retain the original auto partition and report threshold_unmet.
    `partitioner` is an injectable optimisation backend for focused tests.
    """
    for name, value in [('threshold', threshold), ('initial_resolution', initial_resolution),
                        ('min_resolution', min_resolution), ('step', step)]:
        if value is None or not math.isfinite(value):
            raise ValueError(name + ' must be finite')
    if not 0 <= threshold <= 1:
        raise ValueError('threshold must be in [0, 1]')
    if initial_resolution <= 0 or min_resolution <= 0 or step <= 0:
        raise ValueError('Resolution values and step must be positive')
    if not isinstance(max_steps, int) or max_steps < 0:
        raise ValueError('max_steps must be a nonnegative integer')
    names = list(graph.vs['name'])
    if len(set(names)) != len(names):
        raise ValueError('Graph vertex identifiers must be unique')
    if partitioner is None:
        import leidenalg
        def partitioner(g, resolution, rng_seed):
            return leidenalg.find_partition(
                g, leidenalg.RBConfigurationVertexPartition, weights='weight',
                resolution_parameter=resolution, seed=rng_seed)

    def mapping(partition):
        membership = list(partition.membership)
        if len(membership) != len(names):
            raise ValueError('Partition membership does not match graph vertices')
        return dict(zip(names, membership))

    initial_group = mapping(initial_partition)
    selected_partition = initial_partition
    selected_group = initial_group
    selected_resolution = float(initial_resolution)
    trace = []
    partition = initial_partition
    resolution = float(initial_resolution)
    initial_cells = None
    final_cells = None
    threshold_met = False
    reason = None
    reference_support = None
    for iteration in range(max_steps + 1):
        group = mapping(partition)
        metrics, per_cell = partition_conflicts(I_matrix, group, fp_cutoff=fp_cutoff)
        support = (metrics['n_informative_cells'], metrics['total_positive_calls'])
        if reference_support is None:
            reference_support = support
            initial_cells = per_cell
        elif support != reference_support:
            raise RuntimeError('Candidate partition changed the evaluation population')
        record = dict(iteration=iteration, resolution=float(resolution), **metrics)
        trace.append(record)
        if logger_obj:
            logger_obj.info('Resolution backoff: gamma=%.3f groups=%d H=%s budget=%.3f',
                            resolution, metrics['n_groups'], metrics['high_conflict_fraction'], threshold)
        if not metrics['n_informative_cells']:
            reason = 'no_informative_cells'
            break
        if metrics['high_conflict_fraction'] <= threshold:
            selected_partition, selected_group = partition, group
            selected_resolution = resolution
            final_cells = per_cell
            threshold_met = True
            reason = 'initial_budget_met' if iteration == 0 else 'backoff_budget_met'
            break
        if resolution <= min_resolution:
            reason = 'resolution_floor_threshold_unmet'
            break
        if iteration == max_steps:
            reason = 'max_steps_threshold_unmet'
            break
        next_resolution = round(max(min_resolution, initial_resolution-(iteration+1)*step), 10)
        if next_resolution >= resolution:
            raise RuntimeError('Backoff must strictly decrease resolution')
        resolution = next_resolution
        partition = partitioner(graph, resolution, seed)
    if final_cells is None:
        final_cells = initial_cells
    if logger_obj:
        logger_obj.info('Resolution backoff stopped: %s; selected gamma=%.3f; met=%s',
                        reason, selected_resolution, threshold_met)
    return dict(partition=selected_partition, mutation_group=selected_group,
                resolution=float(selected_resolution), initial_resolution=float(initial_resolution),
                threshold=float(threshold), threshold_met=threshold_met, stop_reason=reason,
                trace=trace, initial_cells=initial_cells, final_cells=final_cells)
