# Version history

[Back to README](../README.md)

- **v3.6.0** (2026-10-28): Aligned resolution-backoff and pruning safeguards across modes; cleaned pre-QC/post-QC tree exports with explicit fallback metadata; consistent ROOT mutation scoring; refreshed demos and reorganized documentation with a complete parameter reference ([changelog](https://github.com/douymLab/PhyloSOLID/releases/tag/v3.6.0))
- **v3.5.0** (2026-09-20): scDNA BAM-to-tree end-to-end packaging (bulk BAM optional); self-contained scRNA / scDNA / binary-matrix demos; cleaned vs unpruned export consistency; dedicated visualization environment (`install_vis.sh`, `environment_vis.yml`); English localization of pipeline scripts and docs ([changelog](https://github.com/douymLab/PhyloSOLID/releases/tag/v3.5.0))
- **v3.4.0** (2026-09-03): Adaptive Leiden graph resolution and pruning confidence CV selection: automatic resolution tuning based on data sparsity; new `pruning_confidence` criterion prevents over-pruning by evaluating `pruning_ratio = (Ω_pre-QC - Ω_final) / Ω_final < 10.0`; cleaner output structure with reduced file redundancy ([changelog](https://github.com/douymLab/PhyloSOLID/releases/tag/v3.4.0))
- **v3.3.1** (2026-08-21): Deterministic parallel CV search: fixed non-determinism where same CV threshold produced different results in parallel vs single-run mode ([changelog](https://github.com/douymLab/PhyloSOLID/releases/tag/v3.3.1))
- **v3.3.0** (2026-08-20): Improved optimal tree selection with `Omega_final` criterion and robust tie-breaking mechanism ([changelog](https://github.com/douymLab/PhyloSOLID/releases/tag/v3.3.0))
- **v3.2.1** (2026-08-20): Critical bug fixes for CV threshold search: fixed logger `NameError`, eliminated duplicate mutation nodes in scaffold tree, standardized logger parameter passing ([changelog](https://github.com/douymLab/PhyloSOLID/releases/tag/v3.2.1))
- **v3.2.0** (2026-08-15): Parallel CV traversal, unified optimal selection, structured logging overhaul, full English localization ([changelog](https://github.com/douymLab/PhyloSOLID/releases/tag/v3.2.0))
- **v3.1.1** (2026-08-08): Performance optimization: 4-33x speedup in penalty calculation with caching and vectorization ([changelog](https://github.com/douymLab/PhyloSOLID/releases/tag/v3.1.1))
- **v3.1.0** (2026-07-28): CV threshold auto-search, organized output structure, new `full_tree_builder.py` module ([changelog](https://github.com/douymLab/PhyloSOLID/releases/tag/v3.1.0))
- **v3.0.0** (2026-07-02): Smart conflict resolution, Leiden graph acceleration, multi-threaded data loading ([changelog](https://github.com/douymLab/PhyloSOLID/releases/tag/v3.0.0))
- **v2.0.0** (2026-07-01): Major performance optimization ([changelog](https://github.com/douymLab/PhyloSOLID/releases/tag/v2.0.0))
- **v1.0.0** (2026-02-27): Initial public release ([changelog](https://github.com/douymLab/PhyloSOLID/releases/tag/v1.0.0))

For complete release notes, visit the [Releases page](https://github.com/douymLab/PhyloSOLID/releases).
