# Local mechanical modules embedded in the distributed Ca network

Spatial modules were recomputed from the same current CFU HDF5 inputs used by Fig5/23 global co-occurrence, with ratio<=3, coverage>=0.5, and pattern members>=5.

A module is a local pattern-CFU spatial pair.  For every module pattern, all significant co-occurrence CFUs are joined as candidate network inputs.  A non-local CFU is labelled candidate_upstream only when best_lag<0, because the Fig5/23 convention is CFU onset = motion activation peak + best_lag.  Lag 0 is synchronous and positive lag is candidate_downstream; these are association labels, not causal proof.
