# Compute dynverse Hamming-Ipsen-Mikhailov (HIM) Trajectory Topology Similarity

Quantifies both local edge correspondence (Hamming) and global Laplacian
spectral similarity (Ipsen-Mikhailov) between reference and simulated
trajectory milestone networks.

## Usage

``` r
calc_dynverse_him(ref_traj, sim_traj)
```

## Arguments

- ref_traj:

  Reference trajectory list returned by
  [`infer_dynverse_trajectory`](https://kabilanbio.github.io/scSimEval/reference/infer_dynverse_trajectory.md).

- sim_traj:

  Simulated trajectory list returned by
  [`infer_dynverse_trajectory`](https://kabilanbio.github.io/scSimEval/reference/infer_dynverse_trajectory.md).

## Value

Numeric HIM similarity in `[0, 1]` (higher is better).
