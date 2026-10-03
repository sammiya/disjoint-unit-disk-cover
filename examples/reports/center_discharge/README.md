# Center discharge comparison

Measured on 2026-10-03 using commit `181cfa8`, Apple M5 Pro, 64 GiB RAM,
macOS 26.4, and Apple Clang 21.0.0. The Release build uses `-O3 -DNDEBUG`.
The solver is single-threaded; times include enumeration and search.
Both modes use the same binary. The only switch is `--center-discharge`.
No experimental partition pruning, symmetry reduction, or parallel search is used.

## Results

| Configuration | Original search | With discharge | Result with discharge |
| --- | ---: | ---: | --- |
| 44 points, a = 0.2531 | >180 s (stopped) | 24.42 s | NonCoverable, all 163,137 partitions excluded |
| 44 points, a = 0.2532 | 3.07 s | 2.39 s | Coverable |
| 44 points, a = 0.2500 | 7.15 s | 10.22 s | NonCoverable, all 157,010 partitions excluded |
| 43 points, a = h = 0.2500 | 8.85 s | 7.09 s | NonCoverable, all 90,493 partitions excluded |

Times are medians of three runs, except the 43-point case (two runs) and the
180-second timeout (one run). The three discharged 0.2531 runs took 24.42,
24.03, and 24.65 seconds. The original 0.2531 run did not complete; this
comparison gives a speedup lower bound of about 7.4, not a measured full-run ratio.
The optimization has an overhead on the easier 44-point 0.2500 instance.

The 0.2531 output equals the existing
[44-point report](../offset_grid_44_points.report.toml), including its input.
The other discharged reports are saved here with their exact inputs under
`[input]`. Both modes agreed whenever both completed. All returned positive
witnesses were checked independently using Python integer arithmetic.
For coverable instances, `partitions_evaluated` in the CLI report counts all
enumerated partitions, not just those searched before finding a cover.

## Inputs

All coordinates are exact integers divided by 100,000,000. The 44-point family
consists of `(i*a, j*a)` for `-4 <= i <= 4`, `-2 <= j <= 1`, and
`((i - 1/2)*a, a + h)` for `-3 <= i <= 4`.
For `a = 0.2500, 0.2531, 0.2532`, respectively,
`h = 0.21650635, 0.21919103, 0.21927763`.
These are rational approximations to the equilateral-triangle heights.
The 43-point control instead has a 9-by-3 grid and eight half-offset points
on each outer row; its exact coordinates are in `43_02500_on.report.toml`.

Run the main 44-point instance from the repository root:

```sh
build/release/disjoint_unit_disk_cover \
  --points-path examples/offset_grid_44_points.toml --center-discharge
```

Omit the last flag for the original search. A non-coverable result exits with
code 1, as in the original CLI.

## Why fixing a center is valid

For a candidate rectangle `T_i`, test its four corners for a witness `w` that
covers the assigned points strictly within radius 1. Require also that the
minimum distance from `w` to every other center rectangle is at least 2.
In any feasible choice of centers, replacing center `i` by `w` therefore
preserves coverage and disjointness. Conversely, a solution with center `i`
fixed at `w` was already an allowed solution. Thus feasibility is unchanged.

The implementation replaces `T_i` by the singleton rectangle `{w}`, so it no
longer needs subdivision. Fixed rectangles are restored when backtracking.
Only the choice of corners as witnesses is heuristic; the acceptance test
itself is an exact integer calculation.
