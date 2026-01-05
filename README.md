# disjoint-unit-disk-cover

This repository contains a self-contained research implementation for analyzing whether a planar point set can be covered by pairwise-disjoint unit-radius disks.

## Repository layout

- `src/` – implementation (core + CLI entry point)
- `include/disjoint_unit_disk/` – internal headers used by this project
- `tests/` – GoogleTest-based unit tests
- `examples/` – TOML point sets used for regression / experiments
  - `offset_grid_44_points.toml` – 44 points on an offset grid; NonCoverable
  - `concentric_points_coverable.toml` – 45 points on concentric circles; Coverable
  - `concentric_points_noncoverable.toml` – 45 points on concentric circles; NonCoverable
  - `concentric_points_coverable_with_groups.toml` – same as above, with pre-assigned groups
  - `reports/` – precomputed output reports for the above
- `tools/` – scripts to generate example inputs
- `third_party/` – git submodules (`CLI11`, `toml++`, `googletest`, `cgal`)

## Build & test

```bash
git submodule update --init --recursive

cmake -S . -B build/debug -G Ninja
cmake --build build/debug
ctest --test-dir build/debug

cmake -S . -B build/release -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build/release
ctest --test-dir build/release
```

CGAL requires system packages (Boost, GMP, MPFR, Eigen) available to CMake.

## Run the CLI

```bash
./build/release/disjoint_unit_disk_cover --points-path examples/offset_grid_44_points.toml

# Enable verbose logging
./build/release/disjoint_unit_disk_cover --points-path examples/offset_grid_44_points.toml --verbose

# Write report to TOML
./build/release/disjoint_unit_disk_cover --points-path examples/offset_grid_44_points.toml --report-path examples/reports/offset_grid_44_points.report.toml
```

`--verbose` enables info-level progress logs.

## Input file format (TOML)

A points file contains:

- `scale` (integer, >= 1)
- `[[points]]` entries with integer `x`, `y`
- optional `group` (non-negative integer) for all points (or for none)

Interpretation:

- The solver treats the real coordinates as `(x/scale, y/scale)`.
- Disks have radius `1` in this normalized coordinate system.
- Internally, computations are performed on the integer coordinates with a disk radius of `scale`.

Example:

```toml
scale = 2

[[points]]
x = 0
y = 1

[[points]]
x = 2
y = 3
```

If `group` is provided, points are analyzed with those fixed groups (one disk per group).

## Output status

The CLI reports one of:

- `Coverable` – a disjoint cover was found; sample centers are printed.
- `NonCoverable` – no disjoint cover was found after ruling out all candidates explored by the program.
- `Undecided` – the program could not decide (e.g., when only non-lattice centers are feasible under the current integer subdivision).

## Regenerating example TOMLs

`tools/generate_examples.py` can regenerate some example inputs:

```bash
python tools/generate_examples.py --target concentric-coverable
python tools/generate_examples.py --target concentric-noncoverable
python tools/generate_examples.py --target offset-grid-44-points
```

To regenerate with Arb-based rounding certification, install `python-flint`:

```bash
pip install -r tools/requirements.txt
```

## License

This repository is licensed under GPL-3.0-only because the implementation depends on GPL-licensed CGAL components.
Point configuration data and generated reports are licensed separately (see [DATA_LICENSE.md](./DATA_LICENSE.md)).

| Component               | License                                                            |
| ----------------------- | ------------------------------------------------------------------ |
| Source code             | GPL-3.0-only — See [LICENSE](./LICENSE)                            |
| 44-point data & reports | CC BY 4.0 — See [DATA_LICENSE.md](./DATA_LICENSE.md)               |
| 45-point data & reports | Derived from prior work — See [DATA_LICENSE.md](./DATA_LICENSE.md) |
| Third-party libraries   | See each library's LICENSE in `third_party/`                       |
