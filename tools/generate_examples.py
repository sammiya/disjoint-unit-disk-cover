#!/usr/bin/env python3
"""
Generate example point sets as TOML files.

For each exact coordinate (x_true, y_true), we compute

    sx_true = SCALE * x_true
    sy_true = SCALE * y_true

round them to nearest integers (ix, iy), and use Arb to prove that

    |sx_true - ix| < 1/2  and  |sy_true - iy| < 1/2

for all points.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List, Tuple

from flint import arb, fmpq, ctx

# ----------------------------------------------------------------------
# Global settings
# ----------------------------------------------------------------------

# Arb working precision (in bits)
INITIAL_PREC = 256

# Set initial Arb precision
ctx.prec = INITIAL_PREC

# Fixed integer scale factor
SCALE = 100_000_000  # 1e8

# Absolute error bound in the scaled coordinate system:
# We want to prove |scaled_true - rounded_integer| < 1/2.
ERROR_BOUND = arb(1) / 2


# ----------------------------------------------------------------------
# I/O
# ----------------------------------------------------------------------

def write_points_toml(points_int: Iterable[Tuple[int, int]], output: Path) -> None:
    """
    Write a TOML file with scaled integer coordinates and metadata.

    Each point is represented as integers (ix, iy) such that

        ix ≈ SCALE * x_true
        iy ≈ SCALE * y_true

    and Arb has certified that the rounding error is < 1/2 in absolute value.
    """
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as f:
        f.write(f"scale = {SCALE}\n")
        f.write("\n")
        for ix, iy in points_int:
            f.write("[[points]]\n")
            f.write(f"x = {ix}\n")
            f.write(f"y = {iy}\n\n")


# ----------------------------------------------------------------------
# Point set construction (exact, unscaled, in Arb)
# ----------------------------------------------------------------------

def build_concentric_coverable() -> List[Tuple[arb, arb]]:
    radii = [fmpq(1, 10), fmpq(72, 100), fmpq(10001, 10000)]
    return _build_concentric(radii, [3, 21, 21])


def build_concentric_noncoverable() -> List[Tuple[arb, arb]]:
    radii = [fmpq(1, 10), fmpq(721, 1000), fmpq(10001, 10000)]
    return _build_concentric(radii, [3, 21, 21])


def _build_concentric(radii, counts) -> List[Tuple[arb, arb]]:
    pts: List[Tuple[arb, arb]] = []
    for radius, count in zip(radii, counts):
        r = arb(radius)
        for k in range(count):
            theta = arb(2) * arb.pi() * arb(k) / arb(count)
            x = r * theta.cos()
            y = r * theta.sin()
            pts.append((x, y))
    return pts


def build_offset_grid_44_points() -> List[Tuple[arb, arb]]:
    a = fmpq(2531, 10000)
    pts: List[Tuple[arb, arb]] = []

    for n in range(-4, 5):  # -4..4
        for m in range(-2, 2):  # -2..1
            x = arb(a * n)
            y = arb(a * m)
            pts.append((x, y))

    y_top = arb(a) + (arb(3).sqrt() / arb(2)) * arb(a)
    for n in range(-3, 5):  # -3..4 (8 points)
        x = arb((fmpq(n) - fmpq(1, 2)) * a)
        pts.append((x, y_top))

    return pts


def build_points(target: str) -> List[Tuple[arb, arb]]:
    """Build unscaled exact points (in Arb) for the chosen preset."""
    if target == "concentric-coverable":
        return build_concentric_coverable()
    if target == "concentric-noncoverable":
        return build_concentric_noncoverable()
    if target == "offset-grid-44-points":
        return build_offset_grid_44_points()
    raise ValueError(f"Unknown target {target}")


# ----------------------------------------------------------------------
# Rounding and error checking
# ----------------------------------------------------------------------

def choose_integer_points(points_exact_scaled: Iterable[Tuple[arb, arb]]) -> List[Tuple[int, int]]:
    """
    For each scaled exact point (sx, sy) in Arb, choose an integer pair (ix, iy)
    by rounding the midpoint to the nearest integer.

    The correctness of this choice is *not* assumed; it is verified later
    by Arb via check_error_integer().
    """
    points_int: List[Tuple[int, int]] = []
    for sx, sy in points_exact_scaled:
        # Python float(sx) uses the Arb midpoint; round() gives an integer
        ix = int(round(float(sx)))
        iy = int(round(float(sy)))
        points_int.append((ix, iy))
    return points_int


def check_error_integer(points_int: Iterable[Tuple[int, int]],
                        points_exact_scaled: Iterable[Tuple[arb, arb]]):
    """
    For each coordinate, check with Arb that

        |scaled_true - rounded_integer| < ERROR_BOUND (= 1/2).

    The comparison uses Arb's ball arithmetic, so a successful check is
    a formal proof that the inequality holds for the exact real values.

    Returns:
        ok:        True if all coordinates satisfy the bound.
        max_err:  Approximate maximum error (midpoint of the largest ball).
    """
    ok = True
    max_err_mid = 0.0  # for logging only

    for (ix, iy), (sx, sy) in zip(points_int, points_exact_scaled):
        for n, e in ((ix, sx), (iy, sy)):
            # Arb ball for the absolute error |e_true - n|
            err = abs(e - arb(n))
            # Track the midpoint of the largest error ball, just for diagnostics
            max_err_mid = max(max_err_mid, float(err))
            # Arb's `<` is strict: err < ERROR_BOUND means sup(err) < ERROR_BOUND.
            if not (err < ERROR_BOUND):
                ok = False
    return ok, max_err_mid


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate example TOMLs with scaled integer coordinates and check "
            "|SCALE * coord_true - rounded_int| < 1/2 for all coordinates."
        )
    )
    parser.add_argument(
        "--target",
        choices=[
            "concentric-coverable",
            "concentric-noncoverable",
            "offset-grid-44-points",
        ],
        required=True,
        help="Which preset to generate.",
    )
    parser.add_argument(
        "--output",
        "--points",
        dest="output",
        type=Path,
        default=None,
        help="Output TOML path (defaults depend on target).",
    )
    return parser.parse_args()


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------

def main() -> None:
    args = parse_args()
    default_outputs = {
        "concentric-coverable": Path("examples/concentric_points_coverable.toml"),
        "concentric-noncoverable": Path("examples/concentric_points_noncoverable.toml"),
        "offset-grid-44-points": Path("examples/offset_grid_44_points.toml"),
    }
    output = args.output or default_outputs[args.target]

    prec = INITIAL_PREC
    ctx.prec = prec

    while True:
        print(f"Computing points at precision {prec} bits...")

        # 1. Exact points in Arb (unscaled)
        points_exact = build_points(args.target)

        # 2. Scale all coordinates by SCALE in Arb
        scale_arb = arb(SCALE)
        points_exact_scaled: List[Tuple[arb, arb]] = [
            (scale_arb * px, scale_arb * py) for px, py in points_exact
        ]

        # 3. Choose integer approximations by rounding the midpoint
        points_int: List[Tuple[int, int]] = choose_integer_points(points_exact_scaled)

        # 4. Verify |scaled_true - rounded_int| < 1/2 for all coordinates
        ok, max_err = check_error_integer(points_int, points_exact_scaled)

        print(f"  Approx max |int - Arb_scaled| ≈ {max_err:.3e}")
        print(f"  Error bound checked          : {float(ERROR_BOUND):.3e}")

        if ok:
            print("Check passed at this precision.")
            break

        # If the bound was not proven, increase precision and try again.

        prec *= 2
        ctx.prec = prec
        print("  Error bound not yet proven; increasing precision.\n")

    print("Final result: |int - true_scaled| < 1/2 for all coordinates.")
    write_points_toml(points_int, output)
    print(f"Wrote {len(points_int)} integer points to {output}")


if __name__ == "__main__":
    main()
