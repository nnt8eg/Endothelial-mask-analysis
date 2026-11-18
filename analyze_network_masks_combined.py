#!/usr/bin/env python3
"""
analyze_network_masks_combined.py

Computes tube-formation-assay style metrics from the binary masks produced by
export_network_mask.groovy, aimed at quantifying NETWORK DYSREGULATION
(irregular vessel diameters/lengths, "chunky" branch points, uneven mesh
spacing) rather than just overall network amount.

This is the COMBINED version: it walks masks/ in *every* dataset directory
listed in PROJECT_DIRS below, computes the same metrics for each image, tags
every row with which dataset it came from (a "dataset" column -- some image
names, e.g. "Cntl 1 1", recur across datasets and are NOT guaranteed to be
the same image), and writes ONE combined workbook/csv.

Every spatial scale is reported in real physical units (um / um^2) only --
there is no separate _px column; only counts (n_meshes, n_branch_points,
etc.) and dimensionless ratios (cv_*, network_area_fraction, lacunarity_)
are left as plain numbers.

Two separate mask versions are used deliberately:
  - RAW mask (exactly as exported): used for anything measuring width/size
    (distance-transform diameter, junction area, mesh area, lacunarity,
    network area). Never smoothed -- smoothing would erase the very
    irregularity being studied.
  - CLEANED skeleton (RAW mask -> morphological opening -> skeletonize ->
    iterative short-spur pruning): used only to decide WHERE the real
    branches/junctions/segments are, since manual trace boundaries are never
    perfectly smooth and produce a spray of tiny spurious skeleton branches.
    Opening radius and prune threshold were chosen by visual comparison
    against the study dataset (see manuscript Methods for details);
    raw/uncleaned topology values are kept alongside (still in um) for
    audit.

Junction counting uses a non-destructive node-merge step on top of that
skeleton graph -- see JUNCTION_MERGE_* below -- to correct for skeletonization
splitting one wide/chunky hub into several closely-spaced graph nodes. This
never touches the skeleton image itself, so total_network_length, n_segments,
and n_branch_points are unaffected; only n_junctions and how junction
diameter/area are sampled change.

A dataset whose masks/ folder doesn't exist yet (export not run in QuPath
for that dataset) is skipped with a warning rather than erroring out.

Requirements: numpy, pandas, tifffile, scipy, scikit-image, skan, openpyxl
(see requirements.txt).

Usage:
    python3 analyze_network_masks_combined.py

Before running, edit PROJECT_DIRS below to point at your own exported
masks/ folders (see README.md for the expected directory layout).
"""

import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from skimage.morphology import skeletonize, opening, disk
from skimage.measure import label, regionprops
from skan import Skeleton, summarize

# ---- configuration -----------------------------------------------------------
# Each entry is (dataset_label, path_to_dataset_folder). That folder must
# contain a masks/ subfolder produced by export_network_mask.groovy.
# Results from every listed dataset are combined into ONE output, written
# next to the first dataset in this list.
#
# EXAMPLE PATHS BELOW -- replace with your own dataset folders.
PROJECT_DIRS = [
    ("dataset_1", Path("./dataset_1")),
    ("dataset_2", Path("./dataset_2")),
]

OUTPUT_XLSX = PROJECT_DIRS[0][1] / "network_metrics.xlsx"
OUTPUT_CSV = PROJECT_DIRS[0][1] / "network_metrics.csv"

# Real pixel calibration (physical size per pixel, in micrometers). Set this
# to match your own imaging setup -- QuPath/Bio-Formats can report a generic
# fallback value (e.g. a 72-dpi default) instead of the true calibrated
# pixel size, so confirm via Image > Properties before trusting this.
# The value below (382 px = 250 um) is the value used for this study and is
# carried over to every dataset in PROJECT_DIRS; update per-dataset if your
# imaging setups differ.
PIXEL_SIZE_UM = 250 / 382

# --- topology-cleanup parameters (do NOT affect diameter/area/lacunarity) ---
OPENING_RADIUS_PX = 3
PRUNE_SPUR_LENGTH_PX = 25
PRUNE_MAX_ITERATIONS = 25

# --- dysregulation-metric parameters -----------------------------------------
JUNCTION_LOCAL_WINDOW_PX = 25
LACUNARITY_BOX_SIZES_PX = [4, 8, 16, 32, 64, 128, 256]

# --- junction over-counting fix ------------------------------------------------
# Manual skeletonization of a "chunky"/wide hub can split into several
# closely-spaced degree>=3 graph nodes instead of one. These are merged for
# COUNTING/MEASUREMENT only (via union-find on the graph nodes) -- the
# skeleton image itself, and therefore total_network_length/n_segments/
# n_branch_points, are never touched.
#
# Automated rule: merge two junction nodes connected by a junction-junction
# (type-2) branch if they're closer than JUNCTION_MERGE_K x the smaller
# node's local vessel radius, AND closer than a hard absolute cap (independent
# of vessel width). The absolute cap is required: without it, a large enough K
# to catch tight chunky hubs also falsely merges real, distinct branch points
# spread along a wide vessel trunk (confirmed on a 140px-apart, 3-node case).
#
# This automated rule is intentionally conservative (favors under-merging
# over false merges). A few cases visually confirmed as "one hub" don't meet
# it numerically -- distance and vessel-width proxies could not be made to
# separate them from real wide-trunk trios of distinct branch points, because
# the numbers overlap between the two cases even though the shapes look very
# different by eye (one is a blobby, multi-fingered mass; the other is a
# clean, uniformly-wide corridor with real perpendicular branches). For those,
# specific coordinate pairs confirmed by visual review are listed in
# JUNCTION_MERGE_OVERRIDES below and force-merged regardless of distance.
JUNCTION_MERGE_K = 2.0
JUNCTION_MERGE_ABS_CAP_PX = 40

# Force-merge specific junction pairs confirmed by visual review, keyed by a
# substring of the image base name. Each entry is a list of
# ((r1, c1), (r2, c2)) point pairs in RAW PIXEL coordinates of that image; the
# junction graph node nearest each point (within
# JUNCTION_MERGE_OVERRIDE_SEARCH_PX) is merged. A warning is printed (not an
# error) if no matching node is found within range.
#
# The entry below is specific to this study's dataset; remove it (or replace
# with your own confirmed cases) when analyzing different images.
JUNCTION_MERGE_OVERRIDE_SEARCH_PX = 8
JUNCTION_MERGE_OVERRIDES = {
    # "TGF10 1 5": visually confirmed to be one amorphous, multi-fingered
    # blob (not two separate hubs); no distance/radius/shape rule tested
    # could separate this case from real wide-trunk branch trios without
    # reintroducing false merges elsewhere.
    "TGF10 1 5": [((842, 684), (888, 732))],
}
# -------------------------------------------------------------------------------


class UnionFind:
    def __init__(self, items):
        self.parent = {i: i for i in items}

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


# Order matters: more specific patterns (PEGTGF10, PEGTGF5) must precede the
# substrings they contain (TGF10, TGF5). "PFM" is the condition name; it's
# written "Fc"/"FC" in filenames.
#
# These patterns are specific to this study's filename conventions -- edit
# to match your own treatment-group naming scheme.
GROUP_PATTERNS = [
    ("PEGTGF10", re.compile(r"PEGTGF10", re.IGNORECASE)),
    ("PEGTGF5", re.compile(r"PEGTGF5", re.IGNORECASE)),
    ("TGF10", re.compile(r"TGF10", re.IGNORECASE)),
    ("TGF5", re.compile(r"TGF5", re.IGNORECASE)),
    ("PFM", re.compile(r"(?<![A-Za-z])FC(?![A-Za-z])", re.IGNORECASE)),
    ("Cntl", re.compile(r"Cntl", re.IGNORECASE)),
]


# Explicit per-image overrides for filenames that don't carry a recognizable
# condition string. Checked before GROUP_PATTERNS. Confirmed via visual
# review: "Image018" is PFM in this study's dataset.
GROUP_OVERRIDES = {
    "Image018": "PFM",
}


def parse_group(name: str) -> str:
    for needle, group_label in GROUP_OVERRIDES.items():
        if needle.lower() in name.lower():
            return group_label
    for group_label, pattern in GROUP_PATTERNS:
        if pattern.search(name):
            return group_label
    return "Unknown"


def count_branch_points(skeleton: np.ndarray) -> int:
    neighbor_count = ndi.convolve(
        skeleton.astype(np.uint8),
        np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]]),
        mode="constant",
        cval=0,
    )
    return int(np.sum((neighbor_count >= 3) & skeleton))


def prune_spurs(skel_obj: Skeleton, min_length_px: float, max_iterations: int):
    n_removed_total = 0
    branch_data = summarize(skel_obj, separator="-")
    for _ in range(max_iterations):
        spur_idx = branch_data.index[
            (branch_data["branch-type"] == 1) & (branch_data["branch-distance"] < min_length_px)
        ].tolist()
        if not spur_idx:
            break
        skel_obj = skel_obj.prune_paths(spur_idx)
        n_removed_total += len(spur_idx)
        branch_data = summarize(skel_obj, separator="-")
    return skel_obj, branch_data, n_removed_total


def dist_stats(values: np.ndarray, prefix: str) -> dict:
    if len(values) == 0:
        return {
            f"mean_{prefix}": 0.0, f"median_{prefix}": 0.0,
            f"std_{prefix}": 0.0, f"cv_{prefix}": 0.0,
        }
    mean = float(np.mean(values))
    std = float(np.std(values))
    return {
        f"mean_{prefix}": mean,
        f"median_{prefix}": float(np.median(values)),
        f"std_{prefix}": std,
        f"cv_{prefix}": (std / mean) if mean else 0.0,
    }


def gliding_box_lacunarity(mask: np.ndarray, box_size: int) -> float:
    integral = np.cumsum(np.cumsum(mask.astype(np.float64), axis=0), axis=1)
    integral = np.pad(integral, ((1, 0), (1, 0)), mode="constant")
    br = integral[box_size:, box_size:]
    tl = integral[:-box_size, :-box_size]
    tr = integral[:-box_size, box_size:]
    bl = integral[box_size:, :-box_size]
    box_sums = br - tr - bl + tl
    masses = box_sums.ravel()
    mean_m = masses.mean()
    if mean_m == 0:
        return float("nan")
    var_m = masses.var()
    return float(var_m / (mean_m ** 2) + 1)


def px_to_um(row: dict, px: float) -> dict:
    """Convert every pixel-scaled key in `row` to um, replacing (not
    duplicating) the pixel column. cv_* keys are dimensionless (invariant
    under uniform scaling) so are only renamed, never multiplied. Counts and
    other non-scale keys pass through unchanged."""
    out = {}
    for key, value in row.items():
        is_cv = key.startswith("cv_")
        if "_px2" in key:
            new_key = key.replace("_px2", "_um2")
            out[new_key] = value if is_cv else value * px * px
        elif "_px" in key:
            new_key = key.replace("_px", "_um")
            out[new_key] = value if is_cv else value * px
        else:
            out[key] = value
    return out


def find_nearest_junction(coords, junction_ids, target_r, target_c):
    best_id, best_dist2 = None, None
    for nid in junction_ids:
        r, c = coords[nid]
        d2 = (r - target_r) ** 2 + (c - target_c) ** 2
        if best_dist2 is None or d2 < best_dist2:
            best_id, best_dist2 = nid, d2
    if best_id is None:
        return None, None
    return best_id, best_dist2 ** 0.5


def analyze_image(network_path: Path, negative_path: Path, dataset_label: str) -> dict:
    base_name = network_path.name.replace("_network_mask.tif", "")

    raw_network = tifffile.imread(str(network_path)) > 0
    negative = tifffile.imread(str(negative_path)) > 0
    height, width = raw_network.shape

    row = {
        "dataset": dataset_label,
        "image": base_name,
        "group": parse_group(base_name),
        "width_px": width,
        "height_px": height,
    }

    # --- raw (unpruned, unopened) skeleton, for audit comparison ---
    total_length_px_unpruned = 0.0
    n_branch_points_unpruned = 0
    n_segments_unpruned = 0
    raw_skel = skeletonize(raw_network)
    if raw_skel.sum() > 0:
        raw_skel_obj = Skeleton(raw_skel)
        raw_branch_data = summarize(raw_skel_obj, separator="-")
        total_length_px_unpruned = float(raw_branch_data["branch-distance"].sum())
        n_segments_unpruned = len(raw_branch_data)
        n_branch_points_unpruned = count_branch_points(raw_skel)

    # --- cleaned topology: opening -> skeletonize -> prune short spurs ---
    cleaned_mask = opening(raw_network, disk(OPENING_RADIUS_PX)) if OPENING_RADIUS_PX > 0 else raw_network
    cleaned_skel = skeletonize(cleaned_mask)

    total_length_px = 0.0
    n_branch_points = 0
    n_segments = 0
    n_spurs_pruned = 0
    dt = ndi.distance_transform_edt(raw_network)

    if cleaned_skel.sum() > 0:
        skel_obj = Skeleton(cleaned_skel)
        pruned_skel_obj, branch_data, n_spurs_pruned = prune_spurs(
            skel_obj, PRUNE_SPUR_LENGTH_PX, PRUNE_MAX_ITERATIONS
        )
        final_skel = pruned_skel_obj.skeleton_image > 0
        total_length_px = float(branch_data["branch-distance"].sum())
        n_segments = len(branch_data)
        n_branch_points = count_branch_points(final_skel)

        diam_samples = 2 * dt[final_skel]
        row.update(dist_stats(diam_samples, "vessel_diameter_px"))
        row.update(dist_stats(branch_data["branch-distance"].values, "segment_length_px"))

        degrees = pruned_skel_obj.degrees
        coords = pruned_skel_obj.coordinates
        node_ids = set(branch_data["node-id-src"]).union(branch_data["node-id-dst"])
        junction_ids = set(int(n) for n in node_ids if degrees[int(n)] >= 3)

        def local_radius(nid):
            r, c = coords[int(nid)].astype(int)
            return dt[r, c]

        # --- automated merge: radius-scaled distance + hard absolute cap ---
        uf = UnionFind(junction_ids)
        for _, br in branch_data.iterrows():
            if br["branch-type"] != 2:
                continue
            src, dst = int(br["node-id-src"]), int(br["node-id-dst"])
            if src not in junction_ids or dst not in junction_ids:
                continue
            d = br["branch-distance"]
            thresh = JUNCTION_MERGE_K * min(local_radius(src), local_radius(dst))
            if d < thresh and d < JUNCTION_MERGE_ABS_CAP_PX:
                uf.union(src, dst)

        # --- manual overrides for visually-confirmed cases the rule misses ---
        for needle, pairs in JUNCTION_MERGE_OVERRIDES.items():
            if needle.lower() not in base_name.lower():
                continue
            for (pr, pc), (qr, qc) in pairs:
                n1, dist1 = find_nearest_junction(coords, junction_ids, pr, pc)
                n2, dist2 = find_nearest_junction(coords, junction_ids, qr, qc)
                if (
                    n1 is not None and n2 is not None
                    and dist1 <= JUNCTION_MERGE_OVERRIDE_SEARCH_PX
                    and dist2 <= JUNCTION_MERGE_OVERRIDE_SEARCH_PX
                ):
                    uf.union(n1, n2)
                else:
                    print(
                        f"WARNING [{base_name}]: junction merge override near "
                        f"({pr},{pc})/({qr},{qc}) did not find matching nodes within "
                        f"{JUNCTION_MERGE_OVERRIDE_SEARCH_PX}px -- skipped"
                    )

        clusters = {}
        for nid in junction_ids:
            root = uf.find(nid)
            clusters.setdefault(root, []).append(nid)

        w = JUNCTION_LOCAL_WINDOW_PX
        junction_diams = []
        junction_areas = []
        for members in clusters.values():
            pts = np.array([coords[m] for m in members])
            diam_vals = [2 * dt[int(r), int(c)] for r, c in pts]
            junction_diams.append(float(np.mean(diam_vals)))
            r0 = int(max(0, pts[:, 0].min() - w))
            r1 = int(min(height, pts[:, 0].max() + w + 1))
            c0 = int(max(0, pts[:, 1].min() - w))
            c1 = int(min(width, pts[:, 1].max() + w + 1))
            junction_areas.append(int(raw_network[r0:r1, c0:c1].sum()))

        row["n_junctions"] = len(clusters)
        row["n_junction_merges_applied"] = sum(1 for m in clusters.values() if len(m) > 1)
        row.update(dist_stats(np.array(junction_diams), "junction_diameter_px"))
        row.update(dist_stats(np.array(junction_areas), "junction_local_area_px2"))
    else:
        row["n_junctions"] = 0
        row["n_junction_merges_applied"] = 0
        row.update(dist_stats(np.array([]), "vessel_diameter_px"))
        row.update(dist_stats(np.array([]), "segment_length_px"))
        row.update(dist_stats(np.array([]), "junction_diameter_px"))
        row.update(dist_stats(np.array([]), "junction_local_area_px2"))

    row["total_network_length_px"] = total_length_px
    row["n_branch_points"] = n_branch_points
    row["n_segments"] = n_segments
    row["n_spurs_pruned"] = n_spurs_pruned
    row["total_network_length_px_unpruned"] = total_length_px_unpruned
    row["n_branch_points_unpruned"] = n_branch_points_unpruned
    row["n_segments_unpruned"] = n_segments_unpruned

    # --- mesh metrics (from the traced negative-space regions) ---
    labeled = label(negative, connectivity=2)
    props = regionprops(labeled)
    border_labels = (
        set(np.unique(labeled[0, :]))
        | set(np.unique(labeled[-1, :]))
        | set(np.unique(labeled[:, 0]))
        | set(np.unique(labeled[:, -1]))
    )
    border_labels.discard(0)
    enclosed_areas = np.array([p.area for p in props if p.label not in border_labels], dtype=float)
    all_areas = [p.area for p in props]

    row["n_meshes"] = len(enclosed_areas)
    row["n_meshes_incl_border_touching"] = len(all_areas)
    row.update(dist_stats(enclosed_areas, "mesh_area_px2"))

    # --- network area / area fraction (raw mask); note: area, so _px2 ---
    network_area_px = int(raw_network.sum())
    row["network_area_px2"] = network_area_px
    row["network_area_fraction"] = network_area_px / (height * width)

    row = px_to_um(row, PIXEL_SIZE_UM)

    # --- lacunarity (raw mask, gliding box, multiple scales) ---
    # Box sizes are a spatial scale too -> label columns by their um size.
    lac_values = []
    for box_px in LACUNARITY_BOX_SIZES_PX:
        if box_px >= min(height, width):
            continue
        lac = gliding_box_lacunarity(raw_network, box_px)
        box_um = box_px * PIXEL_SIZE_UM
        row[f"lacunarity_r{box_um:.1f}um"] = lac
        lac_values.append(lac)
    row["lacunarity_mean_across_scales"] = float(np.mean(lac_values)) if lac_values else float("nan")

    return row


def main():
    if not PIXEL_SIZE_UM:
        raise SystemExit("PIXEL_SIZE_UM must be set to produce um-only output.")

    rows = []
    for dataset_label, project_dir in PROJECT_DIRS:
        masks_dir = project_dir / "masks"
        if not masks_dir.is_dir():
            print(f"SKIP dataset '{dataset_label}': no masks/ folder at {masks_dir} "
                  f"(run export_network_mask.groovy for this dataset first)")
            continue
        network_files = sorted(masks_dir.glob("*_network_mask.tif"))
        if not network_files:
            print(f"SKIP dataset '{dataset_label}': no *_network_mask.tif files found in {masks_dir}")
            continue
        for network_path in network_files:
            base_name = network_path.name.replace("_network_mask.tif", "")
            negative_path = masks_dir / f"{base_name}_negative_mask.tif"
            if not negative_path.exists():
                print(f"WARNING [{dataset_label}]: missing negative mask for {base_name}, skipping")
                continue
            try:
                row = analyze_image(network_path, negative_path, dataset_label)
                rows.append(row)
                print(
                    f"OK [{dataset_label}] {base_name}: length={row['total_network_length_um']:.0f}um, "
                    f"branch pts={row['n_branch_points']} ({row['n_spurs_pruned']} spurs pruned), "
                    f"junctions={row['n_junctions']} ({row['n_junction_merges_applied']} merged), "
                    f"meshes={row['n_meshes']}, "
                    f"vessel diam mean={row['mean_vessel_diameter_um']:.1f}um (cv={row['cv_vessel_diameter_um']:.2f}), "
                    f"junction diam mean={row['mean_junction_diameter_um']:.1f}um (cv={row['cv_junction_diameter_um']:.2f})"
                )
            except Exception as e:
                print(f"ERROR [{dataset_label}] {base_name}: {e}")

    if not rows:
        print("No rows produced from any dataset -- nothing written.")
        return

    df = pd.DataFrame(rows)
    df.to_csv(OUTPUT_CSV, index=False)
    try:
        df.to_excel(OUTPUT_XLSX, index=False)
    except Exception as e:
        print(f"(xlsx export skipped: {e})")

    print(f"\nWrote {len(df)} rows to {OUTPUT_XLSX}")
    print(df["dataset"].value_counts().to_string())


if __name__ == "__main__":
    main()
