# Network dysregulation analysis pipeline

This pipeline quantifies **endothelial-fibroblast networks** under TGF-B1 stimulation (irregular vessel diameters, irregular segment lengths, "chunky" branch points, uneven mesh
spacing) from manually-traced negative-space annotations of
endothelial/fibroblast networks in QuPath.

## Pipeline

1. **`export_network_mask.groovy`** — run inside QuPath (`Automate > Show
   script editor > Run > Run for project`) to convert manually-traced
   negative-space annotations into two binary TIFF masks per image:
   - `<image>_negative_mask.tif` — 255 = traced empty space, 0 = elsewhere
   - `<image>_network_mask.tif` — 255 = network, 0 = traced empty space

   Masks are written to a `masks/` subfolder inside the QuPath project
   folder.

2. **`analyze_network_masks_combined.py`** — walks `masks/` in one or more
   dataset folders, computes the dysregulation metrics described below for
   every image, and writes a combined `network_metrics.xlsx` / `.csv`.

## Expected directory layout

```
your_repo/
├── export_network_mask.groovy
├── analyze_network_masks_combined.py
├── requirements.txt
├── dataset_1/
│   └── masks/
│       ├── image_a_negative_mask.tif
│       ├── image_a_network_mask.tif
│       └── ...
└── dataset_2/
    └── masks/
        └── ...
```

Edit `PROJECT_DIRS` at the top of `analyze_network_masks_combined.py` to
point at your own dataset folder(s); a single dataset works the same way
with a one-entry list.

## Setup

```
pip install -r requirements.txt
python3 analyze_network_masks_combined.py
```

## Metrics computed per image

- **Mesh metrics**: number of enclosed meshes, mean/median/std/CV of mesh
  area — from connected components of the negative-space mask.
- **Topology** (from a morphologically opened + skeletonized + pruned
  version of the network mask): total network length, number of branch
  points, number of segments, number of pruned spurs (plus unpruned audit
  versions of each).
- **Vessel diameter**: mean/median/std/CV of local vessel width (2× the
  distance transform of the raw mask), sampled along the cleaned skeleton.
- **Segment length**: mean/median/std/CV of individual branch lengths.
- **Junction metrics**: number of junctions (after a non-destructive
  node-merge correction for skeletonization splitting one wide/chunky hub
  into several nearby graph nodes — see comments in the script),
  mean/median/std/CV of junction diameter and of junction local area.
- **Network area and area fraction** (raw mask).
- **Lacunarity** at multiple box sizes (gliding-box method, raw mask), plus
  the mean across scales — captures spatial heterogeneity of the network
  pattern independently of the size distributions above.

All spatial-scale outputs are reported in physical units (µm / µm²) only,
using the pixel calibration set via `PIXEL_SIZE_UM`. Coefficient-of-variation
(`cv_*`) columns are dimensionless and reported as-is.

Two mask versions are used deliberately throughout: the **raw** mask for
anything measuring width/size (it is never smoothed, since smoothing would
erase the irregularity being studied), and a **cleaned** skeleton (opening →
skeletonize → short-spur pruning) only to determine *where* branches and
junctions are, since manually-traced boundaries produce small spurious
skeleton branches that would otherwise inflate branch/segment counts.

## Notes on study-specific parameters

`GROUP_PATTERNS`, `GROUP_OVERRIDES`, and `JUNCTION_MERGE_OVERRIDES` in
`analyze_network_masks_combined.py` encode filename conventions and specific,
visually-confirmed cases from this study's dataset (treatment-group naming,
one image with no parseable condition string, and one junction pair
confirmed by visual review to be a single fused hub rather than two). These
are left in place for exact reproducibility of the published results; adapt
or remove them when applying this pipeline to a different dataset.

## License

Add a license file appropriate to your institution's and journal's policy
before publishing this repository.
