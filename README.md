# scverse Cell Painting tutorials

Tutorials for **downstream analysis of Cell Painting data with scverse tools** —
`anndata`, `scanpy`, `decoupler`, `pertpy`.

Started at the scverse Cell Painting hackathon. The goal is to take the analysis patterns
that are well established for single-cell transcriptomics — as laid out in
[single-cell best practices](https://www.sc-best-practices.org/) — and work out which of
them transfer to image-based morphological profiling, which need adapting, and which do
not transfer at all. Each notebook states the assumption a step makes, whether it survives
the move to morphological profiles, and what to do instead when it doesn't.

## The idea in one substitution

**A well is an observation, and a CellProfiler feature is a gene.** Where that holds, the
scverse stack works on Cell Painting data unchanged. Where it breaks is where the
interesting work is.

| single-cell concept | Cell Painting equivalent | Transfers? |
| --- | --- | --- |
| Cell | Well (well-level) or segmented cell | yes |
| Gene | Morphological feature (compartment x family x channel) | yes, but features are correlated by construction |
| Count matrix | Continuous feature matrix | **no** — no counts, no sparsity, no overdispersion |
| Size-factor normalisation + `log1p` | Per-plate MAD-robustize against DMSO controls | **no** — controls replace library size |
| Highly variable gene selection | Redundancy pruning (`pycytominer.feature_select`) | replaced — no mean-variance trend to fit |
| Batch (donor, 10x run) | Plate, well position, and imaging site | yes, and **well position is its own confound** |
| PCA -> kNN -> UMAP -> Leiden | identical | yes, after clipping heavy tails |
| Cell-type annotation | Mechanism of action / target annotation | yes, but labels are sparse and noisy |
| Marker genes per cluster | Discriminative features per cluster | yes, with a rank-based test |
| Pathway / GSEA over gene sets | Enrichment over **feature families** from `var` | yes, once `var` is annotated |
| Perturbation modelling (`pertpy`) | Effect size, replicate reproducibility | yes — this is the native question |
| Spatial coordinates | Well row/column on the plate | yes, and underused |

Two of these are what the scverse stack adds over existing Cell Painting tooling:

- **Feature families are gene sets.** Once `var` carries `compartment`, `family` and
  `channel`, thousands of opaque features collapse into ~14 interpretable morphological
  programmes that `decoupler` can score. Notebook 02.
- **Plate layout is a spatial graph.** Well row and column are real 2D coordinates, so
  edge effects are *spatial autocorrelation* — which makes `squidpy.gr.spatial_autocorr` a
  principled plate-artifact detector. Not yet covered; see the roadmap.

## Dataset

The tutorials run on the **EU-OPENSCREEN bioactives** Cell Painting screen
([Wolff et al. 2024](https://doi.org/10.1101/2024.08.27.609964)), `IMTM_HepG2` subset,
converted to AnnData by
[`scverse/cell-painting-io`](https://github.com/scverse/cell-painting-io/blob/cpjump1-anndata-notebook/eu_os_bioactives_to_anndata.ipynb).

**10,668 wells x 2,776 features.** 7 library plates (`B1001`-`B1007`) x 4 replicate plates
(`R1`-`R4`) x 384 wells, HepG2 cells at the IMTM site.

| | |
| --- | --- |
| `.X` | per-plate robust z-score vs DMSO (`mad_robustize`) |
| `layers["aggregated"]` | raw per-well median profiles |
| `obs` | `Plate`, `Replicate`, `Well`, `cell_count`, `EOS` (2,459 compounds), `pert_type`, `compound`, `smiles`, `inchikey`, `mw`, `logp`, `approved_drug`, `tubulin_binder`, `target_genes` |
| `var` | `compartment` (Nuclei/Cells/Cytoplasm), `family` (8), `channel` (DNA/ER/AGP/Mito), `channel_2`, `n_channels` |
| controls | 784 DMSO wells; 112 positive-control wells (nocodazole, tetrandrine) |
| ground truth | 20 annotated tubulin binders |

Note this is a **four-channel** stain — there is no RNA channel — and a **single-dose**
screen at 10 uM, so `concentration_uM` is collinear with `pert_type` and carries no extra
information.

Set `EU_OS_H5AD` if your copy is not at `data/eu_os_imtm_hepg2.h5ad`.

## Tutorials

| Notebook | Chapter it mirrors | Runtime |
| --- | --- | --- |
| [`01_dimensionality_reduction_and_clustering.ipynb`](notebooks/01_dimensionality_reduction_and_clustering.ipynb) | Dimensionality reduction, Clustering | ~30 s |
| [`02_feature_family_enrichment_decoupler.ipynb`](notebooks/02_feature_family_enrichment_decoupler.ipynb) | Pathway / gene-set analysis | ~15 s |
| [`03_perturbation_analysis_pertpy.ipynb`](notebooks/03_perturbation_analysis_pertpy.ipynb) | Perturbation modelling | ~40 s |

These are **backbone notebooks**: every cell runs against the dataset above, and points
where you should make a judgement call are marked **TODO**. They are committed without
outputs — run them to fill in.

The data-preparation notebooks under `notebooks/` (adapted from
[timtreis/2024_broad_hackathon](https://github.com/timtreis/2024_broad_hackathon)) build a
JUMP TARGET2 AnnData from raw profiles; see [`notebooks/README.md`](notebooks/README.md).

## What the tutorials establish

Findings that came out of writing them, each reproducible from the notebooks:

- **Clip the heavy tails or nothing works.** `mad_robustize` divides by the DMSO MAD, so
  near-constant features explode — `.X` spans -2,115 to +42,636. Unclipped, PC1 absorbs
  **93% of the variance** and is driven by `ObjectSkeleton`/`Neighbors` division
  artifacts. Clipping at +/-10 drops it to 31% with sensible `Texture` loadings.
- **Cell count drives PC1** (r = 0.57). The dominant axis of the screen is substantially
  confluence and cytotoxicity, so every cluster needs checking against it before being
  called a mechanism.
- **Well position is a first-class confound.** Per-plate normalisation removes the plate
  offset but not a gradient recurring at the same coordinates on every plate.
- **State your reproducibility null.** Percent replicating is **81%** against a random
  null but **64%** when the null controls for well position — a 17-point swing on a choice
  that often goes unstated. Same-compound wells reach +0.61 mean similarity against +0.12
  for wells merely sharing a plate coordinate, so the signal is real; the null still
  changes the headline.
- **`ulm` equals `zscore` on unweighted sets.** With binary membership and no covariates
  the two correlate at r = 1.00; `ulm` only earns its keep once sets carry weights.
- **Permutation count silently caps significance.** `pt.tl.DistanceTest` cannot return a
  p-value below ~`1/n_perms`, so at `n_perms=100` across 23 groups Holm-Šidák leaves
  **zero** significant results (smallest adjusted p = 0.21) however large the effects.
  At `n_perms=1000`, 22 of 23 are significant. Keep `n_perms` well above `n_groups / alpha`.
- **Feature families recover known biology.** Tubulin binders shift `channel:DNA` by
  -4.2 z and `channel:AGP` by +2.3 — mitotic arrest and cytoskeletal collapse, legible
  because the features were grouped.
- **Positive-control labels are noisy.** Lapatinib, pelitinib and canertinib carry the
  `tubulin_binder` flag and top the E-distance ranking; they are EGFR inhibitors. Inspect
  labels before training on them.

## Install

```bash
mamba env create -f environment.yml     # RAPIDS 26.08 / CUDA 13 / Python 3.13
mamba activate cellpainting
pip install -e ".[tutorials]"
```

CPU-only is fine — nothing in these three notebooks needs a GPU.
`cellpainting-lock.txt` records the exact solved environment used to produce the numbers
above (scanpy 1.12.4, anndata 0.13.3, decoupler 2.2.0, pertpy 1.3.0).

## The `cellpainting_scverse` helper package

Small, installable, and used by the notebooks:

- `parse_feature_names` — rebuilds the `compartment`/`family`/`channel`/`channel_2`/
  `n_channels` schema from raw CellProfiler names. Verified to agree with the curated
  CPJUMP1 and EUbOPEN releases on **all 5,846 features**, so profiles you build yourself
  stay interoperable with them.
- `feature_sets` / `summarize` — turn `var` annotations into a `decoupler` network.
- `profiles_to_anndata`, `sanitize_obs`, `split_features_metadata` — profile table to
  AnnData.
- `drop_image_features` — see the trap below.

## A trap worth knowing about

`pycytominer >= 1.5.0` returns `Image_*` columns from `normalize()` **un-normalized**,
through a passthrough added for OME-Arrow image payloads. But in CellProfiler `Image_` is
the namespace for per-image *measurements*, and because `feature_select` only drops
features it inferred, those raw columns survive every filter and reach `X` on their
original scales. On JUMP TARGET2 this silently turns the documented 591 features into
1,680. Use `drop_image_features` before normalizing, or pass `image_features=True` to both
`normalize` and `feature_select`.

## Roadmap

- Quality control, on a dataset with per-well cell counts (EUbOPEN has them).
- Normalization and batch integration across EU-OS sites (`FMP`, `MEDINA`, `USC`).
- Dose-response with `pertpy` — needs a dose series, so EUbOPEN (10 concentrations)
  rather than EU-OS.
- Plate-layout artifacts with `squidpy.gr.spatial_autocorr`.
- Single-cell morphological profiles, which would unlock `Mixscape`, `Milo` and
  within-well response heterogeneity.

## Credits

Data: [EU-OS bioactives](https://github.com/schmiedc/EU-OS_bioactives)
([Zenodo](https://doi.org/10.5281/zenodo.13309566)), EU-OPENSCREEN. AnnData conversion:
[scverse/cell-painting-io](https://github.com/scverse/cell-painting-io). Data-preparation
notebooks: [timtreis/2024_broad_hackathon](https://github.com/timtreis/2024_broad_hackathon).
Profiling methodology: [pycytominer](https://github.com/cytomining/pycytominer).

## License

BSD 3-Clause. See [LICENSE](LICENSE).
