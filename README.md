# scverse Cell Painting tutorials

Tutorials and helpers for running **downstream analysis on Cell Painting data with
scverse tools** — `anndata`, `scanpy`, `decoupler`, `pertpy`, `squidpy`, and
`rapids-singlecell` for GPU acceleration.

Started at the scverse Cell Painting hackathon. The goal is to take the analysis
patterns that are well established for single-cell transcriptomics — as laid out in
[single-cell best practices](https://www.sc-best-practices.org/) — and work out which
of them transfer to image-based morphological profiling, which need adapting, and
which do not transfer at all.

## The idea in one substitution

**A well (or a segmented cell) is an observation, and a CellProfiler feature is a
gene.** Where that holds, the scverse stack works on Cell Painting data unchanged.
Where it breaks is where the interesting work is.

| single-cell concept | Cell Painting equivalent | Transfers? |
| --- | --- | --- |
| Cell | Well (well-level) or segmented cell | yes |
| Gene | Morphological feature (compartment x family x channel) | yes, but features are correlated by construction |
| Count matrix | Continuous feature matrix | **no** — no counts, no sparsity, no overdispersion |
| Ambient RNA / doublets | Out-of-focus fields, debris, mis-segmentation | different failure modes, same QC role |
| Size-factor normalisation + `log1p` | Per-plate MAD-robustize against DMSO controls | **no** — controls replace library size |
| Highly variable gene selection | `pycytominer.feature_select` | yes |
| Batch (donor, 10x run) | Plate, and imaging source | yes, and sources differ in microscope hardware |
| Harmony / scVI integration | Sphering (ZCA whitening on controls) | yes, and which wins is an open question |
| PCA -> kNN -> UMAP -> Leiden | identical | yes |
| Cell-type annotation | Mechanism of action / target annotation | yes, but labels are sparse and noisy |
| Marker genes per cluster | Discriminative features per cluster | yes |
| Pathway / GSEA over gene sets | Enrichment over **feature families** from `var` | yes, once `var` is annotated |
| Perturbation modelling (`pertpy`) | Dose-response, replicate reproducibility | yes — this is the native question |
| Spatial coordinates | **Well row/column on the plate** | yes, and it is underused |

Two of these are worth spelling out, because they are what the scverse stack adds over
the existing Cell Painting tooling:

- **Feature families are gene sets.** Once `var` carries `compartment`, `family` and
  `channel`, thousands of opaque features collapse into a few dozen interpretable
  morphological programmes that `decoupler` can score. See
  `cellpainting_scverse.feature_sets`.
- **Plate layout is a spatial graph.** Well row and column are real 2D coordinates, so
  edge effects and evaporation gradients are *spatial autocorrelation* — which makes
  `squidpy.gr.spatial_autocorr` a principled plate-artifact detector rather than the
  ad-hoc heuristics the field currently uses.

## Datasets

| Dataset | Shape | What makes it useful |
| --- | --- | --- |
| **CPJUMP1** ([Chandrasekaran et al. 2024](https://www.nature.com/articles/s41592-024-02241-6)) | 19,498 x 903 | Three perturbation modalities (compound, CRISPR, ORF), two cell lines, four timepoints — built for benchmarking |
| **EUbOPEN** Cell Painting | 39,206 x 4,943 | 10 concentration levels for dose-response, and per-well cell counts in `obs` |
| **JUMP `cpg0016` TARGET2** | 64,464 x 591 | Built from raw well-level profiles by the data-preparation notebook |

Both curated releases ship the same `var` schema (`compartment`, `family`, `channel`,
`channel_2`, `n_channels`). `cellpainting_scverse.parse_feature_names` reproduces it
from raw CellProfiler feature names — verified to agree on **100% of all 5,846
features** across both — so profiles you build yourself stay interoperable with them.

## Install

```bash
mamba env create -f environment.yml     # RAPIDS 26.08 / CUDA 13 / Python 3.13
mamba activate cellpainting
pip install -e ".[tutorials]"
```

CPU-only is fine for everything except the `rapids-singlecell` sections:

```bash
pip install -e ".[tutorials]"
```

`cellpainting-lock.txt` records the exact solved environment that produced the
verified results.

## Tutorials

| Notebook | Chapter it mirrors | Status |
| --- | --- | --- |
| Data preparation: JUMP profiles to AnnData | — | working (see `notebooks/`) |
| Quality control | QC | planned |
| Normalization and batch integration | Normalization, Integration | planned |
| Dimensionality reduction and clustering | Dim. reduction, Clustering | planned |
| Feature-family enrichment with `decoupler` | Pathway analysis | planned |
| Perturbation analysis with `pertpy` | Perturbation modelling | planned |
| Plate-layout artifacts with `squidpy` | Spatial statistics | planned |

## A trap worth knowing about

`pycytominer >= 1.5.0` returns `Image_*` columns from `normalize()` **un-normalized**,
through a passthrough added for OME-Arrow image payloads. But in CellProfiler `Image_`
is the namespace for per-image *measurements*, and because `feature_select` only drops
features it inferred, those raw columns survive every filter and reach `X` on their
original scales.

On JUMP TARGET2 this silently turns the documented 591 features into 1,680
(591 morphology plus 1,089 raw `Image_*`). Use
`cellpainting_scverse.drop_image_features` before normalizing, or pass
`image_features=True` to both `normalize` and `feature_select` so they get normalized
and filtered like everything else.

Worth knowing which you want: with proper normalization, `Image_*` features are *not*
more batch-confounded than compartment features (median variance explained by source
0.043 vs 0.037) and carry more compound-associated variance (0.340 vs 0.239). But
`Image_ImageQuality_*` and `Image_Threshold_*` are QC diagnostics and segmentation
internals, so they belong in `obs`, not `X`.

## Credits

The data-preparation notebooks are derived from
[timtreis/2024_broad_hackathon](https://github.com/timtreis/2024_broad_hackathon)
(see `notebooks/README.md`). Data comes from the
[Cell Painting Gallery](https://github.com/broadinstitute/cellpainting-gallery)
(`cpg0016-jump`), produced by the [JUMP-Cell Painting Consortium](https://jump-cellpainting.broadinstitute.org/).
Profiling methodology follows [pycytominer](https://github.com/cytomining/pycytominer).

## License

BSD 3-Clause. See [LICENSE](LICENSE).
