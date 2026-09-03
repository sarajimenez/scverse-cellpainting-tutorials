# Notebooks

## Data preparation (derived from upstream)

These two are adapted from
[timtreis/2024_broad_hackathon](https://github.com/timtreis/2024_broad_hackathon)
and kept close to the originals so the results stay comparable:

| Notebook | What it does |
| --- | --- |
| `1_prepare_wellres_TARGET2_CellProfiler_data.ipynb` | Pulls well-level CellProfiler profiles for the 141 JUMP TARGET2 plates from the public `cellpainting-gallery` S3 bucket, imputes missing values per source, and joins well -> compound -> standardised compound metadata -> JUMP-MOA + CLUE MOA -> microscope config. |
| `2c+3c_pycytonorm_wellres_TARGET2_CellProfiler_by_plate.ipynb` | The [pycytominer profiling workflow](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10516049/): per-plate MAD-robustize against DMSO, feature selection, global sphering, then out to AnnData with UMAP. |

### The one change from upstream

Notebook 2 drops `Image_*` columns immediately after loading:

```python
target2_complete = target2_complete.drop(
    columns=[c for c in target2_complete.columns if c.startswith("Image_")]
)
```

Without this, `pycytominer >= 1.5.0` passes 1,089 raw per-image measurements straight
through `normalize()` and every `feature_select` filter, and the pipeline yields 1,680
features instead of the documented 591. See the trap section in the top-level README.

Verified end to end on pycytominer 1.7.1 / pandas 3.0.5 / anndata 0.13.3:
`normalize` -> 64,464 x 3,699, `feature_select` -> 64,464 x 617,
AnnData -> **64,464 x 591**, whole notebook in 74 s.

### Running them

The JUMP metadata CSVs are not vendored here. Fetch them from upstream first:

```bash
git clone --depth 1 https://github.com/timtreis/2024_broad_hackathon.git /tmp/bh
cp -r /tmp/bh/metadata ..
```

Then create `../data` (the notebooks read and write `../data/...`), pointing it at
wherever you keep large files:

```bash
mkdir -p /path/to/big/disk/cellpainting
ln -s /path/to/big/disk/cellpainting ../data
```

Notebook 1 downloads roughly 1 GB from S3 anonymously; budget about a minute with a
thread pool. Both notebooks together run in well under an hour on CPU.

One caveat if you script them: `jupyter nbconvert --execute` rejects these files with
`NotebookValidationError: 'metadata' is a required property` because of a stored
output inherited from upstream. Clear outputs first. Jupyter itself is unaffected.

## Tutorials (in progress)

The downstream tutorials will run on the curated CPJUMP1 and EUbOPEN AnnData objects
rather than on TARGET2, since those already carry annotated `var` and richer
perturbation metadata. See the tutorial table in the top-level README.
