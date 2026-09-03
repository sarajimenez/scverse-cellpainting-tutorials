"""Get Cell Painting profile tables into an AnnData that scverse tools accept."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import pandas as pd

# Canonical channel names, matching the curated CPJUMP1 and EUbOPEN releases:
# the five Cell Painting stains plus the brightfield planes some sources acquire.
CHANNELS: tuple[str, ...] = (
    "dna",
    "agp",
    "er",
    "mito",
    "rna",
    "brightfield",
    "brightfield_high",
    "brightfield_low",
)

# Raw tokens as they appear in CellProfiler feature names, mapped to the
# canonical names above. CPJUMP1 writes ``HighZBF``/``LowZBF`` while the
# lowercased EUbOPEN names use ``bfhigh``/``bflow``; both mean the same plane.
CHANNEL_ALIASES: dict[str, str] = {
    "dna": "dna",
    "agp": "agp",
    "er": "er",
    "mito": "mito",
    "rna": "rna",
    "brightfield": "brightfield",
    "highzbf": "brightfield_high",
    "lowzbf": "brightfield_low",
    "bfhigh": "brightfield_high",
    "bflow": "brightfield_low",
}

COMPARTMENTS: tuple[str, ...] = ("cells", "cytoplasm", "nuclei", "image")


def parse_feature_names(features: Iterable[str]) -> pd.DataFrame:
    """Decompose CellProfiler feature names into a ``var`` annotation frame.

    Feature names follow ``<compartment>_<family>_<measurement...>[_<channel>...]``,
    for example ``Cells_AreaShape_Compactness``,
    ``Nuclei_Intensity_MeanIntensity_DNA`` or
    ``Cytoplasm_Correlation_Correlation_AGP_Mito``.

    The returned columns match the schema used by the curated CPJUMP1 and EUbOPEN
    AnnData objects (``compartment``, ``family``, ``channel``, ``channel_2``,
    ``n_channels``), so profiles built here can be concatenated with those.
    Values are lowercased; ``channel``/``channel_2`` are NaN where a feature is
    not channel-specific (shape and neighbour measurements).
    """
    features = list(features)
    rows = []
    for name in features:
        tokens = str(name).split("_")
        head = tokens[0].lower() if tokens else "other"
        compartment = head if head in COMPARTMENTS else "other"
        family = tokens[1].lower() if len(tokens) > 1 else "other"
        # Scan the measurement tail in order, so correlation features keep the
        # channel pair in the order CellProfiler wrote it.
        channels = [
            CHANNEL_ALIASES[t.lower()]
            for t in tokens[2:]
            if t.lower() in CHANNEL_ALIASES
        ]
        rows.append(
            (
                compartment,
                family,
                channels[0] if len(channels) > 0 else np.nan,
                channels[1] if len(channels) > 1 else np.nan,
                len(channels),
            )
        )

    var = pd.DataFrame(
        rows,
        columns=["compartment", "family", "channel", "channel_2", "n_channels"],
        index=pd.Index(features, name="feature"),
    )
    for col in ("compartment", "family", "channel", "channel_2"):
        var[col] = var[col].astype("category")
    var["n_channels"] = var["n_channels"].astype("int8")
    return var


def sanitize_obs(obs: pd.DataFrame, max_categories: int = 2000) -> pd.DataFrame:
    """Coerce a metadata frame into dtypes ``AnnData.write_h5ad`` accepts.

    Numeric and boolean columns are left untouched. Anything else becomes
    ``category`` with NaN mapped to the literal ``"NA"``. This is only needed
    when a column genuinely mixes types -- which the JUMP metadata joins do not,
    but other join paths can -- where anndata otherwise raises
    ``TypeError: Can't implicitly convert non-string objects to strings``.
    """
    out = obs.copy()
    for col in out.columns:
        series = out[col]
        if pd.api.types.is_numeric_dtype(series) or pd.api.types.is_bool_dtype(series):
            continue
        filled = series.astype(object).where(series.notna(), "NA").astype(str)
        out[col] = (
            filled.astype("category") if filled.nunique() <= max_categories else filled
        )
    return out


def split_features_metadata(
    df: pd.DataFrame, metadata_prefix: str = "Metadata"
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split a pycytominer-style profile table into ``(features, metadata)``."""
    meta_cols = [c for c in df.columns if metadata_prefix in c]
    feat_cols = [c for c in df.columns if c not in meta_cols]
    return df[feat_cols], df[meta_cols]


def profiles_to_anndata(
    df: pd.DataFrame,
    metadata_prefix: str = "Metadata",
    sanitize: bool = True,
):
    """Build an AnnData from a well-level profile table.

    ``X`` is the feature block (wells x features) as float32, ``obs`` is the
    ``Metadata_*`` block, and ``var`` is parsed from the feature names so the
    feature axis is queryable the way a gene axis is in transcriptomics.
    """
    import anndata as ad

    feats, meta = split_features_metadata(df, metadata_prefix)
    adata = ad.AnnData(
        X=np.ascontiguousarray(feats.to_numpy(dtype=np.float32)),
        obs=sanitize_obs(meta) if sanitize else meta.copy(),
        var=parse_feature_names(feats.columns),
    )
    adata.obs_names = df.index.astype(str)
    return adata


def drop_image_features(df: pd.DataFrame) -> pd.DataFrame:
    """Drop per-image (``Image_*``) CellProfiler measurements.

    pycytominer >= 1.5.0 returns ``Image_*`` columns from ``normalize()``
    *un-normalized*, via a passthrough intended for OME-Arrow image payloads.
    Since ``feature_select`` only drops features it inferred, those raw columns
    then survive every filter and reach ``X`` on their original scales.

    Call this before ``normalize`` to reproduce the pre-1.5.0 behaviour, or pass
    ``image_features=True`` to both ``normalize`` and ``feature_select`` to keep
    them and have them normalized properly.
    """
    return df.drop(columns=[c for c in df.columns if c.startswith("Image_")])
