"""Treat morphological feature families as gene sets.

The single most useful consequence of annotating ``var`` is that the feature
axis gains structure: every CellProfiler feature belongs to a compartment, a
measurement family and (usually) a stain. That is exactly the shape decoupler
consumes for gene sets, which turns thousands of opaque features into a few
dozen interpretable morphological programmes.
"""

from __future__ import annotations

from collections.abc import Sequence

import pandas as pd

DEFAULT_KEYS: tuple[str, ...] = ("compartment", "family", "channel")


def feature_sets(
    adata,
    keys: Sequence[str] = DEFAULT_KEYS,
    combine: bool = False,
    min_size: int = 5,
) -> pd.DataFrame:
    """Build a decoupler-style network from ``var`` annotations.

    Parameters
    ----------
    adata
        AnnData whose ``var`` carries the annotation columns in `keys`.
    keys
        ``var`` columns to turn into feature sets.
    combine
        If True, also emit the crossed sets (e.g. ``cells|intensity|dna``),
        which are more specific but much smaller.
    min_size
        Drop sets with fewer than this many features.

    Returns
    -------
    A long DataFrame with ``source`` (set name) and ``target`` (feature name)
    columns, plus a unit ``weight``, ready to hand to ``decoupler``.
    """
    missing = [k for k in keys if k not in adata.var.columns]
    if missing:
        raise KeyError(
            f"var is missing annotation columns {missing}. "
            "Build var with cellpainting_scverse.io.parse_feature_names, or use a "
            "dataset that ships the compartment/family/channel schema."
        )

    var = adata.var
    records = []
    for key in keys:
        values = var[key].astype(object)
        for value, idx in var.groupby(values, observed=True).groups.items():
            if pd.isna(value):
                continue
            records.append(
                pd.DataFrame({"source": f"{key}:{value}", "target": list(idx)})
            )

    if combine:
        crossed = (
            var[list(keys)]
            .astype(object)
            .agg(
                lambda row: "|".join("NA" if pd.isna(v) else str(v) for v in row),
                axis=1,
            )
        )
        for value, idx in var.groupby(crossed, observed=True).groups.items():
            records.append(
                pd.DataFrame({"source": f"combined:{value}", "target": list(idx)})
            )

    if not records:
        return pd.DataFrame(columns=["source", "target", "weight"])

    net = pd.concat(records, ignore_index=True)
    net["weight"] = 1.0
    sizes = net.groupby("source")["target"].transform("size")
    return net[sizes >= min_size].reset_index(drop=True)


def summarize(net: pd.DataFrame) -> pd.DataFrame:
    """Feature-set sizes, largest first."""
    return (
        net.groupby("source")["target"]
        .size()
        .sort_values(ascending=False)
        .rename("n_features")
        .to_frame()
    )
