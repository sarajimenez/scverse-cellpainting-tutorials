import numpy as np
import pandas as pd
import pytest

from cellpainting_scverse import (
    drop_image_features,
    feature_sets,
    parse_feature_names,
    profiles_to_anndata,
    sanitize_obs,
)

# Expectations taken from the curated CPJUMP1 release, so our parser stays
# schema-compatible with it: (feature, compartment, family, channel, channel_2, n)
CASES = [
    ("Cells_AreaShape_Compactness", "cells", "areashape", None, None, 0),
    (
        "Cells_Intensity_IntegratedIntensityEdge_AGP",
        "cells",
        "intensity",
        "agp",
        None,
        1,
    ),
    ("Cells_Correlation_Correlation_AGP_DNA", "cells", "correlation", "agp", "dna", 2),
    (
        "Cells_Texture_AngularSecondMoment_DNA_10_01_256",
        "cells",
        "texture",
        "dna",
        None,
        1,
    ),
    ("Cells_Granularity_10_AGP", "cells", "granularity", "agp", None, 1),
    (
        "Cells_RadialDistribution_FracAtD_AGP_3of4",
        "cells",
        "radialdistribution",
        "agp",
        None,
        1,
    ),
    (
        "Cells_Neighbors_FirstClosestDistance_Adjacent",
        "cells",
        "neighbors",
        None,
        None,
        0,
    ),
    ("Nuclei_Intensity_MeanIntensity_DNA", "nuclei", "intensity", "dna", None, 1),
    (
        "Cytoplasm_Correlation_Correlation_AGP_Mito",
        "cytoplasm",
        "correlation",
        "agp",
        "mito",
        2,
    ),
    # brightfield planes: CPJUMP1 spells these HighZBF/LowZBF and the curated
    # releases normalize them to brightfield_high / brightfield_low
    (
        "Cells_Correlation_Correlation_DNA_HighZBF",
        "cells",
        "correlation",
        "dna",
        "brightfield_high",
        2,
    ),
    (
        "Cells_Correlation_Correlation_LowZBF_RNA",
        "cells",
        "correlation",
        "brightfield_low",
        "rna",
        2,
    ),
    (
        "image_granularity_10_bfhigh",
        "image",
        "granularity",
        "brightfield_high",
        None,
        1,
    ),
    # lowercased names, as shipped by the EUbOPEN release
    ("image_granularity_10_agp", "image", "granularity", "agp", None, 1),
]


@pytest.mark.parametrize(("name", "compartment", "family", "ch1", "ch2", "n"), CASES)
def test_parse_feature_names(name, compartment, family, ch1, ch2, n):
    var = parse_feature_names([name])
    row = var.loc[name]
    assert row.compartment == compartment
    assert row.family == family
    assert row.n_channels == n
    if ch1 is None:
        assert pd.isna(row.channel)
    else:
        assert row.channel == ch1
    if ch2 is None:
        assert pd.isna(row.channel_2)
    else:
        assert row.channel_2 == ch2


def test_parse_feature_names_schema():
    var = parse_feature_names([c[0] for c in CASES])
    assert list(var.columns) == [
        "compartment",
        "family",
        "channel",
        "channel_2",
        "n_channels",
    ]
    assert var.index.name == "feature"
    assert var.n_channels.dtype == np.int8
    assert str(var.compartment.dtype) == "category"


def test_unknown_compartment_falls_back():
    var = parse_feature_names(["Barcode_Foo_Bar"])
    assert var.loc["Barcode_Foo_Bar"].compartment == "other"


def _toy_profiles(n=40):
    feats = [c[0] for c in CASES]
    df = pd.DataFrame(np.random.default_rng(0).random((n, len(feats))), columns=feats)
    df["Metadata_Plate"] = "p1"
    df["Metadata_moa"] = ["kinase inhibitor"] * (n // 2) + [np.nan] * (n - n // 2)
    df["Metadata_mixed"] = ["a"] * (n // 2) + [1] * (n - n // 2)
    df.index = [f"p1__W{i}" for i in range(n)]
    return df


def test_sanitize_obs_makes_mixed_columns_writable(tmp_path):
    import anndata as ad

    df = _toy_profiles()
    adata = ad.AnnData(np.zeros((len(df), 1), dtype=np.float32))
    adata.obs = sanitize_obs(df[["Metadata_Plate", "Metadata_moa", "Metadata_mixed"]])
    adata.obs_names = df.index
    out = tmp_path / "t.h5ad"
    adata.write(out)
    assert ad.read_h5ad(out).obs.shape == (len(df), 3)


def test_profiles_to_anndata_roundtrip():
    df = _toy_profiles()
    adata = profiles_to_anndata(df)
    assert adata.shape == (len(df), len(CASES))
    assert adata.X.dtype == np.float32
    assert "compartment" in adata.var
    assert list(adata.obs.columns) == [
        "Metadata_Plate",
        "Metadata_moa",
        "Metadata_mixed",
    ]


def test_drop_image_features():
    df = _toy_profiles()
    assert any(c.startswith("Image_") or c.startswith("image_") for c in df.columns)
    out = drop_image_features(
        df.rename(columns={"image_granularity_10_agp": "Image_Granularity_10_AGP"})
    )
    assert not [c for c in out.columns if c.startswith("Image_")]


def test_feature_sets_builds_decoupler_net():
    adata = profiles_to_anndata(_toy_profiles())
    net = feature_sets(adata, min_size=1)
    assert set(net.columns) == {"source", "target", "weight"}
    assert net.source.str.startswith("compartment:").any()
    assert net.source.str.startswith("family:").any()
    assert set(net.target) <= set(adata.var_names)


def test_feature_sets_requires_annotations():
    import anndata as ad

    adata = ad.AnnData(np.zeros((3, 2), dtype=np.float32))
    with pytest.raises(KeyError, match="missing annotation columns"):
        feature_sets(adata)
