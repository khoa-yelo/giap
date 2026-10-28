"""Internal model: Interval, Replicon, Gene, Feature (0-based half-open coordinates)."""

from gippa.model.coords import (
    Interval,
    Replicon,
    RepliconKind,
    Strand,
    Topology,
    check_interval,
    from_one_based_closed,
    interval_sort_key,
    split_at_origin,
    to_one_based_closed,
    wrap_position,
)
from gippa.model.genes import (
    Gene,
    GeneModel,
    GeneType,
    NcFeature,
    NcType,
    gene_model_from_bakta_json,
)

__all__ = [
    "Gene",
    "GeneModel",
    "GeneType",
    "Interval",
    "NcFeature",
    "NcType",
    "Replicon",
    "RepliconKind",
    "Strand",
    "Topology",
    "check_interval",
    "from_one_based_closed",
    "gene_model_from_bakta_json",
    "interval_sort_key",
    "split_at_origin",
    "to_one_based_closed",
    "wrap_position",
]
