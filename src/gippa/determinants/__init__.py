"""Module 8: what explains abundance (promoter, RBS, codon usage)."""

__all__ = ["CAVEATS"]

#: §7.8.2 "Caveats (printed in the report)", verbatim, in order (rendered by the M9 report)
CAVEATS: tuple[str, ...] = (
    "In a synthetic GFP library, mRNA folding near the RBS explained over half the variance; "
    "native genes add regulation and turnover, so expect a lower R².",
    "All features are predictions, and the data come from one condition.",
    "CAI, tAI and GC3 are collinear, so interpret blocks rather than single features.",
)
