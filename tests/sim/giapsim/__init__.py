"""giapsim: synthetic-truth generator for GIAP tests (GIAP_SPEC.md §11.3).

M0 skeleton: one circular chromosome, protein-coding genes in transcription units, a Bakta-like
JSON annotation and a gene truth table.
"""

from giapsim.simulate import SimConfig, SimResult, simulate
from giapsim.strain import StrainConfig

__all__ = ["SimConfig", "SimResult", "StrainConfig", "simulate"]
