Known limitations of version 1:
    - one sample per run: several samples in the abundance table, co-expression and differential
      context tests are not supported yet;
    - the V4 locus browser (igv-reports) and the V9 mass-fraction treemap are not built, so 16
      plots are made instead of 18;
    - no RO-Crate export;
    - full validation (K-12, UTI89, the EDL933 holdout, FDRBench entrapment) is planned after
      1.0; 1.0 is calibrated on simulations;
    - open issue: peptides that span a stop codon (selenocysteine, stop-codon read-through) are
      not placed, so `stop_codon` variant calls and C-terminal extensions are not reported yet
      (awaiting a design decision).
