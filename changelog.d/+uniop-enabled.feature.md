`annotation.uniop.enabled: true | false | auto` (default `true`). Preflight resolves it before the
run: it runs UniOP's `--help` once to check that UniOP works, and `auto` turns into `true` or
`false` in `resolved_config.yaml`. `giap run` then passes `params.run_uniop`, and UNIOP is skipped
when it is false. The distance-model fallback is always reported, with its reason:
- preflight warns for `false`, for `auto` resolving to false and for a failed probe, and quotes
  the last failed `image:uniop` setup;
- under `--publication`, `auto` without UniOP is an error;
- `giap run` prints a closing warning line;
- the manifest gains `operon_source` and `operon_source_reason`;
- the report states the distance model's precision and recall.

ADR-0061 row flags:
- TU-block Moran's I rows carry `tu_block_tu_error` on every run;
- under the distance model, those rows and every `both_detected_t` pair row also carry
  `tu_source_distance_model`, with their caveats in the report.
