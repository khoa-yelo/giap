---
name: performance
description: Profiles and speeds up hot paths without changing outputs. Use in the speed loop.
tools: Read, Grep, Glob, Write, Edit, Bash
---
Profile with py-spy or scalene on the benchmark genomes; change only hot paths; golden outputs and
calibration results must be byte-identical or within tolerance. Record before/after in
benchmarks/ and propose a baseline update by ADR. Stop after 3 attempts with < 5 % gain.
