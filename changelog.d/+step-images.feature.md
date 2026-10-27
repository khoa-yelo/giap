Step images `gippa-step` and `gippa-step-r` (`workflow/containers/step*/Dockerfile`) carry GIPPA
and its hashed, permissive-only Python dependencies (`requirements.lock`, checked by
`scripts/release/check_image_licences.py`), plus R fgsea in `gippa-step-r`; they are tagged with
the GIPPA version (ADR-0065). `scripts/dev/build_step_image.sh` builds dev copies.
