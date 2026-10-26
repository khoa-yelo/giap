The release check now reuses calibration evidence from an earlier commit only when every change
since then lies in documentation, CI, report or plot code (or is the version line alone), and says
so in the row. It also checks that each container image key is current in `images.json`, and the
release workflow verifies that every GHCR image pin is publicly pullable.
