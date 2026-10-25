Releases are built and published by a manual release workflow: it checks the version, the release
readiness record and the changelog, builds and smoke-tests the wheel and sdist, uploads them to PyPI
by trusted publishing and creates the tagged GitHub release with these notes.
