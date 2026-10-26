Tool images on GHCR now carry one immutable tag per build, `<version>-<key>`, where the key hashes
the image's Dockerfile, env file and locks; an existing tag is never rebuilt or overwritten. The
images install from the committed conda and pip locks, and every downloaded source archive is
checked against its sha256. BacTermFinder now has locks too.
