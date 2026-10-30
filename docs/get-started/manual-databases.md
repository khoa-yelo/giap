# Download the databases by hand

`gippa setup` normally downloads every database itself ([Set up the databases](setup-databases.md)).
Use this page when it cannot do that:

- the compute nodes have no internet access, or a firewall or proxy blocks the sources;
- a download failed or timed out (setup does not resume a download; it starts it again);
- your site already holds a shared copy, for example of the Bakta database;
- you installed only some modules and now need another one.

Every URL, version, checksum and path below comes from `src/gippa/data/db_sources.yaml`, the file
`gippa setup` itself reads. Download only from these official sources, never from a mirror or an
unofficial copy. GIPPA's own downloader refuses every other host.

## How `gippa setup` adopts files you placed

`gippa setup` checks each item it plans before it downloads anything (`install` in
`src/gippa/setupdb/run.py`). If a copy is already at the item's path under `--db-dir` and it
verifies, setup records it in `DIR/db_manifest.json` and downloads nothing. It does not
contact the source for that item. The checks run offline (`verify` in `src/gippa/setupdb/items.py`):

| Item | What setup checks in the copy you placed |
|---|---|
| Bakta light / full | `bakta/db-light/version.json` (or `bakta/db/version.json`): `major` 6 and `type` `light` (or `full`) |
| geNOMAD | `genomad_db/version.txt` is `1.9` (or starts with `1.9.`) |
| Promotech | `promotech/RF-HOT.model` and `promotech/RF-HOT.zip` both exist, and the zip's SHA-256 equals the pin |
| antiSMASH | `antismash/` exists and is not empty (presence only) |
| GO | the file's SHA-256 equals the pin, and its `data-version:` header names release `2026-05-19` |
| COG | the file's SHA-256 equals the pin |
| Contaminants | the file's SHA-256 equals the pin, and the file holds 381 FASTA records |
| KEGG | **cannot be adopted**: see [KEGG](#kegg) |
| UniOP archive | the archive's SHA-256 equals the pin (`setup --uniop` adopts it, then builds) |
| Apptainer images | `<cache>/<name>.img` exists (`setup --profile apptainer` skips the pull) |

Three things to know:

- **`--check` only reports. It does not record anything.** Preflight and `gippa run` read
  `db_manifest.json`, not the files. After you place files, run `gippa setup` once *without*
  `--check`, with the opt-outs for the items you did not place. Use `--check` afterwards to
  confirm.
- Every item is on by default, so setup plans all of them unless you opt out. On a node with no
  network, an item that is missing or does not verify still fails with `source <host>
  unreachable`, and setup exits 2. The items that did verify are still recorded. To keep the
  output clean, leave out what you did not place: `--bakta skip`, `--no-genomad`,
  `--no-antismash`, `--no-promotech-models`, `--no-go`, `--no-kegg`, `--no-cog` and
  `--no-contaminants`. Add
  `--no-images` (setup otherwise pulls the images when it finds a container runtime) and
  `--no-nextflow` (setup otherwise fetches the Nextflow launcher when none is on `PATH`).
- For Bakta and geNOMAD, setup reads only the version file. It cannot re-check the data, because
  the archive is gone after unpacking. **Verify the archive checksum before you unpack.**
  Promotech keeps its zip, so setup re-hashes that zip on every check.

When setup adopts a copy, it records the date of that setup run as `retrieved` in the manifest.

## The workflow

1. On a machine with internet access (a transfer node or your workstation; follow your site's
   rules for login nodes), download each archive or file from the table below. Verify its checksum.
2. Put it at its path under the DB directory, and unpack it there.
3. On the cluster, record the items. For example, after placing Bakta light, geNOMAD and
   Promotech (all three are on by default; the opt-outs leave out the rest):

    ```bash
    gippa setup --db-dir /path/to/db --no-antismash --no-go --no-kegg --no-cog \
        --no-contaminants --no-images --no-nextflow
    ```

4. Confirm, and check the run's tools:

    ```bash
    gippa setup --db-dir /path/to/db --no-antismash --no-go --no-kegg --no-cog \
        --no-contaminants --check
    GIPPA_DB_DIR=/path/to/db gippa preflight --genome genome.fna --peptides peptides.tsv \
        --profile apptainer
    ```

    Exit code 0 means every planned item is `ok`. Exit code 2 lists each missing or mismatched item
    and the reason. Preflight takes the DB directory from `GIPPA_DB_DIR`; `gippa run` takes it from
    `--db-dir`.

Setup needs write access to the DB directory: it writes `db_manifest.json` and the lock file
`.giap_setup.lock` there.

## The items

The paths are relative to the DB directory (`--db-dir`, or `GIPPA_DB_DIR`). "Free space" is the
`size_gib` that setup asks for before it downloads an item. It is not the size on disk; for that,
see [Layout and sizes](setup-databases.md#layout-and-sizes).

| Item (`db_sources.yaml`) | `gippa setup` flag | Path | Version | Free space | Needed by |
|---|---|---|---|---|---|
| `bakta_light` | `--bakta light` (default) | `bakta/db-light` | 6.0 | 5 GiB | module `core` (annotation), unless `--annotation` is given |
| `bakta_full` | `--bakta full` | `bakta/db` | 6.0 | 84 GiB | the same, in place of the light DB |
| `genomad` | on; `--no-genomad` | `genomad_db` | 1.9 | 4 GiB | module `mobile` |
| `promotech_models` | on; `--no-promotech-models` | `promotech` | RF-HOT | 5 GiB | module `core` (promoter predictions) |
| `antismash` | on; `--no-antismash` | `antismash` | 8.0.4 | 10 GiB | module `bgc` |
| `go_basic`, `goslim_generic` | on; `--no-go` | `go/go-basic.obo`, `go/goslim_generic.obo` | 2026-05-19 | — | module `function` |
| `kegg_ko_pathway`, `kegg_pathway_names` | on; `--no-kegg` | `kegg/ko_pathway.tsv`, `kegg/pathway_names.tsv` | retrieval date | — | module `function` |
| `cog2020_def`, `cog2020_fun` | on; `--no-cog` | `cog/cog-20.def.tab`, `cog/fun-20.tab` | 2020 | — | module `function`; `cog_category` in `genes.tsv` |
| `contaminants` | on; `--no-contaminants` | `contaminants/contaminants.fasta` | commit `bddf6da9` | — | contaminant labels in mapping; `gippa makedb` |
| `uniop_source` | with the images, or `--uniop` alone | `uniop/UniOP-3ed14dfc2a5d69b993d99a26e50c1f8e9723c087.tar.gz` | commit `3ed14dfc` | — | operon predictions (UniOP; optional, see [Install › UniOP](install.md#uniop)) |

Preflight requires every part of a group to be `ok`: `function` needs both GO files, both KEGG
tables and both COG tables.

In the commands below, set `DB=/path/to/db` first. `md5sum -c` and `sha256sum -c` print `OK` for a
good file. If a check fails, delete the file and download it again from the same URL. Never
replace the pin.

### Bakta database v6.0 (light or full)

Source: Zenodo record 14916843. The archive unpacks into `db-light/` (light) or `db/` (full).

| Type | URL | Checksum (MD5) |
|---|---|---|
| light | `https://zenodo.org/records/14916843/files/db-light.tar.xz` | `4a6e059ded39e9c5537ef4137d2f5648` |
| full | `https://zenodo.org/records/14916843/files/db.tar.xz` | `4c1115e40abfa2b464ae5dd988bdd88e` |

```bash
mkdir -p "$DB/bakta" && cd "$DB/bakta"
curl -fL -C - -o db-light.tar.xz https://zenodo.org/records/14916843/files/db-light.tar.xz
echo "4a6e059ded39e9c5537ef4137d2f5648  db-light.tar.xz" | md5sum -c -
tar -xJf db-light.tar.xz && rm db-light.tar.xz      # creates bakta/db-light/
grep -E '"(major|type)"' db-light/version.json      # "major": 6, "type": "light"
```

For the full DB, use `db.tar.xz`, its MD5, and `bakta/db/`. With `curl -C -`, an interrupted
download resumes where it stopped. That helps with the 84 GiB full archive.

The archive contains an AMRFinderPlus database under `amrfinderplus-db/`. Setup records it as the
item `amrfinderplus`, with the version in `amrfinderplus-db/latest/version.txt`. When setup
downloads Bakta itself, `bakta_db download` also runs `amrfinder_update`. A copy you unpack by hand
keeps the AMRFinderPlus release it shipped with (`2024-12-18.1`, listed in `version.json`).

### geNOMAD database v1.9

Source: Zenodo record 14886553. The archive unpacks into `genomad_db/`.

```bash
cd "$DB"
curl -fL -C - -o genomad_db_v1.9.tar.gz https://zenodo.org/records/14886553/files/genomad_db_v1.9.tar.gz
echo "67244b528bb8bed464d1ca147136d33e  genomad_db_v1.9.tar.gz" | md5sum -c -
tar -xzf genomad_db_v1.9.tar.gz && rm genomad_db_v1.9.tar.gz   # creates genomad_db/
cat genomad_db/version.txt                                      # 1.9
```

### Promotech RF-HOT model

Source: `https://www.cs.mun.ca/~lourdes/public/PromoTech_models/RF-HOT.zip`. Fallback source,
used by `gippa setup` only when `www.cs.mun.ca` cannot be reached: the Promotech authors' GitHub LFS
copy at commit `bdf2f3ce`,
`https://media.githubusercontent.com/media/BioinformaticsLabAtMUN/PromoTech/bdf2f3ce06ffd0e342607eed4bfa32e160816143/models/RF-HOT.zip`
([ADR-0069](../adr/ADR-0069-promotech-model-github-lfs-fallback.md)). The publisher gives no
checksum. The pin is the SHA-256 that GIPPA recorded on its first download
([ADR-0014](../adr/ADR-0014-db-checksum-sources.md)); both sources serve the same bytes. **Keep
the zip**: setup re-hashes it on every check.

When setup reports both sources unreachable, download the zip from either URL on a machine that
reaches it, put `RF-HOT.zip` and its unzipped `RF-HOT.model` in `$DB/promotech/`, then rerun
`gippa setup`; it adopts the pair once the zip's SHA-256 equals the pin.

```bash
mkdir -p "$DB/promotech" && cd "$DB/promotech"
curl -fL -C - -o RF-HOT.zip https://www.cs.mun.ca/~lourdes/public/PromoTech_models/RF-HOT.zip
# or, if www.cs.mun.ca is down:
# curl -fL -C - -o RF-HOT.zip https://media.githubusercontent.com/media/BioinformaticsLabAtMUN/PromoTech/bdf2f3ce06ffd0e342607eed4bfa32e160816143/models/RF-HOT.zip
echo "fcd96ffe5a43b579af965839ee2bb896de94de3601978d88bcabd32c4fb5093f  RF-HOT.zip" | sha256sum -c -
unzip -o RF-HOT.zip            # writes RF-HOT.model; keep RF-HOT.zip next to it
```

### antiSMASH databases (8.0.4)

Source: `https://dl.secondarymetabolites.org/releases/`. GIPPA pins no archive here: antiSMASH's
own downloader, `download-antismash-databases`, fetches the files and checks them. Run it from
the pinned antiSMASH image, on a machine with internet access. The commands below match what
`gippa setup` runs (`ToolRunner` in `src/gippa/setupdb/run.py`):

```bash
mkdir -p "$DB/antismash"
# Apptainer (the image file in your cache; see Container images below)
apptainer exec --bind "$DB/antismash" \
    "$CACHE/antismash-standalone@sha256-d3f390e687cfe02f1387be40256dddbf678e91bc56b96c36da07b56ec8bda29d.img" \
    download-antismash-databases --database-dir "$DB/antismash"
# or Docker
docker run --rm -u "$(id -u):$(id -g)" -v "$DB/antismash:$DB/antismash" \
    --entrypoint download-antismash-databases \
    antismash/standalone@sha256:d3f390e687cfe02f1387be40256dddbf678e91bc56b96c36da07b56ec8bda29d \
    --database-dir "$DB/antismash"
```

Setup checks only that `antismash/` is not empty. Before it downloads, `gippa setup` (antiSMASH
is on by default; `--no-antismash` leaves it out) also runs `antismash --check-prereqs` in the image. If that passes, it records the databases
shipped in the image (manifest URL `image`) and downloads nothing.

### Gene Ontology (release 2026-05-19)

```bash
mkdir -p "$DB/go" && cd "$DB/go"
curl -fL -o go-basic.obo https://purl.obolibrary.org/obo/go/releases/2026-05-19/go-basic.obo
curl -fL -o goslim_generic.obo https://purl.obolibrary.org/obo/go/releases/2026-05-19/subsets/goslim_generic.obo
sha256sum -c - <<'EOF'
923645000afcb1df0e72e252b64950b1c2bcc98c9393fec8c1dd1a0f7c9d4d77  go-basic.obo
dfc49ced033f5f047a9729a2b434e7a70b0c395ba61f950d138beff769555a26  goslim_generic.obo
EOF
```

`purl.obolibrary.org` redirects to `release.geneontology.org`. That is the only redirect GIPPA's
downloader accepts for GO.

### KEGG

KEGG publishes no versions or checksums. Setup downloads `https://rest.kegg.jp/link/pathway/ko` to
`kegg/ko_pathway.tsv` and `https://rest.kegg.jp/list/pathway` to `kegg/pathway_names.tsv`. It
records the SHA-256 of the bytes it received and the retrieval date as the version
([ADR-0017](../adr/ADR-0017-unversioned-sources.md)). On later runs, setup trusts a KEGG file only
when the file still matches the checksum recorded at download. **A KEGG file you place by hand
is reported `absent`** (`no recorded checksum for kegg/ko_pathway.tsv`), and setup tries to
download it again.

What works:

- Run only the KEGG step where `rest.kegg.jp` is reachable, against the same DB directory (for
  example, a transfer node that mounts the cluster filesystem):

    ```bash
    gippa setup --db-dir "$DB" --bakta skip --no-genomad --no-antismash --no-promotech-models \
        --no-go --no-cog --no-contaminants --no-images --no-nextflow
    ```

- Or run `gippa setup` (KEGG is on by default) on a machine with internet access, then copy the whole DB
  directory, `db_manifest.json` included ([Copy to the cluster](#download-on-a-machine-with-internet-copy-to-the-cluster)).
  The manifest records paths relative to the DB directory, so the copy verifies in its new place.

### COG 2020

```bash
mkdir -p "$DB/cog" && cd "$DB/cog"
curl -fL -o cog-20.def.tab https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2020/data/cog-20.def.tab
curl -fL -o fun-20.tab https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2020/data/fun-20.tab
sha256sum -c - <<'EOF'
00e93e42f38f6ada46ef6be5da4a8518829300ea29f2888261660de0ca2f7e77  cog-20.def.tab
fbdce778220f07a6f9701b2bde07fbcb1d3fca9fa7453fa9ab7021b5da18c4f7  fun-20.tab
EOF
```

### Contaminant library (Frankenfield et al. 2022, commit `bddf6da9`)

```bash
mkdir -p "$DB/contaminants" && cd "$DB/contaminants"
curl -fL -o contaminants.fasta \
  'https://raw.githubusercontent.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics/bddf6da9334f4300caed75f7a87c096dd96b6bca/Universal%20protein%20contaminant%20FASTA/0602_Universal%20Contaminants.fasta'
echo "c0782146be6dad6ecf53e23e5a21a7ed4fea9b4802cb01e8d00091404ec66926  contaminants.fasta" | sha256sum -c -
grep -c '^>' contaminants.fasta      # 381
```

### UniOP source archive (commit `3ed14dfc`)

UniOP has no licence file, so GIPPA never distributes it. Setup fetches the authors' archive,
checks it and builds UniOP on your machine ([Install › UniOP](install.md#uniop)). You can place the
archive yourself:

```bash
mkdir -p "$DB/uniop" && chmod 700 "$DB/uniop" && cd "$DB/uniop"
curl -fL -o UniOP-3ed14dfc2a5d69b993d99a26e50c1f8e9723c087.tar.gz \
  https://github.com/hongsua/UniOP/archive/3ed14dfc2a5d69b993d99a26e50c1f8e9723c087.tar.gz
echo "1f6b96f0a8659ad28e4f4b1ed49d538ded82f5f7669af74b265583d837a4108a  UniOP-3ed14dfc2a5d69b993d99a26e50c1f8e9723c087.tar.gz" | sha256sum -c -
chmod 600 UniOP-3ed14dfc2a5d69b993d99a26e50c1f8e9723c087.tar.gz
```

Then run `gippa setup --uniop --profile apptainer --db-dir "$DB" --accept-uniop-terms` (or
`--profile docker`). It adopts the archive and builds the image, but the build itself still
needs network access: the base image and the conda packages come from their registries. Build
where the network is, on a filesystem the run host also sees. Setup keeps `uniop/` at mode `0700`
and the archive at `0600`, so do not share it lab-wide. Without UniOP, the run falls back to the
distance model (`annotation.uniop.enabled` is `auto` by default).

## Container images

The pins are in `workflow/conf/containers.config`, bundled with the package. Pull each image by
its digest:

| Tool | Pinned reference | Version |
|---|---|---|
| antismash | `antismash/standalone@sha256:d3f390e687cfe02f1387be40256dddbf678e91bc56b96c36da07b56ec8bda29d` | 8.0.4 |
| bakta | `oschwengers/bakta@sha256:1e24ae2f435c46f876a322298ab96d994d15aadc2824e1f16a77fa9a1bc07f18` | v1.12.2 |
| blast | `quay.io/biocontainers/blast@sha256:9dfc69f990c0aeb936276ee591ed32919a79f46dfa34060055e05a050a17959c` | 2.16.0--h66d330f_5 |
| genomad | `quay.io/biocontainers/genomad@sha256:d54652ff868af72604ea0e54c541c5c92168c7ac821a422de05dcbad65509e5c` | 1.12.0--pyhdfd78af_0 |
| miniprot | `quay.io/biocontainers/miniprot@sha256:2eb53fea53743b56b1c10ac4b202469b92db959ffbb8bafd83156fdbde80e22f` | 0.18--h577a1d6_0 |
| transtermhp | `quay.io/biocontainers/transtermhp@sha256:fc6fe0404be9340a883ec7a7b90415dfc3826c2bca7333a0314335c7bfee2d03` | 2.09--h9948957_2 |

**Apptainer.** The cache is one directory: `NXF_APPTAINER_CACHEDIR`, else `$GIPPA_ROOT/containers`,
else `DIR/containers` ([Install › Step 3](install.md#the-five-steps)). Inside it, each image
must have the file name Nextflow looks for: the reference with `/` and `:` replaced by `-`, plus
`.img`. For example:

```bash
CACHE=${NXF_APPTAINER_CACHEDIR:-$DB/containers}; mkdir -p "$CACHE"
apptainer pull "$CACHE/oschwengers-bakta@sha256-1e24ae2f435c46f876a322298ab96d994d15aadc2824e1f16a77fa9a1bc07f18.img" \
    docker://oschwengers/bakta@sha256:1e24ae2f435c46f876a322298ab96d994d15aadc2824e1f16a77fa9a1bc07f18
```

Preflight under `--profile apptainer` checks only that the file exists in the cache. It does not
contact the registry, so a cache copied from another machine works offline. To record the images
in `db_manifest.json`, run `gippa setup --db-dir "$DB" --images all --profile apptainer --no-uniop`
(add the opt-out, such as `--no-antismash`, of each database you leave out; the others are
checked, and adopted when present). Setup skips the pull for every image file that is already
there. It does not re-hash an image file, so pull only by digest, as above.

**Docker.** `docker pull <reference>` with the reference from the table. `gippa setup --profile
docker` (images are on by default with a runtime profile) always runs `docker pull` and checks that the pulled digest equals the pin.
Nextflow also refers to each image by digest. Both need the registry, so for compute nodes
without network, use Apptainer with a copied cache.

## Download on a machine with internet, copy to the cluster

The simplest route needs no manual unpacking: run `gippa setup` itself where the network is, then
copy the DB directory. With `NXF_APPTAINER_CACHEDIR` unset, the image cache goes to
`DIR/containers` and is copied along with the databases.

```bash
# on the machine with internet access
unset NXF_APPTAINER_CACHEDIR GIPPA_ROOT
gippa setup --db-dir /scratch/gippa-db --profile apptainer   # every item and image is on by default

# copy (resumable; -a keeps the symlinks and permissions inside the Bakta DB)
rsync -a --partial --info=progress2 /scratch/gippa-db/ cluster:/path/to/db/

# on the cluster
gippa setup --db-dir /path/to/db --check
GIPPA_DB_DIR=/path/to/db gippa preflight --genome genome.fna --peptides peptides.tsv \
    --profile apptainer
```

The databases are recorded by relative paths, so `--check` passes in the new place. The
`image:<tool>` records keep the old machine's cache path; preflight looks the images up by file
name in the cluster's cache. A local UniOP build is tied to the machine and paths it was built on;
build it on the cluster side, as described in [Install › UniOP](install.md#uniop).

## Shared copies at your site

A copy your site already holds counts if it verifies as above. Bakta and geNOMAD copies need
their version files, Promotech needs the zip next to the model, and single files need the pinned
SHA-256. If a copy does not verify, setup downloads its own (and, on a node with no network,
exits 2). GIPPA uses no stand-alone Pfam or other site database: the Bakta and antiSMASH
databases bring their own.

All items live under one DB directory, and setup must be able to write there. To reuse a
read-only shared copy, make your own DB directory and link the item into it, for example
`ln -s /shared/bakta/db-light "$DB/bakta/db-light"`. Inside the containers, only the DB directory
is mounted (`-B DIR:DIR:ro` under Apptainer, `-v DIR:DIR:ro` under Docker; `workflow/nextflow.config`).
A link that points outside it resolves only if the target path is also visible in the container,
for example through your site's Apptainer bind paths. If in doubt, copy the files instead.

## Licences

The terms of each source still apply when you download it yourself. Read
[Licences](../user-guide/licences.md) first. Before any download, `gippa setup` prints the
geNOMAD and KEGG terms when those items are planned. In short:

- **geNOMAD** (database and tool, module `mobile`) is for non-commercial use only. For commercial
  use, skip it (`gippa setup --no-genomad`) and run with `--modules core,bgc,function`.
  `gippa setup --minimal` leaves out antiSMASH too; a later run on that DB directory then
  defaults to `--modules core,function`.
- **KEGG** is free for academic use only; other users need a licence from KEGG, or leave it out
  with `gippa setup --no-kegg`.
- **The contaminant library and UniOP** have no licence file. GIPPA downloads them from their
  authors and never redistributes them; whether your use is permitted is yours to decide.
- **The GPL and AGPL tools** (Bakta, Promotech and others; antiSMASH) run as separate programs in
  their own images or environments. GIPPA never imports their code.
- The Bakta and geNOMAD databases come with the terms on their Zenodo records (14916843, 14886553).
