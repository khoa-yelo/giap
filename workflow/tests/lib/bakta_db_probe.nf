// M11-14 item 3b (14p): a probe named BAKTA, so every `withName: 'BAKTA'` selector applies to it
// (the pinned image from workflow/conf/containers.config, the ps shim, conda off, the per-task
// TMPDIR). Instead of annotating it checks, inside that container, what the real BAKTA needs from
// the docker profile, failing on the first miss: <db_dir>/version.json readable (the profile's
// read-only -v bind: the DB is a val, not an input file, so nothing else mounts it), `ps` for the
// task metrics, and TMPDIR passed in as the task's own tmp/. Used only by docker_db.nf.test.
process BAKTA {
    input:
    val(db)

    output:
    path 'probe.txt', emit: probe

    script:
    """
    test -r '${db}/version.json'
    test "\$TMPDIR" = "\$PWD/tmp" && test -d "\$TMPDIR"
    ps -e -o pid= -o ppid= | grep -q .
    {
        echo "uid=\$(id -u)"
        echo "gid=\$(id -g)"
        echo "version_json=\$(tr -d '\\n' < '${db}/version.json')"
        echo "bakta=\$(bakta --version 2>&1 | head -n 1)"
    } > probe.txt
    """
}
