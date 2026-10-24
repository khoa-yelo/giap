// Shared checks for workflow/tests/*_support.nf.test (M9-17 follow-up: the plot support outputs
// under Nextflow). `published(outDir)`: every published file, OUT-relative, logs/ excluded.
// `supportDirs(outDir)`: the top-level OUT entries named *_support (a support key published as
// its own directory instead of under plots/support/). `argvOuts(dir)`: the --out values the fake
// giap (workflow/tests/bin_giap/giap) recorded in a task directory's giap.argv.
def published = { String outDir ->
    def root = new File(outDir)
    def out = []
    if (root.isDirectory()) root.eachFileRecurse { f ->
        if (f.isFile()) out << root.toPath().relativize(f.toPath()).toString()
    }
    out.findAll { !it.startsWith('logs/') }.sort()
}
def supportDirs = { String outDir ->
    def root = new File(outDir)
    root.isDirectory() ? root.listFiles().findAll { it.name.endsWith('_support') }.collect { it.name }.sort() : []
}
def argvOuts = { File dir ->
    def f = new File(dir, 'giap.argv')
    if (!f.isFile()) return []
    def lines = f.readLines()
    (0..<lines.size()).findAll { lines[it] == '--out' && it + 1 < lines.size() }.collect { lines[it + 1] }
}
return [published: published, supportDirs: supportDirs, argvOuts: argvOuts]
