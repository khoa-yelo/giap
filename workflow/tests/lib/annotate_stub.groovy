// Helpers for annotate.nf.test (M2-07). nf-test 0.9 has no top-level functions, so the test
// loads this file with GroovyShell and gets a map of closures.
// Trace names look like "ANNOTATE:BAKTA (sample)" or "GIPPA_STEP (parse_bakta)".

def procOf = { String name -> name.split(' \\(')[0].tokenize(':')[-1] }
def tagOf = { String name ->
    def m = name =~ /\((.*)\)\s*$/
    m.find() ? m.group(1) : null
}
def names = { trace -> trace.succeeded().collect { it.name.toString() } }
def count = { trace, String proc -> names(trace).count { procOf(it) == proc } }

// A py step ran: GIPPA_STEP tagged with the step name, or a GIPPA_STEP alias named after the step.
def stepCount = { trace, String step ->
    names(trace).count { n ->
        def p = procOf(n)
        (p == 'GIPPA_STEP' && tagOf(n) == step) || p == step.toUpperCase()
    }
}

def taskDirs = { String work ->
    def out = []
    new File(work).eachDir { a -> a.eachDir { t -> if (new File(t, '.command.run').isFile()) out << t } }
    out
}
def nameLine = { File dir ->
    new File(dir, '.command.run').readLines().find { it.startsWith('### name:') } ?: ''
}
// Task directories of a py step (GIPPA_STEP tag or alias) or of a tool process.
def stepDirs = { String work, String step ->
    taskDirs(work).findAll { d ->
        def n = nameLine(d).replaceFirst(/^### name:\s*['"]/, '').replaceFirst(/['"]\s*$/, '')
        def p = procOf(n)
        (p == 'GIPPA_STEP' && tagOf(n) == step) || p == step.toUpperCase()
    }
}
def procDirs = { String work, String proc ->
    taskDirs(work).findAll { d ->
        def n = nameLine(d).replaceFirst(/^### name:\s*['"]/, '').replaceFirst(/['"]\s*$/, '')
        procOf(n) == proc
    }
}
// Ordered GIPPA_STEP inputs as written by its stub to .gippa_inputs ("key=path" per line). Staging
// (`ln -s`) order in .command.run is Nextflow's, not the channel order, so it cannot show input order.
def stepInputs = { File dir ->
    def f = new File(dir, '.gippa_inputs')
    f.isFile() ? f.readLines().findAll { it }.collect { it.split('=', 2) as List } : []
}

// Staging lines of .command.run in order: [[source, target], ...]
def staged = { File dir ->
    new File(dir, '.command.run').readLines().collectMany { line ->
        (line =~ /ln -s\S*\s+(\S+)\s+(\S+)/).collect { m -> [m[1], m[2]] }
    }
}
def files = { File dir ->
    def out = []
    dir.eachFileRecurse { f -> if (f.isFile() && !f.name.startsWith('.')) out << f }
    out
}

// M2-index "Annotation step contract": external tools on by default, and every py step.
def defaultTools = ['BAKTA', 'UNIOP', 'PROMOTECH', 'PROMOTER_CALCULATOR', 'TERMINATORNET',
                    'TRANSTERMHP', 'OSTIR', 'GENOMAD', 'ANTISMASH']
def pySteps = ['parse_bakta', 'uniop_input', 'operon_pairs', 'promoter_windows', 'promoter_layer',
               'tx_layer', 'transterm_input', 'terminator_layer', 'ostir_input', 'rbs_layer',
               'codon_usage', 'region_layer', 'context_layers']

// M2-07 review decision 1: the --in keys of every GIPPA_STEP task whose step is in the registry
// fixture (data/annotate/step_registry.json = gippa.steps registry, kept equal by
// tests/unit/test_nf_annotate_registry.py) are declared inputs, and every required input is given.
def registryProblems = { String work, File regFile ->
    def reg = new groovy.json.JsonSlurper().parse(regFile)
    def problems = []
    taskDirs(work).each { d ->
        def n = nameLine(d).replaceFirst(/^### name:\s*['"]/, '').replaceFirst(/['"]\s*$/, '')
        def step = tagOf(n)
        if (step == null || !reg.containsKey(step)) return
        def spec = reg[step]
        def keys = stepInputs(d).collect { it[0] }.unique()
        def allowed = spec.inputs + spec.optional_inputs
        def extra = keys.findAll { !(it in allowed) }
        def missing = spec.inputs.findAll { !(it in keys) }
        if (extra) problems << "${step}: undeclared input key(s) ${extra}"
        if (missing) problems << "${step}: required input(s) ${missing} not given"
    }
    problems
}

return [registryProblems: registryProblems, defaultTools: defaultTools, pySteps: pySteps, procOf: procOf, tagOf: tagOf, names: names, count: count, stepCount: stepCount,
        taskDirs: taskDirs, stepDirs: stepDirs, procDirs: procDirs, staged: staged, stepInputs: stepInputs, files: files]
