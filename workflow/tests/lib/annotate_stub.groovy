// Helpers for annotate.nf.test (M2-07). nf-test 0.9 has no top-level functions, so the test
// loads this file with GroovyShell and gets a map of closures.
// Trace names look like "ANNOTATE:BAKTA (sample)" or "GIAP_STEP (parse_bakta)".

def procOf = { String name -> name.split(' \\(')[0].tokenize(':')[-1] }
def tagOf = { String name ->
    def m = name =~ /\((.*)\)\s*$/
    m.find() ? m.group(1) : null
}
def names = { trace -> trace.succeeded().collect { it.name.toString() } }
def count = { trace, String proc -> names(trace).count { procOf(it) == proc } }

// A py step ran: GIAP_STEP tagged with the step name, or a GIAP_STEP alias named after the step.
def stepCount = { trace, String step ->
    names(trace).count { n ->
        def p = procOf(n)
        (p == 'GIAP_STEP' && tagOf(n) == step) || p == step.toUpperCase()
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
// Task directories of a py step (GIAP_STEP tag or alias) or of a tool process.
def stepDirs = { String work, String step ->
    taskDirs(work).findAll { d ->
        def n = nameLine(d).replaceFirst(/^### name:\s*'/, '').replaceFirst(/'\s*$/, '')
        def p = procOf(n)
        (p == 'GIAP_STEP' && tagOf(n) == step) || p == step.toUpperCase()
    }
}
def procDirs = { String work, String proc ->
    taskDirs(work).findAll { d ->
        def n = nameLine(d).replaceFirst(/^### name:\s*'/, '').replaceFirst(/'\s*$/, '')
        procOf(n) == proc
    }
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

return [defaultTools: defaultTools, pySteps: pySteps, procOf: procOf, tagOf: tagOf, names: names, count: count, stepCount: stepCount,
        taskDirs: taskDirs, stepDirs: stepDirs, procDirs: procDirs, staged: staged, files: files]
