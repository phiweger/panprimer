nextflow.preview.dsl = 2


process create_graph {
    // echo true
    container 'nanozoo/bifrost:1.0.4--7fae7c4'
    publishDir params.outdir, mode: 'copy', overwrite: true
    cpus 8

    input:
        file(genomes)

    output:
        file('pangenome.gfa')
        file('pangenome.bfg_colors')

    """    
    Bifrost build -v -c -r ${genomes} -t ${task.cpus} -k ${params.k_dbg} -i -d -o pangenome
    """
}


process update_graph {
    // echo true
    container 'nanozoo/bifrost:1.0.4--7fae7c4'
    publishDir params.outdir, mode: 'copy', overwrite: true
    cpus 8

    input:
        file(genomes)
        file(graph)
        file(colors)
        
    output:
        file('pangenome.gfa')
        file('pangenome.bfg_colors')

    """
    Bifrost update -v -g ${graph} -f ${colors} -r ${genomes} -t ${task.cpus} -k ${params.k_dbg} -i -d -o pangenome.update
    
    mv pangenome.update.gfa pangenome.gfa
    mv pangenome.update.bfg_colors pangenome.bfg_colors
    """
}


process annotate_colors {
    container 'nanozoo/blastfrost:0.1--02b3166'
    cpus 8

    input:
        file(graph)
        file(colors)

    output:
        file('pangenome.gfa_colored.gfa')
        file('pangenome.gfa.index')

    """
    BlastFrost -t ${task.cpus} -c -g ${graph} -f ${colors} -o annotation
    """
}


process design_primers {
    echo true
    publishDir params.outdir, mode: 'copy', overwrite: true

    input:
        file(graph)
        file(index)
        file(reference)
        file(include)
        file(exclude)

    output:
        file('primers.csv')        

    shell:

    if( "${params.exclude}" )

    '''
    sed 's/^/1,/' !{include} > pattern.txt
    sed 's/^/0,/' !{exclude} >> pattern.txt

    design_primers.py \
        --graph !{graph} --index !{index} --pattern pattern.txt \
        --minlen !{params.minlen}  --maxlen !{params.maxlen} --maxpairs 100 \
        -n !{params.n_primers} --outfile primers.csv --genome !{reference}
    '''

    else

    '''
    sed 's/^/1,/' !{include} > pattern.txt

    design_primers.py \
        --graph !{graph} --index !{index} --pattern pattern.txt \
        --minlen !{params.minlen}  --maxlen !{params.maxlen} --maxpairs 100 \
        -n !{params.n_primers} --outfile primers.csv --genome !{reference}
    '''
}


process sketch {
    input:
        tuple(val(name), file(genome))

    output:
        tuple(val(name), file(genome), file("${name}.sig"))

    """
    sourmash compute -k ${params.k_sketch} -n ${params.n_sketch} -o ${name}.sig ${genome}
    """
}


process index {
    // echo true

    input:
        file(sketches)

    output:
        file("index.sbt.zip")

    """
    sourmash index index.sbt.zip ${sketches}
    """
}


// Ignore a genome if it is similar to a given index ("db")
process search {
    // echo true
    // errorStrategy 'ignore'

    input:
        tuple(val(name), file(genome), file(sketch))
        file(db)

    output:
        tuple(val(name), file(genome), file("${name}.csv")) optional true

    shell:
    '''
    sourmash search --threshold !{params.maxsim} -o !{name}.csv !{sketch} !{db}
    if [[ $(wc -l < !{name}.csv) -ge 2 ]]
    then
        rm !{name}.csv
    fi
    '''
}


process remove {
    input:
        file(genomes)
        file(remove)

    output:
        file('remaining.txt')

    """
    remove_genomes.py -i ${genomes} --remove ${remove} -o remaining.txt
    """

}


process include {
    // echo true

    input:
        file(include)

    output:
        file('genome_list')

    """
    find -L ${include} -type f -name "*.fasta" > genome_list
    """
}


workflow create_db_wf {
    take:
        genomes

    main:
        index(sketch(genomes).collect {it[2]} )

    emit:
        index.out
}


workflow filter_similarity_wf {
    take:
        genomes
        db

    main:
        search(sketch(genomes), db)
    
    emit:
        search.out
}


workflow primer_wf {
    take:
        genomes

    main:
        create_graph(genomes)
}


workflow {
    // TODO: allow multiple extensions: fa, fa.gz, fna, fna.gz, fasta, fasta.gz
    include_ch = Channel.fromPath(params.include)
                        .map { it -> tuple(it.baseName, it) }
    incl = include_ch.collectFile(name: 'include.txt', newLine: true) { it[1].toString() }

    // We need one random sample on which to orient nodes from the genome graph
    // that are putative targets for primers, such that the primers "match", ie
    // the resulting transcripts aren't synthesized in opposite directions etc.
    reference = include_ch.randomSample(1, 42).map {it[1]}

    // Specify genomes that we want to exclude spcifically
    if (params.exclude) {
        exclude_ch = Channel.fromPath(params.exclude)
                            .map { it -> tuple(it.baseName, it) }

        // If the exclude set contains genomes that are too similar to what we
        // search primers against, we won't find any. So filter those out.
        if (params.filter) {
            db = create_db_wf(include_ch)
            search(sketch(exclude_ch), db)
            excl = search.out.collectFile(name: 'exclude.txt', newLine: true) { it[1].toString() }
        }

        else {
            excl = exclude_ch.collectFile(name: 'exclude.txt', newLine: true) { it[1].toString() }
        }

        create_graph(excl)
        update_graph(incl, create_graph.out)  
        annotate_colors(update_graph.out)
        design_primers(
            annotate_colors.out, reference, incl, excl)
    }


    else {
        create_graph(incl)
        annotate_colors(create_graph.out)
        // Optional input:
        // https://gist.github.com/pditommaso/ff13c333f461ca0d9b839d9e3416b376
        design_primers(
            annotate_colors.out, reference, incl, 'NONE')
    }
}

