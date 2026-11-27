class STRAND {
    static void assertAlignments(workDir, outputDir, forwardOnly) {
        def forward = ['forward_A|reference_A|forward', 'forward_B|reference_B|forward']
        def both = forward + [
            'forward_A|reference_B|reverse',
            'forward_B|reference_A|reverse',
            'reverse_only|reference_C|reverse',
        ]
        def aliases = [
            mature: 'MIRNA_QUANT:BOWTIE_MAP_MATURE',
            hairpin: 'MIRNA_QUANT:BOWTIE_MAP_HAIRPIN',
            seqcluster: 'MIRNA_QUANT:BOWTIE_MAP_SEQCLUSTER',
            genome: 'GENOME_QUANT:BOWTIE_MAP_GENOME',
        ]
        aliases.each { branch, alias ->
            def scripts = []
            new File("${workDir}").eachFileRecurse(groovy.io.FileType.FILES) { file ->
                if (file.name == '.command.run' && file.text.contains("NFCORE_SMRNASEQ:${alias}")) {
                    scripts.add(new File(file.parentFile, '.command.sh').text)
                }
            }
            assert scripts.size() == 1
            def restricted = forwardOnly && branch != 'genome'
            assert scripts[0].contains('--norc') == restricted
            assert !scripts[0].contains('--nofw')

            def observed = new File("${outputDir}/${branch}.tsv").readLines().collect { line ->
                def fields = line.split('\t')
                def orientation = (fields[1].toInteger() & 16) ? 'reverse' : 'forward'
                "${fields[0]}|${fields[2]}|${orientation}".toString()
            }
            assert observed.sort() == (restricted ? forward : both).sort()
        }
    }
}
