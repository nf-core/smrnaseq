// Helpers to check that mirna_quant/mirna_ids.tsv lists every miRNA identifier used in the
// edgeR and mirtop outputs (#357). Identifiers are compared as "<type>:<name>".

class MirnaIds {

    // Rows of mirna_quant/mirna_ids.tsv (without header), as lists of fields
    static List<List<String>> table(String outdir) {
        return new File("${outdir}/mirna_quant/mirna_ids.tsv").readLines().drop(1).collect { line -> line.split('\t', -1) as List }
    }

    // Mature names in the mirtop count tables
    static Set<String> mirtopNames(String outdir) {
        def names = [] as Set
        names.addAll(new File("${outdir}/mirna_quant/mirtop/mirna.tsv").readLines().drop(1).collect { line -> line.split('\t')[0] })
        def joined = new File("${outdir}/mirna_quant/mirtop/joined_samples_mirtop.tsv").readLines()
        def column = joined[0].split('\t').findIndexOf { field -> field == 'miRNA' }
        names.addAll(joined.drop(1).collect { line -> line.split('\t')[column] })
        return names
    }

    // Names in the edgeR count table of the given type (mature or hairpin): samples are rows,
    // miRNAs are the columns of the header line
    static Set<String> edgerNames(String outdir, String type) {
        def header = new File("${outdir}/mirna_quant/edger_qc/${type}_counts.csv").readLines()[0]
        return header.split(',').drop(1).collect { field -> field.replaceAll('"', '') } as Set
    }

    // Identifiers used in the outputs that are missing from mirna_ids.tsv
    static List<String> missing(String outdir) {
        def known = table(outdir).collect { row -> "${row[2]}:${row[0]}".toString() } as Set
        def used = [] as Set
        used.addAll(mirtopNames(outdir).collect { name -> "mature:${name}".toString() })
        used.addAll(edgerNames(outdir, 'mature').collect { name -> "mature:${name}".toString() })
        used.addAll(edgerNames(outdir, 'hairpin').collect { name -> "hairpin:${name}".toString() })
        return (used - known).sort()
    }
}
