nextflow.enable.types = true

//
// Build the reference used for alignment-based UMI deduplication (`--umi_dedup_method fgumi`):
// one "window" per distinct mature miRNA sequence, i.e. the mature sequence plus `flank` nucleotides
// (default 5) on each side, taken from the hairpin.
//
// For each distinct mature sequence (compared in upper case, U as T), in the order of the mature FASTA:
//   - the first hairpin of the hairpin FASTA that contains the sequence is used (first occurrence);
//   - the window is cut from that hairpin and shortened at the ends of the hairpin;
//   - a mature sequence that is in no hairpin is written without flanks.
// Mature miRNAs with the same sequence (paralogs such as hsa-let-7a-5p from let-7a-1/-2/-3) get one
// window, named after the first of them, and identical windows are written once. Reads from paralogs
// therefore align to one record and are deduplicated together.
//
process UMI_DEDUP_WINDOWS {
    label 'process_single'
    tag "$mature"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6b/6b244720eef0bd28a41d4f26e33d3800e75d9fc87f080e81d42d5d676b4960dc/data' :
        'community.wave.seqera.io/library/gawk:5.3.0--180f75ae8b0ce739' }"

    input:
    tuple(meta: Map, mature: Path)
    tuple(meta2: Map, hairpin: Path)

    output:
    fasta   = tuple(meta, file("${task.ext.prefix ?: 'umi_dedup_windows'}.fa"))
    windows = tuple(meta, file("${task.ext.prefix ?: 'umi_dedup_windows'}.tsv"))

    topic:
    tuple(task.process, 'gawk', eval("awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//'")) >> 'versions'

    script:
    // Number of flanking nucleotides on each side of the mature sequence
    def args   = task.ext.args ?: '-v flank=5'
    def prefix = task.ext.prefix ?: 'umi_dedup_windows'
    """
    awk ${args} -v fasta='${prefix}.fa' -v table='${prefix}.tsv' '
        function add_record(    s, len, p, key, n, c, cand, hit, h, start, end, w) {
            if (name == "") return
            s = toupper(seq); gsub(/U/, "T", s); len = length(s)
            if (pass == 1) {
                # Hairpin: index every k-mer by hairpin and position
                nh++; hairpin_name[nh] = name; hairpin_seq[nh] = s
                for (p = 1; p <= len - k + 1; p++) { key = substr(s, p, k); kmers[key] = kmers[key] " " nh ":" p }
            } else if (len > 0) {
                if (s in window_of) {
                    w = window_of[s]
                } else {
                    # First hairpin (and first position in it) that contains the mature sequence
                    h = 0; p = 0
                    if (len >= k) {
                        n = split(kmers[substr(s, 1, k)], cand, " ")
                        for (c = 1; c <= n && !h; c++) {
                            split(cand[c], hit, ":")
                            if (substr(hairpin_seq[hit[1]], hit[2], len) == s) { h = hit[1]; p = hit[2] }
                        }
                    } else {
                        for (c = 1; c <= nh && !h; c++) { p = index(hairpin_seq[c], s); if (p) h = c }
                    }
                    if (h) {
                        start = p - flank < 1 ? 1 : p - flank
                        end   = p + len - 1 + flank > length(hairpin_seq[h]) ? length(hairpin_seq[h]) : p + len - 1 + flank
                        w = substr(hairpin_seq[h], start, end - start + 1)
                    } else {
                        w = s; start = "NA"; end = "NA"
                    }
                    window_of[s] = w
                    if (!(w in window_names)) {
                        nw++; window_order[nw] = w
                        id = name; while (id in used) id = id "_"; used[id] = 1
                        window_id[w] = id
                        window_hairpin[w] = h ? hairpin_name[h] : "NA"; window_start[w] = start; window_end[w] = end
                    }
                }
                window_names[w] = (w in window_names) ? window_names[w] "," name : name
            }
            name = ""; seq = ""
        }
        BEGIN { k = 12; flank += 0 }
        FNR == 1 { add_record(); pass++ }
        /^>/ { add_record(); name = substr(\$0, 2); sub(/[ \\t\\r].*/, "", name); next }
        { gsub(/[ \\t\\r]/, ""); seq = seq \$0 }
        END {
            add_record()
            print "window\\tmature\\thairpin\\tstart\\tend\\tlength" > table
            for (i = 1; i <= nw; i++) {
                w = window_order[i]
                print ">" window_id[w] "\\n" w > fasta
                print window_id[w] "\\t" window_names[w] "\\t" window_hairpin[w] "\\t" window_start[w] "\\t" window_end[w] "\\t" length(w) > table
            }
            if (!nw) { print "No mature miRNA sequences found in ${mature}" > "/dev/stderr"; exit 1 }
        }
    ' ${hairpin} ${mature}
    """

    stub:
    def prefix = task.ext.prefix ?: 'umi_dedup_windows'
    """
    touch ${prefix}.fa
    printf 'window\\tmature\\thairpin\\tstart\\tend\\tlength\\n' > ${prefix}.tsv
    """
}
