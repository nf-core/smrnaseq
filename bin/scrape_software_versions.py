#!/usr/bin/env python
from __future__ import print_function
import os

<<<<<<< HEAD
results = {}
version_files = [x for x in os.listdir(".") if x.endswith(".version.txt")]
for version_file in version_files:
=======
regexes = {
    "nf-core/smrnaseq": ["v_pipeline.txt", r"(\S+)"],
    "R": ["v_R.txt", r"R version (\S+)"],
    "Nextflow": ["v_nextflow.txt", r"(\S+)"],
    "FastQC": ["v_fastqc.txt", r"FastQC v(\S+)"],
    "Trim Galore!": ["v_trim_galore.txt", r"version (\S+)"],
    "Bowtie": ["v_bowtie.txt", r"version (\S+)"],
    "Samtools": ["v_samtools.txt", r"samtools (\S+)"],
    "Htseq": ["v_htseq.txt", r"version (\S+)"],
    "FASTX": ["v_fastx.txt", r"Toolkit (\S+)"],
    "miRTrace": ["v_mirtrace.txt", r"(\S+)"],
    "MultiQC": ["v_multiqc.txt", r"multiqc, version (\S+)"],
    "miRDeep2": ["v_mirdeep2.txt", r"miRDeep(\S+)"],
}
results = OrderedDict()
results["nf-core/smrnaseq"] = '<span style="color:#999999;">N/A</span>'
results["Nextflow"] = '<span style="color:#999999;">N/A</span>'
results["R"] = '<span style="color:#999999;">N/A</span>'
results["FastQC"] = '<span style="color:#999999;">N/A</span>'
results["Trim Galore!"] = '<span style="color:#999999;">N/A</span>'
results["Bowtie"] = '<span style="color:#999999;">N/A</span>'
results["Samtools"] = '<span style="color:#999999;">N/A</span>'
results["Htseq"] = '<span style="color:#999999;">N/A</span>'
results["FASTX"] = '<span style="color:#999999;">N/A</span>'
results["miRTrace"] = '<span style="color:#999999;">N/A</span>'
results["MultiQC"] = '<span style="color:#999999;">N/A</span>'
results["miRDeep2"] = '<span style="color:#999999;">N/A</span>'
>>>>>>> origin/dev

    software = version_file.replace(".version.txt", "")
    if software == "pipeline":
        software = "nf-core/smrnaseq"

    with open(version_file) as fin:
        version = fin.read().strip()
    results[software] = version

# Dump to YAML
print(
    """
id: 'software_versions'
section_name: 'nf-core/smrnaseq Software Versions'
section_href: 'https://github.com/nf-core/smrnaseq'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
"""
)
for k, v in sorted(results.items()):
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k, v))
print("    </dl>")

# Write out as tsv file:
with open("software_versions.tsv", "w") as f:
    for k, v in sorted(results.items()):
        f.write("{}\t{}\n".format(k, v))
