# Troubleshooting

## Generic problems

### Input files not found

If only no file, only one input file , or only read one and not read two is picked up then something is wrong with your input file declaration

1. The path must be enclosed in quotes (`'` or `"`)
2. The path must have at least one `*` wildcard character. This is even if you are only running one paired end sample.

If the pipeline can't find your files then you will get the following error

```
ERROR ~ Cannot find any reads matching: *.fastq.gz
```

Note that if your sample name is "messy" then you have to be very particular with your glob specification. A file name like `L1-1-D-2h_S1_L002_R1_001.fastq.gz` can be difficult enough for a human to read.

The above information is also covered in the [usage README](usage.md#--reads).


### Data organization
The pipeline can't take a list of multiple input files - it takes a glob expression. If your fastq files are scattered in different paths then we recommend that you generate a directory with symlinked files. If running in paired end mode please make sure that your files are sensibly named so that they can be properly paired. See the previous point.

## Extra resources and getting help
If you still have an issue with running the pipeline then feel free to contact us.
Have look at the [issue tracker for our repo](https://github.com/nf-core/smrnaseq/issues). Maybe someone has already had the same problem?

Gitter is a chatt client connected to Github, feel free to come in and chat with us;
[nfcore/smrnaseq Gitter]((https://gitter.im/nf-core/smrnaseq))

If you have problems that are related to Nextflow and not our pipeline then check out the [Nextflow gitter channel](https://gitter.im/nextflow-io/nextflow) or the [google group](https://groups.google.com/forum/#!forum/nextflow).
