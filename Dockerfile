FROM nfcore/base:1.14
LABEL authors="Phil Ewels <phil.ewels@scilifelab.se>, Chuan Wang <chuan.wang@scilifelab.se>, Rickard Hammarén <rickard.hammaren@scilifelab.se>, Lorena Pantano <lorena.pantano@gmail.com>" \
      description="Docker image containing all software requirements for the nf-core/smrnaseq pipeline"

# Install libtbb2 package for bowtie
RUN apt-get update && apt-get install libtbb2 -y

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-smrnaseq-1.1.0/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-smrnaseq-1.1.0 > nf-core-smrnaseq-1.1.0.yml
