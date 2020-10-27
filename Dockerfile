FROM nfcore/base:1.11
LABEL authors="Phil Ewels <phil.ewels@scilifelab.se>, Chuan Wang <chuan.wang@scilifelab.se>, Rickard Hammarén <rickard.hammaren@scilifelab.se>, Lorena Pantano <lorena.pantano@gmail.com>" \
      description="Docker image containing all software requirements for the nf-core/smrnaseq pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-smrnaseq-1.0.1dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-smrnaseq-1.0.1dev > nf-core-smrnaseq-1.0.1dev.yml

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron
