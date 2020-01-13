FROM nfcore/base
LABEL authors="Alexander Peltzer <alex.peltzer@gmail.com>" \
      description="Docker image containing all requirements for nf-core/smrnaseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-smrnaseq-1.0.0/bin:$PATH
