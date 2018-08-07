From:nfcore/base
Bootstrap:docker

%labels
    MAINTAINER Alexander Peltzer <alex.peltzer@gmail.com>
    DESCRIPTION Container image containing all requirements for the nf-core/smrnaseq pipeline
    VERSION 1.0

%files
    environment.yml /

%post
    /opt/conda/bin/conda env update -n root -f /environment.yml
    /opt/conda/bin/conda clean -a
