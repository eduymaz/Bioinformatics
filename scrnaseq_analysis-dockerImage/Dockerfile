#FROM alpine
FROM --platform=linux/amd64 continuumio/miniconda3

WORKDIR /app/conda
RUN conda update conda
COPY env.yaml .
RUN conda env create -f env.yaml
RUN echo "source activate myenv" > ~/.bashrc
ENV PATH=/opt/conda/envs/myenv/bin:$PATH
# RUN conda install bioconductor-dropletutils
# Create the environment:
COPY src /app/src
#WORKDIR /app/lib
#RUN wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
#RUN tar -vxzf sratoolkit.tar.gz
#ENV PATH=$PATH:$PWD/sratoolkit.3.1.0-ubuntu64/bin
#RUN conda install -y bioconda::fastqc
#RUN conda install -y -c bioconda kallisto
#RUN conda install -y bioconda::bustools
#RUN mamba install kallisto
#RUN conda install -y kb-python
#RUN conda install -y conda-forge::r-base
#RUN Rscript /app/Rconf.R


COPY index /app/index

ENV KALLISTO_INDEX="/app/index/humanKallisto_index.idx"
ENV KALLISTO_T2G="/app/index/t2g.txt"
ENV KALLISTO_TRANSCRIPTOME="/app/index/transcriptome.fa"
COPY real_data /app/data
RUN mkdir /app/out
WORKDIR /app/src
#CMD ["bash", "main.sh"]
