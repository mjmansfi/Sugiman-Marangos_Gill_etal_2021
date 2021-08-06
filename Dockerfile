FROM continuumio/miniconda:4.7.12

LABEL authors="michaeljamesmansfield@gmail.com" \
	description="Docker file containing several tools for phylogenetic analysis"

COPY VERSION .

RUN apt-get update \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*
	
COPY environment.yml /

RUN conda env create -f /environment.yml && \
	conda clean -a

ENV PATH /opt/conda/envs/dt/bin:$PATH


