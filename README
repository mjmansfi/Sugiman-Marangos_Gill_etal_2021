# Description
This repository contains all scripts and files needed to reproduce the bioinformatic analysis presented in the forthcoming publication by Sugiman-Marangos, Gill et al. 2021.

## Requirements
In order to run the scripts in this repository, you will need:

* [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html) or [Docker](https://docs.docker.com/get-docker/)
* git
* Patience.

`git` is used to pull this repository.

Singularity or Docker are used in order to run required software in a virtualized environment, making it portable and accessible to different devices.

Patience is needed because the scripts run many combinations of multiple sequence alignments and phylogenetic algorithms, which are computationally demanding and will take time even with a powerful computer.

## Installation
First, clone this repository.
```
git clone https://github.com/mjmansfi/Sugiman-Marangos_Gill_etal_2021
cd Sugiman-Marangos_Gill_etal_2021
```

The `Dockerfile` and `environment.yml` contain the instructions needed to build the virtualized environment. You can either build the image yourself:
```
docker built -t [build tag] .
```

where [build tag] is the desired name of the image. Typically this would be of the form `[Docker Hub username]/[image name]:[version number]`.

Alternatively, you can use my pre-built image by pulling it from Docker Hub:
```
docker pull mjmansfi/dt:0.4.0
```

If this image has been pruned from the Docker Hub, ping me in an issue and I will rebuild it.

There is currently no straightforward way to transmute a `Dockerfile` into a Singularity image directly. So, in order to build the Singularity image, pull it from the Docker Hub as follows:
```
# Note that the scripts expect to find the Singularity image in the ./lib folder.
cd ./lib
singularity pull docker://mjmansfi/dt:0.4.0
cd ../
```
Voila!
