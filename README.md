# NucleosomeLinker
Source Code for the paper 'Nucleosome Spacing Can Fine-Tune Higher Order Chromatin Assembly'

We are delighted to share our model for chromatin phase separation and analysis scripts with the community. Please use it freely and cite our paper: preprint https://www.biorxiv.org/content/10.1101/2024.12.23.627571v1. 

## System requirements
Linux with C++ compilers with MPI. Tested on: CSD3 peta-4 cluster (https://www.hpc.cam.ac.uk/systems/peta-4) with Intel 2017 compliers, Archer HPC (Cray system with GNU and Cray compilers) and Sulis Tier 2 HPC.

## Installation Guide
Please follow the instructions for building our custom LAMMPS code from the multiscale model repository: https://github.com/CollepardoLab/CollepardoLab_Chromatin_Model

## Demo
To run a demo for a slab system:

1. Move to the "Demo" directory

2. run with lammps

| mpirun -np 1 ./lmp_DNA_mpi -in in.run

It will produce two LAMMPS trajectory file "dna.dump" and "cores.dump", a reduced version without the DNA. These can be viewed in Ovito (https://www.ovito.org/), VMD or Pymol. 

## Full simulation setup

These simulations can be used to reproduce the results from this paper. To fully reproduce phase diagrams one needs to vary the parameters E1 and A in the input script. The correct mapping can be found in salt mapping, but this mapping is automated via bash scripting.


## Analysis

Under Analysis, all custom python scripts to analyze the trayectories can be found. 
