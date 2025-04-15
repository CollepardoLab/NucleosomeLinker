# NucleosomeLinker
Source Code for the paper 'Nucleosome Spacing Can Fine-Tune Higher Order Chromatin Assembly'

We are delighted to share our model for chromatin phase separation and analysis scripts with the community. Please use it freely and cite our paper: preprint https://www.biorxiv.org/content/10.1101/2024.12.23.627571v1. 

For questions, contact Julia Maristany in mjm261@cam.ac.uk

## System requirements
- Linux operating system  
- C++ compiler with MPI support (e.g., Intel, GCC, or Cray)  
- MPI library (e.g., Intel MPI, OpenMPI)  
- Python 3.6+ for analysis scripts  
- Optional: Ovito, VMD, or PyMOL for trajectory visualization

### Tested on:
- [CSD3 Cluster](https://docs.hpc.cam.ac.uk/hpc/) with Intel 2017 compilers  
- [ARCHER](https://www.archer2.ac.uk/) HPC (Cray system with GNU and Cray compilers)  
- [Sulis Tier 2 HPC](https://warwick.ac.uk/research/rtp/sc/sulis/)

## Installation Guide
Please follow the instructions for building our custom LAMMPS code from the multiscale model repository: https://github.com/CollepardoLab/CollepardoLab_Chromatin_Model

## Demo
To run a demo for a slab system:

1. Move to the "Demo" directory

2. Run with lammps, on at least 16 cores

```
mpirun -np 16 ./lmp_DNA_mpi -in in.run
```

It will produce two LAMMPS trajectory files, "dna.dump" and "cores.dump", a reduced version without the DNA. These can be viewed in Ovito (https://www.ovito.org/), VMD or Pymol. 

## Full simulation setup

These simulations can be used to reproduce the results from this paper. To fully reproduce phase diagrams one needs to vary the parameters E1 and A in the input script. The correct mapping can be found in salt mapping, but this mapping is automated via bash scripting. More details are in the Simulation folder.


## Analysis

Under Analysis, all custom Python scripts to analyze the trajectories can be found.  More details are in the Analysis folder.
