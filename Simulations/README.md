Automated pipeline to create and run slabs of chromatin fibres of the chosen NRL.

## Single Fiber Replica Exchange

Single fibre creation is set up in 'inputs'. To create a minimal representation of a single fibre:

1 - Modify linker_lenghts to the desired NRLs
2 - Run 
```
python make_variable_nrl_N_from_1kx5.py
```

This creates a data/txt file and a DNA_sequence.txt file corresponding to the desired chromatin fibre. 

In 'sims', we generate via replica exchange a series of snapshots for a single chromatin fibre, in order to speed up sampling in a condensate bulk and avoid arrested configurations.

To do so, simply run 
```
bash initial.sh
```

This creates subfolders of the desired simulation salts and sets up the replica exchange. To run the simulation, submit tremd.in in lammps, in your desired HPC architecture. 


