Automated pipeline to create and run slabs of chromatin fibres of the chosen NRL.

## Single Fiber Replica Exchange

Single fibre creation is set up in 'inputs'. To create a minimal representation of a single fibre:

1. Modify linker_lenghts to the desired NRLs

2. Run 
```
python make_variable_nrl_N_from_1kx5.py
```

This creates a data/txt file and a DNA_sequence.txt file corresponding to the desired chromatin fibre. 

3. After creating a fibre, we move the files to repex, and we generate via replica exchange a series of snapshots for a single chromatin fibre in order to speed up sampling in a condensate bulk and avoid arrested configurations.

To do so, simply run 
```
bash initial.sh
```

This creates subfolders of the desired simulation salts and sets up the replica exchange. 

4. To run the simulation, submit tremd.in in lammps, in your desired HPC architecture, or run, for each salt

```
lmp -partition 16x1 -in tremd.in
```

5. Finally, to obtain representative snapshots of relaxed chromatin fibres, demix the trajectories by running

```
bash demix.sh  
```

and extract frames by 

```
bash create_data.sh
```
which ensures the frames have adequate file structure to be parsed as inputs for the slab conformation.


## Direct Coexistence Simulation

To produce a full-phase diagram:

1. Run ìn the parent directly the following command to move the relevant data into the desired folder in slabs

```
bash move_data.sh
```
2. Run, in slabs,
```
bash create_simulation.sh
```
and follow the promps: simulation name, number of fibers you desire and number of salts you want to sample. Not that these salts have to have sampling in repex before this step. 

3. Inside the simulation folder, run
```
bash submit_all_create.sh
```
This creates the salbs by compressing. Depending on your HPC architecture you may need to modify submission scripts - in this folder,r the working example is for Sulis HPC.
4. Finally, after the creation is over,  run
```
bash submit_all_continue.sh
```
which produces the final trajectory.



