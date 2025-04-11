import numpy as np
import matplotlib
matplotlib.use('Agg')  # Set the backend before importing pyplot
import matplotlib.pyplot as plt
import seaborn as sns
import os
from ovito.io import import_file
from ovito.modifiers import CreateBondsModifier
from ovito.data import *

# Adjust these paths and parameters as necessary
folders = ['172', '173', '174', '175', '176', '177']
data_directory = './'  # Assuming the script is run from a location directly above the folders

def compute_second_neighbor_contacts(data, fibers_in_middle, molecules_per_fiber):
    fiber_contacts = []

    # Iterate over fibers in the middle
    for fiber_id in fibers_in_middle:
        # Determine molecule range for the current fiber
        start_molecule = (fiber_id - 1) * molecules_per_fiber + 1
        end_molecule = fiber_id * molecules_per_fiber

        # Initialize contact count for the current fiber
        fiber_contact_count = 0

        # Iterate over molecule IDs in the current fiber
        for molecule_id in range(start_molecule, end_molecule + 1):
            # Check if the molecule has a bond with its second neighbor
            if has_bond_with_second_neighbor(data, molecule_id):
                # Increment contact count if the second neighbor is bonded
                fiber_contact_count += 1

        # Add contact count for the current fiber to the list
        fiber_contacts.append(fiber_contact_count)

    return fiber_contacts

def has_bond_with_second_neighbor(data, molecule_id):
    # Find the IDs of the two neighbors
    neighbor_id_1 = molecule_id + 2
    #neighbor_id_2 = molecule_id - 2

    # Check if either neighbor has a bond with the molecule
    return is_bonded(data, molecule_id, neighbor_id_1) or is_bonded(data, molecule_id, neighbor_id_2)

def is_bonded(data, molecule_id_1, molecule_id_2):
    # Iterate over bonds in the data
    for bond in data.particles.bonds:
        # Check if the bond involves the two molecules
        if (bond.particle1.index == molecule_id_1 and bond.particle2.index == molecule_id_2) or \
           (bond.particle1.index == molecule_id_2 and bond.particle2.index == molecule_id_1):
            return True  # Return True if the bond is found

    return False  # Return False if no bond is found

def process_folder(folder_path, molecules_per_fiber=12):
    all_fiber_contacts = []
    for filename in sorted(os.listdir(folder_path)):
        if filename.endswith('dna.dump'):
            file_path = os.path.join(folder_path, filename)
            pipeline = import_file(file_path)
            pipeline.modifiers.append(CreateBondsModifier(cutoff=50))
            
            # Get the number of timesteps from the data file
            num_timesteps = pipeline.source.num_frames
    
            # Loop over each timestep
            for frame in range(1):
                data = pipeline.compute(frame*10)

                z_coords = data.particles['Position'][:, 2]
                average_z = np.mean(z_coords)
                std_dev_z = np.std(z_coords)
                middle_start = average_z - std_dev_z / 3
                middle_end = average_z + std_dev_z / 3

                fibers_in_middle = set()
                molecule_ids = data.particles['Molecule Identifier'].array
                for molecule_id in np.unique(molecule_ids):
                    if molecule_id == 0:
                        continue
                    positions = data.particles.positions[data.particles['Molecule Identifier'].array == molecule_id]
                    z_positions = positions[:, 2]
                    if any((z_positions >= middle_start) & (z_positions <= middle_end)):
                        fiber_id = (molecule_id - 1) // molecules_per_fiber
                        fibers_in_middle.add(fiber_id)
                print('Fibers in the middle:')
                print(len(fibers_in_middle))
                # Compute second neighbor contacts for fibers in the middle
                fiber_contacts = compute_second_neighbor_contacts(data, fibers_in_middle, molecules_per_fiber)
                all_fiber_contacts.extend(fiber_contacts)
    return all_fiber_contacts


# Process each folder and collect second neighbor contact distributions
second_neighbor_contact_distributions = [process_folder(os.path.join(data_directory, folder)) for folder in folders]

# Plotting
sns.set(style="whitegrid")
plt.figure(figsize=(12, 8))
sns.violinplot(data=second_neighbor_contact_distributions)
plt.xticks(range(len(folders)), folders)
plt.title('Distribution of Second Neighbor Contacts for Folders 172 to 177')
plt.xlabel('Folder')
plt.ylabel('Number of Contacts')
plt.savefig('Second_neighbor_contacts.png')

