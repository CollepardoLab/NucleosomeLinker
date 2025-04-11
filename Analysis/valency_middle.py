from ovito.io import import_file
from ovito.modifiers import CreateBondsModifier
from ovito.data import *
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Set the backend before importing pyplot
import matplotlib.pyplot as plt

# Load the dataset
pipeline = import_file('172/dna.dump')

# Apply the CreateBondsModifier with an appropriate cutoff
cutoff_distance = 40  # Adjusted cutoff distance
pipeline.modifiers.append(CreateBondsModifier(cutoff=cutoff_distance))

# Function to calculate valency with unique contacts between fibers
def calculate_valency(data, fibers_in_middle, molecules_per_fiber):
    # Access the particle properties, including molecule IDs
    molecule_ids = data.particles['Molecule Identifier'].array
    num_molecules = np.max(molecule_ids)  # Assuming molecule IDs are sequential starting from 1

    # Initialize contact counts for each fiber
    total_fibers = len(fibers_in_middle)
    fiber_contacts = np.zeros(total_fibers, dtype=int)

    # Access computed bonds and bond topology
    bond_topology = data.particles.bonds.topology

    # Dictionary to map molecule ID to its fiber, excluding molecule ID 0
    molecule_to_fiber = {mol_id: (mol_id - 1) // molecules_per_fiber for mol_id in range(1, num_molecules + 1) if mol_id > 0}

    # Adjust mapping for fibers in the middle to index mapping
    fiber_index_mapping = {fiber_id: index for index, fiber_id in enumerate(fibers_in_middle)}

    # Set to store unique contacts between fibers
    unique_contacts = set()

    # Iterate over each bond to identify contacts between fibers
    for bond in bond_topology:
        particle_a, particle_b = bond
        molecule_a, molecule_b = molecule_ids[particle_a], molecule_ids[particle_b]

        # Skip if either molecule ID is 0
        if molecule_a == 0 or molecule_b == 0:
            continue

        fiber_a, fiber_b = molecule_to_fiber.get(molecule_a), molecule_to_fiber.get(molecule_b)

        # Check if molecules belong to different fibers
        if fiber_a is not None and fiber_b is not None and fiber_a != fiber_b:
            # Ensure at least one of the fibers is in the middle
            if fiber_a in fibers_in_middle or fiber_b in fibers_in_middle:
                # Create a unique identifier for the contact
                contact_id = tuple(sorted((fiber_a, fiber_b)))
                if contact_id not in unique_contacts:
                    unique_contacts.add(contact_id)
                    # Update the count for the fiber in the middle
                    if fiber_a in fibers_in_middle:
                        idx_a = fiber_index_mapping[fiber_a]
                        fiber_contacts[idx_a] += 1
                    if fiber_b in fibers_in_middle:
                        idx_b = fiber_index_mapping.get(fiber_b)  # Use get() to avoid KeyError for fibers outside the middle
                        if idx_b is not None:  # Only update if fiber_b is also in the middle
                            fiber_contacts[idx_b] += 1

    return fiber_contacts

# Compute the pipeline to apply modifications
data = pipeline.compute()

# Determine the middle of the slab based on z-coordinates
z_coords = data.particles['Position'][:,2]
average_z = np.mean(z_coords)
std_dev_z = np.std(z_coords)
middle_start = average_z - std_dev_z/2
middle_end = average_z + std_dev_z/2
print(f"Middle region z-direction: {middle_start} to {middle_end}")

# Initialize a set to track fibers in the middle
fibers_in_middle = set()

# Map molecule IDs to fibers
molecule_ids = data.particles['Molecule Identifier']
unique_molecule_ids = np.unique(molecule_ids)
molecules_per_fiber = 12  # Assuming 12 molecules per fiber
total_fibers = len(unique_molecule_ids) // molecules_per_fiber

# Identify fibers with at least one molecule in the middle
for molecule_id in unique_molecule_ids:
    if molecule_id == 0:  # Skip molecule ID 0 if not to be considered
        continue
    positions = data.particles.positions[data.particles['Molecule Identifier'] == molecule_id]
    z_positions = positions[:,2]
    if any((z_positions >= middle_start) & (z_positions <= middle_end)):
        fiber_id = (molecule_id - 1) // molecules_per_fiber
        fibers_in_middle.add(fiber_id)

print(f"Number of fibers recognized in the middle: {len(fibers_in_middle)}")

# Calculate valency using the new function
valencies = calculate_valency(data, fibers_in_middle, 12)  # Adjust molecules_per_fiber if necessary

# Plot the distribution of valencies for fibers in the middle
plt.hist(valencies[valencies > 0], bins=np.arange(1, max(valencies) + 1) - 0.5, edgecolor='black')
plt.title('Distribution of Valencies for Fibers in the Middle')
plt.xlabel('Valency')
plt.ylabel('Frequency')
plt.savefig('valency_distribution.png')
print("Valency distribution plot saved as 'valency_distribution.png'.")

