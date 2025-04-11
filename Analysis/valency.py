from ovito.io import import_file
from ovito.modifiers import CreateBondsModifier
import matplotlib
matplotlib.use('Agg')  # Set the backend before importing pyplot
import matplotlib.pyplot as plt
import numpy as np


# Parameters
molecules_per_fiber = 12  # Number of molecules per fiber
cutoff_distance = 40  # Cutoff distance for contacts, adjust as needed
total_fibers = 125  # Adjust based on your total number of fibers

# Load the dataset
pipeline = import_file('172/dna.dump')

# Apply the CreateBondsModifier with an appropriate cutoff
pipeline.modifiers.append(CreateBondsModifier(cutoff=cutoff_distance))

# Compute the pipeline to apply modifications
data = pipeline.compute()

def calculate_valency(data):
    # Access the particle properties, including molecule IDs
    molecule_ids = data.particles['Molecule Identifier'].array
    num_molecules = np.max(molecule_ids)  # Assuming molecule IDs are sequential starting from 1

    # Initialize contact counts for each fiber
    fiber_contacts = np.zeros(total_fibers, dtype=int)

    # Access computed bonds and bond topology
    bond_topology = data.particles.bonds.topology

    # Dictionary to map molecule ID to its fiber, excluding molecule ID 0
    molecule_to_fiber = {mol_id: (mol_id - 1) // molecules_per_fiber for mol_id in range(1, num_molecules + 1) if mol_id > 0}

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

        # Check if molecules belong to different fibers and are mapped in molecule_to_fiber
        if fiber_a is not None and fiber_b is not None and fiber_a != fiber_b:
            unique_contacts.add((min(fiber_a, fiber_b), max(fiber_a, fiber_b)))

    # Count the unique contacts for each fiber
    for contact in unique_contacts:
        fiber_contacts[contact[0]] += 1
        fiber_contacts[contact[1]] += 1

    return fiber_contacts


# Apply the function and get the output
valency = calculate_valency(data)

# Print or save the valency results
for fiber_id, val in enumerate(valency):
    print(f"Fiber {fiber_id + 1} valency: {val}")

# Plotting function
def plot_valency_distribution(valencies):
    plt.figure(figsize=(10, 6))
    plt.bar(range(1, len(valencies) + 1), valencies, color='skyblue')
    plt.xlabel('Fiber ID')
    plt.ylabel('Valency')
    plt.title('Distribution of Valencies Across Fibers')
    plt.xticks(range(1, len(valencies) + 1))
    plt.savefig('valency_distribution_172.png')  # Save the plot as a PNG file
    plt.close()  # Close the plot to free up memory

# Plot the valency distribution
plot_valency_distribution(valency)

