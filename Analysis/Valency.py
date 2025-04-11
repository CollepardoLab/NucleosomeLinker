from ovito.io import import_file
from ovito.modifiers import CreateBondsModifier
import matplotlib
matplotlib.use('Agg')  # Set the backend before importing pyplot
import matplotlib.pyplot as plt
import numpy as np

molecules_per_fiber = 12  
cutoff_distance = 40  
total_fibers = 125 

pipeline = import_file('172/dna.dump')

pipeline.modifiers.append(CreateBondsModifier(cutoff=cutoff_distance))

data = pipeline.compute()

def calculate_valency(data):
    molecule_ids = data.particles['Molecule Identifier'].array
    num_molecules = np.max(molecule_ids)  

    fiber_contacts = np.zeros(total_fibers, dtype=int)

    bond_topology = data.particles.bonds.topology

    molecule_to_fiber = {mol_id: (mol_id - 1) // molecules_per_fiber for mol_id in range(1, num_molecules + 1) if mol_id > 0}

    unique_contacts = set()

    for bond in bond_topology:
        particle_a, particle_b = bond
        molecule_a, molecule_b = molecule_ids[particle_a], molecule_ids[particle_b]

        if molecule_a == 0 or molecule_b == 0:
            continue

        fiber_a, fiber_b = molecule_to_fiber.get(molecule_a), molecule_to_fiber.get(molecule_b)

        if fiber_a is not None and fiber_b is not None and fiber_a != fiber_b:
            unique_contacts.add((min(fiber_a, fiber_b), max(fiber_a, fiber_b)))

    for contact in unique_contacts:
        fiber_contacts[contact[0]] += 1
        fiber_contacts[contact[1]] += 1

    return fiber_contacts

valency = calculate_valency(data)

for fiber_id, val in enumerate(valency):
    print(f"Fiber {fiber_id + 1} valency: {val}")

def plot_valency_distribution(valencies):
    plt.figure(figsize=(10, 6))
    plt.bar(range(1, len(valencies) + 1), valencies, color='skyblue')
    plt.xlabel('Fiber ID')
    plt.ylabel('Valency')
    plt.title('Distribution of Valencies Across Fibers')
    plt.xticks(range(1, len(valencies) + 1))
    plt.savefig('valency_distribution_172.png')  
    plt.close() 

plot_valency_distribution(valency)

