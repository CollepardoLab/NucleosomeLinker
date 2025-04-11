import numpy as np
import matplotlib.pyplot as plt
import os

# Parameters for the two cases
cases = [
    {
        "file_path": "172/dna.dump",
        "atoms_per_molecule": 424,
        "highlighted_rgs": [177.5962962576766, 198.396908578065],
        "label": "172 experimental data",
        "color": "coral"
    },
    {
        "file_path": "177/dna.dump",
        "atoms_per_molecule": 436,
        "highlighted_rgs": [139.70285927241753, 129.2130024417047, 151.75565155136002, 120.61123109672104, 148.27152270225207, 132.25228416426953, 132.25228416426953],
        "label": "177 data",
        "color": "green"
    }
]

# Function to calculate radius of gyration
def radius_of_gyration(atom_positions):
    positions = np.array(atom_positions)
    com = positions.mean(axis=0)
    rg = np.sqrt(np.mean(np.sum((positions - com) ** 2, axis=1)))
    return rg

# Initialize lists to store Rg values for each case
rg_values_cases = []

# Process each case
for case in cases:
    rg_values = []
    file_path = case["file_path"]
    atoms_per_molecule = case["atoms_per_molecule"]

    # Read the dump file and calculate Rg for each molecule
    with open(file_path, 'r') as file:
        current_atoms = []
        in_atoms_section = False
        for line in file:
            if line.startswith("ITEM: ATOMS"):
                in_atoms_section = True
                current_atoms = []
                continue
            if line.startswith("ITEM:") and in_atoms_section:
                # Calculate Rg for the current molecule
                if current_atoms:
                    molecules = {}
                    for atom_data in current_atoms:
                        atom_id, x, y, z = atom_data
                        mol_id = (atom_id - 1) // atoms_per_molecule + 1
                        if mol_id not in molecules:
                            molecules[mol_id] = []
                        molecules[mol_id].append([x, y, z])

                    # Calculate Rg for each molecule in this frame
                    for atom_positions in molecules.values():
                        rg = radius_of_gyration(atom_positions)
                        rg_values.append(rg)
                in_atoms_section = False
            elif in_atoms_section:
                data = line.strip().split()
                if len(data) >= 4:
                    atom_id = int(data[0])
                    x = float(data[1])
                    y = float(data[2])
                    z = float(data[3])
                    current_atoms.append((atom_id, x, y, z))
    
    # Store the Rg values for this case
    rg_values_cases.append((rg_values, case["highlighted_rgs"], case["label"], case["color"]))

# Plot the histograms as filled areas with shading
plt.figure(figsize=(10, 6))
for rg_values, highlighted_rgs, label, color in rg_values_cases:
    # Normalize histogram and plot as a filled area
    counts, bin_edges = np.histogram(rg_values, bins=50, density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    plt.fill_between(bin_centers, counts, color=color, alpha=0.3, label=label)
    plt.plot(bin_centers, counts, color=color, linewidth=2)

    # Add vertical lines for highlighted Rg values
    for rg_value in highlighted_rgs:
        plt.axvline(rg_value, color=color, linestyle='--', linewidth=1.5, alpha=0.8)

# Customize the plot
plt.xlabel('Radius of Gyration (Rg)')
plt.ylabel('Normalized Frequency')
plt.title('Normalized Rg Distribution for Experimental and Computational Data')
plt.legend()
plt.savefig("Rgs/normalized_rg_distribution.svg")
plt.show()

