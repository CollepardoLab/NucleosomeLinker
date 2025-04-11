import os
from ovito.io import import_file
from ovito.modifiers import CreateBondsModifier
import numpy as np

# Define data directory
data_directory = '172'

# Import data file
pipeline = import_file(os.path.join(data_directory, 'v3/dna.dump'))

# Create bond topology
pipeline.modifiers.append(CreateBondsModifier(cutoff=120))  # Adjust cutoff value as needed

# Compute bonds
data = pipeline.compute(0)

# Get bond information
bond_topology = data.particles.bonds.topology

# Get molecule IDs, types, and quaternion components
molecule_ids = data.particles['Molecule Identifier'].array
types = data.particles['Particle Type'].array

# Extract quaternion components individually
q1 = data.particles['c_q[1]'].array
q2 = data.particles['c_q[2]'].array
q3 = data.particles['c_q[3]'].array
q4 = data.particles['c_q[4]'].array

# Combine quaternion components into a single numpy array
quaternions = np.vstack((q1, q2, q3, q4)).T
positions = data.particles.positions.array

#quaternions = data.particles[['c_q[1]', 'c_q[2]', 'c_q[3]', 'c_q[4]']].array
#positions = data.particles[['Position.X', 'Position.Y', 'Position.Z']].array

# Initialize list to store bond counts per fiber
fiber_bond_counts = {}
fiber_bond_counts_face_to_face = {}
fiber_bond_counts_side_to_side = {}
fiber_bond_counts_face_to_side = {}

def angle_between_vectors(v1, v2):
    unit_vector_1 = v1 / np.linalg.norm(v1)
    unit_vector_2 = v2 / np.linalg.norm(v2)
    dot_product = np.dot(unit_vector_1, unit_vector_2)
    angle = np.arccos(np.clip(dot_product, -1.0, 1.0))
    print(angle*57.2958)
    return angle

def quaternion_to_z_vector(q):
    w, x, y, z = q
    Zx = 2 * (x * z + w * y)
    Zy = 2 * (y * z - w * x)
    Zz = 1 - 2 * (x**2 + y**2)
    return np.array([Zx, Zy, Zz])

def quaternion_to_rotation_matrix(q):
    """Convert a quaternion into a rotation matrix."""
    w, x, y, z = q
    return np.array([
        [1 - 2*y*y - 2*z*z,     2*x*y - 2*z*w,     2*x*z + 2*y*w],
        [2*x*y + 2*z*w,     1 - 2*x*x - 2*z*z,     2*y*z - 2*x*w],
        [2*x*z - 2*y*w,     2*y*z + 2*x*w,     1 - 2*x*x - 2*y*y]
    ])

# Iterate over bonds
for bond in bond_topology:
    particle1, particle2 = bond  # Indices of bonded particles
    
    if types[particle1]==1 and types[particle2]==1:
        # Get molecule IDs of bonded particles
        molecule_id1, molecule_id2 = molecule_ids[particle1], molecule_ids[particle2]
        
        # Calculate bond vector
        bond_vector = positions[particle2] - positions[particle1]
        bond_length = np.linalg.norm(bond_vector)
        bond_vector_normalized = bond_vector / bond_length
            
        # Determine fiber IDs of molecules
        fiber_id1 = (molecule_id1 - 1) // 12 + 1  # Assuming 12 molecules per fiber, indexed from 1
        fiber_id2 = (molecule_id2 - 1) // 12 + 1

        # Check if both molecules belong to the same fiber and are separated by 2 in sequence
        if fiber_id1 != fiber_id2:
            # Ensure molecule_id1 < molecule_id2 to avoid duplicates
            if molecule_id1 < molecule_id2:
                bond_pair = (molecule_id1, molecule_id2)
            else:
                bond_pair = (molecule_id2, molecule_id1)
            
            # Get quaternion and convert to rotation matrices for each particle
            q1 = quaternions[particle1]
            q2 = quaternions[particle2]

            z1 = quaternion_to_z_vector(q1)
            z2 = quaternion_to_z_vector(q2)

            # print(bond_vector_normalized, z1, z2)
            # Determine minimum angle between bond vector and rotated principal axes for both particles
            angles_p1 = [angle_between_vectors(bond_vector_normalized, z1)]
            angles_p2 = [angle_between_vectors(bond_vector_normalized, z2)]
    
            min_angle_p1 = np.min(angles_p1)
            min_angle_p2 = np.min(angles_p2)
            if ( min_angle_p1 < np.pi/4 or min_angle_p1 > 3*np.pi/4 ) and ( min_angle_p2 < np.pi/4 or min_angle_p2 > 3*np.pi/4 ):
                bond_type = 'face-to-face'
                # Correctly initialize or update the specific category dictionary
                if fiber_id1 not in fiber_bond_counts_face_to_face:
                    fiber_bond_counts_face_to_face[fiber_id1] = set()
                fiber_bond_counts_face_to_face[fiber_id1].add(bond_pair)
            elif (min_angle_p1 >= np.pi/4 and min_angle_p1<= 3*np.pi/4) and (min_angle_p2 >= np.pi/4 and min_angle_p2<= 3*np.pi/4):
                bond_type = 'side-to-side'
                # Same for side-to-side
                if fiber_id1 not in fiber_bond_counts_side_to_side:
                    fiber_bond_counts_side_to_side[fiber_id1] = set()
                fiber_bond_counts_side_to_side[fiber_id1].add(bond_pair)
            else:
                bond_type = 'face-to-side'
                # And for face-to-side
                if fiber_id1 not in fiber_bond_counts_face_to_side:
                    fiber_bond_counts_face_to_side[fiber_id1] = set()
                fiber_bond_counts_face_to_side[fiber_id1].add(bond_pair)
            
            # Still increment the total bond count for fiber_id1
            if fiber_id1 not in fiber_bond_counts:
                fiber_bond_counts[fiber_id1] = set()
            fiber_bond_counts[fiber_id1].add(bond_pair)


bond_counts = []  
# Initialize lists to store bond counts for each category
total_bond_counts = []  # Make sure this is defined
ff_counts = []
ss_counts = []
fs_counts = []

# Combined loop to analyze bond counts across all fibers and types
for fiber_id in sorted(set(fiber_bond_counts.keys()) | 
                       set(fiber_bond_counts_face_to_face.keys()) | 
                       set(fiber_bond_counts_side_to_side.keys()) | 
                       set(fiber_bond_counts_face_to_side.keys())):
    total_bonds = fiber_bond_counts.get(fiber_id, set())
    ff_bonds = fiber_bond_counts_face_to_face.get(fiber_id, set())
    ss_bonds = fiber_bond_counts_side_to_side.get(fiber_id, set())
    fs_bonds = fiber_bond_counts_face_to_side.get(fiber_id, set())

    # Append counts to respective lists
    total_bond_counts.append(len(total_bonds))
    ff_counts.append(len(ff_bonds))
    ss_counts.append(len(ss_bonds))
    fs_counts.append(len(fs_bonds))

    print(f"Fiber {fiber_id}: Total Bonds: {len(total_bonds)}, Face-to-Face: {len(ff_bonds)}, Side-to-Side: {len(ss_bonds)}, Face-to-Side: {len(fs_bonds)}")

# Function to calculate and print mean and std deviation, now correctly returning values
def calculate_stats(counts, label):
    mean_counts = np.mean(counts)
    std_counts = np.std(counts)
    print(f"\nAverage {label}: {mean_counts}")
    print(f"Standard deviation of {label}: {std_counts}")
    return mean_counts, std_counts  # Ensure this line is correctly returning values

# Calculate and print stats for each bond type and total bonds
mean_total, std_total = calculate_stats(total_bond_counts, 'total bond counts')
mean_ff, std_ff = calculate_stats(ff_counts, 'face-to-face bond counts')
mean_ss, std_ss = calculate_stats(ss_counts, 'side-to-side bond counts')
mean_fs, std_fs = calculate_stats(fs_counts, 'face-to-side bond counts')

# Define the file to append the results to
output_file_path = 'analysis_results_inter_2.csv'  # Use a specific path for writing in this environment

# Open the file in append mode and write the results
with open(output_file_path, 'a') as file:
    # Check if file is empty and write header accordingly
    file.seek(0, 2)  # Move to the end of the file
    if file.tell() == 0:  # File is empty, write header
        file.write('Data Directory,Mean Total Bonds,Std Total Bonds,Mean FF,Std FF,Mean FS,Std FS,Mean SS,Std SS\n')
    
    # Write the data, ensure variable names match
    file.write(f"{data_directory},{mean_total},{std_total},{mean_ff},{std_ff},{mean_fs},{std_fs},{mean_ss},{std_ss}\n")

print(f"Results for directory {data_directory} appended to {output_file_path}.")

