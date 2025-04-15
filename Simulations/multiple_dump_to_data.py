import numpy as np
import sys

def dump_frame_to_data_file(dump_fname, og_data_fname):
    dump_file = [line.rstrip('\n') for line in open(dump_fname)]
    og_data_file = [line.rstrip('\n') for line in open(og_data_fname)]

    # Get the number of atoms and ellipsoids from the original data file
    num_atoms = int(og_data_file[2].split()[0])
    num_ellipsoids = int(og_data_file[3].split()[0])

    # Extract box dimensions from the dump file
    xx = dump_file[5]
    yy = dump_file[6]
    zz = dump_file[7]

    # Collect atom positions and ellipsoid quaternions
    values = []
    for i in range(9, num_atoms + 9):
        values.append(np.array(dump_file[i].split(), dtype=float))

    # Update box size in the original data file
    og_data_file[13] = xx + " xlo xhi"
    og_data_file[14] = yy + " ylo yhi"
    og_data_file[15] = zz + " zlo zhi"

    # Update atom coordinates
    for i in range(19, num_atoms + 19):
        line = og_data_file[i].split()
        j = int(line[0]) - 1
        new_line = line[:2] + [str(values[j][1]), str(values[j][2]), str(values[j][3])] + line[5:]
        og_data_file[i] = " ".join(new_line)

    # Update ellipsoid quaternions
    for i in range(19 + num_atoms + 3, 19 + num_atoms + 3 + num_ellipsoids):
        line = og_data_file[i].split()
        j = int(line[0]) - 1
        new_line = line[:4] + [str(values[j][4]), str(values[j][5]), str(values[j][6]), str(values[j][7])]
        og_data_file[i] = " ".join(new_line)

    return og_data_file

# Main loop to process each dump file
og_data_fname = sys.argv[1]  # Original data file as the second argument

for frame_num in range(1, 121):
    dump_fname = f"frame_{frame_num}.dump"  # Each frame dump file
    new_data_fname = f"new_data_{frame_num}.txt"  # Output data file

    # Generate data file content from the dump frame
    dat_file_lines = dump_frame_to_data_file(dump_fname, og_data_fname)

    # Write the new data file
    with open(new_data_fname, "w") as dfile:
        for line in dat_file_lines:
            dfile.write(line + "\n")

    print(f"Created {new_data_fname}")


