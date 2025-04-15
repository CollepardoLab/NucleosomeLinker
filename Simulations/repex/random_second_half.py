import random

# Path to the LAMMPS dump file
dump_file_path = 'equil_T0.dump'

# Read and parse the dump file to isolate frames
with open(dump_file_path, 'r') as file:
    lines = file.readlines()

frames = []
frame = []

# Separate each frame
for line in lines:
    if line.startswith("ITEM: TIMESTEP"):
        if frame:  # Add the current frame to the frames list
            frames.append(frame)
        frame = [line]  # Start a new frame
    else:
        frame.append(line)

# Append the last frame if present
if frame:
    frames.append(frame)

# Select frames from the second half only
halfway_index = len(frames) // 2
second_half_frames = frames[halfway_index:]

# Select 120 random frames from the second half
selected_second_half_frames = random.sample(second_half_frames, 120)

# Write selected frames to a new file
output_file_second_half_path = 'random_120_frames_second_half.dump'
with open(output_file_second_half_path, 'w') as file:
    for frame in selected_second_half_frames:
        file.writelines(frame)

print(f"Selected 120 random frames from the second half written to {output_file_second_half_path}")

