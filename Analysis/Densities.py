import numpy as np

def parse_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    timesteps_data = []
    current_data = []
    for line in lines[1:]:  # Skip the first line assuming it's a header
        if len(line.strip()) == 0:
            continue  # Skip any empty lines

        if line.startswith(" "):  # Data line
            parts = line.split()
            if len(parts) == 4:  # Check if it matches the expected format
                current_data.append([float(parts[1]), float(parts[3])])  # x axis, density
        else:  # Header line
            if current_data:  # Save previous timestep if exists
                timesteps_data.append(np.array(current_data))
                current_data = []
    
    if current_data:  # Add the last timestep
        timesteps_data.append(np.array(current_data))

    return timesteps_data

#def find_plateau_midpoint(density, threshold):
#    high_density_indices = np.where(density > threshold)[0]
#    if len(high_density_indices) == 0:
#        return len(density) // 2
#    
#    start = high_density_indices[0]
#    end = high_density_indices[-1]
#    midpoint = (start + end) // 2
#    return midpoint

def center_data(timesteps_data):
    centered_data = []
    for data in timesteps_data:
        densities = data[:, 1]
        threshold = np.mean(densities) + 2 * np.std(densities)
        #midpoint = find_plateau_midpoint(densities, threshold)
        shift_amount = 0 #len(densities) // 2 - midpoint
        centered_densities = np.roll(densities, shift_amount)
        centered_data.append(centered_densities)
    return centered_data

def average_density_profiles(centered_data):
    avg_profile = np.mean(np.array(centered_data), axis=0)
    return avg_profile

def save_density_profile(avg_profile, filename='average_density_profile.txt'):
    np.savetxt(filename, avg_profile, header='Average Density Profile', comments='')

def calculate_high_density_average(avg_profile):
    threshold = np.mean(avg_profile) + 2 * np.std(avg_profile)
    high_density_values = avg_profile[avg_profile > threshold]
    return np.mean(high_density_values) if high_density_values.size > 0 else 0

if __name__ == "__main__":
    filename = "density.profile"
    timesteps_data = parse_file(filename)
    centered_data = center_data(timesteps_data)
    avg_profile = average_density_profiles(centered_data)
    
    x_axis = np.linspace(0, 8000, len(avg_profile))
    
    high_density_avg = calculate_high_density_average(avg_profile)
    save_density_profile(avg_profile, 'average_density_profile.txt')
    print(f"Final Average of High-Density Phase: {high_density_avg:.4f}")
