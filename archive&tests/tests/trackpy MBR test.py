import trackpy as tp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Example: Load or simulate data
# For this example, we'll simulate some simple particle trajectories

# Simulate a DataFrame with particle positions over time
np.random.seed(42)
n_frames = 100
n_particles = 5
data = []

for t in range(n_frames):
    for particle_id in range(n_particles):
        x = np.sin(0.1 * t) + 0.1 * np.random.randn()
        y = np.cos(0.1 * t) + 0.1 * np.random.randn()
        data.append([t, particle_id, x, y])

# Convert to DataFrame
df = pd.DataFrame(data, columns=['frame', 'particle', 'x', 'y'])

# Track the particles
# If you have real data, you would first locate features with tp.locate()
# and then link them with tp.link(), but for this example we already have tracks.

# Calculate displacements
# This assumes that 'particle' and 'frame' columns are present in df.
df['x_disp'] = df.groupby('particle')['x'].diff().fillna(0)
df['y_disp'] = df.groupby('particle')['y'].diff().fillna(0)

# Calculate mean back relaxation
# For simplicity, let's average the displacements over all particles
# and plot them as a function of time.

mean_displacements = df.groupby('frame')[['x_disp', 'y_disp']].mean()

# Plot the results
plt.plot(mean_displacements.index, mean_displacements['x_disp'], label='Mean x displacement')
plt.plot(mean_displacements.index, mean_displacements['y_disp'], label='Mean y displacement')
plt.plot(np.zeros())
plt.xlabel('Frame')
plt.ylabel('Mean displacement')
plt.legend()
plt.title('Mean Back Relaxation')
plt.show()
