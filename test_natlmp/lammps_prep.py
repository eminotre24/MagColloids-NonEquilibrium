# Make the lmpdata and even the script¿
import numpy as np
import os

# Geometry
numof_particles = 1024
packing = 0.3
height = 3.9
radius = 1.4

# Parameters
# ID|type|x|y|z|q|mx|my|mz|d|rho - Use atoms style hybrid: dipole+sphere according to order of atom_style definition
# ID should be counted (start 1), and the type should be 1 if theres no trap. Positions. 0 charge. Momenta in z: 0.0001. Diameter 2.8. Density 0.00514794182.

# Routine from Magcolloids

# Redo Particles Squared
part_at_side = int(np.round(np.sqrt(numof_particles)))
numof_particles = part_at_side**2

area_particles = numof_particles*radius**2*np.pi
area_region = area_particles/packing

length_region = np.sqrt(area_region)
part_separation = length_region/part_at_side # Puntual Position

# Region
region = [np.round(length_region),np.round(length_region),height]

# Positions
x_loc = np.linspace(
        -length_region/2+part_separation/2,
        length_region/2-part_separation/2,part_at_side)
y_loc = np.linspace(
    -length_region/2+part_separation/2,
    length_region/2-part_separation/2,part_at_side)

[X,Y] = np.meshgrid(x_loc,y_loc)
Z = np.zeros(np.shape(X))

initial_positions = np.array([[x,y,z] for (x,y,z) in zip(X.flatten(),Y.flatten(),Z.flatten())])

# Make Document
name = "parts.lmpdata"
xy_lim = round(length_region/2, 2)
z_lim = round(height/2, 2)

with open(name, "w") as f:
    f.write(f"This is the initial atom setup of ./{name}\n")
    f.write(f"{int(numof_particles)} atoms\n")
    f.write("1 atom types\n")
    f.write("0 bonds\n")
    f.write("0 bonds types\n")
    f.write(f"-{xy_lim} {xy_lim} xlo xhi\n")
    f.write(f"-{xy_lim} {xy_lim} ylo yhi\n")
    f.write(f"-{z_lim} {z_lim} zlo zhi\n\n")

    f.write("Atoms\n\n")

    for part in range(numof_particles):
        f.write(f"{part+1} 1 ")
        f.write(f"{initial_positions[part][0]} {initial_positions[part][1]} {initial_positions[part][2]} ")
        f.write("0 0 0 0.0001 ")
        f.write(f"{radius*2} 0.005147941820269123 ")
        f.write("#ID|type|x|y|z|q|mx|my|mz|d|rho\n")

    f.write("\nPairIJ Coeffs\n\n")
    f.write(f"1 1 lj/cut/dipole/cut 0.01 2.49451641079295 {radius*2} 20")
    # 1 over 1 interaction, lennard jones dipole, and parameters...

