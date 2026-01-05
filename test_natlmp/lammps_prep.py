# Make the lmpdata and even the script¿
import numpy as np

# Geometry
numof_particles = 100
packing = 0.3
height = 4
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


# Add this to the end
'''
This is the initial atom setup of ./data.lmpdata
9 atoms 
1 atom types 
0 bonds 
0 bond types 
-7.5 7.5 xlo xhi 
-7.5 7.5 ylo yhi 
-1.95 1.95 zlo zhi 

Atoms

# ID|type|x|y|z|q|mx|my|mz|d|rho

PairIJ Coeffs

1 1 lj/cut/dipole/cut 0.01 2.49451641079295 2.8 20
'''
# 1 over 1 interaction, lennard jones dipole, and parameters...

