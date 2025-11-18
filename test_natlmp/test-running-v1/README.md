Files that generate a sucessful run - generation of pairs due to dipole moment - I run this in Linux using:

`mpirun -np 4 ~/lammpsmgcl/build/lmp -in highfreq.lmpin`

Build generated as:

```
git clone -b release https://github.com/lammps/lammps.git mylammps

cmake ../cmake -D PKG_DIPOLE=yes -D PKG_COLLOID=yes -D PKG_BROWNIAN=yes -D PKG_MOLECULE=yes -D PKG_RIGID=yes -D PKG_KSPACE=yes -D BUILD_MPI=yes -D PKG_MC=yes -D PKG_OPENMP=yes

cmake --build . -j$(nproc)
```
