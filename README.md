GCM Necro
=========

Necrocode from legacy gcm3d branches

Compile
=======

mkdir build/
cd build/
cmake ../gcm/
make

Prepare
=======
cd tasks/
gmsh -3 cover-small.geo
gmsh -3 cube.geo

Run
===
cd build/
./gcm3d --task ../tasks/friction-test.xml --zones ../tasks/basic-contact-task-1cpu-zones.xml --data-dir ../tasks/ --dump-vtk @

Customize
=========
Mesh size is controlled by meshPointDist parameter in geo files. Do not forget to re-create msh file after the parameter was changed!
