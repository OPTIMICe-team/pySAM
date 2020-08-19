#! /bin/bash

# An old compilation tool I was used to use before moving to Makefiles
# Might be worth keeping for 

rm dat/*.dat
rm vtk/*.vtk
rm *.dat
rm *.vtk
g++ -O2 -I ./include/ -fopenmp aggregation.cpp src/aggregate.cpp src/column.cpp src/dendrite.cpp src/geom_lib.cpp src/math_lib.cpp src/hex_prism.cpp src/plate.cpp src/population.cpp src/pristine.cpp src/rosette.cpp -o aggregation
g++ -O2 -I ./include/ -fopenmp collection.cpp src/aggregate.cpp src/column.cpp src/dendrite.cpp src/geom_lib.cpp src/math_lib.cpp src/hex_prism.cpp src/plate.cpp src/population.cpp src/pristine.cpp src/rosette.cpp -o collection