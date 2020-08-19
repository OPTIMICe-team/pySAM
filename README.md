# pySAM

This is mostly a backup repository of the aggregation code I used in Ori et al. (2014).
The code has been enriched with a Makefile for faster compilation and easier usage of the different functions.

To use the code You have to properly setup an application main function. The application function used in Ori et al. (2014) has been implemented in collection.cpp; to compile the code simply type
```
make collection
```
and execute the resulting collection program.

An evolution of the code is aggregation.cpp, this will generate aggregates starting from a distribution (or multiple distributions) of monomers following a differential sedimentation collection kernel (similar to previous works by Westbrook or Leinonen). Make the necessary adjustments in aggregation.cpp and compile with
```
make aggregation
```

## Python tools
Hopefully one day I will make a python interface to all the components so that recompilation won't be necessary anymore

## Caveats
The stocastic engine has a random number generator intiated with some function of the current time so every execution will produce slightly different results.
The output files will likely be overwritten or mixed up, so move the output of previous execution separately if you want to keep them

## Limitations
The melting part is still missing, I have to edit it
The code does not include riming... if I ever implement that we could rename the project Melting and Aggregation of (potentially) Rimed Snowflakes (MARS). Let's se if I manage to get to MARS before NASA
