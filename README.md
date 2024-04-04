# Calculation of the local density of a nuclear system


The purpose of this project is to calculate the local density of a nuclear system and to plot it.
## How to use

Just compile the entire project with make,
Then execute the ./bin/exec with desired option :
    -a : To Compute acceleration values.
    -r : To Compute Density Values and add them to df3 file.
Finally plot the visu.pov file using POV-ray

> make
> ./bin/exec -a
> ./bin/exec -r
> povray visu.pov


