econd
==============

A framework that calculates electrical conductivities from MD trajectories and performs various decomposition methods

### Note

1. With the Fortran program, 'decompose_mpi --nosd' will output a binary file, say, corr-nosd.h5, which contains only the array of C_I(t) and C_IL(t). I call that array as nCorr (the SD array is called sdCorr). Some other information about the system is also written in corr-nosd.h5, such as the charges and the box dimensions. corr-nosd.h5 has to be processes by some Python scripts to actually give us the EC results.

2. Use the Python script 'calCesary.py' to integrate nCorr, and we have a file called 'cesaro.h5' with an array inside called nDCesaro. Note that step 1 and 2 are repeatedly applied to the trajectory of each ns. So for a simulation of 10-ns, we will have 10 corr-nosd.h5 and 10 cesaro.h5.

3. Provide the 10 cesaro.h5 to the Python script 'fitCesaro.py', which will average them, fit them, and does error analysis, according to the fitting ranges input by us. Let's call the output file cesaro.fit.h5.

4. Finally, send cesaro.fit.h5 to the Python script 'plotFit.py', we will then have numerical results (EC, diffusion constants ...) printed out, and the Cesaro plot, as well as g-D-sig figure if in step 1 we did not use --nosd option.
