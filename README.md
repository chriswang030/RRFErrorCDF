# RRFErrorCDF

This repository holds the code used in the paper "Beyond singular value gaps in randomized subspace approximation," coauthored with Alex Townsend and available on the [arXiv](https://arxiv.org/abs/2603.01191). Comments and suggestions welcome.

## Dependencies

Our code requires the `mhg` package of Plamen Koev and Alan Edelman, which is available for download [here](https://sites.google.com/sjsu.edu/plamenkoev/home/software/mhg). This code is the backbone of our algorithm for computing the CDF of the randomized range finder's approximation error; for theory and implementation details, see
> Plamen Koev, Alan Edelman. The efficient evaluation of the hypergeometric function of a matrix argument, _Math. Comp._ 75 (2006), pp. 833-846.

## Usage 

Include the `mhg` package in your MATLAB path. Run `main.m` to generate the plots in Figure 1 of the paper. The main algorithm for computing the CDF can be found in `cdf.m` with arguments described therein. For the plots in Figures 2-4, run `main_bounds.m`, setting the oversampling parameter `p` to be either 1 or 5 (or something of your own choosing).
