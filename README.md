# Phase Relaxed CPMG Excitation
Phase relaxed method for CPMG excitation pulses. The method was presented at the ISMRM conference in 2014 (abstract [here](https://kclpure.kcl.ac.uk/portal/files/38148970/0948.pdf)) and as a [full paper](http://onlinelibrary.wiley.com/doi/10.1002/mrm.25996/abstract) in 2016. This repo contains Matlab code for implementing 2D or 3D pulse designs.

Shaihan Malik, July 2015

### [Releases](https://github.com/mriphysics/phase_relaxed_CPMG_excitation/releases)
Code is available as a release including binary files (Matlab .mat files) containing B1 and B0 field maps. Please see the releases page for more.

A persistent citeable version of release 1.0 is available: [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.19957.svg)](http://dx.doi.org/10.5281/zenodo.19957)

### Usage
`example1_spiral2D.m` and `example2_shells3d.m` are Matlab scripts that run through 2D spiral and 3D 'shells' k-space trajectory versions of the pulse design. They rely on B1 and B0 maps available on the releases page.

The code is dependent on the [reverse-GIRF](https://github.com/mriphysics/reverse-GIRF) project and on [mcgls.m](http://m2matlabdb.ma.tum.de/download.jsp?MC_ID=3&SC_ID=10&MP_ID=126). See next section.

This code is provided freely under the terms of the MIT license. Please cite appropriately.

### Dependencies  

#### Gradient design
The pulse design method uses the time-optimal gradient design method [(Lustig et al, 2008)](http://doi.org/10.1109/TMI.2008.922699) and VERSE based pulse design method [(Lee et al, 2009)](http://doi.org/10.1002/mrm.21950). An implementation of the Lustig 2008 method (`minTimeGradient`) is available  [here](http://www.eecs.berkeley.edu/~mlustig/Software.html).

We have extended [Lustig's code](http://www.eecs.berkeley.edu/~mlustig/Software.html) to implement [(Lee et al, 2009)](http://doi.org/10.1002/mrm.21950) and also a gradient impulse response function (GIRF) based correction for non-ideal gradient systems. **The [reverse-GIRF](https://github.com/mriphysics/reverse-GIRF) project contains all necessary code**  including Lustig's original source (with permission). Please take note of the citations suggested in the [readme](https://github.com/mriphysics/reverse-GIRF/readme).

#### multi-shift CGLS
Least squares minimizations are solved using the multi-shift approach from [van den Eshof et al, 2003](http://doi.org/10.1016/j.apnum.2003.11.010). We used `mcgls.m` that can be downloaded from [here](http://m2matlabdb.ma.tum.de/download.jsp?MC_ID=3&SC_ID=10&MP_ID=126).

If this is unavailable, any other method for least squares solutions of non-square matrices can be used. We have found LSQR (matlab implementation) to work well in the past.
