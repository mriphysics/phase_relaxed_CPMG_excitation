# Phase Relaxed CPMG Excitation
Phase relaxed method for CPMG excitation pulses. The method was presented at the ISMRM conference in 2014 (abstract [here](https://kclpure.kcl.ac.uk/portal/files/38148970/0948.pdf)) and is currently under review for journal publication. A pre-print will be posted on acceptance. This repo contains Matlab code for implementing 2D or 3D pulse designs.

Shaihan Malik, July 2015

### [Releases](link)
Code is available as a release including binary files (Matlab .mat files) containing B1 and B0 field maps. Please see the releases page for more.

Note that a persistent version of release 1.0 has been linked: [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.19957.svg)](http://dx.doi.org/10.5281/zenodo.19957)

### Usage
`example1_spiral2D.m` and `example2_shells3d.m` are Matlab scripts that run through 2D spiral and 3D 'shells' k-space trajectory versions of the pulse design. They rely on B1 and B0 maps available on the releases page.

The code is dependent on the [reverse-GIRF](link) repo and on [mcgls.m](http://m2matlabdb.ma.tum.de/download.jsp?MC_ID=3&SC_ID=10&MP_ID=126). Detailed information is given in the next section.

This code is provided freely under the terms of the MIT license. Please cite appropriately.

### Dependencies  

#### Gradient design
The method uses the time-optimal gradient design method [(Lustig et al, 2008)](http://doi.org/10.1109/TMI.2008.922699) and VERSE guided pulse design method [(Lee et al, 2009)](http://doi.org/10.1002/mrm.21950). An implementation of the Lustig 2008 method (`minTimeGradient`) is available  [here](http://www.eecs.berkeley.edu/~mlustig/Software.html).

Our own implementation of the Lee 2009 method is available as a separate git repo (link) from @mriphysics. The necessary parts of Lusitg's original code are included in this repo with permission. Please cite appropriately.

#### multi-shift CGLS
Least squares minimizations are solved using the multi-shift approach from [van den Eshof et al, 2003](http://doi.org/10.1016/j.apnum.2003.11.010). We used `mcgls.m` that can be downloaded from [here](http://m2matlabdb.ma.tum.de/download.jsp?MC_ID=3&SC_ID=10&MP_ID=126).

If this is unavailable, any other method for least squares solutions of non-square matrices can be used. We have found LSQR (matlab implementation) to work well in the past.
