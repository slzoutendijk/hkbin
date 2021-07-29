# hkbin

hkbin is an alternative binning code for [pyGravSphere][].  It is designed to
give unbiased estimates of the intrinsic velocity dispersion and its
uncertainty, even in the case of large measurement uncertainties and a small
number of velocity measurements.  hkbin was designed for the dynamical analysis
of ultra-faint dwarf galaxies with ~100 velocity measurements in total.

## Changes with respect to pyGravSphere

hkbin is based on the binning code of [pyGravSphere][], but with a few
algorithmic changes.  The main difference is in the estimators of the intrinsic
velocity dispersion and its uncertainty, which are derived from the sample
moments using _h_ and _k_ statistics.  The _h_ and _k_ statistics are
symmetrically unbiased estimators of the population moments and cumulants.  The
bias in the resulting intrinsic velocity dispersion is at most a few per cent
for an ultra-faint dwarf galaxy, which is much smaller than the statistical
uncertainty.

Estimators of the virial shape parameters (VSPs) and their uncertainties are
currently not implemented.  Estimating the uncertainties on the VSPs would
require the calculation of the eighth velocity moment in each bin, which is not
feasible for the number of velocity measurements currently available for
ultra-faint dwarf galaxies.  However, from the statistics side, implementing
estimators for these higher moments should be possible.  If you have a use for
these and are willing to dig into a mathematical statistics paper from 1937,
please file an issue or [contact me](#contact)!

The drawback of the new algorithm is that it is possible to obtain negative,
unphysical, estimates of the intrinsic velocity dispersion.  This is due to the
unbiased nature of the estimators.  It promises that the average results over
many bins or many random realizations of measurements will be correct.  Because
of statistical uncertainties, individual estimates will deviate, and it may
happen that the correction for the measurement uncertainty is larger than the
measured dispersion.  Bins with negative intrinsic dispersion estimates are
discarded from further analysis.

Further changes are: a configurable number of stars per bin, adding remaining
stars to the last bin, using the mean instead of the maximum radius of stars in
a bin as the bin radius, and allowing the fourth moment to be negative.

## Using hkbin

Running hkbin requires a working installation of pyGravSphere, nothing else.
Installation is done by simply copying the source to a directory of your
choosing.  Like pyGravSphere, hkbin currently only works with Python 2.

hkbin can be used in place of the pre-processing step of [pyGravSphere][].
If you have not done so yet, please read the pyGravSphere documentation.
You will have to prepare the `KinPhotDat` directory in the same way.
Optionally, you can configure the number of stars per photometric or kinematic
bin by writing the appropriate number to a configuration file (in this example
for a galaxy named `Galaxy_1`):

```
echo 100 >KinPhotDat/Galaxy_1_PhotBin.txt
echo 11 >KinPhotDat/Galaxy_1_KinBin.txt
```

If you do not configure the number of stars per bin, hkbin will use the same
default as pyGravSphere, which is `floor(N/sqrt(N))`.  This will likely be too
small to give reliable results for an ultra-faint dwarf galaxy.  I recommend
you test a few different values and see what the minimum is to get a stable
density profile.

To perform the binning, in your pyGravSphere working directory, call

```
GravSpherePath=~/software/pyGravSphere/ python2 ~/software/hkbin/hkbin.py Galaxy_1
```

Replace `~/software/pyGravSphere/` and `~/software/hkbin/` with the appropriate
installation directories and `Galaxy_1` with your galaxy name.  You can run
hkbin for multiple galaxies serially by listing multiple galaxy names,
separated by spaces.

hkbin should have generated files in the `GalaxyData` directory.  You can
continue your analysis with pyGravSphere like you would after finishing the
pre-processing step.

## Citing hkbin

hkbin is based on [pyGravSphere][] and re-uses much of its code.  If you
publish a scientific work using hkbin, please cite the paper introducing hkbin
as well as the (py)GravSphere papers:

[Zoutendijk et al. (2021, A&A, 651, A80)](https://ui.adsabs.harvard.edu/abs/2021A%26A...651A..80Z/abstract)

[Read & Steger (2017, MNRAS, 471, 4541)](https://ui.adsabs.harvard.edu/abs/2017MNRAS.471.4541R/abstract)

[Genina et al. (2020, MNRAS, 498, 144)](https://ui.adsabs.harvard.edu/abs/2020MNRAS.498..144G/abstract)

## Contact

In case of bugs or questions, please file an issue on GitHub or send me an
e-mail: `zoutendijk@strw.leidenuniv.nl`.  If your issue has to do with part of
the code inherited from [pyGravSphere][], please report your issue there.  I
will merge improvements to pyGravSphere into hkbin.

## License

This project is licensed under the terms of the [GNU GPL v3+ license](LICENSE).

[pyGravSphere]: https://github.com/AnnaGenina/pyGravSphere
