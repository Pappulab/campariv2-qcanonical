# Overview

This repository implements the q-canonical Monte Carlo sampling method as described in the paper: **Fossat, M. J., & Pappu, R. V.** (2019). [Q-canonical Monte Carlo sampling for modeling the linkage between charge regulation and conformational equilibria of peptides.](https://doi.org/10.1021/acs.jpcb.9b05206) *The Journal of Physical Chemistry B, 123(32)*, 6952-6967.

The method is implemented as an extension to [CAMPARI V2](http://campari.sourceforge.net/), an all-atom Monte Carlo simulation engine for modeling biomolecules, developed by **Andreas Vitalis** in conjunction with the [Pappu lab at Washington University in St. Louis](http://pappu.wustl.edu). The q-canonical extension to CAMPARI V2 is a significant featureset developed by **Martin J. Fossat**, and is thus an unofficial version. It should be noted that CAMPARI is currently at [version 4](http://campari.sourceforge.net/V4) (released in 2020) - consequently, this version of campari should be thought of as a beta-version to a subsequent and updated version of CAMPARI that includes the q-canonical sampling method and the current featureset available in version 4.

# Building & Installation

To build and install this version of CAMPARI, a Fortran compiler is required. CAMPARI uses `make` to build the final executable. However, for ease of use, one can use a local Makefile, i.e. `Makefile.local`, which allows the user to override specific compilation variables within the base `Makefile` such as the compiler flags, optimization flags, etc.

Included in this repository is a `Makefile.local` file which is configured for use with Intel's `ifort` and `mpiifort` Fortran compilers. Alternatively, one can install the requisite GCC version, `gfortran` and related FFTW and XDR library dependencies for compilation. Although `gfortran` produces working binaries, Intel's Fortran compilers produce highly performant binaries compared to other compilers. Where possible, the Intel Fortran compilers should be used.

To build, first navigate to the directory containing all the CAMPARI source files (`*.f90`). Then:

```
mkdir -p ../lib/x86_64/mpi
mkdir -p ../bin/x86_64
```

Next, the user can perform the compilation for `campari_mpi` and `campari` via:

```
make campari_mpi
make campari
```

Once these versions of CAMPARI are built, they can be installed via: `make install campari_mpi` and `make install campari`. If the installation path is to a system directory location, root permissions may be required (i.e. prepend `sudo` to the command).

For more information on building and installation, please consult the official [CAMPARI documentation for V3](http://campari.sourceforge.net/V3/install.html). Although V3 has new features compared to V2, the overall configuration and compilation process is the same.

# Citing

If you use this version of CAMPARI for published scientific research, please cite it. For e.g., a BibTeX entry for the *q-canonical* manuscript is:

```
@article{fossat2019q,
  title={Q-canonical Monte Carlo sampling for modeling the linkage between charge regulation and conformational equilibria of peptides},
  author={Fossat, Martin J and Pappu, Rohit V},
  journal={The Journal of Physical Chemistry B},
  volume={123},
  number={32},
  pages={6952--6967},
  year={2019},
  publisher={ACS Publications}
}
```
