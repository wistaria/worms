worms: a simple worm code
==========================

This is a sample code describing how the non-reversible worm algorithm
(a.k.a directed-loop algorithm) can be implemented in the continuous
imaginary time path integral representation.  It calculatetes an S=1/2
antiferromagnetic Heisenberg chain in a longitudinal magnetic field.

## Prerequisites

 - Boost C++ Libraries - http://www.boost.org
 - CMake - http://www.cmake.org

## Compilation

    rm -rf worms-build
    mkdir worms-build
    cd worms-build
    cmake [absolute-path-to-worms-source-directory]

You have to specify the path to the Boost include files by using -DBoost_ROOT_DIR option, if you have extracted the Boost in a non-standard place.

    make
    make test

## Command line options

 Please see the output of

    ./worms -h

## Release Note

 - Release 1.0 (2013/09/12)
     * Initial version

 - Release 1.1 (2014/02/24)
     * Fixed generation of transition matrix under the magnetic field
     * Updated logic to optimize number of worms per MCS
     * Added energy measurement
     * Added timer
     * Updated documents
     * Added check scripts and tests

## References

 - O.F. Syljuasen and A.W. Sandvik, "Quantum Monte Carlo with directed loops," Phys. Rev. E 66, 046701 (2002).
 - N. Kawashima and K. Harada, "Recent Developments of World-Line Monte Carlo Methods," J. Phys. Soc. Jpn. 73, 1379 (2003).
 - H. Suwa and S. Todo, "Markov Chain Monte Carlo Method without Detailed Balance," Physical Review Letters 105, 120603 (2010).
 - S. Todo, "Loop Algorithm" in Strongly Correlated Systems: Numerical Methods (Springer Series in Solid-State Sciences), ed. A. Avella and F. Mancini (Springer-Verlag, Berlin, 2013), p. 153.

## Developers

 - Synge Todo (University of Tokyo)
