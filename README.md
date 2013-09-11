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

        cmake .
	make
	make test

## Command line options

 Please see the output of

        ./worms -h

## Release Note

 - Release 0.1 (2013/09/12)
     * Initial version

## References

 - O.F. Syljuasen and A.W. Sandvik, "Quantum Monte Carlo with directed loops," Phys. Rev. E 66, 046701 (2002).
 - N. Kawashima and K. Harada, "Recent Developments of World-Line Monte Carlo Methods," J. Phys. Soc. Jpn. 73, 1379 (2003).
 - H. Suwa and S. Todo, "Markov Chain Monte Carlo Method without Detailed Balance," Physical Review Letters 105, 120603 (2010).
 - S. Todo, "Loop Algorithm" in Strongly Correlated Systems: Numerical Methods (Springer Series in Solid-State Sciences), ed. A. Avella and F. Mancini (Springer-Verlag, Berlin, 2013), p. 153.

## Developers

 - Synge Todo (University of Tokyo)
