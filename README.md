# Unfitted HHO

### What is this

This is a prototype for the Unfitted HHO method (see [this paper](https://hal.archives-ouvertes.fr/hal-01625421)).
At the moment it implements the 2D version of the Fictitious Domain problem and fixes bad cuts using the Point Displacement Algorithm.

The code is tested on Mac OS X and on Linux (sorry, I don't have the time to deal with
the non-compliance to standards of Windows). To compile you need, apart of a decent
C++14 compiler that should be already available in your system, the following software:

 * CMake
 	* `brew install cmake` on OS X
 	* `apt-get install cmake` on Linux
 * SILO
 	* `brew install silo` on OS X
 	* `apt-get install libsilo-dev` on Linux
 * Lua
 	* `brew install lua` on OS X
 	* `apt-get install liblua5.3-dev` on Linux
 * Eigen3
 	* `brew install eigen` on OS X
 	* `apt-get install libeigen3-dev` on Linux


