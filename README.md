# Unfitted HHO

## What is this

This is a prototype for the Unfitted HHO method (see [this paper](https://hal.archives-ouvertes.fr/hal-01625421)).
At the moment it implements the 2D version of the Fictitious Domain problem and fixes bad cuts using the Point Displacement Algorithm.

## Building the code

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

If you want to just run the code, on OS X, you can tap my homebrew repo with

	brew tap datafl4sh/code

and install with just
	
	brew install cuthho

## Running the thing
You will end up with two executables, `cuthho_square` and `convergence_test`.
The first is the cutHHO driver, the second does a convergence test of the standard cutHHO method.

The executable `cuthho_square` accepts the following parameters:
 * `-M`: the number of cells in x direction
 * `-N`: the number of cells in y direction
 * `-k`: the degree of the method
 * `-m`: disable Point Displacement

Be sure to run `cuthho_square` in a specific directory, because it produces a bunch of files.

