Hi, this is the readme to testnest, a simple Python demonstration of MultiNest. I hope you find it helpful!

-------------------------------------------------------------------------------------------------

Acknowledgements:

I'd like to point out that this depends on the wonderful Bayesian fitting package MultiNest (and especially the PyMultiNest wrapper) and therefore lots of credit goes to the authors of those packages and you should cite them on any work using their code.


-------------------------------------------------------------------------------------------------

Installation documentation:

This can be a world of pain.

First go to http://ccpforge.cse.rl.ac.uk/gf/project/multinest/, sign up and download multinest v 2.18, which is the version I've used in putting together pymask. I'm sure other versions will work, but there could be dragons.

There's a great blog about how to install MultiNest at http://www.astro.gla.ac.uk/~matthew/blog/?p=342. Personally, I opted for the gfortran install rather than the Intel compilers, which are probably faster but seemed like a lot of effort. For my gfortran install, I changed the top of the Makefile in the main directory to read as:

#FC = mpif90 -DMPI
FC = gfortran
CC = mpicc
CXX = mpiCC
#FFLAGS +=  -w -O3
FFLAGS +=  -w -O3 -ffree-line-length-none -fPIC
CFLAGS += -O3 -DMPI -fPIC

LAPACKLIB = -llapack

NESTLIBDIR = ./

export FC CC CXX FFLAGS CFLAGS LAPACKLIB

Importantly, doing this recquires gfortran to have the appropriate libs so you need to run sudo apt-get install gfortran libblas-dev liblapack-dev. 

Anyway, that should about make it work, and if not, check the blog. 

The next issue is getting PyMultiNest to work. I went to https://github.com/JohannesBuchner/PyMultiNest and downloaded it, and did python setup.py install from the command line. Follow the installation instructions at http://johannesbuchner.github.io/PyMultiNest/install.html#installing-the-python-library as follows:

1) add to .bashrc 
export MULTINEST=/my/multinest/directory

2) go to that directory and type at the command line: make libnest3.so - this creates a dynamic library, so you have to compile it on your own machine as it has relative directory structures.

3) go to the pymultinest directory and type at the command line: make -C multinest_bridge libcnest.so

4) go to .bashrc and add 
export LD_LIBRARY_PATH=$MULTINEST:/my/pymultinest/directory/multinest_bridge

5) test the libraries as 

$ python -c 'import pymultinest'
$ python -c 'import pyapemost'

in the appropriate folders.

6) IMPORTANT: before you run any MultiNest code, you have to create a folder chains/ in your working directory - this is where the results get stored!

I wish you luck!
