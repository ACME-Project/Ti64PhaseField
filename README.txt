README.txt

--- 0. Dependencies ---

This makefile directly depends on the following packages : (unless otherwise specified, the packages are common for both ubuntu and debian) :

c compiler : g++
p compiler : libopenmpi-dev , openmpi-common  
q compiler : Load the xl module after logging into q . (module load xl)

If the error messages still persist,  do a sudo apt-get update and then try make.

=======
This file will explain the use of the Toth XMPF as encoded by David
Crist in C++ using the Mesoscale Microstructure Simulation Package
(MMSP).

--- 0.A. MMSP Installation ---

Download MMSP from its central github repository,
https://github.com/mesoscale/mmsp

Navigate to mmsp/utility/ and build all the utility files by utilizing
the Makefile.

Add a line to your ~/.bashrc file defining a new system variable
$MMSP_PATH by adding an export line to the file, inserting the
appropriate directory where MMSP is stored.

export MMSP_PATH=/home/user_name/direcctory/mmsp

--- 0.B.  OpenMPI Installation ---

You will need to install the packages for OpenMPI, an extensive
library of Message Passing Interface protocols for parallel 
processing.

This can be installed on Debian by using 

sudo apt-get install package

for the following OpenMPI packages:
openmpi-common
openmpi-bin
openmpi-doc
libopenmpi-dev
>>>>>>> b3531ea6a7bebcb17c3e11f72abf406b5a207a28

--- 1. General Use ---

As this simulation software was written using MMSP, it follows the
conventions for generation and execution prescribed by MMSP.

There are two options for building the executable using the Makefile;
sparse and bgq. Sparse uses a sparse representation of field data for
execution on either singular or multiple processors in a generalized
computing environment, while bgq includes libraries specifically
required for building and operation of the code on the AMOS Blue
Gene/Q system. Once the appropriate executable is created, it can be
run using MMSP conventions such as:

./sparse.out --example 2 planar.00000.dat

This will tell the executable sparse.out to generate a 2-dimensional
initial condition microstructure patterned on the planar interface
condition, and to output the data in an MMSP binary data file named
planar.00000.dat

The simulation can be run by executing the following command

./sparse.out planar.00000.dat 25000 250 

which executes the physics modeled in the sp-xmpf.cpp file (from which
sparse.out is built) on the microstructure specified in the MMSP
binary data file split.00000.dat for 25,000 steps, outputing an MMSP
binary data file representation of the microstructure every 250 steps.

This will generate a list of files such as 

planar.00250.dat
planar.00500.dat
planar.00750.dat
planar.01000.dat
planar.01250.dat
...

This process may be sped up by taking advantage of the Message Passing
Interface (MPI) which is natively supported by MMSP and properly
handled in the coding of sp-xmpf.cpp to allow multiprocessing
optimization such as on Ridcully. Thi can be done by executing a
command like

mpirun -np 8 ./sparse.out planar.00000.dat 25000 250

This line performs the same operation as the initial execution line,
but does so using MPI to run it on 8 processors.

--- 2. Data Analysis ---

Interpretation of MMSP binary output datafiles can be done using
tool.cpp included in the utilities directory, and buildable with the
included Makefile.

Use of tool is fairly straightforward: 

./tool file.00000.dat output_directory/

the executable is called, followed by the file to be analyzed, and the
directory to store the analysis output files in. This will generate
two files: output_directory/file.00000_grains.png, showing a grain
boundary outline of the diffuse interfaces and
output_directoty/file.00000_areas.csv, a list of the area of every
grain in the microstructure

tool.cpp can be edited to perform all varieties of operations and
analyses on the MMSP grid data.

This performance can be done for all MMSP files in a directory that
start with the header file.:

for f in file.*; do tool $f output_director; done

--- 3. Operation on AMOS ---

Reminder that execution of commands on AMOS must be done on the
primary front-end node q.

Because of the queue manager on AMOS, execution of code is best done
in an sbatch file. The text of an example sbatch files is shown below:

#!/bin/sh

#SBATCH --job-name=start_ex
#SBATCH --account=GGST
#SBATCH --mail-user=cristd2@rpi.edu
#SBATCH --mail-type=END

#SBATCH --partition medium
#SBATCH --time 720
#SBATCH --nodes 512
#SBATCH --ntasks 8192
#SBATCH --overcommit

#SBATCH -o /gpfs/u/home/GGST/GGSTcrdv/scratch/example/directory/start_ex.log
#SBATCH -D /gpfs/u/home/GGST/GGSTcrdv/scratch/example/directory

#srun --unbuffered /gpfs/u/home/GGST/GGSTcrdv/scratch/example/directory/q_sparse.out --example 2 /gpfs/u/home/GGST/GGSTcrdv/scratch/example/directory/file.00000.dat

Make note to assign directories and files correctly, and to always
output files into the ~/scratch/ directory on AMOS.

As a reminder, AMOS operates in a big-endian binary data format, and
output files generated on AMOS must be converted appropriately to be
analyzed on little-endian x86 systems.

This can be performed using the MMSP default utility wrongendian,
which must be built in the mmsp/utility/ directory, such as

wrongendian file.a.dat file.b.dat

which converts file.a.dat and outputs the endian-swapped file to
file.b.dat. When converting endianness of multiple files, it is often
more convenient to create a directory to send all converted files to
with their original naming convention such as

mkdir converted
wrongendian file.00000.dat converted/

This can be done with a command-line bash script for all MMSP files in
a directory that start with the header file.:

for f in file.*; do wrongendian $f converted/$f; done
