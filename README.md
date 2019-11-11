# Identification and Classification of Network Fibers

Based on the work done by [@makselab](https://github.com/makselab), the codes here are intended to be used for the identification of the minimal fibration building blocks of genetic regulatory networks. Also, we classify the fibers blocks in a given network according to the `fiber numbers` `< n,l|` presented at the paper of Morone et. al. (2019), to be published. The codes are written in `C` for the main algorithm and in `python` for relating properly the edgelist data with the correct gene names. I plan, very soon, to implement the same algorithm fully in `python` to give a general implementation framework for these two widely used languages.  

A more general, and probably more consistent, implementation of the algorithm is presented at the [@makselab](https://github.com/makselab) profile page, while here I aim to show the application of the algorithm on the `Escherichia Coli` transcription regulatory network data in order to reproduce the results obtained by the Professor Makse's Lab. However, in my code I apply a different algorithm approach based on the `Paige & Tarjan (1987)` refinement partitioning algorithm. This approach maintains the runtime complexity of the solution, but it has a simpler implementation and smaller time prefactor than the algorithm presented by `Cardon and Crushemore (1982)`in which the Makse's lab code is based.

## Folder Organization

The folder structure is already given in this repository and all the codes consider this organization to have a proper function. The main file `main.c` for the algorithm implementation uses structures and functions given in the personal modules `fibrationf.h`, `structforfiber.h` and `utilsforfiber.h`, where in these modules I have defined the functions to a proper network construction from an edgelist file together with the vizualization functions and, most important, the functions necessary to implement the coarsest refinement network partitioning algorithm based on `input-tree stability` properties.

In `fibrationf.h` I defined two approaches to calculate the `branching ratio n`: by spectral decomposition using a `GSL` library (I need to fix some implementation details yet) and the other based on direct network verification. Depending on the approaches, generating a executable file is different.

For spectral decomposition, uncomment the function and the line `#include <gsl/gsl_eigen.h>` and use the following command lines:

`gcc -O3 -I/usr/local/include -c main.c`  
`gcc -L/usr/local/include main.o -o executablename.out -lgsl -lgslcblas -lm`

For the direct approach, use the single line:  
`gcc -O3 main.c -o executablename.out`  

and run the code:  
`./executablename.out argc[1] argc[2]`

The code receives TWO commands line arguments. The first one is the string identifier for the edgelist file
`ARG1edgelist.dat` containing all the directed links between nodes (3 columns: `"%d\t%d\t%s\n"` -> Regulator Node/ 
Regulated Node/ Type of regulation). For gene regulatory networks, the type of the regulation can be 'positive', 'negative' 
or 'dual'. The second argument is a flag used to signal the code to properly get the gene names of each node number. For 
that, it is necessary an auxiliary file called `ARG1nameID.dat` containing two columns (formatted as `"%s\t%d\n"` -> Gene 
name/ Gene ID number). Thus, if there is a gene name file, the code will properly link all the node numbers with their 
corresponding name if `ARG2` is passed as `-y`, otherwise just the node ID numbers is stored for each node on the network
structure.

## Code optimization and Documentation

The state of the code written is, of course, not finished. Even though the refinement partitioning algorithm is working correctly and implemented with the correct data structures, I'll do some rewritting sooner to guarantee conciseness, code optimization and good documentation for the purpose of good coding practice. Improvement suggestions are always welcome.

The code is free to use under the common `MIT license`.

