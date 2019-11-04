# Network fibers identification algorithm

Based on the work done by @makselab, the codes here are written in `C` for the main algorithm and in `python` for data analysis and any visualization. The codes are intend to be used for the identification of the minimal fibration building blocks of directed networks. Also, we classify the fibers blocks in a given network according to the `fiber numbers` `< n,l|` presented at the paper of Morone et. al. (2019), to be published.

A more general, and probably more consistent, implementation of the algorithm is presented at the @makselab page, while here I aim to show the application of the algorithm on the `Escherichia Coli` transcription regulatory network data in order to reproduce the results obtained by the Professor Hernan's Lab. However, in my code I apply a different algorithm approach based on the `Paige & Tarjan (1987)` partition refinement algorithm. This approach maintains the runtime complexity of the solution, but it has a simpler implementation and a smaller time prefactors than the algorithm presented by `Cardon and Crushemore (1982)`in which the @makselab code is based.

## Files handling and organization

The folder structure is already given in this repository and all the codes consider this organization in order to have a proper function. The main file for the algorithm implementation is the `minimalfiber.c`. The code uses handmade structures and functions given in the personal modules `fibrationf.h`, `structforfiber.h` and `utilsforfiber.h`, where in these modules I have defined the functions to a proper network construction from an edgelist file together with the vizualization functions and mainly the functions necessary to implement the coarsest refinement network partitioning algorithm. Since we use a `GSL` library, to create an executable file for the `minimalfiber.c` we use the following command lines:

`gcc -O3 -I/usr/local/include -c minimalfiber.c`
`gcc -L/usr/local/include minimalfiber.o -o executablename.out -lgsl -lgslcblas -lm`

and the line
`./executablename.out argc[1] argc[2]`

to run the code.

As we see, the code receives TWO commands line arguments. The first one is the string identifier for the edgelist file 
('ARG1edgelist.dat') containing all the directed links between nodes (3 columns: "%d\t%d\t%s\n" -> Pointing Node/ 
Pointed Node/ Type of regulation). For gene regulatory networks, the type of the regulation can be 'positive', 'negative' 
or 'dual'. The second argument is a flag used to signal the code to properly get the gene names of each node number. For 
that, it is necessary an auxiliary file called 'ARG1nameID.dat' containing two columns (formatted as "%s\t%d\n" -> Gene 
name/ Gene ID number). Thus, if there is a gene name file, the code will properly link all the node numbers with their 
corresponding name if 'ARG2' is passed as '-y', otherwise just the node numbers is stored for each node.

