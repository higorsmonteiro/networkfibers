# Fibration building blocks algorithm

Based on the work done by @makselab, the codes here are written in `C` for the main algorithm and in `python` for data analysis and any visualization. The codes are intend to be used for the identification of the minimal fibration building blocks of directed networks. Also, we classify the fibers blocks in a given network according to the `fiber numbers` `< n,l|` presented at the paper of Morone et. al. (2019), to be published.

A more general, and probably more consistent, implementation of the algorithm is presented at the @makselab page, while here I aim to show the application of the algorithm on the `Escherichia Coli` transcription regulatory network data in order to reproduce the results obtained by the Professor Hernan's Lab. However, in my code I apply a different algorithm approach based on the `Paige & Tarjan (1987)` partition refinement algorithm. This approach maintains the runtime complexity of the solution, but it has a simpler implementation and a smaller time prefactors than the algorithm presented by `Cardon and Crushemore (1982)`in which the @makselab code is based.

## Files handling and organization

The folder structure is already given in this repository and all the codes consider this organization in order to have a proper function. (more after)

