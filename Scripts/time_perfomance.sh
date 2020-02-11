#!/bin/bash
<<COMMENT
    Args:   $1 - Size N of the erdos-renyi network.
            $2 - String to flag which algorithm to run -> 'mbc' or 'ffp'.
            $3 - p parameter of G(n, p) model.
            $4 - Number k of networks to test the algorithm.
            $5 - Expected mean degree c (integer).

    Task:   Considering the two main algorithms of the work (FFP and MBC)
            we get the time perfomance (in seconds) of both methods for
            several random network sizes. For each size $1, we perform the
            desired algorithm in k different random networks with size $1
            and p = $3.

            The k values of the time perfomance are stored in a file with
            path and format as '../Data/time_perfomance/X/k_$5/N$1.dat', where
            X stands either for 'k_MBC' or 'k_FFP'.
COMMENT

NETSIZE=$1
MODE=$2
PROB=$3
NETNUMBER=$4
K_MEAN=$5

MBC_PATH="../Data/time_perfomance/k_MBC/k_$K_MEAN/N$NETSIZE.dat"
FFP_PATH="../Data/time_perfomance/k_FFP/k_$K_MEAN/N$NETSIZE.dat"

for k in $(seq 1 1 $NETNUMBER)
do
    echo sample $k.
    if [ $MODE = 'ffp' ]; then
        python ../PyCode/perfomance_comp.py $NETSIZE $PROB $MODE >> $FFP_PATH
    elif [ $MODE = 'mbc' ]; then
        python ../PyCode/perfomance_comp.py $NETSIZE $PROB $MODE >> $MBC_PATH
    fi
done