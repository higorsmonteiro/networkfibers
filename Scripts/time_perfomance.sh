#!/bin/bash
<<COMMENT
    Args:   $1 - Size N of the erdos-renyi network.
            $2 - String to flag which algorithm to run -> 'mbc' or 'ffp'.
            $3 - Number k of networks to test the algorithm.

    Task:   Considering the two main algorithms of the work (FFP and MBC)
            we get the time perfomance (in seconds) of both methods for
            several random network sizes. For each size $1, we perform the
            desired algorithm in k different random networks with size $1
            and p = 0.1 (default).

            The k values of the time perfomance are stored in a file with
            path and format as '../Data/time_perfomance/X/XN$1.dat', where
            X stands either for 'mbc' or 'ffp'.
COMMENT

NETSIZE=$1
MODE=$2
NETNUMBER=$3

MBC_PATH="../Data/time_perfomance/MBC/mbcN$NETSIZE.dat"
FFP_PATH="../Data/time_perfomance/FFP/ffpN$NETSIZE.dat"

for k in $(seq 1 1 $NETNUMBER)
do
    echo here
    if [ $MODE = 'ffp' ]; then
        python ../PyCode/perfomance_comp.py $NETSIZE $MODE >> $FFP_PATH
    elif [ $MODE = 'mbc' ]; then
        python ../PyCode/perfomance_comp.py $NETSIZE $MODE >> $MBC_PATH
    fi
done