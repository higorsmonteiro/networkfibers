#!bin/bash
<<COMMENT
    Args:   $1 - Largest power of two to the random network size.
            $2 - Expected mean degree <k> of the random network.
            $3 - mode flag -> 'ffp', 'mbc' or 'both'

    Task:   This script is used together with the 'res1_time_perf.sh' 
            script. All the generated networks will have size N equal
            to a power of two, starting with N = 64 = 2^6. The argument
            $1 will give the largest size of q, where N = 2^q.

            Ex: q = 13 -> N_max = 8192; q = 14 -> N_max = 16384.
COMMENT

POWER=$1
K_MEAN=$2
MODE=$3
if [ $POWER -ge 6 ]; then
    for k in $(seq 6 1 $POWER)
    do
        N=$(( 2**$k ))
        if [ $MODE = 'ffp' ]; then
            echo Getting data - size $N and average degree $K_MEAN. mode FFP.
            sh res1_time_perf.sh $N ffp $K_MEAN 30 1
        elif [ $MODE = 'mbc' ]; then
            echo Getting data - size $N and average degree $K_MEAN. mode MBC.
            sh res1_time_perf.sh $N mbc $K_MEAN 30 1
        elif [ $MODE = 'both' ]; then
            echo Getting data - size $N and average degree $K_MEAN. mode FFP.
            sh res1_time_perf.sh $N ffp $K_MEAN 30 1
            echo Getting data - size $N and average degree $K_MEAN. mode MBC.
            sh res1_time_perf.sh $N mbc $K_MEAN 30 1
        fi
    done
fi