#!bin/bash
<<COMMENT
    Args:   $1 - Largest power of two to the random network size.
            $2 - mode flag -> 'ffp', 'mbc' or 'both'
            $3 - Number of edge types.

    Task:   This script is used together with the 'res1_time_perf.sh' 
            script. All the generated networks will have size N equal
            to a power of two, starting with N = 64 = 2^6. The argument
            $1 will give the largest size of q, where N = 2^q.

            Ex: q = 13 -> N_max = 8192; q = 14 -> N_max = 16384.
COMMENT

POWER=$1
MODE=$2
EDGE_TYPE=$3
if [ $POWER -ge 6 ]; then
    for k in $(seq 6 1 $POWER)
    do
        N=$(( 2**$k ))
        if [ $MODE = 'ffp' ]; then
            for k_aver in 1 2 4 8
            do
                echo Getting data - size $N and average degree $k_aver. mode FFP.
                sh res2_time_perf.sh $N ffp $k_aver 30 $EDGE_TYPE
            done
        elif [ $MODE = 'mbc' ]; then
            for k_aver in 1 2 4 8
            do
                echo Getting data - size $N and average degree $k_aver. mode MBC.
                sh res2_time_perf.sh $N mbc $k_aver 30 $EDGE_TYPE
            done
        elif [ $MODE = 'both' ]; then
            for k_aver in 1 2 4 8
            do
                echo Getting data - size $N and average degree $k_aver. mode FFP.
                sh res2_time_perf.sh $N ffp $k_aver 30 $EDGE_TYPE
                echo Getting data - size $N and average degree $k_aver. mode MBC.
                sh res2_time_perf.sh $N mbc $k_aver 30 $EDGE_TYPE
            done
        fi
    done
fi