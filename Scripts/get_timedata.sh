#!bin/bash
<<COMMENT
    Args:   $1 - Largest power of two to the random network size.
            $2 - Expected mean degree <k> of the random network.

    Task:   This script is used together with the 'time_perfomance.sh' 
            script. All the generated networks will have size N equal
            to a power of two, starting with N = 64 = 2^6. The argument
            $1 will give the largest size of q, where N = 2^q.
COMMENT

div ()  # Arguments: dividend and divisor
{
        if [ $2 -eq 0 ]; then echo division by 0; exit; fi
        local p=12                            # precision
        local c=${c:-0}                       # precision counter
        local d=.                             # decimal separator
        local r=$(($1/$2)); echo -n $r        # result of division
        local m=$(($r*$2))
        [ $c -eq 0 ] && [ $m -ne $1 ] && echo -n $d
        [ $1 -eq $m ] || [ $c -eq $p ] && return
        local e=$(($1-$m))
        let c=c+1
        div $(($e*10)) $2
}

POWER=$1
K_MEAN=$2
if [ $POWER -ge 6 ]; then
    for k in $(seq 6 1 $POWER)
    do
        N=$(( 2**$k ))
        PROB=$(div $K_MEAN $N)
        echo Getting data - size $N and p $PROB. mode ffp
        sh time_perfomance.sh $N ffp $PROB 20 $K_MEAN
        echo Getting data - size $N and p $PROB. mode mbc
        sh time_perfomance.sh $N mbc $PROB 20 $K_MEAN
    done
fi