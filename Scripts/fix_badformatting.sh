FOLDER="../Data/time_perfomance/k_MBC/k_8/"

FILE1=$FOLDER$1
while read -r line
do
    count=$(echo $line | wc -w)
    if [ $count -ne 1 ]; then
        sed -i "/$line/d" $FILE1
    fi
done < $FILE1 