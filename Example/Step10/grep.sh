for i in `cat $1`
do
  grep -w -m 1 $i $2
done
