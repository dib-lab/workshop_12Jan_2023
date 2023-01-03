inputBAM=$1
haploid1=$2
haploid2=$3


rm -f mkfifo $haploid1.fifo $haploid2.fifo

mkfifo $haploid1.fifo $haploid2.fifo


cat <(samtools view -H $inputBAM) <(grep "HP:i:1" $haploid1.fifo) |samtools view -b -@8 -  >  $haploid1 &
cat <(samtools view -H $inputBAM) <(grep "HP:i:2" $haploid2.fifo) |samtools view -b -@8 -  >  $haploid2 &


samtools view $inputBAM | tee  $haploid2.fifo >  $haploid1.fifo

samtools index $haploid1 &
samtools index $haploid2 &

wait

rm -f $haploid1.fifo $haploid2.fifo
