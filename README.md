# bio-square-wheels
square wheels used inbioinformatics
for i in $(seq 1 10)
do
 ld_decay.sh -i sample.vcf.gz -p CSR,CSA,CSS -c 8 -t 15 -D -o PopLDdecay${i}
done
