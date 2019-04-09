#!/bin/bash
indir=$(realpath -s $1)
genome=$(realpath -s $2)

novel_patten=$3

[[ -f ${genome}.fai ]] || samtools faidx $genome
genomefai=${genome}.fai
#get N pos
if [ ! -f "ref.N.pos.bga" ]; then 
seqtk seq $genome|~/bin/awk '{if(substr($1,1,1)==">") {chr=substr($1,2);next;} n=split($1,tA,/[nN]+/,tB);nextpos=1;for(i=1;i<=n;++i) { if(length(tA[i])>0) print chr"\t"nextpos-1"\t"nextpos+length(tA[i])-1;  nextpos+=length(tA[i])+length(tB[i]);} }' |sort -k1,1 -k2,2n |bedtools complement -i stdin -g <(sort -k1,1 -k2,2n ${genome}.fai) |awk '{print $0"\tNA"}' > ref.N.pos.bga
fi
ref_Npos=$(realpath -s "ref.N.pos.bga

function stats_bga() { bga=$1;[ $bga=="-" ] && bga="stdin"; Npos="$2";bedtools unionbedg -i $bga $Npos|grep -v 'NA$'|~/bin/awk -v stringf="0.05,1,2.5,5,25,50,75,95,97.5,99,100"  'BEGIN{PROCINFO["sorted_in"]="@ind_num_asc"; split(stringf,f,",");} {tmp=$3-$2;map[$4]+=tmp;s+=tmp;} END{pos=1; for(x in map) { for(i in f) { fi=f[i]*s/100; if(fi>= pos && fi<pos+map[x]) printf("%s%\t%d\n",f[i],x);} pos+=map[x];}}';  }

cd $indir
bedfiles=$(ls *.bed.gz|awk '{s=s" <(unpigz -dc "$1")"} END{print s}')
cmd="sort -k1,1 -k2,2n -m $bedfiles |bedtools genomecov -bga -g $genomefai -i stdin|sed -n '/$novel_patten/!{p;d}; w novel.bga'  |tee >(stats_bga stdin $ref_Npos > GRch38.bga.stats) |pigz -9 >  GRch38.bga.gz && pigz -f9 novel.bga"
echo $cmd
eval $cmd
