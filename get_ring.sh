#!/bin/bash

k=$1;
min_window=$2
cntlist=$3

#####################################################################################
#     ./get_ring.sh 17 23 "3 4 5"
#     ./get_ring.sh 16 20 "4 5"
#
#####################################################################################

################################################commented out########################
((0)) && {
for k in {15,17,19,16};  
do  
  bioawk -c fastx -v k=$k '{l=length($seq)-k+1;for(i=1;i<=l;++i)  {f=r=substr($seq,i,k);r=revcomp(r); if(f>r) o=r;else o=f;print i"\t"o }  }' ref.fa > ref.kmer${k};  
done

for k in {15,16,17}
do
 #filterx -k s -1 'cnt>=1:cut=2' seafood_k17.frq.gz:req=Y:cut=1,2 GRCh38p13_k17.frq.gz bacteria_k17.frq.gz virus_k17.frq.gz|awk -v OFS="\t"  '{for(i=2;i<=NF;++i) if($i=="-") $i=0; $5-=$2;print;  }' > seafood.k17.all.frq 
 cmd="filterx -k s -1 'cnt>=1:cut=2' seafood_k${k}.frq.gz:req=Y:cut=1,2 GRCh38p13_k${k}.frq.gz bacteria_k${k}.frq.gz virus_k${k}.frq.gz |\
 awk -v OFS=\"\\t\"  '{for(i=2;i<=NF;++i) if(\$i==\"-\") \$i=0; \$5-=\$2;print;}' > seafood.k${k}.GRch38.bacteria.virus.frq"
 echo $cmd
done

#just view
k=17;grep -Ff <(awk '{if($3==0 && $4==0 && $5==0) print $1;}' seafood.k${k}.GRch38.bacteria.virus.frq)  \
<(awk '{print "MN908947.3\t"($1-1)"\t"$1"\t"$2"\t"$1}' ref.kmer${k} ) |bedtools merge -i stdin|awk  -v k=$k '{print $0"\t"$3-$2+k-1}'|sort -k4,4nr|less -SN
}

##############################################begin################################
#at leat min_window window to process
grep -Ff <(awk '{if($3==0 && $4==0 && $5==0) print $1;}' seafood.k${k}.GRch38.bacteria.virus.frq) \
 <(awk '{print "MN908947.3\t"($1-1)"\t"$1"\t"$2"\t"$1}' ref.kmer${k} ) |\
bedtools merge -i stdin|awk  -v k=$k '{print $0"\t"$3-$2+k-1}'|\
awk -v OFS="\t" -v min_window=$min_window '$4>=min_window{ 
  cmd="samtools faidx ref.fa "$1":"$2+1"-"$2+$4" 2>/dev/null|seqtk seq|tail -n+2";cmd|getline s;close(cmd); print ++total_cnt,s,$1":"$2"-"$2+$4;
}' > uniq.k${k}.minwindow${min_window}

bioawk -v k=$k  -v OFS="\t" '{windows[$1]=$2;} END{ 
for(x in windows) for(y in windows) 
{ 
  if(x==y) continue;  
  seqs=windows[x]""windows[y]; 
  chkseqs=substr(seqs,length(windows[x])-k+2);len=k-1;
  for(i=1;i<=len;++i) {kmer=rkmer=substr(chkseqs,i,k); rkmer=revcomp(rkmer);if(rkmer<kmer) kmer=rkmer; 
  print kmer,x,y;
} 
} 
}'  uniq.k${k}.minwindow${min_window}  |sortm -k1,1 --temporary-directory=./ --parallel=4 > uniq.k${k}.minwindow${min_window}.chk.tmp

test -s uniq.k${k}.minwindow${min_window}.chk || \
filterx -k s -1 'cnt>=1:cut=2' \
uniq.k${k}.minwindow${min_window}.chk.tmp:grp=1:req=Y:dup=Y:cut=1,2,3 seafood_k${k}.frq.gz GRCh38p13_k${k}.frq.gz  bacteria_k${k}.frq.gz  virus_k${k}.frq.gz|\
awk -v OFS="\t" '{for(i=2;i<=NF;++i) if($i=="-") $i=0; $7-=$4;print;  }' > uniq.k${k}.minwindow${min_window}.chk

awk '{key=$2"\t"$3; if(key=="\t") next; if(! (key in map)) map[key]=1; if(!($5==0 && $6==0 && $7==0)) delete map[key];} END {
 for(x in map) print x"\t"map[x]; 
}' uniq.k${k}.minwindow${min_window}.chk > uniq.k${k}.minwindow${min_window}.pair


for cnt in $cntlist
do
if [ $cnt -eq 3 ]
then
cnt=3
~/bin/awk '{map[$1][$2]=1;} END{
for(i in map) for(j in map[i]) if(j in map) for(k in map[j])  
  if(k in map && i in map[k] && (i!=j && i!=k  && j!=k) ) print i,j,k; 
}' uniq.k${k}.minwindow${min_window}.pair > uniq.k${k}.minwindow${min_window}.pair.out${cnt} &
fi

if [ $cnt -eq 4 ]
then
cnt=4
~/bin/awk '{map[$1][$2]=1;} END{
for(i in map) for(j in map[i]) if(j in map) for(k in map[j]) if(k in map) for(w in map[k]) 
  if(w in map && i in map[w] && (i!=j && i!=k && i!=w && j!=k && j!=w && k!=w) ) print i,j,k,w; 
}' uniq.k${k}.minwindow${min_window}.pair > uniq.k${k}.minwindow${min_window}.pair.out${cnt} &
fi

if [ $cnt -eq 5 ]
then
cnt=5
~/bin/awk '{map[$1][$2]=1;} END{
for(i in map) for(j in map[i]) if(j in map) for(k in map[j]) if(k in map) for(w in map[k]) if(w in map) for(u in map[w]) 
  if(u in map && i in map[u] && (i!=j && i!=k && i!=w && i!=u && j!=k && j!=w && j!=u && k!=w && k!=u && w!=u) ) print i,j,k,w,u; 
}' uniq.k${k}.minwindow${min_window}.pair > uniq.k${k}.minwindow${min_window}.pair.out${cnt} & 
fi

done

wait

for cnt in {3,4,5}
do
file=uniq.k${k}.minwindow${min_window}.pair.out${cnt};
[ ! -f $file ] && continue;
[ $(wc -l $file |awk '{print $1}') -lt 3 ] && continue;
awk -v cnt=$cnt 'NR==FNR{split($3,tA,"[:-]");pos_b[$1]=tA[2];pos_e[$1]=tA[3];seqs[$1]=$2;next;} {
    s=pos_b[$1]"-"pos_e[$1]"\t";for(i=2;i<=cnt;++i) s=s""pos_b[$i]"-"pos_e[$i]"\t"; for(i=1;i<=cnt;++i) s=s""seqs[$i];s=s"\t"$0;print s; 
}' uniq.k${k}.minwindow${min_window} $file | tee ${file}.ovl |awk -v cnt=$cnt '{
    for(x in map) { 
      cir=map[x]; n=split(cir,tA," "); 
      match_cnt=0;for(i=1;i<=n;++i) for(j=1;j<=n;++j) if(tA[i]==$j) match_cnt++;   
      if(match_cnt>=cnt/2) {next;} 
    } 
    print; 
    s=$1;for(i=2;i<=cnt;++i) s=s" "$i; map[++total]=s; 
}' > ${file}.ovl2 & 
done

wait



#############commented out########
((0)) && {
cnt=4;cat uniq.k17.minwindow23.pair.out${cnt}.ovl|awk -v cnt=$cnt '{print length($(cnt+1))"\t"$0;}'|sort -k1,1nr|cut -f2-|awk -v cnt=$cnt '{
    for(x in map) { 
      cir=map[x]; n=split(cir,tA," "); 
      match_cnt=0;for(i=1;i<=n;++i) for(j=1;j<=n;++j) if(tA[i]==$j) match_cnt++;   
      if(match_cnt>=cnt/3) {next;} 
    } 
    print; 
    s=$1;for(i=2;i<=cnt;++i) s=s" "$i; map[++total]=s; 
}'|awk -v cnt=$cnt '{print length($(cnt+1))"\t"$(cnt+1)"\t"$0}'|less -SN
}
