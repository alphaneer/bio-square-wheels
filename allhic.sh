#!/bin/bash

readonly help_msg="\
Usage: $(basename $0) [OPTION]
`basename $0` is used to run allhic
eg `basename $0` -r ref.fa -b bam -t 20 -k 12 -m 90 -d 2 -R 3  -e AAGCTT
eg `basename $0` -r ref.fa -b bam -t 20 -k 15 -m 1131 -d 22 -R 3  -e GATC  
-b FILE 	bam file can be generated with command:
          	   bwa mem -SP5M -t 20 ref.fa r1.fq.gz r2.fq.gz |samtools view -F 0x10c -Shb - -o align.bam && sambamba markdup --tmpdir=. -t 20 align.bam out.bam
                unreliable alignments in bam are filtered out, skip the bam filter by command:
		   ln -sf out.bam out.REduced.bam
-r FILE		reference fasta file
-e STRING	enzyme_sites HindIII: AAGCTT; MboI: GATC
-k NUMBER       number of groups (user defined K value)
-m INT		minimum number of restriction sites (default, 25) CLUSTER_MIN_RE_SITES
-d INT		CLUSTER_MAX_LINK_DENSITY
-R INT		RATIO
-t INT		max threads
"
parse_args () {
    [[ "$#" == 0 ]] && show_msg "missing file argument" && exit -1;
    local OPTIND=1
    while getopts "r:b:t:k:m:e:d:R:h" opt "$@"; do
        case "$opt" in
            r) ref=$(realpath -s "$OPTARG")
               ;;
            b) bam=$(readlink -f "$OPTARG")
               ;;
	    t) ncpu="$OPTARG"
               ;;
	    k) chr_cnt="$OPTARG"
               ;;
	    m) min_restriction_sites="$OPTARG"
	       ;;
	    d) max_link_density="$OPTARG"
               ;;
	    R) non_informative_ratio=$(echo "$OPTARG"|awk '{print int($1)}')
               ;;
            e) enz="$OPTARG"
               ;;
            h) show_msg
               exit 0
               ;;
            *) show_msg "Try -h for more information."
               exit 1
               ;;
        esac
    done
    shift $((OPTIND-1))
}

show_msg () {
    [[ "$1" ]] && echo "`basename $0`: $1" >&2
    echo "$help_msg" >&2
}

parse_args "$@"

pwd=$(pwd)

ncpu=${ncpu:-20}

min_restriction_sites=${min_restriction_sites:-90}
max_link_density=${max_link_density:-2}
non_informative_ratio=${non_informative_ratio:-3}

enz=${enz:-AAGCTT}

#env
export PATH=$SWPATH/ALLHiC/scripts/:$SWPATH/ALLHiC/bin/:$PATH

s=$(basename $bam .bam)
ln -s $ref .
ref=$(basename $ref)

BED_RE_file=${ref}.near_${enz}.500.bed
cmd="`which make_bed_around_RE_site.pl` $ref $enz 500"
cmd="[ ! -s $BED_RE_file ] && $cmd "
eval $cmd

#UNMAP 0x10 ,MUNMAP 0x01,SECONDARY 0x100 ,DUP 0x400
#uniq mapping : mapping_quality >= 1 and not (unmapped or secondary_alignment) and not ([XA] != null or [SA] != null)
#or -q 30 -h -F 0x40c -@ 20 | awk '{if( (match($0,/NM:i:([0-9]+)/,tA) && tA[1]>5) || (match($0,/XA:Z:.*/) )   ) next;print;   }'
#cmd="sambamba view $bam -L $BED_RE_file -t $ncpu -h -f bam -o ${s}.REduced.bam -F 'mapping_quality >= 30 and [NM]<=5 and \
#not (unmapped or mate_is_unmapped or duplicate or secondary_alignment) and not ([XA] != null or [SA] != null)' "

cmd="sambamba view $bam -L $BED_RE_file -t $ncpu -h -F 'mapping_quality >= 30 and [NM]<=5 and not(unmapped or mate_is_unmapped or secondary_alignment)'|\
awk -v OFS=\"\\t\" '/^@/{print;next;}{if(\$1!=pre) {if(r1!=\"\" && r2!=\"\") print r1\"\\n\"r2;r1=\$0;r2=\"\";pre=\$1;} else {r2=\$0;}   }'|samtools view -Shb -@ $ncpu - -o ${s}.REduced.bam "
cmd="[ ! -s  ${s}.REduced.bam ] && $cmd "
echo $cmd
eval $cmd

cmd="samtools flagstat -@ $ncpu ${s}.REduced.bam > ${s}.REduced.bam.flagstat & "
cmd="[ ! -s ${s}.REduced.bam.flagstat ] && $cmd"
eval $cmd

cleanbam=${s}.REduced.bam

if ! (ls ${s}.REduced.counts_${enz}*g*.txt && ls ${s}.REduced.clm)
then
  cmd="`which allhic` extract $cleanbam $ref --RE $enz"
  echo $cmd
  eval $cmd
  cmd="`which allhic` partition ${s}.REduced.counts_${enz}.txt ${s}.REduced.pairs.txt $chr_cnt --minREs ${min_restriction_sites} --maxLinkDensity ${max_link_density} --nonInformativeRatio ${non_informative_ratio}"
  echo $cmd
  eval $cmd
fi

for x in $(ls ${s}.REduced.counts_${enz}*g*.txt)
do
   cmd="`which allhic` optimize $x ${s}.REduced.clm"
   outfile=$(basename $x .txt)".tour"
   echo $cmd > /dev/stderr
   if ! (test -f $outfile && tail -n2  $outfile |grep '^>FLIP' > /dev/null);    then      echo $cmd;    fi;
done | parallel --retries 2 -j $ncpu

`which ALLHiC_build` $ref

seq_len.pl groups.asm.fasta |head -n $chr_cnt > chrn.list
`which ALLHiC_plot` $cleanbam groups.agp chrn.list 500k,1000k,2000k,2500k pdf
`which statAGP.pl` groups.agp > groups.agp.stat
