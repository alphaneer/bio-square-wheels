#!/bin/bash
#bgzip tabix seqtk samtools
readonly help_msg="\
Usage: $(basename $0) [OPTION]
`basename $0` is used to run assembly tools
eg `basename $0` -t 0 -i reads.fa -I '{}' -x pacbio|ont -g genomesize(Mb);  
tools and default options:
0    SMARTdenovo
1    Wtdbg2
2    Canu
3    Flye
4    Ra
5    Unicycler
"
parse_args () {
    [[ "$#" == 0 ]] && show_msg "missing file argument" && exit -1;
    local OPTIND=1
    while getopts "i:I:t:x:p:hg:" opt "$@"; do
        case "$opt" in
            t) tools="$OPTARG"
               ;;
            i) reads=$(readlink -f "$OPTARG")
               ;;
            x) reads_type="$OPTARG"
               ;;
            p) ncpu="$OPTARG"
               ;;
            g) genome_size="$OPTARG"
               ;;
            h) show_msg
               exit 0
               ;;
            I) other_options="$OPTARG"
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

TOOLSBASE=/public/agis/ruanjue_group/lialun/sw/assembly/bin
pwd=$(pwd)

cmd="export PATH=$TOOLSBASE:\$PATH"
echo $cmd >&2
echo "use ~/bin/fastq-dump -Z -M 0 --table SEQUENCE \$SRA 2> err |seqtk seq -A -L 1 | pigz > reads.fa.gz" >&2
#~/bin/awk -v OFS="\t"  '{map[FILENAME][$1]=$2" "$3;} END{for(x in map) {s=x;gsub("/.*","",s);print s,map[x]["Total:"]/1000000,map[x]["Count:"],map[x]["N50:"];}  }' */assembly.fasta.seq_n50


#---------------------------smartdenovo--------------------
mkdir -p $pwd/smartdenovo;cd $pwd/smartdenovo
tmp="smartdenovo.pl -c 1 -t $ncpu $reads > wtasm.mak && make -f wtasm.mak && ln -sf wtasm.dmo.cns assembly.fasta && s50 assembly.fasta";
printf "runit.pl \"%s\" > runit.log 2>&1 \n" "$tmp" > work.sh
chmod u+x work.sh
tmp="cd $pwd/smartdenovo;pstats.sh $(pwd)/work.sh > pstats.log 2>&1 "
[[ $tools -eq 0 ]] && cmd="$tmp"
[[ $tools -eq -1 ]] && allcmd="$tmp"


#---------------------------wtdbg-------------------------
tmp=""
mkdir -p $pwd/wtdbg;cd $pwd/wtdbg

if [[ $reads_type == "pacbio" ]] #-p 19 -AS 2 -s 0.05 -L 5000
then
tmp="wtdbg2 -t $ncpu -x preset1 -i $reads -fo dbg"
tmp="$tmp && wtpoa-cns -t $ncpu -i dbg.ctg.lay.gz -fo dbg.ctg.slf.fa"
tmp="$tmp && minimap2 -t $ncpu -x map-pb -a dbg.ctg.slf.fa $reads | best_sam_hits4longreads | samtools sort -@ $ncpu >dbg.ctg.map.srt.bam";
fi

if [[ $reads_type == "ont" ]] 
then
tmp="wtdbg2 -t 0 -x nanopore -i $reads -fo dbg"
tmp="$tmp && wtpoa-cns -t $ncpu -i dbg.ctg.lay.gz -fo dbg.ctg.slf.fa"
tmp="$tmp && minimap2 -t $ncpu -x map-ont -a dbg.ctg.slf.fa $reads | best_sam_hits4longreads | samtools sort -@ $ncpu >dbg.ctg.map.srt.bam";
fi

tmp="$tmp;samtools view -@ $ncpu dbg.ctg.map.srt.bam |wtpoa-cns -t $ncpu -d dbg.ctg.slf.fa -i - -fo dbg.ctg.lrp.fa && ln -sf dbg.ctg.lrp.fa assembly.fasta && s50 assembly.fasta"

printf "runit.pl \"%s\" > runit.log 2>&1 \n" "$tmp" > work.sh
chmod u+x work.sh
tmp="cd $pwd/wtdbg; pstats.sh $(pwd)/work.sh > pstats.log 2>&1 "
[[ $tools -eq 1 ]] && cmd="$tmp"
[[ $tools -eq -1 ]] && allcmd="$allcmd\n$tmp"

#---------------------------canu-------------------------
mkdir -p $pwd/canu;cd $pwd/canu
tmp="canu -p canu -d out_dir useGrid=false minThreads=$ncpu maxThreads=$ncpu  genomeSize=${genome_size}m "
[[ $reads_type == "pacbio" ]] &&  tmp="$tmp -pacbio-raw $reads "
[[ $reads_type == "ont" ]] &&  tmp="$tmp -nanopore-raw $reads "
tmp="$tmp && ln -sf out_dir/canu.contigs.fasta assembly.fasta && s50 assembly.fasta "
printf "runit.pl \"%s\" > runit.log 2>&1 \n" "$tmp" > work.sh
chmod u+x work.sh
tmp="cd $pwd/canu; pstats.sh $(pwd)/work.sh > pstats.log 2>&1 "
[[ $tools -eq 2 ]] && cmd="$tmp"
[[ $tools -eq -1 ]] && allcmd="$allcmd\n$tmp"

#---------------------------flye-------------------------
mkdir -p $pwd/flye;cd $pwd/flye
tmp="flye  -g ${genome_size}m -o out_dir -t $ncpu "
[[ $reads_type == "pacbio" ]] &&  tmp="$tmp --pacbio-raw $reads "
[[ $reads_type == "ont" ]] &&  tmp="$tmp --nano-raw $reads "
tmp="$tmp && ln -sf out_dir/scaffolds.fasta assembly.fasta && s50 assembly.fasta "
printf "runit.pl \"%s\" > runit.log 2>&1 \n" "$tmp" > work.sh
chmod u+x work.sh
tmp="cd $pwd/flye; pstats.sh $(pwd)/work.sh > pstats.log 2>&1 "
[[ $tools -eq 3 ]] && cmd="$tmp"
[[ $tools -eq -1 ]] && allcmd="$allcmd\n$tmp"

#---------------------------ra-------------------------
mkdir -p $pwd/ra;cd $pwd/ra
[[ $reads_type == "pacbio" ]] &&  tmp="ra -x pb "
[[ $reads_type == "ont" ]] &&  tmp="ra -x ont "
tmp="$tmp  -t $ncpu $reads > assembly.fasta && s50 assembly.fasta"
printf "runit.pl \"%s\" > runit.log 2>&1 \n" "$tmp" > work.sh
chmod u+x work.sh
tmp="cd $pwd/ra; pstats.sh $(pwd)/work.sh > pstats.log 2>&1 "
[[ $tools -eq 4 ]] && cmd="$tmp"
[[ $tools -eq -1 ]] && allcmd="$allcmd\n$tmp"
 s
#---------------------------unicycle-------------------------
#tmp="mkdir -p $pwd/unicycler;cd $pwd/unicycler"
#tmp="$tmp && unicycler -l $reads -o out_dir -t 8"
#[[ $tools -eq 5 ]] && cmd=$tmp
#[[ $tools -eq -1 ]] && allcmd="$allcmd\n$tmp"


[[ $tools -eq -1 ]] && echo -e $allcmd || echo $cmd
