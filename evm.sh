#!/bin/bash
#bgzip tabix samtools
EVMHOME=$SWPATH/EVidenceModeler
readonly help_msg="\
Usage: $(basename $0) [OPTION]
split the genome into contigs and run it one after one.It's always good to print all commands (using -D) and run in parallel
`basename $0` is used to run EVM
`basename $0` -r ref.hardmask.fa -a augustus.gff:1 -p pro.gff3:5 -t pasa_assemblies.gff3:10  -D  > work.sh;
or `basename $0` -r ref.softmask.fa -a augustus.gff:1 -p pro.gff3:5 -t pasa_assemblies.gff3:10  -D -R auto > work.sh;  
sample weight file

PROTEIN nap-nr_minus_rice.fasta 1
PROTEIN genewise-nr_minus_rice.fasta    5
TRANSCRIPT      gap2-plant_gene_index.11282006.fasta    1
TRANSCRIPT      alignAssembly-rice_release4_gmapsim4_02152006   10
ABINITIO_PREDICTION     fgenesh 1
ABINITIO_PREDICTION     genemark        1
ABINITIO_PREDICTION glimmerHMM 1

Options:
  -r FILE       reference genome,soft mask
  -a FILE:INT   ab_init prediction gff3 file:weight of the ab_init
  -p FILE:INT   protein prediction gff3 file:weight of the protein file
  -t FILE:INT   transcripts gff3 file:weight of the transcript
  -D            dry-run
  -w FILE       weight file, use auto to infer from the softmask reference genome
  -h            display this help and exit
"
parse_args () {
    [[ "$#" == 0 ]] && show_msg "missing file argument" && exit -1;
    local OPTIND=1
    while getopts "r:a:p:t:R:w:hD" opt "$@"; do
        case "$opt" in
            r) genome=$(readlink -f "$OPTARG") #mask
               ;;
            R) repeats_file=$(readlink -f "$OPTARG") #mask
               ;;
            a) tmp="$OPTARG";arrtmp=(${tmp//:/ });arrtmp[0]=$(readlink -f ${arrtmp[0]});
               abinit[0]=${arrtmp[0]};abinit[1]=${arrtmp[1]};
               ;;
            p) tmp="$OPTARG";arrtmp=(${tmp//:/ });arrtmp[0]=$(readlink -f ${arrtmp[0]});
               protein[0]=${arrtmp[0]};protein[1]=${arrtmp[1]};
               ;;
            t) tmp="$OPTARG";arrtmp=(${tmp//:/ });arrtmp[0]=$(readlink -f ${arrtmp[0]});
               transcript[0]=${arrtmp[0]};transcript[1]=${arrtmp[1]};
               ;;
            h) show_msg
               exit 0
               ;;
            D) dry_run=1
               ;;
            w) weight_file=$(readlink -f "$OPTARG")
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


#make genome index
[[ -f ${genome}.fai ]] || samtools faidx $genome

#make repeat gff
if [ "$repeats_file" == "auto" ] ;then
  seqtk seq $genome | ~/bin/awk '{if(substr($1,1,1)==">") {chr=substr($1,2);next;}  n=split($1,tA,/[acgtn]+/,tB); nextpos=1;for(i=1;i<=n;++i) { if(length(tA[i])>0) print chr"\t"nextpos-1"\t"n
extpos+length(tA[i])-1;  nextpos+=length(tA[i])+length(tB[i]);} }' |sort -k1,1 -k2,2n |bedtools complement -i stdin -g <(sort -k1,1 -k2,2n ${genome}.fai) > repeat.bed
  awk -v OFS="\t" '{print $1,"repeat",".",$2+1,$3,".","+",".","ID="NR}' repeat.bed > repeat.gff  
  repeats_file="$pwd/repeat.gff"
else 
  repeats_file=$(realpath $repeats_file)
fi


if [ -z "$weight_file" ] ;then
        cmd="$EVMHOME/EvmUtils/create_weights_file.pl -A ${abinit[0]} "
        [[ -z "${protein[0]}" ]] || cmd="$cmd  -P ${protein[0]} "
        [[ -z ${transcript[0]} ]] || cmd="$cmd -T ${transcript[0]} "
        cmd="$cmd > weights.txt"
        eval $cmd
        awk -v w_abinit=${abinit[1]} -v w_protein=${protein[1]}  -v w_transcript=${transcript[1]} -v OFS="\t" '{\
 if($1=="ABINITIO_PREDICTION") $3=w_abinit; else if($1=="PROTEIN") $3=w_protein; else if($1=="TRANSCRIPT") $3=w_transcript;if($3+0<1) $3=1;print;}' weights.txt > tmp.weights.txt
        mv tmp.weights.txt weights.txt
        weight_file=$pwd/weights.txt
fi

cut -f1,2 ${genome}.fai|while read chr len 
do
  tmp=($x)  
  tmpdir="mktemp -dp $pwd";
  cmd="tmpdir=\$($tmpdir);trap 'rm -rf "\$tmpdir"' EXIT; cd \$tmpdir";
  cmd="$cmd ; samtools faidx $genome '$chr' > genome.fasta"
  cmd="$cmd ; $EVMHOME/EvmUtils/gff_range_retriever.pl $chr 1 $len < $repeats_file > repeats.gff3";
  cmd="$cmd ; $EVMHOME/EvmUtils/gff_range_retriever.pl $chr 1 $len < ${abinit[0]} > gene_predictions.gff3";
  [[ -z "${protein[0]}" ]] || cmd="$cmd ; $EVMHOME/EvmUtils/gff_range_retriever.pl $chr 1 $len < ${protein[0]} > protein_alignments.gff3";
  [[ -z "${transcript[0]}" ]] || cmd="$cmd ; $EVMHOME/EvmUtils/gff_range_retriever.pl $chr 1 $len < ${transcript[0]} > transcript_alignments.gff3";
  cmd="$cmd && $EVMHOME/evidence_modeler.pl --genome genome.fasta --gene_predictions gene_predictions.gff3 --weights $weight_file "
  [[ -z ${transcript[0]} ]] || cmd="$cmd --transcript_alignments transcript_alignments.gff3"
  [[ -z ${protein[0]} ]] || cmd="$cmd --protein_alignments protein_alignments.gff3"  
  [[ -z $repeats_file ]] || cmd="$cmd --repeats repeats.gff3"
  cmd="$cmd > evm.out "
  cmd="$cmd; $EVMHOME/EvmUtils/EVM_to_GFF3.pl evm.out $chr > $pwd/${chr}.evm.gff3"
  cmd="$cmd && rm -rf \$tmpdir"
  [[ $dry_run -eq 1 ]] && (echo $cmd) || (eval $cmd)
  #echo $cmd
done
