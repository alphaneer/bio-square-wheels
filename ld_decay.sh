#!/bin/bash
set -eo pipefail
readonly usage_msg="Try ' -h' for more information."
readonly help_msg="\
Usage: `basename $0` [OPTION] eg. `basename $0` -i sample.vcf.gz - pop1,pop2,pop3 -c 8 -t 3 -O outdir_test1

LD decay was calculated using PopLDdecay(v3.40) by measuring the \
correlation coefficients (r2) for all pairs of SNPs within 1000 kb
with the parameters '-MaxDist 1000'
the result was plotted using the Perl scripts provided in PopLDdecay

An equal number of individuals were resampled from each population
The user can disable this option by manually creating a file named sub.pop1/sub.pop2... , 
so that the script will use the user-specified population list

Options:
  -i FILE 		vcf file in bgzip and indexed with tabix
  -c INT		count
  -p FILE1,FILE2...	pop file
  -o DIR		output directory
  -t INT		threads count,use at most number of pops
  -h    		display this help and exit
  -D			dryrun
"

parse_args () {
    [[ "$#" == 0 ]] && die_usage "missing file argument"
    local OPTIND=1
    while getopts "i:p:c:t:o:hD" opt "$@"; do
        case "$opt" in
            i) VCF=$(realpath -s "$OPTARG")
               ;;
            p) POPS="$OPTARG";
	       POPS=(${POPS//,/ });
	       for x in $(seq 1 ${#POPS[@]}); do x=$[$x-1];POPS[$x]=$(realpath -s  ${POPS[$x]}); done
               ;;
	    c) SUBCNT="$OPTARG";
               ;;
            o) OUTDIR=$(realpath -s "$OPTARG")
               ;;
            t) NCPU="$OPTARG";
               ;;
            D) DRY_RUN=1;
               ;;
            h) show_help
               exit 0
	       ;;
            *) die_usage
        esac
    done
    shift $((OPTIND-1))
    return 0
}

die_usage () {
    [[ "$1" ]] && echo "`basename $0`: $1" >&2
    echo "$usage_msg" >&2
    exit 1
}
 
show_help () {
    printf "$help_msg"
}

parse_args "$@"
[ -z $NCPU ] && NCPU=2
[ -z $OUTDIR ] && OUTDIR=$(pwd)
mkdir -p $OUTDIR
DRY_RUN=${DRY_RUN:-0}

[ ! -s ${VCF}.tbi ]  && echo "vcf index not exist ,use tabix -fp vcf $VCF" && exit -1

for pop in ${POPS[@]}
do
    s=$(basename $pop)
    test -f ${OUTDIR}/sub.${s} || shuf -n $SUBCNT $pop > ${OUTDIR}/sub.${s}
done

for chr in $(tabix -l $VCF)
do
    for pop in ${POPS[@]}
    do
        s=$(basename $pop)
	cmd="bcftools view -S sub.${s} $VCF -r $chr | pigz > ${chr}.sub.${s}.vcf.gz"
        cmd="$cmd ; PopLDdecay -MaxDist 1000 -InVCF ${chr}.sub.${s}.vcf.gz -OutStat ld.${chr}.sub.${s}.vcf.gz "
	cmd="$cmd; rm -rf ${chr}.sub.${s}.vcf.gz"
        cmd="cd ${OUTDIR}; $cmd"
        echo $cmd
    done 
done > ${OUTDIR}/work.sh

#plot
cat <<EOF > ${OUTDIR}/plot.sh
cd $OUTDIR
>muti.list
for pop in ${POPS[@]}
do 
   s=$(basename $pop)
   cmd="ls \$(pwd)/ld.*.sub.${s}.vcf.stat.gz > ${s}.list \
 && Plot_OnePop.pl -inList ${s}.list -output ${s}.cat";echo \$cmd;\
 echo ${s}.cat.bin.gz ${s} >> muti.list
done|parallel --will-cite -j $NCPU; 
Plot_MultiPop.pl -inList  muti.list  -output  all.cat
EOF

cmd="parallel --will-cite -j $NCPU < ${OUTDIR}/work.sh && bash ${OUTDIR}/plot.sh"

[ $DRY_RUN -eq 0 ] && eval $cmd || ( echo $cmd > ${OUTDIR}/runme.sh ) 
