#!/bin/bash

readonly help_msg="\
Usage: $(basename $0) [OPTION]
`basename $0` is used to run filterx,wildcard and unsorted(to be completed) files are allowed
 eg filterx.sh -1 'cnt>=1' -k sn */GRch38.bga.gz.stats:1:4
Options:
  -u            unsorted file(to be completed)
  -D            dry-run
  -h            display this help and exit
"
parse_args () {
    [[ "$#" == 0 ]] && show_msg "missing file argument" && exit -1;
    local OPTIND=1
    while getopts "uhDk:b:e:1:2:3:4:5:6:7:8:9:" opt "$@"; do
        case "$opt" in
            h) show_msg
               filterx
               exit 0
               ;;
            D) dry_run=1
               ;;
            u) unsorted=1
               ;;
            k) cmd="$cmd -$opt $OPTARG"
               keys=$OPTARG
               ;;
            [be]) cmd="$cmd -$opt $OPTARG"
               ;;
            [1-9]) cmd="$cmd -$opt '$OPTARG'"
               ;;
        esac
    done
    shift $((OPTIND-1))
    count=1
    while [ -n "$1" ] ;do 
      s=`awk -v para=$1 -v unsorted=$unsorted -v keys=$keys 'BEGIN{\
         if(index(para,":")==0) {s=para":"1;for(i=2;i<=length(keys);++i) s=s":"i; print s;exit 0;}\
         f=substr(para,1,index(para,":")-1); para=substr(para,index(para,":")+1);\
         cmd="ls "f; s="";while((cmd|getline)>0) s=s" "$NF":"para;close(cmd); print s; }'`
      cmd="$cmd $s"
      shift 
    done

    cmd="$cmd $@"
}

show_msg () {
    [[ "$1" ]] && echo "`basename $0`: $1" >&2
    echo "$help_msg" >&2
}

cmd="filterx"
parse_args "$@"

##!/bin/bash

readonly help_msg="\
Usage: $(basename $0) [OPTION]
`basename $0` is used to run filterx,wildcard and unsorted(to be completed) files are allowed
 eg filterx.sh -1 'cnt>=1' -k sn */GRch38.bga.gz.stats:1:4
Options:
  -u            unsorted file(to be completed)
  -D            dry-run
  -h            display this help and exit
"
parse_args () {
    [[ "$#" == 0 ]] && show_msg "missing file argument" && exit -1;
    local OPTIND=1
    while getopts "uhDk:b:e:1:2:3:4:5:6:7:8:9:" opt "$@"; do
        case "$opt" in
            h) show_msg
               filterx
               exit 0
               ;;
            D) dry_run=1
               ;;
            u) unsorted=1
               ;;
            k) cmd="$cmd -$opt $OPTARG"
               keys=$OPTARG
               ;;
            [be]) cmd="$cmd -$opt $OPTARG"
               ;;
            [1-9]) cmd="$cmd -$opt '$OPTARG'"
               ;;
        esac
    done
    shift $((OPTIND-1))
    count=1
    while [ -n "$1" ] ;do 
      s=`awk -v para=$1 -v unsorted=$unsorted -v keys=$keys 'BEGIN{\
         if(index(para,":")==0) {s=para":"1;for(i=2;i<=length(keys);++i) s=s":"i; print s;exit 0;}\
         f=substr(para,1,index(para,":")-1); para=substr(para,index(para,":")+1);\
         cmd="ls "f; s="";while((cmd|getline)>0) s=s" "$NF":"para;close(cmd); print s; }'`
      cmd="$cmd $s"
      shift 
    done

    cmd="$cmd $@"
}

show_msg () {
    [[ "$1" ]] && echo "`basename $0`: $1" >&2
    echo "$help_msg" >&2
}

cmd="filterx"
parse_args "$@"

#commmand substitution is ok
eval $cmd
