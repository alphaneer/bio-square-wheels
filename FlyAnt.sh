#!/bin/bash
readonly help_msg="\
Usage: FlyAnts [-j ncpu*parallel_cnt] [-p 2] [-m mem] [-q small] taskfile task_log_file_prefix task_err_file_prefix
"

readonly usage_msg="Try 'FlyAnts -h' for more information."

not_being_sourced () { [[ "$0" == "$BASH_SOURCE" ]]; }

not_being_sourced && set -eu
#set -o pipefail
                      

main () {
        parse_args "$@" 
        task_cnt=$( wc -l "$taskfile"|cut -f1 -d" " )
        echo "running "$task_cnt" jobs"
        touch ${task_log_file_prefix}.1 ######FOR 'tail -n1 *.log*|awk ' AT LEAST 2 ITEMS EXIST
        touch ${task_log_file_prefix}.2
        cp $taskfile tmp.${taskfile}
        for((i=1;i<=parallel_cnt;++i)) 
        do
          cmd="qqsub '{-l nodes=1:ppn=$ncpu,mem=${memory}gb -q $queue_name}' \
'tasksant.sh $taskfile ${task_log_file_prefix}.$i ${task_err_file_prefix}.$i $runonce' ";
          eval $cmd
          sleep 1
          echo $cmd
        done

        test1=$( wc -l "$taskfile"|cut -f1 -d" " )      
        test2=($(tail -n1  ${task_log_file_prefix}.* |awk -F"\n" -v RS="" \
'{if(NF==1) next; else if(index($2,"EXE_FAIL_TOTAL_RATIO")==0) {run_cmd++;} else {split($2,tA,"\t");s2+=tA[2];s3+=tA[3];} } END{ print run_cmd+0,s2+0,s3+0;  }'));

        while(( ($test1>0)  || (${test2[0]} >0) || (${test2[2]} < $task_cnt)    ))
        do
          sleep 10;
          #echo "1:"$test1
          #echo "2:"${test2[0]}"-"${test2[1]}"-"${test2[2]}
          if (( $test1 < 1 )); then
            sleep 10
            test0="grep -vFf <(cat <(grep BEGIN ${task_log_file_prefix}.* | awk '{gsub(/.*command:/,\"\",\$0);print}') $taskfile) tmp.${taskfile} >> ${taskfile}  "
            echo $test0
            eval $test0 && echo ok || echo err 
            sleep 3
          fi
          echo "-"
          test1=$(wc -l $taskfile|cut -f1 -d' ');
          test2=($(tail -n1  ${task_log_file_prefix}.* |awk -F"\n" -v RS="" '{if(NF==1) next; else if(index($2,"EXE_FAIL_TOTAL_RATIO")==0) {run_cmd++;} else {split($2,tA,"\t");s2+=tA[2];s3+=tA[3];} } END{ print run_cmd+0,s2+0,s3+0;  }'));    
        done
        sleep 10
        
        truncate -s 0 $task_log_file_prefix
        truncate -s 0 $task_err_file_prefix
        for f in ${task_log_file_prefix}.*; do (cat "${f}"; echo) >> $task_log_file_prefix; done
        for f in ${task_err_file_prefix}.*; do (cat "${f}"; echo) >> $task_err_file_prefix; done
        test -s $taskfile || rm -rf $taskfile
        sed -i '/^$/d' $task_err_file_prefix
        test -s $task_err_file_prefix || rm -rf $task_err_file_prefix   
        if(( ${test2[1]} == 0 ));  then
                rm -rf ${task_log_file_prefix}.* ${task_err_file_prefix}.*
                rm -rf tmp.${taskfile}
                exit 0;
        else 
                echo 'error'
                exit -1;
        fi
}


parse_args () {
    dry_run_opt=0
    taskfile=""
    tasklogfile=""
    taskerrfile=""
    local OPTIND=1
    while getopts "Nhj:m:q:p:" opt "$@"; do
        case "$opt" in
            m) memory=$OPTARG
               ;;
            q) queue_name=$OPTARG
               ;;
            p) runonce=$OPTARG
               ;;
            N) dry_run_opt=1
               ;;
            h) show_help
               exit 0
               ;;
            j) str=$OPTARG
               str=${str//[*,:+]/ };
               arr=($str);
               ncpu=${arr[0]}
              parallel_cnt=${arr[1]};
                ;;
            *) die_usage
        esac
    done
    shift $((OPTIND-1))
    [[ "$#" == 0 ]] && die_usage "missing file argument"
    [[ "$#" -gt  4 ]] && die_usage "too many arguments"
    taskfile=$1;task_log_file_prefix=$2;task_err_file_prefix=$3;
}

show_help () {
    printf "$help_msg"
}

die_usage () {
    [[ "$1" ]] && echo "qs2 : $1" >&2
    echo "$usage_msg" >&2
    exit 1
}
 
not_being_sourced && main "$@"
