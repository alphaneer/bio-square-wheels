#!/bin/bash

task_file=$1
task_log_file=$2
task_err_file=$3
parallel=$4
[[ -z $parallel ]] && parallel=1
truncate -s 0 $task_log_file
truncate -s 0 $task_err_file


#read
command_execute () {
        local NP=$(test -n "$PBS_NODEFILE" &&  wc -l $PBS_NODEFILE|awk '{print $1}');
        local CPU=$(grep 'process' /proc/cpuinfo | sort -u | wc -l);CPU=$[$CPU/2];
        NP=${NP:-$CPU};
        local tmpfile=$(mktemp --tmpdir "peeif.XXXXXX");
        trap 'rm "$tmpfile"' EXIT;

        cat|grep -v '^[[:space:]]*$'|sed "s/MAGIC_CPU/${NP}/g"|tee "$tmpfile" |sed "s/^/BEGIN:command:/g" >> $task_log_file;
        #cmd="pstats.sh -o $task_log_file --append bash  $tmpfile";  
        /public/agis/ruanjue_group/lialun/sw/time/bin/time -v  -o $task_log_file --append bash  $tmpfile;
        #eval $cmd;
        local status="$?"  
        #chmod u+x $tmpfile
        #/public/agis/ruanjue_group/ruanjue/bin/runit.pl   $tmpfile >> $task_log_file 2>> $task_log_file        
        #local status=$(tac $task_log_file|grep 'retval'|head -n1|awk '{print $NF}')
        if [[ $status != 0 ]];  then 
                cat $tmpfile >> $task_err_file; 
        fi;
        trap  - EXIT;
        rm "$tmpfile"; 
        echo -e 'FINISHED\n\n\n' >> $task_log_file 
}

#print system info
echo '------------------------------------------------' >> $task_log_file;
egrep '^model name' /proc/cpuinfo |awk -F"\t" '{map[$2]++} END{for(x in map) print map[x],x}' >> $task_log_file;
egrep  'Mem|Cache|Swap' /proc/meminfo >> $task_log_file;
((hostname;uname -a)|tr '\n' ":";echo -e "\n";lsb_release -a;echo -e "\n";) >> $task_log_file;

(echo $PATH;echo -e "\n") >> $task_log_file;

(date;echo -e "\n")  >> $task_log_file;


while [[ -s $task_file ]]
do
 poptail -n $parallel $task_file |for((i=1;i<=${parallel};++i)) ; do command_execute; done
done
date >> $task_log_file;

echo -e "EXE_FAIL_TOTAL_RATIO:\t" $(wc -l $task_err_file|cut -f1 -d" ") "\t" $(grep '^BEGIN:command:' $task_log_file|wc -l|cut -f1 -d" ") >> $task_log_file;
