#!/usr/bin/env python

from fcntl import flock, LOCK_EX, LOCK_UN
from subprocess import PIPE, Popen,check_output
from threading import Timer
import shlex
import argparse
import sys
import time
import os
import tempfile #used for command line
import random

parser = argparse.ArgumentParser(description=
    'An industrious little ant for running serial tasks')
parser.add_argument('-f','--file', type=str,
    help='location of task file', required=True)
parser.add_argument('-l','--log', type=str,
    help='location of task log file', required=True)
parser.add_argument('-j','--jobcnt', type=int,nargs='?', const=1, default=1,
    help='number of jobs taken at a time', required=False)
parser.add_argument('-w','--walltime', type=int,nargs='?', const=300, default=300,
    help='max running time(seconds) of a task ', required=False)

args = parser.parse_args()
#numbers of tasks taken at a time

def runcmd(cmd):
    task=cmd
    #create temp file,text mode
    temp = tempfile.NamedTemporaryFile(prefix="tmpfile")
    #write the task to the temp file
    temp.write(cmd)
    temp.flush()
    cmd=os.environ.get('SWPATH')+"/time/bin/time -f %e-%S-%U-%P/%M-%t-%K/%x  bash "+temp.name
    proc = Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE)
    stdout_data, stderr_data = proc.communicate()

    stderr_data = stderr_data.rstrip(os.linesep).split(os.linesep)
    temp.close()
    rc = proc.returncode
    if(len(stderr_data)!=1):
        task=task+"\t@\t"+' <newline> '.join(stderr_data[0:-1])
    return '{exit_code}\t{stderr_data}\t{command}\n'.format(exit_code=rc,stderr_data=stderr_data[-1],command=task)

def putbackcmd(cmdlist):
    time.sleep(random.randint(0,10))
    with open(args.file, 'ab') as f:
        flock(f, LOCK_EX)
        while len(cmdlist):
           task=cmdlist.pop()
           if task == "":
              continue
           f.write(task+ os.linesep)
        f.flush()
        flock(f, LOCK_UN)
    
#reset log file
f_log = open(args.log, 'w')
f_log.seek(0)
f_log.truncate()
f_log.write(os.linesep.join(["##"+x for x in  (check_output("lscpu && uname -a", shell=True).strip()).decode().split(os.linesep)])+os.linesep*2)

while True:
    # try to open the task file
    with open(args.file, 'rb+') as f:
      flock(f, LOCK_EX)
      f.seek(0, os.SEEK_END)
      pos=f.tell()-1
      offs=100
      tasks=[]
      while True:
        pos=pos-offs if pos>=offs else 0
        f.seek(pos, os.SEEK_SET)
        lines=f.readlines()

        if(len(lines)==0):
           tasks=[]
           break
        elif(len(lines)> args.jobcnt  or pos == 0):
           tasks=lines[-args.jobcnt:]
           sum_len=0
           for x in tasks:
               sum_len+=len(x)           
           pos =f.tell()-sum_len
           break;
        offs*=2
        
      f.seek(pos, os.SEEK_SET)
      f.truncate()
      f.flush()
      flock(f, LOCK_UN)

    #traverse the task list and remove the task after finished
    while len(tasks):
      task=tasks.pop().rstrip(os.linesep)
      if(task ==""):
          continue    
      timer = Timer(args.walltime,putbackcmd,[tasks])
      try:
         timer.start()
         f_log.write(runcmd(task))
      finally:
         timer.cancel()
            
    size = os.path.getsize(args.file)
    if(size==0):
      break;

f_log.close()
