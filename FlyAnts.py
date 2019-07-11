#!/usr/bin/env python

import argparse
import sys,time
import glob

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,description=
   "An industrious little ant for running serial tasks \n\
usage: nohup FlyAnts.py -f work.sh -L work.sh.log.* [-P 100] &   \n\
seq 1 10 |xargs -I {} qsubx +queue7:12:8gb+ taskant.py -f work.sh -l work.sh.log.{} -j 40 &")

parser.add_argument('-f','--file', type=str,
    help='location of task file', required=True)

parser.add_argument('-L','--alllogs', type=str,
    help='match pattern of all the task log files', required=True)

parser.add_argument('-P','--polltime', type=int,nargs='?',const=100, default=100,
    help='Polling interval,default 100', required=False)

args = parser.parse_args()


class MultipleFileManager(object):
    def __init__(self, files):
        self.files = files

    def __enter__(self):
        self.open_files = map(open, self.files)
        return self.open_files

    def __exit__(self,exception_type, exception_value, traceback):
        for f in self.open_files:
            f.close()


if (args.alllogs != None):
    with open(args.file, 'r') as f:
       cmdsets = set( x.strip() for x in f if x.strip()!="")

    alllogs = [f for f in  glob.glob(args.alllogs)]
    with MultipleFileManager(alllogs) as files:
      while len(cmdsets):
        bwait=True
        line=""
        for f in files:
          #get first non-empty body,header exclude
          for line in iter(f.readline,b''):
             line=line.strip()
             if (line !="" and not line.startswith("#")):
               break
          #skip blank rows
          if not line:
             continue
          line=line.split("\t",2)[-1].strip()
          pos=line.find("\t@\t",0)
          if pos !=-1:
            line=line[0:pos]
          cmdsets.remove(line)
          bwait=False

        if bwait==True:  
          time.sleep(args.polltime)
    
    sys.exit(0)
