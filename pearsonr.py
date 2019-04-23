#!/usr/bin/env python

import sys,getopt
opts,args=getopt.getopt(sys.argv[1:],"hi:")
input_file=""
output_file=""
for op,value in opts:
  if op=="-i":
    input_file=value
  elif op=="-h":
    print(sys.argv[0]+" -i infile ")
    sys.exit()

from  scipy.stats  import pearsonr

with open(input_file,"r") as f:
  for line in f.readlines():
    words=line.strip().split("\t",1)
    pearsonr0=list(map(int,words[0].split(",")))
    pearsonr1=list(map(int,words[1].split(",")))
    value=pearsonr(pearsonr0,pearsonr1)
    print(line.strip()+"\t"+str(value[0])+"\t"+str(value[1]))
