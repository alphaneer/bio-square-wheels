#!/usr/bin/env python 

from fcntl import flock, LOCK_EX, LOCK_UN, LOCK_NB
from subprocess import PIPE, Popen,check_output
from threading import Timer
import shlex
import argparse
import sys
import time
import os
import tempfile #used for command line
import random


def putbackcmd(cmdlist):
	with open(args.file, 'ab') as f:
		flock(f, LOCK_EX)
		while len(cmdlist):
			task = cmdlist.pop()
			if task == "":
				continue
			f.write((os.linesep+task).encode())
		flock(f, LOCK_UN)


def runcmd(cmd):
	task=cmd
	#create temp file,text mode
	temp = tempfile.NamedTemporaryFile(prefix="tmpfile")
	#write the task to the temp file
	temp.write(cmd.encode())
	temp.flush()
	cmd="\\time -f %e-%S-%U-%P/%M-%t-%K/%x	bash "+temp.name
	proc = Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE)
	stdout_data, stderr_data = proc.communicate()

	stderr_data = stderr_data.decode().rstrip(os.linesep).split(os.linesep)
	temp.close()
	rc = proc.returncode

	if len(stderr_data) != 1 and  rc != 0:
		task = task+"\t@\t"+' <newline> '.join(stderr_data[0:-1])
	return (rc,'{exit_code}\t{stderr_data}\t{command}'.format(exit_code=rc,
															stderr_data=stderr_data[-1],
															command=task))


def get_last_cmd(file, n):
	res = list()
	with open(args.file, 'r+') as f:
		flock(f, LOCK_EX)
		f.seek(0, os.SEEK_END)
		pos = f.tell() - 1
		for _ in range(0, n):
			while pos > 0 and f.read(len(os.linesep)) != os.linesep:
				pos -= 1
				f.seek(pos, os.SEEK_SET)

		res = f.readlines()
		if pos >= 0:
			f.seek(pos, os.SEEK_SET)
			f.truncate()
		flock(f, LOCK_UN)
	return res


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description=
									 'An industrious little ant for running serial tasks')

	parser.add_argument('-f', '--file', type=str,
						help='location of task file', required=True)

	parser.add_argument('-l', '--log', type=str,
						help='location of task log file', required=False)

	parser.add_argument('-j', '--jobcnt', type=int, nargs='?', const=1, default=1,
						help='number of jobs taken at a time', required=False)
	parser.add_argument('-w', '--walltime', type=int, nargs='?', const=300, default=300,
						help='max running time(seconds) of a task ', required=False)

	args = parser.parse_args()
	retries_dict=dict()

	# reset log file
	f_log = open(args.log, 'w')
	f_log.seek(0)
	f_log.truncate()
	f_log.write(os.linesep.join(["##"+x for x in  (check_output("lscpu && uname -a", shell=True).decode().strip()).split(os.linesep)])+os.linesep*2)
    f_log.flush()
    
	while True:
		# try to open the task file
		tasks = get_last_cmd(args.file, args.jobcnt)
		if len(tasks) == 0:
			f_log.flush()
			break

		# traverse the task list and remove the task after finished
		while len(tasks):
			task = tasks.pop().rstrip(os.linesep)
			if task == "":
				continue
			timer = Timer(args.walltime, putbackcmd, [tasks] )
			try:
				timer.start()
				(rc,ret_str)=runcmd(task)
				if rc != 0:
					if task in retries_dict:
						retries_dict[task]+=1
					else:
						retries_dict[task]=1
					if retries_dict[task] < 3:
						putbackcmd([task])
				f_log.write(ret_str+os.linesep)
			finally:
				timer.cancel()

	f_log.close()
