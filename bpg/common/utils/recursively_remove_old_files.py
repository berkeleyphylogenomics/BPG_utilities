#!/usr/bin/python

import os, sys, stat, time, datetime

if len(sys.argv) < 3:
  print "Usage: %s <dir> <min_num_days_old>" % sys.argv[0]
  sys.exit(0)

dir_to_reap = sys.argv[1]

if not os.path.exists(dir_to_reap):
  print "Couldn't find directory %s" % dir_to_reap

try:
  min_num_days_old = int(sys.argv[2])
except ValueError:
  print "<min_num_days_old> must be an integer"
  sys.exit(0)

os.chdir(dir_to_reap)
all_files = os.listdir(dir_to_reap)
old_files = [file for file in all_files if
              datetime.timedelta(seconds=(time.time() -
              os.stat(file)[stat.ST_MTIME])).days >= min_num_days_old]

for file in old_files:
  os.system("rm -r %s" % file)
