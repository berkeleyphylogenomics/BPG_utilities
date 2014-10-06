#!/usr/bin/python
import os, sys, glob, time
from optparse import OptionParser

def main():
  # parse command line options
  usage = "%prog [options] clustered_seeds_directory"
  opt_parser = OptionParser(usage=usage)
  opt_parser.add_option("-s", "--total_num_shards", dest="num_shards", 
                        default=1,
                        help="Number of shards into which seeds are divided")
  opt_parser.add_option("-n", "--current_shard_number", dest="shard_number",
                        default=0,
                        help="Which shard number the current job is processing")
  (options, args) = opt_parser.parse_args()
  if len(args) > 1:
    opt_parser.error('Incorrect number of arguments')
  if len(args) < 1:
    clustered_seeds_dir = os.getcwd()
  else:
    clustered_seeds_dir = args[0]
  os.chdir(clustered_seeds_dir)
  if not os.path.exists("clusters"):
    opt_parser.error("No file named 'clusters' in %s" % clustered_seeds_dir)
  try:
    num_shards = int(options.num_shards)
  except ValueError:
    opt_parser.error("--total_num_shards must be a number")
  try:
    shard_number = int(options.shard_number)
  except ValueError:
    opt_parser.error("--current_shard_number must be a number")
  f = open("clusters", "rU")
  lines = f.readlines()
  f.close()
  cluster_num = 0
  found_seed = False
  print "%d lines in clusters file" % len(lines)
  for line in lines:
    line = line.rstrip()
    if line[0:9] == 'Cluster #':
      cluster_num = int(line[9:-1])
      found_seed = False
    elif not found_seed and cluster_num % num_shards == shard_number:
      if line.find('|') >= 0:
        seed_id = line.split('|')[1]
      else:
        seed_id = line
      # Replace / in the header by %2f
      seed_id = seed_id.replace('/','%2f')
      # Replace ( in the header by %28
      seed_id = seed_id.replace('(','%28')
      # Replace ) in the header by %29
      seed_id = seed_id.replace(')','%29')
      # Replace - in the header by %2d
      seed_id = seed_id.replace('-','%2d')
      # Replace * in the header by %2a
      seed_id = seed_id.replace('*','%2a')
      os.chdir(seed_id)
      print "Cluster number: %d" % cluster_num
      print "Seed directory: %s" % seed_id
      cmd = "source /etc/profile; " \
          + "flowerpower.pl -i %s.fa -n 7 --tempcheck 1 >& flowerpower.out" \
            % seed_id
      print cmd
      t1 = time.clock()
      os.system(cmd)
      t2 = time.clock()
      print "FlowerPower finished on seed %s in %d seconds" % (seed_id, t2-t1)
      os.chdir('..')
      found_seed = True
  os.system('touch flowerpower_job_%d_of_%d_done' % (shard_number, num_shards))

if __name__ == '__main__':
  main()
