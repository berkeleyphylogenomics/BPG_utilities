#!/usr/bin/env python

import os, sys, glob, string
from matchmaker.shmm_shmm_lib import *
from optparse import OptionParser

def main():
  usage = "%prog [options] <seed_id>"
  opt_parser = OptionParser(usage=usage)
  opt_parser.add_option("-p", "--percent_id", type="int", dest="percent_id", 
                        default=20,
                        help="Minimum %id at which subfamilies were cut")
  opt_parser.add_option("--use_ss", dest="use_ss", action="store_true",
    default=True,
    help="Whether to include secondary structure prediction in the HMM")
  opt_parser.add_option("--no_use_ss", dest="use_ss", action="store_false",
    default=True,
    help="Whether not to include secondary structure prediction in the HMM")
  (options, args) = opt_parser.parse_args()

  if len(args) < 1:
    opt_parser.error("The seed_id argument is required.")

  seed_id = args[0]
  percent_id = int(options.percent_id)
  use_ss = options.use_ss

  subfam_dir = kerf_subfam_dir(seed_id, percent_id)
  os.chdir(subfam_dir)

  kerf_subfamily_files = glob.glob("%s.kerf%d.sf*.fa" % (seed_id, percent_id))
  num_subfams = len(kerf_subfamily_files)

  hmm_dir = kerf_hhalign_hmm_dir(seed_id, percent_id, use_ss)
  make_dir_exist(hmm_dir)
  os.chdir(hmm_dir)
  hmm_file_spec = "%s.hhalign*kerf%d.sf*.hhm" % (seed_id, percent_id)
  num_hmms = len(glob.glob(hmm_file_spec))
  if num_hmms == num_subfams:
    print "HHalign hmms already created for kerf%d on %s" % (percent_id, seed_id)
    return
  if num_hmms > 0:
    cmd = "rm %s" % hmm_file_spec
    n = os.system(cmd)
  trivial_translation = string.maketrans('', '')
  dot = '.'
  if use_ss:
    seed_pred = array('c')
    seed_conf = array('c')
    f = open(os.path.join(single_seqs_dir(), seed_id, 'vanilla_psipred',
              '%s.horiz' % seed_id))
    for line in f.readlines():
      if line.startswith('Pred:'):
        seed_pred.fromstring(line[6:].strip())
      elif line.startswith('Conf:'):
        seed_conf.fromstring(line[6:].strip())
    f.close()
  for kerf_subfam_fa in kerf_subfamily_files:
    subfam_basename = os.path.splitext(kerf_subfam_fa)[0]
    subfam_basename_components = subfam_basename.split('.')
    if use_ss:
      kerf_subfam = "%s.hhalignss%s" % (subfam_basename_components[0],
                                      '.'.join(subfam_basename_components[1:]))
    else:
      kerf_subfam = "%s.hhalign%s" % (subfam_basename_components[0],
                                      '.'.join(subfam_basename_components[1:]))
    inf = open("../%s" % kerf_subfam_fa)
    outf = open("%s_unfiltered.a3m" % kerf_subfam, "w")
    for line in inf.readlines():
      outf.write(line.translate(trivial_translation, dot))
    inf.close()
    if use_ss:
      outf.write(">ss_pred\n")
      outf.write("%s\n" % seed_pred.tostring())
      outf.write(">ss_conf\n")
      outf.write("%s\n" % seed_conf.tostring())
    outf.close()
    cmd = "hhfilter -M a3m -id 90 -cov 20 -diff 0 " \
          + "-i %s_unfiltered.a3m " % kerf_subfam \
          + "-o %s.a3m " % kerf_subfam \
          + "> hhfilter_%s.out " % kerf_subfam
    n = os.system(cmd)
    if n != 0:
      print "hhfilter %s exited with nonzero status %d" % (kerf_subfam, n)
    cmd = "hhmake -i %s.a3m > hhmake_%s.out" % (kerf_subfam, kerf_subfam)
    n = os.system(cmd)
    if n != 0:
      print "hhmake %s exited with nonzero status %d" % (kerf_subfam, n)
    cmd = "hhsearch -cal -i %s.hhm " % kerf_subfam \
          + "-d /clusterfs/ohana/external/cal.hhm " \
          + "> calibrate_%s.out" % kerf_subfam
    n = os.system(cmd)
    if n != 0:
      print "calibrate %s exited with nonzero status %d" % (kerf_subfam, n)


if __name__ == '__main__':
  main()
