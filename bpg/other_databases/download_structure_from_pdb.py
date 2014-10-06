#!/usr/bin/env python

import os
import stat
import ftplib

def download_structure_from_pdb(pdb_id, pdb_dir, 
                                ftp=None, force_download=False):
  pdb_file = os.path.join(pdb_dir, 'pdb%s.ent' % pdb_id)
  if force_download or not os.path.exists(pdb_file) \
      or os.stat(pdb_file)[stat.ST_SIZE] == 0:
    if not ftp:
      ftp = ftplib.FTP('ftp.wwpdb.org')
      ftp.login()
    cur_dir = os.getcwd()
    os.chdir(pdb_dir)
    f = open("%s.gz" % pdb_file, "wb")
    try:
      ftp.retrbinary('RETR pub/pdb/data/structures/divided/pdb/%s/pdb%s.ent.gz'
                      % (pdb_id[1:3], pdb_id), f.write)
    except ftplib.error_perm:
      print "Failed to retrieve pdb%s.ent.gz" % pdb_id
      os.chdir(cur_dir)
      raise AssertionError
    f.close()
    os.system("gunzip %s.gz" % pdb_file)
    os.chdir(cur_dir)
    return ftp
  return None
