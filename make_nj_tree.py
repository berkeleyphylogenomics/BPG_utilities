#!/usr/bin/python

# Make neighbor-joining tree from this alignment.

import os, sys
from matchmaker.shmm_shmm_lib import *

if len(sys.argv) < 2:
  print "Usage: %s <seed_id>" % sys.argv[0]
  sys.exit(0)

seed_id = sys.argv[1]
root = "%s/%s" % (cur_dir(), seed_id)
alignment_path = os.path.join(root, "%s.a2m" % seed_id)
work_dir = os.path.join(root, "NJ_tree_workspace")
make_dir_exist(work_dir)
os.chdir(work_dir)
os.system("ln -s ../%s.a2m %s.a2m" % (seed_id, seed_id))

# Run alignment through pip to simplify identifiers.

cmd = "source /etc/profile; "
cmd += "pip --pre --idfile pip_ids --msa %s.a2m --biglabel" % seed_id

os.system( cmd )

# Prettyalign with option to remove inserts (otherwise clustalw converts
# lowercase to uppercase, and then protdist doesn't see an alignment).
# -f: fasta format.  -m0: no inserts.

cmd = "source /etc/profile; "
cmd += "prettyalign %s.a2m.pip -f -m0 > alignment.afa 2> tree_prettyalign.out" \
        % seed_id

os.system( cmd )

# Use clustalw to convert to phylip format.

cmd = "source /etc/profile; "
cmd += "clustalw -infile=alignment.afa -outfile=alignment.phy " \
               + "-output=phylip -convert >& tree_clustalw.out; "
cmd += "cp alignment.phy infile "

os.system( cmd )

cmd = "source /etc/profile; "
cmd += 'echo "y\n" | protdist >& protdist.out '

os.system( cmd )

# ........................................................................
# Create neighbor-joining tree.

print "Making NJ tree..."

cmd = "source /etc/profile; "
cmd += "mv outfile alignment.dist; "
cmd += "cp alignment.dist infile; "
cmd += 'echo "y\n" | neighbor >& neighbor.out '

os.system( cmd )

# Save as tree.nj.  Re-root the tree.

cmd = "source /etc/profile; "
cmd += "mv outtree tree.nj; "
cmd += "mpreroot tree.nj rooted.nj"

os.system( cmd )

# Use pip to add labels back.  Save result as <seed_id>.nj.

cmd = "source /etc/profile; "
cmd += "pip --post --idfile pip_ids --tree rooted.nj; "
cmd += "mv rooted.nj.pip %s " % os.path.join(root, "%s.nj" % seed_id)

os.system( cmd )

