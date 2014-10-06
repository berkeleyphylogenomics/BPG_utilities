#!/bin/sh

# HTML directory.

html_dir=/var/www/html

# Output/log files in this directory.

log_dir=/usr/local/blastdb/tmp_downloads

# Delete old results directories -- more than 45 days old.

find $html_dir/phylobuild/runs/* \( -name \* \) -mtime +45 -print0 | xargs -0 rm -r -f >& $log_dir/phylobuilder_clean.out

# Also remove results directories in /home/bpg/tmp/phylobuilder more than
# 30 days old.

find /home/bpg/tmp/phylobuild/* -type d -maxdepth 1 -mtime +30 -print0 | xargs -0 rm -r -f >& $log_dir/phylobuilder_clean.out

