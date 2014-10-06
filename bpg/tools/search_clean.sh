#!/bin/sh

# HTML directory.

html_dir=/var/www/html

# Output/log files in this directory.

log_dir=/usr/local/blastdb/tmp_downloads

# Delete old results directories -- more than 60 days old.

(find $html_dir/gdata/* \( -name \* \) -mtime +60 -print0 | xargs -0 rm -r -f) >& $log_dir/search_clean.out

# Also remove tmp files like phylofacts.remote.16226.qjob0.sh -- more
# than 2 days old.

(find /home/bpg/tmp/* -type f -maxdepth 1 -mtime +2 -print0 | xargs -0 rm -f) >& $log_dir/search_clean.out

# Also remove results directories in /home/bpg/tmp/phylofacts more than
# 14 days old.

(find /home/bpg/tmp/phylofacts/* -type d -maxdepth 1 -mtime +14 -print0 | xargs -0 rm -r -f) >& $log_dir/search_clean.out

