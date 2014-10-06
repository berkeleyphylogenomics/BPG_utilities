#!/bin/sh

# Clean up files/directories in /tmp more than 3 days old.
# HTML directory.

find /tmp/* \( -name \* \) -mtime +3 -print0 | xargs -0 rm -rf
