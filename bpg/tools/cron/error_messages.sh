#!/bin/bash
#
# Check output of df command for any volume above 89% capacity

header=$(df -P | head -1)
volumes=$(df -P | grep -v "$header" | awk '$5 > 90')

if [ -n "$volumes" ]; then
    mutt -s "Ohana: Disk Space Over Threshold" bpg.shailen@gmail.com << EOF
    At least one of the filesystems is over its threshold. Please review the
    following disk free command output:


    $(df)
EOF
fi
