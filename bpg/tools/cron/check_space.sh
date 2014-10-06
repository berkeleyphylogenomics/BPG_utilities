#!/bin/bash
#
# Check output of df command for any volume above 89% capacity

v_header=$(df -P | head -1)
i_header=$(df -iP | head -1)
volumes=$(df -P | grep -v "$v_header" | awk '$5 > 90')
i_nodes=$(df -iP | grep -v "$i_header" | awk '$5 > 94')
me=$(uname -n)

if [ -n "$volumes" ]; then
    mutt -s "Ohana ($me): Disk Space Over Threshold" bpg.shailen@gmail.com << EOF
    At least one of the filesystems on $me is over its threshold.
    
    $volumes

    Please review the following disk free command output:

    $(df)
EOF
fi

if [ -n "$i_nodes" ]; then
    mutt -s "Ohana ($me): inodes over threshold" bpg.shailen@gmail.com << EOF

    At least one of the filesystems on $me is over its inode threshold: 

    $i_nodes
    
    Please review the following disk free inode command output (df -i):

    $(df -i)
EOF
fi
