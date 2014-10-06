#!/bin/env python

"""Check Disk space - Ohana module

This module is used as a library (import), or a silent cron script
only. For cron script use assistance, execute this module with --help.

For interactive use, see the manage_ohana.py master script.
"""

import csv
from datetime import datetime
from optparse import OptionParser
import os
from socket import gethostname
import sys
import statvfs
import subprocess

try:
    from bpg.manage_ohana.shared import get_self_executable, logging, press_enter
except ImportError:
    sys.stderr.write("Activate environment first.")
    sys.exit(1)

CAPACITY_DISK_FILE =\
                 '/home/glenjarvis/capacity_planning/capacity_disk.csv'
NOT_APPLICABLE = "N/A"
THRESHOLD = 80

mount_points = [
    '/',
    '/boot',
    '/dev/shm',
    '/mnt',
    '/var/lib/pgsql',
    '/clusterfs/ohana/home',
    '/clusterfs/ohana/software',
    '/clusterfs/ohana/external',
    '/clusterfs/ohana/bpg',
    '/clusterfs/ohana/matchmaker',
    '/clusterfs/ohana/sandbox',
]

nodes = [
    'ohana',
    'makana',
    'db',
    'n0000', 'n0001', 'n0002', 'n0003', 'n0004', 'n0005', 'n0006', 'n0007',
    'n0008', 'n0009', 'n0010', 'n0011', 'n0012',
    'n0018', 'n0019', 'n0020', 'n0021', 'n0022', 'n0023',
    'n0024', 'n0025', 'n0026', 'n0027', 'n0028', 'n0029', 'n0030', 'n0031',
    'n0033', 'n0034', 'n0035', 'n0036', 'n0037', 'n0038', 'n0039',
    'n0040', 'n0041', 'n0042', 'n0043', 'n0044', 'n0045', 'n0046', 'n0047',
    'n0048', 'n0049', 'n0050', 'n0051', 'n0052'
]

binary_prefixes = [
    'B', # bytes
    'KiB', # 1024 bytes (aka 2**10 bytes)
    'MiB', # 1024**2 bytes (aka 2**20 bytes)
    'GiB', # 1024**3 bytes (aka 2**30 bytes)
    'TiB', # 1024**4 bytes (aka 2**40 bytes)
    'PiB', # 1024**5 bytes (aka 2**50 bytes)
    'EiB', # 1024**6 bytes (aka 2**60 bytes)
    'ZiB', # 1024**7 bytes (aka 2**70 bytes)
    'YiB', # 1024**8 bytes (aka 2**80 bytes)
]

def human_readable_units(value):
    for x in binary_prefixes:
        if value < 1024.0:
            return "%3.1f%s" % (value, x)
        value /= 1024.0


def disk_free(print_to_screen=True):
    capacityWriter = csv.writer(open(CAPACITY_DISK_FILE, 'a'), dialect='excel')
    when = datetime.today()

    for point in mount_points:
        try:
            stat = os.statvfs(point)
        except OSError:
            stat = [0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0]
        block_size = stat[statvfs.F_BSIZE]
        blocks_total = stat[statvfs.F_BLOCKS]

        # blocks_free and blocks avail are not the same. Some free blocks are
        # not available because they are reserved for root permissions
        # We will be conservative and calculate blocks_used as total - avail
        blocks_free = stat[statvfs.F_BFREE]
        blocks_avail = stat[statvfs.F_BAVAIL]
        blocks_used = blocks_total - blocks_avail

        size_bytes = blocks_total * block_size # Total size in bytes
        used_bytes = (blocks_total - blocks_avail) * block_size
        avail_bytes = blocks_avail * block_size

        # Percentage is calculated on blocks, not bytes since allocation for use
        # is done atomically by the block
        if blocks_total > 0:
            use_percent = (blocks_used)*100/blocks_total
        else:
            use_percent = NOT_APPLICABLE

        capacityWriter.writerow((point, when.strftime('%Y-%m-%d'),
            gethostname(), size_bytes, used_bytes, avail_bytes, use_percent))
        if print_to_screen:
            if not use_percent == NOT_APPLICABLE:
                if int(use_percent) >= THRESHOLD:
                    print "%s at %s%% (WARNING)" % (point, use_percent)
                    logging.warning("%s:%s at %s capacity" % (gethostname(),
                        point, use_percent))

def check_disk_space(verbosity=False):

    for node in nodes:
        if verbosity:
            print "Checking %s" % node

        self_executable = get_self_executable(__file__)

        retcode = subprocess.call(['/usr/bin/ssh', node, self_executable,
            '-ms'])
        if retcode != 0:
            print "ERROR trying to check on %s..." % node
            logging.error("Error trying to check on the %s..." % node)
        else:
            logging.info("Node %s succesfully checked." % node) 

    if verbosity:
        press_enter()


def main():
    """Options are passed and start_qa is called"""

    parser = OptionParser(version='%prog 0.1')
    parser.add_option('-s', '--scriptable', dest='scriptable',
        action="store_true",
        help="For cron job scripting - no user interaction is needed",
        default=False)
    parser.add_option('-m', '--machine', dest='machine',
        action="store_true",
        help="Run a check on current logged in machine only",
        default=False)

    (options, args) = parser.parse_args()

    if not options.scriptable:
        print __doc__
    else:
        if options.machine:
            disk_free()
        else:
            check_disk_space()

if __name__ == '__main__':
    main()
