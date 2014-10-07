import tempfile
import os
import stat

ROOT_DIR = '/clusterfs/ohana/software/webserver/temp/'

def remove_duplicates(seq): 
  # order preserving
  noDupes = []
  [noDupes.append(i) for i in seq if not noDupes.count(i)]
  return noDupes


def get_temp_dir(prefix='xxx'):

    """Create temporary directory for webserver use using prefix.

    All webserver temporary files are stored in the ROOT_DIR. Further temporary
organization is created by creating directories with the prefix specified
(i.e., 'satchmo', 'phylofacts', 'xxx', etc.). Examples of temporary directories include:

* /clusterfs/ohana/software/webserver/temp/xxx2CcLEH
* /clusterfs/ohana/software/webserver/temp/xxx5CDL3H
* /clusterfs/ohana/software/webserver/temp/satchmo2c99H

    If the directory structure up to the path (e.g., ROOT_DIR isn't complete)
doesn't exist, it is created. The directory is made with all read/write/execute
permissions (777).

    As a final result, the path to the directory is returned.
"""

    target_dir = os.path.join(ROOT_DIR, prefix)
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    tempdir = tempfile.mkdtemp(prefix=target_dir)
    os.chmod(tempdir, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR |\
                      stat.S_IRGRP | stat.S_IWGRP | stat.S_IXGRP |\
                      stat.S_IROTH | stat.S_IWOTH | stat.S_IXOTH)

    return tempdir
