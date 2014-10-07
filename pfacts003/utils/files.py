# Define absolute paths and filename functions
from datetime import datetime
import os
import re
import random

from django.conf import settings
from django.http import HttpResponse, Http404

from id_patterns import bpgid_pat, scopid_pat


pre = settings.QUEUED_INCOMING_DIR

BLAST_DB_DIR = "/clusterfs/ohana/external/"
PHOG_BLAST_DIR = os.path.join(pre, "phog_blast/")
BOOK_ROOT_DIR = "/clusterfs/ohana/bpg/pfacts/"
BOOK_ROOT_DIR_WEB = os.path.join(pre, "pfacts") # so wrong
VIEW_TREE_WORK_DIR = "./adata/aw2WJi.d"


def location_to_dir(location):
  if location == 'phogblast':
    return PHOG_BLAST_DIR
  else:
    return ''   
  
def open_location_file(location, file_name):
  file_name = location_to_dir(location) + file_name
  try:
    return open(file_name, "r")
  except IOError, e:
    return None 
  
def family_dir(family_name, webdir = False):
  if re.match(bpgid_pat, family_name) or re.match(scopid_pat, family_name):
    if webdir:
      dir = BOOK_ROOT_DIR_WEB
    else:
      dir = BOOK_ROOT_DIR + family_name[0:4] + '/'
    return dir + family_name[0:7] + '/' + family_name + '/'
  else:
    return None
  
def view_tree_work_dir():
  return VIEW_TREE_WORK_DIR  
  
def family_nj_file_name(family_name, webdir = False):
  dir = family_dir(family_name, webdir)
  if dir:
    return dir + 'user/' + family_name + '.nj'
  else:
    return None
  
def fasta_input_file_name(location, base_name):
  return location_to_dir(location) + base_name + '_input.fa'
  
def blast_results_file_name(location, base_name):
  return location_to_dir(location) + base_name + '_results.xml'
  
def pfam_color_file_name(family_name):
  if family_dir(family_name):
    return family_dir(family_name) + 'user/' + family_name + '_pfam_colorfile'
  else:
    return None
  
def open_pfam_color_file(family_name):
  file_name = pfam_color_file_name(family_name)
  if file_name:
    try:
      return open(file_name, "r")
    except IOError, e:
      return None 

def datetime_str(append_random = False):
  now = datetime.now()
  result_str = '%d%02d%02d_%02d%02d%02d' \
    % (now.year, now.month, now.day, now.hour, now.minute, now.second)
  if append_random:
    result_str += '_%03d' % random.randint(0, 999)
  return result_str

def create_file(directory, name): 
  handle = None
  count = 0
  not_created = True
  while not_created and count <= 25:
    base_name = datetime_str(True)
    if callable(name):
      full_path_name = name(directory, base_name)
    else:
      full_path_name = directory + base_name + name
    if os.path.exists(full_path_name):
      count += 1   
    else:
      handle = open(full_path_name, "w")
      not_created = False
  if not not_created:
    handle.close()
  return (full_path_name, base_name) 


def serve_file(filename="bpg.txt", content=None, is_download=False,
               mimetype='text/plain'):
    """Return Django method to serve content of file

    The result of this function returns a Django HttpResponse object
    that can be returned directly from a view. The mimetype of that
    object can be specified (or be text by default).

    If content is not None, then that content will be served into the
    HttpRequest object. Additionally, if filename is specified and
    this should be a file download, the filename will be used to serve
    that content.

    If the content is None, the contents of filename is read and served
    via the HttpResponse object. If there is no file at filename, an
    Http404 exception is raised.
    """

    response = HttpResponse(mimetype=mimetype)

    if content is None:
        # Filter out any errors in filename path
        filename = os.path.normpath(filename)
        if os.path.exists(filename):
            f = open(filename, 'r')
            content = f.read()
            f.close()
        else:
            raise Http404

    if is_download:
        response['Content-Disposition'] = \
            'attachment; filename=%s' % filename

    print >>response, content,
    return response
