#!/usr/bin/env python

import os, sys, string, re, cPickle
from Bio.SeqUtils import CheckSum
from Bio import SeqIO

uppercase_translation = string.maketrans(string.lowercase, string.uppercase)
trivial_translation = string.maketrans('', '')
dotdash = '.-'

id_re = re.compile('[^(:),]*')
branch_length_re = re.compile('[0-9]\.[0-9e\-]*')

class node:
  def __init__(self, seguid=None):
    self.children = set()
    self.leftId = None
    self.rightId = None
    if seguid:
      self.contained_seguids = set([seguid])
    else:
      self.contained_seguids = set()

  def addChild(self, child):
    self.children.add(child)
    self.contained_seguids = self.contained_seguids | child.contained_seguids

  def updateLeftId(self, leftId):
    self.leftId = leftId
    rightId = leftId + 1
    for child in self.children:
      rightId = child.updateLeftId(rightId)
    self.rightId = rightId
    return rightId + 1

  def readFromTreeString(self, tree_string, seguids, i):
    if i >= len(tree_string):
      return len(tree_string)
    if tree_string[i] == '(':
      while tree_string[i] != ')':
        child = node()
        i = child.readFromTreeString(tree_string, seguids, i+1)
        self.addChild(child)
      i += 1
    else:
      id = id_re.match(tree_string[i:]).group(0)
      i += len(id)
      self.seguid = seguids[id]
      self.contained_seguids = set([self.seguid])
    if i < len(tree_string) and tree_string[i] == ':':
      i += 1
      branch_length = branch_length_re.match(tree_string[i:]).group(0)
      i += len(branch_length)
    return i

  def updateAlignmentOffsetOfLeftId(self, alignment_offset_dict,
                                    alignment_offset_of_left_id):
    if len(self.contained_seguids) > 0:
      for seguid in self.contained_seguids:
        if seguid in alignment_offset_dict:
          if len(self.contained_seguids) in alignment_offset_dict[seguid]:
            alignment_offset_of_left_id[self.leftId] \
              = alignment_offset_dict[seguid][len(self.contained_seguids)]
          break
    for child in self.children:
      child.updateAlignmentOffsetOfLeftId(alignment_offset_dict,
                                          alignment_offset_of_left_id)

def parse_tree(work_path):
  seguids = {}
  f = open(os.path.join(work_path, "input_unaligned.fasta"))
  for record in SeqIO.parse(f, "fasta"):
    id = record.id.replace(':', '_')
    seguids[id] = CheckSum.seguid(record.seq)
  f.close()
  f = open(os.path.join(work_path, "satchmo_tree.newick"))
  tree_string = f.read()
  f.close()
  tree_string = tree_string.translate(trivial_translation, string.whitespace)
  root = node()
  root.readFromTreeString(tree_string, seguids, 0)
  root.updateLeftId(1)
  return root

def parse_smo(work_path):
  f = open(os.path.join(work_path, "satchmo.smo"))
  alignment_offset_dict = {}
  records = set()
  current_header = ""
  current_sequence = ""
  alignmentOffset = 0
  alignmentNumBytes = 0
  withinAlignment = False
  start_of_line = f.tell()
  line = f.readline()
  while line:
    if line.rstrip() == 'alignment':
      # skip past the line with the curly brace
      f.readline()
      # the next line is the start of the alignment
      alignmentOffset = f.tell()
      withinAlignment = True
    elif line.rstrip() == '//':
      if current_sequence != '':
        seguid = CheckSum.seguid(current_sequence)
        records.add( seguid )
      alignmentNumBytes = start_of_line - alignmentOffset
      for seguid in records:
        if seguid not in alignment_offset_dict:
          alignment_offset_dict[seguid] = {}
        alignment_offset_dict[seguid][len(records)] \
            = (alignmentOffset, alignmentNumBytes)
      records = set()
      current_header = ""
      current_sequence = ""
      withinAlignment = False
    elif withinAlignment:
      if len(line) > 0 and line[0] == '>':
        if current_sequence != '':
          seguid = CheckSum.seguid(current_sequence)
          records.add( seguid )
        current_header = line[1:].rstrip()
        current_sequence = ""
      else:
        current_sequence = current_sequence + \
            line.strip().translate(uppercase_translation, dotdash)
    start_of_line = f.tell()
    line = f.readline()
  f.close()
  return alignment_offset_dict

def find_alignment_offset_of_left_id(work_path):
  alignment_offset_dict = parse_smo(work_path)
  root = parse_tree(work_path)
  alignment_offset_of_left_id = {}
  root.updateAlignmentOffsetOfLeftId(alignment_offset_dict,
                                    alignment_offset_of_left_id)
  return alignment_offset_of_left_id


def main():
  if len(sys.argv) < 2:
    path = os.getcwd()
  else:
    path = sys.argv[1]

  alignment_offset_of_left_id = find_alignment_offset_of_left_id(path)
  f = open(os.path.join(path, 'alignment_offset_of_left_id.pkl'), "w")
  cPickle.dump(alignment_offset_of_left_id, f)
  f.close()

if __name__ == '__main__':
  main()
