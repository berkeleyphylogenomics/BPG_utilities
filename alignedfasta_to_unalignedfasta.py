#!/usr/bin/python

import os, sys, string

if len(sys.argv) >= 2:
  input_handle = open(sys.argv[1], "rU")
else:
  input_handle = sys.stdin

if len(sys.argv) >= 3:
  output_handle = open(sys.argv[2], "w")
else:
  output_handle = sys.stdout

uppercase_translation = string.maketrans(string.lowercase, string.uppercase)
dotdash='.-'
writingSequence = False
for line in input_handle.readlines():
  # ignore blank lines
  if len(line.strip()) > 0:
    if line[0] == '>':
      # This is a header line
      # Terminate the preceding sequence line first
      if writingSequence:
        output_handle.write('\n')
        writingSequence = False;
      output_handle.write(line)
    else:
      writingSequence = True
      translated_line = line.rstrip().translate(uppercase_translation, dotdash)
      output_handle.write(translated_line)
# Terminate the last sequence line
if writingSequence:
  output_handle.write('\n')
