#!/usr/bin/python

import re, sys, commands

# ------------------------------------------------------------------------------
subfamily_re = re.compile("%subfamily N*([0-9]*)")

def identify_subfamily_of_sequence(sequence_id, subfam_file):
  cmd = "grep -E '%s|subfam' %s" % (sequence_id, subfam_file)
  status, output = commands.getstatusoutput(cmd)
  lines = output.split('\n')
  subfam_line_idx = -1
  for i in xrange(len(lines)):
    if lines[i][1:len(sequence_id) + 1] == sequence_id:
      subfam_line_idx = i - 1
      break
  if subfam_line_idx < 1:
    print "Could not find seed %s in subfamily file %s" \
            % (sequence_id, subfam_file)
    return ""
  subfam_match = re.match(subfamily_re, lines[subfam_line_idx])
  if subfam_match == None or len(subfam_match.groups()) < 1:
    print "Could not parse subfamily line"
    print lines[subfam_line_idx]
    return ""
  return subfam_match.group(1)
      
def main():
  if len(sys.argv) < 3:
    print "Usage: %s sequence_id subfamily_file" % sys.argv[0]
    sys.exit(0)
  print identify_subfamily_of_sequence(sys.argv[1], sys.argv[2])

if __name__ == '__main__':
  main()
