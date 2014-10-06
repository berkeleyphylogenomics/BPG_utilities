"""Functions for computing pairwise identity between two sequences

pairwise_identity_KS_1(seq0, seq1) =
  # of columns with same uppercase character in seq0 and seq1
  ------------------------------------------------------------
   # of columns with an uppercase character in seq0 or seq1

pairwise_identity_for_structural_alignments(seq0, seq1) =
  # of columns with same uppercase character in seq0 and seq1
  ------------------------------------------------------------
  # of columns with uppercase characters in both seq0 and seq1

"""

__all__ = ["pairwise_identity_KS_1",
"pairwise_identity_for_structural_alignments"]

def pairwise_identity_KS_1(seq0, seq1):
  i = 0
  j = 0
  num_union = 0
  num_intersection = 0
  def go_past_lowercase_characters_in_seq(seq, k):
    while k < len(seq) and (seq[k].islower() or seq[k] == '.'):
      k = k + 1
    return k
  i = go_past_lowercase_characters_in_seq(seq0, i)
  j = go_past_lowercase_characters_in_seq(seq1, j)
  # At the top of this loop seq0[i] is either uppercase or a dash
  # At the top of this loop seq1[j] is either uppercase or a dash
  while i < len(seq0) or j < len(seq1):
    if i < len(seq0):
      c = seq0[i]
      if c.isupper():
        num_union = num_union + 1
        if j < len(seq1) and seq1[j] == c:
          num_intersection = num_intersection + 1
        i = i + 1
        j = j + 1
      elif c == '-':
        if j < len(seq1):
          if seq1[j].isupper():
            num_union = num_union + 1
          i = i + 1
          j = j + 1
        else:
          print "Error in pairwise_identity_KS_1"
          print "# of uppercase characters and dashes in sequences differs"
          print seq0
          print seq1
          raise AssertionError
      else:
        print "Error in pairwise_identity_KS_1"
        print "Char %d, %s, of sequence is neither uppercase nor a dash" \
              % (i, c)
        print seq0
        raise AssertionError
    else:
      print "Error in pairwise_identity_KS_1"
      print "# of uppercase characters and dashes in sequences differs"
      print seq0
      print seq1
      raise AssertionError
    i = go_past_lowercase_characters_in_seq(seq0, i)
    j = go_past_lowercase_characters_in_seq(seq1, j)
  pwid = float(num_intersection)/float(num_union)
  return pwid

def pairwise_identity_for_structural_alignments(seq0, seq1):
  if len(seq0) != len(seq1):
    print "Error in pairwise_identity_for_structural_alignments"
    print "Lengths %d and %d of sequences differ" % (len(seq0), len(seq1))
    print seq0
    print seq1
    raise AssertionError
  i = 0
  num_superposable = 0
  num_identical = 0
  for i in xrange(0,len(seq0)):
    if seq0[i].isupper(): 
      if  seq1[i].isupper():
        num_superposable = num_superposable + 1
        if seq0[i] == seq1[i]:
          num_identical = num_identical + 1
      elif seq1[i] != '-':
        print "Error in pairwise_identity_for_structural_alignments"
        print "Char #%d, %s, of 2nd sequence is neither uppercase nor a dash" \
          % (i, seq1[i])
        print "Sequences:"
        print seq0
        print seq1
        raise AssertionError
    elif seq0[i] != '-':
      print "Error in pairwise_identity_for_structural_alignments"
      print "Char #%d, %s, of 1st sequence is neither uppercase nor a dash" \
        % (i, seq0[i])
      print "Sequences:"
      print seq0
      print seq1
      raise AssertionError
  pwid = float(num_identical)/float(num_superposable)
  return pwid

def pairwise_identity_belvu(seq0, seq1):
  '''Calculate pairwaise identity according to belvu.
  This is similar to the calculation done in the function pairwise_identity_for_structural_alignments above
  but attempted to be more concise.'''
  if len(seq0) != len(seq1):
    print 'Sequence lengths do not match. Exiting ...'
    return
  match = 0
  identity = 0
  for ind, char0 in enumerate(seq0):
    char1 = seq1[ind]
    if not (char0.isupper() and char1.isupper()):
      continue
    match += 1
    if char0 == char1:
      identity += 1
  try:
    pwid = float(identity)/match
  except:
    pwid = 0
  return pwid
      
