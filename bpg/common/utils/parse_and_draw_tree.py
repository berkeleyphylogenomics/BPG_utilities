#!/usr/bin/env python

import sys, string, re, cairo, math, cPickle
from pfacts003.phylofacts.models import UniProt, UniProtEC, EC, UniProtGO, \
    GO_Term, GO_EvidencePriority, UniProtPDB_Chain, PDB_Chain

id_re = re.compile('[^(:),]*')
branch_length_re = re.compile('[0-9]\.[0-9e\-]*')
EXTENTS_INDEX_OF_WIDTH_ELEMENT = 2
EXTENTS_INDEX_OF_HEIGHT_ELEMENT = 3
species_re = re.compile('OS=')

def rgb_of_hsv(h, s, v):
  h_index = int(math.floor(h/60)) % 6
  f = h / 60.0 - math.floor(h/60)
  p = v * (1 - s)
  q = v * (1 - f * s)
  t = v * (1 - (1 - f) * s)
  if h_index == 0:
    return (v, t, p)
  elif h_index == 1:
    return (q, v, p)
  elif h_index == 2:
    return (p, v, t)
  elif h_index == 3:
    return (p, q, v)
  elif h_index == 4:
    return (t, p, v)
  elif h_index == 5:
    return (v, p, q)

def modified_header(header, includeDBInfo = False):
  escaped_header = header
  escaped_header = escaped_header.replace('(', '%28')
  escaped_header = escaped_header.replace(')', '%29')
  escaped_header = escaped_header.replace(':', '%3A')
  escaped_header = escaped_header.replace(',', '%2C')
  escaped_header = escaped_header.replace(';', '%3B')
  header_words = escaped_header.split()
  id_components = header_words[0].split('|')
  if len(id_components) == 3:
    uniprot_accession = id_components[1]
    for end_description_word_index in range(1,len(header_words)):
      if len(header_words[end_description_word_index]) >= 3 and \
          header_words[end_description_word_index][2] == '=':
        break
    if len(header_words[1].split('|')) == 3:
      description = ' '.join(header_words[2:end_description_word_index])
    else:
      description = ' '.join(header_words[1:end_description_word_index])
    if includeDBInfo:
      uniprot = UniProt.objects.get(accession = uniprot_accession)
      if uniprot.in_swissprot_f:
        description += ' SwissProt'
      uniprot_ecs = UniProtEC.objects.filter(uniprot = uniprot, 
                                          is_in_brenda_f = True)
      if uniprot_ecs:
        description += ' EC:%s' % ';'.join([uniprot_ec.ec.__unicode__() 
                                            for uniprot_ec in uniprot_ecs])
      uniprot_gos = UniProtGO.objects.filter(uniprot = uniprot, 
                                  go_evidence__priority__lte =
                                  GO_EvidencePriority.experimental_priority)
      if uniprot_gos:
        description += ' GO:%s' % ';'.join([uniprot_go.go_term.acc
                                            for uniprot_go in uniprot_gos])
                                            
      uniprot_pdb_chains = UniProtPDB_Chain.objects.filter(uniprot = uniprot)
      if uniprot_pdb_chains:
        description += ' PDB:%s' % ';'.join(['%s%s' %
                                          (uniprot_pdb_chain.pdb_chain.pdb.id,
                                          uniprot_pdb_chain.pdb_chain.chain_id)
                                          for uniprot_pdb_chain 
                                          in uniprot_pdb_chains])
    species = ''
    found_species = False
    for start_species_index in range(end_description_word_index,
                                      len(header_words)):
      if len(header_words[start_species_index]) >= 3 and \
          header_words[start_species_index][0:3] == 'OS=':
        found_species = True
        break
    if found_species:
      end_species_index = len(header_words)
      for end_species_index in range(start_species_index + 1, 
                                      len(header_words)):
        if len(header_words[end_species_index]) >= 3 and \
            header_words[end_species_index][2] == '=':
          break
      species = ' '.join(header_words[start_species_index:
                                        end_species_index])[3:]
      return '%s %s--%s' % (uniprot_accession, species, description)
    return '%s No species--%s' % (uniprot_accession, description)
  else:
    return escaped_header

class node:
  def __init__(self, seqid=None):
    self.parent = None
    self.children = set()
    self.branch_length = -1.0
    self.likelihood = -1.0
    self.seqid = seqid
    if seqid:
      self.contained_seqids = set([seqid])
    else:
      self.contained_seqids = set()

  def addChild(self, child):
    self.children.add(child)
    child.parent = self
    self.contained_seqids = self.contained_seqids | child.contained_seqids

  def readFromTreeString(self, tree_string, i):
    if i >= len(tree_string):
      return len(tree_string)
    if tree_string[i] == '(':
      while tree_string[i] != ')':
        child = node()
        i = child.readFromTreeString(tree_string, i+1)
        self.addChild(child)
      i += 1
      m = branch_length_re.match(tree_string[i:])
      if m:
        likelihood = m.group(0)
        self.likelihood = float(likelihood)
        i += len(likelihood)
    else:
      m = id_re.match(tree_string[i:])
      seqid = m.group(0)
      i += len(seqid)
      self.seqid = seqid
      self.contained_seqids = set([self.seqid])
    if i < len(tree_string) and tree_string[i] == ':':
      i += 1
      branch_length = branch_length_re.match(tree_string[i:]).group(0)
      self.branch_length = float(branch_length)
      i += len(branch_length)
    return i

  def printSelf(self, header_of_seqid, indentLevel = 0, 
                output_file = sys.stdout, write_comma = False):
    if self.seqid:
      for i in xrange(indentLevel):
        output_file.write(' ')
      output_file.write("%s" % header_of_seqid[self.seqid])
    else:
      for i in xrange(indentLevel):
        output_file.write(' ')
      output_file.write("(")
      output_file.write('\n')
      children = list(self.children)
      for i in range(len(children)):
        children[i].printSelf(header_of_seqid, indentLevel = indentLevel + 1, 
                              output_file = output_file, 
                              write_comma = (i < len(children) - 1))
      for i in xrange(indentLevel):
        output_file.write(' ')
      output_file.write(")")
    if self.likelihood >= 0.0:
      output_file.write("%g" % self.likelihood)
    if self.branch_length >= 0.0:
      output_file.write(":%g" % self.branch_length)
    if write_comma:
      output_file.write(",")
    output_file.write('\n')

  def getMaxDepth(self):
    if self.seqid:
      return 0.0
    else:
      maxDepth = 0.0
      for child in self.children:
        childDepth = child.branch_length + child.getMaxDepth()
        if childDepth > maxDepth:
          maxDepth = childDepth
      return maxDepth
        
  def drawSelf(self, ctx, x, y, widthFactor, verticalTick, header_of_seqid):
    if len(self.children) == 0:
      return
    upper_child = list(self.children)[0]
    lower_child = list(self.children)[1]
    xUpper = x + upper_child.branch_length * widthFactor
    xLower = x + lower_child.branch_length * widthFactor
    if len(upper_child.children) == 0:
      yUpper = y - 0.5 * verticalTick
      heightForLowerUpperGrandChild = -1.0
    else:
      lower_upper_grandchild = list(upper_child.children)[1]
      heightForLowerUpperGrandChild = \
          (len(lower_upper_grandchild.contained_seqids)) \
              * verticalTick
      yUpper = y - heightForLowerUpperGrandChild
    if len(lower_child.children) == 0:
      yLower = y + 0.5 * verticalTick
      heightForUpperLowerGrandChild = -1.0
    else:
      upper_lower_grandchild = list(lower_child.children)[0]
      heightForUpperLowerGrandChild = \
          (len(upper_lower_grandchild.contained_seqids)) \
              * verticalTick
      yLower = y + heightForUpperLowerGrandChild
    ctx.move_to(x,y)
    ctx.line_to(x,yUpper)
    ctx.stroke()
    ctx.move_to(x, yUpper)
    ctx.line_to(xUpper,yUpper)
    ctx.stroke()
    if upper_child.seqid:
      ctx.move_to(xUpper,yUpper)
      ctx.show_text(' ' + header_of_seqid[upper_child.seqid])
    ctx.move_to(x,y)
    ctx.line_to(x,yLower)
    ctx.stroke()
    ctx.move_to(x,yLower)
    ctx.line_to(xLower,yLower)
    ctx.stroke()
    if lower_child.seqid:
      ctx.move_to(xLower,yLower)
      ctx.show_text(' ' + header_of_seqid[lower_child.seqid])
    upper_child.drawSelf(ctx, xUpper, yUpper, widthFactor, verticalTick,
                          header_of_seqid)
    lower_child.drawSelf(ctx, xLower, yLower, widthFactor, verticalTick,
                          header_of_seqid)

  def getContainedLeaves(self):
    if self.seqid:
      return set([self])
    else:
      contained_leaves = set()
      for child in self.children:
        contained_leaves = contained_leaves | child.getContainedLeaves()
      return contained_leaves

def main():
  if len(sys.argv) < 4:
    print "Usage: %s <Newick_tree_file> " \
          % sys.argv[0] + " <idmap_file> <output_file> [--include_db_info]"
    sys.exit(0)

  includeDBInfo = False
  if len(sys.argv) >= 5 and sys.argv[4] == "--include_db_info":
    includeDBInfo = True
  trivial_translation = string.maketrans('', '')
  f = open(sys.argv[1])
  tree_str = f.read().translate(trivial_translation, string.whitespace)
  f.close()
  root = node()
  root.readFromTreeString(tree_str, 0)
  f = open(sys.argv[2])
  raw_header_of_seqid = cPickle.load(f)
  f.close()
  totalDepth = root.getMaxDepth()
  if len(root.children) == 0:
    print "Tree has no children, exiting."
    sys.exit(0)
  upper_child = list(root.children)[0]
  verticalTick = 1.0 / (len(root.contained_seqids) + 1.0)
  startY = (len(upper_child.contained_seqids) + 0.5) * verticalTick

  HEIGHT = int(math.ceil(16.0 / verticalTick))
  WIDTH = int(math.ceil(HEIGHT / 1.61803399))

  surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
  ctx = cairo.Context (surface)
  ctx.scale (WIDTH/1.0, HEIGHT/1.0)

  ctx.select_font_face('sans-serif')
  ctx.set_font_size(verticalTick)
  max_header_width = 0.0
  max_header_height = 0.0
  header_of_seqid = {}
  for seqid in raw_header_of_seqid:
    header_of_seqid[seqid] = modified_header(raw_header_of_seqid[seqid],
                                            includeDBInfo = includeDBInfo)
    extents = ctx.text_extents(' ' + header_of_seqid[seqid])
    if extents[EXTENTS_INDEX_OF_WIDTH_ELEMENT] > max_header_width:
      max_header_width = extents[EXTENTS_INDEX_OF_WIDTH_ELEMENT]
    if extents[EXTENTS_INDEX_OF_HEIGHT_ELEMENT] > max_header_height:
      max_header_height = extents[EXTENTS_INDEX_OF_HEIGHT_ELEMENT]
  withinMargin = 1.0 - 3.0 * verticalTick - (max_header_width / WIDTH)
  widthFactor = withinMargin / totalDepth
  ctx.set_source_rgb(1.0, 1.0, 1.0)
  ctx.rectangle (0, 0, 1.0, 1.0)
  ctx.fill()
  ctx.set_source_rgb(0.0, 0.0, 0.0)
  ctx.set_line_width(min(0.005, verticalTick * 0.125))
  root.drawSelf(ctx, verticalTick, startY, widthFactor, verticalTick,
                header_of_seqid)
  outf = open("%s.newick" % sys.argv[3], 'w')
  root.printSelf(header_of_seqid, output_file = outf)
  outf.close()
  surface.write_to_png("%s.png" % sys.argv[3])

if __name__ == '__main__':
  main()
