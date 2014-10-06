#!/usr/bin/env python

import os, sys, glob, re, string

def main():
  cur_dir = os.getcwd()
  dir0, dir1, dir2 = cur_dir.split('/')[-3:]
  if (dir0, dir1, dir2) != ('pfacts003', 'templates', 'kegg'):
    print "Please execute this program in your subversion working directory pfacts003/templates/kegg to generate the KEGG map templates."
    sys.exit(0)
  kegg_map_template_paths \
      = glob.glob('/clusterfs/ohana/external/KEGG/current/map*.html')
  mapdata_start_re = re.compile('<map name="mapdata">')
  mapdata_end_re = re.compile('</map>')
  area_rect_re = re.compile('<area shape=rect')
  coord_re = re.compile('coords=(\d+),(\d+),(\d+),(\d+)')
  ec_re = re.compile('(\d+)\.([0-9\-]+)\.([0-9\-]+)\.([0-9\-]+)')
  dotdash_translation = string.maketrans('.-', '__')
  for path in kegg_map_template_paths:
    inf = open(path)
    inf_lines = inf.readlines()
    inf.close()
    start_index = -1
    end_index = -1
    in_map = False
    coords_of_ec = {}
    for i in xrange(len(inf_lines)):
     if mapdata_start_re.match(inf_lines[i]):
       in_map = True
       start_index = i + 1
     elif mapdata_end_re.match(inf_lines[i]):
       in_map = False
       end_index = i
     elif in_map:
       if area_rect_re.match(inf_lines[i]):
         ec_match = ec_re.search(inf_lines[i])
         if ec_match:
           coord_match = coord_re.search(inf_lines[i])
           ec = ec_match.group(0)
           if ec not in coords_of_ec:
            coords_of_ec[ec] = set()
           coords_of_ec[ec].add((coord_match.group(1),
                                              coord_match.group(2),
                                              coord_match.group(3),
                                              coord_match.group(4)))
    if len(coords_of_ec) > 0:
      outf = open(os.path.split(path)[1], "w")
      outf.write('{% extends "kegg/map_base.html" %}\n\n')
      outf.write('{% block node_overlay %}\n')
      for ec in coords_of_ec:
        ec_string = ec.translate(dotdash_translation, '')
        coords_num = 0
        for coord in coords_of_ec[ec]:
          outf.write('  .node_%s_instance_%d\n' % (ec_string, coords_num))
          outf.write('    { position: absolute;\n')
          outf.write('      left: %spx; top: %spx;\n' % (coord[0], 
                                                        coord[1]))
          outf.write('      width: 46px; height: 17px; ')
          outf.write('z-index:1; opacity:0.5; background: #0000FF;\n')
          outf.write('    }\n')
          coords_num += 1
      outf.write('{% endblock %}\n')
      outf.write('{% block mapdata %}\n')
      for i in xrange(start_index, end_index):
        outf.write(inf_lines[i])
      outf.write('{% endblock %}\n')
      outf.write('{% block color_nodes %}\n')
      for ec in coords_of_ec:
        ec_string = ec.translate(dotdash_translation, '')
        outf.write('  {% if ')
        outf.write('nodes_to_color.%s' % ec_string)
        outf.write(' %}\n')
        coords_num = 0
        for coord in coords_of_ec[ec]:
          outf.write('    <div class="node_%s_instance_%d"></div>\n' 
                      % (ec_string, coords_num))
          coords_num += 1
        outf.write('  {% endif %}\n')
      outf.write('{% endblock %}\n')
      outf.close()

if __name__ == '__main__':
  main()
