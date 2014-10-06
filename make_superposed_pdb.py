#!/usr/bin/python

import os, sys

if len(sys.argv) < 3:
  print "Usage: %s <seed_id> <model_file>" % sys.argv[0]
  sys.exit(0)

seed_id = sys.argv[1]
if len(seed_id) < 4:
  print "<seed_id> must be a pdb identiier"
  sys.exit(0)
chainless_id = seed_id[0:4]
model_file = sys.argv[2]
model_components = os.path.split(model_file)
model_dir = model_components[0]
model_basename = os.path.splitext(model_components[1])[0]
true_structure_path \
  = "/home/ruchira/SHMM-SHMM/single_seqs/%s/pdb%s.ent" % (seed_id, chainless_id)
chain = ' '
if len(seed_id) >= 5:
  chain = seed_id[4]
rotation_file = open(os.path.join(model_dir, "%s_rotation" % seed_id), "r")
rotation_matrix = [[float(x) for x in line.rstrip().split(',')] 
                    for line in rotation_file.readlines()]
num_rows = len(rotation_matrix)
num_cols = len(rotation_matrix[0])
translation_file = open(os.path.join(model_dir, "%s_translation" % seed_id),
                        "r")
translation_vector = [float(x) 
                      for x in translation_file.readline().rstrip().split(',')]
def apply_rotation_matrix_on_left(vector):
  return ([sum(rotation_matrix[i][j] * vector[j] for j in range(len(vector)))
           for i in range(num_rows)])
def get_superposed_model_pos(model_pos):
  rotated_model_pos = apply_rotation_matrix_on_left(model_pos)
  return (rotated_model_pos[0] + translation_vector[0],
          rotated_model_pos[1] + translation_vector[1],
          rotated_model_pos[2] + translation_vector[2])
def superpose_model_line(line):
  model_pos = [float(line[30:38]),float(line[38:46]),float(line[46:54])]
  (x,y,z) = get_superposed_model_pos(model_pos)
  superposed_line = "%sB%s%8.3f%8.3f%8.3f%s" \
                      % (line[0:21], line[22:30], x, y, z, line[54:])
  return(superposed_line)
true_structure = open(true_structure_path, "r")
true_structure_lines = true_structure.readlines()
true_structure.close()
model = open(model_file, "r")
model_lines = model.readlines()
model.close()
superposed_model_lines = [superpose_model_line(line)
                          for line in model_lines if line[0:4] == 'ATOM']
last_superposed_model_line = \
  superposed_model_lines[len(superposed_model_lines)-1]
ter_num = int(last_superposed_model_line[6:11]) + 1
ter_res = last_superposed_model_line[17:20]
ter_res_seq = last_superposed_model_line[22:26]
ter_line = \
 'TER   %5d      %s B%s                                                      ' \
 % (ter_num, ter_res, ter_res_seq) + '\n'
superposed_model_file = os.path.join(model_dir, 
                                      "superposed_%s.pdb" % model_basename)
f = open(superposed_model_file, "w")
seeking_atom = True
printing_chain = False
for line in true_structure_lines:
  if seeking_atom:
    if line[0:4] == 'ATOM':
      seeking_atom = False
      printing_chain = True
      if line[21] == chain:
        f.write("%sA%s" % (line[0:21], line[22:]))
    else:
      f.write(line)
  elif printing_chain:
    if line[0:3] == 'TER':
      f.write(line)
      printing_chain = False
      for superposed_model_line in superposed_model_lines:
        f.write(superposed_model_line)
      f.write(ter_line)
    elif line[21] == chain:
      f.write("%sA%s" % (line[0:21], line[22:]))
  else:
    f.write(line)
f.flush()
f.close()
