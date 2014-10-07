import cPickle
import os
import re, string
import StringIO

from Bio import AlignIO, SubsMat
from Bio.SubsMat import MatrixInfo
from Bio.pairwise2 import dictionary_match
from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import render_to_response
from Bio.SeqUtils import CheckSum

from parse_smo import find_alignment_offset_of_left_id, get_seguids_of_ids

js = {'jquery': True, 'jquery_ui': True, 'expandable': True}

blosum62_of_residues = dictionary_match(SubsMat.SeqMat(MatrixInfo.blosum62))

wrapwidth = 80

def index(request):
    return HttpResponseRedirect('/q/satchmo/')


def results(request, work_path, response_dict):
    pickle_path = os.path.join(work_path,
                              'alignment_offset_of_left_id.pkl')
    if os.path.exists(pickle_path):
      f = open(pickle_path)
      alignment_offset_of_left_id = cPickle.load(f)
      f.close()
    else:
      alignment_offset_of_left_id = find_alignment_offset_of_left_id(work_path)
    pickle_path = os.path.join(work_path,
                                'ids_of_seguid.pkl')
    if os.path.exists(pickle_path):
      f = open(pickle_path)
      ids_of_seguid = cPickle.load(f)
      f.close()
    else:
      seguids, ids_of_seguid = get_seguids_of_ids(work_path)
    left_id = 1
    if 'left_id' in request.GET:
      try:
        left_id = int(request.GET['left_id'].strip())
        if left_id < 1:
          left_id = 1
      except ValueError:
        left_id = 1
    alignments = []
    if left_id != 1:
      if left_id in alignment_offset_of_left_id:
        offset, num_bytes = alignment_offset_of_left_id[left_id]
        f = open(os.path.join(work_path, "satchmo.smo"))
        f.seek(offset)
        fake_f = StringIO.StringIO(f.read(num_bytes))
        f.close()
        alignments = list(AlignIO.parse(fake_f, "fasta"))
        fake_f.close()
      else:
        left_id = 1
    else:
      f = open(os.path.join(work_path, 'satchmo_alignment.fasta'))
      alignments = list(AlignIO.parse(f, "fasta"))
      f.close()
    alignment_blocks = []
    if len(alignments) > 0:
      alignment = alignments[0]
      alignment_length = 0
      aligned_column_indices = set()
      alignment_seqs = {}
      first_pass = True
      i = 0
      k = 0
      prev_seguid = ''
      uppercase_translation = string.maketrans(string.lowercase, 
                                                string.uppercase)
      dotdash = '.-'
      print 
      for row in alignment:
        seq = row.seq.tostring()
        if first_pass:
          alignment_length = len(row.seq)
          for j in range(len(seq)):
            if seq[j] == '-' or seq[j].isupper():
              aligned_column_indices.add(j)
          first_pass = False
        alignment_seqs[i] = seq
        unaligned_seq = seq.translate(uppercase_translation, dotdash)
        seguid = CheckSum.seguid(unaligned_seq)
        if seguid in ids_of_seguid and len(ids_of_seguid[seguid]) >= 1:
          if seguid == prev_seguid:
            if k < len(ids_of_seguid[seguid]) - 1:
              k += 1
          else:
            k = 0
          row.id = ids_of_seguid[seguid][k]
          prev_seguid = seguid
        i += 1
      column_conserved_residue = {}
      column_score = {}
      class_of_column = {}
      for j in aligned_column_indices:
        freq_of_residue = {}
        highest_frequency = 0
        most_frequent_residue = ''
        for i in alignment_seqs.keys():
          residue = alignment_seqs[i][j]
          if residue == '-':
            continue
          if residue not in freq_of_residue:
            freq_of_residue[residue] = 0
          freq_of_residue[residue] += 1
          if freq_of_residue[residue] > highest_frequency:
            highest_frequency = freq_of_residue[residue]
            most_frequent_residue = residue
        column_conserved_residue[j] = most_frequent_residue
        num_pairs = 0
        sum_of_scores = 0.0
        for i0 in range(len(alignment_seqs)):
          residue0 = alignment_seqs[i0][j]
          if residue0 != '-':
            for i1 in range(i0):
              residue1 = alignment_seqs[i1][j]
              if residue1 != '-':
                score = blosum62_of_residues(alignment_seqs[i0][j],
                                              alignment_seqs[i1][j])
                sum_of_scores += score
                num_pairs += 1
        if num_pairs > 0:
          column_score[j] = sum_of_scores / num_pairs
          if column_score[j] >= 3:
            class_of_column[j] = 'align_high'
          elif column_score[j] >= 1.5:
            class_of_column[j] = 'align_moderate'
          elif column_score[j] >= 0.5:
            class_of_column[j] = 'align_low'
      
      num_blocks = alignment_length / wrapwidth
      useless_re = re.compile('^[\.-]*$')
      if alignment_length % wrapwidth > 0:
        num_blocks += 1
      back_count = [0 for i in range(len(alignment))]
      for i in range(num_blocks):
        block = []
        for row_no,row in enumerate(alignment):
          seq = row.seq.tostring()
          seq_piece = seq[(i * wrapwidth):((i + 1) * wrapwidth)]
          if useless_re.match(seq_piece):
            continue
          alignment_row = {}
          alignment_row['id'] = row.id
          alignment_row['seq'] = []
          alignment_row['start'] = back_count[row_no]
          alignment_row['stop'] = back_count[row_no] \
            + len(seq_piece.replace('.','').replace('-',''))
          back_count[row_no] = alignment_row['stop']

          for j in xrange(i*wrapwidth,(i+1)*wrapwidth):
            if j < len(seq):
              residue = seq[j]
              spec = {}
              spec['residue'] = residue
              spec['class'] = ''
              if j in aligned_column_indices and residue != '-' and \
                  j in class_of_column:
                if blosum62_of_residues(residue, column_conserved_residue[j]) \
                    >= column_score[j]:
                  spec['class'] = class_of_column[j]
              alignment_row['seq'] = alignment_row['seq'] + [spec]
            else:
                alignment_row['seq'].append(dict(
                    (('residue',' '), ('class', None))
                ))

          block = block + [alignment_row]
        alignment_blocks = alignment_blocks + [block]
    return render_to_response('satchmo/results.html', dict(response_dict,
        relative_path=os.path.basename(work_path).replace('satchmo', '', 1),
        js=js,
        alignment_blocks=alignment_blocks,
        left_id=left_id,
        left_ids_with_alignments = alignment_offset_of_left_id.keys(),
    ))

def about(request):
    return render_to_response('satchmo/about.html',
        {'main_viewing': 'satchmo', 'sub_viewing': 'about',
         'title': 'About SATCHMO',
        })

def help(request):
    return render_to_response("satchmo/help.html", {
        'main_viewing': 'satchmo', 'sub_viewing': 'help',
        'title': 'SATCHMO Help',
        })
