from glob import glob
import os
import pickle
import re
import subprocess
import shutil
import string
import StringIO

from Bio import AlignIO

from pfacts003.queued.utils import q_function
from pfacts003.satchmo.forms import SatchmoForm as AppForm, AdvancedForm
from pfacts003.satchmo.views import results as general_results

title = 'SATCHMO-JS'
downloadable_files = ['satchmo_alignment.fasta', 'satchmo_tree.newick',
                      'satchmo.smo', 'input_unaligned.fasta',
                      'satchmo_tree_wo_lengths.newick']
js = {'jquery': True, 'jquery_ui': True}

email_subject = 'SATCHMO-JS server results'

@q_function(title='Setting Up',\
    write_field=('satchmo_fasta','input_unaligned.fasta'))
def setup(func_info):
    """Set up beginning of pipeline

    The data entered into the form field 'satchmo_fasta' is written to
    a file that will, by the end of this method, be called
    'input_seequences.fasta'.
         
    Because phylip format is used in several steps along this pipeline,
    and because it strict phylip format limits the defline to be only
    10 characters, the deflines have been saved and replaced with the
    unique pattern SEQ###### (e.g., SEQ000023).  These deflines will be
    placed back into the final results.
    """

    # Convert defline to >SEQ###### format
    in_file = open('./input_unaligned.fasta', 'rU')
    out_file = open('./input_sequences.fasta', 'w')

    adjusted_seqs = []
    seq_dict = {}
    counter = 0
    for line in in_file.readlines():
        if line.startswith('>'):
            key = 'SEQ%06d' % counter
            seq_dict[key] = line
            adjusted_seqs.append(">%s\n" % key)
            counter += 1
        else:
            adjusted_seqs.append(line)
    # Now, the adjusted_msa list contains the adjusted alignment and 
    # seq_dict is the dictionary to translate back

    # Write seq_dict first, in case error occurs before it gets written
    p_file = open('./seq_dict.pkl', 'w')
    p_file.write(pickle.dumps(seq_dict))
    p_file.close()

    # Writing sequences that have deflines replaced
    out_file.writelines(adjusted_seqs)
    out_file.close()
    in_file.close()


@q_function(title='Aligning With MAFFT',\
    command_str='mafft', command_args=['--maxiterate', '1', 'input_sequences.fasta'], \
    command_out='mafft_alignment.fasta', display_outfile='mafft.err')   
def mafft(func_info):
    """Run the MAFFT (Multiple Alignment) algorithm input sequences

    The MAFFT algorithm produces a Multiple Sequence Alignment (MSA)
    that is used as input to the Quick Tree algorithm (next step).
   
    When MAFFT is completed, the MAFFT MSA is examined to find the
    minimum pairwise identity between any two sequences in the MSA. This
    value, plus 10% is used as the cut-off for Kerf.
    
    For example, if the MSA has a minimum pairwise identity of 20%, then
    the Kerf cut (next step) will divide the tree into subtrees such
    that no pair of sequences within any subtree has less than 30%
    identity.

    If minimum pairwise identity of the MSA created by MAFFT is above
    90%, then, the midpoint between the maximum pairwise identity and
    minimum pairwise identity is used to create the Kerf cut.
    
    For example, if the MAFFT MSA has a minimum pairwise identity of
    91%, and a maximum pairwise identity of 99%, then the Keft cut
    (again, next step) will divide the trees into subtrees such that no
    pair of sequences within any subtree has less than 95% identity.

    Because the QuickTree algorithm needs its input in stockholm format,
    we also use Biopython to convert the MAFFT results to stockholm
    format.
 
    This step takes the input_sequences.fasta file previously generated
    and produces the file 'mafft_alignment.stock'.

    The command line used is to produce the fasta file (before
    converting to stockholm format) is:
    mafft input_sequences.fasta > mafft_alignment.fasta
    """
    
    func_info.logger.info('Calculating alignment statistics')

    mafft_fa = os.path.abspath('mafft_alignment.fasta')
    mafft_stockholm = os.path.abspath('mafft_alignment.stock')

    in_file = open(mafft_fa, 'r')
    out_file = open(mafft_stockholm, 'w')

    try:
        align = AlignIO.read(in_file, 'fasta')
        # take the minimum value of the list of pwid's
        pwids = [reduce(lambda m,n: 100.0*sum(m)/sum(n),
            # convert the list of match:align pairs into two lists
            zip(*((
                # do the characters match?
                int(x == y and x != '-'),
                # is the column occupied?
                int(x != '-' or y != '-'),
            # loop over all positions
            ) for x,y in zip(str(a.seq), str(b.seq))))) \
            # loop over all pairs of sequences
            for i,a in enumerate(align) for b in align[i+1:] \
        ]
        minpwid = min(pwids)
        maxpwid = max(pwids)
        func_info.kerf_cut_percent = max((35, minpwid+10))
        func_info.logger.info('Converting MAFFT output to Stockholm format')
        AlignIO.write([align], out_file, 'stockholm')
    except Exception, e:
        func_info.logger.exception(e)
        in_file.close()
        out_file.close()
        return False

    in_file.close()
    out_file.close()
    return func_info.returncode


@q_function(title='Calculating Tree')
def quicktree(func_info):
    """Produce Newick Tree via QuickTree algorithm from MAFFT alignment

    The Quick Tree algorithm takes the previously created
    'mafft_alignment.stock' file, and creates a newick tree file
    'quicktree.newick'.

    The command line used is:
    quicktree mafft_alignment.stock > quicktree.newick

    TODO: Add FastTree to these notes; remove the .stock file creation
        when unnecessary.
    """
    
    treebuilder = func_info.fields['treebuilder']
    func_info.treebuilder = treebuilder

    if treebuilder == 'quicktree':
        treeargs = ['quicktree', 'mafft_alignment.stock']
    else:
        treeargs = ['FastTree', '-nosupport', 'mafft_alignment.fasta']

    fout = open('%s.newick' % treebuilder, 'w')
    ferr = open('%s.err' % treebuilder, 'w')
    func_info.returncode = subprocess.call(treeargs, stdout=fout, stderr=ferr)
    fout.close()
    ferr.close()

    return func_info.returncode


@q_function(title='Cutting With KERF', display_outfile='kerf_output.out')
def kerf(func_info):
    """Create cuts of the tree with minimum pairwise percent identity

    The previously generated MAFFT multiple sequence alignment and
    QuickTree generated tree is taken with the minimum Pairwise Identity
    given in the MAFFT step.  With these inputs, this step outputs cuts
    of the tree into the fewest possible number of subtrees so that no
    pair in any subtree has less than the minimum pairwise identity
    previously stipulated.
       
    The definition for pairwise identity is the ratio of the number of
    aligned pairs that agree exactly over the number of aligned columns
    in the MSA.  When calculating the percent identity of two sequences
    in an MSA, columns in which both of these sequences have a gap
    should not be counted in the numerator nor the denominator of the
    above ratio.
            
    Output: A set of sub-alignments (in the same format as the input),
    and sub-trees: One file for each sub-alignment, and one file for
    each subtree.

    The comand line used is: 
    kerf -1 <optimal cut given by mafft step> quicktree.newick
         mafft_alignment.fasta > kerf_output.out

    TODO: update description to note FastTree.newick
    """

    fout = open('kerf_output.out', 'w')
    ferr = open('kerf.err', 'w')
    func_info.logger.info('Cutting at %i%%' % func_info.kerf_cut_percent)
    func_info.returncode = subprocess.call(
        ['kerf', '-1', str(func_info.kerf_cut_percent),
            '%s.newick' % func_info.treebuilder,
            'mafft_alignment.fasta'],
        stdout=fout, stderr=ferr,
    )
    fout.close()
    ferr.close()

    return func_info.returncode


@q_function(title='Converting KERF Output To SATCHMO Input')
def convert_kerf(func_info):
    """Convert the KERF output into a format that SATCHMO can understand

    The SATCHMO program expects its input in a specific format as
    follows: The input is a directory that contains a file titled
    "summary-partition" and a set of subdirectories. The
    'summary-partition' file contents specify the names of the
    subdirectories. In our context, these subdirectories will contain
    the sub-alignment files previously created from KERF.
                 
    This step converts the previously created create into this format
    that SATCHMO expects.
                      
    Additionally, there is a problem that alignments given to SATCHMO
    were based on a cut of the MAFFT MSA into sub-alignments. This
    could result in many of the sub-alignments containing columns that
    are entirely gapped.  

    The solution is to mask the sub-alignments prior to submitting the
    sub-alignments to SATCHMO. In this step, we mask columns that have
    an excessive fraction of gap characters, by turning the amino acids
    in those columns to lower case letters (and corresponding dashes in
    those positions to dots). SATCHMO should interpret these positions
    as generated in an HMM insert state.
                                
    This is done by two passes of the remove GappyColumns script.
    #1: removeGappyColumns 100 <intermediate_align>
    #2: removeGappyColumns -i 70 <intermediate_align> > <filename>
                                                     
    This will turn any columns which are at least 70% gaps into
    lowercase letters and dots.
                                                          
    Neither version will do anything to columns which are 100% dots,
    however, this is not addressed since the MAFFT algorithm does not
    create files with dots in them.
    """

    satchmo_input_dir = os.path.abspath('satchmo_input')
    os.makedirs(satchmo_input_dir)
    summary_partition = open(os.path.join(satchmo_input_dir,
                             'summary-partition'), 'w')
    summary_partition.write("Header Line\n")

    for f in glob('sub*.fa'):

        subfam_file = os.path.splitext(f)[0]
        summary_partition.write("%s\n" % subfam_file)
        os.makedirs(os.path.join(satchmo_input_dir, subfam_file))

        handle = open('deleteme', 'w')
        subprocess.Popen(['removeGappyColumns', '100', f],
            stdout=handle).wait()
        handle.close()

        handle = open(os.path.join(satchmo_input_dir, subfam_file,
            'acceptedseqs.a2m'), 'w')
        subprocess.Popen(['removeGappyColumns', '-i', '70', 'deleteme'],
            stdout=handle).wait()
        handle.close()
        os.unlink('deleteme')
        
    summary_partition.close()


@q_function(title='Executing SATCHMO', display_outfile='satchmo.err')
def satchmo(func_info):
    """Run the SATCHMO algorithm on the inputs given

    The SATCHMO executable is executed upon the original sequences
    submitted and is using the jump start mode with the directory
    structure previously given by KERF (However, the trees from KERF are
    not currently being used).
     
    The final results that are produced are:
        * A Satchmo alignment ("satchmo_alignment.fasta")
        * A Satchmo Tree ("satchmo_tree.newick"), and
        * A Satchmo SMO file that is used with the SATCHMO Viewer on the
          Windows platform ("satchmo.smo").

    The command line used is:
    satchmo -satchmo input_sequences.fasta -bs satchmo_input
            -out satchmo_alignment.fasta -sv satchmo.smo
            -phy satchmo_tree.newick
    """

    handle_o = open('satchmo.out', 'w')
    handle_e = open('satchmo.err', 'w')
    subprocess.Popen(['satchmo',
        '-satchmo', 'input_sequences.fasta',
        '-bs', 'satchmo_input/',
        '-out', 'satchmo_alignment.fasta',
        '-sv', 'satchmo.smo',
        '-phy', 'satchmo_tree.newick',
        '-minaff', str(func_info.fields['minaff']),
    ], stderr=handle_e, stdout=handle_o).wait()
    handle_e.close()
    handle_o.close()

    f = open('satchmo.err')
    has_errored = not bool(f.read().count(\
                                     'All HMMs have zero length, quitting.\n'))
    f.close()

    return has_errored

@q_function(title='Converting Data For RAxML', command_str='prettyalign',
    command_args=['./satchmo_alignment.fasta', '-f', '-m0'],
    command_out='satchmo_alignment.afa')
def convert_satchmo(func_info):
    """Convert SATCHMO format to a form that RAxML can use

    Because edge lengths are not created with the SATCHMO algorithm, the
    RAxML algorithm is used to create the edge lengths for the SATCHMO
    already generated tree.
     
    Additionally, there were still dot characters in the alignments
    given as inputs in this pipeline. Characters in columns that have
    lowercase letters and dots are not considered aligned to each other.
          
    These characters were masked prior to giving them to the RaXML
    algorithm.  The prettyalign program was used to remove the columns
    with dots and lower case letters.
               
    The output of this step is the file 'satchmo_alignment.phy'
    """
    in_file = open('./satchmo_alignment.afa', 'rU')
    out_file = open('./satchmo_alignment.phy', 'w')

    # Convert alignment from fasta to phylip format and write file
    try:
        alignments = AlignIO.parse(in_file, 'fasta')
        AlignIO.write(alignments, out_file, 'phylip')
    except Exception, e: 
        func_info.logger.exception(e)
        in_file.close()
        return False
    in_file.close()
    out_file.close()
    return func_info.returncode


@q_function(title='Calculating Branch Lengths With RAxML',
    command_str='raxmlHPC', command_args=[
        '-f', 'e', '-t', './satchmo_tree.newick', '-m', 'PROTGAMMAJTT', '-s',
        './satchmo_alignment.phy', '-n', 'post_satchmo_raxml', '-w',
        os.path.abspath('prefix')
    ], command_out='./raxml.out')
def raxml(func_info):
    """Run RAxML algorithm against the SATCHMO tree. 

    The command line issued is:
        raxmlHPC -f e -t ./satchmo_tree.newick -m PROTGAMMAJTT -s
                ./satchmo_alignment.phy -n post_satchmo_raxml
                -w <current working directory>/prefix

    At the end of this step, the file 'satchmo_tree.newick', which was
    the Newick tree output by SATCHMO, is moved to
    'satchmo_tree_wo_lengths.newick', and the Newick tree including branch
    lengths output by RAxML is moved to 'raxml_output_tree.newick'.

    The final result of this step is the file
    'raxml_output_tree.newick', an unrooted tree including branch
    lengths.
    """
    shutil.move(os.path.abspath('satchmo_tree.newick'),
                os.path.abspath('satchmo_tree_wo_lengths.newick'))
    shutil.move(os.path.abspath('prefixRAxML_result.post_satchmo_raxml'),
                os.path.abspath('raxml_output_tree.newick'))
    return func_info.returncode


@q_function(title="Incorporating Branch Lengths",
  command_str="incorporate_branch_lengths_from_unrooted_tree.py", 
  command_args = [os.path.abspath('satchmo_tree_wo_lengths.newick'),
                  os.path.abspath('raxml_output_tree.newick')],
  command_out = './satchmo_tree.newick')
def incorporate_branch_lengths(func_info):
  """Incorporate branch lengths from the RAxML tree into SATCHMO tree.

  The command line issued is:
    incorporate_branch_lengths_from_unrooted_tree.py
    <current_working_directory>/satchmo_tree_wo_lengths.newick 
    <current_working_directory>/raxml_output_tree.newick

  The final result of this step is the file 'satchmo_tree.newick', a
  rooted tree created by SATCHMO that includes branch lengths
  optimized by RAxML.
  """
  return func_info.returncode

@q_function(title='Cleaning Up')
def clean_up(func_info):
    """Convert SEQ###### in deflines and trees back to original format

    The deflines of the original input were modified to be of a format
    similar to this example:

        > SEQ000004

    The original deflines have been retrieved and the sequences in
    question replaced, as they were in the original input.
    """

    # Retrieve original sequences dictionary
    f = open('./seq_dict.pkl', 'r')
    seq_dict = pickle.load(f)
    f.close()

    # Since newick trees use the following symbols, it makes sense that our
    # identifiers should not include these symbols so as to not interrupt the
    # newick format.
    #  : , ) ( ; ] [
    newick_trans = string.maketrans(":,)(;][", '_'*7)


    def replace_seqs(filename, handler):

        in_file = open(filename, 'r')
        modified = handler(in_file.readlines())
        in_file.close()

        out_file = open(filename, 'w')
        out_file.writelines(modified)
        out_file.close()

    def replace_fasta(iterable):
        defline_re = re.compile(r'^>(SEQ\d{6})$')
        results = []

        for line in iterable:
            replacement = defline_re.match(line)
            if replacement:
                results.append(seq_dict[replacement.groups()[0]])
            else:
                # If, for whatever reason, no SEQ -- leave original line
                results.append(line)

        return results

    def replace_inline_seqs(iterable):
        seq_re = re.compile(r'(SEQ\d{6})')

        results = []
        for line in iterable:
            replacement = seq_re.split(line)
            for item in replacement:
                if seq_re.match(item):
                    results.append(seq_dict[item].split(' ')[0][1:].\
                                    translate(newick_trans))
                else:
                    results.append(item)

        return results

    for filename in ['./satchmo_alignment.afa', './satchmo_alignment.fasta',
        './satchmo.smo']:
        replace_seqs(filename, replace_fasta)

    for filename in ['./satchmo_tree.newick',
        './satchmo_tree_wo_lengths.newick']:
        replace_seqs(filename, replace_inline_seqs)

@q_function(title='Parsing SATCHMO Viewer File', command_str='parse_smo.py',
    command_args=[os.path.abspath('.')])
def parse_smo(func_info):
  """Parse the SATCHMO Viewer file to find offsets of the subalignments.

  The command line issued is:
    parse_smo.py <current working directory>

  This command expects the files satchmo_tree.newick and satchmo.smo to
  be in the current working directory.  It parses them both to find the
  offsets within the file satchmo.smo of the subalignments associated
  with nodes of the tree satchmo_tree.newick, and their sizes.  It
  creates a dictionary whose keys are left ids (from modified pre-order
  tree traversal) of nodes in the tree, and whose values are pairs of an
  alignment offset and its size in bytes.

  The final result of this step is the file
  'alignment_offset_of_left_id.pkl', containing the dictionary
  serialized with cPickle.
  """
  return func_info.returncode
    
