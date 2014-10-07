import re
from textwrap import fill

def protein_char_class(alignment=False, nl=True, allow=None, extra_chars=False, stops=False):
    '''Returns the character class string for the sequence described,
    i.e. '[AC-IK-NP-TVWY\\n]'.

Parameters:
    alignment (default=False) Allow '-' character
    nl (default=True) Allow '\\n' characters
    allow (default=None) String of additional characters allowed in the sequence
    extra_chars (default=False) Allow BOUXZ
    stops (default=False) Allow '*' character
'''
    char_class = r'AC-IK-NP-TVWY'

    if nl:
        char_class += r'\n'

    if allow is not None:
        char_class += allow

    if extra_chars:
        char_class += r'XBOUZ'

    if stops:
        char_class += r'*'

    if alignment:
        char_class += r'\-'

    return '[%s]' % char_class

def protein_fasta_re(alignment=False, allow=None, extra_chars=False, stops=False, ignore_case=False, bad_header_chars=None, max_seqs=None, min_seqs=1):
    '''Returns a compiled regex for the described FASTA format.

Parameters:
    alignment (default=False) Allow '-' character
    allow (default=None) String of additional characters allowed in the sequence
    extra_chars (default=False) Allow BOUXZ
    stops (default=False) Allow '*' character
    ignore_case (default=False) Check sequences insensitive to case
    min_seqs (default=1) Minimum number of sequences allowed
    max_seqs (default=None) Maximum number of sequences allowed
'''
    char_class = protein_char_class(nl=False, alignment=alignment, allow=allow,
        extra_chars=extra_chars, stops=stops)
    char_class_and_ws = protein_char_class(nl=True, alignment=alignment,
        allow=(allow or '') + ' ', extra_chars=extra_chars, stops=stops)
    fasta_part = r'(>%s+)\n\s*(%s%s*)' % (
        # either forbidden characters or dot
        bad_header_chars and (r'[^%s\n]' % bad_header_chars) or '.',
        char_class,
        char_class_and_ws,
    )

    if max_seqs is None:
        mult = r'{%i,}' % min_seqs
    else:
        mult = r'{%i,%i}' % (min_seqs,max_seqs)

    pattern = r'^(%s)%s$' % (fasta_part, mult)

    return re.compile(pattern, ignore_case and re.I)

class FastaError(Exception):
    """FASTA format parsing error
    
Attributes:
    message -- explanation of the error
    num -- the number of the offending line in the FASTA file
    line -- the offending line (optional)
"""

    def __init__(self, message, num, line=''):
        self.message = message
        self.line = "line %i: '%s'" % (num, line)

    def __str__(self):
        return self.message

# convenience functions
def format_protein_fasta_and_count(inp, to_upper=False, wrap=False, alignment=False, allow=None, extra_chars=False, stops=False, ignore_case=False, bad_header_chars=None, max_seqs=1, min_seqs=1, max_length=None):
    '''Check FASTA input for format errors and return a tuple with a string with requested changes and the number of sequences.

Parameters:
    to_upper (default=False) Converts sequences to upper case letters
    wrap (default=False) Reset the line wrapping for sequences
    alignment (default=False) Allow '-' character and enforce equal lengths
    allow (default=None) String of additional characters allowed in the sequence
    extra_chars (default=False) Allow BOUXZ
    stops (default=False) Allow '*' character
    ignore_case (default=False) Check sequences insensitive to case
    min_seqs (default=1) Minimum number of sequences allowed
    max_seqs (default=1) Maximum number of sequences allowed
'''

    # Remove unwanted whitespace
    inp = inp.strip().replace('\r','')

    # Handle empty string
    if not inp:
        if min_seqs:
            raise FastaError('input is an empty string; at least %i sequence(s) required' % min_seqs, 0, '')
        else:
            # should be None?
            return ''

    # Convert to_upper and ignore_case to re.I for use as a re.compile argument
    to_upper = to_upper and re.I
    ignore_case = (ignore_case or to_upper) and re.I

    # Handle missing initial header
    if not inp.startswith('>'):
        if max_seqs == 1:
            inp = '>Submission\n' + inp
        else:
            raise FastaError('missing first header',
                0, inp.split('\n',1)[0])

    # If no changes need to be made to the sequences, validate via regex
    if not (wrap or alignment or to_upper or max_length) and \
        protein_fasta_re(alignment=alignment, allow=allow,
        extra_chars=extra_chars, stops=stops, ignore_case=ignore_case,
        bad_header_chars=bad_header_chars, max_seqs=max_seqs,
        min_seqs=min_seqs).match(inp):
        return (inp + '\n', inp.count('\n>')+1)
    
    # Parse FASTA input, line by line
    recs = []
    seq_char_class = protein_char_class(
        alignment=alignment, allow=allow, extra_chars=extra_chars,
        stops=stops, nl=False)
    seq_line_re = re.compile('^%s*$' % seq_char_class, ignore_case)
    for num,line in enumerate(inp.split('\n')):
        line = line.strip()
        if not line:
            continue

        if seq_line_re.match(line.replace(' ', '')): # Catch a sequence line
            recs[-1][1] += line
        elif line.startswith('>'): # Catch a header line
            
            if max_length is not None and recs and len(recs[-1][1]) > max_length:
                raise FastaError("sequence contains %i residues; only %i residue(s) permitted per sequence" % (
                    len(recs[-1][1]), max_length),
                    num, recs[-1][0])

            # Check for empty header
            if not len(line[1:].strip()):
                raise FastaError("illegal empty header", num, line)

            recs.append(['>%s' % line[1:].strip(), '', num])


            # Check for illegal header characters
            bad_chars = set(line).intersection(bad_header_chars)
            if bad_chars:
                raise FastaError("illegal header character(s) '%s'" % (
                    "', '".join(bad_chars)),
                    num, line)
            # Check for too many sequences
            if max_seqs is not None and len(recs) > max_seqs:
                raise FastaError(
                    'contains over %i sequences; up to %i sequences permitted' \
                    % (len(recs)-1, max_seqs), recs[-1][2], recs[-1][0])
        else: # Error
            raise FastaError("illegal sequence character(s) '%s'" % (
                "', '".join(re.compile(seq_char_class,
                ignore_case).sub('', ''.join(set(line.replace(' ', '')))))),
                num, line)
    
    # Check for too few sequences
    if min_seqs is not None and len(recs) < min_seqs:
        num,line = len(recs) and (recs[-1][2], recs[-1][0]) or (0,'')
        raise FastaError('contains %i sequences; at least %i sequences required' % (len(recs), min_seqs), num, line)

    # Rewrap if necessary
    sizer = wrap and fill or (lambda x: x)

    seqs = zip(*recs)[1]

    # Check alignment lengths, if necessary
    if alignment and not reduce(lambda x,y: x == y and x or False,
        map(lambda x: len(x), seqs)):
        len_aln = len(recs[0][1])
        for rec in recs:
            if len(rec[1]) != len_aln:
                raise FastaError("length of sequence '%s' differs from previous sequence(s) in alignment" % rec[0], rec[2], rec[0])

    # Check for empty sequences
    for rec in recs:
        header,seq,num = rec
        if not seq:
            raise FastaError("header '%s' contains no sequence" % header,
                num, header)
        rec[1] = sizer((to_upper and seq.upper()) or seq)

    return ('\n'.join(map(lambda x: '\n'.join(x[0:2]), recs))+'\n',len(rec))

def format_a2m(inp, max_seqs=100, extra_chars=False):
    '''Wrapper function for format_protein_fasta, set to match a2m.
'''
    return format_protein_fasta(inp, ignore_case=True, max_seqs=max_seqs,
        extra_chars=extra_chars, allow='.', alignment=True)

def format_protein_fasta(*args, **kwargs):
    '''Check FASTA input for format errors and return a string with requested changes.

Parameters:
    to_upper (default=False) Converts sequences to upper case letters
    wrap (default=False) Reset the line wrapping for sequences
    alignment (default=False) Allow '-' character and enforce equal lengths
    allow (default=None) String of additional characters allowed in the sequence
    extra_chars (default=False) Allow BOUXZ
    stops (default=False) Allow '*' character
    ignore_case (default=False) Check sequences insensitive to case
    min_seqs (default=1) Minimum number of sequences allowed
    max_seqs (default=1) Maximum number of sequences allowed
'''
    return format_protein_fasta_and_count(*args, **kwargs)[0]
