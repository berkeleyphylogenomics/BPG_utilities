#!/usr/bin/python

###############################################
#
# This script is intended to be a module that
# will check user input to make sure it is 
# appropriate for our servers.
#
# The appropriate input is a single protein
# sequence in fasta format, or a single
# protein sequence with no fasta header
#
###############################################

import sys
import re

###############################################
#
# The inputCheck function takes a single string
# as input and returns a tuple (bool, string).
#
# The first element of the tuple is True if the
# input is valid and False otherwise. The second
# element is a 'validized' user input string (if
# the first element is True) or an error message
# (if the first element is False)
#
# So, you should be able to do something like:
#
#    user_input = getUserInput()
#    (valid, result) = inputCheck(user_input)
#    if valid:
#      processInput(result)
#    else:
#      processError(result)
#
def inputCheck(userinput):
    
    # split user input into lines
    userinput = userinput.replace('\r', '\n')
    userinput = userinput.replace('\n\n', '\n')
    
    lines = userinput.split('\n')

    # remove blank lines
    processed_lines = []

    for line in lines:
        line = line.strip()

        if len(line) > 0:
            processed_lines.append(line)

    #--- error check, is input empty? ---#
    if len(processed_lines) < 1:
        return (False, 'Input is empty, I need a protein sequence')

    # if first line starts with '>', it's a header
    header = lines[0]
    seqstart = 1

    if header[0] == '>':
        # remove any spaces after the '>'
        header = header[1:len(header)]
        header = header.strip()

        # header has characters, so it's okay
        if len(header) > 0:
            header = '>' + header

        # header is empty, replace it
        else:
            header = '>query'

    # if first line contains non-sequence characters, it's a header
    elif not re.match('^[a-zA-Z]+$', header):
        header = '>' + header

    # first line is sequence, create our own header
    else:
        header = '>query'
        seqstart = 0
        
    # parse the sequence data, which should be everything after
    # the header
    sequence = ''
    
    for i in range(seqstart, len(lines)):
        sequence += lines[i].upper() # convert to upper case only

    # replace spaces
    sequence = sequence.replace(' ', '')
    sequence = sequence.replace('\t', '')

    # check that the sequence is valid
    
    #--- error check - all valid characters? ---#
    alphapatt = re.compile('[A-Z]')

    for i in range(len(sequence)):
        if alphapatt.match(sequence[i]) and sequence[i] != 'O':
            pass
        else:
            return (False, 'Illegal residue: \'%s\' at position %d in input sequence' % (sequence[i], i+1))

    #--- error check - is this NT sequence? ---#
    acgtu = 0

    for i in range(len(sequence)):
        if sequence[i]=='A' or sequence[i]=='C' or sequence[i]=='G' or sequence[i]=='T' or sequence[i]=='U':
            acgtu += 1

    acgtuprop = float(acgtu) / float(len(sequence))

    if acgtuprop > 0.95:
        return (False, 'Input appears to be nucleotide sequence, I need a protein sequence')

    # input appears to be valid, accept it #
    myresult = header + '\n' + sequence
    return (True, myresult)
# end inputCheck


#############################################
#
# Tests a single user input sequence, which
# is in a text file, to see if it is valid
# Prints a 'validized' version of the input
#
def main():

    if len(sys.argv) < 2:
        print 'Usage: checkInput.py filename'
        print '  reads filename and checks if it has a valid sequence'
        return 1

    userinput = ''

    handle = open(sys.argv[1])
    
    for line in handle:
        userinput += line

    handle.close()

    print 'checking input:'
    print userinput
    print ''

    (valid, result) = inputCheck(userinput)

    print 'valid: %s' % valid
    print 'result:'
    print result

    return 0
# end main



if __name__ == '__main__':
    main()

