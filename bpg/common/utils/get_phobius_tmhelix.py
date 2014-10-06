#!/usr/bin/env python

# This file contains the code to interact with the phobius webserver.

import mechanize, sys, os, string, glob
from bpg.common.utils.dir_of_family import get_dir_of_family_accession

def get_phobius_results_for_family(acc):

# This function submits a sequence to the phobius webserver, and parses the
# response page. 
#
# INPUT -> family accession.
# OUTPUT -> Returns 3 lists each the same length: (regtype, regstart, regend)  
#
# regtype - list of strings.  each list entry is the regoin type for a phobius 
# predicted region.
# regstart - list of ints.  each list entry is the region start index for a phobius
# predicted region.
# regend - list of ints.   each list entry is the region end index for a phobius
# predicted region.
#
# USAGE:
#
# >>> from bpg.common.utils import get_phobius_tmhelix
# >>> (type, start, end) = get_phobius_tmhelix.get_phobius_results_for_family('bpg0207937')
# >>> type
# ['TOPO_DOM:NON-CYTOPLASMIC.', 'TRANSMEM', 'TOPO_DOM:CYTOPLASMIC.', 'TRANSMEM', 'TOPO_DOM:NON-CYTOPLASMIC.', 'TRANSMEM', 'TOPO_DOM:CYTOPLASMIC.', 'TRANSMEM', 'TOPO_DOM:NON-CYTOPLASMIC.', 'TRANSMEM', 'TOPO_DOM:CYTOPLASMIC.', 'TRANSMEM', 'TOPO_DOM:NON-CYTOPLASMIC.', 'TRANSMEM', 'TOPO_DOM:CYTOPLASMIC.']
# >>> start
# [1, 59, 80, 91, 117, 136, 157, 177, 200, 219, 248, 277, 299, 318, 337]
# >>> end
# [58, 79, 90, 116, 135, 156, 176, 199, 218, 247, 276, 298, 317, 336, 368]
    """
    We had a problem with the phobius server.  Now, instead of doing a web request, we will
    just take whatever info we have in the family directory .phobius file.
    home_response = mechanize.urlopen("http://phobius.sbc.su.se/")
    home_forms = mechanize.ParseResponse(home_response, backwards_compat=False)
    home_forms[0].set_all_readonly(False)
    home_forms[0].set_value(["nog"], name="format")
    home_forms[0].set_value(seq,"protseq")
    result_request = home_forms[0].click(type="submit") 
    result_response = mechanize.urlopen(result_request)
    result_output = result_response.read()
    """
    
    # Get .phobius file from the directory
    directory = get_dir_of_family_accession(acc)
    phobius_file = glob.glob(os.path.join(directory,'*.phobius'))
    if phobius_file:
        read_buffer = open(phobius_file[0],'r')
        result_output = read_buffer.read()
    else:
        return ([],[],[])  
    # Now we have an HTML response containing our TM helix info.  Parse that.

    head = result_output[string.find(result_output, "\n", \
       result_output.find("<h2>Phobius prediction</h2>")):len(result_output)]
    raw_table = head.split("//")[0]
    records = raw_table.split("\n")
    if (len(records) < 3):
        return ([], [], [])
    records = records[1:(len(records)-1)] # remove the first and last empty entries

    # Now we have a list of records, iterate through it and pull out the region type,
    # region start and region end values and store these in lists to return.

    regtype = []
    regstart = []
    regend = []
    for record in records:
        temprec = record.split()
        regstart.append(int(temprec[2]))
        regend.append(int(temprec[3]))
        if len(temprec) == 4:
            regtype.append(temprec[1])
        elif len(temprec) == 5:
            regtype.append(temprec[1]+":"+temprec[4])
        else:
            regtype.append(temprec[1]+":"+temprec[4]+"-"+temprec[5])

    return (regtype, regstart, regend)

def get_transmembrane_regions(regtype, regstart, regend):
# This function takes the three outputs of get_phobius_results_for_seq above, and
# returns the start and end indices of only the transmembrane regions in 2 lists,
# tstart, and tend.
#
# INPUT -> regtype, as described in the comments for get_phobius_results_for_seq.
#          regstart, as described in the comments for get_phobius_results_for_seq,
#          regend, as described in the comments for get_phobius_results_for_seq.
#
# OUTPUT -> tstart -> list of integers corresponding to the start indices of all
#           transmembrane regions.
#           tend -> list of integers corresponding to the end indices of all
#           transmembrane regions.
    tstart = []
    tend = []
    for x in range(len(regtype)):
        if regtype[x] == "TRANSMEM":
            tstart.append(regstart[x])
            tend.append(regend[x])

    return (tstart, tend)

def get_signal_peptide_regions(regtype, regstart, regend):
# This function takes the three outputs of get_phobius_results_for_seq and
# returns .
    tstart = []
    tend = []
    nstart = []
    nend =[]
    hstart = []
    hend =[]
    cstart =[]
    cend = []
    for x in range(len(regtype)):
        if regtype[x] == "SIGNAL":
            tstart.append(regstart[x])
            tend.append(regend[x])
        elif regtype[x] == "REGION:N-REGION.":
            nstart.append(regstart[x])
            nend.append(regend[x])
        elif regtype[x] == "REGION:H-REGION.":
            hstart.append(regstart[x])
            hend.append(regend[x])
        elif regtype[x] == "REGION:C-REGION.":
            cstart.append(regstart[x])
            cend.append(regend[x])
    
    return (tstart, tend, nstart, nend, hstart, hend, cstart, cend)

