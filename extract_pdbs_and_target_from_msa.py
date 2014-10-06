# Input :  An MSA file.  For example, bpg071753.a2m
# Input :  A list of PDB IDs that you want to grab out of your MSA.
#          For example, significant_pdb_against_T0387_ghmm_hits.id
# Output:  An MSA file only containing the sequences you specified in
#          your PDB IDs file.

import string

def main():
    # The two input files
    #msa_input_name = "/home/tfarrah/CASP8/targets/T0387/bpg071753.a2m"
    msa_input_name = "bpg071753.a2m"
    msa_input_handle = open(msa_input_name, 'r')
    #pdb_ids_name = "/home/tfarrah/CASP8/targets/T0387/significant_pdb_against_T0387_ghmm_hits.id"
    pdb_ids_name = "significant_pdb_against_T0387_ghmm_hits.id"
    pdb_ids_handle = open(pdb_ids_name, 'r')

    # Make an output file
    #msa_output_name = "/home/tfarrah/CASP8/targets/T0387/bpg071753_significant.a2m"
    msa_output_name = "bpg071753_significant.a2m"
    msa_output_handle = open(msa_output_name, 'w')

    # Also need the translation table from PDB ID to UniprotID
    #idmap_name = "/usr/local/blastdb/idmapping.tb"
    idmap_name = "idmapping.tb"
    idmap_handle = open(idmap_name, 'r')

    # Make the mapping table from PDB ID to UniprotID
    # Reference "/usr/local/blastdb/idmapping.tb.readme" for columns of idmapping.tb
    print "Making the mapping from PDB IDs to Uniprot IDs..."
    idmap_entry = idmap_handle.readline()
    idmap_dict = {}
    while (idmap_entry != ""):
        idmap_entry_split = idmap_entry.split("\t")
        uniprotID = idmap_entry_split[0]
        pdbID = idmap_entry_split[5]
        if (pdbID != ""):
            # Turn the colon into an underscore
            pdbID = pdbID.replace(":", "_").upper()
            idmap_dict[pdbID] = uniprotID
        idmap_entry = idmap_handle.readline()
    idmap_handle.close()
    print "Done making the mapping."

    # Turn your input MSA into a dict so you can grab stuff from it easily
    # Concatenate the first two lines together of your MSA input
    print "Reading in the input MSA..."
    msa_input_entry = msa_input_handle.readline()+msa_input_handle.readline()
    msa_dict = {}
    while (msa_input_entry != ""):
        # If you're reading the target (it doesn't have a pipe), then
        # immediately write this entry to the output MSA
        if (msa_input_entry.find("|") == -1):
            msa_output_handle.write(msa_input_entry)
        else:
            msa_input_entry_split = msa_input_entry.split("|")
            # Grab the 2nd entry in the split array, because it's the UniprotID
            uniprotID = msa_input_entry_split[1]
            msa_dict[uniprotID] = msa_input_entry
        msa_input_entry = msa_input_handle.readline()+msa_input_handle.readline()
    msa_input_handle.close()
    print "Done reading."

    # Loop through the PDB ID file and grab all the stuff you want
    print "Looping through the PDB IDs for corresponding entries in the MSA..."
    print "Writing results to output..."
    pdbID = pdb_ids_handle.readline().strip().upper()
    num_mapped, num_unmapped = 0, 0
    while (pdbID != ""):
        if pdbID in idmap_dict:
            uniprotID = idmap_dict[pdbID]
            if uniprotID in msa_dict:
                msa_output_entry = msa_dict[uniprotID]
                msa_output_handle.write(msa_output_entry)
                num_mapped += 1
            else:
                num_unmapped += 1
        else:
            num_unmapped += 1
        pdbID = pdb_ids_handle.readline().strip().upper()
    pdb_ids_handle.close()
    msa_output_handle.close()

    print "Done.\n"
    print "Summary statistics:"
    print "There were "+str(num_mapped)+" PDB sequences found in the input MSA that were moved to the output."
    print "There were "+str(num_unmapped)+" PDB sequences not found in the input MSA."
    print "The output file is located at: "+msa_output_name

if __name__ == "__main__":
    main()
