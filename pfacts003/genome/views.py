'''
View that populates templates for the genome page
'''
from django.shortcuts import render_to_response

def genome_page(request, genome_name):
    '''Returns a response with dictionary read from the config file.

    The config file (for E. coli in this example
    should reside at /clusterfs/ohana/bpg/pfacts/genomes/e_coli/config.txt
    '''
    #Create dictionary and render, else render 404
    genome_dict = create_dictionary_for_genome(genome_name)
    # Add an optional variable for the short name of the species.
    # Short name is just the short genome name (without the underscore).
    temp_short_name = genome_name.split('_')
    temp_short_name[0] = temp_short_name[0].capitalize() + '.'
    genome_dict['short_genome_name'] = ' '.join(temp_short_name)
    if genome_dict == {}:
        return render_to_response('404.html')
    return render_to_response('phylofacts/genomes/main.html',
                              {'genome_dict' : genome_dict})
    

def create_dictionary_for_genome(genome_name):
    '''Creates dictionary for specific genome.'''
    #Check if file is present, and get contents.
    #Else return empty dictionary
    try:
        genome_config = [line.strip().split('\t') for line in
                         open('/clusterfs/ohana/bpg/pfacts/genomes/%s/config.txt'
                               % genome_name).readlines()
                         if not (line.startswith('#') or
                                 line.strip() == '')]
    except IOError:
        return {}
    return dict(genome_config)
    
