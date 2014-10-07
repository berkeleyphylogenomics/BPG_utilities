#!/clusterfs/ohana/software/bin/python2.7
'''
This file has methods that are useful for scraping files and putting the information
into tuples.

It does not deal with loading stuff into the db, and it doesn't really know anything
about what our db looks like.  Except for the implicit assumption that the data it
pulls out is the important stuff
'''


'''
http://www.geneontology.org/GO.format.gaf-2_0.shtml
'''
def parse_go_annotation_file(infile):
    for line in infile:
        line = line.strip()
        if len(line) == 0 or line[0] == "!":
            continue
        data = line.split('\t')
        date = "%s-%s-%s" % (data[13][0:4], data[13][4:6], data[13][6:8])
        yield (data[1], data[4], data[6], data[14], date)

def parse_orthodb_annotation_file(infile, date):
    for line in infile:
        line = line.strip()
        if len(line) == 0 or line[0] == "!":
            continue
        data = line.split('\t')
        yield (data[0], data[1], "IEA", "OrthoDB5", date)


'''
We support the OBO XML file because the OBO format is just annoying enough that writing
a parser would take too long, and obo-edit can convert any obo to obo-xml, so we'll
let them deal with their own stupid encapsulation format

We are _not_ trying to pull out all information from the file.  Just the information
that we need to fit the OBO into our idea of what an ontology is
'''
def parse_obo_xml(infile):
    from lxml import etree

    context = etree.iterparse(infile, events=('end',) )

    terms = []
    subsets = []
    subset_membership = []
    connections = []
    xrefs = []
    name = ""

    # As of now, we don't return them, we add them to subsets
    namespaces = set()

    for event, elem in context:
        if elem.tag not in ['header', 'term']:
            continue

        if elem.tag == 'header':
            for subsetdef in elem.findall('subsetdef'):
                subsets.append((subsetdef.find('id').text, subsetdef.find('name').text))
            name = elem.find('ontology').text
        elif elem.tag == 'term':
            terms.append((elem.find('id').text, elem.find('name').text))
            for is_a in elem.findall('is_a'):
                connections.append((
                    elem.find('id').text,
                    'is_a',
                    is_a.text,
                ))
            for relationship in elem.findall('relationship'):
                connections.append((
                    elem.find('id').text,
                    relationship.find('type').text,
                    relationship.find('to').text,
                ))
            for s in elem.findall('subset'):
                subset_membership.append(( elem.find('id').text, s.text ))
            # Namespaces are just a type of subset.
            for s in elem.findall('namespace'):
                subset_membership.append(( elem.find('id').text, s.text ))
                namespaces.add(s.text)
            for xref in elem.findall('xref_analog'):
                xrefs.append( ( elem.find('id').text, xref.find('dbname').text, xref.find('acc').text ) )
        elem.clear()

    for ns in namespaces:
        subsets.append((ns, ''))
    return { 'subsets': subsets, 'terms': terms, 'subset_membership': subset_membership,
             'connections': connections, 'xrefs': xrefs, 'name': name }


if __name__ == '__main__':
    pass
