''' This file conains functions that populate phylo4j tables. '''
from pfacts003.phylo4j.consts import *
from collections import Counter
from py2neo import neo4j

def explorer_node_property_table(node):
    # this function returns the pure html for the node property table 
    # for the phylo4j explorer
    r = '<table><thead><tr><th>Property</th><th>Value</th></tr></thead><tbody>'
    for (p, v) in node.get_properties().items():
        if not isinstance(v, (basestring, int, float)):
            # this is a list or dict, assume list for now
            r += '<tr><td>%s</td><td>' % p
            r += ', '.join([str(s) for s in v])
            r += '</td></tr>'
        else:
            if p == 'sequence':
                r += '<tr><td>%s</td><td class="sequence-text">' % p
            else:
                r += '<tr><td>%s</td><td>' % p
            if isinstance(v, basestring) and len(v) > MAX_TABLE_CHARS:
                # Need to split this row
                for (index, ch) in enumerate(v):
                    if index % MAX_TABLE_CHARS == 0 and not index == 0:
                        r += '<br />%s' % ch
                    else:
                        r += ch
                r += '</td></tr>'
            else:
                r += '%s</td></tr>' % (str(v))
    r += '</tbody></table>'
    return r

def explorer_relationships_table(node):
    # returns the type and number of relatinoships attached to this node
    # again, as the other functions in this file, it returns a json object with
    # the pure html table, and the list of active edges on this node
    c = Counter()
    for rel in node.match():
        c[rel.type] += 1
    r = '<table><thead><tr><th>Relationship</th><th>Number</th></tr></thead><tbody>'
    e = []
    for (p, v) in c.items():
        e.append(p)
        r += '<tr><td>%s</td><td>%s</td></tr>' % (p, str(v))
    r += '</tbody></table>'
    return {'html': r, 'relationships': e}
    

