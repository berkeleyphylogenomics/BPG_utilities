import functions
import dbconnect
import string

#get data(ko file)
f = open('ko')
data = f.read()
data = data.split('///')

#TEST
#data = []
#data.append(f.read())

print 'Updating' 

for val in data:
    #input array will be in form [keggnode, class, sequence_id]
    db_input = []

    #get keggnode ID
    db_input.append(functions.getKeggnode(val))
    
    #get gene class
    db_input.append(functions.getClass(val))

    #get sequence id
    rows = functions.splitGenes(val)
    
    for id in rows:
        org = id.split(':')[0]
        uid = functions.grabUniprot(org, id)
        if uid != 'null':
            dbconnect.cur.execute("INSERT INTO keggnode_keggrefseq(keggnode, type, seq_id) VALUES('%s', '%s', '%s')" % (db_input[0], db_input[1], uid))
    
    print '...'

print 'DONE'


    
    
    





      