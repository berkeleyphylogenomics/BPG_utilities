def getClass(data):
    data = data.split(' ')

    for i in range(1,10):
        for val in data:
            if val == '':
                data.remove(val)

    for val in data:
        if val.find('CLASS')!= -1:
            key = data.index(val)

    key = key + 1

    if data[key] == 'Unclassified;':
        key = key + 1
        u_data = data[key].replace(';', '')
        return u_data

    data = data[key].replace(';', '')
    return data

def splitGenes(data):
    import re
    import string
    
    data = data.split('GENES')
    
    if len(data)<2:
        x = ['null']
        return x
    kegg_gene = data[1].split('\n')
    
    z = []
    for val in kegg_gene:
        val = val.split(':')
        z.append(val)
    
    print '\n'
    
    ids = []
    orgs = []
    for val in z:
        if len(val)>1:
           ids.append(val[1].split(' '))
           orgs.append(val[0])
    
    for val in ids:
        for entry in val:
            if entry == '':
                val.remove(entry)
    
    trimmed_ids = []
    
    for val in ids:
        by_org = []
        for str in val:
            str = re.sub('[(][^)]*[)]', '', str)
            by_org.append(str)
        trimmed_ids.append(by_org)
    
    orgs_trimmed = []
    for val in orgs:
        val = val.replace(' ', '')
        val = string.swapcase(val)
        orgs_trimmed.append(val)
    
    seq_id = []   
    for val in trimmed_ids:
        key = trimmed_ids.index(val)
        for entry in val:
            pp_input = orgs_trimmed[key] + ':' + entry
            seq_id.append(pp_input)
    
    return seq_id       
        
def grabUniprot(organism, seq):
    import ftplib 
    ftp = ftplib.FTP("ftp.genome.jp")
    ftp.login("anonymous", "")
    fp = []
    
    try: 
        ftp.retrlines('RETR ' + "/pub/kegg/genes/organisms/" + organism + "/" + organism + "_uniprot.list", fp.append)
    except:
        return 'null'
    
    ftp.quit()
    for val in fp:
        val = val.split('\tup:')
        if val[0] == seq:
            return val[1]
        
    return 'null'
            

def getKeggnode(data):
    data = data.replace(' ','')
    data = data.split('\n')
    data = data[0].replace('ENTRY', '')
    data = data[:-2]
    return data
    
    

