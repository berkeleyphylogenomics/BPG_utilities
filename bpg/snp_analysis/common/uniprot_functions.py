import urllib,urllib2
import hashlib, base64


def get_acc(name):

    request = 'http://www.uniprot.org/uniprot/'+name+'.list'
    try:
        response = urllib2.urlopen(request)
        result = response.read(20)
    except:
        result = None
    return result


def compute_sequid(seq):
    hashobj = hashlib.sha1()
    hashobj.update(seq)
    seqhash = base64.b64encode(hashobj.digest()).strip('=')
    return seqhash    
