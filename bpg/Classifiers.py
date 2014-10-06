import os, types
from pfacts003.phylofacts.models import Family, UniProt
from time import localtime, strftime
import tempfile, subprocess, shlex


class HMMBLASTSequenceClassifier:
    tempdir = "/clusterfs/ohana/software/webserver/submitted_jobs/"
 
    class ReturnObject:
        """ This class defines a return object for the hmm-blast sequence classifier.  
            The structure of the return object is:
            familyreturnobject = {
                'jobID'       : '' string - has job been created, if so what dir?
                'errors'   : [] list of strings,
                'defline'  : '' string,
                'sequence' : '' string,
                'status'   : '' string
                'families'   : [{  list of dicts,
                        'id'       : '' string,
                        'family'   : '' string,
                        'fdesc'    : '' string
                        'pfam_dom' : [{ list of dicts
                            'name' : '' string,
                            'desc' : '' string
                }]
            }]
        """
        def __init__(self):
            self.stat='' # default for start
            self.families = []
            self.defline = ''
            self.sequence = ''
            self.errors = []
            self.jobID = ''

        def as_dict(self):
            # returns the returnobject as a dictionary
            returndict = \
            {
                'status'    : '', 
                'families'  : [ 
                    {
                    'id'        : '', 
                    'family'    : '', 
                    'fdesc'     : '', 
                    'pfam_dom'  : [{
                        'name'  : '', 
                        'desc'  : ''
                        }]
                    }],
                'errors' : [], 
                'defline' : '', 
                'sequence' : '', 
                'jobID' : ''
            }
            returndict['status'] = self.stat
            returndict['errors'] = [err for err in self.errors]
            returndict['defline'] = self.defline
            returndict['sequence'] = self.sequence
            returndict['families'] = [fam for fam in self.families]
            returndict['jobID'] = self.jobID
            return returndict

        def set_jobid(self, jid):
            if isinstance(jid, types.StringType):
                self.jobID = jid

        def append_error(self, error):
                self.errors.append(str(error))

        def append_pfam(self, pfam, fam):
            for lfam in self.families:
                if lfam['family'] == fam:
                    if isinstance(pfam, types.DictType):
                        lfam['pfam_dom'].append(pfam)

        def get_pfam_list_for_family(self, family):
            ret = []
            for fam in self.families:
                if fam['family'] == family:
                    for pfam in fam['pfam_dom']:
                        if pfam['name'] not in ret:
                            ret.append(pfam['name'])
            return ret

        def set_defline(self, defline):
                self.defline = str(defline)

        def set_sequence(self, sequence):
                self.sequence = str(sequence)

        def append_family(self, family):
            if isinstance(family, types.DictType):
                self.families.append(family) 
    
    def run_HMM_BLAST_classification(self, args=None):
        os.chdir(self.dir)
        st = open(os.path.join(self.dir, 'status'), 'w')
        st.write("\n{'status':'running','stepsdone':'0','totalsteps':'7'}")
        st.close()
        os.chmod("status", 0777)
        script = open(os.path.join(self.dir, 'hmmblastpipe.sh'), 'w')
        os.chmod("hmmblastpipe.sh", 0777)
        script.write("#!/bin/bash\n#PBS -N HMMBLASTClassification\n")
        script.write("#$ -cwd\n#$ -j y\n#$ -S /bin/bash\n")
        script.write("source /clusterfs/ohana/software/ohana_environment.sh\n")
        script.write("development_old\ndevelopment\ndevelopment_old\n")
        #script.write("HMMBLASTPipeline.py -d %s" % self.dir)
        script.write("ls")
        cmd = "qsub -d %s -j oe -o output -lwalltime=288:00:00 hmmblastpipe.sh" % self.dir
        cmd2 = "ls"
        p = subprocess.Popen(shlex.split(cmd)) 
        return

    def status(self):
        f = open(os.path.join(self.dir, 'status'), 'r')
        lines = f.readlines()
        f.close()
        return eval(lines[-1])

    def get_job_id(self):
        return self.id

    def get_dir(self):
        return self.dir

    @classmethod
    def get_job_by_id(cls, id):
        return cls ( "", "", id )
    
    @classmethod
    def create_job( cls, defline, sequence):
        return cls ( defline, sequence, "" )

    def __init__(self, defline, sequence, ID):
        if not ID:
            self.id = "job6"
            self.dir = os.path.join(self.tempdir, self.id)
            os.mkdir(self.dir)
            os.chmod(self.dir, 0777)
            f = open(os.path.join(self.dir, 'query.fa'), 'w')
            f.write(defline + "\n" + sequence)
            f.close()
            os.chmod(os.path.join(self.dir, 'query.fa'), 0777)
        else:
            self.id = ID
            self.dir = os.path.join(self.tempdir, self.id)
        
