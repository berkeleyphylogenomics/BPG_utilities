#!/clusterfs/ohana/software/bin/python2.7
#PBS -z
import sys, time, re, os, smtplib, subprocess
sys.path.insert(0,"/home/cyrus_afrasiabi/ohana_repository")
sys.path.insert(0,"/clusterfs/ohana/software/lib/python2.7/site-packages/")
os.environ["DJANGO_SETTINGS_MODULE"] = "pfacts003.settings"
from pfacts003.phylofacts.models import FxnSitePredictionJob
from pfacts003.utils.hmm import gather_homologs_with_jackhmmer, hmmpress
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from email.MIMEText import MIMEText
from bpg.common.utils.sequenceutils.fasta_utils import get_uniprot_accession_from_header

EMAIL_FROM = "phylofacts.webmaster@gmail.com"
EMAIL_PASSWORD = "lanikai324c"

def get_seed_index(msa, header):
    """
    This function will return the seed sequences row in the MSA as an integer index.
    """
    for (index, seq_record) in enumerate(msa):
        if (seq_record.id in header):
            return index
    return 0


def remove_seed_gaps_and_inserts(msa, seed_index):
    """
    This function identifies all the columns in the msa that are gaps in the seed sequence, and
    also removes any hmm insert states (. or lowercase).  This should return an alignment with 
    length equal to the length of the input query.
    """
    # This is probably what people around here call hacky.  I can see it.
    first = True
    for (index, char) in enumerate(str(msa[seed_index].seq)):
        if (char.isupper()):
            if first:
                return_msa = msa[:,index:index+1]
                first = False
            else:
                return_msa += msa[:,index:index+1]
    return return_msa
 

def remove_seed_duplicates(msa, seed_index):
    sequence = str(msa[seed_index].seq)
    return_msa = MultipleSeqAlignment([])
    for (index, seq_record) in enumerate(msa):
        if (index == seed_index):
            return_msa.append(seq_record)
        else:
            if (str(seq_record.seq) != sequence):
                return_msa.append(seq_record)
    return return_msa
                   
job_id = int(os.environ["job_id"])

job = FxnSitePredictionJob.objects.get(id=job_id, status__id=2)

if not job:
    raise

job.status_id = 3
job.save()

job_directory = os.path.join("/clusterfs/ohana/software/webserver/temp/", "%d" % job_id)

os.mkdir(job_directory)
os.chdir(job_directory)

# load the job parameters
user_email = job.user_email
fasta = str(job.fasta)
treecut_pid = float(job.treecut_pid)
pbs_id = int(job.pbs_job_id)
jackhmmer_iterations = int(job.jackhmmer_iterations)
jackhmmer_evalue = float(job.jackhmmer_evalue)

# create the output file names
full_msa_path = os.path.join(job_directory, "job_%d.msa" % job_id)
full_tree_path = os.path.join(job_directory, "job_%d.tree" % job_id)
summary_msa_path = os.path.join(job_directory, "job_%d_summary.msa" % job_id)
summary_tree_path = os.path.join(job_directory, "job_%d_summary.tree" % job_id)
hmm_path = os.path.join(job_directory, "job_%d.hmm" % job_id)
intrepid_config_path = os.path.join(job_directory, "config.txt")
intrepid_out_path = os.path.join(job_directory, "intrepid.out")

job.status_id = 4
job.save()

# Gather homologs
all_msa = gather_homologs_with_jackhmmer(fasta, E=jackhmmer_evalue, N=jackhmmer_iterations)
jackhmmer_msa = all_msa[0]
del all_msa
job.save()

job.status_id = 5
job.save()

# Process the MSA
seed_index = get_seed_index(jackhmmer_msa, job.fasta_header)
jackhmmer_msa = remove_seed_duplicates(jackhmmer_msa, seed_index)
jackhmmer_msa = remove_seed_gaps_and_inserts(jackhmmer_msa, seed_index)
job.num_homologs = len(jackhmmer_msa)
seed_acc = get_uniprot_accession_from_header(jackhmmer_msa[seed_index].id)
# Now write the MSA
full_msa_file = open(full_msa_path, "w")
AlignIO.write(jackhmmer_msa, full_msa_file, "fasta")
full_msa_file.close()
del jackhmmer_msa

job.status_id = 6
job.save()

# Generate Tree
args = ["/clusterfs/ohana/software/bin/FastTree", "-quiet", "-nopr", full_msa_path]
p = subprocess.Popen(args, shell=False, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
tree = p.communicate()[0]
full_tree_file = open(full_tree_path, "w")
full_tree_file.write(tree)
full_tree_file.close()

job.status_id = 7
job.save()

# Generate HMM
args = ["/clusterfs/ohana/software/bin/hmm3build", "-o", "/dev/null", "--informat", 
        "afa", "--amino", "-n", "PhyloFacts_FSP_Job_%d_HMM" % job_id, hmm_path, full_msa_path] 
p = subprocess.Popen(args, shell=False, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
j = p.communicate()[0]
hmmpress(hmm_path)

job.status_id = 8
job.save()

#Generating Summary Tree and MSA
args = [""]

job.status_id = 9
job.save()

# Running INTREPID
os.chdir(job_directory)
intrepid_config_file = open(intrepid_config_path, "w")
intrepid_config_file.write("msa_file %s\ntree_file %s\nsequence_id %s\n" % (full_msa_path, full_tree_path, seed_acc))
intrepid_config_file.close()
args = ["/clusterfs/ohana/software/test/repo_bin/intrepid.pl", intrepid_config_path]
p = subprocess.Popen(args, shell=False, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
intrepid_out = p.communicate()[0]
intrepid_out_file = open(intrepid_out_path, "w")
intrepid_out_file.write(intrepid_out)
intrepid_out_file.close()

job.status_id = 10
job.save()

# Create email to notify user job is done.
"""
email_body = 'PhyloFacts Functional Site Prediction Job %d has completed.\n\nClick here to view your results page for Job for query %s\n%s\n\n' % (job_id, job.fasta_header, job.fasta_sequence[:10])
msg = MIMEText(email_body)
msg['Subject'] = "PhyloFacts Functional Site Prediction Job %d is complete." % job_id
msg['From'] = EMAIL_FROM
msg['To'] = job.user_email

# Sending email
server = smtplib.SMTP('smtp.gmail.com:587')
server.ehlo()
server.starttls()
server.ehlo()
server.login(EMAIL_FROM, EMAIL_PASSWORD)
server.sendmail(EMAIL_FROM, job.user_email, msg.as_string())
server.close()
"""
job.status_id = 11
job.save()
