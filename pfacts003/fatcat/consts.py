''' File for fatcat constants that aren't user settable, but we still want some kind of control over them.'''

CONSENSUS_UNIPROT_BASE_PARAMETER=0.66
CONSENSUS_UNIPROT_THRESHOLD_FOR_HIGH_CONFIDENCE=0.7
CONSENSUS_UNIPROT_THRESHOLD_FOR_MEDIUM_CONFIDENCE=0.4

STAGE1_MINIMUM_ALIGNMENT_LENGTH=50

NUM_FATCAT_STAGES = 12

HIGH_RECALL_PARAMS = {
    'stage1-eval': 0.001,
    'stage1-mda-family-qcov': 0.6,
    'stage1-mda-family-hcov': 0.6, 
    'stage1-pfam-family-hcov': 0.6,
    'stage2-eval': 0.00001,
    'stage2-mda-family-qcov': 0.6,
    'stage2-mda-family-hcov': 0.6,
    'stage2-pfam-family-hcov': 0.6,
    'query-hmm-pid': 50.0,
    'ec-checkboxes': 'phog kerf oma',
    'ec-kerf-threshold': 50,
    'ec-orthology-number': 1,
    'cluster-similarity': 97,
    'stage2-maxseqs': 200,
    'ortholog-coverage': 70,
    'query-coverage': 70,
    'minimum-pid-for-orthology': 50,
    'stage3-koa': 100,
    'lambda': 0.8,
    'amplitude': 1.0,
    'consensus-uniprot-threshold-high': 70.0,
    'consensus-uniprot-threshold-medium': 40.0
} 
    
HIGH_PRECISION_PARAMS = {
    'stage1-eval': 0.0001,
    'stage1-mda-family-qcov': 0.7,
    'stage1-mda-family-hcov': 0.7, 
    'stage1-pfam-family-hcov': 0.7,
    'stage2-eval': 0.000001,
    'stage2-mda-family-qcov': 0.7,
    'stage2-mda-family-hcov': 0.7,
    'stage2-pfam-family-hcov': 0.7,
    'query-hmm-pid': 60.0,
    'ec-checkboxes': 'phog kerf oma',
    'ec-kerf-threshold': 60,
    'ec-orthology-number': 1,
    'cluster-similarity': 97,
    'stage2-maxseqs': 200,
    'ortholog-coverage': 70,
    'query-coverage': 70,
    'minimum-pid-for-orthology': 73,
    'stage3-koa': 100,
    'lambda': 0.8,
    'amplitude': 1.0,
    'consensus-uniprot-threshold-high': 70.0,
    'consensus-uniprot-threshold-medium': 40.0
} 

REMOTE_HOMOLOG_PARAMS = {
    'stage1-eval': 0.1,
    'stage1-mda-family-qcov': 0.5,
    'stage1-mda-family-hcov': 0.5, 
    'stage1-pfam-family-hcov': 0.5,
    'stage2-eval': 0.001,
    'stage2-mda-family-qcov': 0.5,
    'stage2-mda-family-hcov': 0.5,
    'stage2-pfam-family-hcov': 0.5,
    'query-hmm-pid': 20.0,
    'ec-checkboxes': 'phog kerf oma',
    'ec-kerf-threshold': 20,
    'ec-orthology-number': 1,
    'cluster-similarity': 97,
    'stage2-maxseqs': 200,
    'ortholog-coverage': 30,
    'query-coverage': 30,
    'minimum-pid-for-orthology': 17,
    'stage3-koa': 100,
    'lambda': 0.8,
    'amplitude': 1.0,
    'consensus-uniprot-threshold-high': 70.0,
    'consensus-uniprot-threshold-medium': 40.0
} 

FRAGMENT_PARAMS = {
    'stage1-eval': 1,
    'stage1-mda-family-qcov': 0.0,
    'stage1-mda-family-hcov': 0.0, 
    'stage1-pfam-family-hcov': 0.0,
    'stage2-eval': 1.0,
    'stage2-mda-family-qcov': 0.0,
    'stage2-mda-family-hcov': 0.0,
    'stage2-pfam-family-hcov': 0.0,
    'query-hmm-pid': 30.0,
    'ec-checkboxes': 'phog kerf oma',
    'ec-kerf-threshold': 40,
    'ec-orthology-number': 1,
    'cluster-similarity': 97,
    'stage2-maxseqs': 200,
    'ortholog-coverage': 0,
    'query-coverage': 70,
    'minimum-pid-for-orthology': 50,
    'stage3-koa': 100,
    'lambda': 0.8,
    'amplitude': 1.0,
    'consensus-uniprot-threshold-high': 70.0,
    'consensus-uniprot-threshold-medium': 40.0
}
# numbers of threads for multithreaded stuff
THREADS = 8

# for stage 4 clustering of orthologs/other sequence matches
STAGE4_CLUSTERING_THRESHOLD = 5
