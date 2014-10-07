##################################################################
# TIMEOUT INTERVALS
##################################################################
PROGRESS_PAGE_TIMEOUT_INTERVAL = 3000
MAIN_PAGE_TIMEOUT_INTERVAL = 10000
##################################################################
# STATE MACHINE STUFF
##################################################################
SCAN_FOR_DOMAINS_STAGES = 5
PHYLOBUILDER_DOMAIN_STAGES = 0
INTREPID_DOMAIN_STAGES = 0
DISCERN_DOMAIN_STAGES = 0
##################################################################
# THROTTLING
##################################################################
ANONYMOUS_JOB_LIMIT = 10
AUTHENTICATED_USER_THROTTLE_MESSAGE = "You have submitted %d PhyloFacts Protein Analysis jobs since midnight PST.  You must wait until tomorrow before submitting again.  Please email phylofacts.webmaster@gmail.com to increase your daily job submission quota."
ANONYMOUS_USER_THROTTLE_MESSAGE = 'We have received %d PhyloFacts Protein Analysis submissions from this IP address in the last 24 hours.  This is the maximum allowed for non-registered users.  To become a registered user and increase your quota, please click <a href="/phylofacts/user_create/" target="_blank">here</a>.  Otherwise, you can submit new jobs after midnight PST.'
##################################################################
# FASTA VALIDATION
##################################################################
ACCEPTED_SEQUENCE_CHARS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
NO_FASTA_MESSAGE = 'FASTA input required'
##################################################################
# SEQUENCE LENGTH VALIDATION
##################################################################
MINIMUM_INPUT_LENGTH = 15
MAXIMUM_INPUT_LENGTH = 2000
INPUT_SEQUENCE_TOO_SHORT_MESSAGE = 'Sequence too short.  The PhyloFacts Protein Analysis Webserver does not accept inputs &lt;%d amino acids' % MINIMUM_INPUT_LENGTH
INPUT_SEQUENCE_TOO_LONG_MESSAGE = 'Sequence too long.  The PhyloFacts Protein Analysis Webserver does not accept input &gt;%d amino acids' % MAXIMUM_INPUT_LENGTH
##################################################################
# EMAIL VALIDATION
##################################################################
EMAIL_INVALID_MESSAGE = 'Invalid email'
##################################################################
# JOB COMMENTS VALIDATION
##################################################################
MAXIMUM_COMMENTS_CHARACTERS = 200
JOB_COMMENTS_TOO_LONG_MESSAGE = 'Comments must be &lt;%d characters' % MAXIMUM_COMMENTS_CHARACTERS
##################################################################
# EMAIL SUBJECT VALIDATION
##################################################################
MAXIMUM_EMAIL_SUBJECT_CHARACTERS = 100
EMAIL_SUBJECT_TOO_LONG_MESSAGE = 'Email subject must be &lt;%d characters' % MAXIMUM_EMAIL_SUBJECT_CHARACTERS
##################################################################
# MINIMUM DOMAIN LENGTH VALIDATION
##################################################################
MINIMUM_MINIMUM_DOMAIN_LENGTH = 15
MAXIMUM_MINIMUM_DOMAIN_LENGTH = 200
MINIMUM_DOMAIN_LENGTH_TOO_LONG_MESSAGE = 'Minimum domain length cannot be &lt;%d' % MINIMUM_MINIMUM_DOMAIN_LENGTH
MINIMUM_DOMAIN_LENGTH_TOO_SHORT_MESSAGE = 'Mininmum domain length cannot be &gt;%d' % MAXIMUM_MINIMUM_DOMAIN_LENGTH
MINIMUM_DOMAIN_LENGTH_UNRECOGNIZED_MESSAGE = 'Minimum domain length parameter must be castable to an integer'
MINIMUM_DOMAIN_LENGTH_MISSING = 'Minimum domain length not recognized'
