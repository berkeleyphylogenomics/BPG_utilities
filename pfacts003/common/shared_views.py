from pfacts003.phylofacts.models import *
from pfacts003.common.shared_forms import SequenceInputForm, PDBInputForm
from pfacts003.flowerpower.forms import FlowerPowerForm

from django.http import HttpResponse
from django.shortcuts import render_to_response

import os
import tempfile
import re

# USE django.settings sites instead -- don't put URL_BASE in settings.py file
URL_BASE = 'http://makana.berkeley.edu' # 

work_prefix = '/clusterfs/ohana/software/webserver/temp/'
job_prefix  = '/clusterfs/ohana/software/webserver/submitted_jobs/'

#################################################
#
# handles displaying the initial sequence input 
# form, launching new analysis jobs, and displaying
# the progress of running jobs. Takes a program_name
# argument, which defines which pipeline to run
#
    


def doIndex(request, program_name, command):
    
    FormType = {
        'discern'       : PDBInputForm,
        'flowerpower'   : FlowerPowerForm,
        'intrepid'      : SequenceInputForm,
        'phylobuilder'  : SequenceInputForm,
    }[program_name]
    
    #if request.REQUEST.__contains__('sequence') and request.REQUEST.__contains__('useremail') and request.REQUEST.__contains__('workdir'):
    if request.REQUEST.items():

        form = FormType(request.REQUEST)

        if form.is_valid():
            theseq      = form.cleaned_data['sequence'] + '\n'
            theemail    = form.cleaned_data['useremail']
            theworkdir  = form.cleaned_data['workdir']

            # no working directory, must be a new job  #
            if not theworkdir:
                # create a temporary directory #
                theworkdir = tempfile.mkdtemp( \
                    prefix='/clusterfs/ohana/software/webserver/temp/')
    
                # write seed sequence to file in mytempdir #
                handle = open(theworkdir + '/myseed.fa', 'w')
                handle.write(theseq)
                handle.close()
                os.system('chmod 777 -R %s' % theworkdir)

                # write uploaded file (if it's there)
                for myfile in request.FILES:
                    handle = open(theworkdir + '/MYPDB' + myfile.name, 'w')

                    for chunk in myfile.chunks():
                        print >>handle, chunk

                    handle.close()

                # write qsub script to execute the given pipeline #
                # runjobsd (daemon) will handle running the job   #

                tmpname = theworkdir.rsplit('/',1)[1]

                handle = open('/clusterfs/ohana/software/webserver/submitted_jobs/%s.sh' % tmpname, 'w')

                handle.write( \
'''%s myseed.fa > %s.screen.log
''' % (command, program_name))

                handle.close()
                os.system('chmod 777 %s' % tmpname)

                # render initial submitted job #
                return render_to_response('%s/exec.html' % program_name,
                                          {'sequence': theseq,
                                           'useremail': theemail,
                                           'workdir': theworkdir})


            # working directory given, must be a running job.
            # check if results are finished and pass log, etc
            # to the html form
            #
            # if results are ready, send an email (if email
            # provided)
            else:


                # read log #
                thelog = ''

                if os.path.exists(theworkdir + '/%s.screen.log' % program_name):
                    handle = open(theworkdir + '/%s.screen.log' % program_name, 'r')
                    thelog = handle.read()
                    handle.close()

                # check if program has failed due to some problem

                theerr = ''
                errorpatt = re.compile('FAILED')

                if errorpatt.search(thelog):
                    theerr = 'Results not available: some error has occurred (see below).'

                # check if program has failed due to too few
                # homologs

                errorpatt = re.compile('Sorry')

                if errorpatt.search(thelog):
                    theerr = 'Results not available: Too few homologs found to run analysis (see below).'

                # check if program has finished (log should
                # say "finished.") If it has, look up 
                # book id so that user can look up the results

                itsdone = ''
                finishedpatt = re.compile('finished.')
                #
                if finishedpatt.search(thelog):
                #if thelog.rstrip().endswith('finished.')
                    # check if family id is available #
                    familypatt = re.compile('family: (\d+)')

                    match = familypatt.search(thelog)

                    if match:
                        itsdone = match.group(1)

                    else: theerr = 'Unknown error.'

                rank = ''
                list_dir = [item for item in os.listdir(theworkdir) if item.endswith('.discern')]
                
                if len(list_dir) == 1:
                    handle = open(theworkdir + '/' + list_dir[0])
                    rank = handle.read()
                    handle.close()

                # render results with intrepid log #
                return render_to_response('%s/exec.html' % program_name,
                                          {'sequence': theseq,
                                           'useremail': theemail,
                                           'workdir': theworkdir,
                                           'log': thelog,
                                           'done': itsdone,
                                           'rank': rank,
                                           'error': theerr})
        else:
            # error
            return render_to_response('%s/index.html' % program_name, {'inputform': form})

    # initial form load
    return render_to_response('%s/index.html' % program_name, {'inputform': FormType()})

#
# end index
#
#################################################


#################################################
#
# pulls up a display of existing results in the
# structure viewer
#

class PdbMenuItem:
    value = ''
    displaystr = ''

class PdbAlignedCol:
    querychar = ''
    querypos  = 0
    pdbchar   = ''
    pdbpos    = 0
    score     = 0.0
    match     = False

def viewStructResults(request, familyid, pdbid, application):

    myfam = Family.objects.get(pk=familyid)

    familydict = {'urlbase': URL_BASE,
                  'id': familyid,
                  'nseq': myfam.n_sequences,
                  'date': myfam.build_date,
                  'seed': myfam.seed_sequence_header.header,
                  'avec': myfam.average_chars}

    # check for browser compatibility
    browserstr = request.META['HTTP_USER_AGENT']

    if browserstr.find('Safari') > -1:
        familydict['safari'] = browserstr

    ## get list of all pdb objects available ##
    allpdblist = []

    for pdbobj in FamilyPdb.objects.filter(family=familyid):
        pdbstr    = pdbobj.pdb_id + pdbobj.chain_id
        pdbheader = '%s [e-value: %.2e]' % (pdbstr, pdbobj.e_val)

        x = PdbMenuItem()
        x.value = pdbstr
        x.displaystr = pdbheader

        allpdblist.append(x)

    familydict['allpdb'] = allpdblist

    ## get query sequence ##
    familydict['querysequence'] = myfam.seed_sequence_header.sequence.chars

    ## get the best pdb id ##
    if pdbid == '':
        best_evalue = 100000.0

        for pdbobj in FamilyPdb.objects.filter(family=familyid):
            if pdbobj.e_val < best_evalue:
                best_evalue = pdbobj.e_val
                pdbid       = pdbobj.pdb_id + pdbobj.chain_id

    if pdbid == '':
        familydict['pdbid']    = ''
        familydict['pdbdir']   = ''
        familydict['pdb']      = ''
        familydict['pdbchain'] = ''

        # create 'alignment' of query seq #
        alignmentarr = []

        for i in range(len(familydict['querysequence'])):
            myindex = i+1
            x = PdbAlignedCol()
            x.querychar = familydict['querysequence'][i]
            x.querypos  = myindex
            alignmentarr.append(x)

        familydict['alignment'] = alignmentarr

    else:
        familydict['pdbid']    = pdbid
        familydict['pdbdir']   = pdbid[1:3]
        familydict['pdb']      = pdbid[0:len(pdbid)-1]
        familydict['pdbchain'] = pdbid[len(pdbid)-1]

        ## align query sequence to pdb structure ##
        
        # create temporary directory
        theworkdir = tempfile.mkdtemp(prefix='/clusterfs/ohana/software/webserver/temp/')
        os.system('chmod 777 %s' % theworkdir)

        # print query sequence to file
        handle = open('%s/seqs.fa' % theworkdir, 'w')
        print >>handle, '>query'
        print >>handle, familydict['querysequence']
        handle.close()

        # print pdb sequence to file
        cmd = '/clusterfs/ohana/software/bin/fastacmd -d /clusterfs/ohana/external/pdb_rcsb_full -s %s_%s >> %s/seqs.fa' % (familydict['pdb'], familydict['pdbchain'], theworkdir)
        os.system(cmd)

        # print model to model file
        mymodel = FamilyHmm.objects.get(family=familyid).hmm

        handle = open('%s/model.mod' % theworkdir, 'w')
        print >>handle, mymodel
        handle.close()

        # align to model
        cmd = '/clusterfs/ohana/software/bin/align2model %s/seqs -modelfile %s/model.mod -db %s/seqs.fa -sw 2 &> /dev/null' % (theworkdir, theworkdir, theworkdir)
        os.system(cmd)

        # read alignment
        handle = open('%s/seqs.a2m' % theworkdir, 'r')

        line = handle.readline()
        line = handle.readline()

        queryseq = ''

        while line and line[0] != '>':
            line = line.strip()
            queryseq += line
            line = handle.readline()

        line = handle.readline()

        pdbseq = ''

        while line:
            line = line.strip()
            pdbseq += line
            line = handle.readline()

        handle.close()

        # get residue scores #
        residuescores = {}       

        res_score_objs = ResidueScore.objects.filter(family=familyid, method=application)

        for i in res_score_objs:
            residuescores[i.residue_number] = i.score

        # build alignment array #

        alignmentarr = []

        pdbinquery = ''

        curr_query_pos = 0
        curr_pdb_pos   = 0

        for i in range(len(queryseq)):
            querychar = queryseq[i]
            pdbchar   = pdbseq[i]

            x = PdbAlignedCol()
            x.querychar = querychar
            x.pdbchar   = pdbchar

            if querychar != '-' and querychar != '.':
                curr_query_pos += 1
                x.querypos = curr_query_pos
                
                if curr_query_pos in residuescores:
                    x.score = residuescores[curr_query_pos]


            if pdbchar != '-' and pdbchar != '.':
                curr_pdb_pos += 1
                x.pdbpos = curr_pdb_pos

            if x.querypos != 0 and x.pdbpos != 0 and x.querychar.isupper() and x.pdbchar.isupper():
                pdbinquery += '%d:%s,' % (x.pdbpos, familydict['pdbchain'])
                x.match = True

            alignmentarr.append(x)

        familydict['alignment'] = alignmentarr

        if len(pdbinquery) > 0:
            pdbinquery = pdbinquery[0:len(pdbinquery)-1]
            familydict['pdbinquery'] = 'select %s; color yellow;' % pdbinquery

        else:
            familydict['pdbinquery'] = ''
        # delete temporary directory
        os.system('rm -r %s' % theworkdir)

    return render_to_response('%s/results.html' % application, familydict)

#
# end viewStructResults
#
#################################################

#################################################
#
# viewSequenceScores
#

def viewSequenceScores(request, familyid, application):

    thedict = {}    

    # get residue scores
    scores_objs = ResidueScore.objects.filter(family = familyid, method=application).order_by('residue_number') 

    if scores_objs:
        myscores = ''

        for scoreobj in scores_objs:
            myscores += '%s\t%s\n' % (scoreobj.residue_number, scoreobj.score)

        thedict['scores'] = myscores

    else:
        thedict['scores'] = 'no %s residue scores are available for protein family %s' % (application, familyid)

    return render_to_response('common/scores.html', thedict)

#
# end viewSequenceScores
#
#################################################

