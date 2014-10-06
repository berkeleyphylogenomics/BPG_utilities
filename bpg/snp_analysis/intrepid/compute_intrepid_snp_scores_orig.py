#!/clusterfs/ohana/software/bin/python2.6 -Wignore::DeprecationWarning
'''
This script takes a bpg accession and a uniprot accession, and run interpid to
score the residues in the sequence corresponding to the uniprot accession

Inputs: "snp_rost_data_unique" table in pfacts003_test database
Outputs: "snp_rost_intrepid_header" and "snp_rost_intrepid_detail" tables in pfacts003_test database 

Author: Curt Hansen
Created: 20Apr2012
Modified: 
'''


import sys
sys.path.append('../RostLab') #Need this for PFConnect
from bpg.common.utils.dir_of_family import get_dir_of_family_accession as get_dir #Need this for Intrepid
import glob, os, shutil, time
import PFConnect as db
from datetime import datetime

class NullDevice: #Need this to capture print messages.
    def write(self,s):
        pass

def get_file(family_dir,extension):
    '''get the msa and tree file'''
    file=glob.glob(os.path.join(family_dir,extension))
    if len(file) != 1:
        return 'error'
    else:
        return file[0]

def CheckInputs(protid,familyid):
    OKstatus=['Y','Y','Y'] #Initialize to all OK.
    data=[]
    family_dir = get_dir(familyid)
    msa_file = get_file(family_dir, '*mafft.afa')
    if msa_file=='error':
        OKstatus[0]='N'
    else:
        data.append(msa_file)
    tree_file = get_file(family_dir, '*ml.rooted.tre')
    if tree_file=='error':
        OKstatus[1]='N'
    else:
        data.append(tree_file)
    idmap_file = get_file(family_dir, '*ungapped.idmap')
    if idmap_file=='error':
        OKstatus[2]='N'
    else:
        data.append(idmap_file)
    return [OKstatus,data]

def ScoreProtein(protid,data):
    #data is list with msa_file,tree_file,idmap_file.
    cmd = "python "
    cmd+="/home/curth/ohana_repository/bpg/intrepid/run_intrepid_against_manuscripts.py "
    cmd+="%s %s %s %s" % (protid,data[0],data[1],data[2])  
    os.system(cmd)

def main():
    if len(sys.argv) != 4:
        print "Usage: %s <start num> <end num> <printstep>" % sys.argv[0]
        sys.exit(1)

    startprot=int(sys.argv[1])
    endprot=int(sys.argv[2])
    printstep=int(sys.argv[3])
    sys.stdout=NullDevice() #Turn off print messages.
    startTime=datetime.now()
    pgmName='ComputeIntrepidSNPScores'
    fout=open(pgmName+".out."+startTime.strftime("%Y%m%d.%Hh%M"),"w")
    fout.write("\nStarting program...\n")
    fout.write("Turned off deprecation warnings.\n")
   
    dbConn=db.Conn('pfacts003_test','webuser')
    dbCur=db.Cursor(dbConn)
    dbConn.conn.set_isolation_level(0)

    try:
        dbCur.execute("DROP TABLE snp_rost_intrepid_header;")
        fout.write("Previous version of 'snp_rost_intrepid_header' deleted.\n")
    except:
        fout.write("Table 'snp_rost_intrepid_header' does not exist and will be created.\n")
    try:
        dbCur.execute("DROP TABLE snp_rost_intrepid_detail;")
        fout.write("Previous version of 'snp_rost_intrepid_detail' deleted.\n")
    except:
        fout.write("Table 'snp_rost_intrepid_detail' does not exist and will be created.\n")

    dbCur.execute("CREATE TABLE snp_rost_intrepid_header \
          (id serial PRIMARY KEY, pfactsid integer, protaccid varchar, familyid integer, \
          familytype char(1), msaok char(1), treeok char(1), idmapok char(1), familyscore integer, starttime varchar);")
    fout.write("Created new or replacement 'snp_rost_intrepid_header' table in database.\n")

    dbCur.execute("CREATE TABLE snp_rost_intrepid_detail \
      (id serial PRIMARY KEY, pfactsid integer, protaccid varchar, familyid integer, position integer,resseq integer,\
      chain varchar,residue char(1),icode char(1),scoreimpmap float,scorejsmap float,scoreremap float,\
      scoreglobalmap float,scoreglobalksmap float,scoreglobalremap float);")
    fout.write("Created new or replacement 'snp_rost_intrepid_detail' table in database.\n")

    proteinsList=dbCur.sqlQueryList("SELECT u.id, u.accession FROM uniprot as u INNER JOIN snp_rost_unique as s \
    ON s.pfactsid=u.id WHERE s.matchstatus='M';")

    totCounter=0
    totResCounter=0

    l1=65
    l2=15
    fout.write("\nProgress Details\n")
    fout.write("-"*(l1+l2)+"\n")
    for prot in proteinsList:
        if totCounter+1>endprot or totCounter+1<startprot:
            break
        
        protPfactsID=prot[0]
        protAccID=prot[1]

        famsList=dbCur.sqlQueryList("SELECT f.id, f.family_type_id, f.score FROM family_sequence_taxa_2012_02_17 as fs \
        LEFT JOIN family as f ON fs.family_id=f.id WHERE fs.uniprot_accession=%s \
        ORDER BY f.family_type_id DESC, f.score DESC;",[protAccID])
        if len(famsList)==0:
            famid=0
            famtype='N'
            famscore=0
            msaok='N'
            treeok='N'
            idmapok='N'
            dbCur.insert("INSERT INTO snp_rost_intrepid_header \
              (pfactsid, protaccid, familyid, familytype, msaok, treeok, idmapok, familyscore, starttime) \
              VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s);",\
              [protPfactsID,protAccID,famid,famtype,msaok,treeok,idmapok,famscore,startTime])
        else:
            flgGHGFound=False
            for fam in famsList:
                famid=fam[0]
                famtype=fam[1]
                famscore=fam[2]
                if famtype=='C' or (famtype=='G' and flgGHGFound==False):
                    if famtype=='G':
                        flgGHGFound=True
                    reldata=CheckInputs(protAccID,"bpg0"+str(famid))
                    msaok=reldata[0][0]
                    treeok=reldata[0][1]
                    idmapok=reldata[0][2]
                    dbCur.insert("INSERT INTO snp_rost_intrepid_header \
                    (pfactsid, protaccid, familyid, familytype, msaok, treeok, idmapok, familyscore, starttime) \
                    VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s);",\
                      [protPfactsID,protAccID,famid,famtype,msaok,treeok,idmapok,famscore,startTime])
        
                    if reldata[0]==['Y','Y','Y']: #Means inputs found for protein, so run Intrepid.
                        ScoreProtein(protAccID,reldata[1]) #Run intrepid for protein / fmaily combo.
                        path='./'+protAccID+'_intrepid/'
                        filepath=path+'output.aux'
                        score_file=open(filepath,'r')
                        linenum=0
                        for line in score_file:
                            linenum+=1
                            if linenum==1:
                                pass #Skip the header row.
                            else:
                                totResCounter+=1 #Increment counter of total number of detail lines added.
                                details=line.split('|') #Parse line.
                                details[10]=details[10][0:-1] #Drop trailing carriage return for last field.
                                dbCur.insert("INSERT INTO snp_rost_intrepid_detail \
                                (pfactsid, protaccid, familyid, position, resseq, chain, residue, icode, scoreimpmap, \
                                scorejsmap, scoreremap, scoreglobalmap, scoreglobalksmap, scoreglobalremap) \
                                VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s);", \
                                  [protPfactsID,protAccID,famid,details[0],details[1],details[2],details[3],\
                                   details[4],details[5],details[6],details[7],details[8],details[9],details[10]])
                        score_file.close()
                        shutil.rmtree(path,ignore_errors=True) #Delete directory of Intrepid results for this protein.

        totCounter+=1
        if totCounter % printstep==0:
            msg="Completed processing entry numnber "+str(totCounter)+" in protein list.\n"
            fout.write(msg)

    fout.write("\nResults\n")
    fout.write("-"*(l1+l2)+"\n")
    fout.write("Number of proteins scored:".ljust(l1)+str(totCounter).rjust(l2)+"\n")
    fout.write("Number of protein positions scored:".ljust(l1)+str(totResCounter).rjust(l2)+"\n") 
    endTime=datetime.now()
    diff=endTime-startTime
    fout.write("Elapsed time (HH:MM:SS):".ljust(l1)+str(diff).rjust(l2)+"\n\n")
    fout.close()
    dbConn.Commit()
    dbConn.Close()


if __name__=='__main__':
  main()
