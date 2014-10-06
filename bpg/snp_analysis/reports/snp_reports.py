'''
This is a module of commands used to produce reports.
Author: Curt Hansen
Created: May 22, 2012
Modified:
'''


import os,sys
import pf_connect as db
import latex_functions as l
import file_functions as f
from datetime import datetime


def run_query(query):
    
    dbConn=db.connect_to_server('pfacts003_test','webuser')
    dbCur=db.get_cursor(dbConn)
    dbConn.set_isolation_level(0)
    return db.get_query_list(dbCur,query)
    

def structure_table(content):

    rowheads = []
    colheads = []
    for rowtag,coltag,value in content:
        if rowtag not in rowheads:
            rowheads.append(rowtag)
        if coltag not in colheads:
            colheads.append(coltag)
    rowheads.sort()
    colheads.sort()
    data = [[' ']*len(rowheads) for x in xrange(len(colheads))]
    for rowtag,coltag,value in content:
        data[rowheads.index(rowtag)][colheads.index(coltag)] = value
    return [rowheads,colheads,data]


def create_rost_ressub_table(dir='~/snp_analysis_work/reports'):

    startTime=datetime.now()
    query = 'SELECT resref, resmut, cast(cast(SUM(label) as float)/cast(COUNT(*) as float) as decimal(8,2)) as pct '+\
            'FROM snp_rost_detail GROUP BY resref, resmut ORDER BY resref, resmut;'
    result = run_query(query)
    [rowheads,colheads,data] = structure_table(result)
    caption = 'Each entry indicates the percentage of substitutions from the reference residue '+\
              'to the mutation that are non-neutral.'
    l.create_latex_table('rostsub.tbl',dir,rowheads,colheads,data,'tbl:rostsubs',caption,\
                         subcaption='Rost dataset Substitutions',comment=startTime.strftime("%Y%m%d.%Hh%M"))


def create_rost_intrepid_summary_data(dir='~/snp_analysis_work/reports'):

    startTime = datetime.now()
    comments = startTime.strftime("%Y%m%d.%Hh%M")
    results = []
    query = "SELECT i.scoreimpmap,i.scorejsmap,i.scoreremap,i.scoreglobalmap,i.scoreglobalksmap,i.scoreglobalremap "+\
            "FROM intrepid_detail as i, intrepid_header as ih, snp_rost_master as m, snp_rost_detail as d "+\
            "WHERE i.accession=m.accession AND m.rost_prot_name=d.rost_prot_name AND m.status='M' " +\
            "AND i.accession=ih.accession AND i.fam_id=ih.fam_id "+\
            "AND i.position+ih.start_pos-1=d.respos AND d.label=0 AND m.prependm='false';"
    result1 = run_query(query)
    results.append(result1)
    query = "SELECT i.scoreimpmap,i.scorejsmap,i.scoreremap,i.scoreglobalmap,i.scoreglobalksmap,i.scoreglobalremap "+\
            "FROM intrepid_detail as i, intrepid_header as ih, snp_rost_master as m, snp_rost_detail as d "+\
            "WHERE i.accession=m.accession AND m.rost_prot_name=d.rost_prot_name AND m.status='M' " +\
            "AND i.accession=ih.accession AND i.fam_id=ih.fam_id "+\
            "AND i.position+1+ih.start_pos-1=d.respos AND d.label=0 AND m.prependm='true';"
    result2 = run_query(query)
    results.append(result2)

    f.write_csv_file('rost_intrepid_summary_neutral',dir, \
                     "imp map,js map,re map,global map,global ks map,global re map",\
                     results,comments)
                     
    results = []
    query = "SELECT i.scoreimpmap,i.scorejsmap,i.scoreremap,i.scoreglobalmap,i.scoreglobalksmap,i.scoreglobalremap "+\
            "FROM intrepid_detail as i, intrepid_header as ih, snp_rost_master as m, snp_rost_detail as d "+\
            "WHERE i.accession=m.accession AND m.rost_prot_name=d.rost_prot_name AND m.status='M' " +\
            "AND i.accession=ih.accession AND i.fam_id=ih.fam_id "+\
            "AND i.position+ih.start_pos-1=d.respos AND d.label=1 AND m.prependm='false';"
    result1 = run_query(query)
    results.append(result1)
    query = "SELECT i.scoreimpmap,i.scorejsmap,i.scoreremap,i.scoreglobalmap,i.scoreglobalksmap,i.scoreglobalremap "+\
            "FROM intrepid_detail as i, intrepid_header as ih, snp_rost_master as m, snp_rost_detail as d "+\
            "WHERE i.accession=m.accession AND m.rost_prot_name=d.rost_prot_name AND m.status='M' " +\
            "AND i.accession=ih.accession AND i.fam_id=ih.fam_id "+\
            "AND i.position+1+ih.start_pos-1=d.respos AND d.label=1 AND m.prependm='true';"
    result2 = run_query(query)
    results.append(result2)

    f.write_csv_file('rost_intrepid_summary_nonneutral',dir, \
                     "imp map,js map,re map,global map,global ks map,global re map",\
                     results,comments)
                     
