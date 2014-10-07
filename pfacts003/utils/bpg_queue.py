"""

Some functions for dealing with our queue.  We should migrate these to use
PyPBS when that is figured out (the pbs_statque currently gives a seg fault).

"""

import os, sys, subprocess

POSSIBLE_QUEUES = ['web','research','library'] 

def qstat_full_parse(rec):
    lines = rec.split("\n")
    ret = {}
    for l in lines:
        line = l.lstrip()
        if line.startswith("Job Id: "):
            ret['jobId'] = line.lstrip("Job Id: ").split(".")[0]
        elif line.startswith("Job_Name = "):
            ret['job_name'] = line.lstrip("Job_Name = ")
        elif line.startswith("Job_Owner = "):
            ret['job_owner'] = line.split()[-1].split("@")[0]
        elif line.startswith("job_state = "):
            code = line.lstrip("job_state = ")
            if code == 'Q':
                ret['job_status'] = "Queued"
            elif code == 'R':
                ret['job_status'] = "Running"
            elif code == 'C':
                ret['job_status'] = "Completed"
            elif code == 'E':
                ret['job_status'] = "Exiting"
            elif code == 'H':
                ret['job_status'] = "Held"
            elif code == "T":
                ret['job_status'] = "Moving"
            elif code == "W":
                ret['job_status'] = "Waiting"
            elif code == "S":
                ret['job_status'] = "Suspended"
        elif line.startswith("qtime = "):
            ret['qtime'] = line.lstrip("qtime = ")
        elif line.startswith("Walltime.Remaining = "):
            ret['walltime_remaining'] = int(line.lstrip('Walltime.Remaining = '))
        elif line.startswith("resources_used.cput = "):
            ret['cpu_time'] = line.lstrip("resources_used.cput = ")
        elif line.startswith("resources_used.mem = "):
            ret['mem_used'] = int(line.lstrip("resources_used.mem = ").rstrip("kb"))*1024
        elif line.startswith("resources_used.vmem = "):
            ret['virtual_mem_used'] = int(line.lstrip("resources_used.vmem = ").rstrip("kb"))*1024
        elif line.startswith("resources_used.walltime = "):
            ret['walltime_used'] = line.lstrip("resources_used.walltime = ")
        elif line.startswith("queue = "):
            ret['queue'] = line.strip("queue = ")
    return ret

def ohana_job_status(id):
    # gets the status of a single job
    args = ["qstat", "-f", "%d" % id]
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
    (out, err) = p.communicate()
    return qstat_full_parse(out)

def ohana_single_queue_status(queue):
    args = ["qstat", "-f", queue]
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
    (out, err) = p.communicate()
    ret = {}
    if not out:
        return ret
    records = out.split("\n\n")
    for record in records:
        r = qstat_full_parse(record)
        if r:
            ret[r['jobId']] = r
    return ret

def ohana_queue_status(queue = None, jobId = None):
    # This function returns the status of a queue on ohana
    if jobId:
        return ohana_job_status(jobId) 
    
    if queue:
        return ohana_single_queue_status(queue)

    ret = {}

    for q in POSSIBLE_QUEUES:
        ret[q] = ohana_single_queue_status(q)

    return ret
