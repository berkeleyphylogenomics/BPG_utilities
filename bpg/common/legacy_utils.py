
# Misc generic utilities

import string, os, sys, fnmatch, stat, time, tempfile, socket, pwd

RunCmdVerbose = 0
ErrorTraceback = 0
PrintCmdTimes = 0
PBSNoSubmit = 0
PBSMaxRunning = 0

# -----------------------------------------------------------------

def Error(message) :
    print "Directory:", os.getcwd()
    print "Host:", socket.gethostname()
    if ErrorTraceback :
        raise message
    else :
        print "Error:", message
        sys.exit(1)

# -----------------------------------------------------------------

# If if_error = 1 then throw an exception if there is a nonzero exit
# status.  If it's 0 then print an error message but continue.  If it's
# -1 then ignore it entirely.

def RunCmd(cmd, out = "", verbose = "", if_error = 1) :
    if PrintCmdTimes :
        verbose = 1
    elif verbose == "" :
        verbose = RunCmdVerbose
    if out == "stdout" or cmd.find('>') >= 0 :
        cmd1 = cmd
    elif not out :
        cmd1 = cmd + " >/dev/null 2>&1"
    else :
        cmd1 = cmd + " >%s 2>&1" % out
        cmd = cmd1 # looks better when echoing
    if PrintCmdTimes :
        startt = time.time()
    if verbose : print "*RUN: ", cmd
    ret = os.system( "source /etc/profile; " + cmd1 )
    if ret :
        if if_error == 1 :
            Error("%s exited with status %d" % (cmd1, ret))
        elif if_error == 0 :
            print "Warning: %s exited with status %d" % (cmd1, ret)

    if PrintCmdTimes :
        endt = time.time()
        Print("*TIME: %g seconds" % (endt - startt))

    return ret

# -----------------------------------------------------------------

def Chdir(dir, verbose = 0) :
    if verbose : print "Chdir:", dir
    os.chdir(dir)
    
# -----------------------------------------------------------------

def Plural(base, num) :
    if num == 1 :
        return "1 %s" % base
    elif base[-1] == 's' :
        return "%d %ses" % (num, base)
    else :
        return "%d %ss" % (num, base)

# -----------------------------------------------------------------

def Plural(base, num) :
    if num == 1 :
        return "1 %s" % base
    elif base[-1] == 's' :
        return "%d %ses" % (num, base)
    else :
        return "%d %ss" % (num, base)

# -----------------------------------------------------------------

def CapFirst(str) :
    return str[0].upper() + str[1:]

# -----------------------------------------------------------------

def Whoami() :
    return pwd.getpwuid(os.getuid())[0]

# -----------------------------------------------------------------

def Print(*args) :
    str = ""
    for a in args :
        if str : str = str + " "
        str = str + "%s" % a
    print str
    sys.stdout.flush()

# -----------------------------------------------------------------

def PbsSubmitCmds(cmds, outfile, tag, priority = 0) :
    tag = tag[:15]
    if not PBSNoSubmit and PbsRunning(tag) :
        print tag, "is already running"
        return
    if PBSMaxRunning and PbsNumRunning() >= PBSMaxRunning :
        secs = 2
        while 1 :
            time.sleep(secs)
            if secs < 64 : secs *= 2
            nr = PbsNumRunning()
            if nr < PBSMaxRunning : break
            if secs == 4 : print nr, "running, sleeping"
            
    if not os.access("qsub", os.F_OK) :
        os.system("mkdir qsub")

    job = "qsub/" + tag
    fp = file(job, "w")
    fp.write("#!/bin/csh -f\n")
    fp.write("cd %s\n" % os.getcwd())
    nc = 0
    if PBSNoSubmit :
        into_out = ""
    else :
        into_out = ">>& " + outfile
        
    for cmd in cmds :
        if not nc and not PBSNoSubmit : 
            fp.write("rm -f %s\n" % outfile)

        written = 0
        for pat in ["foreach ", "end", "cd ", "set "] :
            if cmd.strip().startswith(pat) :
                fp.write("%s\n" % (cmd))
                written = 1
                break

        if not written :
            fp.write("echo \"COMMAND: %s\" %s\n" % (cmd, into_out))
            if cmd.find('>') >= 0 :
                fp.write("%s\n" % (cmd))
            else :
                fp.write("%s %s\n" % (cmd, into_out))

            fp.write("echo ===================================== %s\n" \
                                                                   % (into_out))
            if not PBSNoSubmit :
                fp.write("echo >> %s\n" % (outfile))

        nc = nc + 1

    fp.write("exit 0\n")
    fp.close()

    os.system("chmod 0755 %s" % job)
    if PBSNoSubmit :
        os.system(job)
    else :
        qcmd = "qsub -cwd -p %d -N %s -e %s.stderr -o %s.stdout %s" % (priority, tag, outfile, outfile, job)
        print qcmd
        os.system(qcmd)
    
# -----------------------------------------------------------------

def PbsSubmit(cmd, outfile, tag, priority = 0) :
    if PBSNoSubmit :
        # Just run it locally
        print "RUNNING:", cmd
        os.system(cmd)
        return

    PbsSubmitCmds([cmd], outfile, tag, priority)
    
# -----------------------------------------------------------------

def PbsRunning(tag) :
    for line in os.popen("qstat", "r") :
        if line.find(tag) >= 0 : return 1
    return 0
    
# -----------------------------------------------------------------

def PbsNumRunning(user = Whoami()) :
    num = 0
    for line in os.popen("qstat", "r") :
        if line.find(user) >= 0 : num += 1
    return num
    
# -----------------------------------------------------------------

def SameFilesystem(file1, file2) :
    return os.stat(file1)[stat.ST_DEV] == os.stat(file2)[stat.ST_DEV]

# -----------------------------------------------------------------

def FileNonEmpty(fname) :
    if not os.access(fname, os.F_OK) : return 0
    if not os.stat(fname)[6] : return 0
    return 1

# -----------------------------------------------------------------

def SendEmail(recip, subject, text) :
    mfile = tempfile.mktemp()
    fp = file(mfile, "w")
    fp.write(text + "\n")
    fp.close()
    cmd = 'mail -s "%s" %s < %s' % (subject, recip, mfile)
    if RunCmdVerbose : print cmd
    os.system(cmd)
    os.system("rm %s" % mfile)
