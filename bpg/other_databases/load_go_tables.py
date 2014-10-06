#!/usr/bin/env python

import os, glob, datetime, commands

def do_command(logf, cmd):
  logf.write(cmd)
  logf.write('\n')
  status, output = commands.getstatusoutput(cmd)
  if status != 0:
    logf.write('Command exited with nonzero status %d\n'
                % status)
  logf.write("Command output:\n")
  logf.write(output)
  logf.write('\n')

def main():
  logf = open('/clusterfs/ohana/external/tmp_downloads/load_go_tables.log', 'a')
  logf.write('---')
  logf.write(datetime.date.today().strftime('%Y-%m-%d\n'))
  os.chdir('/clusterfs/ohana/external/go/current')
  dir = glob.glob('go_*-assocdb-tables')[0]
  os.chdir(dir)
  all_table_files = set(glob.glob('*.sql'))
  previous_table_files = set(glob.glob('go*.sql'))
  table_files = all_table_files - previous_table_files
  table_names = [os.path.splitext(name)[0] for name in table_files]
  for name in table_names:
    cmd = "cat %s.sql | sed -e 's/CREATE TABLE `/CREATE TABLE `go_/'" % name \
          + " | sed -e 's/`/\"/g' " \
          + "> go_%s.mysql" % name
    do_command(logf, cmd)
    cmd = 'echo "SET ROLE bpg_user;" > go_%s.sql' % name
    do_command(logf, cmd)
    cmd = "/clusterfs/ohana/software/bin/my2pg.pl go_%s.mysql " % name \
          + " | sed -e 's/INT4 NOT NULL auto_increment/SERIAL NOT NULL/' " \
          + ">> go_%s.sql" % name
    do_command(logf, cmd)
    cmd = "psql -h db -f go_%s.sql pfacts003_test" % name
    do_command(logf, cmd)
    cmd = "cat %s.txt | sed -e 's//_/g' > go_%s.txt" % (name, name)
    do_command(logf, cmd)
  of = open("go_load_tables.sql", "w")
  of.write("SET ROLE bpg_user;\n")
  for name in table_names:
    of.write("\COPY go_%s FROM go_%s.txt\n" % (name, name))
  of.close()
  cmd = "psql -h db -f go_load_tables.sql pfacts003_test"
  do_command(logf, cmd)

if __name__ == '__main__':
  main()
