#!/bin/env python

import sys, os
import csv

reader = csv.reader(sys.stdin)

print "<html>"
print "<table border=1 cellpadding=2>"


for row in reader:
  print "<tr>"
  for item in row:
    print "<td>%s</td>" % item
  print "</tr>"
 


print "</table>"
print "</html>"
