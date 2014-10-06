#!/usr/bin/python
'''
Updates the static coverage html file.

Input arguments are the csv file (generated using outputs from simple_genome_coverage.py) and the html file to be updated.
(usually at ohana_repository/pfacts003/templates/phylofacts/coverage.html)
output will be called coverage2.html in the directory where coverage.html resides. Rename this to coverage.html when done and delete the older copy.
'''
import os
import sys
import re

class StaticCoverageUpdate:
    def __init__(self, *args):
        '''Initialize arguments. Exits if the order of arguments is mismatched.'''
        if not args[0].endswith('.csv'):
            print 'The first argument should have a csv extension.'
            sys.exit()
        if not args[1].endswith('.html'):
            print 'The second argument should have a html extension.'
            sys.exit()                            


        self.coverage_csv = args[0]
        self.coverage_html = args[1]
        self.output_filename = self.coverage_html[:-5] + '1.html'

    def read_csv_file(self):
        '''Reads csv file into a dictionary'''
        import csv
        self.coverage_dict = {}
        csv_reader = csv.reader(open(self.coverage_csv, 'rU'))
        # Ignore the header line
        csv_reader.next() 
        self.csv_dict = {}
        for row in csv_reader:
            self.csv_dict[row[0].strip()] = row[1:]
        
    def read_html_file(self):
        '''Reads html file into array'''
        self.html_array = open(self.coverage_html, 'rU').readlines()


    def is_species_line(self, line):
        '''Is this html line one with species info'''
        return ('<td class="species">' in line)
    
    def map_csv_onto_html(self):
        '''Maps csv dictionary info onto the html file.'''
        #Ideally sgmlparser should be used. Here I do simple regular expression
        #substitutions
        count = 0
        unset = set()
        orig_re_val = 'class\s*=\s*\"species\"\s*>(.*?)</td><td>(.*)</td><td>(.*)</td><td>(.*?)</td>'
        re_val = re.compile(orig_re_val)
        replace_re_val = 'class="species">\\1</td><td>%s</td><td>%s</td><td>%s</td>'
        for num, line in enumerate(self.html_array):
            if self.is_species_line(line):
                species_values = re_val.findall(line)
                for species_value in species_values:
                    for csv_key, csv_item in self.csv_dict.items():
                        if species_value[0] in csv_key:
                            new_line = re.sub(orig_re_val,
                                              replace_re_val % (csv_item[0], csv_item[-2],
                                                                csv_item[-1]),
                                              line)
                            self.html_array[num] = new_line

    def write_new_html_file(self):
        '''Write self.html_array to a new file'''
        f = open(self.output_filename, 'w')
        f.write('\n'.join(self.html_array))
        f.close()
        print "Wrote new information into %s." % self.output_filename
    
    def main(self):
        '''Main function that glues together other commands'''
        self.read_html_file()
        self.read_csv_file()
        self.map_csv_onto_html()
        self.write_new_html_file()
        
if __name__ == "__main__":
    if not len(sys.argv) == 3:
        print 'Usage: python %s csv_file html_file' % sys.argv[0]
        sys.exit()
    scu = StaticCoverageUpdate(sys.argv[1], sys.argv[2])
    scu.main()
    
