#!/usr/bin/python

from bookends import Bookend
import os
import re
import subprocess
import pickle

newliner = re.compile(r'(.{70})')
fasta_seq_trash = re.compile(r'[\s\.-]')

def parse_astats(text):
    astats_dict = {}
    for num, section in enumerate(re.compile('^-{2,}$', re.MULTILINE).split(text)):
        sub_dict = {}
        for line in re.findall(r'(\w+( [()\w]+)+ {2,}\d+(\.\d+)?%?)\n', section):
            key, data = re.split(r' {2,}', line[0])[0:2]
            sub_dict[re.sub(r'\([\w ]+\)', '', key).strip().replace(' ','_')] = data
        astats_dict[('general','all','aligned','general')[num]] = sub_dict
    return astats_dict


class PostFlowerPower(Bookend):

    def __init__(self):
        super(PostFlowerPower,self).__init__('flowerpower', r'FlowerPower retrieved \d+ sequences\.\n$')

    def __run__(self):

        print >>self.out, "prettyaligning the multiple sequence alignment ...",
        os.system('prettyalign final.a2m -f > flowerpower_alignment.fa 2> /dev/null')
        print >>self.out, " done"

        print >>self.out, "Calculating astats for FlowerPower alignment ...",
        handle = open('flowerpower_alignment_astats.p', 'w')
        pickle.dump(parse_astats(subprocess.Popen(
            ['astats', 'flowerpower_alignment.fa'],
            stdout=subprocess.PIPE,
        ).stdout.read()),handle)
        handle.close()
        print >>self.out, " done"

        print >>self.out, "Writing unaligned sequences ...",
        handle = open('final.a2m', 'r')
        final_a2m = handle.read()
        handle.close()

        # Surely there is a better way of doing this...
        handle = open('flowerpower_sequences_unaligned.fa', 'w')
        handle.writelines(reduce(
            # Flatten the 2-tuple
            lambda x,y: x+y,
            # Returns a list of FASTA header,sequence,header,sequence...etc
            zip(
                # FASTA headers
                re.compile(r'>.*\n').findall(final_a2m),
                # FASTA sequences
                (
                    # Add newlines to the fasta sequence
                    newliner.sub(
                        r'\1\n',
                        # Remove junk characters; sequence becomes a single line
                        fasta_seq_trash.sub('', seq)
                    )+'\n' for seq in \
                    re.findall(r'\n[AC-IK-NP-TVWY\s.-]+', final_a2m)
                ),
            )
        ))
        handle.close()
        print >>self.out, " done"
        
        print >>self.out, "Realigning sequences with MUSCLE ...",
        if self.cleaned_data['mode'] == 'global':
            os.system('muscle -in final.a2m -out muscle_realignment.fa -maxiters 2 -verbose >& muscle.log')
            print >>self.out, " done"
            print >>self.out, "Generating astats for MUSCLE alignment ...",
            handle = open('muscle_realignment_astats.p', 'w')
            pickle.dump(parse_astats(subprocess.Popen(
                ['astats', 'muscle_realignment.fa'],
                stdout=subprocess.PIPE,
            ).stdout.read()),handle)
            handle.close()
            print >>self.out, " done"
        else:
            print >>self.out, " skipping"
        return 1

if __name__ == '__main__':
    PostFlowerPower().run()
