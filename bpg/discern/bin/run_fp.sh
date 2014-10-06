#! /bin/sh -l

input=$1
flowerpower.pl -i $input -n 10 --mode glocal --tempcheck 1 >& flowerpower.out
removeGappyColumns -i -m 70 final.a2m > makebook_input_alignment
make_nj_tree.py makebook_input_alignment
