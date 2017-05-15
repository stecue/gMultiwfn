#!/bin/bash
# Written by Xing (stecue@gmail) to do some "porting" work for gMultiwfn.
srcDir=$1
for srcF in ${srcDir}/*f90
do
##Fix tab using blanks
    dos2unix "${srcF}"
#    ntab=`grep -P '^ *\t' "${srcF}"|wc -l`
#    while [ $ntab -gt 0 ]
#    do
#        sed -i -e 's/^\( *\)\t/\1    /g' "${srcF}"
#        ntab=`grep -P '^ *\t' "${srcF}"|wc -l`
#    done
    sed -i -e 's/\t/    /g' "${srcF}"
    sed -i -e 's/\(^ *\)read *(\*,\*)\( *\)$/\1pause\2/g' "${srcF}"
done
