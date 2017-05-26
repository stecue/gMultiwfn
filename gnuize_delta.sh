#!/bin/bash
# Written by Xing (stecue@gmail) to do some "porting" work for gMultiwfn.
srcDir=$1
dos2unix ${srcDir}/*F
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
#Note that following will add too many "getNThreads". It should only work for _delta_ update!
    sed -i -e 's/^[\t ]*\!\$/!$/' "${srcF}"
    sed -i -e '/^..OMP PARALLEL.*nthreads/i\
nthreads=getNThreads()' "${srcF}"
    sed -i -e 's/\([^=]\)==\.false\./\1.eqv. .false./g' "${srcF}"
    sed -i -e 's/\([^=]\)==\.true\./\1.eqv. .true./g' "${srcF}"
    sed -i -e 's/^\(.*call KMP_SET_STACKSIZE_S.*\)/!\1/' "${srcF}"
    # Delete the unsupported domaingui for otherfunction2.f90
    sed -i -e  '/call drawdomaingui/d' "${srcF}"
done
