#!/bin/bash
# Written by Xing (stecue@gmail) to do some "porting" work for gMultiwfn.
srcDir=src.orig
for srcF in ${srcDir}/*f90
do
##Fix tab using blanks
    sed -i -e 's/\t/    /g' "${srcF}"
    sed -i -e 's/\(^ *\)read *(\*,\*)\( *\)$/\1pause\2/g' "${srcF}"
#Note that following will add too many "getNThreads". It should only work for _delta_ update!
    sed -i -e 's/^[\t ]*\!\$/!$/' "${srcF}"
    sed -i -e '/^..OMP PARALLEL.*nthreads/i\
nthreads=getNThreads()' "${srcF}"
    sed -i -e '/^..OMP parallel.*nthreads/i\
nthreads=getNThreads()' "${srcF}"
    sed -i -e 's/num_threads( nthreads  )/num_threads(nthreads)/g' "${srcF}"
    sed -i -e 's/\([^=]\)==\.false\./\1.eqv. .false./g' "${srcF}"
    sed -i -e 's/\(^ *write.*\)nthreads/\1 getNThreads()/g' "${srcF}"
#The following will remove usefull lines.
#    uniq "${srcF}" > "${srcF}.1"
#    mv "${srcF}.1" "${srcF}"
done
