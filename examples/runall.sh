#!/bin/bash
for inf in *.gjf
do
echo Running ${inf} ...
time g09 < ${inf} > ${inf//gjf/out}
echo ${inf} is finished
echo
done
