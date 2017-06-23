#!/bin/bash
cat << EOF > allNTO.txt
18
6
examples/NTO/uracil.out
EOF
for ((i=1;i<=3;i=i+1))
do
cat << EOF >> allNTO.txt
$i
2
S$i.fch
6
EOF
done
./Multiwfn examples/NTO/uracil.fch < allNTO.txt
rm ./allNTO.txt
