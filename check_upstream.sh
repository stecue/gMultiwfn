#!/bin/bash
LATEST=latest_src.zip
CURRENT=current_src.zip
> $LATEST
ifwrong=`( unzip -l ${LATEST} > /dev/null ) 2>&1 |wc -l`
while [ $ifwrong -gt 0 ]
do
    echo $ifwrong
    wget 'http://sobereva.com/multiwfn/misc/Multiwfn_3.5(dev)_src_Linux.zip' -O $LATEST
#    wget 'http://sobereva.com/multiwfn/misc/Multiwfn_3.4(dev)_src_Linux.zip' -O $LATEST
    ifwrong=`( unzip -l ${LATEST} > /dev/null ) 2>&1 |wc -l`
done
md5 () {
    this_md5=`md5sum -b "$1"|cut -d' ' -f1`
    echo ${this_md5}
}
if [ x`md5 $LATEST` = x`md5 $CURRENT` ]
then
    echo "Noting to do..."
else
    echo "You need to update the build."
    parallel --lb <<EOF
wget 'http://sobereva.com/multiwfn/misc/Manual_3.5(dev).pdf' -O doc/Manual_dev.pdf
wget 'http://sobereva.com/multiwfn/misc/Multiwfn_3.5(dev)_bin_Win64.rar' -O examples.rar
EOF
    unrar x examples.rar
    cp -r Multi*/examples/* examples/
    rm -r Multi*
fi
