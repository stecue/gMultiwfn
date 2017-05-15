#!/bin/bash
dirdiff=vdiff
#First let's generate the patch
if [ -d $dirdiff ];then
    rm -rf $dirdiff
fi
mkdir $dirdiff
cd $dirdiff
pushd .
#Get the old src
ln -s ../current_src.zip ./
unzip current_src.zip
cd Mult*
find . -name '* *' -print0|xargs -0 rm -f
rm GUI.f90 && rm *.a
sed -i -e 's/grep/grep -a/g' noGUI.sh
sh noGUI.sh && mv noGUI ../src_old
popd && mv Multi*src_Linux src_full_old
pushd .
#Get the new src
ln -s ../latest_src.zip ./
unzip latest_src.zip
cd Mult*
find . -name '* *' -print0|xargs -0 rm -f
rm GUI.f90 && rm *.a
sed -i -e 's/grep/grep -a/g' noGUI.sh
sh noGUI.sh && mv noGUI ../src_new
popd && mv Multi*src_Linux src_full_new
../gnuize.sh src_old
../gnuize.sh src_new
diff -uNr src_old src_new > src_update.diff
#backup the original src
cd ..
if [ -d src.orig ];then
    rm -r src
    cp -pr src.orig src
else
    cp -pr src src.orig
fi
#apply the patch
cd src
patch -p1 < ../$dirdiff/src_update.diff
cd ..
# rename .zip #should be down manually
#mv current_src.zip previous_src.zip
#mv latest_src.zip current_src.zip
