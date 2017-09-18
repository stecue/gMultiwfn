#!/bin/bash
#Try update_src_delta.sh first!
#This will apply the "gnu" patch to the new src.
dirdiff=vdiff
function clean_noGUI {
    rm GUI.f90 plot.f90 && rm *.a
}
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
#rm file not mentioned in the noGUI Makefile
clean_noGUI
sed -i -e 's/grep/grep -a/g' noGUI.sh
sh noGUI.sh && mv noGUI ../src_old
popd && mv Multi*src_Linux src_full_old
../gnuize.sh src_old
pushd .
#Get my patch for current src
if [ -d ../src.orig ]; then
    cp -pr ../src.orig src_gnu
else
    cp -pr ../src src_gnu
fi
diff -uNr src_old src_gnu > patch_gnu.diff
#Get the new src and get the delta diff from the old to new ifort src
ln -s ../latest_src.zip ./
unzip latest_src.zip
cd Mult*
find . -name '* *' -print0|xargs -0 rm -f
clean_noGUI
sed -i -e 's/grep/grep -a/g' noGUI.sh
sh noGUI.sh && mv noGUI ../src_new
popd && mv Multi*src_Linux src_full_new
../gnuize.sh src_new
diff -uNr src_old src_new > src_update.diff
#backup the original src
# cd ..
# if [ -d src.orig ];then
#     rm -r src
#     cp -pr src.orig src
# else
#     cp -pr src src.orig
# fi
# #apply the patch
# cd src
# patch -p1 < ../$dirdiff/src_update.diff
cd src_new && patch -p1 --no-backup-if-mismatch < ../patch_gnu.diff
cd .. #we are in $dirdiff
rm -r ../src
cp -r src_new ../src
# rename .zip #should be down manually
#mv current_src.zip previous_src.zip
#mv latest_src.zip current_src.zip
cd .. #go back to the project home folder
echo "Number of rejected patches is: `ls -l src/*rej 2> /dev/null | wc -l`"
