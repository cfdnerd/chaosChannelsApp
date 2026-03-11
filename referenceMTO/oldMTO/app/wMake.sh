git pull
cd ../src_TF/
rm -rf Make/linux*
wmake -j
cd ../app/
git pull
# for copying purposes: wmake -j 2>&1 | tee build.log
