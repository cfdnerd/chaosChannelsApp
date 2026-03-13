git pull
cd ../src/
rm -rf Make/linux*
wmake -j
cd ../app/
git pull
# for copying purposes: wmake -j 2>&1 | tee build.log
