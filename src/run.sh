rm mesh.* log
#make -f makefile
echo ' '
echo 'Running...'


./cons1d_cons.exe >& log_consoverset
mkdir ConsOverset
mv log* mesh* ConsOverset
echo 'Done Cons Overset'

./cons1d_baseline.exe >& log_baseline
mkdir Baseline
mv log* mesh* Baseline
echo 'Done Baseline '
