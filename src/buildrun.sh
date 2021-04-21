rm *.o *.mod *exe mesh.* log
make -f makefile
echo ' '
echo 'Running...'

./cons1d.exe > log
echo 'Done'
