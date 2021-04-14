rm *.o *.mod fort.* mesh.* log
make -f makefile
echo ' '
echo 'Running...'

./cons1d.exe > log
