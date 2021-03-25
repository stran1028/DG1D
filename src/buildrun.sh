rm *.o *.mod fort.* mesh.* log
make -f makefile
./cons1d.exe > log
