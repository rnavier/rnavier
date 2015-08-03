#Test 1d sod problem in x direction
python ${RNAVIERDATA}/rnavier-configini.py sod1d.cfg
./test_sod.exe sod1d.ini
gnuplot -p test_sod_1d.gpi
