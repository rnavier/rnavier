#Test 2d sod problem in all directions
python ${RNAVIERDATA}/rnavier-configini.py sod2d_xy.cfg
python ${RNAVIERDATA}/rnavier-configini.py sod2d_yz.cfg
python ${RNAVIERDATA}/rnavier-configini.py sod2d_xz.cfg
./test_sod_2d.exe sod2d_xy.ini
./test_sod_2d.exe sod2d_yz.ini
./test_sod_2d.exe sod2d_xz.ini
gnuplot -p test_sod_2d.gpi
