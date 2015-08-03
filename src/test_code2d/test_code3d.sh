python ${RNAVIERDATA}/rnavier-configini.py ./code3d.cfg
./test_code3d.exe code3d.ini
gnuplot -p ./test_code3d.gpi
