python ${RNAVIERDATA}/rnavier-configini.py gaus_x.cfg
python ${RNAVIERDATA}/rnavier-configini.py gaus_y.cfg

./test_gaus.exe gaus_x.ini
./test_gaus.exe gaus_y.ini
gnuplot -p test_gaus.gpi

