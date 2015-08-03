python ${RNAVIERDATA}/rnavier-configini.py gubserfo.cfg
./test_fo.exe gubserfo.ini 2.0 1.0 0.001 
gnuplot -p test_gubserfo.gpi
