python ${RNAVIERDATA}/rnavier-configini.py ./code2d.cfg
./test_code2d.exe code2d.ini
gnuplot -p ./test_code2d.gpi 
