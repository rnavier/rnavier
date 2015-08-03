python ${RNAVIERDATA}/rnavier-configini.py icbj.cfg
./test_icbj.exe icbj.ini
./test_icbj.exe --no-entropy-rewrite icbj.ini
./test_icbj_hydro.exe icbj.ini
gnuplot -p test_icbj.gpi


