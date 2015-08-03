#Test Bjorken expansion
python ${RNAVIERDATA}/rnavier-configini.py icgb.cfg
python ${RNAVIERDATA}/rnavier-configini.py icgb2.cfg
./test_icgb_hydro.exe icgb.ini
./test_icgb_hydro.exe icgb2.ini
gnuplot -p test_icgb.gpi &
gnuplot -p test_icgb_ideal.gpi
