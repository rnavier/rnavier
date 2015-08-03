
python ${RNAVIERDATA}/rnavier-configini.py perturb_x.cfg
python ${RNAVIERDATA}/rnavier-configini.py perturb_y.cfg
python ${RNAVIERDATA}/rnavier-configini.py perturb_z.cfg
./test_perturb.exe perturb_x.ini
./test_perturb.exe perturb_y.ini
./test_perturb.exe perturb_z.ini
gnuplot -p -e "xyz='x'; idir=0;" test_perturb.gpi
gnuplot -p -e "xyz='y'; idir=1;" test_perturb.gpi
gnuplot -p -e "xyz='z'; idir=2;" test_perturb.gpi

