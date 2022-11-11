gfortran -O3 -fPIC -march=native -std=legacy -fopenmp -c stokes_sherman_lauricella_solver.f -o stokes_sherman_lauricella_solver.o -L/usr/local/lib -lfmm2d
gfortran -O3 -fPIC -march=native -std=legacy -fopenmp bhgmres_dr.f stokes_sherman_lauricella_solver.o -o int2 -L/usr/local/lib -lfmm2d
./int2
