Compilation and Run commands are as follows

gcc -g -O4 -fopt-info-optall-optimized -ftree-vectorize -march=native -ftree-vectorizer-verbose=2 -fopt-info-vec-optimized -mtune=native -Wunused-variable matmat_auto.c util.h util.c -o matmat_auto

./matmat_auto -n 100
./matmat_auto -n 200
./matmat_auto -n 400
./matmat_auto -n 800
./matmat_auto -n 1600

gcc -Wall -march=native -mtune=native -Wno-incompatible-pointer-types -Wno-unused-variable main_drivers.c main_drivers.h matvec.c util.h util.c -o matvec
./matvec -n 100
./matvec -n 200
./matvec -n 400
./matvec -n 800
./matvec -n 1600

gcc -Wall -march=native -mtune=native -Wno-incompatible-pointer-types -Wno-unused-variable main_drivers.c main_drivers.h matmat.c util.h util.c -o matmat
./matmat -n 100
./matmat -n 200
./matmat -n 400
./matmat -n 800
./matmat -n 1600