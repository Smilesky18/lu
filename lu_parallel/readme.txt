gcc -O2 demos.c lu_gp_sparse.c parallel.c -L. nicslu.a nicslu_util.a -lrt -lpthread -lm -o demos -fopenmp
./demos G2_circuit.mtx
