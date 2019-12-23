all: main.cpp
	~/intel/bin/icpc -O2 main.cpp Source/lu_gp_sparse.cpp Lib/libcxsparse.a Lib/libsuitesparseconfig.a -fpermissive -fPIC -lsuperlu -lblas -lm -o sparse_lu
main: main.cpp
	g++ -g -O2 main.cpp Source/lu_gp_sparse.cpp -fpermissive Lib/libcxsparse.a Lib/libsuitesparseconfig.a -fPIC -lsuperlu -lblas -lm -o main
test: main.cpp
	gcc -g -O2 -fdump-rtl-expand -fexceptions -fPIC -fopenmp main.cpp Source/lu_gp_sparse.cpp -fpermissive Lib/libcxsparse.a Lib/libsuitesparseconfig.a -lsuperlu -lblas -lm -lrt -W -o sparse_lu
sn_detect: supernode_detect.cpp
	g++ -g -O2 supernode_detect.cpp -fpermissive Lib/libcxsparse.a Lib/libsuitesparseconfig.a -fPIC -lblas -lm -o sn_detect
sn_sort: supernode_sort.cpp
	icpc -O3 supernode_sort.cpp -fpermissive -fPIC -lblas -lm -o sn_sort
block: block.cpp
	g++ -g -O2 block.cpp -fpermissive Lib/libcxsparse.a Lib/libsuitesparseconfig.a -fPIC -lblas -lm -o block
row_computing: row_computing.cpp
	g++ -g -O2 row_computing.cpp Source/lu_gp_sparse_row_computing.cpp Source/lu_gp_sparse.cpp -fpermissive Lib/libcxsparse.a Lib/libsuitesparseconfig.a -fPIC -lblas -lm -o row_computing
sn_computing: supernode_computing.cpp
	icpc -O3 -g -mavx -mavx2 supernode_computing.cpp Source/lu_gp_sparse_avx2.cpp Source/lu_gp_sparse_sn_u.cpp Source/lu_gp_sparse_supernode_computing.cpp Source/lu_gp_sparse.cpp Source/lu_gp_sparse_supernode_computing_blas.cpp Include/libcblas.a Include/liblapack.a Include/librefblas.a -fpermissive -fPIC -lblas -lm -lgfortran -o sn_computing
clean:
	rm -rf  sparse_lu sparse_lu.g row_computing supernode_computing sn_computing
