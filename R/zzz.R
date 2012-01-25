.First.lib <-function (lib, pkg) {
	library.dynam("iBMQ", pkg, lib)
}



#.First.lib<-function(lib, pkg){


#system("cp ./src/main_parallel_sparse.c main_parallel_sparse.c")
#system("cp ./src/sparse.c sparse.c")
#system("cp ./src/sparse.h sparse.h")
#system("cp ./src/norm_gamma_generation.c norm_gamma_generation.c")
#system("cp ./src/norm_gamma_generation.h norm_gamma_generation.h")
#system("cp ./src/RngStream.c RngStream.c")
#system("cp ./src/RngStream.h RngStream.h")
#system("cp ./src/Makevars Makevars")
#system("R CMD SHLIB main_parallel_sparse.c RngStream.c sparse.c norm_gamma_generation.c")
#dyn.load("main_parallel_sparse.so")
#library.dynam("main_parallel_sparse.so")
#}
