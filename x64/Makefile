#To compile with gcc: "make runtests gcc=1"
CLANGFLAGS = -O3 -march=native -DMAIN
GCCFLAGS = -O3 -march=native -mtune=native -DMAIN

all:
ifdef gcc
	gcc -o run_gcc $(GCCFLAGS) dh.c
	./run_gcc
else
	clang -o run_clang $(CLANGFLAGS) dh.c
	./run_clang
endif
