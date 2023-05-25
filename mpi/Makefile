CFLAGS=-Wextra
NUM_OF_PROC=4

run:build
	mpiexec -n $(NUM_OF_PROC) mpiq

build:mpiq.c
	mpicc mpiq.c -o mpiq $(CFLAGS)
