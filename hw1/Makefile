FLAGS=O0 O2 O3 
GCC_FLAGS=$(addprefix -,$(FLAGS))
ADD_TEST_SIZE=3000 10000 30000 100000
.PHONY=add test-add vec
#TODO one by one
add:
	gcc -O0 add.c -o add_O0
	gcc -O2 add.c -o add_O2
	gcc -O3 add.c -o add_O3

test-add:add
	add_O0 3000
	add_O0 10000
	add_O0 30000
	add_O0 100000
	add_O2 3000
	add_O2 10000
	add_O2 30000
	add_O2 100000

vec:
	gcc -O0 vector.c -o vector_O0
	gcc -O2 vector.c -o vector_O2
	gcc -O3 vector.c -o vector_O3