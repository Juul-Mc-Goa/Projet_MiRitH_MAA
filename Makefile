##
# Projet MiRitH Master Algèbre Appliquée
#
# @file
# @version 0.1

test_dir := tests
all_tests := $(wildcard ${test_dir}/*.c)

export matrix_deps := matrix.o field_arithmetics.o constants.o
export key_gen_deps := key_generation.o matrix.o field_arithmetics.o constants.o

export key_gen_flags := -lgmp

prog: main.o key_generation.o matrix.o field_arithmetics.o constants.o mpc.o
	gcc -o prog main.o \
key_generation.o \
field_arithmetics.o \
matrix.o \
constants.o \
mpc.o \
-lgmp -lcrypto

main.o: main.c
	gcc -c -Wall main.c

key_generation.o: key_generation.c
	gcc -c -Wall key_generation.c

mpc.o: mpc.c
	gcc -c -Wall mpc.c

matrix.o: matrix.c
	gcc -c -Wall matrix.c

field_arithmetics.o: field_arithmetics.c
	gcc -c -Wall field_arithmetics.c

constants.o: constants.c
	gcc -c -Wall constants.c

test:
	cd ${test_dir} && $(MAKE)

clean:
	rm -f prog *.o

# end
