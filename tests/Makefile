##
# MiRitH tests
#
# @file
# @version 0.1

local_matrix_deps := $(matrix_deps:%.o=../%.o)
local_key_gen_deps := $(key_gen_deps:%.o=../%.o)
local_mpc_deps := $(mpc_deps:%.o=../%.o)

# Linked binaries
all: share_alpha key_gen_m0 key_gen random_matrix matrix_sum matrix_prod

share_alpha: share_alpha.o $(local_mpc_deps)
	gcc -o $@ $^ $$mpc_flags

key_gen_m0: key_gen_m0.o $(local_key_gen_deps)
	gcc -o key_gen_m0 key_gen_m0.o $(local_key_gen_deps) $$key_gen_flags

key_gen: key_gen.o $(local_key_gen_deps)
	gcc -o key_gen key_gen.o $(local_key_gen_deps) $$key_gen_flags

random_matrix: random_matrix.o $(local_key_gen_deps)
	gcc -o random_matrix random_matrix.o $(local_key_gen_deps) $$key_gen_flags

matrix_sum: matrix_sum.o matrix_common.o $(local_matrix_deps)
	gcc -o matrix_sum matrix_sum.o matrix_common.o $(local_matrix_deps)

matrix_prod: matrix_prod.o matrix_common.o
	gcc -o matrix_prod matrix_prod.o matrix_common.o $(local_matrix_deps)


# Unlinked binaries

key_gen_m0.o: key_gen_m0.c
	gcc -c -Wall key_gen_m0.c

matrix_common.o: matrix_common.c
	gcc -c -Wall matrix_common.c

matrix_sum.o: matrix_sum.c
	gcc -c -Wall matrix_sum.c

matrix_prod.o: matrix_prod.c
	gcc -c -Wall matrix_prod.c

random_matrix.o: random_matrix.c
	gcc -c -Wall random_matrix.c

# end
