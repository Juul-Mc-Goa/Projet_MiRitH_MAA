##
# Projet MiRitH Master Algèbre Appliquée
#
# @file
# @version 0.1

.PHONY: clean

test_dir := tests
all_tests := $(wildcard $(test_dir)/*.c)

# define dependencies for various files, used in `tests/Makefile`
field_deps := field_arithmetics.o constants.o
export matrix_deps := matrix.o field_arithmetics.o constants.o
random_matrix_deps := random.o $(matrix_deps)
export key_gen_deps := key_generation.o $(random_matrix_deps)
export mpc_deps := packing.o mpc.o $(random_matrix_deps)
export all_deps := seed_tree.o packing.o mpc.o key_generation.o $(random_matrix_deps)

export key_gen_flags := -lgmp
export mpc_flags := -lgmp
export all_flags := -lgmp -lcrypto -lm

objects := main.o sign.o verif.o key_generation.o random.o matrix.o \
field_arithmetics.o constants.o mpc.o seed_tree.o packing.o

# Linked binary
prog: $(objects)
	gcc -o $@ $^ $(all_flags)

# Objects
$(objects): %.o: %.c
	gcc -c -Wall $^

# Other
test:
	cd $(test_dir) && $(MAKE)

clean:
	rm -f prog $(objects) && cd $(test_dir) && $(MAKE) clean

# end
