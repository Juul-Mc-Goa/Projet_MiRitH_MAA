##
# Projet MiRitH Master Algèbre Appliquée
#
# @file
# @version 0.1

.PHONY: clean

test_dir := tests
all_tests := $(wildcard $(test_dir)/*.c)

export matrix_deps := matrix.o field_arithmetics.o constants.o
export key_gen_deps := key_generation.o random.o matrix.o field_arithmetics.o constants.o
export mpc_deps := mpc.o key_generation.o random.o matrix.o field_arithmetics.o constants.o

export key_gen_flags := -lgmp
export mpc_flags := -lgmp

objects := main.o key_generation.o random.o matrix.o field_arithmetics.o constants.o mpc.o

prog: $(objects)
	gcc -o $@ $^ -lgmp -lcrypto

$(objects): %.o: %.c
	gcc -c -Wall $^

test:
	cd $(test_dir) && $(MAKE)

clean:
	rm -f prog $(objects)

# end
