##
# MiRitH tests
#
# @file
# @version 0.1

local_matrix_deps := $(matrix_deps:%.o=../%.o) matrix_common.o
local_key_gen_deps := $(key_gen_deps:%.o=../%.o)
local_mpc_deps := $(mpc_deps:%.o=../%.o)

.PHONY: all clean

# Linked binaries
bin_with_mpc := share_alpha_a share_alpha
bin_with_key_gen := key_gen random_matrix
bin_with_matrix := matrix_sum matrix_prod
bins := $(bin_with_mpc) $(bin_with_key_gen) $(bin_with_matrix)

all: $(bins)

$(bin_with_mpc): %: %.o $(local_mpc_deps)
	gcc -o $@ $^ $$mpc_flags

$(bin_with_key_gen): %: %.o $(local_key_gen_deps)
	gcc -o $@ $^  $$key_gen_flags

$(bin_with_matrix): %: %.o $(local_matrix_deps)
	gcc -o $@ $^

# Unlinked binaries
objects := $(foreach bin,$(bins),$(bin).o) matrix_common.o

$(objects): %.o: %.c
	gcc -c -Wall $^

clean:
	rm -f $(objects) $(bins)
# end
