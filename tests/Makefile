##
# MiRitH tests
#
# @file
# @version 0.1

.PHONY: all clean

# inherit dependencies exported from `../Makefile`
local_matrix_deps := $(matrix_deps:%.o=../%.o) matrix_common.o
local_key_gen_deps := $(key_gen_deps:%.o=../%.o)
local_mpc_deps := $(mpc_deps:%.o=../%.o) mpc_common.o

# Linked binaries
bin_with_matrix := matrix_sum matrix_prod
bin_with_key_gen := random_matrix key_gen
bin_with_mpc := share_alpha share_alpha_a share_alpha_a_c_k
bins := $(bin_with_matrix) $(bin_with_key_gen) $(bin_with_mpc)

# Run all tests
all: $(bins)
	$(foreach bin,$(bins),./$(bin) && ) echo "All tests done."

$(bin_with_mpc): %: %.o $(local_mpc_deps)
	gcc -o $@ $^ $$mpc_flags

$(bin_with_key_gen): %: %.o $(local_key_gen_deps)
	gcc -o $@ $^  $$key_gen_flags

$(bin_with_matrix): %: %.o $(local_matrix_deps)
	gcc -o $@ $^

# Objects
objects := $(foreach bin,$(bins),$(bin).o) matrix_common.o mpc_common.o

$(objects): %.o: %.c
	gcc -c -Wall $^

# Other
clean:
	rm -f $(objects) $(bins)

# end
