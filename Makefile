##
# Projet MiRitH Master Algèbre Appliquée
#
# @file
# @version 0.1

prog: main.o key_generation.o matrix.o field_arithmetics.o
	gcc -o prog main.o key_generation.o field_arithmetics.o matrix.o -lgmp -lcrypto

main.o: main.c
	gcc -c -Wall main.c

key_generation.o: key_generation.c
	gcc -c -Wall key_generation.c

field_arithmetics.o: field_arithmetics.c
	gcc -c -Wall field_arithmetics.c

matrix.o: matrix.c
	gcc -c -Wall matrix.c

clean:
	rm -f prog *.o

# end
