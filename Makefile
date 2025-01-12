##
# Projet MiRitH Master Algèbre Appliquée
#
# @file
# @version 0.1

prog: main.o key_generation.o matrix.o
	gcc -o prog main.o key_generation.o matrix.o -lgmp -lcrypto

main.o: main.c
	gcc -c -Wall main.c

key_generation.o: key_generation.c
	gcc -c -Wall key_generation.c

matrix.o: matrix.c
	gcc -c -Wall matrix.c

clean:
	rm -f prog *.o

# end
