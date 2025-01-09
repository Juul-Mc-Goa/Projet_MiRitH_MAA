##
# Projet MiRitH Master Algèbre Appliquée
#
# @file
# @version 0.1

prog: main.o
	gcc -o prog main.o -lgmp -lcrypto

main.o: main.c
	gcc -c -Wall main.c

clean:
	rm -f prog *.o

# end
