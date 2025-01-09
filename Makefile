##
# Projet MiRitH Master Algèbre Appliquée
#
# @file
# @version 0.1

prog: main.o
	gcc -o prog main.o

main.o: main.c
	gcc -c -Wall main.c -lgmp

clean:
	rm -f prog *.o

# end
