# Makefile, versao 02
# Sistemas Operativos, DEI/IST/ULisboa 2017-18

CFLAGS= -g -Wall -pedantic -std=gnu99
CC=gcc
PTHREADS = -pthread

heatSim: main.o matrix2d.o 
	$(CC) $(CFLAGS) -o heatSim main.o matrix2d.o $(PTHREADS) 

main.o: main.c matrix2d.h 
	$(CC) $(CFLAGS) -c main.c

matrix2d.o: matrix2d.c matrix2d.h
	$(CC) $(CFLAGS) -c matrix2d.c

clean:
	rm -f *.o heatSim

zip:
	zip  3EntregaProjetoSOGrupo100.zip main.c matrix2d.c matrix2d.h readme.txt Makefile

run:
	./heatSim 5 10.0 10.0 0.0 0.0 10 5 0.3

