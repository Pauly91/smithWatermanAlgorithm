all:
	gcc smithWaterManCPU.c -Wall -O3 -g -o dna
clean:
	rm dna