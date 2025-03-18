CC           = gcc
#CC          = icx  #for intel compiler
CCFLAGS      = -O3 -lgsl -lgslcblas -lm -pthread


gendynGssEALF : gendynGssEALF.c
	$(CC) gendynGssEALF.c $(CCFLAGS) -o gendynGssEALF

clean:; 	@rm -f gendynGssEALF


