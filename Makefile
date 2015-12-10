     # Diretorio base onde estarão os diretórios de biblioteca
    PREFIX = ./

    FONTES=$(wildcard *.c)
    OBJECTS=$(FONTES:.c=.o)

     # Arquivo final
    FILE = pdeSolver

    CC = gcc
    #FLAGS = -DLIKWID_PERFMON -O3 -mavx -march=native -lm
    FLAGS = -DLIKWID_PERFMON -O3 -lm -march=native
    INCLUDE = -I/home/soft/likwid/include
    LIKWID = -llikwid -L/usr/lib -L/home/soft/likwid/lib -DLIKWID_PERFMON

.PHONY: all clean

%.o:  %.c
	$(CC) -c -o $< -I/home/soft/likwid/include $(FLAGS)

all: pdeSolver

$(FILE): $(OBJECTS)
	$(CC) -o $@ $^ $(LIKWID) $(FLAGS)

%.o: %.c
	$(CC) -c $< -o $@ $(INCLUDE) $(FLAGS)

clean:
	@rm -f *% *.bak *~ *.o $(FILE) core *.swp
