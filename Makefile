    FONTES=$(wildcard *.c)
    OBJECTS=$(FONTES:.c=.o)

    # Arquivo final
    FILE = pdeSolver

    CC = gcc

    #FLAGS = -DLIKWID_PERFMON -O3 -mavx -march=native -lm
    # Sem Likwid
    #FLAGS = -O3 -lm -march=native

    # com Likwid
    FLAGS = -DLIKWID_PERFMON -O3 -lm -mavx -march=native

    INCLUDE = -I/home/soft/likwid/include
    LIKWID = -llikwid -L/usr/lib -L/home/soft/likwid/lib -DLIKWID_PERFMON

.PHONY: all clean

%.o:  %.c
	$(CC) -c -o $< -I/home/soft/likwid/include $(FLAGS)

all: pdeSolver

pdeSolver: $(OBJECTS)
	$(CC) -o $@ $^ $(LIKWID) $(FLAGS)

%.o: %.c
	$(CC) -c $< -o $@ $(INCLUDE) $(FLAGS)

clean:
	@rm -f *% *.bak *~ *.o $(FILE) core *.swp
