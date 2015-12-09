     # Diretorio base onde estarão os diretórios de biblioteca
    PREFIX = ./

     # Arquivo final
    FILE = pdeSolver

    CC = gcc
    AR = ar -rcu
    INSTALL = install
    FLAGS=-DLIKWID_PERFMON -O3 -mavx -march=native
#    LIKWID = -llikwid -L/usr/lib -L/home/soft/likwid/lib -DLIKWID_PERFMON

.PHONY: clean install all

%.o:  %.c
	$(CC) -c $(CFLAGS) $< -I/home/soft/likwid/include $(FLAGS)

all: install $(FILE).o

$(FILE).o: $(FILE).c
	$(CC) -static -o $(PREFIX)$(FILE) $(FILE).c -llikwid -L/usr/lib -L/home/soft/likwid/lib -lm $(FLAGS)

clean:
	@rm -f *% *.bak *~ *.o $(FILE) core *.swp
