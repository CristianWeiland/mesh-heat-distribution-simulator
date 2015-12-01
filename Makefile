     # Diretorio base onde estarão os diretórios de biblioteca
    PREFIX = ./

     # Arquivo final
    FILE = pdeSolver

    CC = gcc -g
    AR = ar -rcu
    INSTALL = install
    OPTIMIZATION = -O2

.PHONY: clean install all

%.o:  %.c
	$(CC) -c $(CFLAGS) $<

all: install $(FILE).o

$(FILE).o: $(FILE).c
	$(CC) -static -o $(PREFIX)$(FILE) $(OPTIMIZATION) $(FILE).c -lm

clean:
	@rm -f *% *.bak *~ *.o $(FILE) core *.swp