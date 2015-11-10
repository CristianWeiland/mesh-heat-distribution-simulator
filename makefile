     # Nome da biblioteca, usado para gerar o arquivo libnome.a
#    LIBNOME = 
#    LIBNOME2 = 

     # Nome do arquivo header com as declarações da biblioteca
#    INCFILES = 

     # Modulos que contém as funções da biblioteca
#    OBJECTS = 

     # Diretorio base onde estarão os diretórios de biblioteca
    PREFIX = ./

     # Arquivo de teste
    TEST = pdeSolver

    CC = gcc -g
    AR = ar -rcu
    INSTALL = install

.PHONY: clean distclean install all

#%.o:  %.c $(INCFILES)
%.o:  %.c
	$(CC) -c $(CFLAGS) $<

#all: install lib$(LIBNOME).a lib$(LIBNOME2).a $(TEST).o
all: install $(TEST).o

#lib$(LIBNOME).a: $(OBJECTS)
#	$(AR) $@ $?
#	ranlib $@

#lib$(LIBNOME2).a: $(OBJECTS)
#	$(AR) $@ $?
#	ranlib $@

$(TEST).o: $(TEST).c
	$(CC) -static -o $(PREFIX)$(TEST) $(TEST).c -lm
#	$(CC) -static -o $(PREFIX)$(TEST) -L -l$(LIBNOME) -L -l$(LIBNOME2) $(TEST).c


clean:
	@rm -f *% *.bak *~

distclean:   limpa
	@rm -rf *.o lib$(LIBNOME).* lib$(LIBNOME2).* $(TEST)
