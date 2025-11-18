CC=gcc

FFTWLFLAGS= -lfftw3 -lfftw3f

SRC = $(wildcard *.c)
OBJSRL = $(addprefix .o.serial/,$(SRC:.c=.o))

all:	serial

serial:	$(OBJSRL) | bin
	$(CC) -o ../bin/rism3d $(OBJSRL) \
	    ../argtable2-13/lib/libargtable2.a \
            $(FFTWLFLAGS) -lm -lopenblas

.o.serial/%.o:	%.c
	$(CC) -c -I../argtable2-13/include -o $@ $<

$(OBJSRL):	| .o.serial

.o.serial:
	mkdir .o.serial

bin:
	mkdir -p ../bin

.PHONY:	clean
clean:
	-rm -rf ../bin
	-rm -rf .o.serial
