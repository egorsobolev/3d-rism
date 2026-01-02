CC=gcc

FFTFLAGS= ../cpocketfft/src/libcpocketfft.a -lstdc++

SRC = $(wildcard *.c)
OBJSRL = $(addprefix .o.serial/,$(SRC:.c=.o))

all:	serial

serial:	$(OBJSRL) | bin
	$(CC) -o ../bin/rism3d $(OBJSRL) \
	    ../argtable2-13/lib/libargtable2.a \
            $(FFTFLAGS) -lm -lopenblas -fopenmp

.o.serial/%.o:	%.c
	$(CC) -c -O3 -fopenmp -I../argtable2-13/include -I../cpocketfft/src -o $@ $<

$(OBJSRL):	| .o.serial

.o.serial:
	mkdir .o.serial

bin:
	mkdir -p ../bin

.PHONY:	clean
clean:
	-rm -rf ../bin
	-rm -rf .o.serial
