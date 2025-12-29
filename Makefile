
all:	serial

serial:	argtable cpocketfft/src/libcpocketfft.a
	cd src && $(MAKE) serial

argtable:       | argtable2-13

argtable2-13: argtable2-13.tar.gz
	( \
            tar -xzf argtable2-13.tar.gz; \
            cd argtable2-13; \
            CFLAGS='--std=c90' ./configure --prefix=`pwd`; \
            make install; \
        )

argtable2-13.tar.gz:
	wget http://prdownloads.sourceforge.net/argtable/argtable2-13.tar.gz;

cpocketfft/src/libcpocketfft.a:
	cd cpocketfft/src && $(MAKE) libcpocketfft.a

.PHONY: clean test

test:
	cd test && $(MAKE)

clean:
	rm -rf argtable2-13
	cd cpocketfft/src && $(MAKE) clean
	cd src && $(MAKE) clean
	cd test && $(MAKE) clean
