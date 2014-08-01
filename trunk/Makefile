.PHONY: all install clean app tests

all: app

lib:
	make -C src install

app: lib
	make -C bin
	make -C bin install

test: lib
	make -C tests
	
clean:
	make -C src clean
	make -C bin clean
	make -C tests clean
