.PHONY: all install clean app tests

all: app

lib:
	make -C src install

app: lib
	make -C bin install

test: lib
	make -C tests
    
tool:
	make -C tools install
    
matlab:
	make -C mex install
	
clean:
	make -C src clean
	make -C bin clean
	make -C tests clean
	make -C tools clean
	make -C mex clean
