include ../makevars.in

SUBDIR = ewald steric spme lanczos

all:
	@for dir in $(SUBDIR); \
          do \
            make -C $${dir} || exit 1; \
          done
	@for dir in $(SUBDIR); \
          do \
            (cd $${dir} && LD_LIBRARY_PATH=$${LD_LIBRARY_PATH}:../../install/lib sh test.sh) || exit 1; \
          done

clean:
	@for dir in $(SUBDIR); \
          do \
            make -C $${dir} clean || exit 1; \
          done
