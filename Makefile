PREFIX = /usr/local
NAME = EventGenerator

# ensure lib directory exists
$(shell mkdir -p lib)

.PHONY: all examples clean install

all: lib$(NAME).a
	make examples

lib$(NAME).a:
	make -C src
	ar crus lib/$@ `ls src/*.o`

install: lib$(NAME).a
	rsync -a lib/lib$(NAME).a $(PREFIX)/lib
	rsync -a include/*.hh $(PREFIX)/include

examples:
	make -C examples

clean:
	rm -fv lib/*
	make -C examples clean
	make -C src clean
