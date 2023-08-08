
SOURCES = nnparams.c 

SRCS=$(SOURCES)
         
OBJECTS= $(patsubst %.c,%.o,$(SOURCES))

LIBFILE= libthermo.a
RANLIB= ranlib


include ../global.mk


all: $(LIBFILE)

clean:
	rm -rf $(OBJECTS) $(LIBFILE)

$(LIBFILE): $(OBJECTS)
	ar -cr $@ $?
	$(RANLIB) $@
