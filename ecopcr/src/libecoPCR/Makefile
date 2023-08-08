
SOURCES = ecoapat.c \
         ecodna.c \
         ecoError.c \
         ecoIOUtils.c \
         ecoMalloc.c \
         ecorank.c \
         ecoseq.c \
         ecotax.c \
         ecofilter.c \
         econame.c     

SRCS=$(SOURCES)
         
OBJECTS= $(patsubst %.c,%.o,$(SOURCES))

LIBFILE= libecoPCR.a
RANLIB= ranlib


include ../global.mk


all: $(LIBFILE)

clean:
	rm -rf $(OBJECTS) $(LIBFILE)

$(LIBFILE): $(OBJECTS)
	ar -cr $@ $?
	$(RANLIB) $@
