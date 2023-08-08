EXEC=ecoPCR ecofind ecogrep

PCR_SRC= ecopcr.c
PCR_OBJ= $(patsubst %.c,%.o,$(PCR_SRC))

FIND_SRC= ecofind.c
FIND_OBJ= $(patsubst %.c,%.o,$(FIND_SRC))

GREP_SRC= ecogrep.c
GREP_OBJ= $(patsubst %.c,%.o,$(GREP_SRC))

IUT_SRC= ecoisundertaxon.c
IUT_OBJ= $(patsubst %.c,%.o,$(IUT_SRC))

SRCS= $(PCR_SRC) $(FIND_SRC) $(IUT_SRC)

LIB= -lecoPCR -lthermo -lapat -lz -lm

LIBFILE= libapat/libapat.a \
         libecoPCR/libecoPCR.a \
         libthermo/libthermo.a


include global.mk

all: $(EXEC)


########
#
# ecoPCR compilation
#
########

# executable compilation and link

ecoPCR: $(PCR_OBJ) $(LIBFILE) 
	$(CC) $(LDFLAGS) -o $@ $< $(LIBPATH) $(LIB)
	
########
#
# ecofind compilation
#
########
	
# executable compilation and link

ecofind: $(FIND_OBJ) $(LIBFILE)
	$(CC) $(LDFLAGS) -o $@ $< $(LIBPATH) $(LIB)
	
########
#
# ecogrep compilation
#
########
	
# executable compilation and link

ecogrep: $(GREP_OBJ) $(LIBFILE)
	$(CC) $(LDFLAGS) -o $@ $< $(LIBPATH) $(LIB)
	
########
#
# IsUnderTaxon compilation
#
########
	
# executable compilation and link

ecoisundertaxon: $(IUT_OBJ) $(LIBFILE)
	$(CC) $(LDFLAGS) -o $@ $< $(LIBPATH) $(LIB)

########
#
# library compilation
#
########

libapat/libapat.a:
	$(MAKE) -C libapat
	
libecoPCR/libecoPCR.a:
	$(MAKE) -C libecoPCR
	
libthermo/libthermo.a:
	$(MAKE) -C libthermo
	
	
########
#
# project management
#
########

clean:
	rm -f *.o	
	rm -f $(EXEC)
	$(MAKE) -C libapat clean
	$(MAKE) -C libecoPCR clean
	$(MAKE) -C libthermo clean
		

		