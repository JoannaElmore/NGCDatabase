LIBPATH= -Llibapat -LlibecoPCR -Llibthermo
MAKEDEPEND = gcc -M $(CPPFLAGS) -o $*.d $<

CC=gcc
CFLAGS= -O3 -w

default: all

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

%.P : %.c
	$(MAKEDEPEND)
	@sed 's/\($*\)\.o[ :]*/\1.o $@ : /g' < $*.d > $@; \
	rm -f $*.d; [ -s $@ ] || rm -f $@

include $(SRCS:.c=.P)