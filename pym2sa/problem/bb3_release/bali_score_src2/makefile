install: bali_score 

clean:
	rm *.o

HEADERS = general.h clustalw.h

CC	= cc
CFLAGS  = -c -O -I$(EXPAT_INC)
LFLAGS	= -O -lm -L$(EXPAT_LIB) -lexpat
EXPAT_LIB	= /usr/local/lib
EXPAT_INC	= /usr/local/include


bali_score : readxml.o init.o util.o bali_score.o
	$(CC) -o $@ readxml.o init.o util.o bali_score.o $(LFLAGS)
	cp bali_score /usr/local/bin/baliscore 

.c.o :	$(HEADERS)
	$(CC) $(CFLAGS) $?

