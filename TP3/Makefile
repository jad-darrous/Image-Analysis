#######################################
# Makefile Vision                     #  
#                                     #
#                                     #
#######################################
 

PROG = myprogram
all : $(PROG) 




# Variables pour la compilation des fichiers
CC        =  cc
CFLAGS    =  -g -Wall
CPPFLAGS  =  -DDEBUG 
LDFLAGS   =  -g 


# gmake
%.o : %.c utils.h pixmapIO.h
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $< -o $@ 

$(PROG): % : %.o utils.o pixmapIO.o
	$(CC) $(LDFLAGS) -o $@ $< utils.o pixmapIO.o 
#	$(CC) $(LDFLAGS) -o $@ $< utils.o pixmapIO.o  -lm -ll
 

clean:
	/bin/rm -f *.o
	/bin/rm -f myprogram