#######################################
# Makefile Vision                     #  
#                                     #
#                                     #
#######################################
 

PROG = main
all : $(PROG) 


# Variables pour la compilation des fichiers
CC        =  cc
CFLAGS    =  -g -Wall
CPPFLAGS  =  -DDEBUG 
LDFLAGS   =  -g 


# gmake
%.o : %.c utils.h pixmapIO.h kmeans.h
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $< -o $@ 

$(PROG): % : %.o utils.o pixmapIO.o kmeans.o
	$(CC) $(LDFLAGS) -o $@ $< utils.o pixmapIO.o kmeans.o
 

clean:
	/bin/rm -f *.o
	/bin/rm -f main
	/bin/rm -f main.exe