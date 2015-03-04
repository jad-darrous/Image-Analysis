#######################################
# Makefile PBM                        #
#                                     #
# E.B.                                #
#######################################

# Variables for file compilation
CC        =  gcc
CFLAGS    =  -g -Wall
CPPFLAGS  =  -DDEBUG
LDFLAGS   =  -g -lm

main: main.c Util.h Util.c ImgUtils.h ImgUtils.c 
	$(CC) $(CFLAGS) $^ -o $@

clean :
	@rm -f *.o main.exe

