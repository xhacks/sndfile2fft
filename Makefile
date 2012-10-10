CC = gcc
CFLAGS= -std=gnu99
LIBS = -lfftw3 -lsndfile
INCLUDES = -I .

all: 
	     $(CC) $(CFLAGS) $(INCLUDES) $(LIBS) sndfile2fft.c -o sndfile2fft 

clean:
		rm -rf sndfile2fft
