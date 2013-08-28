CC = gcc
CFLAGS= -std=gnu99
LDFLAGS = -L/opt/local/lib
LIBS = -lfftw3 -lsndfile
INCLUDES = -I/opt/local/include

all: 
	     $(CC) $(CFLAGS) $(LDFLAGS) $(INCLUDES) $(LIBS) sndfile2fft.c -o sndfile2fft 

clean:
		rm -rf sndfile2fft
