sndfile2fft
===========

This is a simple application I hacked together while testing USB Audio Interfaces.

It takes a wav file, runs an FFT (using FFTW), plots the output (using gnuplot), just some simple checks for glitches (expects sine waves)

Requirements
-------------
  * Gnuplot for generating plots
  * libsndfie (Version used 1.0.25) - http://www.mega-nerd.com/libsndfile/#Download
    * (I came accross this issue when building from source: https://trac.macports.org/ticket/35358)
  * fftw3 (build from source or use macports)

Building
---------
Something like the following should do it..

gcc -lfftw3 sndfile2fft.c -lsndfile -std=gnu99 -o sndfile2fft

Usage
------

sndfile2fft <input file> <output file>

Where the output file will contain a line for each freq bin. First column freq value, then column for each channel.

Testing
--------
Only tested on OSX and a few recorded wav files (some with injected errors).


TODO
----
  * Only for stereo really
  * Could do with command line options for the various options
  * Needs tidy up - Very hacky! 
  * Check detected peaks dont move
