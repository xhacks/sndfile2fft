/*
** Copyright (C) 2008-2011 Erik de Castro Lopo <erikd@mega-nerd.com>
**
** All rights reserved.
**
** Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**
**     * Redistributions of source code must retain the above copyright
**       notice, this list of conditions and the following disclaimer.
**     * Redistributions in binary form must reproduce the above copyright
**       notice, this list of conditions and the following disclaimer in
**       the documentation and/or other materials provided with the
**       distribution.
**     * Neither the author nor the names of any contributors may be used
**       to endorse or promote products derived from this software without
**       specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
** TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
** PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
** CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
** EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
** PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
** OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
** WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
** OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
** ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include <sndfile.h>
#include "fftw3.h"
char 		*progname, *infilename, *outfilename ;

// Find enerygy in peak - save +/- 1 bin.
// find energy out of peak.

// TODO run with larges and small window?

#define MAX_BLOCK_SIZE  102400
#define MAX_CHANS       10
#define BLOCK_SIZE      (1024*8)
//#define BLOCK_SIZE 128
double g_div = 0;

static double getwindowcoef(int i)
{
    double mult = 1;
    double PI = 3.1415926535;
#if 1
    // Blackman-harris
    double a0 = 0.35875;
    double a1 = 0.48829;
    double a2 = 0.14128;
    double a3 = 0.01168;
    mult = a0 - a1*cos(2*PI*i/(BLOCK_SIZE-1)) + a2 * cos(4*PI*i/(BLOCK_SIZE-1)) - a3*cos(6*PI*i/(BLOCK_SIZE-1));

#endif
#if 0
            //Hann
            double x = sin(PI * i / (BLOCK_SIZE-1));
            double hann = x*x;
            mult = hann;
#endif
#if 0        
            // Blackman
            double a0 = 7938.0/18608.0;
            double a1 = 9240.0/18608.0;
            double a2 = 1430.0/18608.0;


            double black = a0 -(a1*cos(2*PI*i/(BLOCK_SIZE-1))) + (a2*cos(4*PI*i)/(BLOCK_SIZE-1));
            mult = black;
#endif
#if 0
            double hamming = 0.54 -0.46*cos(2*PI*i/(BLOCK_SIZE-1));
            mult = hamming;

#endif

    return mult;

}


fftw_complex out[MAX_BLOCK_SIZE];
double input[MAX_BLOCK_SIZE];
double mag[MAX_BLOCK_SIZE];
double peakhold[MAX_CHANS][MAX_BLOCK_SIZE];

double otherEnergy[MAX_CHANS];
double peakEnergy[MAX_CHANS];

int g_peakIndex[MAX_CHANS];

int g_channels = 0;
int g_samplerate = 0;
int g_expectedFreq = 0;

int g_blockCount = 0;

// N set to sample rate for 1 second blocks

static void
print_usage (char *progname)
{	fprintf (stderr, "\nUsage : %s <input file> <output file> <expected freq>\n", progname) ;
	puts ("\n"
		"    Where the output file will contain a line for each frame\n"
		"    and a column for each channel.\n"
		) ;

} /* print_usage */


int properfft(int nreal, int channel)
{
    fftw_plan p;
    int i;
    double max = 0;
    double max2 = 0;
    int maxIndex = 0;
 
    g_blockCount++;

    p = fftw_plan_dft_r2c_1d(BLOCK_SIZE, input, out, FFTW_ESTIMATE);
    
    fftw_execute(p);
    
    fftw_destroy_plan(p);
    for(i=0;i<nreal;i++)
    {
        mag[i] = hypot(out[i][0], out[i][1]);
        if (mag[i] > max2 && i != 0) 
        {
            max2 = mag[i];
        }
    }
    
    for(i=0;i<nreal;i++)
    {
        mag[i] = 20*log10(mag[i]/max2+1e-10);

        if(mag[i] > max)
        {
            max = mag[i];
            maxIndex = i;
        }
        
        if(mag[i]>peakhold[channel][i])
        {
            peakhold[channel][i] = mag[i];
        }
    }
    
    //fprintf(stderr, "i: %d\n", maxIndex);
    return maxIndex;
}

static int doBlock(float buf[], int readcount, int lastIndex)
{

    int i = 0;
    double mult = 0;
  //  int m = 0;
    double val;
    int index;

    for(int c = 0; c < g_channels; c++)
    {
        i = 0;
        for (int k = 0 ; k < readcount ; k++)
	    {
            mult = getwindowcoef(i);

            val = buf[k*g_channels+c];
            input[i] = val*mult;// /(float)((unsigned)(1<<31)); // real
            
            i++;
            if (i >= BLOCK_SIZE) break;
        } ;
        
        /* Do the fft over the block */
        g_peakIndex[c] = properfft(i, c);

#if 0 
        if(lastIndex != -1)
        {
            if(lastIndex!=index)
            {
                //fprintf(stderr, "PEAK MOVED\n");
            //exit(1);
            }
        }
#endif
    }
    
    /* Return index of peak */
    return index;
}



static void
convert_to_text (SNDFILE * infile, FILE * outfile, int channels)
{	
    float buf1 [channels * BLOCK_SIZE] ;
	float buf2 [channels * BLOCK_SIZE] ;
	float block [channels * BLOCK_SIZE] ;
	int k, m, readcount ;
    
    double edge, val, time = 0;
    int index = -1;
    int seconds = 0;
    int firstRun = 1;
    double mult;


    // Init peak hold
    for(int i = 0; i < BLOCK_SIZE; i++)
        for(int j = 0; j < g_channels; j++)
        {
            peakhold[j][i] = -100000;
            otherEnergy[j] = 0.0;
            peakEnergy[j] = 0.0;
            
        }
    /* Ignore first few blocks - hack to avoid 0's in start of wav */
    while(1)
    {
        if((readcount = sf_readf_float (infile, buf1, 1)) <= 0)
        {
            exit(1);
        }

        if(buf1[0] != 0.0)
            break;
    }   
    
    if((readcount = sf_readf_float (infile, buf1, BLOCK_SIZE)) <= 0)
        exit(1);    
	
    while (1)
    {
        readcount = sf_readf_float (infile, buf2, BLOCK_SIZE);
        if(readcount < 0)
        {
            break;
        }
        else if(readcount != BLOCK_SIZE)
        {
            fprintf(stderr, "Info: Dropped last block\n");
            break;
        }

        for(int i = 0; i < (BLOCK_SIZE); i++)
        {
            block[i] = buf1[i+(BLOCK_SIZE)];
        } 
        for(int i = 0; i < (BLOCK_SIZE); i++)
        {
            block[i+(BLOCK_SIZE)] = buf2[i];
        }

        index = doBlock(block, readcount, index);
        index = doBlock(buf2, readcount, index);
        
        readcount = sf_readf_float(infile, buf1, BLOCK_SIZE);
        if(readcount < 0)
        {
            break;
        }
        else if(readcount != BLOCK_SIZE)
        {
            fprintf(stderr, "Info: Dropped last block\n");
            break;
        }

        for(int i = 0; i < (BLOCK_SIZE); i++)
        {
            block[i] = buf2[i+(BLOCK_SIZE)];
        } 
        for(int i = 0; i < (BLOCK_SIZE); i++)
        {
            block[i+(BLOCK_SIZE)] = buf1[i];
        }

        index = doBlock(block, readcount, index);
        index = doBlock(buf1, readcount, index);

    } ;
 
   FILE *fftout = fopen ("fftout.txt", "w");

    /* Print out the output */
    for(int i = 0; i < BLOCK_SIZE; i++)
    {
        fprintf(fftout, "%10d %10g %10g ", (int)(i*g_div), mag[i], input[i]);
        for(int c = 0; c<g_channels; c++)
        {
            fprintf(fftout, "%10g ", peakhold[c][i]);
        }
        fprintf(fftout, "\n");
    }

#define PEAK_RANGE 5
    for(int c = 0; c < g_channels; c++)
    { 
        int othercount = 0;  
        int peakIndex = g_peakIndex[c];
       
        for(int i = (-PEAK_RANGE); i< PEAK_RANGE;i++)
        {
            peakEnergy[c] += peakhold[c][g_peakIndex[c]+i];
        }
        
        for(int i = 0; i < BLOCK_SIZE; i++)
        {
            if((i < (peakIndex-PEAK_RANGE)) || (i > (peakIndex+PEAK_RANGE)))
            {
                if(peakhold[c][i] > -200)
                {
                    otherEnergy[c] -= peakhold[c][i];
                    othercount++;
                }
            } 
        }

        peakEnergy[c] /= (PEAK_RANGE*2+1);
        otherEnergy[c] /= (othercount-(PEAK_RANGE*2+1));
        otherEnergy[c] *= -1;

        /* Convert to freq */
        g_peakIndex[c] = (int)((double)g_peakIndex[c]*g_div); 
 
        fprintf(stderr, "Channel %d : Peak found at approx %d Hz, Peak avg: %10g, Other Freqs: %10g\n", c, g_peakIndex[c], (peakEnergy[c]), otherEnergy[c]);

        if(otherEnergy[c] > -85)
        {
            fprintf(stderr, "Suspect glitch on this channel!\n\n");
        }
        else
        {
            fprintf(stderr, "Looks okay - manual inspect to be sure!\n\n");
        }

    }
    fprintf(stderr, "\n");

    FILE *pipe = popen("gnuplot ","w");
    fprintf(pipe, "set xtic auto; set ytic auto; set title 'Spectrum (%s)'; set ylabel 'Relative Amp (dB)'; \
        set xlabel 'Freq (Hz)'; set xr[0:20000]; set yr [-140:1];",infilename );

// set xr [0:25000]

    fprintf(pipe, "set arrow from %d,-140 to %d,-120 lc rgb 'red';", g_peakIndex[0], g_peakIndex[0]);
    fprintf(pipe, "set label '%d Hz' at %d,-138;", g_peakIndex[0], g_peakIndex[0]);

    fprintf(pipe, "set arrow from %d,-140 to %d,-120 lc rgb 'navy';", g_peakIndex[1], g_peakIndex[1]);
    fprintf(pipe, "set label '%d Hz' at %d,-138;", g_peakIndex[1], g_peakIndex[1]);

    fprintf(pipe, "set arrow from 1,%d to 20000, %d lc rgb 'red' nohead ;", (int)otherEnergy[0], (int)otherEnergy[0]);
    fprintf(pipe, "set arrow from 1,%d to 20000, %d lc rgb 'navy' nohead ;", (int)otherEnergy[1], (int)otherEnergy[1]);

    fprintf(pipe, "plot 'fftout.txt' using 1:4 title 'Left (peak hold)' with lines lc rgb 'red', \
            'fftout.txt' using 1:5 title 'Right (peak hold)' with lines lc rgb 'navy'");

    close(pipe);     

    return ;
} /* convert_to_text */

int
main (int argc, char * argv [])
{		SNDFILE	 	*infile = NULL ;
	FILE		*outfile = NULL ;
	SF_INFO	 	sfinfo ;

	progname = strrchr (argv [0], '/') ;
	progname = progname ? progname + 1 : argv [0] ;

	if (argc != 3)
	{	print_usage (progname) ;
		return 1 ;
		} ;

	infilename = argv [1] ;
	outfilename = argv [2] ;


	if (strcmp (infilename, outfilename) == 0)
	{	printf ("Error : Input and output filenames are the same.\n\n") ;
		print_usage (progname) ;
		return 1 ;
		} ;

	if (infilename [0] == '-')
	{	printf ("Error : Input filename (%s) looks like an option.\n\n", infilename) ;
		print_usage (progname) ;
		return 1 ;
		} ;

	if (outfilename [0] == '-')
	{	printf ("Error : Output filename (%s) looks like an option.\n\n", outfilename) ;
		print_usage (progname) ;
		return 1 ;
		} ;

	if ((infile = sf_open (infilename, SFM_READ, &sfinfo)) == NULL)
	{	printf ("Not able to open input file %s.\n", infilename) ;
		puts (sf_strerror (NULL)) ;
		return 1 ;
		} ;

	/* Open the output file. */
	if ((outfile = fopen (outfilename, "w")) == NULL)
	{	printf ("Not able to open output file %s : %s\n", outfilename, sf_strerror (NULL)) ;
		return 1 ;
		} ;

	fprintf (outfile, "# Converted from file %s.\n", infilename) ;
	fprintf (outfile, "# Channels %d, Sample rate %d\n", sfinfo.channels, sfinfo.samplerate) ;

	fprintf (stderr, "Running from input file %s.\n", infilename) ;
	fprintf (stderr, "Channels %d, Sample rate %d\n", sfinfo.channels, sfinfo.samplerate) ;
 
    g_samplerate = sfinfo.samplerate;
    g_channels = sfinfo.channels; 

    //BLOCK_SIZE = g_samplerate/DIV; 
    g_div = (double)g_samplerate/(double)(BLOCK_SIZE);
     

	convert_to_text (infile, outfile, sfinfo.channels) ;

	sf_close (infile) ;
	fclose (outfile) ;

	return 0 ;
} /* main */

