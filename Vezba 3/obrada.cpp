#include <math.h>
#include <stdlib.h>
#include "obrada.h"
#include "sr_fft.h"

double time_buffer[FFT_SIZE];
double fft_buffer[FFT_SIZE];
double in_delay[FFT_SIZE/2];
double out_delay[FFT_SIZE/2];

extern double window[FFT_SIZE];

#define PI 3.14

void obrada(double *in, double *out, int N)
{
  int i;
  int k1 = 1024 * 300 / 44100;
  int k2 = 1024 * 4000 / 44100; 

  for(i=0; i<N; i++)
  {
	  time_buffer[i] = in_delay[i] * window[i];
	  time_buffer[N+i] = in[i] * window[N+i];
	  in_delay[i] = in[i];
  }

  fft(time_buffer, fft_buffer, FFT_ORDER);
  
  /*for(i=0; i<2*N; i+=2)
  {
	  if(i < 2*k1 || i > 2*k2 )
	  {
		  fft_buffer[i] = 0;
		  fft_buffer[i + 1] = 0;
	  }
  }*/

  for(i=1; i<2*N; i+=2)
  {
	  double tmp = sqrt(pow(fft_buffer[i],2) + pow(fft_buffer[i+1],2));
	  fft_buffer[i] = tmp * cos(0.25 * PI * (0.25*i*i));
	  fft_buffer[i + 1] = tmp * sin(0.25 * PI * (0.25*i*i));

  }

  ifft(fft_buffer, time_buffer, FFT_ORDER);

  for(i=0; i<N; i++)
  {

	  out[i] = out_delay[i] * window[N+i] + time_buffer[i] * window[i];
	  out_delay[i] = time_buffer[N+i];
  }

}
