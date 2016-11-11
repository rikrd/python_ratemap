/*
 *=========================================================================
 * An efficient C implementation of MPC's makeRateMap matlab code
 *-------------------------------------------------------------------------
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *-------------------------------------------------------------------------
%
%  ratemap = ratemap(x,fs,lowcf,highcf,numchans,frameshift,ti,compression)
% 
%  x           input signal
%  fs          sampling frequency in Hz (8000)
%  lowcf       centre frequency of lowest filter in Hz (50)
%  highcf      centre frequency of highest filter in Hz (3500)
%  numchans    number of channels in filterbank (32)
%  frameshift  interval between successive frames in ms (10)
%  ti          temporal integration in ms (8)
%  compression type of compression ['cuberoot','log','none'] ('cuberoot')
%
%  e.g. ratemap = ratemap(x,8000,50,3850,32,10);
%
%  For more detail on this implementation, see
%  http://www.dcs.shef.ac.uk/~ning/resources/ratemap/
%
%  Ning Ma, University of Sheffield
%  n.ma@dcs.shef.ac.uk, 08 Dec 2005
%
 * $ JPB: 7 August 2014 $
 *  pure C implementation stripped from Ning's MATLAB mex code
 *  for integration with Python. See ratemap.py and test_ratemap.py
 *  compile using: cc -fPIC -shared -o libratemap.so ratemap.c
 *
 * $Test: 22 Dec 2005 $
 *  Tested using a Linux PC Pentium4 2.0G Hz. This implementation
 *  takes 0.29 secs to compute a 64-channel ratemap for a 
 *  signal with a duration of 3.7654 secs and a sampling rate 
 *  of 20K. Matlab code takes 9.56 secs for the same signal.
 *
 *=============================================================
 */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

/*=======================
 * Useful Const
 *=======================
 */
#define BW_CORRECTION   1.019
#define VERY_SMALL_NUMBER  1e-200

#ifndef M_PI
#define M_PI            3.14159265358979323846
#endif

/*=======================
 * Utility functions
 *=======================
 */
#define getMax(x,y)     ((x)>(y)?(x):(y))
#define getRound(x)     ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

#define erb(x)          (24.7*(4.37e-3*(x)+1.0))
#define HzToErbRate(x)  (21.4*log10(4.37e-3*(x)+1.0))
#define ErbRateToHz(x)  ((pow(10.0,((x)/21.4))-1.0)/4.37e-3)

/*=======================
 * Main Function
 *=======================
 */


void ratemap(double *x, int nsamples, int fs, double lowcf, double highcf, int numchans, double frameshift, double ti, char* compression, double* result)
{
  double *senv;
  int i, j, chan;
  int nsamples_padded, frameshift_samples, numframes;
  double lowErb, highErb, spaceErb, cf;
  double a, tpt, tptbw, gain, intdecay, intgain, sumEnv;
  double p0r, p1r, p2r, p3r, p4r, p0i, p1i, p2i, p3i, p4i;
  double a1, a2, a3, a4, a5, cs, sn, u0r, u0i;
  double senv1, oldcs, oldsn, coscf, sincf;

  
  /*=========================================
   * input arguments
   *=========================================
   */
  
  frameshift_samples = getRound(frameshift*fs/1000);
  numframes = (int)ceil((double)nsamples / (double)frameshift_samples);
  nsamples_padded = numframes*frameshift_samples;

  /*=========================================
   * output arguments
   *=========================================
   */
  
  lowErb = HzToErbRate(lowcf);
  highErb = HzToErbRate(highcf);
  
  if (numchans > 1)  spaceErb = (highErb-lowErb)/(numchans-1);
  else  spaceErb = 0.0;

  /* Allocate memory for ratemap and smoothed envelope */
  senv = (double*) malloc(nsamples_padded * sizeof(double));

  tpt = 2 * M_PI / fs;
  intdecay = exp(-(1000.0/(fs*ti)));
  intgain = 1-intdecay;

  for (chan=0; chan<numchans; chan++)
  {
    cf = ErbRateToHz(lowErb+chan*spaceErb);
    tptbw = tpt * erb(cf) * BW_CORRECTION;
    a = exp(-tptbw);
    gain = (tptbw*tptbw*tptbw*tptbw)/3;

    /* Update filter coefficiants */
    a1 = 4.0*a; a2 = -6.0*a*a; a3 = 4.0*a*a*a; a4 = -a*a*a*a; a5 = a*a;

    p0r=0.0; p1r=0.0; p2r=0.0; p3r=0.0; p4r=0.0;
    p0i=0.0; p1i=0.0; p2i=0.0; p3i=0.0; p4i=0.0;
    senv1=0.0;

    /*====================================================================================
     * complex z=x+j*y, exp(z) = exp(x)*(cos(y)+j*sin(y)) = exp(x)*cos(x)+j*exp(x)*sin(y).
     * z = -j * tpti * cf, exp(z) = cos(tpti*cf) - j * sin(tpti*cf)
     *====================================================================================
     */
    coscf = cos(tpt*cf);
    sincf = sin(tpt*cf);
    cs = 1; sn = 0;
    for (i=0; i<nsamples; i++)
    {      
      p0r = cs*x[i] + a1*p1r + a2*p2r + a3*p3r + a4*p4r;
      p0i = sn*x[i] + a1*p1i + a2*p2i + a3*p3i + a4*p4i;
     
      /* Clip coefficients to stop them from becoming too close to zero */
      if (fabs(p0r) < VERY_SMALL_NUMBER)
        p0r = 0.0F;
      if (fabs(p0i) < VERY_SMALL_NUMBER)
        p0i = 0.0F;

      u0r = p0r + a1*p1r + a5*p2r;
      u0i = p0i + a1*p1i + a5*p2i;

      p4r = p3r; p3r = p2r; p2r = p1r; p1r = p0r;
      p4i = p3i; p3i = p2i; p2i = p1i; p1i = p0i;

     /*==========================================
      * Smoothed envelope by temporal integration 
      *==========================================
      */
      senv1 = senv[i] = sqrt(u0r*u0r+u0i*u0i) * gain + intdecay*senv1;

     /*=========================================
      * cs = cos(tpt*i*cf);
      * sn = -sin(tpt*i*cf);
      *=========================================
      */
      cs = (oldcs=cs)*coscf + (oldsn=sn)*sincf;
      sn = oldsn*coscf - oldcs*sincf;
    }
    for (i=nsamples; i<nsamples_padded; i++)
    {
      p0r = a1*p1r + a2*p2r + a3*p3r + a4*p4r;
      p0i = a1*p1i + a2*p2i + a3*p3i + a4*p4i;

      u0r = p0r + a1*p1r + a5*p2r;
      u0i = p0i + a1*p1i + a5*p2i;

      p4r = p3r; p3r = p2r; p2r = p1r; p1r = p0r;
      p4i = p3i; p3i = p2i; p2i = p1i; p1i = p0i;

     /*==========================================
      * Envelope
      *==========================================
      * env0 = sqrt(u0r*u0r+u0i*u0i) * gain;
      */

     /*==========================================
      * Smoothed envelope by temporal integration 
      *==========================================
      */
      senv1 = senv[i] = sqrt(u0r*u0r+u0i*u0i) * gain + intdecay*senv1;
    }
    
   /*==================================================================================
    * we take the mean of the smoothed envelope as the energy value in each frame
    * rather than simply sampling it.
    * ratemap(c,:) = intgain.*mean(reshape(smoothed_env,frameshift_samples,numframes));
    *==================================================================================
    */
    for (j=0; j<numframes; j++)
    {
      sumEnv = 0.0;
      for (i=j*frameshift_samples; i<(j+1)*frameshift_samples; i++)
      {
        sumEnv += senv[i];
      }
      result[chan+numchans*j] = intgain * sumEnv / frameshift_samples;
    }
  }

  if (strcmp(compression, "cuberoot") == 0)
  {
    for (i=0; i<numchans*numframes; i++)
      result[i] = pow(result[i], 0.3);
  }
  else if (strcmp(compression, "log") == 0)
  {
    for (i=0; i<numchans*numframes; i++)
      result[i] = 20 * log10(result[i]);
  }

  free(senv);
};

