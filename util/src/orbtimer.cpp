/*
   This file is part of SystemConfidence.

   Copyright (C) 2012, UT-Battelle, LLC.

   This product includes software produced by UT-Battelle, LLC under Contract No.
   DE-AC05-00OR22725 with the Department of Energy.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the New BSD 3-clause software license (LICENSE).

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   LICENSE for more details.

   For more information please contact the SystemConfidence developers at:
   systemconfidence-info@googlegroups.com

 */

/**************************************************************************
* ORBtimer -- The Oak Ridge Benchmarking Timer Library
*
* Author: Jeffery A. Kuehn
* Copyright (C) 2009 UT-Battelle, LLC.
* All Rights Reserved.
*
**************************************************************************/
#if !WIN32
  #if !CYGWIN

    #include <unistd.h>
    #include <sys/time.h>
    #include "orbconfig.h"

    #define  ORBTIMER_LIBRARY
    #include "orbtimer.h"

static ORB_tick_t Csum = 0;
static ORB_tick_t Gsum = 0;
static ORB_tick_t nsamples = 0;
static ORB_tick_t ndummy = 0;
  #endif
#endif

void ORB_calibrate(){
    #if WIN32 || CYGWIN || _CYGWIN
    return;
    #else
    int i, j;
    double seconds;
    struct timeval tv1, tv2;
    ORB_t t1, t2, t3, cycles;
    ORB_tick_t nsam;
    ORB_tick_t cmin, gmin, csum, gsum, c21, c32;

    #if defined(ORB_IS_FIXEDFREQUENCY)
    ORB_ref_freq = ORB_IS_FIXEDFREQUENCY;
    #else
    #if defined(ORB_IS_FLOATINGPOINT)
    #   error .........................................................
    #   error .... Calibration algorithm does not support floating ....
    #   error .... point timers without a known fixed frequency    ....
    #   error .........................................................
    #endif /* ORB_IS_FLOATINGPOINT  */
    #endif /* ORB_IS_FIXEDFREQUENCY */
    cmin = ORB_min_lat_cyc;
    gmin = GTD_min_lat_cyc + ORB_avg_lat_cyc;
    for(j = 0; j < 4; j++){
        nsam = csum = gsum = 0;         /* keep only last */
        for(i = 0; i < 1000000; i++){           /* sample */
            ORB_read(t1);
            ORB_read(t2);
            gettimeofday(&tv1, 0);
            ORB_read(t3);
            c21 = ORB_cycles_u(t2, t1);
            c32 = ORB_cycles_u(t3, t2);
            if((c21 >= 0) && (c32 >= 0)){               /* GTD timers aren't monotonic */
                if(c21 < cmin){
                    cmin = c21;
                }
                if(c32 < gmin){
                    gmin = c32;
                }
                csum += c21;
                gsum += c32;
                nsam++;
            }
        }
        ndummy += nsam + csum + gsum;
    }
    Csum += csum;
    Gsum += gsum;
    nsamples += nsam;
    ORB_avg_lat_cyc = (Csum + (nsamples >> 1)) / nsamples;
    ORB_min_lat_cyc = cmin;
    GTD_avg_lat_cyc = (Gsum - Csum + (nsamples >> 1)) / nsamples;
    GTD_min_lat_cyc = gmin - ORB_avg_lat_cyc;
    #if !defined(ORB_IS_FIXEDFREQUENCY) /* discover frequency */
    gettimeofday(&tv1, 0);
    ORB_read(t1);
    sleep(5);                   /* region = sleep()+2*ORB()+(0.5+0.5)*gtd() */
    ORB_read(t2);
    gettimeofday(&tv2, 0);
    seconds = ((double)(tv2.tv_sec - tv1.tv_sec)) + ((double)(tv2.tv_usec - tv1.tv_usec)) / 1.0e+6;
    cycles = ORB_cycles_u(t2, t1) + ORB_avg_lat_cyc + GTD_avg_lat_cyc;
    ORB_ref_freq = ((double)(cycles)) / seconds;
    #endif                      /* ORB_IS_FIXEDFREQUENCY */
    ORB_avg_lat_sec = ORB_avg_lat_cyc / ORB_ref_freq;
    ORB_min_lat_sec = ORB_min_lat_cyc / ORB_ref_freq;
    GTD_avg_lat_sec = GTD_avg_lat_cyc / ORB_ref_freq;
    GTD_min_lat_sec = GTD_min_lat_cyc / ORB_ref_freq;
    #endif // if WIN32 || CYGWIN || _CYGWIN
} // ORB_calibrate

