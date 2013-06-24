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
 **************************************************************************
 * API TYPES:
 * ORB_t
 *         An opaque type compatible with the underlying native timer
 *	   this is the internal/opaque type & may be int, double, or struct
 * ORB_tick_t
 *         Syntactic shorthand for 'unsigned long long', 64bit unsigned int
 *	   this is the exposed 'tick' type
 **************************************************************************
 * API FUNCTIONS/MACROS:
 * void
 * ORB_read(T)
 *         sets (ORB_t)T to the value of the current counter. ORB_t is
 *         treated as an opaque type for portability across systems
 *         with different timer implementations
 * void
 * ORB_calibrate()
 *         Initialize ORB() and estimate overheads and ref freq
 *         subsequent calls to ORB_calibrate() further refine estimates
 * ORB_tick_t
 * ORB_cycles(ORB_t end, ORB_t start)
 * ORB_cycles_m(ORB_t end, ORB_t start)
 *         Return the number of cycles between end and start
 *         adjusting for the minimum overhead of ORB_RDTIMER().
 * ORBtick_t
 * ORB_cycles_a(ORB_t end, ORB_t start)
 *         Return the number of cycles between end and start
 *         adjusting for the average overhead of ORB_RDTIMER().
 * ORBtick_t
 * ORB_cycles_u(ORB_t end, ORB_t start)
 *         Return the number of cycles between end and start
 *         unadjusted for the overhead of ORB_RDTIMER().
 * double
 * ORB_seconds(ORB_t end, ORB_t start)
 * ORB_seconds_m(ORB_t end, ORB_t start)
 *         Return the number of seconds between end and start
 *         adjusting for the minimum overhead of ORB_RDTIMER().
 * double
 * ORB_seconds_a(ORB_t end, ORB_t start)
 *         Return the number of seconds between end and start
 *         adjusting for the average overhead of ORB_RDTIMER().
 * double
 * ORB_seconds_u(ORB_t end, ORB_t start)
 *         Return the number of seconds between end and start
 *         unadjusted for the overhead of ORB_RDTIMER().
 **************************************************************************
 * API-EXPOSED USEFUL VALUES:
 * (ORB_tick_t) GTD_MINLAT        gettimeofday() minimum latency in ticks
 * (double)     GTD_MINLATSEC     gettimeofday() minimum latency in seconds
 * (ORB_tick_t) GTD_AVGLAT        gettimeofday() average latency in ticks
 * (double)     GTD_AVGLATSEC     gettimeofday() average latency in seconds
 * (ORB_tick_t) ORB_MINLAT        ORB_read() minimum latency in ticks
 * (double)     ORB_MINLATSEC     ORB_read() minimum latency in seconds
 * (ORB_tick_t) ORB_AVGLAT        ORB_read() average latency in ticks
 * (double)     ORB_AVGLATSEC     ORB_read() average latency in seconds
 * (ORB_tick_t) ORB_IREFFREQ      timer ticks per second
 * (double)     ORB_REFFREQ       timer ticks per second (floating point)
 * (int)        GTD_REFFREQ       GTD timer ticks per second (integer)
 *************************************************************************/

#ifndef HAVE_ORBTIMER
#    define HAVE_ORBTIMER
#    ifdef ORBTIMER_LIBRARY
#        define ORBEXTERN(x,y) x=y
#    else
#        define ORBEXTERN(x,y) extern x
#    endif

void ORB_calibrate();
typedef unsigned long long ORB_tick_t;
ORBEXTERN(int GTD_ref_freq, 1000000);
ORBEXTERN(double ORB_ref_freq, 0.0);
ORBEXTERN(int ORB_FreqTest, 0);
ORBEXTERN(ORB_tick_t GTD_avg_lat_cyc, 0);
ORBEXTERN(ORB_tick_t GTD_min_lat_cyc, 1000000);
ORBEXTERN(ORB_tick_t ORB_avg_lat_cyc, 0);
ORBEXTERN(ORB_tick_t ORB_min_lat_cyc, 1000000);
ORBEXTERN(double GTD_avg_lat_sec, 0.0);
ORBEXTERN(double GTD_min_lat_sec, 1000000.0);
ORBEXTERN(double ORB_avg_lat_sec, 0.0);
ORBEXTERN(double ORB_min_lat_sec, 1000000.0);
#    define GTD_AVGLAT            ( (ORB_tick_t) (GTD_avg_lat_cyc))
#    define GTD_AVGLATSEC         ( (double)    (GTD_avg_lat_sec) )
#    define GTD_MINLAT            ( (ORB_tick_t) (GTD_min_lat_cyc))
#    define GTD_MINLATSEC         ( (double)    (GTD_min_lat_sec) )
#    define GTD_REFFREQ           ( (int   )    (1000000))
#    define ORB_AVGLAT            ( (ORB_tick_t) (ORB_avg_lat_cyc))
#    define ORB_AVGLATSEC         ( (double)    (ORB_avg_lat_sec) )
#    define ORB_MINLAT            ( (ORB_tick_t) (ORB_min_lat_cyc))
#    define ORB_MINLATSEC         ( (double)    (ORB_min_lat_sec) )
#    define ORB_REFFREQ           ( (double)    (ORB_ref_freq))

#    if defined(TIMER_X86_64)
       /**********************************************************
        * NOTES:
        * A hardware cycle counter is preferred if it is 
        * predictable -- ie. the cost is relatively constant.
	*
	* If multiple hardware counters are available, prefer
	* one located in-core as it will have lower latency and
	* less variability.
        *
        * This module is a compromise for Intel and AMD
        * processors running in 64 bit mode. Hardware Performance
	* Counters can have less latency and variability, but they
	* are available on all systems.
        * 
        * RDTSC alone is faster, but without a serializing
        * instruction, its result can be impacted by
        * out-of-order execution.  While the documentation
        * suggests using CPUID to serialize the instruction
        * stream, the CPUID instruction timing is
        * enormously variable, and thus unreliable. MFENCE
        * on the other hand is predictable. Newer Intel
        * and AMD processors provide a serializing
        * version, RDTSCP, however, its timing is not
        * more predictable in its timing than RDTSC, and
        * since we can estimate and subtract the avg cost
        * of timing, the constant cost difference between
        * RDTSC and RDTSCP will not impact our accuracy.
	*
	* ORB_t		-- CPU Cycles
	* ORB_tick_t	-- CPU Cycles
	* ORB_REFFREQ	-- CPU core Frequency
        *********************************************************/
#       ifndef HAVE_ORBTIMER_NATIVE
#           define HAVE_ORBTIMER_NATIVE "X86_64"
typedef unsigned long long ORB_t;
#            define ORB_cycles(T2,T1)     ( (ORB_tick_t) (T2-T1-ORB_min_lat_cyc))
#            define ORB_cycles_m(T2,T1)   ( (ORB_tick_t) (T2-T1-ORB_min_lat_cyc))
#            define ORB_cycles_a(T2,T1)   ( (ORB_tick_t) (T2-T1-ORB_avg_lat_cyc))
#            define ORB_cycles_u(T2,T1)   ( (ORB_tick_t) (T2-T1)                )
#            define ORB_seconds(T2,T1)    ( ((double)   (T2-T1-ORB_min_lat_cyc)) / ORB_ref_freq)
#            define ORB_seconds_m(T2,T1)  ( ((double)   (T2-T1-ORB_min_lat_cyc)) / ORB_ref_freq)
#            define ORB_seconds_a(T2,T1)  ( ((double)   (T2-T1-ORB_avg_lat_cyc)) / ORB_ref_freq)
#            define ORB_seconds_u(T2,T1)  ( ((double)   (T2-T1)                ) / ORB_ref_freq)
#            define ORB_read(T)        __asm__ __volatile__ ( "  \n\t" \
                                                  "mfence           \n\t" \
                                                  "rdtsc            \n\t" \
                                                  "movl %%eax,%%eax \n\t" \
                                                  "salq $32,%%rdx   \n\t" \
                                                  "orq %%rdx,%%rax  \n\t" : "=a" (T) : : "%rdx")
#        else			/* HAVE_NATIVE_TIMER */
#            error "Multiple native timers. " __FILE__ " detected previous value " HAVE_ORBTIMER_NATIVE
#        endif			/* HAVE_ORBTIMER_NATIVE */
#    elif defined(TIMER_PPC64)
#        ifndef HAVE_ORBTIMER_NATIVE
typedef unsigned long ORB_t;
#            define HAVE_ORBTIMER_NATIVE
#            define ORB_cycles(T2,T1)     ( (ORB_tick_t) (T2-T1-ORB_min_lat_cyc))
#            define ORB_cycles_m(T2,T1)   ( (ORB_tick_t) (T2-T1-ORB_min_lat_cyc))
#            define ORB_cycles_a(T2,T1)   ( (ORB_tick_t) (T2-T1-ORB_avg_lat_cyc))
#            define ORB_cycles_u(T2,T1)   ( (ORB_tick_t) (T2-T1)                )
#            define ORB_seconds(T2,T1)    ( ((double)   (T2-T1-ORB_min_lat_cyc)) / ORB_ref_freq)
#            define ORB_seconds_m(T2,T1)  ( ((double)   (T2-T1-ORB_min_lat_cyc)) / ORB_ref_freq)
#            define ORB_seconds_a(T2,T1)  ( ((double)   (T2-T1-ORB_avg_lat_cyc)) / ORB_ref_freq)
#            define ORB_seconds_u(T2,T1)  ( ((double)   (T2-T1)                ) / ORB_ref_freq)
#            define ORB_read(T)       __asm__ __volatile__ ( "  \n\t" \
                                                 "mftb %0           \n\t" : "=r" (T) )
#        else			/* HAVE_ORBTIMER_NATIVE */
#            error "Multiple native timers. " __FILE__ " detected previous value " HAVE_ORBTIMER_NATIVE
#        endif			/* HAVE_ORBTIMER_NATIVE */
#    elif defined(TIMER_MPI_WTIME)
       /**********************************************************
        * NOTES:
        * On some systems MPI_Wtime() is a very predictable on 
        * others, less so. Compare against gettimeofday() for
        * reliability. Prefer CPU Cycle counter if available.
	*
	* ORB_t		-- floating point seconds
	* ORB_tick_t	-- number of MPI_Wtick() intervals
	* ORB_REFFREQ	-- (double)1.0/MPI_Wtick()
        *********************************************************/
#        ifndef HAVE_ORBTIMER_NATIVE
#            define HAVE_ORBTIMER_NATIVE "MPI_WTIME"
#            define ORB_IS_FLOATINGPOINT
#            define ORB_IS_FIXEDFREQUENCY (1.0/MPI_Wtick())
typedef double ORB_t;
extern double MPI_Wtime();
extern double MPI_Wtick();
#            define ORB_cycles(T2,T1)     ( (ORB_tick_t) ((T2-T1-ORB_min_lat_sec)*ORB_ref_freq) )
#            define ORB_cycles_m(T2,T1)   ( (ORB_tick_t) ((T2-T1-ORB_min_lat_sec)*ORB_ref_freq) )
#            define ORB_cycles_a(T2,T1)   ( (ORB_tick_t) ((T2-T1-ORB_avg_lat_sec)*ORB_ref_freq) )
#            define ORB_cycles_u(T2,T1)   ( (ORB_tick_t) ((T2-T1                )*ORB_ref_freq) )
#            define ORB_seconds(T2,T1)    ( (double)    ((T2-T1-ORB_min_lat_sec)             ) )
#            define ORB_seconds_m(T2,T1)  ( (double)    ((T2-T1-ORB_min_lat_sec)             ) )
#            define ORB_seconds_a(T2,T1)  ( (double)    ((T2-T1-ORB_avg_lat_sec)             ) )
#            define ORB_seconds_u(T2,T1)  ( (double)    ((T2-T1                )             ) )
#            define ORB_read(T)        T = MPI_Wtime()
#        else			/* HAVE_ORBTIMER_NATIVE */
#            error "Multiple native timers. " __FILE__ " detected previous value " HAVE_ORBTIMER_NATIVE
#        endif			/* HAVE_ORBTIMER_NATIVE */
#    elif defined(TIMER_GTD)
       /**********************************************************
        * NOTES:
        * gettimeofday() is a standard fallback. At best it is
        * incremented at 1 MHz, however on some systems it could
        * be incremented at 100 Hz... generally not good enough.
        * ORB_t 	-- struct timeval
        * ORB_tick_t 	-- microseconds
        * ORB_REFFREQ 	-- 1 MHz
        *********************************************************/
#        ifndef HAVE_ORBTIMER_NATIVE
#            define HAVE_ORBTIMER_NATIVE "GETTIMEOFDAY"
#            define ORB_IS_FIXEDFREQUENCY (1000000)
typedef struct timeval ORB_t;
#            define ORB_USEC(T)           ( (ORB_tick_t)((ORB_t)(T)).tv_sec * (GTD_REFFREQ) + (ORB_tick_t)((ORB_t)T).tv_usec )
#            define ORB_cycles(T2,T1)     ( (ORB_tick_t)( ORB_USEC(T2) - ORB_USEC(T1) - ORB_min_lat_cyc) )
#            define ORB_cycles_m(T2,T1)   ( (ORB_tick_t)( ORB_USEC(T2) - ORB_USEC(T1) - ORB_min_lat_cyc) )
#            define ORB_cycles_a(T2,T1)   ( (ORB_tick_t)( ORB_USEC(T2) - ORB_USEC(T1) - ORB_avg_lat_cyc) )
#            define ORB_cycles_u(T2,T1)   ( (ORB_tick_t)( ORB_USEC(T2) - ORB_USEC(T1) ) )
#            define ORB_seconds(T2,T1)    ( ( (double)ORB_cycles(T2,T1)   ) / ORB_REFFREQ )
#            define ORB_seconds_m(T2,T1)  ( ( (double)ORB_cycles_m(T2,T1) ) / ORB_REFFREQ )
#            define ORB_seconds_a(T2,T1)  ( ( (double)ORB_cycles_a(T2,T1) ) / ORB_REFFREQ )
#            define ORB_seconds_u(T2,T1)  ( ( (double)ORB_cycles_u(T2,T1) ) / ORB_REFFREQ )
#            define ORB_read(T)        gettimeofday(&T,NULL)
#        else			/* HAVE_ORBTIMER_NATIVE */
#            error "Multiple native timers. " __FILE__ " detected previous value " HAVE_ORBTIMER_NATIVE
#        endif			/* HAVE_ORBTIMER_NATIVE */
#    else			/* TIMER_<type> ... no timer type declared */
#          error ....................................................................
#          error ....................................................................
#          error .... No timer selected in __FILE__
#          error .... Add one of the following to the compile:
#          error ................ -DTIMER_X86_64 ....... uses X86-64 RDTSC instruction
#          error ................ -DTIMER_MPI_WTIME .... uses MPI_Wtime()
#          error ................ -DTIMER_GTD ......... uses gettimeofday()
#          error .... Undeclared variable complaints below are a side effect of this
#          error ....................................................................
#          error ....................................................................
#    endif			/* TIMER_<type> */
#    undef ORBEXTERN
#endif				/* HAVE_ORBTIMER */
