/*
This file is part of INDDGO.

Copyright (C) 2012, Oak Ridge National Laboratory 

This product includes software produced by UT-Battelle, LLC under Contract No. 
DE-AC05-00OR22725 with the Department of Energy. 

This program is free software; you can redistribute it and/or modify
it under the terms of the New BSD 3-clause software license (LICENSE). 

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
LICENSE for more details.

For more information please contact the INDDGO developers at: 
inddgo-info@googlegroups.com

*/

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <ctime>

#include <Log.h>

#if WIN32
// Windows doesn't have timezone and gettimeofday()
// Wrapper taken from
// http://suacommunity.com/dictionary/gettimeofday-entry.php
#include <time.h>
#include <windows.h>

#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64

struct timezone
{
	int  tz_minuteswest; /* minutes W of Greenwich */
	int  tz_dsttime;     /* type of dst correction */
};

int gettimeofday(struct timeval *tv, struct timezone *tz)
{
	// Define a structure to receive the current Windows filetime
	FILETIME ft;

	// Initialize the present time to 0 and the timezone to UTC
	unsigned __int64 tmpres = 0;
	static int tzflag = 0;

	if (NULL != tv)
	{
		GetSystemTimeAsFileTime(&ft);

		// The GetSystemTimeAsFileTime returns the number of 100 nanosecond 
		// intervals since Jan 1, 1601 in a structure. Copy the high bits to 
		// the 64 bit tmpres, shift it left by 32 then or in the low 32 bits.
		tmpres |= ft.dwHighDateTime;
		tmpres <<= 32;
		tmpres |= ft.dwLowDateTime;

		// Convert to microseconds by dividing by 10
		tmpres /= 10;

		// The Unix epoch starts on Jan 1 1970.  Need to subtract the difference 
		// in seconds from Jan 1 1601.
		tmpres -= DELTA_EPOCH_IN_MICROSECS;

		// Finally change microseconds to seconds and place in the seconds value. 
		// The modulus picks up the microseconds.
		tv->tv_sec = (long)(tmpres / 1000000UL);
		tv->tv_usec = (long)(tmpres % 1000000UL);
	}

	if (NULL != tz)
	{
		if (!tzflag)
		{
			_tzset();
			tzflag++;
		}

		// Adjust for the timezone west of Greenwich
		tz->tz_minuteswest = _timezone / 60;
		tz->tz_dsttime = _daylight;
	}

	return 0;
};

#endif


// global variables defined for log
double __log_stime__;
FILE *__log_file__;       
FILE *__error_file__;
int __log_level__;

#ifndef __PARALLEL__
struct timeval tval;
struct timezone tzval;
#endif  // __PARALLEL__

void log_init (const char *filename, const char *errorfile, int loglevel)
{
#ifndef __PARALLEL__
	gettimeofday(&tval, &tzval);
#endif
	__log_stime__ = STIME;
	__log_level__ = 0;
	if (!(loglevel < 0))
		__log_level__ = loglevel;

#ifdef __LOG_ENABLE__

	__log_file__ = fopen (filename, "a+");
	if (!__log_file__)
	{
		fprintf (stderr, "Unable to open __log_file__ for writing\n");
	}

	if (!errorfile)
	{
		errorfile = filename;
	}

	__error_file__ = fopen (errorfile, "a+");
	if (!__error_file__)
	{
		fprintf (stderr, "Unable to open __error_file__ for writing\n");
	}

#endif /* __LOG_ENABLE__ */
}


void log_close ()
{
#ifdef __LOG_ENABLE__
	GEN("############################################################\n");
	fclose (__log_file__);
	fclose (__error_file__);
#endif /* __LOG_ENABLE__ */
}
