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

#ifndef _BINTK_T_
#define _BINTK_T_

//This is necessary for newer versions of gcc where ULLONG_MAX is not defined
//in limits.h 
#include <climits> 
#include <stdio.h>
#include <stdlib.h>
#ifdef __MADNESS__
#include <world/archive.h>
#endif /* __MADNESS__ */

// set to 0 to use uint128_t if you have it
#define NOUINT128 1

#if NOUINT128 || __CYGWIN__ || WIN32 || _WIN32
// Assume we don't have a uint128_t type on windows or cygwin
#define BIGINT_WORD unsigned long long
// Below is cast as int to avoid compiler unsigned vs. signed warning...
#define BIGINT_WORD_SIZE ((int)((sizeof(unsigned long long))<<3))
#define BIGINT_MAX_WORD_VAL ULLONG_MAX
#else
// Assume we have the uint128_t type
#define BIGINT_WORD  __uint128_t
// Below is cast as int to avoid compiler unsigned vs. signed warning...
#define BIGINT_WORD_SIZE ((int)((sizeof(__uint128_t))<<3))
#define BIGINT_MAX_WORD_VAL  (__uint128_t)((__int128_t)(-1))
#endif

#define BIGINT_ONE (BIGINT_WORD)1

// Represents a BIGINT_WORD_SIZE*S-bit integer as S BIGINT_WORD_SIZE-bit words
// word[x-1]word[x-2]...word[3]word[2]word[1]word[0]
// Should be independent of a machine's endian-ness.
class bigint_t
{
public:
	friend bool operator==(bigint_t const& val1, bigint_t const& val2);

	// The # of words
	int S;
	// The data is here
	BIGINT_WORD *words;

	bigint_t()
	{
		this->S = 1;
		this->words = new BIGINT_WORD[this->S];
		for (int i = 0; i < this->S; i++)
			this->words[i] = 0;
	}
	;

	// Constructor for a bigint_t
	inline bigint_t(int num_words)
	{
		this->S = num_words;
		this->words = new BIGINT_WORD[this->S];
		// Set to zeros
		for (int i = 0; i < this->S; i++)
			this->words[i] = 0;
	}
	;

	// Copy constructor - allocates memory and fills with values from source.
	inline bigint_t(const bigint_t &source)
	{
		this->S = source.S;
		this->words = new BIGINT_WORD[this->S];
		for (int i = 0; i < this->S; i++)
			if (source.words)
				this->words[i] = source.words[i];
	}
	;

	// Destructor
	inline ~bigint_t()
	{
		if (this->words != NULL)
			delete[] this->words;
		this->words = NULL;
	}
	;

	inline void sprint(char *buff)
	{
		// Assumes buff is big enough!!
		int nchars = 0;
#if NOUINT128 || __CYGWIN__ || WIN32 || _WIN32
		for (int i = this->S - 1; i >= 1; i--)
			nchars += sprintf(buff + nchars, "%016llx", words[i]);
		sprintf(buff + nchars, "%016llx", words[0]);
#else //for use with uint_128_t
		for(int i=this->S-1;i>=1;i--)
		{
			nchars+=sprintf(buff+nchars,"%16llx",(unsigned long long)(words[i]>>64));
			nchars+=sprintf(buff+nchars,"%16llx",(unsigned long long)(words[i]&ULLONG_MAX));
		}
		nchars+=sprintf(buff+nchars,"%16llx|",(unsigned long long)(words[0]>>64));
		sprintf(buff+nchars,"%16llx\n",(unsigned long long)(words[0]&ULLONG_MAX));
#endif

	}
	;

	// Prints out the words in hex
	void print(FILE *stream)
	{
#if NOUINT128 || __CYGWIN__ || WIN32 || _WIN32
		for (int i = S - 1; i >= 1; i--)
			fprintf(stream, "%016llx", words[i]);
		fprintf(stream, "%016llx\n", words[0]);
		//printing if we use uint_128_t
#else
		for(int i=S-1;i>=1;i--)
		{
			fprintf(stream,"%16llx|",(unsigned long long)(words[i]>>64));
			fprintf(stream,"%16llx|",(unsigned long long)(words[i]&ULLONG_MAX));
		}
		fprintf(stream,"%16llx|",(unsigned long long)(words[0]>>64));
		fprintf(stream,"%16llx\n",(unsigned long long)(words[0]&ULLONG_MAX));
#endif

		return;
	}
	;

	// Prints out the words in hex with a newline b/w words
	void dprint(FILE *stream)
	{
#if NOUINT128 || __CYGWIN__ || WIN32 || _WIN32
		for (int i = S - 1; i >= 1; i--)
			fprintf(stream, "%016llx\n", words[i]);
		fprintf(stream, "%016llx\n", words[0]);
		//printing if we use uint_128_t
#else
		for(int i=S-1;i>=1;i--)
		{
			fprintf(stream,"%16llx",(unsigned long long)(words[i]>>64));
			fprintf(stream,"%16llx\n",(unsigned long long)(words[i]&ULLONG_MAX));
		}
		fprintf(stream,"%16llx",(unsigned long long)(words[0]>>64));
		fprintf(stream,"%16llx\n",(unsigned long long)(words[0]&ULLONG_MAX));
#endif

		return;
	}
	;

	// The &= operator (AND)
	inline bigint_t& operator&=(const bigint_t& val)
	{
		for (int i = 0; i < this->S; i++)
			this->words[i] &= val.words[i];
		return *this;
	}
	;

	// The |= operator (OR)
	inline bigint_t& operator|=(const bigint_t& val)
	{
		for (int i = 0; i < this->S; i++)
			this->words[i] |= val.words[i];
		return *this;
	}
	;

	// The ^= operator (XOR)
	inline bigint_t& operator^=(const bigint_t& val)
	{
		for (int i = 0; i < this->S; i++)
			this->words[i] ^= val.words[i];
		return *this;
	}
	;

	// The ++ operator
	inline bigint_t operator++()
	{
		// Special case with 1 word - assume we don't have BIGINT_MAX_WORD_VAL w/
		// size 1!!
		if (this->S == 1)
		{
			this->words[0]++;
			return *this;
		}

		int k = 0;
#if NOUINT128 || __CYGWIN__ || WIN32 || _WIN32
		while (k <= this->S && this->words[k] == BIGINT_MAX_WORD_VAL)
			k++;
		//using uint_128_t
#else
		while(k<=this->S &&
			((this->words[k]>>64) == ULLONG_MAX) &&
			((this->words[k]&ULLONG_MAX) == ULLONG_MAX)
			)
			k++;
#endif
		// Words 0,1,...,k-1 are 0xff...ff - set them to zero and increment word k by 1
		for (int i = 0; i <= k - 1; i++)
			this->words[i] = (BIGINT_WORD) 0;
		this->words[k]++;

		return *this;
	}
	;

	//The << operator - note that we cross word boundaries
	inline bigint_t operator<<(int shift)
	{
		BIGINT_WORD temp, carry = 0, shift_mask;
		// Is the copy really necessary?
		bigint_t ans(this->S);//(this->num_words);
		ans = *this;
		int i, j = 0, r;

		// If shift>=BIGINT_WORD_SIZE, write shift=BIGINT_WORD_SIZE*j+r where
		// 0<=r<=BIGINT_WORD_SIZE-1
		if (shift >= BIGINT_WORD_SIZE)
		{
			j = 1;
			while (shift >= BIGINT_WORD_SIZE * j)
				j++;
			j--;
			// write shift=BIGINT_WORD_SIZE*j+r
			r = shift - BIGINT_WORD_SIZE * j;

			// Shift the entire words
			for (i = this->S - 1; i >= j; i--)
				ans.words[i] = ans.words[i - j];
			// The bottom j words must be 0
			for (i = 0; i < j; i++)
				ans.words[i] = (BIGINT_WORD) 0;
		}
		else
			r = shift;

		shift_mask = (((BIGINT_ONE) << r) - 1) << (BIGINT_WORD_SIZE - r);
		for (i = j; i < this->S; i++)
		{
			// Remember its orig value
			temp = ans.words[i];
			// Shift the word
			ans.words[i] <<= r;
			// Add in whatever necessary from the lower word
			if (carry)
				ans.words[i] |= carry;
			// Now look at temp to see if we need to modify words[i+1]
			carry = (temp & shift_mask) >> (BIGINT_WORD_SIZE - r);
		}

		//ans.print(stderr);
		return ans;
	}
	;

	// The = assignment operator - assumes already allocated and of right size!
	inline bigint_t& operator=(const bigint_t& val)
	{

		for (int i = 0; i < this->S; i++)
			this->words[i] = val.words[i];

		return *this;
	}
	;

	// The = operator for a ULL - also assumes allocated
	inline bigint_t& operator=(BIGINT_WORD &val)
	{
		this->words[0] = val;
		for (int i = 1; i < this->S; i++)
			this->words[i] = 0;
		return *this;
	}
	;

	// The = operator for a const ULL - also assumes allocated
	inline bigint_t& operator=(const BIGINT_WORD &val)
	{
		this->words[0] = val;
		for (int i = 1; i < this->S; i++)
			this->words[i] = 0;
		return *this;
	}
	;

	// Function to set all words to 0
	inline void zeroize()
	{
		// Set to zeros
		for (int i = 0; i < this->S; i++)
			this->words[i] = (BIGINT_WORD) 0;
	}

	inline void all_ones(int k)
	{
		if (this->S == 1)
		{
			this->words[0] = ((BIGINT_ONE) << k);
			this->words[0]--;
			return;
		}
		if (k > BIGINT_WORD_SIZE * this->S)
		{
			fprintf(stderr, "Asking for too big of an all_ones mask! k=%d\n", k);
			exit(-1);
		}
		int i = 1, j;
		while (k >= BIGINT_WORD_SIZE * i)
			i++;
		// write k=BIGINT_WORD_SIZE*a+r
		int r = k - BIGINT_WORD_SIZE * (i - 1);

		this->zeroize();// ans is all zeros
		//word_0, word_1, word_{i-2} are all ones
		for (j = 0; j <= i - 2; j++)
			this->words[j] = BIGINT_MAX_WORD_VAL;
		// In words[i-1], set the first r bits to 1
		for (j = 0; j < r; j++)
			this->words[i - 1] |= (BIGINT_ONE << j);

		return;
	}

	// Set this to 2^x
	inline void two_x(int x)
	{
		if (this->S == 1)
		{
			this->words[0] = ((BIGINT_ONE) << x);
			return;
		}

		if (x > BIGINT_WORD_SIZE * this->S)
		{
			fprintf(stderr, "%s:  too big a number!\n", __FUNCTION__);
			exit(-1);
		}

		this->zeroize();
		int i, r;

		// Find which word x is in
		i = 1;
		while (x >= BIGINT_WORD_SIZE * i)
			i++;
		r = x - BIGINT_WORD_SIZE * (i - 1);
		this->words[i - 1] |= (BIGINT_ONE << r);

		return;
	}

	// Set this=(1<<p1)|(1<<p2)
	inline void two_ones(int p1, int p2)
	{
		if (this->S == 1)
		{
			this->words[0] = ((BIGINT_ONE) << p1 | (BIGINT_ONE) << p2);
			return;
		}

		if (p1 > BIGINT_WORD_SIZE * this->S || p2 > BIGINT_WORD_SIZE * this->S)
		{
			fprintf(stderr, "%s:  too big a number! p1=%d, p2=%d; S=%d\n",
				__FUNCTION__, p1, p2, this->S);
			exit(-1);
		}
		this->zeroize();
		int i, r;

		// Find which word p1 is in
		i = 1;
		while (p1 >= BIGINT_WORD_SIZE * i)
			i++;
		r = p1 - BIGINT_WORD_SIZE * (i - 1);
		this->words[i - 1] |= (BIGINT_ONE << r);

		// Find which word p2 is in
		i = 1;
		while (p2 >= BIGINT_WORD_SIZE * i)
			i++;
		r = p2 - BIGINT_WORD_SIZE * (i - 1);
		this->words[i - 1] |= (BIGINT_ONE << r);

		return;
	}

	// Sets bit p in this
	inline void set_bit(int p)
	{
		if (this->S == 1)
		{
			this->words[0] |= ((BIGINT_ONE) << p);
			return;
		}
		// Find which word p is in
		int i = 1;
		while (p >= BIGINT_WORD_SIZE * i)
			i++;
		int r = p - BIGINT_WORD_SIZE * (i - 1);

		// Set position r in the correct word and return
		this->words[i - 1] |= (BIGINT_ONE << r);
		return;
	}

	// XORs bit p in this
	inline void xor_bit(int p)
	{
		if (this->S == 1)
		{
			this->words[0] ^= ((BIGINT_ONE) << p);
			return;
		}
		// Find which word p is in
		int i = 1;
		while (p >= BIGINT_WORD_SIZE * i)
			i++;
		int r = p - BIGINT_WORD_SIZE * (i - 1);

		// Set position r in the correct word and return
		this->words[i - 1] ^= (BIGINT_ONE << r);
		return;
	}

	// ORs bit p in this
	inline void or_bit(int p)
	{
		if (this->S == 1)
		{
			this->words[0] |= ((BIGINT_ONE) << p);
			return;
		}
		// Find which word p is in
		int i = 1;
		while (p >= BIGINT_WORD_SIZE * i)
			i++;
		int r = p - BIGINT_WORD_SIZE * (i - 1);

		// Set position r in the correct word and return
		this->words[i - 1] |= (BIGINT_ONE << r);
		return;
	}

	// Return true if bit p is set, false otherwise.
	inline bool test_bit(int p)
	{
		if (this->S == 1)
		{
			if (((this->words[0]) & ((BIGINT_ONE) << p)))
				return true;
			else
				return false;
		}

		// Find which word p is in
		int i = 1;
		while (p >= BIGINT_WORD_SIZE * i)
			i++;
		int r = p - BIGINT_WORD_SIZE * (i - 1);

		// Check position r in the correct word and return true/false
		if ((this->words[i - 1]) & (BIGINT_ONE << r))
			return true;
		return false;

	}

	inline bool is_zero()
	{
		for (int i = 0; i < this->S; i++)
			if (this->words[i] != 0)
				return false;
		return true;
	}

	// Inserts a zero in bit position p (counting from the right, starting
	// at 0).
	inline void insert_zero_bit(int p)
	{
		if (this->is_zero())
		{
			//print_message(0,"%s returning 0\n",__FUNCTION__);
			return;
		}

		if (this->S == 1)
		{
			BIGINT_WORD mask = (BIGINT_WORD) 0;
			for (int i = BIGINT_WORD_SIZE - 1; i >= p; i--)
				mask |= ((BIGINT_ONE) << i);
			this->words[0] = ((this->words[0] & mask) << 1) | (this->words[0]
			& (((BIGINT_ONE) << p) - 1));
			return;
		}

		bigint_t temp(this->S);
		bigint_t mask(this->S);
		int i = 1;
		while (p >= BIGINT_WORD_SIZE * i)
			i++;

		// Create an all ones mask and then just shift over p spots
		// to get the high (-p) bits as we can just shift them into
		// oblivion at the left end
		for (i = 0; i < this->S; i++)
			mask.words[i] = BIGINT_MAX_WORD_VAL;
		mask = (mask << p);

		// Seems like we need a copy of *this
		temp = (*this);
		//this->and_with(&mask);
		(*this) &= mask;
		(*this) = (*this) << 1;
		mask.all_ones(p);
		temp &= mask;
		(*this) |= temp;
		return;
	}
	;

	// Removes the bit in position p(counting from the right, starting at 0).
	// Then squeezes the two halves together:
	// remove_bit(0110,2)=010
	// This is currently slow...
	inline void remove_bit(int p)
	{
		if (this->is_zero())
		{
			//print_message(0,"%s returning 0\n",__FUNCTION__);
			return;
		}

		int i;
		bigint_t temp(this->S);
		// bits 0,1,...,p-1 are unchanged
		for (i = 0; i < p; i++)
		{
			if (this->test_bit(i))
				temp.set_bit(i);
		}
		// bits p+1 are slid over to the right one spot
		for (i = p + 1; i < this->S * BIGINT_WORD_SIZE; i++)
		{
			if (this->test_bit(i))
				temp.set_bit(i - 1);
		}
		// Now copy temp over to this
		for (i = 0; i < this->S; i++)
			this->words[i] = temp.words[i];

		return;
	}
	;

	inline bool operator<(const bigint_t& val)
	{
		for (int i = this->S - 1; i >= 0; i--)
		{
			if (this->words[i] < val.words[i])
				return true;
			if (this->words[i] > val.words[i])
				return false;
		}

		// They must be equal if we get here, so return false
		return false;
	}
	;

	// The < operator
	inline bool operator<(const bigint_t val) const
	{
		for (int i = this->S - 1; i >= 0; i--)
		{
			if (this->words[i] < val.words[i])
				return true;
			if (this->words[i] > val.words[i])
				return false;
		}

		// They must be equal if we get here, so return false
		return false;
	}
	;

#ifdef __MADNESS__
	template<class Archive>
	void serialize(Archive& ar)
	{
		int i = 0;
		ar & this->S;
		for (i = 0; i < this->S; i++)
			ar & words[i];
	}
#endif /* __MADNESS__ */
};

class int_bigint
{
public:
	bigint_t *w;
	int k;

	int_bigint(int S)
	{
		this->w = new bigint_t(S);
	}
	;
	~int_bigint()
	{
		if (this->w)
			delete this->w;
	}
	;

#ifdef __MADNESS__
	template<class Archive>
	void serialize(Archive& ar)
	{
		ar & this->k;
		ar & *this->w;
	}
#endif /* __MADNESS__ */
};

#endif

