// Fingerprint.h
//
// Copyright (c) 2013
// Universitaet Duisburg-Essen
// Campus Duisburg
// Institut fuer Soziologie
// Prof. Dr. Rainer Schnell
// Lotharstr. 65
// 47057 Duisburg 
//
// This file is part of the R-Package "canopyClustering".
//
// "canopyClustering" is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// "canopyClustering" is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with "canopyClustering". If not, see <http://www.gnu.org/licenses/>.

#ifndef FINGERPRINT_H
#define FINGERPRINT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Misc.h"
#include "QueryResult.h"

typedef unsigned int WORDTYPE;			// 32bit-words
#define WORD_LEN 32				// word-length in bits
#define BIT1 1u					// unsigned int literal 1

#define FOLDED_WORDS (128/WORD_LEN)		// array-length for hash-key

// Objects of class Fingerprint hold a single bit-vector of arbitrary length
// The class also provides functions for modifying and checking single bits of
// the bit-vector and to compute its cardinality.
// Furthermore there are functions to compute the Tanimoto coefficient of
// two fingerprints and the faster 128-bit folding-hash estimation of the Tanimoto
// coefficient.

class Fingerprint {
	protected:

	char *mId;				// pointer to ID string
	WORDTYPE *mArray;			// array for stored bits
	int mLength;				// length of fingerprint in bits
	int mOrigin;				// origin, eg fileA=1, fileB=2
	int mCardinality;			// cardinality buffer of fingerprint
	bool mInP;				// indicator for set P
	static int sCardinalityMap[0x10000];	// static 16-bit cardinality-map
	
	// calculate length of word-array
	
	inline int arrayLength() {
		return (mLength - 1) / WORD_LEN + 1;
	}
	
	// initialize word-array
	
	inline void allocate() {
		mLength = MAX(mLength, 128);
		mArray = new WORDTYPE[arrayLength()];
	}

	// calculate cardinality of 32bit-word
	
	inline int cardWord(WORDTYPE word) {
		return sCardinalityMap[word & 0xFFFF] + sCardinalityMap[(word >> 16) & 0xFFFF];
	}
	
	public:

	// constructor for empty fingerprint of given bit-length
	
	inline Fingerprint(int length, int origin) {
		mId = NULL;
		mLength = length;
		mOrigin = origin;
		mInP = true;
		allocate();
		clear();
	}
	
	// constructor for fingerprint based on given ascii-string

	inline Fingerprint(char * id, const char *str, int origin) {
		int l = 0;
		
		mId = id;
		mOrigin = origin;
		
		while ((str[l] == '0') || (str[l] == '1')) {
			l++;
		}

		mLength = l;
		mInP = true;
		
		allocate();
		clear();
	
		for (int i = 0; i < l; i++) {
			if (str[i] != '0') {
				setBit(i);
				mCardinality++;
			}
		}
	}
	
	// destructor
	
	inline ~Fingerprint() {
		if (mId != NULL) {
			delete[] mId;
		}
		delete[] mArray;
	}

	// get id
	
	inline char *getId() {
		return mId;
	}

	inline int getOrigin() {
		return mOrigin;
	}

	// clear all bits
	
	inline void clear() {
		for (int i = 0; i < arrayLength(); i++) {
			mArray[i] = 0;
		}
		mCardinality = 0;
	}
	
	// check if all bits are cleared
	
	inline int isEmpty() {
		for (int i = 0; i < arrayLength(); i++) {
			if (mArray[i] != 0) {
				return 0;
			}
		}
	
		return 1;
	}
	
	// get length in bits
	
	inline int getLength() {
		return mLength;
	}
	
	inline bool isInP() {
		return mInP;
	}

	inline void removeFromP() {
		mInP = false;
	}

	private:
	// get bit at position n
	
	inline WORDTYPE getBit(int n) {
		return (mArray[n / WORD_LEN] & (BIT1 << (n % WORD_LEN)));
	}
	
	// set bit at position n
	
	inline void setBit(int n) {
		mArray[n / WORD_LEN] |= (BIT1 << (n % WORD_LEN));
	}
	
	// unset bit at position n
	
	inline void unsetBit(int n) {
		mArray[n / WORD_LEN] &= 0xFFFF - (BIT1 << (n % WORD_LEN));
	}

	public:
	// compute tanimoto-index

	inline double tanimoto(Fingerprint *print) {
		int count_and = 0;
		int count_or = 0;
		int min, len;

		len = arrayLength();
		min = print->arrayLength();
		
		// handle different bit-length
		if (min <= len) {
			for (int i = min; i < len; i++) {
				count_or += cardWord(mArray[i]);
			}
		} else {
			for (int i = len; i < min; i++) {
				count_or += cardWord(print->mArray[i]);
			}
			min = len;
		}
		
		for (int i = 0; i < min; i++) {
			int a = mArray[i] & print->mArray[i];
			int o = mArray[i] | print->mArray[i];
			count_and += cardWord(a);
			count_or  += cardWord(o);
		}
		
		return ((double) count_and) / count_or;
	}
	
	// generic similarity estimation

	inline double estimatedSimilarity(QueryResult *result, Fingerprint *print) {
		result->mStat3++;
		return tanimoto(print);
	}

	// generic similarity check

	inline void match(QueryResult *result, Fingerprint *print, double minSimilarity) {
		result->mStat1++;
		double sim = tanimoto(print);
		if (sim >= minSimilarity) {
			if (mOrigin == 1) {
				result->add(mId, print->getId(), sim);
			} else {
				result->add(print->getId(), mId, sim);
			}
		}
	}

	// static class function to initialise cardinality-map
	
	static void init();
	
	// convert fingerprint to ascii-string of size n
	
	void copyToString(char *str, int n);
};

class FingerprintXOR : public Fingerprint {
	protected:

	WORDTYPE mHashArray[FOLDED_WORDS];	// 128 Bit folded Hash-Key

	public:

	inline FingerprintXOR(int length, int origin) : Fingerprint(length, origin) {
		fold();
	}

	inline FingerprintXOR(char * id, const char *str, int origin) : Fingerprint(id, str, origin) {
		fold();
	}
	
	// compute hash-key
	
	inline void fold() {
		int len;
		
		len = arrayLength();
		
		for (int i = 0; i < FOLDED_WORDS; i++) {
			mHashArray[i] = mArray[i];
		}
		
		for (int i = FOLDED_WORDS; i < len; i++) {
			mHashArray[i % FOLDED_WORDS] ^= mArray[i];
		}
	}
	
	// compute tanimoto estimation on hash-keys
	
	inline double tanimotoXOR(FingerprintXOR *print, int AB) {

		int xorCount  = cardWord(mHashArray[0] ^ print->mHashArray[0])
			      + cardWord(mHashArray[1] ^ print->mHashArray[1])
			      + cardWord(mHashArray[2] ^ print->mHashArray[2])
			      + cardWord(mHashArray[3] ^ print->mHashArray[3]);
		
		return ((double) (AB - xorCount)) / (AB + xorCount);
	}

	// generic similarity estimation

	inline double estimatedSimilarity(QueryResult *result, FingerprintXOR *print) {
		result->mStat3++;
		return tanimotoXOR(print, mCardinality + print->mCardinality);
	}

	// generic similarity check

	inline void match(QueryResult *result, FingerprintXOR *print, double minSimilarity) {
		double sim;

		result->mStat2++;
		sim = tanimotoXOR(print, mCardinality + print->mCardinality);

		if (sim >= minSimilarity) {
			Fingerprint::match(result, print, minSimilarity);
		}
	}
};
#endif
