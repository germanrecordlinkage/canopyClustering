// Fingerprint.cpp
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <limits>

#include "Fingerprint.h"

// static class member for fast cardinality calculation
// array contains pre-computed cardinalities for all possible 16bit-words

int Fingerprint::sCardinalityMap[0x10000];

// static class function to initialise cardinality-map

void Fingerprint::init() {
	for (int i = 0; i < 0x10000; i++) {
		sCardinalityMap[i] = 0;
		for (int j = 0; j < 16; j++) {
			if ((i & (1 << j)) != 0) {
				sCardinalityMap[i]++;
			}
		}
	}
}

// convert fingerprint to ascii-string of size n
	
void Fingerprint::copyToString(char *str, int n) {
	int l;
	
	l = MIN(mLength, n - 1);
	
	for (int i = 0; i < l; i++) {
		if (getBit(i) == 0) {
			str[i] = '0';
		} else {
			str[i] = '1';
		}
	}
	
	str[l] = 0;
}
