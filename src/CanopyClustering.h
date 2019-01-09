// Canopy.h
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

#ifndef CANOPY_CLUSTERING_H
#define CANOPY_CLUSTERING_H

#include <math.h>
#include "Fingerprint.h"
#include "QueryResult.h"

// class CanopyClustering reveives an array of instances of the class T.
// The array will be clustered into canopies by a given minimal similarity,
// loose threshold and tight threshold. The objects of the same
// clusters will be matched by the given minimal similarity.

template<class T>
class CanopyClustering {
	private:

	long long mSizeA;		// size of object-array part A
	long long mSizeB;		// size of object-array part B

	T **mObjects;			// array of objects
	
	public:

	// constructor
	inline CanopyClustering(T **objects, long long sizeA, long long sizeB) {
		mSizeA = sizeA;
		mSizeB = sizeB;
		mObjects = objects;
	}

	// destructor
	inline ~CanopyClustering() {
		for (int i = 0; i < (mSizeA + mSizeB); i++) {
			if (mObjects[i]) {
				delete mObjects[i];
			}
		}

		delete[] mObjects;
	}

	// perform a search for <minSimilarity> and <looseThreshold>, <tightThreshold> and add the result to <result>
	inline void search(QueryResult *result, double minSimilarity, double looseThreshold, double tightThreshold) {
		double sim;
		long long clusterColSize = 0;
		long long *clusterCollection = new long long[mSizeA + mSizeB];

		for (int center = 0; center < mSizeA; center++) {
			clusterColSize = 0;
			if (mObjects[center]->isInP()) {
				mObjects[center]->removeFromP();
				clusterCollection[clusterColSize] = center;
				clusterColSize++;
				
				for (int i = (center + 1); i < (mSizeA + mSizeB); i++) {
					if (mObjects[i]->isInP()) {    				
						// check estimated similarity
						sim = mObjects[center]->estimatedSimilarity(result, mObjects[i]); 
						if (sim >= looseThreshold) {
							clusterCollection[clusterColSize] = i;
							clusterColSize++;
							if (sim >= tightThreshold) {
								mObjects[i]->removeFromP();
							}
						}
					}
				}

				for (int i = 0; i < clusterColSize; i++) {
					for (int j = (i + 1); j < clusterColSize; j++) {
						if (mObjects[clusterCollection[i]]->getOrigin() != mObjects[clusterCollection[j]]->getOrigin()) {
							mObjects[clusterCollection[i]]->match(result, mObjects[clusterCollection[j]], minSimilarity);
						}
					}
				}
			}
		}
		delete[] clusterCollection;
	}
};
#endif
