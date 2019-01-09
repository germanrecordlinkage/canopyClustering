// QueryResult.h
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

#ifndef QUERYRESULT_H
#define QUERYRESULT_H

// Instances of QueryResultNode store one match and
// the corresponding similarity coefficient within a sorted
// tree. If sorting is disabled, QueryResultNodes are simply
// used as nodes of a linked list.
 
typedef struct QueryResultNodeStruct {
	char *mQueryId;				// ID of query object
	char *mMatchId;				// ID of match object
	double mSimilarity;			// corresponding similarity coefficient
	struct QueryResultNodeStruct *mLeft;	// left subtree (higher similarity coeffs)
	struct QueryResultNodeStruct *mRight;	// right subtree (lower similarity coeffs)
} QueryResultNode;

// Objects of class QueryResult store search results,
// which consist of two ids and the computed
// similarity coefficient. The QueryResult is designed as
// a sorted tree. On construction the QueryResult can be
// flagged to sort the results while collecting them by
// the similarity coefficient or not. 

class QueryResult {
	private:
	
	QueryResultNode *mRootNode;		// root node of sorted tree or linked list
	int mSize;				// size of result set
	int mSort;				// flag if the query results have ro be sorted

	void deleteNode(QueryResultNode *node);	// delete a node and all its subnodes

	public:

	long long mStat1;
	long long mStat2;
	long long mStat3;

	// constructor
	QueryResult(int sort);

	// destructor
	~QueryResult();
	
	// add a Fingerprint and the corresponding Tanimoto coefficient to the query result
	void add(char *queryId, char *machtId, double similarity);
	
	// return root node for reading
	inline QueryResultNode *getRootNode() {
		return mRootNode;
	}

	// return size of result set
	inline int getSize() {
		return mSize;
	}
};
#endif
