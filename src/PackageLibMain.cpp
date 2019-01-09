// PackageLibMain.cpp
//
// Copyright (c) 2014
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

#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <string.h>

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

#include "Misc.h"
#include "CanopyClustering.h"

// This file contains the pure c-functions for the R-library-interface.
// R_init_useCall		register .Call-Methods
// ccMatchFilesCall		match files
// ccStatistics			print statistics 

// maximal line size = length of ascii representation of fingerprint
#define STRSIZE 4000

static int gMethod = 0;                         // static method indicator

static long long gStat1 = 0;
static long long gStat2 = 0;
static long long gStat3 = 0;
static long long gStatSizeLastSearch = 0;

// check, if a character is considered to be a white-space or seperator
int isWS(char c) {
	return (
		(c == '"') ||
		(c == '\'') ||
		(c == ',') ||
		(c == ';') ||
		(c == ' ') ||
		(c == '\t')
	);
}

// check, if a character is considered to be end of line
int isEOL(char c) {
	return (
		(c == 10) ||
		(c == 13) ||
		(c == 0)
	);
}

// parse prints from file, return number of parsed fields
int parseLine(FILE *in, char *str, long long *idx1, long long *end1, long long *idx2, long long *end2) {
	// read line from file
	if (fgets(str, STRSIZE, in) == NULL) {
		return 0;
	}
	
	// parse line
	
	// find start of first field
	*idx1 = 0;
	while (isWS(str[*idx1]) && (*idx1 < STRSIZE-1)) {
		(*idx1)++;
	}

	// find end of first field
	*end1 = *idx1;
	while (!isWS(str[*end1]) && (*end1 < STRSIZE-1) && !isEOL(str[*end1])) {
		(*end1)++;
	}

	if (!isEOL(str[*end1])) {
		// find start of second field
		*idx2 = *end1;
		while (isWS(str[*idx2]) && (*idx2 < STRSIZE-1)) {
			(*idx2)++;
		}
		
		// find end of second field
		*end2 = *idx2;
		while (!isWS(str[*end2]) && (*end2 < STRSIZE-1) && !isEOL(str[*end2])) {
			(*end2)++;
		}
	} else {
		*idx2 = 0;
		*end2 = 0;
	}

	// cut field 1
	str[*end1] = 0;

	if (*idx2 != *end2) {
		// if there are two fields, cut field 2
		str[*end2] = 0;
		return 2;
	}

	return 1;
}

// copy QueryResults into R vectors for queries, matches and similarity coefficients
void insertQueryResultNodesWithId(long long *idx, SEXP queries, SEXP matches, double *similaritiesPtr, QueryResultNode *node, long long sizeResult) {
	while ((node != NULL) && ((*idx) < sizeResult)) {
		SET_STRING_ELT(queries, *idx, mkChar(node->mQueryId));
		SET_STRING_ELT(matches, *idx, mkChar(node->mMatchId));
		similaritiesPtr[*idx] = node->mSimilarity;		// copy similarity coefficent

		(*idx)++;
		node = node->mRight;
	}
}

// count number of lines in file
long long countFile(const char * filename) {
        FILE *in;
        char str[STRSIZE];
        long long size = 0;

        in = fopen(filename, "r");

        if (in != NULL) {
                while (fgets(str, STRSIZE, in)) {
                        size++;
                }
                fclose(in);
        }

        return(size);
}

// read prints from file
template <class T>
void readFile(const char *filename, T **prints, long long start, long long& size, int& nBits, int origin) {
        FILE *in;
        char str[STRSIZE];
        long long fields;
        long long idx1, end1, idx2, end2;
        char *idStr;
        int len;

        in = fopen(filename, "r");

        if (in != NULL) {
                for (long long i = start; i < (start + size); i++) {
                        // parse line from file
                        fields = parseLine(in, str, &idx1, &end1, &idx2, &end2);

                        if (fields == 0) {
                                size = i - start;
                                break;
                        }

                        if (fields == 1) {
                                // if there is only one field, use line as id
                                idStr = new char[13];
                                sprintf(idStr, "%012" PRId64, i+1);
                                // create fingerprint
                                prints[i] = new T(idStr, str+idx1, origin);
                        } else {
                                // if there are two fields, use first string as id
                                idStr = new char[end1-idx1+1];
                                strcpy(idStr, str+idx1);
                                // create fingerprint
                                prints[i] = new T(idStr, str+idx2, origin);
                        }

                        // compute maximal length of Fingerprints
                        len = prints[i]->getLength();
                        if (len > nBits) {
                                nBits = len;
                        }
                }

                fclose(in);
        }
}

// load files into internal tanimoto fingerprint data structure
template <class T>
void matchTanimoto(QueryResult *queryResult, const char *filenameA, const char *filenameB, double minSimilarity, double looseThreshold, double tightThreshold) {
        T **prints;
        CanopyClustering<T> *cc;
        long long sizePrintsA;
        long long sizePrintsB;
        int nBits;

        // initialize Fingerprint data structure (cardinality-map)
        T::init();

        sizePrintsA = countFile(filenameA);

        if (sizePrintsA == 0) {
                return;
        }

        sizePrintsB = countFile(filenameB);

        if (sizePrintsB == 0) {
                return;
        }

        // create array of Fingerprints
        prints = new T*[sizePrintsA + sizePrintsB];

        nBits = 0;

        readFile(filenameA, prints, 0, sizePrintsA, nBits, 1);
        readFile(filenameB, prints, sizePrintsA, sizePrintsB, nBits, 2);

        // build sn data structure
        cc = new CanopyClustering<T>(prints, sizePrintsA, sizePrintsB);
        cc->search(queryResult, minSimilarity, looseThreshold, tightThreshold);
        delete cc;

        // total count of comparisons for statistics
        gStatSizeLastSearch = sizePrintsA * sizePrintsB;
}

// searching for similar fingerprints
SEXP ccMatchFiles(const char *filenameA, const char *filenameB, double minSimilarity, double looseThreshold, double tightThreshold, int method) {
	SEXP result;
	SEXP names;
	SEXP queries;
	SEXP matches;
	SEXP similarities;
	double *similaritiesPtr;
	long long idx;
	long long sizeResult;
	FILE *in;
	char str[STRSIZE];
	long long fields;
	long long idx1, end1, idx2, end2;
	long long i;
	char *idStr;
	char *strfinger;
	QueryResult queryResult(0);

        gMethod = method;
        gStat1 = 0;
        gStat2 = 0;
        gStat3 = 0;
        gStatSizeLastSearch = 0;

        if (gMethod == 1) {
                // method "tanimoto"

                matchTanimoto<Fingerprint>(&queryResult, filenameA, filenameB, minSimilarity, looseThreshold, tightThreshold);

                gStat1 = queryResult.mStat1;
                gStat3 = queryResult.mStat3;

        } else if (gMethod == 2) {
                // method "tanimotoXOR"

                matchTanimoto<FingerprintXOR>(&queryResult, filenameA, filenameB, minSimilarity, looseThreshold, tightThreshold);

                gStat1 = queryResult.mStat1;
                gStat2 = queryResult.mStat2;
                gStat3 = queryResult.mStat3;
	}

	sizeResult = queryResult.getSize();

	// allocate R data structures for result
	PROTECT(queries = allocVector(STRSXP, sizeResult));
	PROTECT(matches = allocVector(STRSXP, sizeResult));
	PROTECT(similarities = allocVector(REALSXP, sizeResult));
	similaritiesPtr = REAL(similarities);

	// copy result into R data structures
	idx = 0;
	insertQueryResultNodesWithId(&idx, queries, matches, similaritiesPtr, queryResult.getRootNode(), sizeResult);

	// allocate vector for the three result vectors
	PROTECT(result = allocVector(VECSXP, 3));

	SET_VECTOR_ELT(result, 0, queries);
	SET_VECTOR_ELT(result, 1, matches);
	SET_VECTOR_ELT(result, 2, similarities);

	// set name attributes for the two result vectors
	PROTECT(names = allocVector(STRSXP, 3));

	SET_STRING_ELT(names, 0, mkChar("query"));
	SET_STRING_ELT(names, 1, mkChar("match"));
	SET_STRING_ELT(names, 2, mkChar("similarity"));
	setAttrib(result, R_NamesSymbol, names);

	UNPROTECT(5);

	return(result);
}

// Output of statistic values
SEXP ccStatistics() {
	SEXP result;
	SEXP names;
	SEXP values;
	SEXP params;
	SEXP percents;

        double *valuesPtr;
        double *percentsPtr;

        if (gMethod == 1) {
                PROTECT(params = allocVector(STRSXP, 3));
                PROTECT(values = allocVector(REALSXP, 3));
                PROTECT(percents = allocVector(REALSXP, 3));

                valuesPtr = REAL(values);
                percentsPtr = REAL(percents);

                valuesPtr[0] = (double)gStat3;
                valuesPtr[1] = (double)gStat1;
                valuesPtr[2] = (double)gStatSizeLastSearch;

                percentsPtr[0] = (double)gStat3 / gStatSizeLastSearch * 100;
                percentsPtr[1] = (double)gStat1 / gStatSizeLastSearch * 100;
                percentsPtr[2] = 100.0;

                SET_STRING_ELT(params, 0, mkChar("Tanimoto while Clustering"));
                SET_STRING_ELT(params, 1, mkChar("Tanimoto while Matching"));
                SET_STRING_ELT(params, 2, mkChar("Total while Matching"));
        } else if (gMethod == 2) {
                PROTECT(params = allocVector(STRSXP, 4));
                PROTECT(values = allocVector(REALSXP, 4));
                PROTECT(percents = allocVector(REALSXP, 4));

                valuesPtr = REAL(values);
                percentsPtr = REAL(percents);

                valuesPtr[0] = (double)gStat3;
                valuesPtr[1] = (double)gStat2;
                valuesPtr[2] = (double)gStat1;
                valuesPtr[3] = (double)gStatSizeLastSearch;

                percentsPtr[0] = (double)gStat3 / gStatSizeLastSearch * 100;
                percentsPtr[1] = (double)gStat2 / gStatSizeLastSearch * 100;
                percentsPtr[2] = (double)gStat1 / gStatSizeLastSearch * 100;
                percentsPtr[3] = 100.0;

                SET_STRING_ELT(params, 0, mkChar("XOR while Clustering"));
                SET_STRING_ELT(params, 1, mkChar("XOR while Matching"));
                SET_STRING_ELT(params, 2, mkChar("Tanimoto while Matching"));
                SET_STRING_ELT(params, 3, mkChar("Total while Matching"));
        } else {
                return(NULL);
        }

        PROTECT(result = allocVector(VECSXP, 3));

        SET_VECTOR_ELT(result, 0, params);
        SET_VECTOR_ELT(result, 1, values);
        SET_VECTOR_ELT(result, 2, percents);

        // set name attributes for the two result vectors
        PROTECT(names = allocVector(STRSXP, 3));

        SET_STRING_ELT(names, 0, mkChar("Checkpoint"));
        SET_STRING_ELT(names, 1, mkChar("Count"));
        SET_STRING_ELT(names, 2, mkChar("Percentage"));
        setAttrib(result, R_NamesSymbol, names);

        UNPROTECT(5);

	return(result);

}

#ifdef __cplusplus
extern "C" {
#endif

// wrapper for R-function ccMatchFilesCall
SEXP ccMatchFilesCall(SEXP filenameA, SEXP filenameB, SEXP minSimilarity, SEXP looseThreshold, SEXP tightThreshold, SEXP method) {
	SEXP result;

	PROTECT(filenameA = AS_CHARACTER(filenameA));
	PROTECT(filenameB = AS_CHARACTER(filenameB));
	PROTECT(minSimilarity = AS_NUMERIC(minSimilarity));
	PROTECT(looseThreshold = AS_NUMERIC(looseThreshold));
	PROTECT(tightThreshold = AS_NUMERIC(tightThreshold));
        PROTECT(method = AS_INTEGER(method));

	result = ccMatchFiles(
		CHAR(STRING_ELT(filenameA, 0)),
		CHAR(STRING_ELT(filenameB, 0)),
		REAL(minSimilarity)[0],
		REAL(looseThreshold)[0],
		REAL(tightThreshold)[0],
		INTEGER(method)[0]
	);
	
	UNPROTECT(6);

	return(result);
}


SEXP ccStatisticsCall() {
	SEXP result;

	result = ccStatistics();
  
	return(result);
}

// register wrapper-functions
void R_init_useCall(DllInfo *info) {
	R_CallMethodDef callMethods[]  = {
	  {"ccMatchFilesCall", (DL_FUNC) &ccMatchFilesCall, 4},
	  {"ccStatisticsCall", (DL_FUNC) &ccStatisticsCall, 0},
	  {NULL, NULL, 0}
	};
	
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

#ifdef __cplusplus
}
#endif
