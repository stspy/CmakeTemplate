/**
 * @file miiMinPath2D.h
 * @brief
 * 
 *
 * @see miiMinPath2D
 *
 * @author L. L. 
 * @version 0.0.1
 * @date 2013/08/08
 */

#ifndef _miiMinPath2D_h
#define _miiMinPath2D_h

// the definition of fast-marching-method
#define FMM_ALIVE 0
#define FMM_TRIAL 1
#define FMM_FAR 2

// define a value for infinity
#define FMM_INF 800000000

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>


#include "miiMinHeap.h"

using namespace std;

/** 
 * @class miiMinPath2D
 * @brief The parent class for 2D minimal path. 
 * 
 *
 * @ingroup miiMinPath2D
 */
class miiMinPath2D
{
public:
	///
	miiMinPath2D(){};
	miiMinPath2D(int nImgWX, int nImgWY);

	///
	virtual ~miiMinPath2D();	 

	/// 
	virtual int FastMarchingEvolution(int nMaxItrNum) = 0;

	///
	bool FindMinPath(int gnStartPoint[2], int gnEndPoint[2], int nMaxPath = 2000);

	///
	bool FindMinPath(int gnStartPoint[2], int nMaxPath = 2000);

	///
	bool FindMinPath(int nMaxPath = 2000);


	/// generate the result image
	bool GenMinPathGraph(short *gsImgData, short *gsOutImgData, int nMethod = 0);

protected:
	/// the image size
	int m_nImgWX, m_nImgWY;

	/// define the starting and end point of the coronary 
	miiCNode<> m_iStartPoint, m_iEndPoint;

	/// define the distance function U in Upwind
	double *m_gfU;

	/// define the map of the fast-marching-method for 
	/// saving 'alive', 'trial(Narrow Band)', and 'far'
	int *m_nFmMap;

	/// potential of minimal path
	double *m_gfP2;

	/// define the min-heap for NarrowBand
	miiMinHeap<> *m_iMinHeap;

	/// define the vector of Narrow Band for 'trail' in fast-marching-method
	vector<miiCNode<>> m_vNarrowBand;

	/// define a vector for saving the points of minimal path
	vector<miiCNode<>> m_vMinPath;

	/// define the variable for recording the number of FMM iteration
	int m_nFmItr;

	//// the initiation of fast-marching-method(FMM)
	void FastMarchingInitBase(short *gsImgData, int gnStartPoint[], int gnEndPoint[]);
	
	/// 
	void UpWind(int x, int y);

	///
	double UpWind(double *gfU, int x, int y, double fP2);

	bool UpdateNarrowBandVal(miiCNode<>);

	/// solve the quadratic function
	bool QuadraticRoots(double a, double b, double c, double gfRoots[2]);	


};

#endif

/* --------------------------------- ENDING LINE ------------------------------------- */
