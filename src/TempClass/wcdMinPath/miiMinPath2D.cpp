/**
 * @file miiMinPath2D.cpp
 * @brief
 * 
 *
 * @see miiMinPath2D
 *
 * @author L. L. 
 * @version 0.0.1
 * @date 2013/08/08
 */

/************************************
* upwind   
-------------------%%%%%++++++
-----------------%%%%%++++++++
---------------%%%%%++++++++++
-------------%%%%%++++++++++++
-----------%%%%%++++++++++++++
---------%%%%%++++++++++++++++
-------%%%%%++++++++++++++++++
-----%%%%%++++++++++++++++++++
* notes: 
* '-' --> 'alive', set 0
* '%' --> 'trial', set 1, Narrow Band
* '+' --> 'far', set 2
*
************************************/

/** include files **/
#include "miiMinPath2D.h"

/**
 * @brief miiMinPath2D()
 *
 * 
 * @author L. L.
 * 
 * @param nImgWX, nImgWY: the image size. 
 *
 * @return None.
 */
miiMinPath2D::miiMinPath2D(int nImgWX, int nImgWY):
					m_nImgWX(nImgWX), m_nImgWY(nImgWY)
{
	m_gfU = new double[nImgWX * nImgWY];
	m_nFmMap = new int[nImgWX * nImgWY];
	m_gfP2 = new double[nImgWX * nImgWY];

	// create a min-heap for Narrow Band
	m_iMinHeap = new miiMinHeap<>();

}

/**
 * @brief ~miiMinPath2D()
 *
 * 
 * @author L. L.
 * 
 * @param None. 
 *
 * @return None.
 */
miiMinPath2D::~miiMinPath2D()
{
	if (m_gfU != NULL)
		delete[] m_gfU;

	if (m_nFmMap != NULL)	
		delete[] m_nFmMap;

	if (m_gfP2 != NULL)
		delete[] m_gfP2;

	if (m_iMinHeap != NULL)
		delete m_iMinHeap;
}

/**
 * @brief FastMarchingInit()
 *
 * Initiate the FMM's data, including Narrow Band, distance function, and FMM map.
 * @author L. L.
 * 
 * @param *sImgData: the pointer for the data of 2D image.
 * @param gnStartPoint[]: the coordinate of the starting point.  
 * @param gnEndPoint[]: the coordinate of the end point.
 *						
 * @return None
 */
void miiMinPath2D::FastMarchingInitBase(short *gsImgData, int gnStartPoint[], int gnEndPoint[])
{
	if (m_iStartPoint.x >= m_nImgWX || m_iStartPoint.y >= m_nImgWY || 
		m_iEndPoint.x >= m_nImgWX || m_iEndPoint.y >= m_nImgWY)
	{
		cout << "The starting or end point is out of boundary!" << endl;
	}
	// set starting point
	m_iStartPoint.x = gnStartPoint[0]; 
	m_iStartPoint.y = gnStartPoint[1]; 
	m_iStartPoint.val = 0;

	// set end point
	m_iEndPoint.x = gnEndPoint[0];
	m_iEndPoint.y = gnEndPoint[1];
	m_iEndPoint.val = 0;

	// initiate the number of the iteration for FMM
	m_nFmItr = 0;

	// check narrow band
	if (m_vNarrowBand.size() > 0)
	{
		m_vNarrowBand.clear();
	}

	// add the starting point to Narrow Band
	m_vNarrowBand.push_back(m_iStartPoint);
	m_vNarrowBand.push_back(m_iStartPoint);

	// build the min-heap
	m_iMinHeap->BuildMinHeap(m_vNarrowBand);

	// initiate the distance function U and FMM map
	for (int i = 0; i < (m_nImgWX * m_nImgWY); i++)
	{
		m_gfU[i] = FMM_INF;
		m_nFmMap[i] = FMM_FAR;
	}

	// set the distance from starting point to starting point as zero 
	m_gfU[m_iStartPoint.x * m_nImgWY + m_iStartPoint.y] = 0;

	// set starting point as 'trail' in FMM map
	m_nFmMap[m_iStartPoint.x * m_nImgWY + m_iStartPoint.y] = FMM_TRIAL;

//	PotentialFunction(gsImgData, nMethod);
}

/**
 * @brief UpdateNarrowBandVal()
 *
 * Update the node value in the narrow band and re-organize the min-heap of the narrow band.
 * @author L. L.
 * 
 * @param iNode: the pixel needed to update in narrow band.
 *
 * @return Successful or failed.
 */
bool miiMinPath2D::UpdateNarrowBandVal(miiCNode<> iNode)
{
	for (int i = 1; i < m_vNarrowBand.size(); i++)
	{
		if (m_vNarrowBand[i].x == iNode.x && m_vNarrowBand[i].y == iNode.y)
		{
			if (iNode.val > m_vNarrowBand[i].val)
			{
				m_vNarrowBand[i].val = iNode.val;
				m_iMinHeap->MinHeapify(m_vNarrowBand, i);
			}
			else
			{
				m_iMinHeap->HeapDecreaseKey(m_vNarrowBand, i, iNode);
			}			

			return true;
		}		
	}

	return false;
}

/**
 * @brief UpWind()
 *
 * 
 * @author L. L.
 * 
 * @param x, y: the coordinate.
 *
 * @return None.
 */
void miiMinPath2D::UpWind(int x, int y)
{
	int u_l_x = x - 1, u_l_y = y;  
	int u_r_x = x + 1, u_r_y = y;
	int u_u_x = x,     u_u_y = y - 1;
	int u_d_x = x,     u_d_y = y + 1;
	int u_b_x = x,     u_b_y = y;
	int u_f_x = x,     u_f_y = y;
	double u_l, u_r, u_u, u_d;

	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (double)m_gfU[x * m_nImgWY + y];
	else
		u_l = (double)m_gfU[u_l_x * m_nImgWY + u_l_y];
	
	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (double)m_gfU[x * m_nImgWY + y];
	else
		u_r = (double)m_gfU[u_r_x * m_nImgWY + u_r_y];
	
	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (double)m_gfU[x * m_nImgWY + y];
	else
		u_u = (double)m_gfU[u_u_x * m_nImgWY + u_u_y];
	
	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (double)m_gfU[x * m_nImgWY + y];
	else
		u_d = (double)m_gfU[u_d_x * m_nImgWY + u_d_y];

	double a = min(u_l, u_r);
	double b = min(u_u, u_d);

	double dMax, dMin;
	if (a > b) 
	{
		dMax = a;
		dMin = b;
	}
	else
	{
		dMax = b;
		dMin = a;
	}

	a = dMax;
	b = dMin;

	double gfRoots[2] = {0};
	bool bSolv = QuadraticRoots(2, -2*(a+b), \
		a*a+b*b+m_gfP2[x * m_nImgWY + y], gfRoots);

	double fRtMax = max(gfRoots[0], gfRoots[1]);
	
	if (fRtMax > a && bSolv == true)
	{
		m_gfU[x * m_nImgWY + y] = fRtMax;
	}
	else
	{
		m_gfU[x * m_nImgWY + y] = \
			sqrt(m_gfP2[x * m_nImgWY + y]) + b;
	}
}

/**
 * @brief UpWind()
 *
 * 
 * @author L. L.
 * 
 * @param x, y: the coordinate.
 *
 * @return None.
 */
double miiMinPath2D::UpWind(double *gfU, int x, int y, double fP2)
{
	int u_l_x = x - 1, u_l_y = y;  
	int u_r_x = x + 1, u_r_y = y;
	int u_u_x = x,     u_u_y = y - 1;
	int u_d_x = x,     u_d_y = y + 1;
	int u_b_x = x,     u_b_y = y;
	int u_f_x = x,     u_f_y = y;
	double u_l, u_r, u_u, u_d;

	if (u_l_x < 0 || u_l_x >= m_nImgWX)
		u_l = (double)gfU[x * m_nImgWY + y];
	else
		u_l = (double)gfU[u_l_x * m_nImgWY + u_l_y];
	
	if (u_r_x < 0 || u_r_x >= m_nImgWX)
		u_r = (double)gfU[x * m_nImgWY + y];
	else
		u_r = (double)gfU[u_r_x * m_nImgWY + u_r_y];
	
	if (u_u_y < 0 || u_u_y >= m_nImgWY)
		u_u = (double)gfU[x * m_nImgWY + y];
	else
		u_u = (double)gfU[u_u_x * m_nImgWY + u_u_y];
	
	if (u_d_y < 0 || u_d_y >= m_nImgWY)
		u_d = (double)gfU[x * m_nImgWY + y];
	else
		u_d = (double)gfU[u_d_x * m_nImgWY + u_d_y];

	double a = min(u_l, u_r);
	double b = min(u_u, u_d);

	double dMax, dMin;
	if (a > b) 
	{
		dMax = a;
		dMin = b;
	}
	else
	{
		dMax = b;
		dMin = a;
	}

	a = dMax;
	b = dMin;

	double gfRoots[2] = {0};
	bool bSolv = QuadraticRoots(2, -2*(a+b), \
		a*a+b*b+fP2, gfRoots);

	double fRtMax = max(gfRoots[0], gfRoots[1]);
	
	double fU;
	if (fRtMax > a && bSolv == true)
	{
		fU = fRtMax;
	}
	else
	{
		fU = sqrt(m_gfP2[x * m_nImgWY + y]) + b;
	}

	return fU;	
}

/**
 * @brief QuadraticRoots()
 *
 * Calculate the solution of the quadratic function. 
 * @author L. L.
 * 
 * @param gfRoots[2]: two rooters. {a*(x^2) + b*x + c =0}		   
 *
 * @return Successful or failed.
 */
bool miiMinPath2D::QuadraticRoots(double a, double b, double c, double gfRoots[2])
{
	double delta = b * b - 4 * a * c;

	if (delta < 0)
	{
		return false;
	} 
	else
	{
		delta = sqrt(delta);
		gfRoots[0] = (-b + delta) / (2 * a + 0.0000001);
		gfRoots[1] = (-b - delta) / (2 * a + 0.0000001);
	}

	return true;
}

/**
 * @brief FindMinPath()
 *
 * Find the minimal path by using back-propagation in the distant map. 
 * @author L. L.
 * 
 * @param gnStartPoint[2]: the coordinate of the starting point.
 * @param nMaxPath: the maximum number of the iteration.
 *
 * @return Successful or failed.
 */
bool miiMinPath2D::FindMinPath(int gnStartPoint[2], int gnEndPoint[2], int nMaxPath)
{
	if (m_gfU == NULL || m_nFmItr == 0)
	{
		cout << "Cannot find minimal path!" << endl;
		return false;
	}

	// define a array for neighbor searching
	int gNbr[4][2] = { {-1, 0}, \
				       { 1, 0}, \
				       { 0,-1}, \
				       { 0, 1} };

	// define the end point as the source point
	int nMinX = gnEndPoint[0];
	int nMinY = gnEndPoint[1];

	int nCentX = gnEndPoint[0];
	int nCentY = gnEndPoint[0];

	int nx, ny;

	vector<miiCNode<>> vSgmtMinPath;

	if (m_vMinPath.size() > 0)
	{
		m_vMinPath.clear();
	}
	// save the source point
	miiCNode<> iTempPoint;
	iTempPoint.x = gnEndPoint[0];
	iTempPoint.y = gnEndPoint[1];
	iTempPoint.val = 0;
	vSgmtMinPath.push_back(iTempPoint);

	miiCNode<> iOrgStartPoint;

	int nPathNum = 0;

	// set starting point
	iOrgStartPoint.x = gnStartPoint[0]; 
	iOrgStartPoint.y = gnStartPoint[1]; 

	if (iOrgStartPoint.x >= m_nImgWX)
		iOrgStartPoint.x = m_nImgWX - 1;
	if (iOrgStartPoint.y >= m_nImgWY)
		iOrgStartPoint.y = m_nImgWY - 1;

	while (1)
	{
		nPathNum++;
		if (nPathNum > nMaxPath)
		{
			return false;
		}

		for (int i = 0; i < 4; i++)
		{
			nx = nCentX + gNbr[i][0];
			ny = nCentY + gNbr[i][1];

			// reach the starting point
			if (nx == iOrgStartPoint.x && ny == iOrgStartPoint.y)
			{
				iTempPoint.x = nx;
				iTempPoint.y = ny;
				iTempPoint.val = m_gfU[nx * m_nImgWY + ny];

				vSgmtMinPath.push_back(iTempPoint);

				// save last point (starting point)
				vSgmtMinPath.push_back(iOrgStartPoint);

				for (int j = (int)vSgmtMinPath.size() - 1; j >= 0; j--)
				{
					m_vMinPath.push_back(vSgmtMinPath[j]);
				}

				return true;
			}

			// prevent out of boundary
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY)
			{
				// find the minimal value around the center point 
				if (m_gfU[nx * m_nImgWY + ny] <= \
					m_gfU[nMinX * m_nImgWY + nMinY])
				{
					nMinX = nx;
					nMinY = ny;
				}
			}
		}

		nCentX = nMinX;
		nCentY = nMinY;

		// save the minimal value
		iTempPoint.x = nMinX;
		iTempPoint.y = nMinY;
		iTempPoint.val = m_gfU[nMinX * m_nImgWY + nMinY];

		vSgmtMinPath.push_back(iTempPoint);

		int nDist = (nMinX - iOrgStartPoint.x) * (nMinX - iOrgStartPoint.x) \
			+ (nMinY - iOrgStartPoint.y) * (nMinY - iOrgStartPoint.y);

		nDist = (int)(sqrt((double)nDist) + 0.5);

		if (nDist < 3)
		{
			// save last point (starting point)
			vSgmtMinPath.push_back(iOrgStartPoint);

			for (int j = (int)vSgmtMinPath.size() - 1; j >= 0; j--)
			{
				m_vMinPath.push_back(vSgmtMinPath[j]);
			}

			return true;
		}		
	}
}

/**
 * @brief FindMinPath()
 *
 * Find the minimal path by using back-propagation in the distant map. 
 * @author L. L.
 * 
 * @param gnStartPoint[2]: the coordinate of the starting point.
 * @param nMaxPath: the maximum number of the iteration.
 *
 * @return Successful or failed.
 */
bool miiMinPath2D::FindMinPath(int gnStartPoint[2], int nMaxPath)
{
	if (m_gfU == NULL || m_nFmItr == 0)
	{
		cout << "Cannot find minimal path!" << endl;
		return false;
	}

	// define a array for neighbor searching
	int gNbr[4][2] = { {-1, 0}, \
				       { 1, 0}, \
				       { 0,-1}, \
				       { 0, 1} };

	// define the end point as the source point
	int nMinX = m_iEndPoint.x;
	int nMinY = m_iEndPoint.y;

	int nCentX = m_iEndPoint.x;
	int nCentY = m_iEndPoint.y;

	int nx, ny;

	vector<miiCNode<>> vSgmtMinPath;

	if (m_vMinPath.size() > 0)
	{
		m_vMinPath.clear();
	}
	// save the source point
	vSgmtMinPath.push_back(m_iEndPoint);

	miiCNode<> iTempPoint, iOrgStartPoint;

	int nPathNum = 0;

	// set starting point
	iOrgStartPoint.x = gnStartPoint[0]; 
	iOrgStartPoint.y = gnStartPoint[1]; 

	if (iOrgStartPoint.x >= m_nImgWX)
		iOrgStartPoint.x = m_nImgWX - 1;
	if (iOrgStartPoint.y >= m_nImgWY)
		iOrgStartPoint.y = m_nImgWY - 1;

	while (1)
	{
		nPathNum++;
		if (nPathNum > nMaxPath)
		{
			return false;
		}

		for (int i = 0; i < 4; i++)
		{
			nx = nCentX + gNbr[i][0];
			ny = nCentY + gNbr[i][1];

			// reach the starting point
			if (nx == iOrgStartPoint.x && ny == iOrgStartPoint.y)
			{
				iTempPoint.x = nx;
				iTempPoint.y = ny;
				iTempPoint.val = m_gfU[nx * m_nImgWY + ny];

				vSgmtMinPath.push_back(iTempPoint);

				// save last point (starting point)
				vSgmtMinPath.push_back(iOrgStartPoint);

				for (int j = (int)vSgmtMinPath.size() - 1; j >= 0; j--)
				{
					m_vMinPath.push_back(vSgmtMinPath[j]);
				}

				return true;
			}

			// prevent out of boundary
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY)
			{
				// find the minimal value around the center point 
				if (m_gfU[nx * m_nImgWY + ny] <= \
					m_gfU[nMinX * m_nImgWY + nMinY])
				{
					nMinX = nx;
					nMinY = ny;
				}
			}
		}

		nCentX = nMinX;
		nCentY = nMinY;

		// save the minimal value
		iTempPoint.x = nMinX;
		iTempPoint.y = nMinY;
		iTempPoint.val = m_gfU[nMinX * m_nImgWY + nMinY];

		vSgmtMinPath.push_back(iTempPoint);

		int nDist = (nMinX - iOrgStartPoint.x) * (nMinX - iOrgStartPoint.x) \
			+ (nMinY - iOrgStartPoint.y) * (nMinY - iOrgStartPoint.y);

		nDist = (int)(sqrt((double)nDist) + 0.5);

		if (nDist < 3)
		{
			// save last point (starting point)
			vSgmtMinPath.push_back(iOrgStartPoint);

			for (int j = (int)vSgmtMinPath.size() - 1; j >= 0; j--)
			{
				m_vMinPath.push_back(vSgmtMinPath[j]);
			}

			return true;
		}		
	}
}

/**
 * @brief FindMinPath()
 *
 * Find the minimal path by using back-propagation in the distant map. 
 * @author L. L.
 * 
 * @param nMaxPath: the max length of the path for preventing overtime.
 *
 * @return Successful or failed.
 */
bool miiMinPath2D::FindMinPath(int nMaxPath)
{
	if (m_gfU == NULL || m_nFmItr == 0)
	{
		cout << "Cannot find minimal path!" << endl;
		return false;
	}

	// define a array for neighbor searching
	int gNbr[4][2] = { {-1, 0}, \
	                   { 1, 0}, \
	                   { 0,-1}, \
	                   { 0, 1}};

	// define the end point as the source point
	int nMinX = m_iEndPoint.x;
	int nMinY = m_iEndPoint.y;
	int nCentX = m_iEndPoint.x;
	int nCentY = m_iEndPoint.y;

	int nx, ny;

	vector<miiCNode<double>> vSgmtMinPath;

	// 	if (m_vMinPath.size() > 0)
	// 	{
	// 		m_vMinPath.clear();
	// 	}
	// save the source point
	vSgmtMinPath.push_back(m_iEndPoint);

	miiCNode<double> iTempPoint;
	int nPathNum = 0;

	while (1)
	{
		nPathNum++;
		if (nPathNum > nMaxPath)
		{
			return false;
		}

		for (int i = 0; i < 4; i++)
		{
			nx = nCentX + gNbr[i][0];
			ny = nCentY + gNbr[i][1];

			// reach the starting point
			if (nx == m_iStartPoint.x && ny == m_iStartPoint.y)
			{
				iTempPoint.x = nx;
				iTempPoint.y = ny;
				iTempPoint.val = m_gfU[nx * m_nImgWY + ny];

				vSgmtMinPath.push_back(iTempPoint);

				// save last point (starting point)
				vSgmtMinPath.push_back(m_iStartPoint);

				for (int j = (int)vSgmtMinPath.size() - 1; j >= 0; j--)
				{
					m_vMinPath.push_back(vSgmtMinPath[j]);
				}

				return true;
			}

			// prevent out of boundary
			if (nx >= 0 && nx < m_nImgWX && ny >= 0 && ny < m_nImgWY)
			{
				// find the minimal value around the center point 
				if (m_gfU[nx * m_nImgWY + ny] <= \
					m_gfU[nMinX * m_nImgWY + nMinY])
				{
					nMinX = nx;
					nMinY = ny;
				}
			}
		}

		nCentX = nMinX;
		nCentY = nMinY;

		// save the minimal value
		iTempPoint.x = nMinX;
		iTempPoint.y = nMinY;
		iTempPoint.val = m_gfU[nMinX * m_nImgWY + nMinY];

		vSgmtMinPath.push_back(iTempPoint);

		int nDist = (nMinX - m_iStartPoint.x) * (nMinX - m_iStartPoint.x) \
			+ (nMinY - m_iStartPoint.y) * (nMinY - m_iStartPoint.y);

		nDist = (int)(sqrt((double)nDist) + 0.5);

		if (nDist < 3)
		{
			// save last point (starting point)
			vSgmtMinPath.push_back(m_iStartPoint);

			for (int j = (int)vSgmtMinPath.size() - 1; j >= 0; j--)
			{
				m_vMinPath.push_back(vSgmtMinPath[j]);
			}

			return true;
		}		
	}
}

/**
 * @brief GenMinPathGraph()
 *
 * Generate the result image by overlapping the minimal path on the image.
 * @author L. L.
 * 
 * @param *gsImgData: the original image data.
 * @param *gsOutData: the result image data.
 * @param nMethod: the result image type.
 *                 '0' - "ImageLine"
 *				   '1' - "Line"
 *
 * @return Successful or failed.
 */
bool miiMinPath2D::GenMinPathGraph(short *gsImgData, short *gsOutImgData, int nMethod)
{
	if (m_vMinPath.size() == 0)
	{
		cout << "The path is null in the stage of the saving image!" << endl;
		return false;
	}

	if (gsOutImgData == NULL)
	{
		cout << "The pointer of the result image is null!" << endl;
		return false;
	}
	
	int gnCord[2] = {0};

	for (int iy = 0; iy < m_nImgWY; iy++)
	{
		for (int ix = 0; ix < m_nImgWX; ix++)
		{
			if (nMethod == 0)
			{
				gnCord[0] = ix; 
				gnCord[1] = iy; 

				if (gnCord[0] >= m_nImgWX)
					gnCord[0] = m_nImgWX - 1;
				if (gnCord[1] >= m_nImgWY)
					gnCord[1] = m_nImgWY - 1;

				gsOutImgData[ix * m_nImgWY + iy] = \
						gsImgData[gnCord[0] * m_nImgWY + gnCord[1]];
			}
			else
			{
				gsOutImgData[ix * m_nImgWY+ iy] = 0; 
			}	
		}
	}


	for (int i = 0; i < m_vMinPath.size(); i++)
	{
		gsOutImgData[m_vMinPath[i].x * m_nImgWY + m_vMinPath[i].y] = 1000;
	}

	// define a array for neighbor searching
	int gNbr[4][2] = { {-1, 0}, \
	                   { 1, 0}, \
	                   { 0,-1}, \
	                   { 0, 1} };

	// label the starting and end points
	for (int i = 0; i < 4; i++)
	{
		int nx_1 = m_vMinPath[0].x + gNbr[i][0];
		int ny_1 = m_vMinPath[0].y + gNbr[i][1];

		int nx_2 = m_vMinPath[m_vMinPath.size()-1].x + gNbr[i][0];
		int ny_2 = m_vMinPath[m_vMinPath.size()-1].y + gNbr[i][1];

		gsOutImgData[nx_1 * m_nImgWY + ny_1] = 2000;

		gsOutImgData[nx_2 * m_nImgWY + ny_2] = 2000;
	}

	// save the image data
//	SaveImage(m_pBaseImgInfo, sTrsmData, strFileName);

	return true;
}

/* --------------------------------- ENDING LINE ------------------------------------- */
