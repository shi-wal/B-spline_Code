#pragma once
#include "BasicTypes.h"
#include "Pnt2D.h"
#include "Pnt3D.h"
#include <array>
#include <vector>

namespace KernelBridgeNS {

	inline SBool apprxEq(SReal x, 
						 SReal y)
	{
		return (fabs(x - y) < SEps);
	}



	inline SBool apprxEqEps(SReal x, 
							SReal y, 
							SReal Eps)
	{
		return (fabs(x - y) < Eps);
	}



	inline SReal sSign2(SReal a, 
						SReal b)
	{
		return (b >= 0.0 ? fabs(a) : -fabs(a));
	}



	template<class T>
	inline T boundMinMax(T x, 
						 T tMin, 
						 T tMax)
	{
		return max(min(x, tMax), tMin);
	}



	inline SReal randomReal(SReal min, 
							SReal max)
	{
		int
			r = rand();

		return min + (max - min) * ((SReal)(r & RAND_MAX)) / ((SReal)RAND_MAX);
	}



	inline SReal determinant(const Pnt2D & v0, 
							 const Pnt2D & v1)
	{
		return v0.coord(0) * v1.coord(1) - v0.coord(1) * v1.coord(0);
	}


	inline SReal determinant(const vector<Pnt2D> & vecs)
	{
		return determinant(vecs[0], vecs[1]);
	}


	inline SReal determinant(const Pnt3D & v0,
							 const Pnt3D & v1,
							 const Pnt3D & v2)
	{
		return v0[0] * v1[1] * v2[2] + v0[1] * v1[2] * v2[0] + v0[2] * v1[0] * v2[1] -
			   v2[0] * v1[1] * v0[2] - v2[1] * v1[2] * v0[0] - v2[2] * v1[0] * v0[1];
	}



	inline SReal determinant(const vector<Pnt3D> & vecs)
	{
		return determinant(vecs[0], vecs[1], vecs[2]);
	}



	inline void elimSimilarReals(vector<SReal> & tVals,
								 SReal eps)
	{
		SInt i = 0,
			nPts = (SInt)tVals.size();

		while (i < nPts) {
			SBool
				purgePt = SFalse;

			for (SInt j = i + 1; purgePt == SFalse && j < nPts; j++)
				if (abs(tVals[i] - tVals[j]) < eps)
					purgePt = STrue;

			if (purgePt) {
				tVals.erase(tVals.begin() + i);
				nPts--;
			}
			else
				i++;
		}
	}



	template<class T>
	inline void elimSimilarPts(vector<T> & pts, 
							   SReal eps)
	{
		SInt i = 0,
			nPts = (SInt)pts.size();

		while (i < nPts) {
			SBool
				purgePt = SFalse;

			for (SInt j = i + 1; purgePt == SFalse && j < nPts; j++) 
				if (pts[i].isEqual(pts[j], eps))
					purgePt = STrue;
			
			if (purgePt) {
				pts.erase(pts.begin() + i);
				nPts--;
			}
			else
				i++;
		}
	}
	



	template<class T>
	inline void elimSimilarPts2(vector<T> & pts,
								SReal eps)
	{
		SInt i = 0,
			nPts = (SInt)pts.size();

		while (i < nPts) {
			SBool
				purgePt = SFalse;

			for (SInt j = i + 1; purgePt == SFalse && j < nPts; j++)
				if (pts[i].isEqual2(pts[j], eps))
					purgePt = STrue;

			if (purgePt) {
				pts.erase(pts.begin() + i);
				nPts--;
			}
			else
				i++;
		}
	}
	   



	template<class T>
	SInt id2Indx(const vector<T> & items, 
				 SInt id)
	{
		SInt indx = 0;

		for (auto itr = items.begin(); itr != items.end(); itr++, indx++) {
			if ((*itr)._id == id)
				return indx;
		}
		return -1;
	}



	template<class T>
	typename vector<T>::const_iterator id2ItrBinary(const vector<T> & items, 
													SInt id)
	{
		auto itr = std::lower_bound(items.begin(), items.end(), id,
			[](const T & itm, SInt id) {
			return itm._id < id;
		});

		return itr;
	}



	template<class T>
	SInt id2IndxBinary(const vector<T> & items, 
					   SInt id)
	{
		auto itr = id2ItrBinary(items, id);

		if (itr == items.end())
			return -1;
		else
			return itr - items.begin();
	}



	template<class T>
	typename vector<T>::const_iterator id2Itr(const vector<T> & items, 
											  SInt id)
	{
		for (auto itr = items.begin(); itr != items.end(); itr++) {
			if ((*itr)._id == id)
				return itr;
		}
		return items.end();
	}

	

	template<class T>
	SInt id2Indx(const vector<T> & items, 
				 SInt id1, 
				 SInt id2)
	{
		SInt indx = 0;

		for (auto itr = items.begin(); itr != items.end(); itr++, indx++) {
			if ((*itr)._id1 == id1 && (*itr)._id2 == id2)
				return indx;
		}
		return -1;
	}



	template<class T>
	typename vector<T>::const_iterator id2Itr(const vector<T> & items, 
											  SInt id1, 
											  SInt id2)
	{
		for (auto itr = items.begin(); itr != items.end(); itr++) {
			if ((*itr)._id1 == id1 && (*itr)._id2 == id2)
				return itr;
		}
		return items.end();
	}



	/*	-1:	first segment
		0:	p1
		1:	second segment
	*/
	template<class PType>
	SInt closestSegment(const PType & p0,
						const PType & p1,
						const PType & p2,
						const PType & pt)
	{
		SReal
			a = (pt - p0).magnitudeL1() - (p1 - p0).magnitudeL1(),
			b = (pt - p2).magnitudeL1() - (p1 - p2).magnitudeL1();

		if (apprxEqEps(a, b, SEEps))
			return 0;
		else if (a < b)
			return -1;
		else
			return 1;
	}

	/* projPt = ptA + (ptB - ptA) * projPar */
	template<class PType>
	SBool projectPtOnLine(const PType & ptA,
						  const PType & ptB,
						  const PType & pt,
						  SReal & projPar,
						  PType & projPt)
	{
		PType
			diffBA = ptB - ptA,
			diffPA = pt - ptA;
		SReal
			magSqrBA = diffBA.magnitudeSquare();

		if (magSqrBA < SEEps) {
			projPar = -SInfnty;
			return SFalse;
		}

		projPar = diffPA.dot(diffBA) / magSqrBA;
		projPt = ptA + diffBA * projPar;
		return STrue;	
	}



	inline SBool linearCombination(const Pnt3D & a, 
								   const Pnt3D & b, 
								   const Pnt3D & c, 
								   Pnt2D & coeffs)
	{
		SReal det,
			eps = SEps * 10;

		det = a[0] * b[1] - a[1] * b[0];		

		if (abs(det) > eps) {
			coeffs.setCoord(0, (c[0] * b[1] - c[1] * b[0]) / det);
			coeffs.setCoord(1, (a[0] * c[1] - a[1] * c[0]) / det);
			return STrue;
		}

		det = a[1] * b[2] - a[2] * b[1];

		if (abs(det) > eps) {
			coeffs.setCoord(0, (c[1] * b[2] - c[2] * b[1]) / det);
			coeffs.setCoord(1, (a[1] * c[2] - a[2] * c[1]) / det);
			return STrue;
		}

		det = a[0] * b[2] - a[2] * b[0];

		if (abs(det) > eps) {
			coeffs.setCoord(0, (c[0] * b[2] - c[2] * b[0]) / det);
			coeffs.setCoord(1, (a[0] * c[2] - a[2] * c[0]) / det);
			return STrue;
		}

		return SFalse;
	}



	inline void intrvlUnite(array<SReal, 2> & a, 
						    const array<SReal, 2> & b)
	{
		a[0] = std::min(a[0], b[0]);
		a[1] = std::max(a[1], b[1]);
	}



	inline SInt lowerBndIndx(const vector<SReal> & tVals, 
							  SReal t)
	{
		vector<SReal>::const_iterator
			itr = std::lower_bound(tVals.begin(), tVals.end(), t);

		return (SInt)(itr - tVals.begin());
	}


	inline SBool realPresentInSortedVec(const vector<SReal> & rVals,
										SReal r,
										SReal tol)
	{
		SInt
			indx = lowerBndIndx(rVals, r - tol);

		if (indx >= rVals.size())
			return SFalse;

		if (std::abs(rVals[indx] - r) < tol || (indx > 0 && std::abs(rVals[indx - 1] - r) < tol))
			return STrue;
		else
			return SFalse;
	}


	inline SBool lineSegIntrsct2D(const Pnt2D & p0,
								  const Pnt2D & p1,
								  const Pnt2D & q0,
								  const Pnt2D & q1,
								  Pnt2D & par,
								  SReal eps)	
	{
		SBool
			segsIntersect = SFalse;
		Pnt2D
			r = p1 - p0,
			s = q1 - q0,
			pq = q0 - p0;
		SReal
			rxs = determinant(r, s),
			pqxr = determinant(pq, r);

		if (std::fabs(rxs) < eps) {	/* no transveral intersection */
			if (std::fabs(pqxr) < eps) { /* colinear line-segments */
				SReal
					rr = r.dot(r),
					sr = s.dot(r),
					t0 = pq.dot(r) / rr,
					t1 = t0 + sr / rr;

				if ((-eps < t0 && t0 < 1.0 + eps) || (-eps < t1 && t1 < 1.0 + eps)) { /* overlapping line-segments */
					par.fillAllCoord(0.5);
					segsIntersect = STrue;
				}
			}
		}
		else {	/* transveral intersection */
			SReal
				pqxs = determinant(pq, s),
				par0 = pqxs / rxs,
				par1 = pqxr / rxs;

			par.setAllCoord(par0, par1);
			if (-eps < par0 && par0 < 1.0 + eps && -eps < par1 && par1 < 1.0 + eps) /* line-segments intersect */
				segsIntersect = STrue;
		}

		return segsIntersect;
	}


	inline SInt pairOf(const vector<array<SInt, 2>> & pairIds,
					   SInt id)
	{
		for (SInt i = 0; i < pairIds.size(); i++) {
			if (pairIds[i][0] == id)
				return pairIds[i][1];
			else if (pairIds[i][1] == id)
				return pairIds[i][0];
		}
		return -1;
	}
	   

}


