#pragma once

#include "BasicTypes.h"
#include "Pnt1D.h"
#include "Pnt2D.h"
#include "Pnt3D.h"
#include "Pnt4D.h"
#include "Pnt5D.h"
#include "Pnt6D.h"
#include <vector>
#include <array>
#include <assert.h>
#include "Utils.h"

/*  TODO:
	replace asserts with exceptions
	optimizations in evaluation, making basisVec, indexVec, indexFirst as members and allocating only once in constructor */

using namespace KernelBridgeNS;
using namespace std;

typedef enum {			  /* List of all possible tokens enumerated. */
	TOKEN_NONE,

	TOKEN_OPEN_PAREN,
	TOKEN_CLOSE_PAREN,

	TOKEN_E1,
	TOKEN_E2,
	TOKEN_E3,
	TOKEN_E4,
	TOKEN_E5,
	TOKEN_E6,
	TOKEN_E7,
	TOKEN_E8,
	TOKEN_E9,

	IP_TOKEN_OBJECT,
	TOKEN_BEZIER,
	TOKEN_BSPLINE,
	TOKEN_PTYPE,
	TOKEN_NUM_PTS,
	TOKEN_ORDER,
	TOKEN_KV,
	TOKEN_MULTIVAR,

	TOKEN_OTHER = 100,		/* Probably names & numbers. */
	TOKEN_QUOTED,			/* A quoted string. */

	TOKEN_EOF = -1
} ITDTokenType;


namespace SweepNS {
	template <class PTypeD , class PTypeR>
	class BspMap
	{
	private:
		template<class PTypeD, class PTypeR>
		friend class BspMap;

		SInt _domainDim;
		SInt _rangeDim;
		vector<SInt> _orders;
		vector<SInt> _lengths;
		vector<SInt> _subSpaces;
		vector<vector<SReal>> _knotVectors;
		vector<PTypeR> _points;
		SInt _id1;
		SInt _id2;

		static const SInt _nCkMaxOrder = 99;
		static const SInt _nCkMaxOrder2 = 15;
		static const array<array<SReal, _nCkMaxOrder + 1>, _nCkMaxOrder + 1> _NCKTable;
		static const array<array<array<array<SReal, _nCkMaxOrder2 + 1>, _nCkMaxOrder2 + 1>, _nCkMaxOrder2 + 1>, _nCkMaxOrder2 + 1> _IcKJcMIJcKMTable;
		static const SReal _Eps;
		static const SReal _dmnEps;
		static const SReal _numGradEps;
		static const SReal _epsRoundKnot;

		struct AlphaCoeff {
			SInt _order;
			SInt _length;
			SInt _refLength;
			vector<vector<SReal>> _mat;
			vector<vector<SReal>> _matTransp;
			vector<SInt> _colInd;
			vector<SInt> _colLen;

			AlphaCoeff();
			AlphaCoeff(SInt order, SInt lenKVT, SInt lenKVt);
		};
		AlphaCoeff evalAlphaCoeff(SInt order, SInt lenKVT, SInt lenKVt, const vector<SReal> & KVT, const vector<SReal> & KVt)	const;
		void alphaBlendStep(AlphaCoeff & alphaC, SInt iMin, SInt iMax, typename vector<PTypeR>::const_iterator & origPts, SInt origPtsStep, vector<PTypeR> & refPtsVec, SInt refPtsSt, SInt refPtsStep)    const;

		struct BlsmAlphaCoeff {
			SInt _order;
			SInt _length;
			SInt _newOrder;
			SInt _newLength;
			vector<vector<SReal>> _mat;
			vector<vector<SReal>> _matTransp;
			vector<SInt> _colInd;
			vector<SInt> _colLen;
			vector<SReal> _KV;
			vector<SReal> _newKV;

			BlsmAlphaCoeff() {}
			BlsmAlphaCoeff(SInt order, SInt length, SInt newOrder, SInt newLength);
			void addRow(const vector<SReal> & coefs, SInt aRow, SInt colIndex, SInt colLength);
			void scale(SReal scl);
			void setDomain();
		};
		struct BlsmSymb {
			vector<SReal> _coefs;
			SInt _min;
			SInt _max;

			BlsmSymb()
			{
			}
			BlsmSymb(SInt len) : _coefs(len)
			{
			}
		};
		struct BlsmEvalCache {
			SInt _idxFirst;
			SInt _idxLength;
			SInt _order;
			vector<SReal> _coefs;
			vector<SReal> _blsmVals;
			vector<SReal> _knots;
			vector<BlsmSymb> _symbVec;
		};

		BlsmAlphaCoeff blsmDegRaiseMat(const vector<SReal> & KV, SInt order, SInt len)	const;
		BlsmAlphaCoeff blsmDegRaiseNMat(const vector<SReal> & KV, SInt order, SInt newOrder, SInt len)	const;
		BlsmAlphaCoeff blsmDegRaiseMatProd(const BlsmAlphaCoeff &a1, const BlsmAlphaCoeff &a2)	const;
		void blsmEvalSymb(SInt order, const vector<SReal> & knots, const vector<SReal> & blsmVals, SInt blsmLen, SInt & retIdxFirst, SInt & retLength, BlsmEvalCache & blsmCache)   const;
		SReal validateMinMaxDomain(SReal t, SReal tMin, SReal tMax) const;

		SReal nChooseK(SInt n, SInt k)	const;
		SBool paramInDomain(SInt dim, SReal t)	const;
		vector<SReal> coxDeBoorBasis(SInt dim, SReal t, SInt &indexFirst)   const;
		void bezierBasis(SInt order, SReal t, vector<SReal> & basisVec)	const;

		BspMap bspMultiplyAux(const BspMap &mv)	const;

		SBool hasBezierKnotVector(const vector<SReal> & knotVector, SInt order, SInt length)	const;
		SInt knotsLastIndexLE(const vector<SReal> & knotVector, SReal t)    const;
		SInt knotsLastIndexL(const vector<SReal> & knotVector, SReal t)	    const;
		SInt knotsFirstIndexG(const vector<SReal> & knotVector, SReal t)    const;
		vector<SReal> knotsDegreeRaise(const vector<SReal> & kv, SInt len, SInt order, SInt newOrder)	const;
		SBool knotC0Discont(const vector<SReal> & kv, SInt order, SInt length, SReal & t)   const;
		SBool knotC1Discont(const vector<SReal> & kv, SInt order, SInt length, SReal & t)   const;
		void uniformOpenKnotVector(SInt dir);
		SInt interiorKnots(SReal & knot)    const;
		vector<SReal> mergeKnotVectors(const vector<SReal> & knotVector1, const vector<SReal> & knotVector2, SInt mult)	const;
		SBool knotsMakeRobust(vector<SReal> &kv);
		vector<SReal> knotSubtract(const vector<SReal> & kv1, const vector<SReal> & kv2);
		void knotAffineTransform(vector<SReal> & kv, SReal translate, SReal scale)  const;
		void knotAffineTransform2(vector<SReal> & kv, SReal minVal, SReal maxVal)  const;
		void knotAffineTransformOrder2(vector<SReal> & kv, SInt order, SReal minVal, SReal maxVal)	const;
		SBool discontKnot(const vector<SReal> &kv, SReal t, SInt order)	const;
		SBool knotsHasOpenEnd(const vector<SReal> &kv, SInt len, SInt order)	const;
		void reverseKnotVec(SInt dir);

		void extractBspMap(const BspMap<PTypeD,PTypeR> & mv, SReal t, SInt dir);
		void extendBspMap(const BspMap<PTypeD,PTypeR> & mv, SInt axis);
		void shiftAxes(SInt axis);
		PTypeR evaluateFromPoints(SInt inc, SInt order, SInt length, const vector<SReal> & basisVec, SInt ptsInd, SInt indexFirst) const;
		SInt linearizePointsIndex(const vector<SInt> & indexVec)    const;
		SInt incremetPointsIndex(vector<SInt> & indexVec, SInt & index)	const;
		SInt incremetPointsOrderIndex(vector<SInt> & indexVec, SInt & index)	const;
		SInt incremetPointsBoundIndex(vector<SInt> & indexVec, vector<SInt> & lowerBnd, vector<SInt> & upperBnd, SInt & index)	const;
		SInt incrementPointsSkipIndex(vector<SInt> & indexVec, SInt skipDir, SInt & index)	const;
		SBool verifyTDomain(SReal tMin, SReal tMax, SReal t)	const;
		void subdivideBezierCtrlMesh(BspMap & lmv, BspMap & rmv, SReal t, SInt dir);
		void subdivideBezierCtrlMeshOneSide(BspMap & lmv, BspMap & rmv, SReal t, SInt dir, SBool leftSide)  const;


		int getRealNumber(const string & strNum, double & realNum);
		SBool getStringToken(ifstream & file, string & fileBuffer, SInt & bufferInd, string & stringToken, SBool & quoted);
		ITDTokenType getToken(ifstream & file, string & fileBuffer, SInt & bufferInd, string & stringToken);

	public:
		BspMap();		

		BspMap(
			const vector<SInt> &orders,
			const vector<SInt> &lengths);

		BspMap(
			const vector<SInt> &orders,
			const vector<SInt> &lengths,
			const vector<vector<SReal>> &knotVectors);

		BspMap(
			const vector<SInt> &orders,
			const vector<SInt> &lengths,
			const vector<PTypeR> &points,
			SInt id1 = -1,
			SInt id2 = -1);

		BspMap(
			const vector<SInt> &orders,
			const vector<SInt> &lengths,
			const vector<vector<SReal>> &knotVectors,
			const vector<PTypeR> &points);

		BspMap(const string & FileName);

		BspMap(ifstream & inFS);

		BspMap(const BspMap &mv);
		
		BspMap(const BspMap<PTypeD,PTypeR> & mv, SReal t, SInt dir);

		template <class PTypeD2>
		BspMap(const BspMap<PTypeD2,PTypeR> & mv, SReal t, SInt dir);
		
		BspMap(const BspMap<PTypeD,PTypeR> &mv, SInt axis);
		
		template <class PTypeD2>
		BspMap(const BspMap<PTypeD2,PTypeR> &mv, SInt axis);

		BspMap(const BspMap<PTypeD,PTypeR> &mv, SInt newDim, SInt startAxis);

		BspMap(const vector<BspMap<Pnt1D, Pnt1D>> &scalarMVs);

		BspMap & operator = (const BspMap &mv);

		SInt rangeDim()	const;

		SInt domainDim()    const;

		void domain(SReal & tMin, SReal & tMax, SInt axis)    const;

		const vector<PTypeR> & controlPoints()	const;

		void knotsMultiplicities(SInt dir, vector<SReal> & knotVals, vector<SInt> & multiplicities, SReal eps)	const;

		void domain(PTypeD & tMin, PTypeD & tMax)	const;

		void setDomain(SReal tMin, SReal tMax, SInt axis);

		SBool lengthEqOrder(SInt axis)	const;

		SReal midKnotVal(SInt axis) const;

		void transform(PTypeR translate, SReal scale);

		const vector<SInt> & orders()	const;

		const vector<SInt> & lengths()	const;

		const vector<vector<SReal>> & knotVectors()	const;

		PTypeR evaluate(const PTypeD & params)	const;

		BspMap derive(SInt dir)	const;

		void boundingBox(PTypeR & min, PTypeR & max)	const;

		SBool domainCheck(const BspMap &mv)	const;

		void subdivide(SReal t, SInt dir, BspMap<PTypeD,PTypeR> & lmv, BspMap<PTypeD,PTypeR> & rmv)   const;

		BspMap subdivideOneSide(SReal t, SInt dir, SBool leftSide)   const;

		BspMap extractRegion(SReal t1, SReal t2, SInt dir)	const;

		BspMap extractRegion(const PTypeD & par1, const PTypeD & par2)	const;

		BspMap degreeRaise(vector<SInt> & newOrders)  const;

		BspMap degreeRaise2(vector<SInt> & newOrders)  const;

		BspMap bezierMultiply(const BspMap &mv) const;

		BspMap refineAtParams(SInt dir, SBool replace, vector<SReal> & knots)  const;

		SBool makeCompatible(BspMap & mv, SBool sameOrders, SBool sameKVs);

		SBool makeOneDimCompatible(BspMap & mv, SInt dim, SBool sameOrders, SBool sameKVs);

		BspMap merge(const BspMap & mv, SInt dir, SBool discont)	const;

		vector<BspMap<Pnt1D, Pnt1D>> splitScalar()  const;

		BspMap operator + (const BspMap &mv)	const;

		BspMap operator + (const PTypeR & pt)	const;

		BspMap operator - (const BspMap &mv)	const;

		BspMap operator * (const BspMap &mv)    const;

		BspMap reverseDir(SInt axis)	const;

		void preconditionScale(SReal scale);

		SReal maxScalarPtMagnitude()	const;

		BspMap scale(SReal scale)	const;
		
		BspMap<Pnt1D, Pnt1D> dotProd(const BspMap &mv)	const;

		BspMap<Pnt1D, Pnt1D> dotProd(const PTypeR & pt)	const;

		BspMap crossProd(const BspMap &mv)	const;

		BspMap scalarProd(const BspMap<Pnt1D, Pnt1D> &mv)	const;

		SBool noZeroCrossing(SReal eps)	const;

		SBool noZeroCrossingSub(const BspMap &mv, SReal eps)	const;

		SBool allPtNegative()	const;

		SBool allPtNegativeSub(const BspMap &mv)    const;

		SBool allPtPositive()	const;

		SBool allPtPositiveSub(const BspMap &mv)    const;

		PTypeD evalGradientNumeric(const PTypeD & param)    const;

		vector<Pnt2D> deriveAllBounds()	const;

		SBool knotC1Discont(SInt dir, SReal & t)    const;

		void perturb(SReal eps);

		void readFromStream(ifstream & inFS);
		void writeToFile(const string & FileName, SBool append)	const;
		void writeToStream(ofstream & outFS)	const;
	};
}