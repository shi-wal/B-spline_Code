#include "BspMap.h"
#include "Utils.h"
#include <assert.h>

using namespace KernelBridgeNS;
using namespace SweepNS;
using namespace std;



template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR> BspMap<PTypeD, PTypeR>::refineAtParams(SInt dir,
	SBool replace,
	vector<SReal> & knots)	const
{
	SInt lengthKVt,
		length = _lengths[dir],
		order = _orders[dir],
		n = (SInt)knots.size();

	BspMap<PTypeD, PTypeR> refMV;

	assert(dir >= 0 && dir < _domainDim);

	if (replace) {
		for (SInt i = 1; i < n; i++)
			assert(knots[i] >= knots[i - 1]);

		assert(order + length == n);

		refMV = *this;

		refMV._knotVectors[dir] = knots;
	}
	else if (n == 0) {
		refMV = *this;
	}
	else {
		vector<SInt>
			lenVec = _lengths;

		lenVec[dir] += n;
		refMV = BspMap<PTypeD, PTypeR>( _orders, lenVec);

		vector<SReal>
			mergedKV = mergeKnotVectors(_knotVectors[dir], knots, 0);

		lengthKVt = (SInt)mergedKV.size();

		AlphaCoeff
			A = evalAlphaCoeff(order, length, lengthKVt - order, _knotVectors[dir], mergedKV);

		for (SInt i = 0; i < _domainDim; i++) {
			if (i == dir)
				refMV._knotVectors[i] = mergedKV;
			else
				refMV._knotVectors[i] = _knotVectors[i];
		}

		SInt
			rIndex = 0,
			rLength = refMV._lengths[dir];
		vector<SInt>
			indices(_domainDim, 0);

		do {
			SInt
				index = linearizePointsIndex(indices),
				rStep = refMV._subSpaces[dir],
				step = _subSpaces[dir];
			typename vector<PTypeR>::const_iterator
				pts = _points.begin() + index;

			alphaBlendStep(A, 0, rLength, pts, step, refMV._points, rIndex, rStep);
		} while (refMV.incrementPointsSkipIndex(indices, dir, rIndex));
	}

	for (SInt i = 0; i < _domainDim; i++)
		refMV.knotsMakeRobust(refMV._knotVectors[i]);

	return refMV;
}


template <class PTypeD, class PTypeR>
SBool BspMap<PTypeD, PTypeR>::makeCompatible(BspMap<PTypeD, PTypeR> & mv,
	SBool sameOrders,
	SBool sameKVs)
{
	if (sameOrders) {
		SBool doRaise;
		SInt l;
		vector<SInt> newOrders(_domainDim);

		for (l = 0; l < _domainDim; l++)
			newOrders[l] = max(_orders[l], mv._orders[l]);

		for (l = 0, doRaise = SFalse; l < _domainDim; l++) {
			if (_orders[l] != newOrders[l])
				doRaise = STrue;
		}

		if (doRaise)
			*this = degreeRaise(newOrders);

		for (l = 0, doRaise = SFalse; l < _domainDim; l++) {
			if (mv._orders[l] != newOrders[l])
				doRaise = STrue;
		}

		if (doRaise)
			mv = mv.degreeRaise(newOrders);

		if (sameKVs) {
			for (SInt l = 0; l < _domainDim; l++) {
				vector<SReal>
					& kv1 = _knotVectors[l],
					&kv2 = mv._knotVectors[l];
				SInt
					order = _orders[l],
					kv1Len = (SInt)kv1.size(),
					kv2Len = (SInt)kv2.size();

				knotAffineTransform(kv2, kv1[order - 1] - kv2[order - 1], (kv1[kv1Len - order] - kv1[order - 1]) / (kv2[kv2Len - order] - kv2[order - 1]));

				vector<SReal>
					kv1Tmp(kv1.begin() + order, kv1.end() - order + 1),
					kv2Tmp(kv2.begin() + order, kv2.end() - order),
					refKV = knotSubtract(kv2Tmp, kv1Tmp);

				SInt
					refLen = (SInt)refKV.size();

				if (refLen > 0) {
					*this = refineAtParams(l, SFalse, refKV);
					kv1 = _knotVectors[l];
					kv1Len = (SInt)kv1.size();
				}
				kv1Tmp = vector<SReal>(kv1.begin() + order, kv1.end() - order);
				kv2Tmp = vector<SReal>(kv2.begin() + order, kv2.end() - order + 1);
				refKV = knotSubtract(kv1Tmp, kv2Tmp);
				refLen = (SInt)refKV.size();

				if (refLen > 0) {
					mv = mv.refineAtParams(l, SFalse, refKV);
				}
			}
		}
	}

	return STrue;
}



template <class PTypeD, class PTypeR>
SBool BspMap<PTypeD, PTypeR>::makeOneDimCompatible(BspMap<PTypeD, PTypeR> & mv,
	SInt dim,
	SBool sameOrders,
	SBool sameKVs)
{
	if (sameOrders && _orders[dim] != mv._orders[dim]) {
		vector<SInt> newOrders(_domainDim);

		if (_orders[dim] < mv._orders[dim]) {
			newOrders = _orders;
			newOrders[dim] = mv._orders[dim];
			*this = degreeRaise(newOrders);
		}
		else {
			newOrders = mv._orders;
			newOrders[dim] = _orders[dim];
			mv = mv.degreeRaise(newOrders);
		}
	}

	if (sameOrders && sameKVs) {
		vector<SReal>
			&kv1 = _knotVectors[dim],
			&kv2 = mv._knotVectors[dim];
		SInt
			order = _orders[dim],
			kv1Len = (SInt)kv1.size(),
			kv2Len = (SInt)kv2.size();

		knotAffineTransform(kv2, kv1[order - 1] - kv2[order - 1], (kv1[kv1Len - order] - kv1[order - 1]) / (kv2[kv2Len - order] - kv2[order - 1]));

		vector<SReal>
			kv1Tmp(kv1.begin() + order, kv1.end() - order),
			kv2Tmp(kv2.begin() + order, kv2.end() - order),
			refKV = knotSubtract(kv2Tmp, kv1Tmp);

		SInt
			refLen = (SInt)refKV.size();

		if (refLen > 0) {
			*this = refineAtParams(dim, SFalse, refKV);
			kv1 = _knotVectors[dim];
			kv1Len = (SInt)kv1.size();
			kv1Tmp = vector<SReal>(kv1.begin() + order, kv1.end() - order);
		}

		refKV = knotSubtract(kv1Tmp, kv2Tmp);
		refLen = (SInt)refKV.size();

		if (refLen > 0) {
			mv = mv.refineAtParams(dim, SFalse, refKV);
		}
	}

	return STrue;
}




template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR> BspMap<PTypeD, PTypeR>::merge(const BspMap<PTypeD, PTypeR> & mv,
	SInt dir,
	SBool discont)	const
{
	SBool
		mergedEdgeDiscont = discont || _orders[dir] == 1 || mv._orders[dir] == 1;
	SReal min1, max1, min2, max2;

	domain(min1, max1, dir);
	mv.domain(min2, max2, dir);

	BspMap<PTypeD, PTypeR>
		cpMV1 = *this,
		cpMV2 = mv,
		&mv1 = min1 > min2 ? cpMV2 : cpMV1,
		&mv2 = min1 > min2 ? cpMV1 : cpMV2;

	vector<SInt>
		lengths = mv1._lengths;

	for (SInt i = 0; i < _domainDim; i++) {
		if (i == dir)
			lengths[i] = mv1._lengths[i] + mv2._lengths[i] - 1 + mergedEdgeDiscont;
		else if (mv1._lengths[i] != mv2._lengths[i])
			mv1.makeOneDimCompatible(mv2, i, STrue, STrue);
	}

	BspMap<PTypeD, PTypeR>
		mergedMV( mv1._orders, lengths);

	for (SInt i = 0; i < _domainDim; i++) {
		if (i == dir) {
			vector<SReal> tmpKV(mv2._knotVectors[i].begin() + mv2._orders[i], mv2._knotVectors[i].end());

			knotAffineTransform(tmpKV, mv1._knotVectors[i][mv1._lengths[i] + mv1._orders[i] - 2 + mergedEdgeDiscont] - mv2._knotVectors[i][0], 1.0);
			std::copy(mv1._knotVectors[i].begin(), mv1._knotVectors[i].end(), mergedMV._knotVectors[i].begin());
			std::copy(tmpKV.begin(), tmpKV.end(), mergedMV._knotVectors[i].begin() + (mv1._lengths[i] + mv1._orders[i] - 1 + mergedEdgeDiscont));
		}
		else {
			mergedMV._knotVectors[i] = mv1._knotVectors[i];
		}
	}

	vector<SInt> mergedIndices(_domainDim, 0), lowerBound(_domainDim, 0),
		upperBound = mergedMV._lengths;
	SInt
		index1 = 0,
		mergedIndex = 0;

	upperBound[dir] = mv1._lengths[dir];

	do {
		mergedMV._points[mergedIndex] = mv1._points[index1];
		mergedMV.incremetPointsBoundIndex(mergedIndices, lowerBound, upperBound, mergedIndex);
	} while (++index1 < mv1._subSpaces[_domainDim]);

	lowerBound[dir] = mv1._lengths[dir] - 1 + mergedEdgeDiscont;
	upperBound[dir] = mergedMV._lengths[dir];

	SInt
		index2 = 0;

	mergedIndices = lowerBound;
	mergedIndex = mergedMV.linearizePointsIndex(mergedIndices);

	do {
		mergedMV._points[mergedIndex] = mv2._points[index2];
		mergedMV.incremetPointsBoundIndex(mergedIndices, lowerBound, upperBound, mergedIndex);
	} while (++index2 < mv2._subSpaces[_domainDim]);

	return mergedMV;
}


template class BspMap<Pnt1D, Pnt1D>;
template class BspMap<Pnt1D, Pnt2D>;
template class BspMap<Pnt1D, Pnt3D>;
template class BspMap<Pnt1D, Pnt4D>;
template class BspMap<Pnt1D, Pnt5D>;
template class BspMap<Pnt1D, Pnt6D>;

template class BspMap<Pnt2D, Pnt1D>;
template class BspMap<Pnt2D, Pnt2D>;
template class BspMap<Pnt2D, Pnt3D>;
template class BspMap<Pnt2D, Pnt4D>;
template class BspMap<Pnt2D, Pnt5D>;
template class BspMap<Pnt2D, Pnt6D>;

template class BspMap<Pnt3D, Pnt1D>;
template class BspMap<Pnt3D, Pnt2D>;
template class BspMap<Pnt3D, Pnt3D>;
template class BspMap<Pnt3D, Pnt4D>;
template class BspMap<Pnt3D, Pnt5D>;
template class BspMap<Pnt3D, Pnt6D>;

template class BspMap<Pnt4D, Pnt1D>;
template class BspMap<Pnt4D, Pnt2D>;
template class BspMap<Pnt4D, Pnt3D>;
template class BspMap<Pnt4D, Pnt4D>;
template class BspMap<Pnt4D, Pnt5D>;
template class BspMap<Pnt4D, Pnt6D>;

template class BspMap<Pnt5D, Pnt1D>;
template class BspMap<Pnt5D, Pnt2D>;
template class BspMap<Pnt5D, Pnt3D>;
template class BspMap<Pnt5D, Pnt4D>;
template class BspMap<Pnt5D, Pnt5D>;
template class BspMap<Pnt5D, Pnt6D>;

template class BspMap<Pnt6D, Pnt1D>;
template class BspMap<Pnt6D, Pnt2D>;
template class BspMap<Pnt6D, Pnt3D>;
template class BspMap<Pnt6D, Pnt4D>;
template class BspMap<Pnt6D, Pnt5D>;
template class BspMap<Pnt6D, Pnt6D>;