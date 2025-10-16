#include "BspMap.h"
#include <assert.h>
#include "Utils.h"

using namespace KernelBridgeNS;
using namespace SweepNS;


template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR> BspMap<PTypeD, PTypeR>::degreeRaise(vector<SInt> & newOrders)    const
{
	SReal t;
	BspMap<PTypeD, PTypeR> retMV;

	for (SInt i = 0; i < _domainDim; i++) {
		if (knotC0Discont(_knotVectors[i], _orders[i], _lengths[i], t)) {
			BspMap<PTypeD, PTypeR> lmv, rmv;

			subdivide(t, i, lmv, rmv);

			BspMap<PTypeD, PTypeR>
				retMV1 = lmv.degreeRaise(newOrders),
				retMV2 = rmv.degreeRaise(newOrders);

			retMV = retMV1.merge(retMV2, i, STrue);
			return retMV;
		}
	}

	/* If the B-spline MV is actually a Bezier, treat as such. */
	if (interiorKnots(t) == -1) {
		retMV = degreeRaise2(newOrders);
		return retMV;
	}

	retMV = *this;

	for (SInt k = 0; k < _domainDim; k++) {
		if (retMV._orders[k] >= newOrders[k])
			continue;

		BlsmAlphaCoeff
			A = blsmDegRaiseNMat(retMV._knotVectors[k], retMV._orders[k], newOrders[k], retMV._lengths[k]);
		SInt
			newLen = A._newLength;

		std::swap(retMV._lengths[k], newLen);
		std::swap(retMV._orders[k], newOrders[k]);

		BspMap<PTypeD, PTypeR> rmv( retMV._orders, retMV._lengths);

		std::swap(retMV._lengths[k], newLen);
		std::swap(retMV._orders[k], newOrders[k]);

		for (SInt i = 0; i < _domainDim; i++) {
			if (i == k)
				rmv._knotVectors[i] = A._newKV;
			else
				rmv._knotVectors[i] = retMV._knotVectors[i];
		}

		SInt
			rStepSize = rmv._subSpaces[k],
			stepSize = retMV._subSpaces[k],
			index = 0,
			rIndex = 0;

		vector<SInt> indices(_domainDim, 0), rIndices(_domainDim, 0);

		do {
			SInt
				rmvInd = rIndex;

			for (SInt l = 0; l < newLen; l++, rmvInd += rStepSize) {
				vector<SReal>::iterator
					blendVals = A._mat[l].begin() + A._colInd[l];
				SInt
					pInd = index + A._colInd[l] * stepSize;

				for (SInt j = 0; j < A._colLen[l]; j++, pInd += stepSize)
					rmv._points[rmvInd] += retMV._points[pInd] * *blendVals++;

			}

			rmv.incrementPointsSkipIndex(rIndices, k, rIndex);
		} while (retMV.incrementPointsSkipIndex(indices, k, index));

		retMV = rmv;
	}

	return retMV;
}


template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR> BspMap<PTypeD, PTypeR>::degreeRaise2(vector<SInt> & newOrders)    const
{
	vector<SInt> orders(_domainDim);

	for (SInt i = 0; i < _domainDim; i++) {
		assert(_orders[i] <= newOrders[i]);
		orders[i] = newOrders[i] - _orders[i] + 1;
	}

	BspMap<PTypeD, PTypeR> unitMV( orders, orders);

	for (SInt i = 0; i < _domainDim; i++) {
		SReal min, max;
		vector<SReal>::iterator
			kv = unitMV._knotVectors[i].begin();

		domain(min, max, i);
		for (SInt j = 0; j < orders[i]; j++)
			*kv++ = min;
		for (SInt j = 0; j < orders[i]; j++)
			*kv++ = max;
	}

	for (SInt i = 0; i < unitMV._subSpaces[_domainDim]; i++)
		for (SInt j = 0; j < _rangeDim; j++)
			unitMV._points[i].setCoord(j, 1.0);

	BspMap<PTypeD, PTypeR>
		raisedMV = *this * unitMV;

	return raisedMV;
}


template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR> BspMap<PTypeD, PTypeR>::bezierMultiply(const BspMap<PTypeD, PTypeR> &mv)   const
{
	for (SInt i = 0; i < _domainDim; i++) {
		assert(_orders[i] == _lengths[i]);
		assert(mv._orders[i] == mv._lengths[i]);
	}

	assert(_rangeDim == mv._rangeDim);

	vector<SInt> lengths(_domainDim);

	for (SInt i = 0; i < _domainDim; i++)
		lengths[i] = _lengths[i] + mv._lengths[i] - 1;

	BspMap<PTypeD, PTypeR> prodMV( lengths, lengths);
	SBool
		useIChooseKTable = STrue;

	for (SInt i = 0; i < _domainDim; i++) {
		if (prodMV._orders[i] >= _nCkMaxOrder) {
			useIChooseKTable = SFalse;
			break;
		}
	}

	SBool
		useMultiIChooseKTable = STrue;

	for (SInt i = 0; i < _domainDim; i++) {
		if (_orders[i] >= _nCkMaxOrder2 || mv._orders[i] >= _nCkMaxOrder2) {
			useMultiIChooseKTable = SFalse;
			break;
		}
	}

	vector<SInt> indices1(_domainDim, 0), indices2(_domainDim, 0);
	SInt
		index1 = 0,
		index2 = 0,
		prodIndex = 0;

	if (useIChooseKTable) {
		do {
			SReal
				coef0 = 1.0;

			do {
				SReal coef;

				if (indices2[0] == 0) {
					vector<SInt>::iterator
						i1 = indices1.begin() + 1,
						i2 = indices2.begin() + 1,
						prodSubSpcs = prodMV._subSpaces.begin();
					vector<SInt>::const_iterator
						o1 = _orders.begin(),
						o2 = mv._orders.begin();

					coef0 = 1.0;

					if (useMultiIChooseKTable) {
						for (SInt i = 1; i < _domainDim; i++, i1++, i2++) {
							coef0 *= _IcKJcMIJcKMTable[*++o1 - 1][*++o2 - 1][*i1][*i2];
						}
					}
					else {
						vector<SInt>::iterator
							oP = prodMV._orders.begin();

						for (SInt i = 1; i < _domainDim; i++, i1++, i2++) {
							coef0 *= _NCKTable[*++o1 - 1][*i1] * _NCKTable[*++o2 - 1][*i2] / _NCKTable[*++oP - 1][*i1 + *i2];
						}
					}

					SInt i;
					for (i = prodIndex = 0; i < _domainDim; i++)
						prodIndex += *prodSubSpcs++ * (indices1[i] + indices2[i]);
				}
				else {
					prodIndex++;
				}

				if (useMultiIChooseKTable) {
					coef = coef0 * _IcKJcMIJcKMTable[_orders[0] - 1][mv._orders[0] - 1][indices1[0]][indices2[0]];
				}
				else {
					coef = coef0 * _NCKTable[_orders[0] - 1][indices1[0]] * _NCKTable[mv._orders[0] - 1][indices2[0]] /
						_NCKTable[prodMV._orders[0] - 1][indices1[0] + indices2[0]];
				}

				prodMV._points[prodIndex] += _points[index1].mult(mv._points[index2]) * coef;
			} while (mv.incremetPointsIndex(indices2, index2));
		} while (incremetPointsIndex(indices1, index1));
	}
	else {
		do {
			SReal
				coef0 = 1.0;

			do {
				SReal coef;

				if (indices2[0] == 0) {
					vector<SInt>::iterator
						i1 = indices1.begin() + 1,
						i2 = indices2.begin() + 1,
						oP = prodMV._orders.begin(),
						prodSubSpcs = prodMV._subSpaces.begin();
					vector<SInt>::const_iterator
						o1 = _orders.begin(),
						o2 = mv._orders.begin();

					coef0 = 1.0;

					SInt i;
					for (i = 1; i < _domainDim; i++, i1++, i2++)
						coef0 *= nChooseK(*i1, *++o1 - 1) * nChooseK(*i2, *++o2 - 1) / nChooseK(*i1 + *i2, *++oP - 1);

					for (i = prodIndex = 0; i < _domainDim; i++)
						prodIndex += *prodSubSpcs++ * (indices1[i] + indices2[i]);

				}
				else {
					prodIndex++;
				}

				coef = coef0 * nChooseK(indices1[0], _orders[0] - 1) * nChooseK(indices2[0], mv._orders[0] - 1) /
					nChooseK(indices1[0] + indices2[0], prodMV._orders[0] - 1);

				prodMV._points[prodIndex] += _points[index1].mult(mv._points[index2]) * coef;

			} while (mv.incremetPointsIndex(indices2, index2));
		} while (incremetPointsIndex(indices1, index1));
	}

	return prodMV;
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