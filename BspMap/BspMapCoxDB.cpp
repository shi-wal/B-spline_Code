#include "BspMap.h"
#include "Utils.h"
#include <assert.h>

using namespace KernelBridgeNS;
using namespace SweepNS;
using namespace std;


template <class PTypeD, class PTypeR>
SReal BspMap<PTypeD, PTypeR>::nChooseK(SInt n,
	SInt k) const
{
	SInt j;
	SReal
		c = 1.0;

	if (k < _nCkMaxOrder)
		return _NCKTable[k][n];

	if ((k >> 1) > n) {
		for (j = k - n + 1; j <= k; j++)
			c *= j;
		for (j = 2; j <= n; j++)
			c /= j;
	}
	else {
		for (j = n + 1; j <= k; j++)
			c *= j;
		for (j = 2; j <= k - n; j++)
			c /= j;
	}

	return c;

}

template <class PTypeD, class PTypeR>
void BspMap<PTypeD, PTypeR>::bezierBasis(SInt order,
	SReal t,
	vector<SReal> & basisVec)	const
{
	SInt i;
	SReal r,
		t1 = 1.0 - t;

	if (order < _nCkMaxOrder) {
		basisVec[0] = _NCKTable[order - 1][0];

		for (i = 0, r = t; ++i < order; r *= t)
			basisVec[i] = r * _NCKTable[order - 1][i];
	}
	else {
		basisVec[0] = nChooseK(0, order - 1);

		for (i = 0, r = t; ++i < order; r *= t)
			basisVec[i] = r * nChooseK(i, order - 1);
	}

	for (i = order - 1, r = t1; --i >= 0; r *= t1)
		basisVec[i] *= r;

}




template <class PTypeD, class PTypeR>
SBool BspMap<PTypeD, PTypeR>::paramInDomain(SInt dim,
	SReal t)	const
{
	SInt
		order = _orders[dim],
		length = _lengths[dim];
	const vector<SReal> &
		knotVector = _knotVectors[dim];
	SReal t1, t2;

	t1 = knotVector[order - 1];
	t2 = knotVector[length];

	return (t1 < t || apprxEqEps(t1, t, SEEps))
		&& (t2 > t || apprxEqEps(t2, t, SEEps));
}




template <class PTypeD, class PTypeR>
vector<SReal> BspMap<PTypeD, PTypeR>::coxDeBoorBasis(SInt dim,
	SReal t,
	SInt &indexFirst)	const
{
	SInt
		i, l, index, kVLen,
		length = _lengths[dim],
		order = _orders[dim],
		origLen = length;
	vector<SReal> basisVec(order, 0);
	const vector<SReal> &
		knotVector = _knotVectors[dim];

	assert(paramInDomain(dim, t));

	if (hasBezierKnotVector(knotVector, order, length)) {
		indexFirst = 0;
		bezierBasis(order, (t - knotVector[order - 1]) / (knotVector[order] - knotVector[order - 1]), basisVec);
		return basisVec;
	}

	kVLen = order + length;
		
	index = knotsLastIndexLE(knotVector, t);

	if (index >= kVLen - 1) {
		basisVec[order - 1] = 1.0;
		indexFirst = length - order;
		return basisVec;
	}
	else
		basisVec[0] = 1.0;

	for (i = 2; i <= order; i++) {
		SInt
			indxKV = index + i - 1,
			indxKV1 = indxKV + 1,
			indxKV2 = index,
			indxBasis = i - 1;
		SReal s1, s2, s2Inv;

		if ((s2 = knotVector[indxKV1] - knotVector[indxKV2 + 1]) >= _Eps)
			s2Inv = 1.0 / s2;
		else
			s2Inv = 0.0;

		for (l = i - 1; l >= 0; l--) {
			if (s2Inv == 0.0) {
				basisVec[indxBasis--] = 0.0;
				indxKV1--;
			}
			else
				basisVec[indxBasis--] *= (knotVector[indxKV1--] - t) * s2Inv;

			if (l > 0 && (s1 = knotVector[indxKV--] - knotVector[indxKV2--]) >= _Eps) {
				s2Inv = 1.0 / s1;
				basisVec[indxBasis + 1] += basisVec[indxBasis] * (t - knotVector[indxKV2 + 1]) * s2Inv;
			}
			else
				s2Inv = 0.0;

		}

	}

	if ((indexFirst = index - order + 1) >= origLen)
		indexFirst -= origLen;

	return basisVec;
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