#include "BspMap.h"
#include <assert.h>
#include "Utils.h"

using namespace KernelBridgeNS;
using namespace SweepNS;


template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR>::BlsmAlphaCoeff::BlsmAlphaCoeff(SInt order,
	SInt length,
	SInt newOrder,
	SInt newLength) : _order(order), _length(length), _newOrder(newOrder), _newLength(newLength),
	_mat(newLength + 1, vector<SReal>(length + 1)), _colInd(newLength), _colLen(newLength),
	_KV(length + order + 1), _newKV(newLength + newOrder + 1)
{
	assert(order > 0 && length >= order);
}


template <class PTypeD, class PTypeR>
void BspMap<PTypeD, PTypeR>::BlsmAlphaCoeff::addRow(const vector<SReal> & coefs,
	SInt aRow,
	SInt colIndex,
	SInt colLength)
{
	vector<SReal>::iterator
		colPtr = _mat[aRow].begin() + colIndex;
	vector<SReal>::const_iterator
		coefsPtr = coefs.begin() + colIndex;

	assert(aRow >= 0 && aRow < _newLength && colIndex >= 0 && colLength >= 1 && colIndex + colLength <= _length);

	while (colLength--)
		*colPtr++ += *coefsPtr++;
}

template <class PTypeD, class PTypeR>
void BspMap<PTypeD, PTypeR>::BlsmAlphaCoeff::scale(SReal scl)
{
	for (SInt j = 0; j < _newLength; j++)
		for (SInt i = 0; i < _length; i++)
			_mat[j][i] *= scl;

}

template <class PTypeD, class PTypeR>
void BspMap<PTypeD, PTypeR>::BlsmAlphaCoeff::setDomain()
{
	SInt j, i;

	for (j = 0; j < _newLength; j++) {
		for (i = 0; i < _length; i++) {
			if (_mat[j][i] != 0.0)
				break;
		}
		_colInd[j] = i;

		for (i = _length - 1; i >= 0; i--) {
			if (_mat[j][i] != 0.0)
				break;
		}
		_colLen[j] = i - _colInd[j] + 1;
	}
}


template <class PTypeD, class PTypeR>
typename BspMap<PTypeD, PTypeR>::BlsmAlphaCoeff BspMap<PTypeD, PTypeR>::blsmDegRaiseNMat(const vector<SReal> & KV,
	SInt order,
	SInt newOrder,
	SInt len)	const
{
	BlsmAlphaCoeff A, A2;
	BlsmAlphaCoeff & A1 = A;
	SInt
		i = 0;

	while (order < newOrder) {
		if (i == 0) {
			A = blsmDegRaiseMat(KV, order, len);
			A1 = A;
		}
		else {
			A2 = blsmDegRaiseMat(A1._newKV, A1._newOrder, A1._newLength);
			A = blsmDegRaiseMatProd(A1, A2);
			A1 = A;
		}
		order++;
		i++;
	}
	return A;
}



template <class PTypeD, class PTypeR>
typename BspMap<PTypeD, PTypeR>::BlsmAlphaCoeff BspMap<PTypeD, PTypeR>::blsmDegRaiseMat(const vector<SReal> & KV,
	SInt order,
	SInt len)	const
{
	SReal t;

	assert(order > 0 && order <= len);

	if (knotC0Discont(KV, order, len, t))
		assert(0);

	vector<SReal>
		newKV = knotsDegreeRaise(KV, len, order, order + 1);

	SInt
		newLen = (SInt)newKV.size();

	newLen -= order + 1;

	BlsmAlphaCoeff bAlphaC(order, len, order + 1, newLen);

	if (order > 1) {
		BlsmEvalCache blsmCache;
		vector<SReal> blossomValues(order);

		for (SInt l = 0; l < newLen; l++) {
			std::copy(newKV.begin() + (l + 2), newKV.begin() + (l + 2 + order), blossomValues.begin());

			for (SInt j = 0; j < order; j++) {
				SInt index, length;

				blsmEvalSymb(order, KV, blossomValues, order - 1, index, length, blsmCache);
				bAlphaC.addRow(blsmCache._coefs, l, index, length);
				blossomValues[j] = newKV[l + 1 + j];
			}
		}
		bAlphaC.scale(1.0 / order);
	}
	else {
		for (SInt l = 0; l < newLen; l++)
			bAlphaC._mat[l][knotsLastIndexLE(KV, newKV[l])] = 1.0;
	}
	bAlphaC.setDomain();
	bAlphaC._KV = KV;
	bAlphaC._newKV = newKV;

	return bAlphaC;
}


/* returns answer in blsmCache*/
template <class PTypeD, class PTypeR>
void BspMap<PTypeD, PTypeR>::blsmEvalSymb(SInt order,
	const vector<SReal> & knots,
	const vector<SReal> & blsmVals,
	SInt blsmLen,
	SInt & retIdxFirst,
	SInt & retLength,
	BlsmEvalCache & blsmCache)  const
{
	SInt i, k, m, symbVecLen,
		knotsLen = (SInt)knots.size(),
		ptLen = knotsLen - order,
		J = 0,
		lastJ = -1;
	SReal
		tMin = knots[0],
		tMax = knots[knotsLen - 1];

	symbVecLen = knotsLen + 1;	    //TODO: revisit (+1)

	if (knotsLen > (SInt)blsmCache._knots.size())
		blsmCache._symbVec = vector<BlsmSymb>(symbVecLen, BlsmSymb(symbVecLen));

	vector<BlsmSymb> &
		symbVec = blsmCache._symbVec;

	if (blsmCache._order == order && blsmCache._blsmVals == blsmVals && blsmCache._knots == knots) {
		retIdxFirst = blsmCache._idxFirst;
		retLength = blsmCache._idxLength;
		return;
	}
	blsmCache._order = order;
	blsmCache._blsmVals = blsmVals;
	blsmCache._knots = knots;

	assert(order > blsmLen);

	for (i = 0; i < knotsLen; i++) {
		SInt
			minIdx = max(0, i - order * 2),
			maxIdx = min(order * 4, knotsLen - minIdx);

		std::fill_n(symbVec[i]._coefs.begin() + minIdx, maxIdx, 0.0);
		symbVec[i]._coefs[i] = 1.0;
		symbVec[i]._min = symbVec[i]._max = i;
	}

	for (k = 1; k < blsmLen; k++)
		if (blsmVals[k - 1] > blsmVals[k])
			break;

	SBool
		blossomValsNotInOrder = k < blsmLen;

	for (k = 0; k < blsmLen; k++) {
		typename vector<BlsmSymb>::iterator S;
		vector<SReal>::const_iterator K;
		SReal
			blsmVal = blsmVals[k],
			clippedBlsmVal = blsmVal;
		SInt
			orderK = order - 1 - k;

		assert(blsmVal > tMin - SEEps && blsmVal < tMax + SEEps);

		if (k == 0 || blossomValsNotInOrder) {
			clippedBlsmVal = validateMinMaxDomain(clippedBlsmVal, tMin, tMax);
			J = knotsLastIndexLE(knots, clippedBlsmVal);
		}
		else {
			for (J = lastJ; J < knotsLen && (knots[J] <= clippedBlsmVal || apprxEqEps(knots[J], clippedBlsmVal, SEEps)); J++);
			J--;
		}

		if (J - lastJ > 1) {
			if (k == 0)
				lastJ = J;
			else
				J = lastJ + 1;
		}
		lastJ = J;

		for (i = J, S = symbVec.begin() + J, K = knots.begin() + i;
			i > J - orderK && i > 0;
			i--, S--, K--) {
			SInt
				minIdx = min(S[0]._min, S[-1]._min),
				maxIdx = max(S[0]._max, S[-1]._max);
			vector<SReal>::iterator
				s0Coefs = S[0]._coefs.begin() + minIdx,
				s1Coefs = S[-1]._coefs.begin() + minIdx;
			SReal
				u = orderK + i < knotsLen && K[orderK] - *K > 0 ? (blsmVal - *K) / (K[orderK] - *K) : 0.0,
				u1 = 1.0 - u;

			for (m = maxIdx - minIdx; m-- >= 0; s0Coefs++)
				*s0Coefs = u * *s0Coefs + u1 * *s1Coefs++;

			S[0]._min = minIdx;
			S[0]._max = maxIdx;
		}
	}

	for (i = symbVec[J]._min; symbVec[J]._coefs[i] == 0 && i < knotsLen; i++);
	retIdxFirst = i;

	for (i = symbVec[J]._max; symbVec[J]._coefs[i] == 0 && i > 0; i--);
	retLength = i - retIdxFirst + 1;

	blsmCache._coefs = symbVec[J]._coefs;
	blsmCache._idxFirst = retIdxFirst;
	blsmCache._idxLength = retLength;
}


template <class PTypeD, class PTypeR>
SReal BspMap<PTypeD, PTypeR>::validateMinMaxDomain(SReal t, SReal tMin, SReal tMax)   const
{
	if (t >= tMax - SDmnEps) {
		t = tMax - SDmnEps;
		if (t >= tMax) {
			t = tMax - SDmnEps * tMax;
			assert(t < tMax);
		}
	}
	if (t < tMin)
		t = tMin;
	return t;
}


template <class PTypeD, class PTypeR>
typename BspMap<PTypeD, PTypeR>::BlsmAlphaCoeff BspMap<PTypeD, PTypeR>::blsmDegRaiseMatProd(
	const BlsmAlphaCoeff &a1,
	const BlsmAlphaCoeff &a2)   const
{
	SInt
		a2NewLength = a2._newLength;

	assert(a1._newOrder == a2._order && a1._newLength == a2._length);

	BlsmAlphaCoeff a = a2;

	for (SInt i = 0; i < a1._length; i++) {
		for (SInt j = 0; j < a2NewLength; j++) {
			SInt
				colIndex = a2._colInd[j],
				colLength = colIndex + a2._colLen[j];
			SReal
				R = 0.0;
			vector<SReal>::const_iterator
				a2r = a2._mat[j].begin() + colIndex;

			for (SInt k = colIndex; k < colLength; k++)
				R += a1._mat[k][i] * *a2r++;

			a._mat[j][i] = R;
		}
	}

	a._order = a1._order;
	a._length = a1._length;
	a._KV = a1._KV;
	a.setDomain();
	return a;
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