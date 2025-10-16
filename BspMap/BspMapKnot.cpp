#include "BspMap.h"
#include "Utils.h"
#include <assert.h>

using namespace KernelBridgeNS;
using namespace SweepNS;
using namespace std;


template <class PTypeD, class PTypeR>
const SReal BspMap<PTypeD, PTypeR>::_epsRoundKnot = 1e-9;


template <class PTypeD, class PTypeR>
SReal BspMap<PTypeD, PTypeR>::midKnotVal(SInt axis) const
{
	SInt
		nKnots = (SInt)_knotVectors[axis].size();

	return _knotVectors[axis][nKnots >> 1];
}


template <class PTypeD, class PTypeR>
void BspMap<PTypeD, PTypeR>::knotsMultiplicities(SInt dir, 
										vector<SReal> & knotVals, 
										vector<SInt> & multiplicities,
										SReal eps)	const
{
	knotVals.clear();
	multiplicities.clear();

	assert(dir < _domainDim);
	if (_knotVectors[dir].size() < 1)
		return;

	knotVals.push_back(_knotVectors[dir][0]);
	multiplicities.push_back(1);

	for (SInt i = 1, j = 0; i < _knotVectors[dir].size(); i++) {
		if (apprxEqEps(_knotVectors[dir][i], knotVals.back(), eps)) {
			multiplicities.back()++;
		}
		else {
			knotVals.push_back(_knotVectors[dir][i]);			
			multiplicities.push_back(1);
		}
	}
}


template <class PTypeD, class PTypeR>
vector<SReal> BspMap<PTypeD, PTypeR>::knotsDegreeRaise(
	const vector<SReal> & kv,
	SInt len,
	SInt order,
	SInt newOrder)  const
{
	SInt i, j, k, mult,
		dOrder = newOrder - order + 1;
	vector<SReal> newKV((len * dOrder + order) * 2);

	for (i = j = mult = 0; i < len + order - 1; i++) {
		if (apprxEqEps(kv[i], kv[i + 1], SEEps))
			mult++;
		else {
			for (k = 0; k < mult + dOrder; k++)
				newKV[j++] = kv[i];
			mult = 0;
		}
	}
	if (mult > 0 || order == 1) {
		for (k = 0; k < mult + dOrder; k++)
			newKV[j++] = kv[i];
	}

	newKV.resize(j);
	return newKV;
}


template <class PTypeD, class PTypeR>
vector<SReal> BspMap<PTypeD, PTypeR>::mergeKnotVectors(const vector<SReal> & knotVector1,
	const vector<SReal> & knotVector2,
	SInt mult)	const
{
	SInt
		i = 0,
		j = 0,
		k = 0,
		m = 0,
		len1 = (SInt)knotVector1.size(),
		len2 = (SInt)knotVector2.size(),
		len = len1 + len2;
	SReal t;
	vector<SReal>
		newKnotVector = vector<SReal>(len);

	if (mult == 0)
		mult = len + 1;

	while (i < len1 && j < len2) {
		if (knotVector1[i] < knotVector2[j]) {
			t = knotVector1[i++];
		}
		else {
			t = knotVector2[j++];
		}

		if (k == 0 || (k > 0 && !apprxEqEps(newKnotVector[k - 1], t, SEEps))) {
			newKnotVector[k++] = t;
			m = 1;
		}
		else if (m < mult) {
			newKnotVector[k++] = t;
			m++;
		}
	}

	while (i < len1)
		newKnotVector[k++] = knotVector1[i++];
	while (j < len2)
		newKnotVector[k++] = knotVector2[j++];

	newKnotVector.resize(k);

	return newKnotVector;
}



template <class PTypeD, class PTypeR>
SBool BspMap<PTypeD, PTypeR>::hasBezierKnotVector(const vector<SReal> & knotVector,
	SInt order,
	SInt length)	const
{
	SInt
		l = length + order - 1;

	if (order != length)
		return SFalse;

	switch (order) {

	case 1:
		return STrue;

	case 2:
		return apprxEqEps(knotVector[0], knotVector[1], SEEps) &&
			apprxEqEps(knotVector[l], knotVector[l - 1], SEEps);

	case 3:
		return apprxEqEps(knotVector[0], knotVector[2], SEEps) &&
			apprxEqEps(knotVector[0], knotVector[1], SEEps) &&
			apprxEqEps(knotVector[l], knotVector[l - 2], SEEps) &&
			apprxEqEps(knotVector[l], knotVector[l - 1], SEEps);

	case 4:
		return apprxEqEps(knotVector[0], knotVector[3], SEEps) &&
			apprxEqEps(knotVector[0], knotVector[2], SEEps) &&
			apprxEqEps(knotVector[0], knotVector[1], SEEps) &&
			apprxEqEps(knotVector[l], knotVector[l - 3], SEEps) &&
			apprxEqEps(knotVector[l], knotVector[l - 2], SEEps) &&
			apprxEqEps(knotVector[l], knotVector[l - 1], SEEps);

	case 5:
		return apprxEqEps(knotVector[0], knotVector[4], SEEps) &&
			apprxEqEps(knotVector[0], knotVector[3], SEEps) &&
			apprxEqEps(knotVector[0], knotVector[2], SEEps) &&
			apprxEqEps(knotVector[0], knotVector[1], SEEps) &&
			apprxEqEps(knotVector[l], knotVector[l - 4], SEEps) &&
			apprxEqEps(knotVector[l], knotVector[l - 3], SEEps) &&
			apprxEqEps(knotVector[l], knotVector[l - 2], SEEps) &&
			apprxEqEps(knotVector[l], knotVector[l - 1], SEEps);

	default:
		for (SInt i = 0; ++i < order; )
			if (!apprxEqEps(knotVector[0], knotVector[i], SEEps))
				return SFalse;
		for (SInt i = l; --i >= length; )
			if (!apprxEqEps(knotVector[l], knotVector[i], SEEps))
				return SFalse;
	}

	return STrue;
}



template <class PTypeD, class PTypeR>
SInt BspMap<PTypeD, PTypeR>::knotsLastIndexL(const vector<SReal> & knotVector,
	SReal t)	const
{
	SInt i,
		kvLength = (SInt)knotVector.size(),
		step = kvLength >> 1,
		startIndex = 0;
	SReal tmpR;

	while (step > 2) {
		tmpR = knotVector[startIndex + step];
		if (tmpR < t && !apprxEqEps(tmpR, t, SEEps))
			startIndex += step;
		step >>= 1;
	}

	for (i = startIndex; i < kvLength; i++) {
		tmpR = knotVector[i];
		if (tmpR >= t || apprxEqEps(tmpR, t, SEEps))
			break;
	}

	return i - 1;
}


template <class PTypeD, class PTypeR>
SInt BspMap<PTypeD, PTypeR>::knotsLastIndexLE(const vector<SReal> & knotVector,
	SReal t)	const
{
	SInt i,
		kvLength = (SInt)knotVector.size(),
		step = kvLength >> 1,
		startIndex = 0;
	SReal tmpR;

	while (step > 2) {
		tmpR = knotVector[startIndex + step];
		if (tmpR <= t || apprxEqEps(tmpR, t, SEEps))
			startIndex += step;
		step >>= 1;
	}

	for (i = startIndex; i < kvLength; i++) {
		tmpR = knotVector[i];
		if (tmpR > t && !apprxEqEps(tmpR, t, SEEps))
			break;
	}

	return i - 1;
}



template <class PTypeD, class PTypeR>
SInt BspMap<PTypeD, PTypeR>::knotsFirstIndexG(const vector<SReal> & knotVector,
	SReal t)	const
{
	SInt i,
		kvLength = (SInt)knotVector.size(),
		step = kvLength >> 1,
		startIndex = kvLength - 1;
	SReal tmpR;

	while (step > 2) {
		tmpR = knotVector[startIndex - step];
		if (tmpR > t && !apprxEqEps(tmpR, t, SEEps))
			startIndex -= step;
		step >>= 1;
	}

	for (i = startIndex; i >= 0; i--) {
		tmpR = knotVector[i];
		if (tmpR <= t || apprxEqEps(tmpR, t, SEEps))
			break;
	}

	return i + 1;
}


template <class PTypeD, class PTypeR>
void BspMap<PTypeD, PTypeR>::uniformOpenKnotVector(SInt dir)
{
	assert(dir < _domainDim);

	SInt i, j,
		len = _lengths[dir],
		order = _orders[dir];
	SReal
		interiorKnots = 1 + len - order;

	_knotVectors[dir] = vector<SReal>(order + len);

	vector<SReal>::iterator
		kvptr = _knotVectors[dir].begin();

	for (i = 0; i < order; i++)
		*kvptr++ = 0.0;
	for (i = 1, j = len - order; i <= j; )
		*kvptr++ = i++ / interiorKnots;
	for (j = 0; j < order; j++)
		*kvptr++ = i / interiorKnots;
}


template <class PTypeD, class PTypeR>
SBool BspMap<PTypeD, PTypeR>::knotC0Discont(const vector<SReal> & kv,
	SInt order,
	SInt length,
	SReal & t)	const
{
	SInt i, count;
	SReal
		lastT = kv[0] - 1.0;

	for (i = max(order, (length >> 1)), count = 0; i < length; i++) {
		if (apprxEqEps(lastT, kv[i], SEEps))
			count++;
		else {
			count = 1;
			lastT = kv[i];
		}

		if (count >= order) {
			t = lastT;
			return STrue;
		}
	}

	for (i = min(max(order, (length >> 1)) + order - 2, length - 1), count = 0; i >= order; i--) {
		if (apprxEqEps(lastT, kv[i], SEEps))
			count++;
		else {
			count = 1;
			lastT = kv[i];
		}

		if (count >= order) {
			t = lastT;
			return STrue;
		}
	}
	return SFalse;
}


template <class PTypeD, class PTypeR>
SBool BspMap<PTypeD, PTypeR>::knotC1Discont(const vector<SReal> & kv,
	SInt order,
	SInt length,
	SReal & t)	const
{
	SInt i, count;
	SReal
		lastT = kv[0] - 1.0;

	for (i = max(order, (length >> 1)), count = 0; i < length; i++) {
		if (apprxEqEps(lastT, kv[i], SEEps))
			count++;
		else {
			count = 1;
			lastT = kv[i];
		}

		if (count >= order - 1) {
			t = lastT;
			return STrue;
		}
	}

	for (i = min(max(order, (length >> 1)) + order - 2, length - 1), count = 0; i >= order; i--) {
		if (apprxEqEps(lastT, kv[i], SEEps))
			count++;
		else {
			count = 1;
			lastT = kv[i];
		}

		if (count >= order - 1) {
			t = lastT;
			return STrue;
		}
	}
	return SFalse;
}


template <class PTypeD, class PTypeR>
SBool BspMap<PTypeD, PTypeR>::knotC1Discont(SInt dir,
	SReal & t)	const
{
	return knotC1Discont(_knotVectors[dir], _orders[dir], _lengths[dir], t);
}


template <class PTypeD, class PTypeR>
SInt BspMap<PTypeD, PTypeR>::interiorKnots(SReal & knot)	const
{
	SInt
		dirWithKnot = -1,
		numIntrKnots = -1;

	for (SInt j = 0; j < _domainDim; j++) {
		if (_lengths[j] != _orders[j] && _lengths[j] - _orders[j] > numIntrKnots) {
			knot = _knotVectors[j][(_lengths[j] + _orders[j]) / 2];
			numIntrKnots = _lengths[j] - _orders[j];
			dirWithKnot = j;
		}
	}
	return dirWithKnot;
}


template <class PTypeD, class PTypeR>
SBool BspMap<PTypeD, PTypeR>::knotsMakeRobust(vector<SReal> &kv)
{
	SBool
		retVal = SFalse;
	vector<SReal>::iterator
		kvIter = kv.begin(),
		kvLast = kv.end();

	while (++kvIter < kvLast) {
		if (*kvIter < kvIter[-1]) {
			SReal
				d = kvIter[-1] - *kvIter;

			*kvIter = kvIter[-1];
			retVal = d > SEEps;
		}
	}
	return retVal;
}

template <class PTypeD, class PTypeR>
vector<SReal> BspMap<PTypeD, PTypeR>::knotSubtract(const vector<SReal> & kv1,
	const vector<SReal> & kv2)
{
	SInt
		len1 = (SInt)kv1.size(),
		len2 = (SInt)kv2.size();

	if (len1 <= 0)
		return vector<SReal>(0);

	SReal
		eps = SKnotEps * (kv1[len1 - 1] - kv1[0] + SEps);
	vector<SReal> newKV(len1);
	vector<SReal>::iterator
		t = newKV.begin();
	SInt
		i = 0,
		j = 0,
		newLen = 0;

	while (i < len1 && j < len2) {
		if (apprxEqEps(kv1[i], kv2[j], eps)) {
			i++;
			if (j < len2)
				j++;
		}
		else if (kv1[i] > kv2[j]) {
			j++;
		}
		else {
			*t++ = kv1[i++];
			newLen++;
		}
	}
	newKV.resize(newLen);
	return newKV;
}

template <class PTypeD, class PTypeR>
void BspMap<PTypeD, PTypeR>::knotAffineTransform(vector<SReal> & kv,
	SReal translate,
	SReal scale) const
{
	SReal
		kv0 = kv[0];
	SInt
		len = (SInt)kv.size();

	for (SInt i = 0; i < len; i++)
		kv[i] = (kv[i] - kv0) * scale + kv0 + translate;
}


template <class PTypeD, class PTypeR>
void BspMap<PTypeD, PTypeR>::knotAffineTransform2(vector<SReal> & kv,
	SReal minVal,
	SReal maxVal) const
{
	SInt
		len = (SInt)kv.size();
	SReal
		kv0 = kv[0],
		kvN = kv[len - 1],
		scale = (maxVal - minVal) / (kvN - kv0);

	for (SInt i = 0; i < len; i++)
		kv[i] = (kv[i] - kv0) * scale + minVal;
}

template <class PTypeD, class PTypeR>
void BspMap<PTypeD, PTypeR>::knotAffineTransformOrder2(vector<SReal> & kv,
	SInt order,
	SReal minVal,
	SReal maxVal)	const
{
	SInt
		len = (SInt)kv.size();
	SReal
		kv0 = kv[order - 1],
		kvN = kv[len - order],
		scale = (maxVal - minVal) / (kvN - kv0);	

	for (SInt i = 0; i < len; i++)
		kv[i] = (kv[i] - kv0) * scale + minVal;
}



template <class PTypeD, class PTypeR>
SBool BspMap<PTypeD, PTypeR>::discontKnot(const vector<SReal> &kv,
	SReal t,
	SInt order) const
{
	return (knotsFirstIndexG(kv, t) - knotsLastIndexL(kv, t)) > order;
}


template <class PTypeD, class PTypeR>
SBool BspMap<PTypeD, PTypeR>::knotsHasOpenEnd(const vector<SReal> &kv,
	SInt len,
	SInt order)	const
{
	SInt
		l = len + order - 1;

	switch (order) {
	case 1:
		return STrue;
	case 2:
		return apprxEqEps(kv[0], kv[1], SEEps)
			&& apprxEqEps(kv[l], kv[l - 1], SEEps);
	case 3:
		return apprxEqEps(kv[0], kv[1], SEEps) && apprxEqEps(kv[0], kv[2], SEEps)
			&& apprxEqEps(kv[l], kv[l - 1], SEEps) && apprxEqEps(kv[l], kv[l - 2], SEEps);
	case 4:
		return apprxEqEps(kv[0], kv[1], SEEps) && apprxEqEps(kv[0], kv[2], SEEps) && apprxEqEps(kv[0], kv[3], SEEps)
			&& apprxEqEps(kv[l], kv[l - 1], SEEps) && apprxEqEps(kv[l], kv[l - 2], SEEps) && apprxEqEps(kv[l], kv[l - 3], SEEps);

	case 5:
		return apprxEqEps(kv[0], kv[1], SEEps) && apprxEqEps(kv[0], kv[2], SEEps) && apprxEqEps(kv[0], kv[3], SEEps) && apprxEqEps(kv[0], kv[4], SEEps)
			&& apprxEqEps(kv[l], kv[l - 1], SEEps) && apprxEqEps(kv[l], kv[l - 2], SEEps) && apprxEqEps(kv[l], kv[l - 3], SEEps) && apprxEqEps(kv[l], kv[l - 4], SEEps);

	default:
		for (SInt i = 1; i < order; ++i)
			if (!apprxEqEps(kv[0], kv[i], SEEps))
				return SFalse;
		for (SInt i = l - 1; i >= len; --i)
			if (!apprxEqEps(kv[l], kv[i], SEEps))
				return SFalse;
		break;
	}
	return STrue;
}


template <class PTypeD, class PTypeR>
void BspMap<PTypeD, PTypeR>::reverseKnotVec(SInt dir)
{
	SInt
		len = (SInt)_knotVectors[dir].size();
	vector<SReal>
		newKnotVec(len);
	SReal
		t = _knotVectors[dir][len - 1] + _knotVectors[dir][0];

	for (SInt i = 0; i < len; i++)
		newKnotVec[i] = t - _knotVectors[dir][len - i - 1];

	_knotVectors[dir] = newKnotVec;
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