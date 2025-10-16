#include "BspMap.h"
#include <assert.h>
#include "Utils.h"

using namespace KernelBridgeNS;
using namespace SweepNS;

template <class PTypeD, class PTypeR>
SBool BspMap<PTypeD, PTypeR>::domainCheck(const BspMap &mv)	const
{
	if (_domainDim != mv._domainDim)
		return SFalse;

	for (SInt i = 0; i < _domainDim; i++) {
		SReal min1, max1, min2, max2;

		domain(min1, max1, i);
		mv.domain(min2, max2, i);
		if (!apprxEq(min1, min2) || !apprxEq(max1, max2))
			return SFalse;
	}
	return STrue;
}


template <class PTypeD, class PTypeR>
vector<BspMap<Pnt1D, Pnt1D>> BspMap<PTypeD, PTypeR>::splitScalar()    const
{
	SInt
		nCoords = PTypeR::dim(),
		nPts = (SInt)_points.size();
	vector<BspMap<Pnt1D, Pnt1D>> retMVs(nCoords);

	for (SInt i = 0; i < nCoords; i++) {
		retMVs[i] = BspMap<Pnt1D, Pnt1D>( _orders, _lengths);
		retMVs[i]._orders = _orders;
		retMVs[i]._lengths = _lengths;

		for (SInt j = 0; j < _domainDim; j++)
			retMVs[i]._knotVectors[j] = _knotVectors[j];

		for (SInt j = 0; j < nPts; j++)
			retMVs[i]._points[j].setCoord(0, _points[j].coord(i));
	}
	return retMVs;
}


template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR>::BspMap(const vector<BspMap<Pnt1D, Pnt1D>> &scalarMVs) :
	_domainDim(scalarMVs[0]._domainDim), _rangeDim(PTypeR::dim())
{
	assert(_rangeDim == scalarMVs.size());

	vector<BspMap<Pnt1D, Pnt1D>>
		copyMVs = scalarMVs;

	for (SInt i = 0; i < copyMVs.size(); i++) {
		for (SInt j = i + 1; j < copyMVs.size(); j++) {
			copyMVs[i].makeCompatible(copyMVs[j], STrue, STrue);
		}
	}

	for (SInt i = 0; i < _domainDim; i++) {
		for (SInt j = 1; j < _rangeDim; j++) {
			assert(copyMVs[0]._orders[i] == copyMVs[j]._orders[i]);
			assert(copyMVs[0]._lengths[i] == copyMVs[j]._lengths[i]);
		}
	}

	_orders = copyMVs[0]._orders;
	_lengths = copyMVs[0]._lengths;	
	_subSpaces = copyMVs[0]._subSpaces;
	_knotVectors = copyMVs[0]._knotVectors;
	_points = vector<PTypeR>(copyMVs[0]._points.size());

	for (SInt i = 0; i < copyMVs[0]._points.size(); i++) {
		for (SInt j = 0; j < _rangeDim; j++) {
			_points[i].setCoord(j, copyMVs[j]._points[i].coord(0));
		}
	}
}


template <class PTypeD, class PTypeR>
void BspMap<PTypeD, PTypeR>::preconditionScale(SReal scale)
{
	
	for (auto ptItr = _points.begin(); ptItr != _points.end(); ++ptItr)
		(*ptItr) *= scale;
}


template <class PTypeD, class PTypeR>
SReal BspMap<PTypeD, PTypeR>::maxScalarPtMagnitude()	const
{
	SReal
		rMax = 0.0;

	for (auto ptItr = _points.begin(); ptItr != _points.end(); ++ptItr)
		rMax = std::max(rMax, fabs(ptItr->coord(0)));

	return rMax;
}



template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR> BspMap<PTypeD, PTypeR>::scale(SReal scale)	const
{
	BspMap<PTypeD, PTypeR>
		retBsp = *this;

	for (auto ptItr = retBsp._points.begin(); ptItr != retBsp._points.end(); ++ptItr)
		(*ptItr) *= scale;

	return retBsp;
}


template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR> BspMap<PTypeD, PTypeR>::operator + (const BspMap<PTypeD, PTypeR> &mv) const
{
	assert(domainCheck(mv));

	BspMap<PTypeD, PTypeR>
		mv1 = *this,
		mv2 = mv;

	mv1.makeCompatible(mv2, STrue, STrue);

	SInt
		nPts = (SInt)mv1._points.size();

	for (SInt i = 0; i < nPts; i++)
		mv1._points[i] += mv2._points[i];

	return mv1;
}


template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR> BspMap<PTypeD, PTypeR>::operator + (const PTypeR &pt) const
{
	BspMap<PTypeD, PTypeR>
		retBsp = *this;

	for (auto ptItr = retBsp._points.begin(); ptItr != retBsp._points.end(); ++ptItr)
		(*ptItr) += pt;

	return retBsp;

}


template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR> BspMap<PTypeD, PTypeR>::operator - (const BspMap<PTypeD, PTypeR> &mv) const
{
	assert(domainCheck(mv));

	BspMap<PTypeD, PTypeR>
		mv1 = *this,
		mv2 = mv;

	mv1.makeCompatible(mv2, STrue, STrue);

	SInt
		nPts = (SInt)mv1._points.size();

	for (SInt i = 0; i < nPts; i++)
		mv1._points[i] -= mv2._points[i];

	return mv1;
}


template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR> BspMap<PTypeD, PTypeR>::operator *(const BspMap<PTypeD, PTypeR> &mv)   const
{
	assert(domainCheck(mv));

	BspMap<PTypeD, PTypeR>
		retMV = bspMultiplyAux(mv);

	return retMV;
}

template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR> BspMap<PTypeD, PTypeR>::bspMultiplyAux(const BspMap<PTypeD, PTypeR> &mv) const
{
	BspMap<PTypeD, PTypeR> prodMV;
	SBool
		isBezier = STrue;

	for (SInt i = 0; i < _domainDim; i++) {
		if (_lengths[i] != _orders[i] || mv._lengths[i] != mv._orders[i]) {
			SReal
				t = _lengths[i] != _orders[i] ? _knotVectors[i][(_lengths[i] + _orders[i]) >> 1] :
				mv._knotVectors[i][(mv._lengths[i] + mv._orders[i]) >> 1];

			BspMap<PTypeD, PTypeR> subMV1LMV, subMV1RMV, subMV2LMV, subMV2RMV;

			subdivide(t, i, subMV1LMV, subMV1RMV);
			mv.subdivide(t, i, subMV2LMV, subMV2RMV);

			SBool
				discontMerge = discontKnot(_knotVectors[i], t, _orders[i]) || discontKnot(mv._knotVectors[i], t, mv._orders[i]);

			BspMap<PTypeD, PTypeR>
				mv1Prod = subMV1LMV.bspMultiplyAux(subMV2LMV),
				mv2Prod = subMV1RMV.bspMultiplyAux(subMV2RMV);

			prodMV = mv1Prod.merge(mv2Prod, i, discontMerge);

			isBezier = SFalse;
			break;
		}
	}

	if (isBezier) {
		prodMV = bezierMultiply(mv);
		for (SInt i = 0; i < _domainDim; i++)
			prodMV.uniformOpenKnotVector(i);
	}

	for (SInt i = 0; i < _domainDim; i++) {
		SReal min, max;

		domain(min, max, i);
		knotAffineTransform2(prodMV._knotVectors[i], min, max);
	}

	return prodMV;
}

template <class PTypeD, class PTypeR>
BspMap<Pnt1D, Pnt1D> BspMap<PTypeD, PTypeR>::dotProd(const BspMap<PTypeD, PTypeR> &mv)   const
{
	BspMap<PTypeD, PTypeR>
		prodMV = *this * mv;

	BspMap<Pnt1D, Pnt1D> retMV( prodMV._orders, prodMV._lengths, prodMV._knotVectors);

	SInt
		nPts = (SInt)prodMV._points.size();

	for (SInt i = 0; i < nPts; i++) {
		SReal
			sum = 0.0;

		for (SInt j = 0; j < prodMV._rangeDim; j++) {
			sum += prodMV._points[i].coord(j);
		}
		retMV._points[i].setAllCoord(sum);
	}
	return retMV;
}


template <class PTypeD, class PTypeR>
BspMap<Pnt1D, Pnt1D> BspMap<PTypeD, PTypeR>::dotProd(const PTypeR & pt)	const
{
	vector<BspMap<Pnt1D, Pnt1D>>
		scalarBsps = splitScalar();

	BspMap<Pnt1D, Pnt1D>
		retBsp = scalarBsps[0].scale(pt.coord(0));

	for (SInt i = 1; i < _rangeDim; i++) 
		retBsp = retBsp + scalarBsps[i].scale(pt.coord(i));
	
	return retBsp;
}

template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR> BspMap<PTypeD, PTypeR>::crossProd(const BspMap<PTypeD, PTypeR> &mv)	const
{
	assert(_rangeDim == 3 && mv._rangeDim == 3);

	vector<BspMap<Pnt1D, Pnt1D>>
		scalarsA = splitScalar(),
		scalarsB = mv.splitScalar(),
		scalarsAxB = {
				(scalarsA[1] * scalarsB[2]) - (scalarsA[2] * scalarsB[1]),
				(scalarsA[2] * scalarsB[0]) - (scalarsA[0] * scalarsB[2]),
				(scalarsA[0] * scalarsB[1]) - (scalarsA[1] * scalarsB[0]) };

	return BspMap<PTypeD, PTypeR>(scalarsAxB);
}


template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR> BspMap<PTypeD, PTypeR>::scalarProd(const BspMap<Pnt1D, Pnt1D> &mv)	const
{
	assert(_rangeDim == 3);

	vector<BspMap<Pnt1D, Pnt1D>> 
		scalars = splitScalar(),
		scalarsPr = {
				scalars[0] * mv,
				scalars[1] * mv,
				scalars[2] * mv	};

	return BspMap<PTypeD, PTypeR>(scalarsPr);
}


template <class PTypeD, class PTypeR>
SBool BspMap<PTypeD, PTypeR>::noZeroCrossing(SReal eps)	const
{
	SInt
		len = _subSpaces[_domainDim],
		numZeros = 0;
	SBool
		pos = SFalse,
		neg = SFalse;

	for (SInt i = 0; i < len; i++) {
		SReal
			x = _points[i].coord(0);

		numZeros += (x > -eps && x < eps);
		pos |= (x > 0.0);
		neg |= (x < 0.0);
		if (pos && neg)
			return SFalse;
	}

	if (numZeros == len)
		return STrue;

	if (numZeros > 0 && (pos || neg))
		return SFalse;

	return STrue;
}

template <class PTypeD, class PTypeR>
SBool BspMap<PTypeD, PTypeR>::noZeroCrossingSub(const BspMap<PTypeD, PTypeR> & mv,
	SReal eps)	const
{
	SInt
		len1 = _subSpaces[_domainDim],
		len2 = mv._subSpaces[_domainDim],
		numZeros = 0;
	SBool
		pos = SFalse,
		neg = SFalse;

	for (SInt i = 0; i < len1; i++) {
		for (SInt j = 0; j < len2; j++) {
			SReal
				x = _points[i].coord(0) - mv._points[j].coord(0);

			numZeros += (x > -eps && x < eps);
			pos |= (x > 0.0);
			neg |= (x < 0.0);
			if (pos && neg)
				return SFalse;
		}
	}

	if (numZeros == len1 * len2)
		return STrue;

	if (numZeros > 0 && (pos || neg))
		return SFalse;

	return STrue;
}

template <class PTypeD, class PTypeR>
SBool BspMap<PTypeD, PTypeR>::allPtNegative()	const
{
	SInt
		len = _subSpaces[_domainDim];

	for (SInt i = 0; i < len; i++) {
		if (_points[i].coord(0) >= 0.0)
			return SFalse;
	}

	return STrue;
}


template <class PTypeD, class PTypeR>
SBool BspMap<PTypeD, PTypeR>::allPtNegativeSub(const BspMap<PTypeD, PTypeR> &mv)	const
{
	SInt
		len1 = _subSpaces[_domainDim],
		len2 = mv._subSpaces[_domainDim];

	for (SInt i = 0; i < len1; i++) {
		for (SInt j = 0; j < len2; j++) {
			if (_points[i].coord(0) - mv._points[j].coord(0) >= 0.0)
				return SFalse;
		}
	}

	return STrue;
}


template <class PTypeD, class PTypeR>
SBool BspMap<PTypeD, PTypeR>::allPtPositive()	const
{
	SInt
		len = _subSpaces[_domainDim];

	for (SInt i = 0; i < len; i++) {
		if (_points[i].coord(0) <= 0.0)
			return SFalse;
	}
	return STrue;
}


template <class PTypeD, class PTypeR>
SBool BspMap<PTypeD, PTypeR>::allPtPositiveSub(const BspMap<PTypeD, PTypeR> &mv)	const
{
	SInt
		len1 = _subSpaces[_domainDim],
		len2 = mv._subSpaces[_domainDim];

	for (SInt i = 0; i < len1; i++) {
		for (SInt j = 0; j < len2; j++) {
			if (_points[i].coord(0) - mv._points[j].coord(0) <= 0.0)
				return SFalse;
		}
	}

	return STrue;
}



template <class PTypeD, class PTypeR>
vector<Pnt2D> BspMap<PTypeD, PTypeR>::deriveAllBounds()	const
{
	assert(PTypeR::dim() == 1);

	vector<Pnt2D> minMax(_domainDim);

	for (SInt i = 0; i < _domainDim; i++)
		minMax[i].setAllCoord(SInfnty, -SInfnty);

	SInt
		maxLength = 0;

	for (SInt i = 0; i < _domainDim; i++)
		maxLength = max(maxLength, _lengths[i]);

	vector<vector<SReal>> denoms(_domainDim, vector<SReal>(maxLength));

	for (SInt dir = 0; dir < _domainDim; dir++) {
		SInt
			order = _orders[dir],
			length = _lengths[dir];
		vector<SReal>
			& denomDir = denoms[dir];

		for (SInt i = 0; i < length; i++) {
			denomDir[i] = _knotVectors[dir][i + order] - _knotVectors[dir][i + 1];
			denomDir[i] = denomDir[i] < SEEps ? 1.0 / SEEps : 1.0 / denomDir[i];
		}
	}

	vector<SInt> indexVec(_domainDim, 0);
	SInt
		index = 0;

	do {
		SReal
			ptIndex = _points[index].coord(0);

		for (SInt dir = 0; dir < _domainDim; dir++) {
			SReal der;
			Pnt2D
				& mm = minMax[dir];
			SInt order;

			if ((order = _orders[dir]) > 1) {
				if (indexVec[dir] >= _lengths[dir] - 1) {
					continue;
				}
				der = (_points[index + _subSpaces[dir]].coord(0) - ptIndex) * denoms[dir][indexVec[dir]];
			}
			else
				der = 0.0;

			if (mm.coord(0) > der)
				mm.setCoord(0, der);
			if (mm.coord(1) < der)
				mm.setCoord(1, der);
		}

	} while (incremetPointsIndex(indexVec, index));

	for (SInt i = 0; i < _domainDim; i++)
		minMax[i] *= (_orders[i] - 1);

	return minMax;
}


template <class PTypeD, class PTypeR>
PTypeD BspMap<PTypeD, PTypeR>::evalGradientNumeric(const PTypeD & param)    const
{
	assert(PTypeR::dim() == 1);

	SReal
		pos = evaluate(param).coord(0);
	PTypeD grad,
		prmEps = param;

	for (SInt i = 0; i < _domainDim; i++) {
		SReal tMin, tMax, eps,
			oldParam = prmEps.coord(i);

		domain(tMin, tMax, i);

		SReal
			dtEps = max((tMax - tMin) * _numGradEps, SEEps),
			iCoord = prmEps.coord(i);

		if (iCoord + dtEps > tMax) {
			eps = -dtEps;
			prmEps.setCoord(i, iCoord - dtEps);
		}
		else {
			eps = dtEps;
			prmEps.setCoord(i, iCoord + dtEps);
		}

		SReal
			pos2 = evaluate(prmEps).coord(0);

		grad.setCoord(i, (pos2 - pos) / eps);

		prmEps.setCoord(i, oldParam);
	}
	return grad;
}


template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR> BspMap<PTypeD, PTypeR>::reverseDir(SInt axis)	const
{
	assert(axis >= 0 && axis < _domainDim);
	SInt
		lengthAxis = _lengths[axis];
	BspMap<PTypeD, PTypeR>
		newBsp = *this;
	vector<SInt> indices(_domainDim, 0);
	SInt
		index = 0;

	newBsp.reverseKnotVec(axis);

	do {
		for (SInt k = 0; k < lengthAxis / 2; k++) {
			SInt
				indx1 = index + k * _subSpaces[axis],
				indx2 = index + (lengthAxis - k - 1) * _subSpaces[axis];

			std::swap(newBsp._points[indx1], newBsp._points[indx2]);
		}
	} while (newBsp.incrementPointsSkipIndex(indices, axis, index));

	return newBsp;
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

template Pnt1D BspMap<Pnt1D,Pnt1D>::evalGradientNumeric(const Pnt1D & param) const;
template Pnt2D BspMap<Pnt2D,Pnt1D>::evalGradientNumeric(const Pnt2D & param) const;
template Pnt3D BspMap<Pnt3D,Pnt1D>::evalGradientNumeric(const Pnt3D & param) const;
template Pnt4D BspMap<Pnt4D,Pnt1D>::evalGradientNumeric(const Pnt4D & param) const;
template Pnt5D BspMap<Pnt5D,Pnt1D>::evalGradientNumeric(const Pnt5D & param) const;
template Pnt6D BspMap<Pnt6D,Pnt1D>::evalGradientNumeric(const Pnt6D & param) const;