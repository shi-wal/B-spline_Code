#include "BspMap.h"
#include <assert.h>
#include "Utils.h"
#include <iostream>
#include <omp.h>

using namespace KernelBridgeNS;
using namespace SweepNS;



template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR>::BspMap() :
	_domainDim(0), _rangeDim(PTypeR::dim())
{
}


template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR>::BspMap(
	const vector<SInt> &orders,
	const vector<SInt> &lengths) :
	_domainDim(PTypeD::dim()), _rangeDim(PTypeR::dim()), _orders(orders), _lengths(lengths), _subSpaces(_domainDim + 1), _knotVectors(_domainDim), _id1(-1), _id2(-1)
{
	_subSpaces[0] = 1;
	for (SInt i = 0; i < _domainDim; i++) {
		assert(_orders[i] <= _lengths[i]);
		_subSpaces[i + 1] = _subSpaces[i] * _lengths[i];
	}

	for (SInt i = 0; i < _domainDim; i++)
		_knotVectors[i] = vector<SReal>(_orders[i] + _lengths[i]);

	_points = vector<PTypeR>(_subSpaces[_domainDim]);
}


template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR>::BspMap(
	const vector<SInt> &orders,
	const vector<SInt> &lengths,
	const vector<vector<SReal>> &knotVectors) :
	_domainDim(PTypeD::dim()), _rangeDim(PTypeR::dim()), _orders(orders), _lengths(lengths), _subSpaces(_domainDim + 1), _knotVectors(knotVectors), _id1(-1), _id2(-1)
{
	_subSpaces[0] = 1;
	for (SInt i = 0; i < _domainDim; i++) {
		assert(_orders[i] <= _lengths[i]);
		_subSpaces[i + 1] = _subSpaces[i] * _lengths[i];
	}

	_points = vector<PTypeR>(_subSpaces[_domainDim]);
}


template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR>::BspMap(
	const vector<SInt> &orders,
	const vector<SInt> &lengths,
	const vector<PTypeR> &points,
	SInt id1,
	SInt id2) :
	_domainDim(PTypeD::dim()), _rangeDim(PTypeR::dim()), _orders(orders), _lengths(lengths), _subSpaces(_domainDim + 1), _knotVectors(_domainDim), _points(points), _id1(id1), _id2(id2)
{
	_subSpaces[0] = 1;
	for (SInt i = 0; i < _domainDim; i++) {
		assert(_orders[i] <= _lengths[i]);
		_subSpaces[i + 1] = _subSpaces[i] * _lengths[i];
		uniformOpenKnotVector(i);
	}
}


template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR>::BspMap(
	const vector<SInt> &orders,
	const vector<SInt> &lengths,
	const vector<vector<SReal>> &knotVectors,
	const vector<PTypeR> &points) :
	_domainDim(PTypeD::dim()), _rangeDim(PTypeR::dim()), _orders(orders), _lengths(lengths), _subSpaces(_domainDim + 1), _knotVectors(knotVectors), _points(points), _id1(-1), _id2(-1)
{
	assert(knotVectors.size() == _domainDim);

	_subSpaces[0] = 1;
	for (SInt i = 0; i < _domainDim; i++) {
		assert(_orders[i] <= _lengths[i]);
		_subSpaces[i + 1] = _subSpaces[i] * _lengths[i];
	}

}

template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR>::BspMap(const BspMap &mv) :
	_domainDim(mv._domainDim), _rangeDim(mv._rangeDim), _orders(mv._orders), _lengths(mv._lengths), _subSpaces(mv._subSpaces), _knotVectors(mv._knotVectors), _points(mv._points), _id1(mv._id1), _id2(mv._id2)
{
}

template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR> & BspMap<PTypeD, PTypeR>::operator = (const BspMap &mv)
{
	if (this != &mv) {
		_domainDim = mv._domainDim;
		_rangeDim = mv._rangeDim;
		_orders = mv._orders;
		_lengths = mv._lengths;
		_subSpaces = mv._subSpaces;
		_knotVectors = mv._knotVectors;
		_points = mv._points;
		_id1 = mv._id1;
		_id2 = mv._id2;
	}
	return *this;
}


template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR>::BspMap(const BspMap<PTypeD, PTypeR>& mv, SReal t, SInt dir)
    : _domainDim(mv._domainDim - 1), _rangeDim(PTypeR::dim()), _orders(_domainDim), _lengths(_domainDim),
      _subSpaces(_domainDim + 1), _knotVectors(_domainDim), _id1(-1), _id2(-1)
{
    extractBspMap(mv, t, dir);
}

template <class PTypeD, class PTypeR>
template <class PTypeD2>
BspMap<PTypeD, PTypeR>::BspMap(const BspMap<PTypeD2, PTypeR>& mv, SReal t, SInt dir)
    : _domainDim(PTypeD::dim()), _rangeDim(PTypeR::dim()), _orders(_domainDim), _lengths(_domainDim),
      _subSpaces(_domainDim + 1), _knotVectors(_domainDim), _id1(-1), _id2(-1)
{
    assert(PTypeD2::dim() == PTypeD::dim() + 1);
    assert(dir >= 0 && dir < PTypeD2::dim());

	BspMap<PTypeD2, PTypeR> sliced(mv, t, dir);

    _orders = sliced._orders;
    _lengths = sliced._lengths;
    _knotVectors = sliced._knotVectors;
    _subSpaces = sliced._subSpaces;
    _points = sliced._points;
}


template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR>::BspMap(const BspMap<PTypeD, PTypeR> & mv,
	SInt axis) :
	_domainDim(mv._domainDim + 1), _rangeDim(PTypeR::dim()), _orders(_domainDim), _lengths(_domainDim), _subSpaces(_domainDim + 1), _knotVectors(_domainDim), _id1(-1), _id2(-1)
{
	extendBspMap(mv, axis);
}


template <class PTypeD, class PTypeR>
template <class PTypeD2>
BspMap<PTypeD, PTypeR>::BspMap(const BspMap<PTypeD2, PTypeR>& mv, SInt axis)
    : _domainDim(PTypeD::dim()), _rangeDim(PTypeR::dim()), _orders(_domainDim), _lengths(_domainDim),
      _subSpaces(_domainDim + 1), _knotVectors(_domainDim), _id1(-1), _id2(-1)
{
	
    assert(PTypeD2::dim() + 1 == PTypeD::dim());
    assert(axis >= 0 && axis <= PTypeD::dim());

	BspMap<PTypeD2, PTypeR> lifted(mv, axis);

    _orders = lifted._orders;
    _lengths = lifted._lengths;
    _knotVectors = lifted._knotVectors;
    _subSpaces = lifted._subSpaces;
    _points = lifted._points;
}


template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR>::BspMap(const BspMap<PTypeD, PTypeR> &mv, SInt newDim, SInt startAxis)
{
	assert(startAxis + mv._domainDim <= newDim);

	BspMap<PTypeD, PTypeR> tmpMap(mv);

	for (SInt i = 0; i < startAxis; i++) 
		tmpMap = BspMap<PTypeD, PTypeR>(tmpMap, 0);
	
	for (SInt i = tmpMap._domainDim; i < newDim; i++)
		tmpMap = BspMap<PTypeD, PTypeR>(tmpMap, tmpMap._domainDim);

	*this = tmpMap;
}


template <class PTypeD, class PTypeR>
void BspMap<PTypeD, PTypeR>::extendBspMap(const BspMap<PTypeD, PTypeR> & mv,
	SInt axis)
{
	assert(axis >= 0 && axis <= _domainDim);

	assert(_domainDim == mv._domainDim + 1);

	extractBspMap(mv, 0.0, -1);

	if (axis != _domainDim - 1)
		shiftAxes(axis);
}



template <class PTypeD, class PTypeR>
void BspMap<PTypeD, PTypeR>::extractBspMap(const BspMap<PTypeD, PTypeR> & mv,
	SReal t,
	SInt dir)
{
	assert(_domainDim == (dir >= 0 ? mv._domainDim - 1 : mv._domainDim + 1));

	assert(dir < mv._domainDim && _domainDim >= 1);

	assert(dir < 0 || mv.paramInDomain(dir, t));

	if (dir >= 0) {
		for (SInt i = 0; i < _domainDim; i++) {
			SInt l = i < dir ? i : i + 1;

			_lengths[i] = mv._lengths[l];
			_orders[i] = mv._orders[l];
			_knotVectors[i] = mv._knotVectors[l];
		}
	}
	else {
		std::copy(mv._lengths.begin(), mv._lengths.end(), _lengths.begin());
		_lengths[_domainDim - 1] = -dir;
		std::copy(mv._orders.begin(), mv._orders.end(), _orders.begin());
		_orders[_domainDim - 1] = -dir;
		std::copy(mv._knotVectors.begin(), mv._knotVectors.end(), _knotVectors.begin());
		uniformOpenKnotVector(_domainDim - 1);
	}

	_subSpaces[0] = 1;
	for (SInt i = 0; i < _domainDim; i++)
		_subSpaces[i + 1] = _subSpaces[i] * _lengths[i];

	_points = vector<PTypeR>(_subSpaces[_domainDim]);

	if (dir >= 0) {
		SInt indexFirst, index, newIndex;
		vector<SInt> lowerBnd(mv._domainDim, 0), upperBnd(mv._lengths), indexVec(mv._domainDim, 0), newIndexVec(_domainDim, 0);
		vector<SReal>
			basisVec = mv.coxDeBoorBasis(dir, t, indexFirst);

		upperBnd[dir] = 0;
		index = newIndex = 0;

		do {
			_points[newIndex] = mv.evaluateFromPoints(mv._subSpaces[dir], mv._orders[dir], mv._lengths[dir], basisVec, index, indexFirst);
			incremetPointsIndex(newIndexVec, newIndex);
		} while (mv.incremetPointsBoundIndex(indexVec, lowerBnd, upperBnd, index));
	}
	else {
		SInt
			oldSize = (SInt)mv._points.size();

		for (SInt j = 0; j < -dir; j++)
			std::copy(mv._points.begin(), mv._points.end(), _points.begin() + (oldSize * j));
	}
}


template <class PTypeD, class PTypeR>
PTypeR BspMap<PTypeD, PTypeR>::evaluateFromPoints(SInt inc,
	SInt order,
	SInt length,
	const vector<SReal> & basisVec,
	SInt ptsInd,
	SInt indexFirst)    const
{
	vector<SReal>::const_iterator
		basis = basisVec.begin();
	PTypeR retPt;

	if (inc == 1) {
		for (SInt i = 0; i < order; i++) {
			retPt += _points[ptsInd + indexFirst] * *basis++;
			if (++indexFirst >= length)
				indexFirst -= length;
		}
	}
	else {
		SInt
			indexFirstInc = indexFirst;

		indexFirstInc *= inc;
		for (SInt i = 0; i < order; i++) {
			retPt += _points[ptsInd + indexFirstInc] * *basis++;
			indexFirstInc += inc;
			if (++indexFirst >= length) {
				indexFirst -= length;
				indexFirstInc -= length * inc;
			}
		}
	}
	return retPt;
}


template <class PTypeD, class PTypeR>
void BspMap<PTypeD, PTypeR>::shiftAxes(SInt axis)
{
	SInt i, saveIndex;
	vector<SReal> saveKV;

	if (axis == _domainDim - 1)
		return;

	assert(axis >= 0 && axis < _domainDim);

	BspMap<PTypeD, PTypeR>
		origMV = *this;

	saveIndex = _lengths[_domainDim - 1];
	for (i = _domainDim - 1; i > axis; i--)
		_lengths[i] = _lengths[i - 1];
	_lengths[axis] = saveIndex;

	saveIndex = _orders[_domainDim - 1];
	for (i = _domainDim - 1; i > axis; i--)
		_orders[i] = _orders[i - 1];
	_orders[axis] = saveIndex;

	saveKV = _knotVectors[_domainDim - 1];
	for (i = _domainDim - 1; i > axis; i--)
		_knotVectors[i] = _knotVectors[i - 1];
	_knotVectors[axis] = saveKV;

	for (i = 1; i < _domainDim; i++)
		_subSpaces[i] = _subSpaces[i - 1] * _lengths[i - 1];

	SInt
		index = 0;
	vector<SInt> indexVec(_domainDim, 0);
	
	do {
		int revIndex;

		saveIndex = indexVec[_domainDim - 1];
		for (i = _domainDim - 1; i > axis; i--)
			indexVec[i] = indexVec[i - 1];
		indexVec[axis] = saveIndex;

		revIndex = linearizePointsIndex(indexVec);
		saveIndex = indexVec[axis];
		for (i = axis; i < _domainDim - 1; i++)
			indexVec[i] = indexVec[i + 1];
		indexVec[_domainDim - 1] = saveIndex;
		
		_points[revIndex] = origMV._points[index];

	} while (origMV.incremetPointsIndex(indexVec, index));
	
}




template <class PTypeD, class PTypeR>
SInt BspMap<PTypeD, PTypeR>::linearizePointsIndex(const vector<SInt> & IndexVec)	const
{
	SInt Index;

	switch (_domainDim) {
	case 1:
		return _subSpaces[0] * IndexVec[0];
	case 2:
		return _subSpaces[0] * IndexVec[0] +
			_subSpaces[1] * IndexVec[1];
	case 3:
		return _subSpaces[0] * IndexVec[0] +
			_subSpaces[1] * IndexVec[1] +
			_subSpaces[2] * IndexVec[2];
	case 4:
		return _subSpaces[0] * IndexVec[0] +
			_subSpaces[1] * IndexVec[1] +
			_subSpaces[2] * IndexVec[2] +
			_subSpaces[3] * IndexVec[3];
	case 5:
		return _subSpaces[0] * IndexVec[0] +
			_subSpaces[1] * IndexVec[1] +
			_subSpaces[2] * IndexVec[2] +
			_subSpaces[3] * IndexVec[3] +
			_subSpaces[4] * IndexVec[4];
	default:
		Index = 0;
		for (SInt i = 0; i < _domainDim; i++)
			Index += _subSpaces[i] * IndexVec[i];
		return Index;
	}
}


template <class PTypeD, class PTypeR>
SInt BspMap<PTypeD, PTypeR>::incremetPointsIndex(vector<SInt> & indexVec,
	SInt & index)	const
{
	SInt i;

	if (++(indexVec[0]) < _lengths[0])
		return ++index;

	indexVec[0] = 0;

	for (i = 1; i < _domainDim; i++) {
		if (++(indexVec[i]) < _lengths[i]) {
			return ++index;
		}
		indexVec[i] = 0;
	}
	return index = 0;
}


template <class PTypeD, class PTypeR>
SInt BspMap<PTypeD, PTypeR>::incremetPointsOrderIndex(vector<SInt> & indexVec,
	SInt & index)	const
{
	SInt i;

	if (++(indexVec[0]) < _orders[0])
		return ++index;

	indexVec[0] = 0;
	index -= _orders[0];

	for (i = 1; i < _domainDim; i++) {
		if (++(indexVec[i]) < _orders[i]) {
			index += _subSpaces[i];
			return ++index;
		}
		indexVec[i] = 0;
		index -= (_orders[i] - 1) * _subSpaces[i];
	}
	return index = 0;
}


template <class PTypeD, class PTypeR>
SInt BspMap<PTypeD, PTypeR>::incremetPointsBoundIndex(vector<SInt> & indexVec,
	vector<SInt> & lowerBnd,
	vector<SInt> & upperBnd,
	SInt & index)   const
{
	SInt i;

	if (++(indexVec[0]) < upperBnd[0])
		return ++index;

	indexVec[0] = lowerBnd[0];
	index -= (upperBnd[0] == 0 ? upperBnd[0] - lowerBnd[0]
		: upperBnd[0] - lowerBnd[0] - 1);

	for (i = 1; i < _domainDim; i++) {
		if (++(indexVec[i]) < upperBnd[i]) {
			return (index += _subSpaces[i]);
		}
		indexVec[i] = lowerBnd[i];
		index -= (upperBnd[i] == 0 ? upperBnd[i] - lowerBnd[i]
			: upperBnd[i] - lowerBnd[i] - 1) * _subSpaces[i];
	}
	return index = 0;
}


template <class PTypeD, class PTypeR>
SInt BspMap<PTypeD, PTypeR>::incrementPointsSkipIndex(vector<SInt> & indexVec,
	SInt skipDir,
	SInt & index)   const
{
	if (_domainDim <= 1)
		return (index = 0);

	if (++(indexVec[skipDir == 0]) < _lengths[skipDir == 0])
		return (index += _subSpaces[skipDir == 0]);
	else {
		switch (skipDir) {
		case 0:
			indexVec[1] = 0;
			index -= (_lengths[1] - 1) * _subSpaces[1];

			for (SInt i = 2; i < _domainDim; i++) {
				if (++(indexVec[i]) < _lengths[i])
					return (index += _subSpaces[i]);
				indexVec[i] = 0;
				index -= (_lengths[i] - 1) * _subSpaces[i];
			}
			break;

		case 1:
			indexVec[0] = 0;
			index -= (_lengths[0] - 1) * _subSpaces[0];

			for (SInt i = 2; i < _domainDim; i++) {
				if (++(indexVec[i]) < _lengths[i])
					return (index += _subSpaces[i]);
				indexVec[i] = 0;
				index -= (_lengths[i] - 1) * _subSpaces[i];
			}
			break;

		default:
			indexVec[0] = 0;
			index -= (_lengths[0] - 1) * _subSpaces[0];

			for (SInt i = 1; i < _domainDim; ++i == skipDir ? i++ : i) {
				if (++(indexVec[i]) < _lengths[i])
					return (index += _subSpaces[i]);
				indexVec[i] = 0;
				index -= (_lengths[i] - 1) * _subSpaces[i];
			}

		}
		return (index = 0);
	}
}



template <class PTypeD, class PTypeR>
SInt BspMap<PTypeD, PTypeR>::domainDim()    const
{
	return _domainDim;
}

template <class PTypeD, class PTypeR>
SInt BspMap<PTypeD, PTypeR>::rangeDim()	const
{
	return _rangeDim;
}

template <class PTypeD, class PTypeR>
void BspMap<PTypeD, PTypeR>::domain(SReal & tMin, SReal & tMax, SInt axis)	const
{
	assert(axis >= 0 && axis < _domainDim);
	SInt
		order = _orders[axis],
		length = _lengths[axis];

	tMin = _knotVectors[axis][order - 1];
	tMax = _knotVectors[axis][length];
}



template <class PTypeD, class PTypeR>
void BspMap<PTypeD, PTypeR>::domain(PTypeD & tMin, PTypeD & tMax) const
{
	for (SInt i = 0; i < _domainDim; i++) {
		SInt
			order = _orders[i],
			length = _lengths[i];

		tMin.setCoord(i, _knotVectors[i][order - 1]);
		tMax.setCoord(i, _knotVectors[i][length]);
	}
}


template <class PTypeD, class PTypeR>
const vector<SInt> & BspMap<PTypeD, PTypeR>::orders(void) const
{
	return _orders;
}


template <class PTypeD, class PTypeR>
const vector<SInt> & BspMap<PTypeD, PTypeR>::lengths()	const
{
	return _lengths;
}


template <class PTypeD, class PTypeR>
const vector<vector<SReal>> & BspMap<PTypeD, PTypeR>::knotVectors()	const
{
	return _knotVectors;
}


template <class PTypeD, class PTypeR>
const vector<PTypeR> & BspMap<PTypeD, PTypeR>::controlPoints()	const
{
	return _points;
}



template <class PTypeD, class PTypeR>
void  BspMap<PTypeD, PTypeR>::setDomain(SReal tMin, SReal tMax, SInt axis)
{
	knotAffineTransformOrder2(_knotVectors[axis], _orders[axis], tMin, tMax);
}


template <class PTypeD, class PTypeR>
SBool BspMap<PTypeD, PTypeR>::lengthEqOrder(SInt axis)	const
{
	return (_lengths[axis] == _orders[axis]);
}


template <class PTypeD, class PTypeR>
void BspMap<PTypeD, PTypeR>::transform(PTypeR translate, SReal scale)
{
	for (auto pt = _points.begin(); pt != _points.end(); ++pt) {
		(*pt) += translate;
		(*pt) *= scale;
	}
}

/*
template <class PTypeD, class PTypeR>
PTypeR BspMap<PTypeD, PTypeR>::evaluate(const PTypeD & params)	const
{
    SInt i, index;
    vector<SInt> indexVec(_domainDim, 0);
    vector<SInt> indexVecFirst(_domainDim);
    vector<vector<SReal>> basisVec(_domainDim);
    PTypeR retPt;
    SReal basisProd;

    assert(params.dim() == _domainDim);

    for (i = 0; i < _domainDim; i++)
        basisVec[i] = coxDeBoorBasis(i, params.coord(i), indexVecFirst[i]);

    index = linearizePointsIndex(indexVecFirst);

	
	 do {
		for (i = 0, basisProd = 1.0; i < _domainDim; i++)
			basisProd *= basisVec[i][indexVec[i]];
		retPt += _points[index] * basisProd;
	} while (incremetPointsOrderIndex(indexVec, index));
	
    return retPt;
}

*/



template <class PTypeD, class PTypeR>
PTypeR BspMap<PTypeD, PTypeR>::evaluate(const PTypeD& params) const
{
	SInt i;
	vector<SInt> indexVecFirst(_domainDim);
	vector<vector<SReal>> basisVec(_domainDim);

	assert(params.dim() == _domainDim);

	// Compute basis functions for each dimension
	for (i = 0; i < _domainDim; i++)
		basisVec[i] = coxDeBoorBasis(i, params.coord(i), indexVecFirst[i]);
	SInt total = 1;
	for (i = 0; i < _domainDim; i++)
		total *= _orders[i];
	PTypeR retPt;
	// Parallel evaluation only beneficial for larger problems
	if (total < 64) {
		// Serial version for small problems
		vector<SInt> indexVec(_domainDim, 0);
		SInt index = linearizePointsIndex(indexVecFirst);
		SReal basisProd;

		do {
			for (i = 0, basisProd = 1.0; i < _domainDim; i++)
				basisProd *= basisVec[i][indexVec[i]];
			retPt += _points[index] * basisProd;
		} while (incremetPointsOrderIndex(indexVec, index));

		return retPt;
	}
	else {
		// Parallel version for larger problems
		omp_set_num_threads(8);

#pragma omp parallel
		{
			PTypeR localRetPt;

#pragma omp for schedule(static) nowait

			for (SInt flatIdx = 0; flatIdx < total; ++flatIdx) {
				// Decompose flat index into multi-dimensional indices
				SInt remaining = flatIdx;
				SInt localIndex = 0;
				SReal localBasisProd = 1.0;

				for (SInt d = 0; d < _domainDim; ++d) {
					SInt digit = remaining % _orders[d];
					remaining /= _orders[d];

					localIndex += (indexVecFirst[d] + digit) * _subSpaces[d];
					localBasisProd *= basisVec[d][digit];
				}

				localRetPt += _points[localIndex] * localBasisProd;
			}

#pragma omp critical
			{
				retPt += localRetPt;
			}
		}
	}
	return retPt;
}

/*
template <class PTypeD, class PTypeR>
PTypeR BspMap<PTypeD, PTypeR>::evaluate(const PTypeD& params) const
{
	SInt i;
	vector<SInt> indexVecFirst(_domainDim);
	vector<vector<SReal>> basisVec(_domainDim);

	assert(params.dim() == _domainDim);

	// Compute basis vectors for each dimension
	for (i = 0; i < _domainDim; i++)
		basisVec[i] = coxDeBoorBasis(i, params.coord(i), indexVecFirst[i]);

	SInt baseIndex = linearizePointsIndex(indexVecFirst);

	// Proper initialization (works for all PntXD types)
	PTypeR retPt = PTypeR();   // default constructor, all zeros

	// ---- 1D case ----
	if (_domainDim == 1)
	{
		omp_set_num_threads(4);
#pragma omp parallel
		{
			PTypeR localRet = PTypeR();
#pragma omp for schedule(static) nowait
			for (SInt i0 = 0; i0 < _orders[0]; ++i0)
			{
				SReal b = basisVec[0][i0];
				SInt idx = baseIndex + _subSpaces[0] * i0;
				localRet += _points[idx] * b;
			}
#pragma omp critical
			retPt += localRet;
		}
	}
	// ---- 2D case ----
	else if (_domainDim == 2)
	{
		omp_set_num_threads(4);
#pragma omp parallel
		{
			PTypeR localRet = PTypeR();
#pragma omp for schedule(static) nowait
			for (SInt i0 = 0; i0 < _orders[0]; ++i0)
				for (SInt i1 = 0; i1 < _orders[1]; ++i1)
				{
					SReal b = basisVec[0][i0] * basisVec[1][i1];
					SInt idx = baseIndex + _subSpaces[0] * i0 + _subSpaces[1] * i1;
					localRet += _points[idx] * b;
				}
#pragma omp critical
			retPt += localRet;
		}
	}
	// ---- 3D case ----
	else if (_domainDim == 3)
	{
		omp_set_num_threads(8);
#pragma omp parallel
		{
			PTypeR localRet = PTypeR();
#pragma omp for schedule(static) nowait
			for (SInt i0 = 0; i0 < _orders[0]; ++i0)
				for (SInt i1 = 0; i1 < _orders[1]; ++i1)
					for (SInt i2 = 0; i2 < _orders[2]; ++i2)
					{
						SReal b = basisVec[0][i0] * basisVec[1][i1] * basisVec[2][i2];
						SInt idx = baseIndex
							+ _subSpaces[0] * i0
							+ _subSpaces[1] * i1
							+ _subSpaces[2] * i2;
						localRet += _points[idx] * b;
					}
#pragma omp critical
			retPt += localRet;
		}
	}
	// ---- Higher dimensions: fallback ----
	else
	{
		SInt totalComb = 1;
		for (i = 0; i < _domainDim; i++)
			totalComb *= _orders[i];
		omp_set_num_threads(8);
#pragma omp parallel
		{
			PTypeR localRet = PTypeR();
#pragma omp for schedule(static) nowait
			for (SInt flatIdx = 0; flatIdx < totalComb; ++flatIdx)
			{
				vector<SInt> localIdx(_domainDim);
				SInt tmp = flatIdx;
				SInt idx = baseIndex;

				for (SInt d = 0; d < _domainDim; d++)
				{
					localIdx[d] = tmp % _orders[d];
					tmp /= _orders[d];
					idx += _subSpaces[d] * localIdx[d];
				}

				SReal bProd = 1.0;
				for (SInt d = 0; d < _domainDim; d++)
					bProd *= basisVec[d][localIdx[d]];

				localRet += _points[idx] * bProd;
			}
#pragma omp critical
			retPt += localRet;
		}
	}

	return retPt;
}
*/



template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR> BspMap<PTypeD, PTypeR>::derive(SInt dir)	const
{
	SInt newLength, newOrder, i, newMVIndex,
		order = _orders[dir],
		length = _lengths[dir];
	vector<SInt> indexVec(_domainDim, 0);

	newLength = order < 2 ? length : length - 1;
	newOrder = max(order - 1, 1);

	vector<SInt> lengthsNew = _lengths;
	lengthsNew[dir] = newLength;

	vector<SInt> ordersNew = _orders;
	ordersNew[dir] = newOrder;

	BspMap<PTypeD, PTypeR>
		newMV = BspMap<PTypeD, PTypeR>( ordersNew, lengthsNew);

	for (i = 0; i < _domainDim; i++) {
		if (i == dir) {
			SInt
				offset = order < 2 ? 0 : 1;

			std::copy(_knotVectors[i].begin() + offset, _knotVectors[i].begin() + (offset + newLength + newOrder), newMV._knotVectors[i].begin());
		}
		else {			
			newMV._knotVectors[i] = _knotVectors[i];
		}
	}

	newMVIndex = 0;

	do {
		SInt
			index = linearizePointsIndex(indexVec),
			nextIndex = index + _subSpaces[dir];
		SReal denomInv,
			denom = _knotVectors[dir][indexVec[dir] + order] - _knotVectors[dir][indexVec[dir] + 1];

		if (apprxEqEps(denom, 0.0, SEEps))
			denom = SEEps;

		denomInv = 1.0 / denom;

		if (order < 2)
			newMV._points[newMVIndex] = PTypeR();
		else
			newMV._points[newMVIndex] = (_points[nextIndex] - _points[index]) * ((order - 1) * denomInv);

	} while (newMV.incremetPointsIndex(indexVec, newMVIndex));

	return newMV;
}


template <class PTypeD, class PTypeR>
void BspMap<PTypeD, PTypeR>::boundingBox(PTypeR & min, 
						PTypeR & max)	const
{
	min.fillAllCoord(SInfnty);
	max.fillAllCoord(-SInfnty);

	for (auto pt = _points.begin(); pt != _points.end(); ++pt) {
		min.min(*pt);
		max.max(*pt);
	}
}


template <class PTypeD, class PTypeR>
SBool BspMap<PTypeD, PTypeR>::verifyTDomain(SReal tMin,
	SReal tMax,
	SReal t)  const
{
	if (t < tMin)
		t += _dmnEps;
	if (t > tMax)
		t -= _dmnEps;
	if (t < tMin || t > tMax)
		return SFalse;
	else
		return STrue;
}


template <class PTypeD, class PTypeR>
void BspMap<PTypeD, PTypeR>::subdivideBezierCtrlMesh(BspMap<PTypeD, PTypeR> & lmv,
	BspMap<PTypeD, PTypeR> & rmv,
	SReal t,
	SInt dir)
{
	SReal
		t1 = 1.0 - t;

	lmv._points = _points;
	rmv._points = _points;

	vector<SInt> rIndexVec(_domainDim, 0);
	SInt
		rIndex0 = 0,
		step = _subSpaces[dir],
		length = _lengths[dir];

	do {
		SInt rIndex, rIndex1,
			lIndex = rIndex0;

		for (SInt i = 1; i < length; i++) {
			rIndex = rIndex0;

			for (SInt l = length - i; l-- > 0; ) {
				rIndex1 = rIndex + step;
				rmv._points[rIndex] = rmv._points[rIndex] * t1 + rmv._points[rIndex1] * t;
				rIndex = rIndex1;
			}
			lmv._points[lIndex] = rmv._points[rIndex0];
		}

	} while (incrementPointsSkipIndex(rIndexVec, dir, rIndex0));

}



template <class PTypeD, class PTypeR>
void BspMap<PTypeD, PTypeR>::subdivide(SReal t,
	SInt dir,
	BspMap<PTypeD, PTypeR> & lmv,
	BspMap<PTypeD, PTypeR> & rmv)	const
{
	const vector<SReal>
		& refKV = _knotVectors[dir];
	SInt i, j, index1, index2, mult,
		order = _orders[dir],
		length = _lengths[dir],
		kvLen = order + length;
	AlphaCoeff A;

	assert(dir >= 0 && dir < _domainDim);

	index1 = knotsLastIndexL(refKV, t);
	if (index1 + 1 < order)
		index1 = order - 1;
	index2 = knotsFirstIndexG(refKV, t);
	if (index2 > length)
		index2 = length;

	mult = order - 1 - (index2 - index1 - 1);

	vector<SInt>
		lengths(_lengths);

	lengths[dir] = index1 + 1;

	lmv = BspMap<PTypeD, PTypeR>(_orders, lengths);

	lengths[dir] = length - index2 + order;

	rmv = BspMap<PTypeD, PTypeR>( _orders, lengths);

	for (i = 0; i < _domainDim; i++) {
		if (i == dir) {
			std::copy(_knotVectors[i].begin(), _knotVectors[i].begin() + (index1 + 1), lmv._knotVectors[i].begin());
			for (j = index1 + 1; j <= index1 + order; j++)
				lmv._knotVectors[i][j] = t;

			std::copy(_knotVectors[i].begin() + index2, _knotVectors[i].begin() + (length + order), rmv._knotVectors[i].begin() + order);
			for (j = 0; j < order; j++)
				rmv._knotVectors[i][j] = t;
		}
		else {
			lmv._knotVectors[i] = _knotVectors[i];
			rmv._knotVectors[i] = _knotVectors[i];
		}
	}

	if (hasBezierKnotVector(refKV, _orders[dir], _lengths[dir])) {
		SReal tMin, tMax;
		SInt
			rIndex0 = 0,
			step = _subSpaces[dir];
		vector<SInt> rIndexVec(_domainDim, 0);

		domain(tMin, tMax, dir);
		t = (t - tMin) / (tMax - tMin);

		SReal t1 = 1.0 - t;

		rmv._points = _points;
		lmv._points = _points;

		typename vector<PTypeR>::iterator
			l1Points = lmv._points.begin(),
			r1Points = rmv._points.begin();

		do {
			SInt
				lIndex = rIndex0;

			for (i = 1; i < length; i++) {
				typename vector<PTypeR>::iterator
					rPts = r1Points + rIndex0;

				for (j = length - i; j-- > 0; ) {
					*rPts = *rPts * t1 + rPts[step] * t;
					rPts += step;
				}
				l1Points[lIndex += step] = r1Points[rIndex0];
			}

		} while (rmv.incrementPointsSkipIndex(rIndexVec, dir, rIndex0));
	}
	else {
		if (mult > 0) {
			SReal tMin, tMax;
			vector<SReal> newKV(mult);

			domain(tMin, tMax, dir);
			assert(verifyTDomain(tMin, tMax, t));
			if (t == tMax)
				t -= _dmnEps;
			for (i = 0; i < mult; i++)
				newKV[i] = t;

			vector<SReal>
				tmpKV = mergeKnotVectors(refKV, newKV, 0);

			A = evalAlphaCoeff(order, length, (SInt)tmpKV.size() - order, refKV, tmpKV);
		}
		else {
			A = evalAlphaCoeff(order, length, length, refKV, refKV);
		}
		mult = mult >= 0 ? 0 : -mult;

		vector<SInt> indexVec(_domainDim, 0);
		SInt
			lLength = lmv._lengths[dir],
			rLength = rmv._lengths[dir],
			index = 0;

		do {
			SInt
				lIndex = lmv.linearizePointsIndex(indexVec);
			typename vector<PTypeR>::const_iterator
				origPts = _points.begin() + index;

			alphaBlendStep(A, 0, lLength, origPts, _subSpaces[dir], lmv._points, lIndex, lmv._subSpaces[dir]);

		} while (incrementPointsSkipIndex(indexVec, dir, index));

		indexVec.assign(_domainDim, 0);
		index = 0;

		do {
			SInt
				rIndex = rmv.linearizePointsIndex(indexVec),
				offset = lLength - 1;
			typename vector<PTypeR>::const_iterator
				origPts = _points.begin() + index;

			alphaBlendStep(A, offset + mult, rLength + offset + mult, origPts, _subSpaces[dir], rmv._points, rIndex, rmv._subSpaces[dir]);

		} while (incrementPointsSkipIndex(indexVec, dir, index));
	}
}



template <class PTypeD, class PTypeR>
void BspMap<PTypeD, PTypeR>::subdivideBezierCtrlMeshOneSide(BspMap<PTypeD, PTypeR> & lmv,
	BspMap<PTypeD, PTypeR> & rmv,
	SReal t,
	SInt dir,
	SBool leftSide) const
{
	SReal
		t1 = 1.0 - t;
	SInt
		step = _subSpaces[dir],
		length = _lengths[dir];

	if (leftSide) {
		vector<SInt> lIndexVec(_domainDim, 0);
		SInt
			lIndex0 = 0;

		lmv._points = _points;

		do {
			SInt lIndex,
				lIndex00 = lIndex0;

			for (SInt i = 1; i < length; i++) {
				lIndex = lIndex00 + step * (length - i - 1);
				for (SInt l = length - i; l-- > 0; ) {
					lmv._points[lIndex + step] = lmv._points[lIndex] * t1 + lmv._points[lIndex + step] * t;
					lIndex -= step;
				}
				lIndex00 += step;
			}

		} while (incrementPointsSkipIndex(lIndexVec, dir, lIndex0));
	}
	else {
		vector<SInt> rIndexVec(_domainDim, 0);
		SInt
			rIndex0 = 0;

		rmv._points = _points;

		do {
			SInt rIndex;

			for (SInt i = 1; i < length; i++) {
				rIndex = rIndex0;
				for (SInt l = length - i; l-- > 0; ) {
					rmv._points[rIndex] = rmv._points[rIndex] * t1 + rmv._points[rIndex + step] * t;
					rIndex += step;
				}
			}

		} while (incrementPointsSkipIndex(rIndexVec, dir, rIndex0));
	}
}


template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR> BspMap<PTypeD, PTypeR>::subdivideOneSide(SReal t,
	SInt dir,
	SBool leftSide)   const
{
	const vector<SReal>
		& refKV = _knotVectors[dir];
	SInt i, j, index1, index2, mult,
		order = _orders[dir],
		length = _lengths[dir],
		kvLen = order + length;
	AlphaCoeff A;

	BspMap<PTypeD, PTypeR> lmv, rmv;

	assert(dir >= 0 && dir < _domainDim);

	index1 = knotsLastIndexL(refKV, t);
	if (index1 + 1 < order)
		index1 = order - 1;
	index2 = knotsFirstIndexG(refKV, t);
	if (index2 > length)
		index2 = length;

	mult = order - 1 - (index2 - index1 - 1);

	vector<SInt> lengths(_lengths);

	if (leftSide) {
		lengths[dir] = index1 + 1;
		lmv = BspMap<PTypeD, PTypeR>(_orders, lengths);
	}
	else {
		lengths[dir] = length - index2 + order;
		rmv = BspMap<PTypeD, PTypeR>(_orders, lengths);
	}

	for (i = 0; i < _domainDim; i++) {
		if (i == dir) {
			if (leftSide) {
				std::copy(_knotVectors[i].begin(), _knotVectors[i].begin() + (index1 + 1), lmv._knotVectors[i].begin());
				for (j = index1 + 1; j <= index1 + order; j++)
					lmv._knotVectors[i][j] = t;
			}
			else {
				std::copy(_knotVectors[i].begin() + index2, _knotVectors[i].begin() + (length + order), rmv._knotVectors[i].begin() + order);
				for (j = 0; j < order; j++)
					rmv._knotVectors[i][j] = t;
			}
		}
		else {
			if (leftSide)
				lmv._knotVectors[i] = _knotVectors[i];
			else
				rmv._knotVectors[i] = _knotVectors[i];
		}
	}

	if (hasBezierKnotVector(refKV, _orders[dir], _lengths[dir])) {
		SReal
			tMin, tMax;

		domain(tMin, tMax, dir);
		subdivideBezierCtrlMeshOneSide(lmv, rmv, (t - tMin) / (tMax - tMin), dir, leftSide);
	}
	else {
		if (mult > 0) {
			SReal tMin, tMax;
			vector<SReal> newKV(mult);

			domain(tMin, tMax, dir);
			assert(verifyTDomain(tMin, tMax, t));
			if (t == tMax)
				t -= _dmnEps;
			for (i = 0; i < mult; i++)
				newKV[i] = t;

			vector<SReal>
				tmpKV = mergeKnotVectors(refKV, newKV, 0);

			A = evalAlphaCoeff(order, length, (SInt)tmpKV.size() - order, refKV, tmpKV);
		}
		else {
			A = evalAlphaCoeff(order, length, length, refKV, refKV);
		}
		mult = mult >= 0 ? 0 : -mult;

		if (leftSide) {
			vector<SInt> indexVec(_domainDim, 0);
			SInt
				lLength = lmv._lengths[dir],
				index = 0;
			do {
				SInt
					lIndex = lmv.linearizePointsIndex(indexVec);
				typename vector<PTypeR>::const_iterator
					origPts = _points.begin() + index;

				alphaBlendStep(A, 0, lLength, origPts, _subSpaces[dir], lmv._points, lIndex, lmv._subSpaces[dir]);
			} while (incrementPointsSkipIndex(indexVec, dir, index));
		}
		else {
			vector<SInt> indexVec(_domainDim, 0);
			SInt
				rLength = rmv._lengths[dir],
				lLength = index1 + 1,
				index = 0;
			do {
				SInt
					rIndex = rmv.linearizePointsIndex(indexVec),
					offset = lLength - 1;
				typename vector<PTypeR>::const_iterator
					origPts = _points.begin() + index;

				alphaBlendStep(A, offset + mult, rLength + offset + mult, origPts, _subSpaces[dir], rmv._points, rIndex, rmv._subSpaces[dir]);
			} while (incrementPointsSkipIndex(indexVec, dir, index));
		}
	}

	if (leftSide)
		return lmv;
	else
		return rmv;
}


template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR> BspMap<PTypeD, PTypeR>::extractRegion(SReal t1, 
											 SReal t2, 
											 SInt dir)	const
{
	SReal tMin, tMax;
	SBool
		openEnd = knotsHasOpenEnd(_knotVectors[dir], _lengths[dir], _orders[dir]);

	domain(tMin, tMax, dir);

	if (t1 > t2)
		std::swap(t1, t2);

	BspMap<PTypeD, PTypeR> tmpBsp;

	if (!apprxEqEps(t1, tMin, _epsRoundKnot) || !openEnd)
		tmpBsp = subdivideOneSide(t1, dir, SFalse);
	else
		tmpBsp = *this;

	if (apprxEqEps(t2, tMax, _epsRoundKnot) && openEnd)
		return tmpBsp;
	else
		return tmpBsp.subdivideOneSide(t2, dir, STrue);
}


template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR> BspMap<PTypeD, PTypeR>::extractRegion(const PTypeD & par1, 
											 const PTypeD & par2)	const
{
	assert(par1.dim() == _domainDim);

	BspMap<PTypeD, PTypeR>
		tmpBsp = *this;

	for (SInt i = 0; i < _domainDim; i++) 
		tmpBsp = tmpBsp.extractRegion(par1.coord(i), par2.coord(i), i);
	
	return tmpBsp;
}

template <class PTypeD, class PTypeR>
void BspMap<PTypeD, PTypeR>::perturb(SReal eps)
{
	for (SInt i = 0; i < _points.size(); i++) {
		for (SInt j = 0; j < _domainDim; j++) {
			_points[i].setCoord(j, _points[i][j] + randomReal(-eps, eps));
		}
	}
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


template Pnt1D BspMap<Pnt1D,Pnt1D>::evaluate(const Pnt1D & params)	const;
template Pnt2D BspMap<Pnt1D,Pnt2D>::evaluate(const Pnt1D & params)	const;
template Pnt3D BspMap<Pnt1D,Pnt3D>::evaluate(const Pnt1D & params)	const;
template Pnt4D BspMap<Pnt1D,Pnt4D>::evaluate(const Pnt1D & params)	const;
template Pnt5D BspMap<Pnt1D,Pnt5D>::evaluate(const Pnt1D & params)	const;
template Pnt6D BspMap<Pnt1D,Pnt6D>::evaluate(const Pnt1D & params)	const;


template Pnt1D BspMap<Pnt2D,Pnt1D>::evaluate(const Pnt2D & params)	const;
template Pnt2D BspMap<Pnt2D,Pnt2D>::evaluate(const Pnt2D & params)	const;
template Pnt3D BspMap<Pnt2D,Pnt3D>::evaluate(const Pnt2D & params)	const;
template Pnt4D BspMap<Pnt2D,Pnt4D>::evaluate(const Pnt2D & params)	const;
template Pnt5D BspMap<Pnt2D,Pnt5D>::evaluate(const Pnt2D & params)	const;
template Pnt6D BspMap<Pnt2D,Pnt6D>::evaluate(const Pnt2D & params)	const;


template Pnt1D BspMap<Pnt3D,Pnt1D>::evaluate(const Pnt3D & params)	const;
template Pnt2D BspMap<Pnt3D,Pnt2D>::evaluate(const Pnt3D & params)	const;
template Pnt3D BspMap<Pnt3D,Pnt3D>::evaluate(const Pnt3D & params)	const;
template Pnt4D BspMap<Pnt3D,Pnt4D>::evaluate(const Pnt3D & params)	const;
template Pnt5D BspMap<Pnt3D,Pnt5D>::evaluate(const Pnt3D & params)	const;
template Pnt6D BspMap<Pnt3D,Pnt6D>::evaluate(const Pnt3D & params)	const;


template Pnt1D BspMap<Pnt4D,Pnt1D>::evaluate(const Pnt4D & params)	const;
template Pnt2D BspMap<Pnt4D,Pnt2D>::evaluate(const Pnt4D & params)	const;
template Pnt3D BspMap<Pnt4D,Pnt3D>::evaluate(const Pnt4D & params)	const;
template Pnt4D BspMap<Pnt4D,Pnt4D>::evaluate(const Pnt4D & params)	const;
template Pnt5D BspMap<Pnt4D,Pnt5D>::evaluate(const Pnt4D & params)	const;
template Pnt6D BspMap<Pnt4D,Pnt6D>::evaluate(const Pnt4D & params)	const;


template Pnt1D BspMap<Pnt5D,Pnt1D>::evaluate(const Pnt5D & params)	const;
template Pnt2D BspMap<Pnt5D,Pnt2D>::evaluate(const Pnt5D & params)	const;
template Pnt3D BspMap<Pnt5D,Pnt3D>::evaluate(const Pnt5D & params)	const;
template Pnt4D BspMap<Pnt5D,Pnt4D>::evaluate(const Pnt5D & params)	const;
template Pnt5D BspMap<Pnt5D,Pnt5D>::evaluate(const Pnt5D & params)	const;
template Pnt6D BspMap<Pnt5D,Pnt6D>::evaluate(const Pnt5D & params)	const;


template Pnt1D BspMap<Pnt6D,Pnt1D>::evaluate(const Pnt6D & params)	const;
template Pnt2D BspMap<Pnt6D,Pnt2D>::evaluate(const Pnt6D & params)	const;
template Pnt3D BspMap<Pnt6D,Pnt3D>::evaluate(const Pnt6D & params)	const;
template Pnt4D BspMap<Pnt6D,Pnt4D>::evaluate(const Pnt6D & params)	const;
template Pnt5D BspMap<Pnt6D,Pnt5D>::evaluate(const Pnt6D & params)	const;
template Pnt6D BspMap<Pnt6D,Pnt6D>::evaluate(const Pnt6D & params)	const;


template void BspMap<Pnt1D, Pnt1D>::domain(Pnt1D & tMin, Pnt1D & tMax)	const;
template void BspMap<Pnt1D, Pnt2D>::domain(Pnt1D & tMin, Pnt1D & tMax)	const;
template void BspMap<Pnt1D, Pnt3D>::domain(Pnt1D & tMin, Pnt1D & tMax)	const;
template void BspMap<Pnt1D, Pnt4D>::domain(Pnt1D & tMin, Pnt1D & tMax)	const;
template void BspMap<Pnt1D, Pnt5D>::domain(Pnt1D & tMin, Pnt1D & tMax)	const;
template void BspMap<Pnt1D, Pnt6D>::domain(Pnt1D & tMin, Pnt1D & tMax)	const;

template void BspMap<Pnt2D, Pnt1D>::domain(Pnt2D & tMin, Pnt2D & tMax)	const;
template void BspMap<Pnt2D, Pnt2D>::domain(Pnt2D & tMin, Pnt2D & tMax)	const;
template void BspMap<Pnt2D, Pnt3D>::domain(Pnt2D & tMin, Pnt2D & tMax)	const;
template void BspMap<Pnt2D, Pnt4D>::domain(Pnt2D & tMin, Pnt2D & tMax)	const;
template void BspMap<Pnt2D, Pnt5D>::domain(Pnt2D & tMin, Pnt2D & tMax)	const;
template void BspMap<Pnt2D, Pnt6D>::domain(Pnt2D & tMin, Pnt2D & tMax)	const;

template void BspMap<Pnt3D, Pnt1D>::domain(Pnt3D & tMin, Pnt3D & tMax)	const;
template void BspMap<Pnt3D, Pnt2D>::domain(Pnt3D & tMin, Pnt3D & tMax)	const;
template void BspMap<Pnt3D, Pnt3D>::domain(Pnt3D & tMin, Pnt3D & tMax)	const;
template void BspMap<Pnt3D, Pnt4D>::domain(Pnt3D & tMin, Pnt3D & tMax)	const;
template void BspMap<Pnt3D, Pnt5D>::domain(Pnt3D & tMin, Pnt3D & tMax)	const;
template void BspMap<Pnt3D, Pnt6D>::domain(Pnt3D & tMin, Pnt3D & tMax)	const;

template void BspMap<Pnt4D, Pnt1D>::domain(Pnt4D & tMin, Pnt4D & tMax)	const;
template void BspMap<Pnt4D, Pnt2D>::domain(Pnt4D & tMin, Pnt4D & tMax)	const;
template void BspMap<Pnt4D, Pnt3D>::domain(Pnt4D & tMin, Pnt4D & tMax)	const;
template void BspMap<Pnt4D, Pnt4D>::domain(Pnt4D & tMin, Pnt4D & tMax)	const;
template void BspMap<Pnt4D, Pnt5D>::domain(Pnt4D & tMin, Pnt4D & tMax)	const;
template void BspMap<Pnt4D, Pnt6D>::domain(Pnt4D & tMin, Pnt4D & tMax)	const;

template void BspMap<Pnt5D, Pnt1D>::domain(Pnt5D & tMin, Pnt5D & tMax)	const;
template void BspMap<Pnt5D, Pnt2D>::domain(Pnt5D & tMin, Pnt5D & tMax)	const;
template void BspMap<Pnt5D, Pnt3D>::domain(Pnt5D & tMin, Pnt5D & tMax)	const;
template void BspMap<Pnt5D, Pnt4D>::domain(Pnt5D & tMin, Pnt5D & tMax)	const;
template void BspMap<Pnt5D, Pnt5D>::domain(Pnt5D & tMin, Pnt5D & tMax)	const;
template void BspMap<Pnt5D, Pnt6D>::domain(Pnt5D & tMin, Pnt5D & tMax)	const;

template void BspMap<Pnt6D, Pnt1D>::domain(Pnt6D & tMin, Pnt6D & tMax)	const;
template void BspMap<Pnt6D, Pnt2D>::domain(Pnt6D & tMin, Pnt6D & tMax)	const;
template void BspMap<Pnt6D, Pnt3D>::domain(Pnt6D & tMin, Pnt6D & tMax)	const;
template void BspMap<Pnt6D, Pnt4D>::domain(Pnt6D & tMin, Pnt6D & tMax)	const;
template void BspMap<Pnt6D, Pnt5D>::domain(Pnt6D & tMin, Pnt6D & tMax)	const;
template void BspMap<Pnt6D, Pnt6D>::domain(Pnt6D & tMin, Pnt6D & tMax)	const;


template BspMap<Pnt1D,Pnt1D> BspMap<Pnt1D,Pnt1D>::extractRegion(const Pnt1D & par1, const Pnt1D & par2)	const;
template BspMap<Pnt1D,Pnt2D> BspMap<Pnt1D,Pnt2D>::extractRegion(const Pnt1D & par1, const Pnt1D & par2)	const;
template BspMap<Pnt1D,Pnt3D> BspMap<Pnt1D,Pnt3D>::extractRegion(const Pnt1D & par1, const Pnt1D & par2)	const;
template BspMap<Pnt1D,Pnt4D> BspMap<Pnt1D,Pnt4D>::extractRegion(const Pnt1D & par1, const Pnt1D & par2)	const;
template BspMap<Pnt1D,Pnt5D> BspMap<Pnt1D,Pnt5D>::extractRegion(const Pnt1D & par1, const Pnt1D & par2)	const;
template BspMap<Pnt1D,Pnt6D> BspMap<Pnt1D,Pnt6D>::extractRegion(const Pnt1D & par1, const Pnt1D & par2)	const;

template BspMap<Pnt2D,Pnt1D> BspMap<Pnt2D,Pnt1D>::extractRegion(const Pnt2D & par1, const Pnt2D & par2)	const;
template BspMap<Pnt2D,Pnt2D> BspMap<Pnt2D,Pnt2D>::extractRegion(const Pnt2D & par1, const Pnt2D & par2)	const;
template BspMap<Pnt2D,Pnt3D> BspMap<Pnt2D,Pnt3D>::extractRegion(const Pnt2D & par1, const Pnt2D & par2)	const;
template BspMap<Pnt2D,Pnt4D> BspMap<Pnt2D,Pnt4D>::extractRegion(const Pnt2D & par1, const Pnt2D & par2)	const;
template BspMap<Pnt2D,Pnt5D> BspMap<Pnt2D,Pnt5D>::extractRegion(const Pnt2D & par1, const Pnt2D & par2)	const;
template BspMap<Pnt2D,Pnt6D> BspMap<Pnt2D,Pnt6D>::extractRegion(const Pnt2D & par1, const Pnt2D & par2)	const;

template BspMap<Pnt3D,Pnt1D> BspMap<Pnt3D,Pnt1D>::extractRegion(const Pnt3D & par1, const Pnt3D & par2)	const;
template BspMap<Pnt3D,Pnt2D> BspMap<Pnt3D,Pnt2D>::extractRegion(const Pnt3D & par1, const Pnt3D & par2)	const;
template BspMap<Pnt3D,Pnt3D> BspMap<Pnt3D,Pnt3D>::extractRegion(const Pnt3D & par1, const Pnt3D & par2)	const;
template BspMap<Pnt3D,Pnt4D> BspMap<Pnt3D,Pnt4D>::extractRegion(const Pnt3D & par1, const Pnt3D & par2)	const;
template BspMap<Pnt3D,Pnt5D> BspMap<Pnt3D,Pnt5D>::extractRegion(const Pnt3D & par1, const Pnt3D & par2)	const;
template BspMap<Pnt3D,Pnt6D> BspMap<Pnt3D,Pnt6D>::extractRegion(const Pnt3D & par1, const Pnt3D & par2)	const;

template BspMap<Pnt4D,Pnt1D> BspMap<Pnt4D,Pnt1D>::extractRegion(const Pnt4D & par1, const Pnt4D & par2)	const;
template BspMap<Pnt4D,Pnt2D> BspMap<Pnt4D,Pnt2D>::extractRegion(const Pnt4D & par1, const Pnt4D & par2)	const;
template BspMap<Pnt4D,Pnt3D> BspMap<Pnt4D,Pnt3D>::extractRegion(const Pnt4D & par1, const Pnt4D & par2)	const;
template BspMap<Pnt4D,Pnt4D> BspMap<Pnt4D,Pnt4D>::extractRegion(const Pnt4D & par1, const Pnt4D & par2)	const;
template BspMap<Pnt4D,Pnt5D> BspMap<Pnt4D,Pnt5D>::extractRegion(const Pnt4D & par1, const Pnt4D & par2)	const;
template BspMap<Pnt4D,Pnt6D> BspMap<Pnt4D,Pnt6D>::extractRegion(const Pnt4D & par1, const Pnt4D & par2)	const;

template BspMap<Pnt5D,Pnt1D> BspMap<Pnt5D,Pnt1D>::extractRegion(const Pnt5D & par1, const Pnt5D & par2)	const;
template BspMap<Pnt5D,Pnt2D> BspMap<Pnt5D,Pnt2D>::extractRegion(const Pnt5D & par1, const Pnt5D & par2)	const;
template BspMap<Pnt5D,Pnt3D> BspMap<Pnt5D,Pnt3D>::extractRegion(const Pnt5D & par1, const Pnt5D & par2)	const;
template BspMap<Pnt5D,Pnt4D> BspMap<Pnt5D,Pnt4D>::extractRegion(const Pnt5D & par1, const Pnt5D & par2)	const;
template BspMap<Pnt5D,Pnt5D> BspMap<Pnt5D,Pnt5D>::extractRegion(const Pnt5D & par1, const Pnt5D & par2)	const;
template BspMap<Pnt5D,Pnt6D> BspMap<Pnt5D,Pnt6D>::extractRegion(const Pnt5D & par1, const Pnt5D & par2)	const;

template BspMap<Pnt6D,Pnt1D> BspMap<Pnt6D,Pnt1D>::extractRegion(const Pnt6D & par1, const Pnt6D & par2)	const;
template BspMap<Pnt6D,Pnt2D> BspMap<Pnt6D,Pnt2D>::extractRegion(const Pnt6D & par1, const Pnt6D & par2)	const;
template BspMap<Pnt6D,Pnt3D> BspMap<Pnt6D,Pnt3D>::extractRegion(const Pnt6D & par1, const Pnt6D & par2)	const;
template BspMap<Pnt6D,Pnt4D> BspMap<Pnt6D,Pnt4D>::extractRegion(const Pnt6D & par1, const Pnt6D & par2)	const;
template BspMap<Pnt6D,Pnt5D> BspMap<Pnt6D,Pnt5D>::extractRegion(const Pnt6D & par1, const Pnt6D & par2)	const;
template BspMap<Pnt6D,Pnt6D> BspMap<Pnt6D,Pnt6D>::extractRegion(const Pnt6D & par1, const Pnt6D & par2)	const;


template BspMap<Pnt1D, Pnt1D>::BspMap(const BspMap<Pnt2D, Pnt1D>&, SReal, SInt);
template BspMap<Pnt2D, Pnt1D>::BspMap(const BspMap<Pnt3D, Pnt1D>&, SReal, SInt);
template BspMap<Pnt3D, Pnt1D>::BspMap(const BspMap<Pnt4D, Pnt1D>&, SReal, SInt);
template BspMap<Pnt4D, Pnt1D>::BspMap(const BspMap<Pnt5D, Pnt1D>&, SReal, SInt);
template BspMap<Pnt5D, Pnt1D>::BspMap(const BspMap<Pnt6D, Pnt1D>&, SReal, SInt);

template BspMap<Pnt1D, Pnt2D>::BspMap(const BspMap<Pnt2D, Pnt2D>&, SReal, SInt);
template BspMap<Pnt2D, Pnt2D>::BspMap(const BspMap<Pnt3D, Pnt2D>&, SReal, SInt);
template BspMap<Pnt3D, Pnt2D>::BspMap(const BspMap<Pnt4D, Pnt2D>&, SReal, SInt);
template BspMap<Pnt4D, Pnt2D>::BspMap(const BspMap<Pnt5D, Pnt2D>&, SReal, SInt);
template BspMap<Pnt5D, Pnt2D>::BspMap(const BspMap<Pnt6D, Pnt2D>&, SReal, SInt);

template BspMap<Pnt1D, Pnt3D>::BspMap(const BspMap<Pnt2D, Pnt3D>&, SReal, SInt);
template BspMap<Pnt2D, Pnt3D>::BspMap(const BspMap<Pnt3D, Pnt3D>&, SReal, SInt);
template BspMap<Pnt3D, Pnt3D>::BspMap(const BspMap<Pnt4D, Pnt3D>&, SReal, SInt);
template BspMap<Pnt4D, Pnt3D>::BspMap(const BspMap<Pnt5D, Pnt3D>&, SReal, SInt);
template BspMap<Pnt5D, Pnt3D>::BspMap(const BspMap<Pnt6D, Pnt3D>&, SReal, SInt);

template BspMap<Pnt1D, Pnt4D>::BspMap(const BspMap<Pnt2D, Pnt4D>&, SReal, SInt);
template BspMap<Pnt2D, Pnt4D>::BspMap(const BspMap<Pnt3D, Pnt4D>&, SReal, SInt);
template BspMap<Pnt3D, Pnt4D>::BspMap(const BspMap<Pnt4D, Pnt4D>&, SReal, SInt);
template BspMap<Pnt4D, Pnt4D>::BspMap(const BspMap<Pnt5D, Pnt4D>&, SReal, SInt);
template BspMap<Pnt5D, Pnt4D>::BspMap(const BspMap<Pnt6D, Pnt4D>&, SReal, SInt);

template BspMap<Pnt1D, Pnt5D>::BspMap(const BspMap<Pnt2D, Pnt5D>&, SReal, SInt);
template BspMap<Pnt2D, Pnt5D>::BspMap(const BspMap<Pnt3D, Pnt5D>&, SReal, SInt);
template BspMap<Pnt3D, Pnt5D>::BspMap(const BspMap<Pnt4D, Pnt5D>&, SReal, SInt);
template BspMap<Pnt4D, Pnt5D>::BspMap(const BspMap<Pnt5D, Pnt5D>&, SReal, SInt);
template BspMap<Pnt5D, Pnt5D>::BspMap(const BspMap<Pnt6D, Pnt5D>&, SReal, SInt);

template BspMap<Pnt1D, Pnt6D>::BspMap(const BspMap<Pnt2D, Pnt6D>&, SReal, SInt);
template BspMap<Pnt2D, Pnt6D>::BspMap(const BspMap<Pnt3D, Pnt6D>&, SReal, SInt);
template BspMap<Pnt3D, Pnt6D>::BspMap(const BspMap<Pnt4D, Pnt6D>&, SReal, SInt);
template BspMap<Pnt4D, Pnt6D>::BspMap(const BspMap<Pnt5D, Pnt6D>&, SReal, SInt);
template BspMap<Pnt5D, Pnt6D>::BspMap(const BspMap<Pnt6D, Pnt6D>&, SReal, SInt);

template BspMap<Pnt2D, Pnt1D>::BspMap(const BspMap<Pnt1D, Pnt1D>&, SInt);
template BspMap<Pnt3D, Pnt1D>::BspMap(const BspMap<Pnt2D, Pnt1D>&, SInt);
template BspMap<Pnt4D, Pnt1D>::BspMap(const BspMap<Pnt3D, Pnt1D>&, SInt);
template BspMap<Pnt5D, Pnt1D>::BspMap(const BspMap<Pnt4D, Pnt1D>&, SInt);
template BspMap<Pnt6D, Pnt1D>::BspMap(const BspMap<Pnt5D, Pnt1D>&, SInt);

template BspMap<Pnt2D, Pnt2D>::BspMap(const BspMap<Pnt1D, Pnt2D>&, SInt);
template BspMap<Pnt3D, Pnt2D>::BspMap(const BspMap<Pnt2D, Pnt2D>&, SInt);
template BspMap<Pnt4D, Pnt2D>::BspMap(const BspMap<Pnt3D, Pnt2D>&, SInt);
template BspMap<Pnt5D, Pnt2D>::BspMap(const BspMap<Pnt4D, Pnt2D>&, SInt);
template BspMap<Pnt6D, Pnt2D>::BspMap(const BspMap<Pnt5D, Pnt2D>&, SInt);

template BspMap<Pnt2D, Pnt3D>::BspMap(const BspMap<Pnt1D, Pnt3D>&, SInt);
template BspMap<Pnt3D, Pnt3D>::BspMap(const BspMap<Pnt2D, Pnt3D>&, SInt);
template BspMap<Pnt4D, Pnt3D>::BspMap(const BspMap<Pnt3D, Pnt3D>&, SInt);
template BspMap<Pnt5D, Pnt3D>::BspMap(const BspMap<Pnt4D, Pnt3D>&, SInt);
template BspMap<Pnt6D, Pnt3D>::BspMap(const BspMap<Pnt5D, Pnt3D>&, SInt);

template BspMap<Pnt2D, Pnt4D>::BspMap(const BspMap<Pnt1D, Pnt4D>&, SInt);
template BspMap<Pnt3D, Pnt4D>::BspMap(const BspMap<Pnt2D, Pnt4D>&, SInt);
template BspMap<Pnt4D, Pnt4D>::BspMap(const BspMap<Pnt3D, Pnt4D>&, SInt);
template BspMap<Pnt5D, Pnt4D>::BspMap(const BspMap<Pnt4D, Pnt4D>&, SInt);
template BspMap<Pnt6D, Pnt4D>::BspMap(const BspMap<Pnt5D, Pnt4D>&, SInt);

template BspMap<Pnt2D, Pnt5D>::BspMap(const BspMap<Pnt1D, Pnt5D>&, SInt);
template BspMap<Pnt3D, Pnt5D>::BspMap(const BspMap<Pnt2D, Pnt5D>&, SInt);
template BspMap<Pnt4D, Pnt5D>::BspMap(const BspMap<Pnt3D, Pnt5D>&, SInt);
template BspMap<Pnt5D, Pnt5D>::BspMap(const BspMap<Pnt4D, Pnt5D>&, SInt);
template BspMap<Pnt6D, Pnt5D>::BspMap(const BspMap<Pnt5D, Pnt5D>&, SInt);

template BspMap<Pnt2D, Pnt6D>::BspMap(const BspMap<Pnt1D, Pnt6D>&, SInt);
template BspMap<Pnt3D, Pnt6D>::BspMap(const BspMap<Pnt2D, Pnt6D>&, SInt);
template BspMap<Pnt4D, Pnt6D>::BspMap(const BspMap<Pnt3D, Pnt6D>&, SInt);
template BspMap<Pnt5D, Pnt6D>::BspMap(const BspMap<Pnt4D, Pnt6D>&, SInt);
template BspMap<Pnt6D, Pnt6D>::BspMap(const BspMap<Pnt5D, Pnt6D>&, SInt);



