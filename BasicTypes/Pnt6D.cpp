#include "Pnt6D.h"
#include <assert.h>
#include "Utils.h"

using namespace KernelBridgeNS;
using namespace std;


const SReal Pnt6D::_numTol = 1e-10;

SBool Pnt6D::isEqual(const Pnt6D & Pnt,
					 const SReal Tol)	const
{
	SReal Diff;
	Diff = _c0 - Pnt._c0;
	if (Diff < 0)
		Diff = -Diff;
	if (Diff > Tol)
		return SFalse;
	Diff = _c1 - Pnt._c1;
	if (Diff < 0)
		Diff = -Diff;
	if (Diff > Tol)
		return SFalse;
	Diff = _c2 - Pnt._c2;
	if (Diff < 0)
		Diff = -Diff;
	if (Diff > Tol)
		return SFalse;
	Diff = _c3 - Pnt._c3;
	if (Diff < 0)
		Diff = -Diff;
	if (Diff > Tol)
		return SFalse;
	Diff = _c4 - Pnt._c4;
	if (Diff < 0)
		Diff = -Diff;
	if (Diff > Tol)
		return SFalse;
	Diff = _c5 - Pnt._c5;
	if (Diff < 0)
		Diff = -Diff;
	if (Diff > Tol)
		return SFalse;
	return STrue;
}


SBool Pnt6D::orthogonalize(vector<Pnt6D> & vecs)
{
	SInt i;

	for (i = 0; i < _maxOrthoIters; i++) {
		if (!vecs[0].normalizeSafe())
			vecs[0].fillAllCoord(0.0);

		vecs[1] -= vecs[0] * vecs[1].dot(vecs[0]);
		if (!vecs[1].normalizeSafe())
			vecs[1].fillAllCoord(0.0);

		vecs[2] -= vecs[0] * vecs[2].dot(vecs[0]);
		vecs[2] -= vecs[1] * vecs[2].dot(vecs[1]);
		if (!vecs[2].normalizeSafe())
			vecs[2].fillAllCoord(0.0);

		vecs[3] -= vecs[0] * vecs[3].dot(vecs[0]);
		vecs[3] -= vecs[1] * vecs[3].dot(vecs[1]);
		vecs[3] -= vecs[2] * vecs[3].dot(vecs[2]);
		if (!vecs[3].normalizeSafe())
			vecs[3].fillAllCoord(0.0);

		vecs[4] -= vecs[0] * vecs[4].dot(vecs[0]);
		vecs[4] -= vecs[1] * vecs[4].dot(vecs[1]);
		vecs[4] -= vecs[2] * vecs[4].dot(vecs[2]);
		vecs[4] -= vecs[3] * vecs[4].dot(vecs[3]);
		if (!vecs[4].normalizeSafe())
			vecs[4].fillAllCoord(0.0);

		if (fabs(vecs[0].dot(vecs[1])) <= _numTol
			&& fabs(vecs[0].dot(vecs[2])) <= _numTol
			&& fabs(vecs[0].dot(vecs[3])) <= _numTol
			&& fabs(vecs[0].dot(vecs[4])) <= _numTol
			&& fabs(vecs[1].dot(vecs[2])) <= _numTol
			&& fabs(vecs[1].dot(vecs[3])) <= _numTol
			&& fabs(vecs[1].dot(vecs[4])) <= _numTol
			&& fabs(vecs[2].dot(vecs[3])) <= _numTol
			&& fabs(vecs[2].dot(vecs[4])) <= _numTol
			&& fabs(vecs[3].dot(vecs[4])) <= _numTol)
			break;
	}

	if (i >= _maxOrthoIters)
		return SFalse;

	return STrue;
}


Pnt6D::Pnt6D(const vector<Pnt6D> &vecs)
{
	vector<Pnt6D>
		vecsCopy = vecs;
	SBool
		tryAgain = STrue;

	orthogonalize(vecsCopy);

	while (tryAgain) {
		_c0 = randomReal(-1.0, 1.0);
		_c1 = randomReal(-1.0, 1.0);
		_c2 = randomReal(-1.0, 1.0);
		_c3 = randomReal(-1.0, 1.0);
		_c4 = randomReal(-1.0, 1.0);
		_c5 = randomReal(-1.0, 1.0);

		tryAgain = SFalse;

		(*this) -= vecsCopy[0] * dot(vecsCopy[0]);
		(*this) -= vecsCopy[1] * dot(vecsCopy[1]);
		(*this) -= vecsCopy[2] * dot(vecsCopy[2]);
		(*this) -= vecsCopy[3] * dot(vecsCopy[3]);
		(*this) -= vecsCopy[4] * dot(vecsCopy[4]);

		if (!normalizeSafe())
			tryAgain = STrue;

		if (fabs(dot(vecsCopy[0])) >= SEps 
			|| fabs(dot(vecsCopy[1])) >= SEps
			|| fabs(dot(vecsCopy[2])) >= SEps
			|| fabs(dot(vecsCopy[3])) >= SEps
			|| fabs(dot(vecsCopy[4])) >= SEps)
			tryAgain = STrue;
	}
}