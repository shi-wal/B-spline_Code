#include "Pnt4D.h"
#include <assert.h>
#include "Utils.h"

using namespace KernelBridgeNS;

const SReal Pnt4D::_numTol = 1e-10;

SBool Pnt4D::isEqual(const Pnt4D & Pnt,
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
	return STrue;
}


SBool Pnt4D::orthogonalize2(vector<Pnt4D> & vecs)		//TODO: hardcode
{
	SInt nVec = (SInt)vecs.size();
	assert(nVec == 3);
	SInt n;

	for (n = 0; n < _maxOrthoIters; n++) {
		SBool
			ortho = STrue;

		for (SInt j = 0; j < nVec; j++) {
			for (SInt i = 0; i < j; i++)
				vecs[j] -= vecs[i] * vecs[j].dot(vecs[i]);

			if (!vecs[j].normalizeSafe())
				vecs[j].fillAllCoord(0.0);
		}

		for (SInt i = 1; i < nVec; i++) {
			for (SInt j = 0; j < i; j++) {
				SReal
					r = vecs[i].dot(vecs[j]);

				if (fabs(r) > _numTol)
					ortho = SFalse;
			}
		}

		if (ortho)
			break;
	}

	if (n >= _maxOrthoIters)
		return SFalse;

	return STrue;
}


SBool Pnt4D::orthogonalize(vector<Pnt4D> & vecs)
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

		if (fabs(vecs[0].dot(vecs[1])) <= _numTol && fabs(vecs[0].dot(vecs[2])) <= _numTol && fabs(vecs[1].dot(vecs[2])) <= _numTol)
			break;
	}

	if (i >= _maxOrthoIters)
		return SFalse;

	return STrue;

}


Pnt4D::Pnt4D(const vector<Pnt4D> &vecs)
{
	std::vector<Pnt4D>
		vecsCopy = vecs;
	SBool
		tryAgain = STrue;

	orthogonalize(vecsCopy);

	while (tryAgain) {
		_c0 = randomReal(-1.0, 1.0);
		_c1 = randomReal(-1.0, 1.0);
		_c2 = randomReal(-1.0, 1.0);
		_c3 = randomReal(-1.0, 1.0);

		tryAgain = SFalse;

		(*this) -= vecsCopy[0] * dot(vecsCopy[0]);
		(*this) -= vecsCopy[1] * dot(vecsCopy[1]);
		(*this) -= vecsCopy[2] * dot(vecsCopy[2]);

		if (!normalizeSafe())
			tryAgain = STrue;

		if (std::fabs(dot(vecsCopy[0])) >= SEps ||
			std::fabs(dot(vecsCopy[1])) >= SEps ||
			std::fabs(dot(vecsCopy[2])) >= SEps)
			tryAgain = STrue;
	}
}