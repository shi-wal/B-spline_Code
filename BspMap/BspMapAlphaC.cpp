#include "BspMap.h"
#include "Utils.h"
#include <assert.h>

using namespace KernelBridgeNS;
using namespace SweepNS;
using namespace std;

template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR>::AlphaCoeff::AlphaCoeff()
{
}


template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR>::AlphaCoeff::AlphaCoeff(SInt order,
	SInt lenKVT,
	SInt lenKVt) : _order(order), _length(lenKVT), _refLength(lenKVt),
	_mat(_length + 1, vector<SReal>(_refLength + 1, 0)), _matTransp(_refLength + 1, vector<SReal>(_length + 1, 0)),
	_colInd(_refLength), _colLen(_refLength)
{
}


template <class PTypeD, class PTypeR>
typename BspMap<PTypeD, PTypeR>::AlphaCoeff BspMap<PTypeD, PTypeR>::evalAlphaCoeff(SInt order,
	SInt lenKVT,
	SInt lenKVt,
	const vector<SReal> & KVT,
	const vector<SReal> & KVt)    const
{
	AlphaCoeff alphaC(order, lenKVT, lenKVt);
	SInt i, j, o, nextStart;
	vector<SInt> noZeroKVtMin(alphaC._length + 2), noZeroKVtMax(alphaC._length + 2);

	for (i = nextStart = 0; i < alphaC._length; i++) {
		SInt
			last = -1,
			first = -1;
		vector<SReal>::iterator
			rowI = alphaC._mat[i].begin() + nextStart;
		vector<SReal>::const_iterator
			KVtp = KVt.begin() + nextStart;
		SReal
			KVTI = KVT[i],
			KVTI1 = KVT[i + 1];

		for (j = nextStart - 1; ++j < alphaC._refLength; rowI++, KVtp++) {
			if (*KVtp >= KVTI1) {
				if (first >= 0)
					last = j - 1;
				break;
			}
			else if (KVTI <= *KVtp) {
				*rowI = 1.0;
				if (first < 0) {
					first = j;
					if (nextStart < first)
						nextStart = first;
				}
			}
		}
		if (j >= alphaC._refLength)
			last = alphaC._refLength - 1;

		noZeroKVtMin[i] = first < 0 ? alphaC._refLength + 2 : first;
		noZeroKVtMax[i] = last;
	}

	for (o = 2; o <= alphaC._order; o++) {
		for (i = 0; i < alphaC._length; i++) {
			SInt
				jMin = noZeroKVtMin[i] = min(noZeroKVtMin[i], noZeroKVtMin[i + 1]),
				jMax = noZeroKVtMax[i] = max(noZeroKVtMax[i], noZeroKVtMax[i + 1]);

			if (jMax >= jMin) {
				vector<SReal>::iterator
					rowI = alphaC._mat[i].begin() + jMin,
					rowI1 = alphaC._mat[i + 1].begin() + jMin;
				vector<SReal>::const_iterator
					KVtp = KVt.begin() + jMin + o - 1;
				SReal
					KVTI = KVT[i],
					KVTIO = KVT[i + o],
					t1 = KVT[i + o - 1] - KVTI,
					t2 = KVTIO - KVT[i + 1];

				t1 = t1 < _Eps ? 0.0 : 1.0 / t1;
				t2 = t2 < _Eps ? 0.0 : 1.0 / t2;

				for (j = jMin; j <= jMax; j++, rowI++, rowI1++, KVtp++) {
					*rowI = *rowI * (*KVtp - KVTI) * t1 + *rowI1 * (KVTIO - *KVtp) * t2;
				}
			}
		}
	}

	nextStart = alphaC._length - 1;

	for (i = alphaC._refLength; --i >= 0; ) {
		SInt
			first = 0,
			last = -1;

		for (j = nextStart + 1; --j >= 0; ) {
			if ((alphaC._matTransp[i][j] = alphaC._mat[j][i]) > _Eps) {
				first = last < 0 ? last = j : j;
			}
			else if (last >= 0)
				break;
		}
		alphaC._colInd[i] = first;
		alphaC._colLen[i] = last - first + 1;
		if (nextStart > last && last > -1)
			nextStart = last;
	}
	return alphaC;
}


template <class PTypeD, class PTypeR>
void BspMap<PTypeD, PTypeR>::alphaBlendStep(AlphaCoeff & alphaC,
	SInt iMin,
	SInt iMax,
	typename vector<PTypeR>::const_iterator & origPts,
	SInt origPtsStep,
	vector<PTypeR> & refPtsVec,
	SInt refPtsSt,
	SInt refPtsStep)  const
{
	vector<SInt>::const_iterator
		colLen = alphaC._colLen.begin() + iMin,
		colInd = alphaC._colInd.begin() + iMin;
	typename vector<PTypeR>::iterator
		refPts = refPtsVec.begin() + refPtsSt;
	SInt
		refPtsLen = (SInt)refPtsVec.size();
	vector<SReal>::const_iterator r;
	typename vector<PTypeR>::const_iterator p;

	for (SInt i = iMin; i < iMax; i++) {
		switch (*colLen++) {
		case 1:
			*refPts = origPts[*colInd++ * origPtsStep];
			break;
		case 2:
			p = origPts + *colInd * origPtsStep;
			r = alphaC._matTransp[i].begin() + *colInd++;
			*refPts = *p * *r +
				*(p + origPtsStep) * *(r + 1);
			break;
		case 3:
			p = origPts + *colInd * origPtsStep;
			r = alphaC._matTransp[i].begin() + *colInd++;
			*refPts = *p * *r +
				*(p + origPtsStep) * *(r + 1) +
				*(p + origPtsStep * 2) * *(r + 2);
			break;
		case 4:
			p = origPts + *colInd * origPtsStep;
			r = alphaC._matTransp[i].begin() + *colInd++;
			*refPts = *p * *r +
				*(p + origPtsStep) * *(r + 1) +
				*(p + origPtsStep * 2) * *(r + 2) +
				*(p + origPtsStep * 3) * *(r + 3);
			break;
		default:
			SInt
				len = *(colLen - 1),
				rangeDim = origPts[0].dim();

			p = origPts + *colInd * origPtsStep;
			r = alphaC._matTransp[i].begin() + *colInd++;

			PTypeR newPt;
			for (SInt j = 0; len-- > 0; j++/*p += origPtsStep*/)
				newPt += *(p + origPtsStep * j) * *r++;

			*refPts = newPt;
			break;
		}
		refPts += min(refPtsStep, refPtsLen - refPtsSt);
		refPtsSt += min(refPtsStep, refPtsLen - refPtsSt);
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