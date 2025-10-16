#pragma once

#include "BasicTypes.h"
#include "Pnt3D.h"
#include <assert.h>
#include <algorithm>
#include <vector>
#include <fstream>

/* TODO: replace asserts with exceptions */
using namespace std;

namespace KernelBridgeNS {

	class Pnt4D
	{
	private:
		SReal _c0;
		SReal _c1;
		SReal _c2;
		SReal _c3;
		static const SInt
			_dim = 4;
		static const SInt
			_maxOrthoIters = 10;
		static const SReal _numTol;

		SBool orthogonalize2(vector<Pnt4D> & vecs);
		SBool orthogonalize(vector<Pnt4D> & vecs);

	public:
		Pnt4D() : _c0(0.), _c1(0.), _c2(0.), _c3(0.)
		{
		}

		Pnt4D(const SReal c0,
			const SReal c1,
			const SReal c2,
			const SReal c3) : _c0(c0), _c1(c1), _c2(c2), _c3(c3)
		{
		}

		Pnt4D(const Pnt4D& Pt) : _c0(Pt._c0), _c1(Pt._c1), _c2(Pt._c2), _c3(Pt._c3)
		{
		}

		Pnt4D(const vector<Pnt4D> &vecs);


		Pnt4D(const Pnt3D & pt, SInt dimPad, SReal valuPad)
		{
			switch (dimPad) {
			case 0:
				_c0 = valuPad;
				_c1 = pt.coord(0);
				_c2 = pt.coord(1);
				_c3 = pt.coord(2);
				break;

			case 1:
				_c0 = pt.coord(0);
				_c1 = valuPad;
				_c2 = pt.coord(1);
				_c3 = pt.coord(2);
				break;

			case 2:
				_c0 = pt.coord(0);
				_c1 = pt.coord(1);
				_c2 = valuPad;
				_c3 = pt.coord(2);
				break;

			case 3:
				_c0 = pt.coord(0);
				_c1 = pt.coord(1);
				_c2 = pt.coord(2);
				_c3 = valuPad;
				break;

			default:
				assert(0);
			}
		}


		Pnt4D & operator = (const Pnt4D & Pt)
		{
			if (this == &Pt)
				return *this;

			_c0 = Pt._c0;
			_c1 = Pt._c1;
			_c2 = Pt._c2;
			_c3 = Pt._c3;
			return *this;
		}

		static SInt dim()
		{
			return _dim;
		}

		void setCoord(const SInt i,
			const SReal x)
		{
			(&_c0)[i] = x;
		}


		void incrCoord(const SInt i,
			const SReal dx)
		{
			(&_c0)[i] += dx;
		}

		void setAllCoord(const SReal c0,
			const SReal c1,
			const SReal c2,
			const SReal c3)
		{
			_c0 = c0;
			_c1 = c1;
			_c2 = c2;
			_c3 = c3;
		}

		void fillAllCoord(const SReal x)
		{
			_c0 = _c1 = _c2 = _c3 = x;
		}

		SReal coord(const SInt i)   const
		{
			return (&_c0)[i];
		}

		SReal operator [] (const SInt i)	const
		{
			return (&_c0)[i];
		}

		SReal magnitude()   const
		{
			return sqrt(_c0 * _c0
				+ _c1 * _c1
				+ _c2 * _c2
				+ _c3 * _c3);
		}

		SReal magnitudeSquare()	const
		{
			return _c0 * _c0
				+ _c1 * _c1
				+ _c2 * _c2
				+ _c3 * _c3;
		}

		SReal magnitudeL1() const
		{
			return fabs(_c0) + fabs(_c1) + fabs(_c2) + fabs(_c3);
		}

		SReal magnitudeLInf()	const
		{
			return std::max(std::max(std::max(fabs(_c0), fabs(_c1)), fabs(_c2)), fabs(_c3));
		}

		SBool isEqual(const Pnt4D & Pnt,
			const SReal Tol)	const;

		static SBool linDependent(const vector<Pnt4D> & vecs)
		{
			assert(vecs.size() == 3);

			Pnt4D corr, ortho1, ortho2, tmpVec,
				ortho0 = vecs[0];

			if (ortho0.magnitude() < SEps)
				return STrue;

			ortho1 = vecs[1] - (ortho0 * (ortho0.dot(vecs[1]) / ortho0.magnitudeSquare()));

			if (ortho1.magnitude() < SEps)
				return STrue;

			ortho2 = vecs[2] - (ortho0 * (ortho0.dot(vecs[2]) / ortho0.magnitudeSquare()))
				- (ortho1 * (ortho1.dot(vecs[2]) / ortho1.magnitudeSquare()));

			if (ortho2.magnitude() < SEps)
				return STrue;

			return SFalse;
		}


		SBool min(const Pnt4D & pt)
		{
			SBool
				updated = SFalse;

			if (pt._c0 < _c0) {
				_c0 = pt._c0;
				updated = STrue;
			}
			if (pt._c1 < _c1) {
				_c1 = pt._c1;
				updated = STrue;
			}
			if (pt._c2 < _c2) {
				_c2 = pt._c2;
				updated = STrue;
			}
			if (pt._c3 < _c3) {
				_c3 = pt._c3;
				updated = STrue;
			}
			return updated;
		}

		SBool max(const Pnt4D & pt)
		{
			SBool
				updated = SFalse;

			if (pt._c0 > _c0) {
				_c0 = pt._c0;
				updated = STrue;
			}
			if (pt._c1 > _c1) {
				_c1 = pt._c1;
				updated = STrue;
			}
			if (pt._c2 > _c2) {
				_c2 = pt._c2;
				updated = STrue;
			}
			if (pt._c3 > _c3) {
				_c3 = pt._c3;
				updated = STrue;
			}
			return updated;
		}


		void clamp(const Pnt4D & bMin, const Pnt4D & bMax)
		{
			_c0 = std::max(std::min(_c0, bMax._c0), bMin._c0);
			_c1 = std::max(std::min(_c1, bMax._c1), bMin._c1);
			_c2 = std::max(std::min(_c2, bMax._c2), bMin._c2);
			_c3 = std::max(std::min(_c3, bMax._c3), bMin._c3);
		}

		Pnt4D ceil()	const
		{
			return Pnt4D(std::ceil(_c0), std::ceil(_c1), std::ceil(_c2), std::ceil(_c3));
		}

		Pnt4D round()	const
		{
			return Pnt4D(std::round(_c0), std::round(_c1), std::round(_c2), std::round(_c3));
		}

		Pnt4D floor()	const
		{
			return Pnt4D(std::floor(_c0), std::floor(_c1), std::floor(_c2), std::floor(_c3));
		}

		Pnt4D extractDigits(SLongInt l)	const
		{
			SReal c3, c2, c1, c0;

			c3 = std::floor(l / _c3 + SEEps);
			l = l % ((SInt)_c3);
			c2 = std::floor(l / _c2 + SEEps);
			l = l % ((SInt)_c2);
			c1 = std::floor(l / _c1 + SEEps);
			c0 = (SReal)(l % ((SInt)_c1));

			return Pnt4D(c0, c1, c2, c3);
		}

		void gridNbrs(vector<Pnt4D> & nbrs)	const
		{
			assert(0);
		}

		SBool operator < (const Pnt4D & Pnt)	const
		{
			return (_c0 < Pnt._c0 && _c1 < Pnt._c1 && _c2 < Pnt._c2 && _c3 < Pnt._c3);
		}
		
		SBool operator <= (const Pnt4D & Pnt)	const
		{
			return (_c0 <= Pnt._c0 && _c1 <= Pnt._c1 && _c2 <= Pnt._c2 && _c3 <= Pnt._c3);
		}

		SBool operator > (const Pnt4D & Pnt)	const
		{
			return (_c0 > Pnt._c0 && _c1 > Pnt._c1 && _c2 > Pnt._c2 && _c3 > Pnt._c3);
		}

		SBool operator >= (const Pnt4D & Pnt)	const
		{
			return (_c0 >= Pnt._c0 && _c1 >= Pnt._c1 && _c2 >= Pnt._c2 && _c3 >= Pnt._c3);
		}


		void operator += (const Pnt4D & Pnt)
		{
			_c0 += Pnt._c0;
			_c1 += Pnt._c1;
			_c2 += Pnt._c2;
			_c3 += Pnt._c3;
		}

		Pnt4D operator + (const Pnt4D & Pnt)	const
		{
			return Pnt4D(_c0 + Pnt._c0,
				_c1 + Pnt._c1,
				_c2 + Pnt._c2,
				_c3 + Pnt._c3);
		}

		void operator -= (const Pnt4D & Pnt)
		{
			_c0 -= Pnt._c0;
			_c1 -= Pnt._c1;
			_c2 -= Pnt._c2;
			_c3 -= Pnt._c3;
		}

		Pnt4D operator - (const Pnt4D & Pnt)	const
		{
			return Pnt4D(_c0 - Pnt._c0,
				_c1 - Pnt._c1,
				_c2 - Pnt._c2,
				_c3 - Pnt._c3);
		}

		void operator *= (const SReal Scalar)
		{
			_c0 *= Scalar;
			_c1 *= Scalar;
			_c2 *= Scalar;
			_c3 *= Scalar;
		}

		Pnt4D operator * (const SReal Scalar)	const
		{
			return Pnt4D(_c0 * Scalar,
				_c1 * Scalar,
				_c2 * Scalar,
				_c3 * Scalar);
		}

		SReal operator * (const Pnt4D & Pnt)	const
		{
			return  _c0 * Pnt._c0
				+ _c1 * Pnt._c1
				+ _c2 * Pnt._c2
				+ _c3 * Pnt._c3;
		}


		SReal dot(const Pnt4D & Pnt)	const
		{
			return _c0 * Pnt._c0
				+ _c1 * Pnt._c1
				+ _c2 * Pnt._c2
				+ _c3 * Pnt._c3;
		}

		Pnt4D mult(const Pnt4D & Pnt)	const
		{
			return Pnt4D(_c0 * Pnt._c0,
				_c1 * Pnt._c1,
				_c2 * Pnt._c2,
				_c3 * Pnt._c3);
		}

		void normalize()
		{
			SReal
				mag = magnitude();

			assert(mag > SEEps);

			_c0 /= mag;
			_c1 /= mag;
			_c2 /= mag;
			_c3 /= mag;
		}

		SBool normalizeSafe()
		{
			SReal
				mag = magnitude();

			if (mag < SEEps)
				return SFalse;
			_c0 /= mag;
			_c1 /= mag;
			_c2 /= mag;
			_c3 /= mag;

			return STrue;
		}


		Pnt4D normalized()  const
		{
			SReal
				mag = magnitude();

			assert(mag > SEEps);

			return Pnt4D(_c0 / mag,
				_c1 / mag,
				_c2 / mag,
				_c3 / mag);
		}

		Pnt4D normalizedSafe()  const
		{
			SReal
				Mag = magnitude();

			if (Mag < SEEps)
				return *this;
			else
				return Pnt4D(_c0 / Mag,
					_c1 / Mag,
					_c2 / Mag,
					_c3 / Mag);
		}

		void writeToStream(std::ofstream & ofStr, const std::string & prefix, const std::string & name)
		{
			char buffer[200];

			sprintf_s(buffer, 200, "%s[OBJECT %s\n%s\t[POINT %.15f %.15f %.15f %.15f]\n]\n", prefix.c_str(), name.c_str(), prefix.c_str(), _c0, _c1, _c2, _c3);
			ofStr << buffer;
		}

	};

}



