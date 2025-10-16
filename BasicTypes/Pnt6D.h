#pragma once

#include "BasicTypes.h"
#include "Pnt5D.h"
#include <assert.h>
#include <algorithm>
#include <vector>
#include <fstream>

using namespace std;

namespace KernelBridgeNS {

	class Pnt6D
	{
	private:
		SReal _c0;
		SReal _c1;
		SReal _c2;
		SReal _c3;
		SReal _c4;
		SReal _c5;
		static const SInt
			_dim = 6;
		static const SInt
			_maxOrthoIters = 10;
		static const SReal _numTol;
				
		SBool orthogonalize(vector<Pnt6D> & vecs);

	public:
		Pnt6D() : _c0(0.), _c1(0.), _c2(0.), _c3(0.), _c4(0.), _c5(0.)
		{
		}


		Pnt6D(const SReal c0,
			  const SReal c1,
			  const SReal c2,
			  const SReal c3,
			  const SReal c4,
			  const SReal c5) : _c0(c0), _c1(c1), _c2(c2), _c3(c3), _c4(c4), _c5(c5)
		{
		}


		Pnt6D(const Pnt6D& Pt) : _c0(Pt._c0), _c1(Pt._c1), _c2(Pt._c2), _c3(Pt._c3), _c4(Pt._c4), _c5(Pt._c5)
		{
		}


		Pnt6D(const vector<Pnt6D> &vecs);


		Pnt6D(const Pnt5D & pt, 
			  SInt dimPad, 
			  SReal valuPad)
		{
			switch (dimPad) {
			case 0:
				_c0 = valuPad;
				_c1 = pt.coord(0);
				_c2 = pt.coord(1);
				_c3 = pt.coord(2);
				_c4 = pt.coord(3);
				_c5 = pt.coord(4);
				break;

			case 1:
				_c0 = pt.coord(0);
				_c1 = valuPad;
				_c2 = pt.coord(1);
				_c3 = pt.coord(2);
				_c4 = pt.coord(3);
				_c5 = pt.coord(4);
				break;

			case 2:
				_c0 = pt.coord(0);
				_c1 = pt.coord(1);
				_c2 = valuPad;
				_c3 = pt.coord(2);
				_c4 = pt.coord(3);
				_c5 = pt.coord(4);
				break;

			case 3:
				_c0 = pt.coord(0);
				_c1 = pt.coord(1);
				_c2 = pt.coord(2);
				_c3 = valuPad;
				_c4 = pt.coord(3);
				_c5 = pt.coord(4);
				break;

			case 4:
				_c0 = pt.coord(0);
				_c1 = pt.coord(1);
				_c2 = pt.coord(2);
				_c3 = pt.coord(3);
				_c4 = valuPad;
				_c5 = pt.coord(4);
				break;

			case 5:
				_c0 = pt.coord(0);
				_c1 = pt.coord(1);
				_c2 = pt.coord(2);
				_c3 = pt.coord(3);
				_c4 = pt.coord(4);
				_c5 = valuPad;
				break;

			default:
				assert(0);
			}
		}


		Pnt6D & operator = (const Pnt6D & Pt)
		{
			if (this == &Pt)
				return *this;

			_c0 = Pt._c0;
			_c1 = Pt._c1;
			_c2 = Pt._c2;
			_c3 = Pt._c3;
			_c4 = Pt._c4;
			_c5 = Pt._c5;
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
						 const SReal c3,
						 const SReal c4,
						 const SReal c5)
		{
			_c0 = c0;
			_c1 = c1;
			_c2 = c2;
			_c3 = c3;
			_c4 = c4;
			_c5 = c5;
		}


		void fillAllCoord(const SReal x)
		{
			_c0 = _c1 = _c2 = _c3 = _c4 = _c5 = x;
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
				+ _c3 * _c3
				+ _c4 * _c4
				+ _c5 * _c5);
		}

		SReal magnitudeSquare()	const
		{
			return _c0 * _c0
				+ _c1 * _c1
				+ _c2 * _c2
				+ _c3 * _c3
				+ _c4 * _c4
				+ _c5 * _c5;
		}

		SReal magnitudeL1() const
		{
			return fabs(_c0) + fabs(_c1) + fabs(_c2) + fabs(_c3) + fabs(_c4) + fabs(_c5);
		}


		SReal magnitudeLInf()	const
		{
			return std::max(std::max(std::max(std::max(std::max(fabs(_c0), fabs(_c1)), fabs(_c2)), fabs(_c3)), fabs(_c4)), fabs(_c5));
		}


		SBool isEqual(const Pnt6D & Pnt,
			const SReal Tol)	const;


		static SBool linDependent(const vector<Pnt6D> & vecs)
		{
			assert(vecs.size() == 5);

			Pnt6D corr, ortho1, ortho2, ortho3, ortho4, tmpVec,
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

			ortho3 = vecs[3] - (ortho0 * (ortho0.dot(vecs[3]) / ortho0.magnitudeSquare()))
				- (ortho1 * (ortho1.dot(vecs[3]) / ortho1.magnitudeSquare()))
				- (ortho2 * (ortho2.dot(vecs[3]) / ortho2.magnitudeSquare()));

			if (ortho3.magnitude() < SEps)
				return STrue;

			ortho4 = vecs[4] - (ortho0 * (ortho0.dot(vecs[4]) / ortho0.magnitudeSquare()))
				- (ortho1 * (ortho1.dot(vecs[4]) / ortho1.magnitudeSquare()))
				- (ortho2 * (ortho2.dot(vecs[4]) / ortho2.magnitudeSquare()))
				- (ortho3 * (ortho3.dot(vecs[4]) / ortho3.magnitudeSquare()));

			if (ortho4.magnitude() < SEps)
				return STrue;

			return SFalse;
		}


		SBool min(const Pnt6D & pt)
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
			if (pt._c4 < _c4) {
				_c4 = pt._c4;
				updated = STrue;
			}
			if (pt._c5 < _c5) {
				_c5 = pt._c5;
				updated = STrue;
			}
			return updated;
		}

		SBool max(const Pnt6D & pt)
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
			if (pt._c4 > _c4) {
				_c4 = pt._c4;
				updated = STrue;
			}
			if (pt._c5 > _c5) {
				_c5 = pt._c5;
				updated = STrue;
			}
			return updated;
		}


		void clamp(const Pnt6D & bMin, const Pnt6D & bMax)
		{
			_c0 = std::max(std::min(_c0, bMax._c0), bMin._c0);
			_c1 = std::max(std::min(_c1, bMax._c1), bMin._c1);
			_c2 = std::max(std::min(_c2, bMax._c2), bMin._c2);
			_c3 = std::max(std::min(_c3, bMax._c3), bMin._c3);
			_c4 = std::max(std::min(_c4, bMax._c4), bMin._c4);
			_c5 = std::max(std::min(_c5, bMax._c5), bMin._c5);
		}

		Pnt6D ceil()	const
		{
			return Pnt6D(std::ceil(_c0), std::ceil(_c1), std::ceil(_c2), std::ceil(_c3), std::ceil(_c4), std::ceil(_c5));
		}

		Pnt6D round()	const
		{
			return Pnt6D(std::round(_c0), std::round(_c1), std::round(_c2), std::round(_c3), std::round(_c4), std::round(_c5));
		}

		Pnt6D floor()	const
		{
			return Pnt6D(std::floor(_c0), std::floor(_c1), std::floor(_c2), std::floor(_c3), std::floor(_c4), std::floor(_c5));
		}


		Pnt6D extractDigits(SLongInt l)	const
		{
			SReal c5, c4, c3, c2, c1, c0;

			c5 = std::floor(l / _c5 + SEEps);
			l = l % ((SInt)_c5);
			c4 = std::floor(l / _c4 + SEEps);
			l = l % ((SInt)_c4);
			c3 = std::floor(l / _c3 + SEEps);
			l = l % ((SInt)_c3);
			c2 = std::floor(l / _c2 + SEEps);
			l = l % ((SInt)_c2);
			c1 = std::floor(l / _c1 + SEEps);
			c0 = (SReal)(l % ((SInt)_c1));

			return Pnt6D(c0, c1, c2, c3, c4, c5);
		}


		void gridNbrs(vector<Pnt5D> & nbrs)	const
		{
			assert(0);
		}


		SBool operator < (const Pnt6D & Pnt)	const
		{
			return (_c0 < Pnt._c0 
				&& _c1 < Pnt._c1 
				&& _c2 < Pnt._c2 
				&& _c3 < Pnt._c3 
				&& _c4 < Pnt._c4 
				&& _c5 < Pnt._c5);
		}


		SBool operator <= (const Pnt6D & Pnt)	const
		{
			return (_c0 <= Pnt._c0 
				&& _c1 <= Pnt._c1 
				&& _c2 <= Pnt._c2 
				&& _c3 <= Pnt._c3 
				&& _c4 <= Pnt._c4 
				&& _c5 <= Pnt._c5);
		}


		SBool operator > (const Pnt6D & Pnt)	const
		{
			return (_c0 > Pnt._c0 
				&& _c1 > Pnt._c1 
				&& _c2 > Pnt._c2 
				&& _c3 > Pnt._c3 
				&& _c4 > Pnt._c4 
				&& _c5 > Pnt._c5);
		}


		SBool operator >= (const Pnt6D & Pnt)	const
		{
			return (_c0 >= Pnt._c0 
				&& _c1 >= Pnt._c1 
				&& _c2 >= Pnt._c2 
				&& _c3 >= Pnt._c3 
				&& _c4 >= Pnt._c4 
				&& _c5 >= Pnt._c5);
		}

		void operator += (const Pnt6D & Pnt)
		{
			_c0 += Pnt._c0;
			_c1 += Pnt._c1;
			_c2 += Pnt._c2;
			_c3 += Pnt._c3;
			_c4 += Pnt._c4;
			_c5 += Pnt._c5;
		}


		Pnt6D operator + (const Pnt6D & Pnt)	const
		{
			return Pnt6D(_c0 + Pnt._c0,
				_c1 + Pnt._c1,
				_c2 + Pnt._c2,
				_c3 + Pnt._c3,
				_c4 + Pnt._c4,
				_c5 + Pnt._c5);
		}


		void operator -= (const Pnt6D & Pnt)
		{
			_c0 -= Pnt._c0;
			_c1 -= Pnt._c1;
			_c2 -= Pnt._c2;
			_c3 -= Pnt._c3;
			_c4 -= Pnt._c4;
			_c5 -= Pnt._c5;
		}


		Pnt6D operator - (const Pnt6D & Pnt)	const
		{
			return Pnt6D(_c0 - Pnt._c0,
				_c1 - Pnt._c1,
				_c2 - Pnt._c2,
				_c3 - Pnt._c3,
				_c4 - Pnt._c4,
				_c5 - Pnt._c5);
		}


		void operator *= (const SReal Scalar)
		{
			_c0 *= Scalar;
			_c1 *= Scalar;
			_c2 *= Scalar;
			_c3 *= Scalar;
			_c4 *= Scalar;
			_c5 *= Scalar;
		}


		Pnt6D operator * (const SReal Scalar)	const
		{
			return Pnt6D(_c0 * Scalar,
				_c1 * Scalar,
				_c2 * Scalar,
				_c3 * Scalar,
				_c4 * Scalar,
				_c5 * Scalar);
		}


		SReal operator * (const Pnt6D & Pnt)	const
		{
			return  _c0 * Pnt._c0
				+ _c1 * Pnt._c1
				+ _c2 * Pnt._c2
				+ _c3 * Pnt._c3
				+ _c4 * Pnt._c4
				+ _c5 * Pnt._c5;
		}


		SReal dot(const Pnt6D & Pnt)	const
		{
			return _c0 * Pnt._c0
				+ _c1 * Pnt._c1
				+ _c2 * Pnt._c2
				+ _c3 * Pnt._c3
				+ _c4 * Pnt._c4
				+ _c5 * Pnt._c5;
		}


		Pnt6D mult(const Pnt6D & Pnt)	const
		{
			return Pnt6D(_c0 * Pnt._c0,
				_c1 * Pnt._c1,
				_c2 * Pnt._c2,
				_c3 * Pnt._c3,
				_c4 * Pnt._c4,
				_c5 * Pnt._c5);
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
			_c4 /= mag;
			_c5 /= mag;
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
			_c4 /= mag;
			_c5 /= mag;

			return STrue;
		}


		Pnt6D normalized()  const
		{
			SReal
				mag = magnitude();

			assert(mag > SEEps);

			return Pnt6D(_c0 / mag,
				_c1 / mag,
				_c2 / mag,
				_c3 / mag,
				_c4 / mag,
				_c5 / mag);
		}

		Pnt6D normalizedSafe()  const
		{
			SReal
				Mag = magnitude();

			if (Mag < SEEps)
				return *this;
			else
				return Pnt6D(_c0 / Mag,
					_c1 / Mag,
					_c2 / Mag,
					_c3 / Mag,
					_c4 / Mag,
					_c5 / Mag);
		}

		void writeToStream(std::ofstream & ofStr, const std::string & prefix, const std::string & name)
		{
			char buffer[200];

			sprintf_s(buffer, 200, "%s[OBJECT %s\n%s\t[POINT %.15f %.15f %.15f %.15f %.15f %.15f]\n]\n", prefix.c_str(), name.c_str(), prefix.c_str(), _c0, _c1, _c2, _c3, _c4, _c5);
			ofStr << buffer;
		}


	};
}