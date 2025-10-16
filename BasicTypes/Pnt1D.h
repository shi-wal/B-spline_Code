#pragma once

#include "BasicTypes.h"
#include <assert.h>
#include <algorithm>
#include <vector>
#include <fstream>

/* TODO: replace asserts with exceptions */

using namespace std;

namespace KernelBridgeNS {

	class Pnt1D
	{
	private:
		SReal _x;
		const static SInt
			_dim = 1;

	public:
		Pnt1D() : _x(0.)
		{
		}

		Pnt1D(const SReal x) : _x(x)
		{
		}

		Pnt1D(const Pnt1D& Pt) : _x(Pt._x)
		{
		}

		Pnt1D(const std::vector<Pnt1D> &pts) : _x(0.)
		{
		}

		Pnt1D & operator = (const Pnt1D & Pt)
		{
			if (this == &Pt)
				return *this;

			_x = Pt._x;
			return *this;
		}

		static SInt dim()
		{
			return _dim;
		}

		void setCoord(const SInt i,
			const SReal x)
		{
			_x = x;
		}

		void incrCoord(const SInt i,
			const SReal dx)
		{
			(&_x)[i] += dx;
		}

		void setAllCoord(const SReal x)
		{
			_x = x;
		}

		void fillAllCoord(const SReal x)
		{
			_x = x;
		}

		SReal coord(const SInt i)   const
		{
			return (&_x)[i];
		}

		SReal operator [] (const SInt i)	const
		{
			return (&_x)[i];
		}

		SReal magnitude()   const
		{
			return fabs(_x);
		}

		SReal magnitudeSquare()	const
		{
			return _x * _x;
		}

		SReal magnitudeL1() const
		{
			return fabs(_x);
		}

		SReal magnitudeLInf()	const
		{
			return fabs(_x);
		}

		SBool isEqual(const Pnt1D & Pnt,
			const SReal Tol)	const;

		static SBool linDependent(const std::vector<Pnt1D> & vecs)
		{
			return SFalse;
		}

		SBool min(const Pnt1D & pt)
		{
			SBool
				updated = SFalse;

			if (pt._x < _x) {
				_x = pt._x;
				updated = STrue;
			}
			return updated;
		}

		SBool max(const Pnt1D & pt)
		{
			SBool
				updated = SFalse;

			if (pt._x > _x) {
				_x = pt._x;
				updated = STrue;
			}
			return updated;
		}


		void clamp(const Pnt1D & bMin, const Pnt1D & bMax)
		{
			_x = std::max(std::min(_x, bMax._x), bMin._x);
		}

		Pnt1D ceil()	const
		{
			return Pnt1D(std::ceil(_x));
		}

		Pnt1D round()	const
		{
			return Pnt1D(std::round(_x));
		}

		Pnt1D floor()	const
		{
			return Pnt1D(std::floor(_x));
		}

		Pnt1D extractDigits(SLongInt l)	const
		{
			return Pnt1D((SReal)l);
		}

		void gridNbrs(std::vector<Pnt1D> & nbrs)	const
		{
			nbrs[0] = Pnt1D(_x + 1);
		}

		SBool operator < (const Pnt1D & Pnt)	const
		{
			return (_x < Pnt._x);
		}

		SBool operator <= (const Pnt1D & Pnt)	const
		{
			return (_x <= Pnt._x);
		}

		SBool operator > (const Pnt1D & Pnt)	const
		{
			return (_x > Pnt._x);
		}

		SBool operator >= (const Pnt1D & Pnt)	const
		{
			return (_x >= Pnt._x);
		}

		void operator += (const Pnt1D & Pnt)
		{
			_x += Pnt._x;
		}

		Pnt1D operator + (const Pnt1D & Pnt)	const
		{
			return Pnt1D(_x + Pnt._x);
		}

		void operator -= (const Pnt1D & Pnt)
		{
			_x -= Pnt._x;
		}

		Pnt1D operator - (const Pnt1D & Pnt)	const
		{
			return Pnt1D(_x - Pnt._x);
		}

		void operator *= (const SReal Scalar)
		{
			_x *= Scalar;
		}

		SReal operator * (const Pnt1D & Pnt)	const
		{
			return  _x * Pnt._x;
		}

		Pnt1D operator * (const SReal Scalar)	const
		{
			return Pnt1D(_x * Scalar);
		}


		SReal dot(const Pnt1D & Pnt)	const
		{
			return _x * Pnt._x;
		}

		Pnt1D mult(const Pnt1D & Pnt)	const
		{
			return Pnt1D(_x * Pnt._x);
		}

		void normalize()
		{
			SReal
				mag = magnitude();

			assert(mag > SEEps);

			_x /= mag;
		}

		SBool normalizeSafe()
		{
			SReal
				mag = magnitude();

			if (mag < SEEps)
				return SFalse;
			_x /= mag;

			return STrue;
		}

		Pnt1D normalized()  const
		{
			SReal
				Mag = magnitude();

			assert(Mag > SEEps);

			return Pnt1D(_x / Mag);
		}

		Pnt1D normalizedSafe()  const
		{
			SReal
				Mag = magnitude();

			if (Mag < SEEps)
				return *this;
			else
				return Pnt1D(_x / Mag);
		}

		void writeToStream(std::ofstream & ofStr, const std::string & prefix, const std::string & name)
		{
			char buffer[200];

			sprintf_s(buffer, 200, "%s[OBJECT %s\n%s\t[POINT %.15f 0]\n]\n", prefix.c_str(), name.c_str(), prefix.c_str(), _x);			
			ofStr << buffer;
		}

	};

}


