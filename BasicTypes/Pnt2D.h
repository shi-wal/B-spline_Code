#pragma once

#include "BasicTypes.h"
#include "Pnt1D.h"
#include <assert.h>
#include <algorithm>
#include <vector>
#include <fstream>

/* TODO: replace asserts with exceptions */

using namespace std;

namespace KernelBridgeNS {
	
	class Pnt2D
	{
	private:
		SReal _x;
		SReal _y;
		const static SInt
			_dim = 2;

	public:
		Pnt2D() : _x(0.), _y(0.)
		{
		}

		Pnt2D(const SReal x,
			const SReal y) : _x(x), _y(y)
		{
		}

		Pnt2D(const Pnt2D& Pt) : _x(Pt._x), _y(Pt._y)
		{
		}

		Pnt2D(const std::vector<Pnt2D> &pts) : _x(-pts[0]._y), _y(pts[0]._x)
		{
			assert(pts.size() == 1);
		}

		Pnt2D(const Pnt1D & pt, SInt dimPad, SReal valuPad)
		{
			switch (dimPad) {
			case 0:
				_x = valuPad;
				_y = pt.coord(0);
				break;
			
			case 1:
				_x = pt.coord(0);
				_y = valuPad;
				break;

			default:
				assert(0);
			}
		}

		Pnt2D & operator = (const Pnt2D & Pt)
		{
			if (this == &Pt)
				return *this;

			_x = Pt._x;
			_y = Pt._y;
			return *this;
		}

		static SInt dim()
		{
			return _dim;
		}

		void setCoord(const SInt i,
			const SReal x)
		{
			(&_x)[i] = x;
		}


		void incrCoord(const SInt i,
			const SReal dx)
		{
			(&_x)[i] += dx;
		}

		void setAllCoord(const SReal x,
			const SReal y)
		{
			_x = x;
			_y = y;
		}

		void fillAllCoord(const SReal x)
		{
			_x = _y = x;
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
			return sqrt(_x * _x
				+ _y * _y);
		}

		SReal magnitudeSquare()	const
		{
			return _x * _x
				+ _y * _y;
		}

		SReal magnitudeL1() const
		{
			return fabs(_x) + fabs(_y);
		}

		SReal magnitudeLInf()	const
		{
			return std::max(fabs(_x), fabs(_y));
		}

		SBool isEqual(const Pnt2D & Pnt,
			const SReal Tol)	const;

		static SBool linDependent(const std::vector<Pnt2D> & vecs)
		{
			assert(vecs.size() == 1);

			if (vecs[0].magnitude() < SEps)
				return STrue;
			else
				return SFalse;
		}

		SBool min(const Pnt2D & pt)
		{
			SBool
				updated = SFalse;

			if (pt._x < _x) {
				_x = pt._x;
				updated = STrue;
			}
			if (pt._y < _y) {
				_y = pt._y;
				updated = STrue;
			}
			return updated;
		}

		SBool max(const Pnt2D & pt)
		{
			SBool
				updated = SFalse;

			if (pt._x > _x) {
				_x = pt._x;
				updated = STrue;
			}
			if (pt._y > _y) {
				_y = pt._y;
				updated = STrue;
			}
			return updated;
		}


		void clamp(const Pnt2D & bMin, const Pnt2D & bMax)
		{
			_x = std::max(std::min(_x, bMax._x), bMin._x);
			_y = std::max(std::min(_y, bMax._y), bMin._y);
		}

		Pnt2D ceil()	const
		{
			return Pnt2D(std::ceil(_x), std::ceil(_y));
		}

		Pnt2D round()	const
		{
			return Pnt2D(std::round(_x), std::round(_y));
		}

		Pnt2D floor()	const
		{
			return Pnt2D(std::floor(_x), std::floor(_y));
		}

		Pnt2D extractDigits(SLongInt l)	const
		{
			SReal rY, rX;

			rY = std::floor(l / _y + SEEps);
			rX = (SReal)(l % ((SInt)_y));

			return Pnt2D(rX, rY);
		}

		void gridNbrs(std::vector<Pnt2D> & nbrs)	const
		{
			nbrs[0] = Pnt2D(_x + 1, _y);
			nbrs[1] = Pnt2D(_x, _y + 1);
			nbrs[2] = Pnt2D(_x + 1, _y + 1);
			nbrs[3] = Pnt2D(_x - 1, _y + 1);
		}

		SBool operator < (const Pnt2D & Pnt)	const
		{
			return (_x < Pnt._x && _y < Pnt._y);
		}

		SBool operator <= (const Pnt2D & Pnt)	const
		{
			return (_x <= Pnt._x && _y <= Pnt._y);
		}

		SBool operator > (const Pnt2D & Pnt)	const
		{
			return (_x > Pnt._x && _y > Pnt._y);
		}

		SBool operator >= (const Pnt2D & Pnt)	const
		{
			return (_x >= Pnt._x && _y >= Pnt._y);
		}

		void operator += (const Pnt2D & Pnt)
		{
			_x += Pnt._x;
			_y += Pnt._y;
		}

		Pnt2D operator + (const Pnt2D & Pnt)	const
		{
			return Pnt2D(_x + Pnt._x,
				_y + Pnt._y);
		}

		void operator -= (const Pnt2D & Pnt)
		{
			_x -= Pnt._x;
			_y -= Pnt._y;
		}

		Pnt2D operator - (const Pnt2D & Pnt)	const
		{
			return Pnt2D(_x - Pnt._x,
				_y - Pnt._y);
		}

		void operator *= (const SReal Scalar)
		{
			_x *= Scalar;
			_y *= Scalar;
		}

		SReal operator * (const Pnt2D & Pnt)	const
		{
			return  _x * Pnt._x +
				_y * Pnt._y;
		}

		Pnt2D operator * (const SReal Scalar)	const
		{
			return Pnt2D(_x * Scalar,
				_y * Scalar);
		}

		SReal operator ^ (const Pnt2D & Pnt)	const
		{
			return _x * Pnt._y
				- _y * Pnt._x;
		}

		SReal dot(const Pnt2D & Pnt)	const
		{
			return _x * Pnt._x
				+ _y * Pnt._y;
		}

		Pnt2D mult(const Pnt2D & Pnt)	const
		{
			return Pnt2D(_x * Pnt._x,
				_y * Pnt._y);
		}

		void normalize()
		{
			SReal
				mag = magnitude();

			assert(mag > SEEps);

			_x /= mag;
			_y /= mag;
		}

		SBool normalizeSafe()
		{
			SReal
				mag = magnitude();

			if (mag < SEEps)
				return SFalse;
			_x /= mag;
			_y /= mag;

			return STrue;
		}


		Pnt2D normalized()  const
		{
			SReal
				Mag = magnitude();

			assert(Mag > SEEps);

			return Pnt2D(_x / Mag,
				_y / Mag);
		}

		Pnt2D normalizedSafe()  const
		{
			SReal
				Mag = magnitude();

			if (Mag < SEEps)
				return *this;
			else
				return Pnt2D(_x / Mag,
					_y / Mag);
		}

		void writeToStream(std::ofstream & ofStr, const std::string & prefix, const std::string & name)
		{
			char buffer[200];

			sprintf_s(buffer, 200, "%s[OBJECT %s\n%s\t[POINT %.15f %.15f 0]\n]\n", prefix.c_str(), name.c_str(), prefix.c_str(), _x, _y);
			ofStr << buffer;
		}

	};

}

