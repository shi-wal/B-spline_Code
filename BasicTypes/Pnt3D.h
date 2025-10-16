#pragma once

#include "BasicTypes.h"
#include "Pnt2D.h"
#include <assert.h>
#include <algorithm>
#include <vector>
#include <fstream>
#include <string>

/* TODO: replace asserts with exceptions */
using namespace std;

namespace KernelBridgeNS {

	class Pnt3D
	{
	private:
		SReal _x;
		SReal _y;
		SReal _z;
		const static SInt
			_dim = 3;

	public:
		Pnt3D() : _x(0.), _y(0.), _z(0.)
		{
		}

		Pnt3D(const SReal x,
			const SReal y,
			const SReal z) : _x(x), _y(y), _z(z)
		{
		}

		Pnt3D(const Pnt3D& Pt) : _x(Pt._x), _y(Pt._y), _z(Pt._z)
		{
		}

		Pnt3D(const vector<Pnt3D> &pts) : _x(pts[0]._y * pts[1]._z - pts[0]._z * pts[1]._y),
										  _y(pts[0]._z * pts[1]._x - pts[0]._x * pts[1]._z),
										  _z(pts[0]._x * pts[1]._y - pts[0]._y * pts[1]._x)
		{
			assert(pts.size() == 2);
		}

		Pnt3D(const Pnt2D & pt, SInt dimPad, SReal valuPad)
		{
			switch (dimPad) {
			case 0:
				_x = valuPad;
				_y = pt.coord(0);
				_z = pt.coord(1);
				break;

			case 1:
				_x = pt.coord(0);
				_y = valuPad;
				_z = pt.coord(1);
				break;

			case 2:
				_x = pt.coord(0);
				_y = pt.coord(1);
				_z = valuPad;
				break;

			default:
				assert(0);
			}
		}

		Pnt3D & operator = (const Pnt3D & Pt)
		{
			if (this == &Pt)
				return *this;

			_x = Pt._x;
			_y = Pt._y;
			_z = Pt._z;
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
			const SReal y,
			const SReal z)
		{
			_x = x;
			_y = y;
			_z = z;
		}

		void fillAllCoord(const SReal x)
		{
			_x = _y = _z = x;
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
				+ _y * _y
				+ _z * _z);
		}

		SReal magnitudeSquare()	const
		{
			return _x * _x
				+ _y * _y
				+ _z * _z;
		}


		SReal magnitudeL1() const
		{
			return fabs(_x) + fabs(_y) + fabs(_z);
		}

		SReal magnitudeLInf()	const
		{
			return std::max(std::max(fabs(_x), fabs(_y)), fabs(_z));
		}


		SBool isEqual(const Pnt3D & Pnt,
			const SReal Tol)	const;
		
		
		static SBool linDependent(const vector<Pnt3D> & vecs)
		{
			assert(vecs.size() == 2);

			Pnt3D crossProd(vecs);

			if (crossProd.magnitude() < SEps)
				return STrue;
			else
				return SFalse;
		}


		SBool min(const Pnt3D & pt)
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
			if (pt._z < _z) {
				_z = pt._z;
				updated = STrue;
			}
			return updated;
		}

		SBool max(const Pnt3D & pt)
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
			if (pt._z > _z) {
				_z = pt._z;
				updated = STrue;
			}
			return updated;
		}

		void clamp(const Pnt3D & bMin, const Pnt3D & bMax)
		{
			_x = std::max(std::min(_x, bMax._x), bMin._x);
			_y = std::max(std::min(_y, bMax._y), bMin._y);
			_z = std::max(std::min(_z, bMax._z), bMin._z);
		}

		Pnt3D ceil()	const
		{
			return Pnt3D(std::ceil(_x), std::ceil(_y), std::ceil(_z));
		}

		Pnt3D round()	const
		{
			return Pnt3D(std::round(_x), std::round(_y), std::round(_z));
		}

		Pnt3D floor()	const
		{
			return Pnt3D(std::floor(_x), std::floor(_y), std::floor(_z));
		}

		Pnt3D extractDigits(SLongInt l)	const
		{
			SReal rZ, rY, rX;

			rZ = std::floor(l / _z + SEEps);
			l = l % ((SInt)_z);
			rY = std::floor(l / _y + SEEps);
			rX = (SReal)(l % ((SInt)_y));

			return Pnt3D(rX, rY, rZ);
		}

		void gridNbrs(vector<Pnt3D> & nbrs)	const
		{
			nbrs[0] = Pnt3D(_x + 1, _y, _z);
			nbrs[1] = Pnt3D(_x, _y + 1, _z);
			nbrs[2] = Pnt3D(_x, _y, _z + 1);
			nbrs[3] = Pnt3D(_x + 1, _y + 1, _z);
			nbrs[4] = Pnt3D(_x, _y + 1, _z + 1);
			nbrs[5] = Pnt3D(_x + 1, _y, _z + 1);
			nbrs[6] = Pnt3D(_x + 1, _y + 1, _z + 1);
			nbrs[7] = Pnt3D(_x + 1, _y - 1, _z);
			nbrs[8] = Pnt3D(_x, _y + 1, _z - 1);
			nbrs[9] = Pnt3D(_x - 1, _y, _z + 1);
			nbrs[10] = Pnt3D(_x + 1, _y + 1, _z - 1);
			nbrs[11] = Pnt3D(_x + 1, _y - 1, _z + 1);
			nbrs[12] = Pnt3D(_x - 1, _y + 1, _z + 1);
		}

		SBool operator < (const Pnt3D & Pnt)	const
		{
			return (_x < Pnt._x && _y < Pnt._y && _z < Pnt._z);
		}

		SBool operator <= (const Pnt3D & Pnt)	const
		{
			return (_x <= Pnt._x && _y <= Pnt._y && _z <= Pnt._z);
		}

		SBool operator > (const Pnt3D & Pnt)	const
		{
			return (_x > Pnt._x && _y > Pnt._y && _z > Pnt._z);
		}

		SBool operator >= (const Pnt3D & Pnt)	const
		{
			return (_x >= Pnt._x && _y >= Pnt._y && _z >= Pnt._z);
		}

		void operator += (const Pnt3D & Pnt)
		{
			_x += Pnt._x;
			_y += Pnt._y;
			_z += Pnt._z;
		}

		Pnt3D operator + (const Pnt3D & Pnt)	const
		{
			return Pnt3D(_x + Pnt._x,
				_y + Pnt._y,
				_z + Pnt._z);
		}

		void operator -= (const Pnt3D & Pnt)
		{
			_x -= Pnt._x;
			_y -= Pnt._y;
			_z -= Pnt._z;
		}

		Pnt3D operator - (const Pnt3D & Pnt)	const
		{
			return Pnt3D(_x - Pnt._x,
				_y - Pnt._y,
				_z - Pnt._z);
		}

		void operator *= (const SReal Scalar)
		{
			_x *= Scalar;
			_y *= Scalar;
			_z *= Scalar;
		}

		Pnt3D operator * (const SReal Scalar)	const
		{
			return Pnt3D(_x * Scalar,
				_y * Scalar,
				_z * Scalar);
		}

		SReal operator * (const Pnt3D & Pnt)	const
		{
			return  _x * Pnt._x
				+ _y * Pnt._y
				+ _z * Pnt._z;
		}


		SReal dot(const Pnt3D & Pnt)	const
		{
			return _x * Pnt._x
				+ _y * Pnt._y
				+ _z * Pnt._z;
		}

		Pnt3D mult(const Pnt3D & Pnt)	const
		{
			return Pnt3D(_x * Pnt._x,
				_y * Pnt._y,
				_z * Pnt._z);
		}

		void normalize()
		{
			SReal
				mag = magnitude();

			assert(mag > SEEps);
				
			_x /= mag;
			_y /= mag;
			_z /= mag;
		}

		SBool normalizeSafe()
		{
			SReal
				mag = magnitude();

			if (mag < SEEps)
				return SFalse;
			_x /= mag;
			_y /= mag;
			_z /= mag;

			return STrue;
		}

		Pnt3D normalized()  const
		{
			SReal
				mag = magnitude();

			assert(mag > SEEps);

			return Pnt3D(_x / mag,
				_y / mag,
				_z / mag);
		}

		Pnt3D normalizedSafe()  const
		{
			SReal
				Mag = magnitude();

			if (Mag < SEEps)
				return *this;
			else
				return Pnt3D(_x / Mag,
					_y / Mag,
					_z / Mag);
		}

		void writeToStream(std::ofstream & ofStr, const std::string & prefix, const std::string & name)
		{
			char buffer[200];

			sprintf_s(buffer, 200, "%s[OBJECT %s\n%s\t[POINT %.15f %.15f %.15f]\n]\n", prefix.c_str(), name.c_str(), prefix.c_str(), _x, _y, _z);
			ofStr << buffer;
		}

	};

}

