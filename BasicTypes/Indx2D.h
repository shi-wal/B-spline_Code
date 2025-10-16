#pragma once

#include "BasicTypes.h"
#include "Pnt2D.h"
#include <assert.h>
#include <vector>

using namespace std;

namespace KernelBridgeNS {

	class Indx2D
	{
	private:
		SInt _i;
		SInt _j;
		const static SInt
			_dim = 2;


	public:
		Indx2D() : _i(0), _j(0)
		{
		}

		Indx2D(const SInt i,
			const SInt j) : _i(i), _j(j)
		{
		}

		Indx2D(const Indx2D& idx) : _i(idx._i), _j(idx._j)
		{
		}	

		Indx2D(const Pnt2D & pt) : _i((SInt)pt.coord(0)), _j((SInt)pt.coord(1))
		{
		}

		Indx2D & operator = (const Indx2D & idx)
		{
			if (this == &idx)
				return *this;

			_i = idx._i;
			_j = idx._j;
			return *this;
		}

		static SInt dim()
		{
			return _dim;
		}

		void setCoord(const SInt i,
			const SInt val)
		{
			(&_i)[i] = val;
		}

		void fillAllCoord(const SInt i)
		{
			_i = _j = i;
		}

		SInt coord(const SInt i)   const
		{
			return (&_i)[i];
		}

		SInt operator [] (const SInt i)	const
		{
			return (&_i)[i];
		}

		Indx2D operator + (const Indx2D & idx)	const
		{
			return Indx2D(_i + idx._i,
				_j + idx._j);
		}

		SInt dot(const Indx2D & idx)	const
		{
			return _i * idx._i
				+ _j * idx._j;
		}

		Pnt2D operator * (const SReal Scalar)	const
		{
			return Pnt2D(_i * Scalar,
				_j * Scalar);
		}

		SBool operator < (const Indx2D & idx)	const
		{
			return (_i < idx._i && _j < idx._j);
		}

		SBool operator <= (const Indx2D & idx)	const
		{
			return (_i <= idx._i && _j <= idx._j);
		}

		SBool operator > (const Indx2D & idx)	const
		{
			return (_i > idx._i && _j > idx._j);
		}

		SBool operator >= (const Indx2D & idx)	const
		{
			return (_i >= idx._i && _j >= idx._j);
		}


		Indx2D extractDigits(SLongInt l)	const
		{			
			assert(l >= 0);

			//intY = (SInt)std::floor(l / (SReal)_j + SEEps);
			SInt
				intY = SInt(l / _j),
				intX = l % _j;

			return Indx2D(intX, intY);
		}


		void gridNbrs(vector<Indx2D> & nbrs)	const
		{
			assert(nbrs.size() >= 4);

			nbrs[0] = Indx2D(_i + 1, _j);
			nbrs[1] = Indx2D(_i, _j + 1);
			nbrs[2] = Indx2D(_i + 1, _j + 1);
			nbrs[3] = Indx2D(_i - 1, _j + 1);
		}

	};

}