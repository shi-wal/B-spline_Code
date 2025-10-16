#pragma once

#include "BasicTypes.h"
#include "Pnt2D.h"


using namespace std;

namespace KernelBridgeNS {

	class Cone2D
	{
		SReal _stAngle;
		SReal _enAngle;
		SReal _totalAngle;
		SBool _valid;
		SBool _obtuse;

		void computeTotalAngle();

	public:
		Cone2D();
		Cone2D(SReal stAngle, SReal enAngle);
		Cone2D(const Pnt2D & stVec, const Pnt2D & enVec);

		SBool isValid()		const;
		SBool isObtuse()	const;
		SBool vecInCone(const Pnt2D & vec)	const;
	};

}
