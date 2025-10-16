#include "Cone2D.h"


using namespace KernelBridgeNS;


Cone2D::Cone2D() : _stAngle(SInfnty), _enAngle(-SInfnty), _valid(SFalse)
{
}


Cone2D::Cone2D(SReal stAngle, 
			   SReal enAngle) : _stAngle(stAngle), _enAngle(enAngle), _valid(STrue)
{
	computeTotalAngle();
}



Cone2D::Cone2D(const Pnt2D & stVec, 
			   const Pnt2D & enVec) : _stAngle(std::atan2(stVec[1], stVec[0])), _enAngle(std::atan2(enVec[1], enVec[0])), _valid(STrue)
{	
	computeTotalAngle();
}


void Cone2D::computeTotalAngle()
{
	_totalAngle = _enAngle - _stAngle;

	if (_totalAngle < 0)
		_totalAngle += 2.0 * M_PI;
	
	_obtuse = (_totalAngle > M_PI);
}



SBool Cone2D::isValid()		const
{
	return _valid;
}


SBool Cone2D::isObtuse()	const
{
	return _obtuse;
}


SBool Cone2D::vecInCone(const Pnt2D & vec)	const
{
	SReal
		vecAngle = std::atan2(vec[1], vec[0]);

	if (_stAngle < _enAngle)
		return (_stAngle <= vecAngle && vecAngle <= _enAngle);
	else
		return (_stAngle <= vecAngle || vecAngle <= _enAngle);	
}
