#include "Pnt1D.h"
#include <assert.h>

using namespace KernelBridgeNS;



SBool Pnt1D::isEqual(const Pnt1D & Pnt,
    const SReal Tol)	const
{
    SReal Diff;
    Diff = _x - Pnt._x;
    if (Diff < 0) Diff = -Diff;
    if (Diff > Tol) return SFalse;    
    return STrue;
}