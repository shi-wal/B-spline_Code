#include "Pnt2D.h"
#include <assert.h>

using namespace KernelBridgeNS;



SBool Pnt2D::isEqual(const Pnt2D & Pnt,
		     const SReal Tol)	const
{    
    SReal Diff;
    Diff = _x - Pnt._x;
    if (Diff < 0) Diff = -Diff;
    if (Diff > Tol) return SFalse;
    Diff = _y - Pnt._y;
    if (Diff < 0) Diff = -Diff;
    if (Diff > Tol) return SFalse;
    return STrue;    
}
