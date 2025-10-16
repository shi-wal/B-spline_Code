#include "Pnt3D.h"
#include <assert.h>

using namespace KernelBridgeNS;



SBool Pnt3D::isEqual(const Pnt3D & Pnt,
    const SReal Tol)	const
{
    SReal Diff;
    Diff = _x - Pnt._x;
    if (Diff < 0) Diff = -Diff;
    if (Diff > Tol) return SFalse;
    Diff = _y - Pnt._y;
    if (Diff < 0) Diff = -Diff;
    if (Diff > Tol) return SFalse;
    Diff = _z - Pnt._z;
    if (Diff < 0) Diff = -Diff;
    if (Diff > Tol) return SFalse;
    return STrue;
}
