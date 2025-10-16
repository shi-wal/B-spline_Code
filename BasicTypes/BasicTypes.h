#pragma once

#include <atomic>

namespace KernelBridgeNS {

    #define SFalse	false
    #define STrue	true

    typedef int    SInt;
    typedef double	SReal;
    typedef float  SShortReal;
    typedef bool    SBool;
	typedef char  SChar;
	typedef std::atomic_int	SAtomicInt;
	typedef int64_t	SLongInt;

	const SReal SPi = 3.14159265358979323846;
    const SReal	SEps = 1e-5;
    const SReal SEEps = 1e-14;
    const SReal SDmnEps = SEEps * 10;
    const SReal SKnotEps = 1e-12;
    const SReal SInfnty = 2.3197171528332553e+25;

    enum ConstraintType {
		CNSTRNT_TYPE_ZERO,
		CNSTRNT_TYPE_POSITIVE,
		CNSTRNT_TYPE_NEGATIVE,
		CNSTRNT_TYPE_ZERO_AUX
    };

	enum ConstraintClass {
		CNSTRNT_BSP,
		CNSTRNT_BSPSUB
	};

	#define M_PI       3.14159265358979323846   // pi

}
