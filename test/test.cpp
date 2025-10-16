/*
#include <iostream>
#include <vector>
#include <iomanip>
#include "BspMap.h"
#include <chrono>

using namespace SweepNS;
using namespace std;

inline double genVal(int idx, int comp) {
    return 0.1 * (idx + comp + 1);
}

void test_deterministic_point_eval() {
    cout << fixed << setprecision(3);

    for (int dim = 1; dim <= 3; ++dim) {  // test for 1D–3D
        vector<SInt> orders(dim, 4);       // order = 4
        vector<SInt> lens(dim, 8);         // length = 8
        int total_points = 1;
        for (auto l : lens) total_points *= l;

        cout << "\n=== " << dim << "D B-spline Evaluation ===\n";
        cout << "Orders: [";
        for (int i = 0; i < dim; ++i)
            cout << (i ? ", " : "") << orders[i];
        cout << "], Lengths: [";
        for (int i = 0; i < dim; ++i)
            cout << (i ? ", " : "") << lens[i];
        cout << "], Total control points: " << total_points << "\n";

        if (dim == 1) {
            vector<Pnt1D> pts;
            for (int i = 0; i < total_points; ++i)
                pts.emplace_back(genVal(i, 0));

            BspMap<Pnt1D, Pnt1D> curve(orders, lens, pts);

            // evaluate at midpoint
            const auto& kv = curve.knotVectors()[0];
            double t = 0.5 * (kv[curve.orders()[0] - 1] + kv[curve.lengths()[0]]);
            Pnt1D eval = curve.evaluate(Pnt1D(t));

            cout << "Parameter t = " << t << "\n";
            cout << "Evaluated point: (" << eval.coord(0) << ")\n";
        }

        else if (dim == 2) {
            vector<Pnt2D> pts;
            for (int i = 0; i < total_points; ++i)
                pts.emplace_back(genVal(i, 0), genVal(i, 1));

            BspMap<Pnt2D, Pnt2D> surface(orders, lens, pts);

            vector<double> t(2);
            for (int d = 0; d < 2; ++d) {
                const auto& kv = surface.knotVectors()[d];
                t[d] = 0.5 * (kv[surface.orders()[d] - 1] + kv[surface.lengths()[d]]);
            }

            Pnt2D eval = surface.evaluate(Pnt2D(t[0], t[1]));

            cout << "Parameters (u,v) = (" << t[0] << ", " << t[1] << ")\n";
            cout << "Evaluated point: (" << eval.coord(0)
                << ", " << eval.coord(1) << ")\n";
        }

        else if (dim == 3) {
            vector<Pnt3D> pts;
            for (int i = 0; i < total_points; ++i)
                pts.emplace_back(genVal(i, 0), genVal(i, 1), genVal(i, 2));

            BspMap<Pnt3D, Pnt3D> volume(orders, lens, pts);

            vector<double> t(3);
            for (int d = 0; d < 3; ++d) {
                const auto& kv = volume.knotVectors()[d];
                t[d] = 0.5 * (kv[volume.orders()[d] - 1] + kv[volume.lengths()[d]]);
            }

            Pnt3D eval = volume.evaluate(Pnt3D(t[0], t[1], t[2]));

            cout << "Parameters (u,v,w) = (" << t[0] << ", " << t[1]
                << ", " << t[2] << ")\n";
            cout << "Evaluated point: (" << eval.coord(0)
                << ", " << eval.coord(1)
                << ", " << eval.coord(2) << ")\n";
        }
    }
}

int main() {
    auto start = chrono::high_resolution_clock::now();

    test_deterministic_point_eval();

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    cout << "\n Serial Version..." << endl;
    cout << "\nTotal execution time: " << elapsed.count() << " seconds\n";
}
*/

#include <iostream>
#include <vector>
#include <iomanip>
#include "BspMap.h"
#include <chrono>
#include <string>

using namespace SweepNS;
using namespace std;

// Deterministic value generator
inline double genVal(int idx, int comp) {
    return 0.1 * (idx + comp + 1);
}

// Helper to get parameter values deterministically
vector<double> getParameterValues(int domainDim, const vector<vector<double>>& knotVectors,
    const vector<SInt>& orders, const vector<SInt>& lengths,
    int testIndex) {
    vector<double> params(domainDim);
    for (int d = 0; d < domainDim; ++d) {
        const auto& kv = knotVectors[d];
        double start = kv[orders[d] - 1];
        double end = kv[lengths[d]];
        // Distribute test points evenly across parameter space
        double t = start + (end - start) * (testIndex + 1) / 11.0;
        params[d] = t;
    }
    return params;
}

// Template function to create point of any dimension
template<typename T>
T createPoint(const vector<double>& coords);

template<>
Pnt1D createPoint<Pnt1D>(const vector<double>& coords) {
    return Pnt1D(coords[0]);
}

template<>
Pnt2D createPoint<Pnt2D>(const vector<double>& coords) {
    return Pnt2D(coords[0], coords[1]);
}

template<>
Pnt3D createPoint<Pnt3D>(const vector<double>& coords) {
    return Pnt3D(coords[0], coords[1], coords[2]);
}

// Test function for domain dimension 1
template<typename RangeType>
void testMapping1D(int rangeDim, int testCount) {
    vector<SInt> orders(1, 4);
    vector<SInt> lens(1, 512);
    int total_points = lens[0];

    vector<RangeType> pts;
    for (int i = 0; i < total_points; ++i) {
        vector<double> coords(rangeDim);
        for (int c = 0; c < rangeDim; ++c) {
            coords[c] = genVal(i, c);
        }
        pts.push_back(createPoint<RangeType>(coords));
    }

    BspMap<Pnt1D, RangeType> bspline(orders, lens, pts);

// #pragma omp parallel for
    for (int t = 0; t < testCount; ++t) {
        vector<double> params = getParameterValues(1, bspline.knotVectors(),
            bspline.orders(), bspline.lengths(), t);
        Pnt1D param = createPoint<Pnt1D>(params);
        RangeType result = bspline.evaluate(param);
    }
}

// Test function for domain dimension 2
template<typename RangeType>
void testMapping2D(int rangeDim, int testCount) {
    vector<SInt> orders(2, 4);
    vector<SInt> lens(2, 256);
    int total_points = lens[0] * lens[1];

    vector<RangeType> pts;
    for (int i = 0; i < total_points; ++i) {
        vector<double> coords(rangeDim);
        for (int c = 0; c < rangeDim; ++c) {
            coords[c] = genVal(i, c);
        }
        pts.push_back(createPoint<RangeType>(coords));
    }

    BspMap<Pnt2D, RangeType> bspline(orders, lens, pts);

    for (int t = 0; t < testCount; ++t) {
        vector<double> params = getParameterValues(2, bspline.knotVectors(),
            bspline.orders(), bspline.lengths(), t);
        Pnt2D param = createPoint<Pnt2D>(params);
        RangeType result = bspline.evaluate(param);
    }
}

// Test function for domain dimension 3
template<typename RangeType>
void testMapping3D(int rangeDim, int testCount) {
    vector<SInt> orders(3, 4);
    vector<SInt> lens(3, 128);
    int total_points = lens[0] * lens[1] * lens[2];

    vector<RangeType> pts;
    for (int i = 0; i < total_points; ++i) {
        vector<double> coords(rangeDim);
        for (int c = 0; c < rangeDim; ++c) {
            coords[c] = genVal(i, c);
        }
        pts.push_back(createPoint<RangeType>(coords));
    }

    BspMap<Pnt3D, RangeType> bspline(orders, lens, pts);

    for (int t = 0; t < testCount; ++t) {
        vector<double> params = getParameterValues(3, bspline.knotVectors(),
            bspline.orders(), bspline.lengths(), t);
        Pnt3D param = createPoint<Pnt3D>(params);
        RangeType result = bspline.evaluate(param);
    }
}

void runAllTests() {
    cout << fixed << setprecision(6);
    const int testsPerMapping = 10;
    int totalTests = 0;

    cout << "\n========================================\n";
    cout << "B-spline Evaluation Test Suite\n";
    cout << "Testing all 9 domain->range combinations\n";
    cout << "Tests per mapping: " << testsPerMapping << "\n";
    cout << "========================================\n\n";

    // Mapping descriptions
    vector<string> descriptions = {
        "1D->1D: Curve in 1D space (scalar function)",
        "1D->2D: Curve in 2D space (planar curve)",
        "1D->3D: Curve in 3D space (space curve)",
        "2D->1D: Surface to scalar (height field)",
        "2D->2D: Surface in 2D space",
        "2D->3D: Surface in 3D space",
        "3D->1D: Volume to scalar (density field)",
        "3D->2D: Volume to 2D (projection-like)",
        "3D->3D: Volume in 3D space"
    };

    int mappingIdx = 0;

    // Domain dimension 1
    for (int rangeDim = 1; rangeDim <= 3; ++rangeDim) {
        cout << "Testing " << descriptions[mappingIdx++] << "\n";
        auto start = chrono::high_resolution_clock::now();

        if (rangeDim == 1) {
            testMapping1D<Pnt1D>(1, testsPerMapping);
        }
        else if (rangeDim == 2) {
            testMapping1D<Pnt2D>(2, testsPerMapping);
        }
        else {
            testMapping1D<Pnt3D>(3, testsPerMapping);
        }

        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = end - start;
        totalTests += testsPerMapping;
        cout << "  Completed " << testsPerMapping << " tests in "
            << elapsed.count() << " seconds\n\n";
    }

    // Domain dimension 2
    for (int rangeDim = 1; rangeDim <= 3; ++rangeDim) {
        cout << "Testing " << descriptions[mappingIdx++] << "\n";
        auto start = chrono::high_resolution_clock::now();

        if (rangeDim == 1) {
            testMapping2D<Pnt1D>(1, testsPerMapping);
        }
        else if (rangeDim == 2) {
            testMapping2D<Pnt2D>(2, testsPerMapping);
        }
        else {
            testMapping2D<Pnt3D>(3, testsPerMapping);
        }

        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = end - start;
        totalTests += testsPerMapping;
        cout << "  Completed " << testsPerMapping << " tests in "
            << elapsed.count() << " seconds\n\n";
    }

    // Domain dimension 3
    for (int rangeDim = 1; rangeDim <= 3; ++rangeDim) {
        cout << "Testing " << descriptions[mappingIdx++] << "\n";
        auto start = chrono::high_resolution_clock::now();

        if (rangeDim == 1) {
            testMapping3D<Pnt1D>(1, testsPerMapping);
        }
        else if (rangeDim == 2) {
            testMapping3D<Pnt2D>(2, testsPerMapping);
        }
        else {
            testMapping3D<Pnt3D>(3, testsPerMapping);
        }

        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = end - start;
        totalTests += testsPerMapping;
        cout << "  Completed " << testsPerMapping << " tests in "
            << elapsed.count() << " seconds\n\n";
    }

    cout << "========================================\n";
    cout << "Total tests executed: " << totalTests << "\n";
    cout << "========================================\n";
}

int main() {
   cout << "\n=== SERIAL VERSION ===\n";
   // cout << "\n=== OpenMP VERSION ===\n";

    auto start = chrono::high_resolution_clock::now();
    runAllTests();
    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> elapsed = end - start;
    cout << "\nTotal execution time: " << elapsed.count() << " seconds\n\n";

    return 0;
}

