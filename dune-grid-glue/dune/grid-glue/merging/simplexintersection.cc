// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

namespace Dune {
namespace GridGlue {

template <int dimworld, typename T>
inline void simplexSubdivision(std::integral_constant<int,0>,const std::vector<Dune::FieldVector<T, dimworld> >& elementCorners,
                                                                     std::vector<std::vector<unsigned int> >& subElements,
                                                                     std::vector<std::vector<int> >& faceIds);
template <int dimworld, typename T>
inline void simplexSubdivision(std::integral_constant<int,1>,const std::vector<Dune::FieldVector<T, dimworld> >& elementCorners,
                                                                     std::vector<std::vector<unsigned int> >& subElements,
                                                                     std::vector<std::vector<int> >& faceIds);
template <int dimworld, typename T>
inline void simplexSubdivision(std::integral_constant<int,2>,const std::vector<Dune::FieldVector<T, dimworld> >& elementCorners,
                                                                     std::vector<std::vector<unsigned int> >& subElements,
                                                                     std::vector<std::vector<int> >& faceIds);
template <int dimworld, typename T>
inline void simplexSubdivision(std::integral_constant<int,3>,const std::vector<Dune::FieldVector<T, dimworld> >& elementCorners,
                                                                     std::vector<std::vector<unsigned int> >& subElements,
                                                                     std::vector<std::vector<int> >& faceIds);



// *****************SIMPLEX INTERSECTION COMPUTATION METHODS ***************************
template<int dimWorld,int dim1, int dim2, typename T>
class SimplexMethod : public ComputationMethod<dimWorld,dim1,dim2,T>{
    static_assert(dim1 > dim2, "Specialization missing");
    friend class ComputationMethod<dimWorld,dim1,dim2,T>;
public:
    typedef FieldVector<T, dimWorld> Vector;
    static const int grid1Dimension = dim1;
    static const int grid2Dimension = dim2;
    static const int intersectionDimension = dim2;

    static bool computeIntersectionPoints(const std::vector<FieldVector<T,dimWorld> >&   X,
                 const std::vector<FieldVector<T,dimWorld> >&   Y,
                 std::vector<std::vector<int> >         & SX,
                 std::vector<std::vector<int> >         & SY,
                 std::vector<FieldVector<T,dimWorld> > & P)
    {
        return SimplexMethod<dimWorld,dim2,dim1,T>::computeIntersectionPoints(Y, X, SY, SX, P);
    }

    static void grid1_subdivisions(const std::vector<Vector>& elementCorners,
                                   std::vector<std::vector<unsigned int> >& subElements,
                                   std::vector<std::vector<int> >& faceIds)
    {
        simplexSubdivision<dimWorld,T>(std::integral_constant<int,grid1Dimension>(),
                                       elementCorners,subElements, faceIds);
    }

    static void grid2_subdivisions(const std::vector<Vector>& elementCorners,
                                   std::vector<std::vector<unsigned int> >& subElements,
                                   std::vector<std::vector<int> >& faceIds)
    {
        simplexSubdivision<dimWorld,T>(std::integral_constant<int,grid2Dimension>(),
                                       elementCorners, subElements, faceIds);
    }
};



// POINTS ARE EQUAL
template<int dimWorld,typename T>
class SimplexMethod<dimWorld,0,0,T> : public ComputationMethod<dimWorld,0,0,T>{
    friend class ComputationMethod<dimWorld,0,0,T>;

public:
    typedef FieldVector<T, dimWorld> Vector;
    static const int grid1Dimension = 0;
    static const int grid2Dimension = 0;
    static const int intersectionDimension = 0;

    static bool computeIntersectionPoints(
            const std::vector<FieldVector<T,dimWorld> >&   X,
            const std::vector<FieldVector<T,dimWorld> >&   Y,
            std::vector<std::vector<int> >         & SX,
            std::vector<std::vector<int> >         & SY,
            std::vector<FieldVector<T,dimWorld> > & P)
    {
        assert(X.size() == 1 && Y.size() == 1);

        P.clear(); SX.clear(); SY.clear();
        int k;

        T eps = 1e-8;
        T a = X[0].infinity_norm();
        T b =  Y[0].infinity_norm();
        T c = (X[0] - Y[0]).infinity_norm();

        if (c <= eps*a || c <= eps*b ||
                (a<eps && b< eps && c < 0.5*eps) ) {
            k = insertPoint(X[0],P);

            return true;
        }
        return false;
    }

    static void grid1_subdivisions(const std::vector<Vector>& elementCorners,
                                   std::vector<std::vector<unsigned int> >& subElements,
                                   std::vector<std::vector<int> >& faceIds)
    {
        simplexSubdivision<dimWorld,T>(std::integral_constant<int,grid1Dimension>(),
                                       elementCorners,subElements, faceIds);
    }

    static void grid2_subdivisions(const std::vector<Vector>& elementCorners,
                                   std::vector<std::vector<unsigned int> >& subElements,
                                   std::vector<std::vector<int> >& faceIds)
    {
        simplexSubdivision<dimWorld,T>(std::integral_constant<int,grid2Dimension>(),
                                       elementCorners, subElements, faceIds);
    }

};


// POINT ON LINE SEGMENT - :)
template<int dimWorld,typename T>
class SimplexMethod<dimWorld,0,1,T> : public ComputationMethod<dimWorld,0,1,T>{
    friend class ComputationMethod<dimWorld,0,1,T>;

public:
    typedef FieldVector<T, dimWorld> Vector;
    static const int grid1Dimension = 0;
    static const int grid2Dimension = 1;
    static const int intersectionDimension = 0;

    static bool computeIntersectionPoints(const std::vector<FieldVector<T,dimWorld> >&   X,
                                          const std::vector<FieldVector<T,dimWorld> >&   Y,
                                          std::vector<std::vector<int> >         & SX,
                                          std::vector<std::vector<int> >         & SY,
                                          std::vector<FieldVector<T,dimWorld> > & P)
    {
        assert(X.size() == 1 && Y.size() == 2);

        P.clear(); SX.clear(); SY.clear();
        SY.resize(2);

        if (dimWorld == 1) {
            T lowerBound = std::max(X[0][0], std::min(Y[0][0],Y[1][0]));
            T upperBound = std::min(X[0][0], std::max(Y[0][0],Y[1][0]));

            if (lowerBound <= upperBound) {      // Intersection is non-empty
                insertPoint(X[0],P);
                return true;
            }
        } else {

            T eps = 1e-8;

            // check whether the point is on the segment
            FieldVector<T,dimWorld> v0 = X[0] - Y[0];
            FieldVector<T,dimWorld> v1 = X[0] - Y[1];
            FieldVector<T,dimWorld> v2 = Y[1] - Y[0];

            T s = v0.dot(v1);
            T t = v0.two_norm()/v2.two_norm();
            v2*=t;
            v2+=Y[0];
            v2-=X[0];

            if (v2.infinity_norm() < eps && s<=eps && t<=1+eps) {

                int k = insertPoint(X[0],P);
                if (s  < eps && t < eps)
                    SY[0].push_back(k);
                else if (s < eps && t>1-eps)
                    SY[1].push_back(k);
                return true;
            }
        }
        return false;
    }

    static void grid1_subdivisions(const std::vector<Vector>& elementCorners,
                                   std::vector<std::vector<unsigned int> >& subElements,
                                   std::vector<std::vector<int> >& faceIds)
    {
        simplexSubdivision<dimWorld,T>(std::integral_constant<int,grid1Dimension>(),
                                       elementCorners,subElements, faceIds);
    }

    static void grid2_subdivisions(const std::vector<Vector>& elementCorners,
                                   std::vector<std::vector<unsigned int> >& subElements,
                                   std::vector<std::vector<int> >& faceIds)
    {
        simplexSubdivision<dimWorld,T>(std::integral_constant<int,grid2Dimension>(),
                                       elementCorners, subElements, faceIds);
    }

};


// POINT IN TRIANGLE - :)
template<int dimWorld,typename T>
class SimplexMethod<dimWorld,0,2,T> : public ComputationMethod<dimWorld,0,2,T>{
    friend class ComputationMethod<dimWorld,0,2,T>;

public:
    typedef FieldVector<T, dimWorld> Vector;
    static const int grid1Dimension = 0;
    static const int grid2Dimension = 2;
    static const int intersectionDimension = 0;

    static bool computeIntersectionPoints(const std::vector<FieldVector<T,dimWorld> >&   X,
                 const std::vector<FieldVector<T,dimWorld> >&   Y,
                 std::vector<std::vector<int> >         & SX,
                 std::vector<std::vector<int> >         & SY,
                 std::vector<FieldVector<T,dimWorld> > & P)
    {
        assert(X.size() == 1 && Y.size() == 3 && dimWorld > 1);

        P.clear(); SX.clear(); SY.clear();
        SY.resize(3);
        int k;

        // If not, check whether it is inside the triangle
        double eps= 1e-8 ;     // tolerance for relative error

        FieldVector<T,dimWorld> v0,v1,v2,r;

        v0 = Y[1] - Y[0];
        v1 = Y[2] - Y[0];
        v2 = X[0] - Y[0];

        T s,t,d;

        d = ((v0.dot(v0))*(v1.dot(v1)) - (v0.dot(v1))*(v0.dot(v1)));

        s = ((v1.dot(v1))*(v0.dot(v2)) - (v0.dot(v1))*(v1.dot(v2))) / d;
        t = ((v0.dot(v0))*(v1.dot(v2)) - (v0.dot(v1))*(v0.dot(v2))) / d;

        v0*=s;
        v1*=t;
        r = Y[0] + v0 + v1;

        if (s > -eps && t > -eps && (s+t)< 1+eps && (r-X[0]).infinity_norm() < eps) {
            k = insertPoint(X[0],P);

            if (t < eps)  { // t ~ 0, s varies -> edge 0
                SY[0].push_back(k);
            }
            if (s < eps)  { // s ~ 0, t varies -> edge 1
                SY[1].push_back(k);
            }
            if (s+t > 1-eps) { // s+t ~ 1 -> edge 2
                SY[2].push_back(k);
            }

            return true;
        }

        return false;
    }

    static void grid1_subdivisions(const std::vector<Vector>& elementCorners,
                                   std::vector<std::vector<unsigned int> >& subElements,
                                   std::vector<std::vector<int> >& faceIds)
    {
        simplexSubdivision<dimWorld,T>(std::integral_constant<int,grid1Dimension>(),
                                       elementCorners,subElements, faceIds);
    }

    static void grid2_subdivisions(const std::vector<Vector>& elementCorners,
                                   std::vector<std::vector<unsigned int> >& subElements,
                                   std::vector<std::vector<int> >& faceIds)
    {
        simplexSubdivision<dimWorld,T>(std::integral_constant<int,grid2Dimension>(),
                                       elementCorners, subElements, faceIds);
    }
};


// POINT IN TETRAHEDRON - : )
template<int dimWorld,typename T>
class SimplexMethod<dimWorld,0,3,T> : public ComputationMethod<dimWorld,0,3,T>{
    friend class ComputationMethod<dimWorld,0,3,T>;

public:
    typedef FieldVector<T, dimWorld> Vector;
    static const int grid1Dimension = 0;
    static const int grid2Dimension = 3;
    static const int intersectionDimension = 0;

    static bool computeIntersectionPoints(const std::vector<FieldVector<T,dimWorld> >&   X,
                 const std::vector<FieldVector<T,dimWorld> >&   Y,
                 std::vector<std::vector<int> >         & SX,
                 std::vector<std::vector<int> >         & SY,
                 std::vector<FieldVector<T,dimWorld> > & P)
    {
        assert(X.size() == 1 && Y.size() == 4 && dimWorld == 3);

        P.clear(); SX.clear(); SY.clear();
        SY.resize(4);

        T eps = 1e-8;
        // if not, check whether its inside the tetrahedron
        FieldMatrix<T,dimWorld+1,dimWorld+1>  D,DD ;

        D[0][0] =  Y[0][0] ;  D[0][1] =  Y[1][0] ;  D[0][2] =  Y[2][0] ;  D[0][3] =  Y[3][0] ;
        D[1][0] =  Y[0][1] ;  D[1][1] =  Y[1][1] ;  D[1][2] =  Y[2][1] ;  D[1][3] =  Y[3][1] ;
        D[2][0] =  Y[0][2] ;  D[2][1] =  Y[1][2] ;  D[2][2] =  Y[2][2] ;  D[2][3] =  Y[3][2] ;
        D[3][0] =        1 ;  D[3][1] =        1 ;  D[3][2] =        1 ;  D[3][3] =        1 ;

        std::array<T, 5> detD;
        detD[0] = D.determinant();

        for(unsigned i = 1; i < detD.size(); ++i) {
            DD = D;
            for (unsigned d = 0; d < dimWorld; ++d)
                DD[d][i-1] = X[0][d];
            detD[i] = DD.determinant();
            if (std::abs(detD[i]) > eps && std::signbit(detD[0]) != std::signbit(detD[i]))
                return false; // We are outside.
        }

        int k = insertPoint(X[0],P);
        unsigned int faces[4] = {3,2,1,0};
        for (unsigned i = 1; i < detD.size(); ++i)
            if(std::abs(detD[i]) < eps)
                SY[faces[i-1]].push_back(k);   // on triangle not containing node i-1

        return true;
    }

    static void grid1_subdivisions(const std::vector<Vector>& elementCorners,
                                   std::vector<std::vector<unsigned int> >& subElements,
                                   std::vector<std::vector<int> >& faceIds)
    {
        simplexSubdivision<dimWorld,T>(std::integral_constant<int,grid1Dimension>(),
                                       elementCorners,subElements, faceIds);
    }

    static void grid2_subdivisions(const std::vector<Vector>& elementCorners,
                                   std::vector<std::vector<unsigned int> >& subElements,
                                   std::vector<std::vector<int> >& faceIds)
    {
        simplexSubdivision<dimWorld,T>(std::integral_constant<int,grid2Dimension>(),
                                       elementCorners, subElements, faceIds);
    }
};

// SEGEMENT-SEGMENT INTERSECTION - :)
template<int dimWorld,typename T>
class SimplexMethod<dimWorld,1,1,T> : public ComputationMethod<dimWorld,1,1,T>{
    friend class ComputationMethod<dimWorld,1,1,T>;

public:
    typedef FieldVector<T, dimWorld> Vector;
    static const int grid1Dimension = 1;
    static const int grid2Dimension = 1;
    static const int intersectionDimension = 1;

    static bool computeIntersectionPoints(const std::vector<FieldVector<T,dimWorld> >&   X,
                 const std::vector<FieldVector<T,dimWorld> >&   Y,
                 std::vector<std::vector<int> >         & SX,
                 std::vector<std::vector<int> >         & SY,
                 std::vector<FieldVector<T,dimWorld> > & P)
    {
        assert(X.size() == 2 && Y.size() == 2);

        P.clear(); SX.clear(); SY.clear();
        SX.resize(2);
        SY.resize(2);
        T eps = 1e-8;
        int k,idX_min=-1,idX_max=-1,idY_min=-1,idY_max=-1;

        // compute intersections
        switch (dimWorld) {
            case 1:  // d
            {
                FieldVector<T,dimWorld> lowerbound(std::max(std::min(X[0][0],X[1][0]), std::min(Y[0][0],X[1][0])));
                FieldVector<T,dimWorld> upperbound(std::min(std::max(X[0][0],X[1][0]), std::max(Y[0][0],Y[1][0])));

                if (lowerbound[0] < upperbound[0]) {      // Intersection is non-empty

                    idX_min = (std::min(X[0][0],X[1][0]) < std::min(Y[0][0],Y[1][0]))?(-1):((X[0][0]<X[1][0])?(0):(1));
                    if (idX_min < 0)
                        idY_min = ((Y[0][0]<Y[1][0])?(0):(1));

                    idX_max = (std::max(X[0][0],X[1][0]) > std::max(Y[0][0],Y[1][0]))?(-1):((X[0][0]>X[1][0])?(0):(1));
                    if (idX_max < 0)
                        idY_max = ((Y[0][0]>Y[1][0])?(0):(1));

                    k = insertPoint(lowerbound,P);
                    if (idX_min >= 0)
                        SX[idX_min].push_back(k);
                    else
                        SY[idY_min].push_back(k);

                    k = insertPoint(upperbound,P);
                    if (idX_max >= 0)
                        SX[idX_max].push_back(k);
                    else
                        SY[idY_max].push_back(k);

                    return true;
                }
                return false;
            }
            case 2: // solve X0 + r_0 * (X1 - X0) = Y0 + r_1 * (Y1 - Y0)
            {    // get size_type for all the vectors we are using

                FieldMatrix<T,dimWorld, dimWorld> A;
                A[0][0] =  X[1][0] - X[0][0]; A[0][1] =  Y[0][0] - Y[1][0];
                A[1][0] =  X[1][1] - X[0][1]; A[1][1] =  Y[0][1] - Y[1][1];

                if (std::abs(A.determinant())>eps) {
                    // lines are non parallel and not degenerated
                    FieldVector<T,dimWorld>  p,r,b = Y[0] - X[0];
                    A.solve(r,b) ;

                    if ((r[0]>-eps)&&(r[0]<=1+eps)&&(r[1]>-eps)&&(r[1]<1+eps)) {
                        p = X[1] - X[0];
                        p *= r[0] ;
                        p += X[0] ;
                        k = insertPoint(p,P);
                        if(r[0] < eps) {            // X = X_0 + r_0 (X_1 - X_0) = X_0
                            SX[0].push_back(k);
                            P[k] = X[0];
                        }
                        else if(r[0] > 1-eps) {     // X = X_0 + r_0 (X_1 - X_0) = X_1
                            SX[1].push_back(k);
                            P[k] = X[1];
                        }
                        if(r[1] < eps){             // Y = Y_0 + r_1 (Y_1 - Y_0) = Y_0
                            SY[0].push_back(k);
                            P[k] = Y[0];
                        }
                        else if(r[1] > 1-eps) {     // Y = Y_0 + r_1 (Y_1 - Y_0) = Y_1
                            SY[1].push_back(k);
                            P[k] = Y[1];
                        }
                        return true;
                    }
                } else if ((X[1]-X[0]).infinity_norm() > eps && (Y[1]-Y[0]).infinity_norm() > eps) {
                    // lines are paralles, but non degenerated
                    bool found = false;

                    // use triangle equality ||a - b||_2 = || a -c ||_2 + || c - b ||_2 for non perpendicular lines
                    for (unsigned i = 0; i < 2; ++i) {
                        if (std::abs((Y[i]-X[0]).two_norm() + std::abs((Y[i]-X[1]).two_norm())
                                     - std::abs((X[1]-X[0]).two_norm())) < eps) {
                            k = insertPoint(Y[i],P);
                            SY[i].push_back(k);
                            found = true;
                        }
                        if (std::abs((X[i]-Y[0]).two_norm() + std::abs((X[i]-Y[1]).two_norm())
                                     - std::abs((Y[1]-Y[0]).two_norm())) < eps) {
                            k = insertPoint(X[i],P);
                            SX[i].push_back(k);
                            found = true;
                        }
                    }
                    return found;
                }
                return false;
            }
            case 3: // solve X0 + r_0 * (X1 - X0) = Y0 + r_1 * (Y1 - Y0)
            {    FieldVector<T,dimWorld>  dX, dY, dZ, cXY, cYZ;

                dX = X[1]-X[0];
                dY = Y[1]-Y[0];
                dZ = Y[0]-X[0];

                cXY[0] = dX[1]* dY[2] - dX[2]* dY[1];
                cXY[1] = dX[2]* dY[0] - dX[0]* dY[2];
                cXY[2] = dX[0]* dY[1] - dX[1]* dY[0];

                if (fabs(dZ.dot(cXY)) < eps*1e+3 && cXY.infinity_norm()>eps) { // coplanar, but not aligned

                    cYZ[0] = dY[1]* dZ[2] - dY[2]* dZ[1];
                    cYZ[1] = dY[2]* dZ[0] - dY[0]* dZ[2];
                    cYZ[2] = dY[0]* dZ[1] - dY[1]* dZ[0];

                    T s = -cYZ.dot(cXY) / cXY.two_norm2();

                    if (s > -eps && s < 1+eps) {
                        dX*= s;
                        dX+= X[0];
                        T o = (dX - Y[0]).two_norm() + (dX- Y[1]).two_norm();

                        if (std::abs(o-dY.two_norm()) < eps) {
                            k = insertPoint(dX,P);
    \
                            if (s<eps) {
                                P[k] = X[0];
                                SX[0].push_back(k);
                            } else if(s > 1-eps) {
                                P[k] = X[1];
                                SX[1].push_back(k);
                            } else if ((dX - Y[0]).two_norm() < eps) {
                                P[k] = Y[0];
                                SY[0].push_back(k);
                            } else if((dX - Y[1]).two_norm() < eps) {
                                P[k] = Y[1];
                                SY[1].push_back(k);
                            }

                            return true;
                        }
                    }
                } else if (cXY.infinity_norm() <= eps) {// lines are parallel

                    bool found = false;

                    // use triangle equality ||a - b||_2 = || a -c ||_2 + || c - b ||_2,
                    // under the assumption (a-c)*(c-b) > 0 or (a-c) = 0 or (c-b) = 0
                    for (unsigned i = 0; i < 2; ++i) {
                        if ((std::abs((Y[i]-X[0]).two_norm() + std::abs((Y[i]-X[1]).two_norm()) // triangle equality
                                     - std::abs((X[1]-X[0]).two_norm())) < eps) &&
                                (std::abs((Y[i]-X[0]).dot((Y[i]-X[1]))) > eps                   // assumption
                                 || (Y[i]-X[0]).infinity_norm() < eps || (Y[i]-X[1]).infinity_norm() < eps)) {
                            k = insertPoint(Y[i],P);
                            SY[i].push_back(k);
                            found = true;
                        }
                        if (std::abs((X[i]-Y[0]).two_norm() + std::abs((X[i]-Y[1]).two_norm())  // triangle equality
                                     - std::abs((Y[1]-Y[0]).two_norm())) < eps &&
                                (std::abs((X[i]-Y[0]).dot((X[i]-Y[1]))) > eps                   // assumption
                                    || (X[i]-Y[0]).infinity_norm() < eps || (X[i]-Y[1]).infinity_norm() < eps)){
                            k = insertPoint(X[i],P);
                            SX[i].push_back(k);
                            found = true;
                        }
                    }
                    return found;
                }

                return false;
            }
        }
        return false;
    }

    static void grid1_subdivisions(const std::vector<Vector>& elementCorners,
                                   std::vector<std::vector<unsigned int> >& subElements,
                                   std::vector<std::vector<int> >& faceIds)
    {
        simplexSubdivision<dimWorld,T>(std::integral_constant<int,grid1Dimension>(),
                                       elementCorners,subElements, faceIds);
    }

    static void grid2_subdivisions(const std::vector<Vector>& elementCorners,
                                   std::vector<std::vector<unsigned int> >& subElements,
                                   std::vector<std::vector<int> >& faceIds)
    {
        simplexSubdivision<dimWorld,T>(std::integral_constant<int,grid2Dimension>(),
                                       elementCorners, subElements, faceIds);
    }
};

// SEGEMENT-TRIANGLE INTERSECTION - : )
template<int dimWorld,typename T>
class SimplexMethod<dimWorld,1,2,T> : public ComputationMethod<dimWorld,1,2,T>{
    friend class ComputationMethod<dimWorld,1,2,T>;

public:
    typedef FieldVector<T, dimWorld> Vector;
    static const int grid1Dimension = 1;
    static const int grid2Dimension = 2;
    static const int intersectionDimension = 1;

    static bool computeIntersectionPoints(const std::vector<FieldVector<T,dimWorld> >&   X,
                 const std::vector<FieldVector<T,dimWorld> >&   Y,
                 std::vector<std::vector<int> >         & SX,
                 std::vector<std::vector<int> >         & SY,
                 std::vector<FieldVector<T,dimWorld> > & P)
    {

        assert(X.size() == 2 && Y.size() == 3 && dimWorld > 1);
        P.clear(); SX.clear(); SY.clear();
        SX.resize(2);
        SY.resize(3);

        int k;
        std::vector<FieldVector<T,dimWorld> > surfPts, edge(2), pni(1);
        std::vector<std::vector<int> > hSX, hSY;

        // is any segment point inside the triangle?
        for (unsigned ni = 0; ni < 2; ++ni) {
            pni[0] = X[ni];

            if (SimplexMethod<dimWorld,0,2,T>::computeIntersectionPoints(pni,Y,hSX,hSY,surfPts)) {
                    k = insertPoint(X[ni],P);
                    SX[ni].push_back(k);
                    for (unsigned e=0; e < 3; ++e)
                        if (hSY[e].size() > 0)
                            SY[e].push_back(k);

            }
            surfPts.clear(); hSX.clear(); hSY.clear();
        }

        if (P.size() >= 2)  // we cannot have more than two intersection points
            return true;

        unsigned int faces[3] = {0,2,1};
        // do triangle faces intersect with the segment?
        for (unsigned ni = 0; ni < 3; ++ni) {
            edge[0] = Y[ni];
            edge[1] = Y[(ni+1)%3];

            if (SimplexMethod<dimWorld,1,1,T>::computeIntersectionPoints(X,edge,hSX,hSY,surfPts)) {
                for (unsigned ne = 0; ne < surfPts.size(); ++ne) {
                    k = insertPoint(surfPts[ne],P);
                    SY[faces[ni]].push_back(k);
                    if (hSX[0].size() > 0)
                        SX[0].push_back(k);
                    if (hSX[1].size() > 0)
                        SX[1].push_back(k);
                }

                if (P.size() >= 2)  // we cannot have more than two intersection points
                    return true;

                surfPts.clear(); hSX.clear(); hSY.clear();
            }
        }

        if (P.size() >= 2)  // we cannot have more than two intersection points
            return true;

        // if line and triangle are not coplanar in 3d and do not intersect at boundaries
        if (dimWorld == 3) {
            T eps = 1e-8;

            Dune::FieldVector<T,dimWorld>      B,r,p ;
            Dune::FieldMatrix<T,dimWorld,dimWorld>    A ;

            B = Y[0] - X[0] ;

            for (unsigned i = 0; i < dimWorld; ++i) {
                A[i][0] = (X[1][i] - X[0][i]);
                A[i][1] = (Y[0][i] - Y[1][i]);
                A[i][dimWorld-1] = (Y[0][i] - Y[dimWorld-1][i]);
            }

            if (std::abs(A.determinant())>eps) {

                A.solve(r,B) ;

                if ((r[0]>=-eps)&&(r[0]<=1+eps)
                        &&(r[1]>=-eps)&&(r[1]<=1+eps)
                        &&(r[2]>=-eps)&&(r[2]<=1+eps)
                        &&(r[1]+r[2]>=-eps) &&(r[1]+r[2]<=1+eps)) {
                    p =  X[1] - X[0] ;
                    p *= r[0] ;
                    p += X[0] ;
                    k = insertPoint(p,P);

                    if (std::abs(r[0]) < eps) // we prefer exact locations
                        P[k] = X[0];
                    else if (std::abs(r[0]) > 1-eps)
                        P[k] = X[1];
                    else if (std::abs(r[1]) < eps && std::abs(r[2]) < eps)
                        P[k] = Y[0];
                    else if (std::abs(r[1]) < eps && std::abs(r[2]) > 1-eps)
                        P[k] = Y[2];
                    else if (std::abs(r[1]) > 1-eps && std::abs(r[2]) < eps)
                        P[k] = Y[1];

                    if (std::abs(r[1])  < eps)
                        SY[1].push_back(k);
                    if (std::fabs(r[dimWorld-1])  < eps)
                        SY[0].push_back(k);
                    if (std::fabs(r[1]+r[dimWorld-1] - 1)  < eps)
                        SY[2].push_back(k);


                    return true ;
                }
            }
        }
        return false ;
    }

    static void grid1_subdivisions(const std::vector<Vector>& elementCorners,
                                   std::vector<std::vector<unsigned int> >& subElements,
                                   std::vector<std::vector<int> >& faceIds)
    {
        simplexSubdivision<dimWorld,T>(std::integral_constant<int,grid1Dimension>(),
                                       elementCorners,subElements, faceIds);
    }

    static void grid2_subdivisions(const std::vector<Vector>& elementCorners,
                                   std::vector<std::vector<unsigned int> >& subElements,
                                   std::vector<std::vector<int> >& faceIds)
    {
        simplexSubdivision<dimWorld,T>(std::integral_constant<int,grid2Dimension>(),
                                       elementCorners, subElements, faceIds);
    }
};

// SEGEMENT-TETRAHEDRON INTERSECTION - : )
template<int dimWorld,typename T>
class SimplexMethod<dimWorld,1,3,T> : public ComputationMethod<dimWorld,1,3,T>{
    friend class ComputationMethod<dimWorld,1,3,T>;

public:
    typedef FieldVector<T, dimWorld> Vector;
    static const int grid1Dimension = 1;
    static const int grid2Dimension = 3;
    static const int intersectionDimension = 1;

    static bool computeIntersectionPoints(const std::vector<FieldVector<T,dimWorld> >&   X,
                 const std::vector<FieldVector<T,dimWorld> >&   Y,
                 std::vector<std::vector<int> >         & SX,
                 std::vector<std::vector<int> >         & SY,
                 std::vector<FieldVector<T,dimWorld> > & P)
    {
        assert(X.size() == 2 && Y.size() == 4 && dimWorld > 2);

        std::vector<int> indices;
        std::vector<FieldVector<T,dimWorld> > surfPts;
        std::vector<std::vector<int> > hSX, hSY;

        P.clear(); SX.clear(); SY.clear();
        SX.resize(2);
        SY.resize(4);

        std::vector<FieldVector<T,dimWorld> > pni(1);

        bool b = false;

        // check whether the corners of the segment are contained in the tetrahedron
        for (unsigned int ci = 0; ci < 2; ++ ci) {

            pni[0] = X[ci];
            if(SimplexMethod<dimWorld,0,3,T>::computeIntersectionPoints(pni,Y,hSX,hSY,surfPts)) {
                int k = insertPoint(X[ci],P);
                SX[ci].push_back(k);
                b=true;
            }
            surfPts.clear();
            hSX.clear(); hSY.clear();
        }

        if (P.size() == 2)
            return true;

        unsigned int faces[4] = {0,3,2,1};
        // check whether tetrahedron faces intersect with segment
        std::vector<FieldVector<T,dimWorld> > triangle(3);
        for (unsigned int ci = 0; ci < 4; ++ci) { // iterate over all faces
            triangle[0] = Y[ci];
            triangle[1] = Y[(ci+1)%4];
            triangle[2] = Y[(ci+2)%4];

            if (SimplexMethod<dimWorld,1,2,T>::computeIntersectionPoints(X,triangle,hSX,hSY,surfPts)) { // seg - triangle intersection
                std::vector<int> indices(surfPts.size());
                for (unsigned int np = 0; np < surfPts.size(); ++np) {
                    int k = insertPoint(surfPts[np],P);
                    indices[np]=k;
                    SY[faces[ci]].push_back(k);
                }

                // hSX[*] is nonempty if the intersection point is on an edge of the current face of Y
                for (unsigned int np = 0; np < hSX[0].size(); ++np)
                    SX[0].push_back(indices[hSX[0][np]]);
                for (unsigned int np = 0; np < hSX[1].size(); ++np)
                    SX[0].push_back(indices[hSX[1][np]]);

                b = true;
            }
            surfPts.clear();
            hSX.clear(); hSY.clear();
        }

        return b;
    }

    static void grid1_subdivisions(const std::vector<Vector>& elementCorners,
                                   std::vector<std::vector<unsigned int> >& subElements,
                                   std::vector<std::vector<int> >& faceIds)
    {
        simplexSubdivision<dimWorld,T>(std::integral_constant<int,grid1Dimension>(),
                                       elementCorners,subElements, faceIds);
    }

    static void grid2_subdivisions(const std::vector<Vector>& elementCorners,
                                   std::vector<std::vector<unsigned int> >& subElements,
                                   std::vector<std::vector<int> >& faceIds)
    {
        simplexSubdivision<dimWorld,T>(std::integral_constant<int,grid2Dimension>(),
                                       elementCorners, subElements, faceIds);
    }
};

// TRIANGLE -TRIANGLE INTERSECTION
template<int dimWorld,typename T>
class SimplexMethod<dimWorld,2,2,T> : public ComputationMethod<dimWorld,2,2,T>{
    friend class ComputationMethod<dimWorld,2,2,T>;

public:
    typedef FieldVector<T, dimWorld> Vector;
    static const int grid1Dimension = 2;
    static const int grid2Dimension = 2;
    static const int intersectionDimension = 2;

    static bool computeIntersectionPoints(const std::vector<FieldVector<T,dimWorld> >&   X,
                 const std::vector<FieldVector<T,dimWorld> >&   Y,
                 std::vector<std::vector<int> >         & SX,
                 std::vector<std::vector<int> >         & SY,
                 std::vector<FieldVector<T,dimWorld> > & P)
    {
        assert(X.size() == 3 && Y.size() == 3 && dimWorld > 1);

        P.clear(); SX.clear(); SY.clear();

        SX.resize(3);
        SY.resize(3);

        bool b = false;

        std::vector<FieldVector<T,dimWorld> > edge(2);
        std::vector<FieldVector<T,dimWorld> > surfPts;
        std::vector<std::vector<int> > hSX, hSY;
        std::vector<int> indices;

        unsigned int faces[3] = {0,2,1};

        for (unsigned int ni = 0; ni < 3; ++ni) {
            // check whether the faces of triangle Y intersect the triangle X
            edge[0] = Y[ni];
            edge[1] = Y[(ni+1)%3];

            if(SimplexMethod<dimWorld,2,1,T>::computeIntersectionPoints(X,edge,hSX,hSY,surfPts)) {

                indices.resize(surfPts.size());
                // add intersections of edges of Y with triangle X
                for (unsigned int np = 0; np < surfPts.size(); ++np) {
                    int k = insertPoint(surfPts[np],P);
                    indices[np] = k;
                    SY[faces[ni]].push_back(k); // add edge data
                }

                b=true;
            }
            if (P.size() >= 6)
                return true;

            surfPts.clear(); hSX.clear(); hSY.clear();
            // check whether the faces of triangle X intersect the triangle Y
            edge[0] = X[ni];
            edge[1] = X[(ni+1)%3];

            if(SimplexMethod<dimWorld,1,2,T>::computeIntersectionPoints(edge,Y,hSX,hSY,surfPts)) {

                indices.resize(surfPts.size());
                // add intersections of edges of X with triangle Y
                for (unsigned int np = 0; np < surfPts.size(); ++np) {
                    int k = insertPoint(surfPts[np],P);
                    indices[np] = k;
                    SX[faces[ni]].push_back(k); // add edge data
                }

                b=true;
            }
            if (P.size() >= 6)
                return true;

            surfPts.clear(); hSX.clear(); hSY.clear();
        }

        return b;
    }

    static void grid1_subdivisions(const std::vector<Vector>& elementCorners,
                                   std::vector<std::vector<unsigned int> >& subElements,
                                   std::vector<std::vector<int> >& faceIds)
    {
        simplexSubdivision<dimWorld,T>(std::integral_constant<int,grid1Dimension>(),
                                       elementCorners,subElements, faceIds);
    }

    static void grid2_subdivisions(const std::vector<Vector>& elementCorners,
                                   std::vector<std::vector<unsigned int> >& subElements,
                                   std::vector<std::vector<int> >& faceIds)
    {
        simplexSubdivision<dimWorld,T>(std::integral_constant<int,grid2Dimension>(),
                                       elementCorners, subElements, faceIds);
    }
};

// TRIANGLE -TETRAHEDRON INTERSECTION -: )
template<int dimWorld,typename T>
class SimplexMethod<dimWorld,2,3,T> : public ComputationMethod<dimWorld,2,3,T>{
    friend class ComputationMethod<dimWorld,2,3,T>;

public:
    typedef FieldVector<T, dimWorld> Vector;
    static const int grid1Dimension = 2;
    static const int grid2Dimension = 3;
    static const int intersectionDimension = 2;

    static bool computeIntersectionPoints(const std::vector<FieldVector<T,dimWorld> >&   X,
                 const std::vector<FieldVector<T,dimWorld> >&   Y,
                 std::vector<std::vector<int> >         & SX,
                 std::vector<std::vector<int> >         & SY,
                 std::vector<FieldVector<T,dimWorld> > & P)
    {
        assert(X.size() == 3 && Y.size() == 4 && dimWorld > 2);
        P.clear(); SX.clear(); SY.clear();
        SX.resize(3);
        SY.resize(4);

        bool b = false;
        int k,ni,np, ep;
        std::vector<FieldVector<T,dimWorld> > surfPts, xni(1);
        std::vector<std::vector<int> > hSX, hSY;

        unsigned int fiX[3][2];

        fiX[0][0] = 0; fiX[1][0] = 0; fiX[2][0] = 1; // faces to node
        fiX[0][1] = 1; fiX[1][1] = 2; fiX[2][1] = 2;
        // 1st step: check whether the points of the triangle are contained in the tetrahedron
        for (ni = 0; ni < 3; ++ni) {

            xni[0] = X[ni];
            if (SimplexMethod<dimWorld,0,3,T>::computeIntersectionPoints(xni,Y,hSX,hSY,surfPts)) {
                std::vector<int> indices(surfPts.size());
                for (np = 0; np < surfPts.size(); ++np) {
                    k = insertPoint(X[ni],P);
                    indices[np] = k;
                    SX[fiX[ni][0]].push_back(k); // the corresponding edges to the node are ni and (ni+2)%3
                    SX[fiX[ni][1]].push_back(k);
                }

                for (ep = 0;  ep < 4; ++ep) {
                    for (np = 0; np < hSY[ep].size();++np) {
                        SY[ep].push_back(indices[hSY[ep][np]]);
                    }
                }
                b = true;
            }
            hSX.clear(); hSY.clear(); surfPts.clear();
        }

        if (P.size() == 3) // intersection is given by all three corners of the triangle
            return true;

        // note: points of Y in X is called indirectly via triangle-triangle intesection

        // check whether the triangles of the one tetrahedron intersect the triangle
        unsigned int facesY[4] = {0,3,2,1}; // face ordering

        std::vector<FieldVector<T,dimWorld> > triangle(3);
        for (ni = 0; ni < 4; ++ni) {

            triangle[0] = Y[ni];
            triangle[1] = Y[(ni+1)%4];
            triangle[2] = Y[(ni+2)%4];

            if (SimplexMethod<dimWorld,2,2,T>::computeIntersectionPoints(X,triangle,hSX,hSY,surfPts)) {
                std::vector<int> indices(surfPts.size());
                for (np = 0; np < surfPts.size(); ++np) {
                    k = insertPoint(surfPts[np],P);
                    indices[np]=k;
                    SY[facesY[ni]].push_back(k);
                }

                // SX[*] is nonempty if the face * of X is intersected
                for (np = 0; np < 3; ++np) {
                    for (ep = 0; ep < hSX[np].size(); ++ep) {
                        SX[np].push_back(indices[hSX[np][ep]]);
                    }
                }

                b = true;
            }
            hSX.clear(); hSY.clear(); surfPts.clear();
        }
        return b;
    }

    static void grid1_subdivisions(const std::vector<Vector>& elementCorners,
                                   std::vector<std::vector<unsigned int> >& subElements,
                                   std::vector<std::vector<int> >& faceIds)
    {
        simplexSubdivision<dimWorld,T>(std::integral_constant<int,grid1Dimension>(),
                                       elementCorners,subElements, faceIds);
    }

    static void grid2_subdivisions(const std::vector<Vector>& elementCorners,
                                   std::vector<std::vector<unsigned int> >& subElements,
                                   std::vector<std::vector<int> >& faceIds)
    {
        simplexSubdivision<dimWorld,T>(std::integral_constant<int,grid2Dimension>(),
                                       elementCorners, subElements, faceIds);
    }
};

template<int dimWorld,typename T>
class SimplexMethod<dimWorld,3,3,T> : public ComputationMethod<dimWorld,3,3,T>{
    friend class ComputationMethod<dimWorld,3,3,T>;

public:
    typedef FieldVector<T, dimWorld> Vector;
    static const int grid1Dimension = 3;
    static const int grid2Dimension = 3;
    static const int intersectionDimension = 3;

    static bool computeIntersectionPoints(const std::vector<FieldVector<T,dimWorld> >&   X,
                 const std::vector<FieldVector<T,dimWorld> >&   Y,
                 std::vector<std::vector<int> >         & SX,
                 std::vector<std::vector<int> >         & SY,
                 std::vector<FieldVector<T,dimWorld> > & P)
    {
        assert(X.size() == 4 && Y.size() == 4 && dimWorld > 2);
        P.clear(); SX.clear(); SY.clear();

        SX.resize(4);
        SY.resize(4);

        bool b = false;
        int ci,np,ne,k,ni[4][3];

        ni[0][0]= 0 ;   ni[0][1]= 1 ;   ni[0][2]= 2 ;   // faces touching each node
        ni[1][0]= 0 ;   ni[1][1]= 1 ;   ni[1][2]= 3 ;
        ni[2][0]= 0 ;   ni[2][1]= 2 ;   ni[2][2]= 3 ;
        ni[3][0]= 1 ;   ni[3][1]= 2 ;   ni[3][2]= 3 ;

        // temporal data container
        std::vector<FieldVector<T,dimWorld> > surfPts, pci(1);
        std::vector<std::vector<int> > hSX, hSY;

        // check whether the points of the one tetrahedron are contained in the other tetrahedron
        for (ci = 0; ci < 3; ++ci) {
            pci[0] = X[ci];
            // point of X is inside Y
            if (SimplexMethod<dimWorld,0,3,T>::computeIntersectionPoints(pci,Y,hSX,hSY,surfPts)) {
                std::vector<int> indices(surfPts.size());

                for (np = 0; np < surfPts.size(); ++np) {

                    k = insertPoint(X[ci],P);
                    indices[np]=k;
                    SX[ni[ci][0]].push_back(k); // add face data
                    SX[ni[ci][1]].push_back(k);
                    SX[ni[ci][2]].push_back(k);
                }

                // hSY[*] is nonempty if X[ci] is on the face * of Y
                for (ne = 0; ne < 4; ++ne) {
                    for (np = 0; np < hSY[ne].size(); ++np) {
                        SY[ne].push_back(indices[hSY[ne][np]]);
                    }
                }

                b = true;
            }
            hSX.clear(); hSY.clear(); surfPts.clear();

            // probably one point of Y is inside X
            surfPts.resize(0);
            pci[0]=Y[ci];
            if (SimplexMethod<dimWorld,0,3,T>::computeIntersectionPoints(pci,X,hSY,hSX,surfPts)) {
                std::vector<int> indices(surfPts.size());

                for (np = 0; np < surfPts.size(); ++np) {
                    k = insertPoint(Y[ci],P);
                    indices[np]=k;
                    SY[ni[ci][0]].push_back(k); // add face data
                    SY[ni[ci][1]].push_back(k);
                    SY[ni[ci][2]].push_back(k);
                }

                // hSX[*] is nonempty if the point Y[ci] is on the face * of X
                for (ne = 0; ne < 4; ++ne) {
                    for (np = 0; np < hSX[ne].size(); ++np) {
                        SX[ne].push_back(indices[hSX[ne][np]]);
                    }
                }
                b = true;
            }
            hSX.clear(); hSY.clear(); surfPts.clear();
        }

        // check whether the triangles of the one tetrahedron intersect the triangles
        // of the other tetrahedron
        unsigned int faces[4] = {0,3,2,1}; // face ordering

        std::vector<FieldVector<T,dimWorld> > triangle(3);
        for (ci = 0; ci < 4; ++ci) { // iterate over faces of Y

            triangle[0] = Y[ci];
            triangle[1] = Y[(ci+1)%4];
            triangle[2] = Y[(ci+2)%4];

            if(SimplexMethod<dimWorld,3,2,T>::computeIntersectionPoints(X, triangle,hSX,hSY,surfPts)) {

                // add Triangle of Y intersects tetrahedron Y data
                for (np = 0; np < surfPts.size(); ++np) {
                    k = insertPoint(surfPts[np],P);
                    SY[faces[ci]].push_back(k); // add face data
                }
                b = true;
            }
            hSX.clear(); hSY.clear(); surfPts.clear();

            triangle[0] = X[ci];
            triangle[1] = X[(ci+1)%4];
            triangle[2] = X[(ci+2)%4];

            if(SimplexMethod<dimWorld,3,2,T>::computeIntersectionPoints(Y, triangle,hSY,hSX,surfPts)) {

                // add Triangle of Y intersects tetrahedron Y data
                for (np = 0; np < surfPts.size(); ++np) {
                    k = insertPoint(surfPts[np],P);
                    SX[faces[ci]].push_back(k); // add face data
                }
                b = true;
            }
            hSX.clear(); hSY.clear(); surfPts.clear();
        }

        return b;
    }

    static void grid1_subdivisions(const std::vector<Vector>& elementCorners,
                                   std::vector<std::vector<unsigned int> >& subElements,
                                   std::vector<std::vector<int> >& faceIds)
    {
        simplexSubdivision<dimWorld,T>(std::integral_constant<int,grid1Dimension>(),
                                       elementCorners,subElements, faceIds);
    }

    static void grid2_subdivisions(const std::vector<Vector>& elementCorners,
                                   std::vector<std::vector<unsigned int> >& subElements,
                                   std::vector<std::vector<int> >& faceIds)
    {
        simplexSubdivision<dimWorld,T>(std::integral_constant<int,grid2Dimension>(),
                                       elementCorners, subElements, faceIds);
    }
};

template <int dimworld, typename T>
inline void simplexSubdivision(std::integral_constant<int,0>,
                               const std::vector<Dune::FieldVector<T, dimworld> >& elementCorners,
                               std::vector<std::vector<unsigned int> >& subElements,
                               std::vector<std::vector<int> >& faceIds)
{
    subElements.resize(1);
    faceIds.resize(0);

    subElements[0].push_back(0);
}

template <int dimworld, typename T>
inline void simplexSubdivision(std::integral_constant<int,1>,
                               const std::vector<Dune::FieldVector<T, dimworld> >& elementCorners,
                               std::vector<std::vector<unsigned int> >& subElements,
                               std::vector<std::vector<int> >& faceIds)
{
    subElements.resize(1);
    faceIds.resize(1);

    subElements[0].push_back(0);
    subElements[0].push_back(1);

    faceIds[0].push_back(0);
    faceIds[0].push_back(1);
}

template <int dimworld, typename T>
inline void simplexSubdivision(std::integral_constant<int,2>,
                               const std::vector<Dune::FieldVector<T, dimworld> >& elementCorners,
                               std::vector<std::vector<unsigned int> >& subElements,
                               std::vector<std::vector<int> >& faceIds)
{
    subElements.clear();
    faceIds.clear();

    if (elementCorners.size() == 3) { // triangle
        subElements.resize(1);
        faceIds.resize(1);

        subElements[0].push_back(0);
        subElements[0].push_back(1);
        subElements[0].push_back(2);

        faceIds[0].push_back(0);
        faceIds[0].push_back(1);
        faceIds[0].push_back(2);
    } else if (elementCorners.size() == 4) { // quadrilateral => 2 triangles
        subElements.resize(2);
        faceIds.resize(2);

        subElements[0].push_back(0);
        subElements[0].push_back(1);
        subElements[0].push_back(2);

        subElements[1].push_back(1);
        subElements[1].push_back(2);
        subElements[1].push_back(3);

        faceIds[0].push_back(2);
        faceIds[0].push_back(0);
        faceIds[0].push_back(-1);

        faceIds[1].push_back(-1);
        faceIds[1].push_back(1);
        faceIds[1].push_back(3);
    }
}

template <int dimworld, typename T>
inline void simplexSubdivision(std::integral_constant<int,3>,
                               const std::vector<Dune::FieldVector<T, dimworld> >& elementCorners,
                               std::vector<std::vector<unsigned int> >& subElements,
                               std::vector<std::vector<int> >& faceIds)
{
    subElements.clear();
    faceIds.clear();

    if (elementCorners.size() == 4) { // tetrahedron
        subElements.resize(1);
        faceIds.resize(1);

        subElements[0].push_back(0);
        subElements[0].push_back(1);
        subElements[0].push_back(2);
        subElements[0].push_back(3);

        faceIds[0].push_back(0);
        faceIds[0].push_back(1);
        faceIds[0].push_back(2);
        faceIds[0].push_back(3);

    } else if (elementCorners.size() == 8) { // cube => 5 tetrahedra
        subElements.resize(5);
        faceIds.resize(5);

        subElements[0].push_back(0);
        subElements[0].push_back(2);
        subElements[0].push_back(3);
        subElements[0].push_back(6);

        subElements[1].push_back(0);
        subElements[1].push_back(1);
        subElements[1].push_back(3);
        subElements[1].push_back(5);

        subElements[2].push_back(0);
        subElements[2].push_back(3);
        subElements[2].push_back(5);
        subElements[2].push_back(6);

        subElements[3].push_back(0);
        subElements[3].push_back(4);
        subElements[3].push_back(5);
        subElements[3].push_back(6);

        subElements[4].push_back(3);
        subElements[4].push_back(5);
        subElements[4].push_back(6);
        subElements[4].push_back(7);

        faceIds[0].push_back(4);
        faceIds[0].push_back(0);
        faceIds[0].push_back(-1);
        faceIds[0].push_back(3);

        faceIds[1].push_back(4);
        faceIds[1].push_back(2);
        faceIds[1].push_back(-1);
        faceIds[1].push_back(1);

        faceIds[2].push_back(-1);
        faceIds[2].push_back(-1);
        faceIds[2].push_back(-1);
        faceIds[2].push_back(-1);

        faceIds[3].push_back(2);
        faceIds[3].push_back(0);
        faceIds[3].push_back(-1);
        faceIds[3].push_back(5);

        faceIds[4].push_back(-1);
        faceIds[4].push_back(1);
        faceIds[4].push_back(3);
        faceIds[4].push_back(5);
    }
}

} /* namespace Dune::GridGlue */
} /* namespace Dune */
