// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

namespace Dune {
namespace GridGlue {

//****************************************************************************************
// PUBLIC
//****************************************************************************************

template<class CM>
bool IntersectionComputation<CM>::computeIntersection(const std::vector<V>& X,
                                                      const std::vector<V>& Y,
                                                      std::vector<std::vector<int> >& SX,
                                                      std::vector<std::vector<int> >& SY,
                                                      std::vector<V>& P) {

    std::vector<std::vector<unsigned int> > subElementsX, subElementsY;
    std::vector<std::vector<int> > faceIdsX, faceIdsY;
    std::vector<V> subElementX(CM::grid1Dimension+1), subElementY(CM::grid2Dimension+1), sP;
    std::vector<std::vector<int> > sSX, sSY;

    CM::grid1_subdivisions(X,subElementsX,faceIdsX);
    CM::grid2_subdivisions(Y,subElementsY,faceIdsY);

    bool intersectionFound = false;

    for (unsigned int i = 0; i < subElementsX.size(); ++i) { // iterate over all X subelements
        for (unsigned int ki = 0; ki < subElementsX[i].size(); ++ki) // define the X subelement
            subElementX[ki] = X[subElementsX[i][ki]];
        for (unsigned int j = 0; j < subElementsY.size(); ++j) { // iterate over all Y subelemetns
            for (unsigned int kj = 0; kj < subElementsY[j].size(); ++kj) // define the Y subleement
                subElementY[kj] = Y[subElementsY[j][kj]];

            sP.clear();

            // compute the intersection
            bool b = CM::computeIntersectionPoints(subElementX,subElementY,sSX,sSY,sP);
            intersectionFound = intersectionFound || b;

            // only insert points on outer faces
            for (unsigned int ki = 0; ki < sSX.size(); ++ki) { // iterate over all faces
                if (faceIdsX[i][ki] >= 0) {
                    for (unsigned int kii = 0; kii < sSX[ki].size(); ++kii) {
                        int k = insertPoint(sP[sSX[ki][kii]],P);  // determine index in P
                        SX[faceIdsX[i][ki]].push_back(k);
                    }
                }
            }
            for (unsigned int kj = 0; kj < sSY.size(); ++kj) { // iterate over all faces
                if (faceIdsY[j][kj] >= 0) {
                    for (unsigned int kjj = 0; kjj < sSY[kj].size(); ++kjj) {
                        int k = insertPoint(sP[sSY[kj][kjj]],P);  // determine index in P
                        SY[faceIdsY[j][kj]].push_back(k);
                    }
                }
            }
        }
    }

    return intersectionFound;
}

//****************************************************************************************
// PRIVATE
//****************************************************************************************

template<class CM>
void IntersectionComputation<CM>::orderPoints_(std::integral_constant<int,3>,
                                               std::integral_constant<int,3>,
                                               const V& centroid,
                                               const std::vector<std::vector<int> >& SX,
                                               const std::vector<std::vector<int> >& SY,
                                               const std::vector<V>&  P,
                                               std::vector<std::vector<int> >& H)
{
    int n_facesX = SX.size();
    int n_facesY = SY.size();
    int m;

    std::vector<int> no,id,temp ;
    std::vector<V> p ;
    std::vector<std::vector<int> > tempH;

    std::vector<int> faceOrderingX(n_facesX);
    std::vector<int> faceOrderingY(n_facesY);

    if (n_facesX==3) {
        faceOrderingX[0] = 0; faceOrderingX[1] = 2; faceOrderingX[2] = 1;
    } else {
        faceOrderingX[0] = 0; faceOrderingX[1] = 3; faceOrderingX[2] = 2; faceOrderingX[3] = 1;
    }
    if (n_facesY==3) {
        faceOrderingY[0] = 0; faceOrderingY[1] = 2; faceOrderingY[2] = 1;
    } else {
        faceOrderingY[0] = 0; faceOrderingY[1] = 3; faceOrderingY[2] = 2; faceOrderingY[3] = 1;
    }

    if (P.size() > 3) {
        for (int i = 0; i < n_facesX; ++i) { // loop on faces of X
            if (SX[i].size() > 0) {
                no = SX[faceOrderingX[i]];
                removeDuplicates(no);
                m = no.size() ;
                if ((m>=3) && newFace3D(no,tempH))  // don't compute degenerate polygons and check if face is new
                {
                    for ( int l=0; l<m; ++l)
                        p.push_back(P[no[l]]);
                    orderPointsCC(std::integral_constant<int,3>(), centroid,id,p); // order points counter-clock-wise
                    for ( int l=0; l<m; ++l)
                        temp.push_back(no[id[l]]) ;
                    tempH.push_back(temp) ;
                    temp.clear();
                    p.clear();
                    id.clear();                 // clean
                }
                no.clear() ;             // clean
            }
        }
        for (int i = 0; i < n_facesY; ++i) { // loop on faces of Y
            if (SY[i].size() > 0) {
                no = SY[faceOrderingY[i]];
                removeDuplicates(no);
                m = no.size() ;
                if ((m>=3) && newFace3D(no,tempH))  // don't compute degenerate polygons and check if face is new
                {
                    for ( int l=0; l<m; ++l)
                        p.push_back(P[no[l]]) ;
                    orderPointsCC(std::integral_constant<int,3>(),centroid,id,p); // order points counter-clock-wise
                    for ( int l=0; l<m; ++l)
                        temp.push_back(no[id[l]]) ;
                    tempH.push_back(temp) ;
                    temp.clear();
                    p.clear();
                    id.clear();                 // clean
                }
                no.clear() ;             // clean
            }
        }
    }

    for (int i = 0; i < tempH.size(); ++i) {
        int hs = tempH[i].size();
        if (hs >= 3) {
            for (int j = 1; j <= hs-2;++j) {
                temp.clear();
                temp.push_back(tempH[i][0]);
                for (int k = 0; k < 2; ++k)
                    temp.push_back(tempH[i][j+k]);
                H.push_back(temp);
            }
        }
    }
}

template<class CM>
void IntersectionComputation<CM>::orderPoints_(std::integral_constant<int,2>,
                                               std::integral_constant<int,2>,
                                               const V& centroid,
                                               const std::vector<std::vector<int> >& SX,
                                               const std::vector<std::vector<int> >& SY,
                                               const std::vector<V>&  P,
                                               std::vector<std::vector<int> >& H)
{
    H.clear();
    std::vector<int> id, temp(2);

    orderPointsCC(std::integral_constant<int,2>(),centroid,id,P);

    for (std::size_t i = 0; i < id.size();++i) {
        temp[0] = id[i];
        temp[1] = id[(i+1)%(id.size())];
        H.push_back(temp);
    }
}

template<class CM>
void IntersectionComputation<CM>::orderPoints_(std::integral_constant<int,2>,
                                               std::integral_constant<int,3>,
                                               const V& centroid,
                                               const std::vector<std::vector<int> >& SX,
                                               const std::vector<std::vector<int> >& SY,
                                               const std::vector<V>&  P,
                                               std::vector<std::vector<int> >& H)
{
    H.clear();
    std::vector<int> id, temp(2);

    orderPointsCC(std::integral_constant<int,3>(),centroid,id,P);

    for (int i = 0; i < id.size();++i) {
        temp[0] = id[i];
        temp[1] = id[(i+1)%(id.size())];
        H.push_back(temp);
    }
}

template<class CM>
void IntersectionComputation<CM>::removeDuplicates(std::vector<int> & p)
{
    sort(p.begin(),p.end());
    std::vector<int>::iterator it = std::unique(p.begin(),p.end());
    p.erase(it,p.end());
}

template<class CM>
bool IntersectionComputation<CM>::newFace3D(const std::vector<int>& id,
                                            const std::vector<std::vector<int> >& H)
{
    // get size_type for all the vectors we are using
    typedef typename std::vector<Empty>::size_type size_type;

    int n = H.size() ;
    int m = id.size() ;
    std::vector<int> A ;
    std::vector<int> B = id ;
    sort(B.begin(),B.end()) ;
    int i = 0 ;
    bool b = true ;
    double tp ;

    while ( b && (i<n) )
    {
        if ((H[i].size())>=m)
        {
            A=H[i] ;
            sort(A.begin(),A.end());
            tp = 0 ;
            for ( size_type j=0 ; j < m; j++)
                tp += std::fabs(A[j]-B[j]) ;
            b = (tp>0) ;
        }

        i += 1 ;
    }

    return b ;
}


template<class CM>
void IntersectionComputation<CM>::orderPointsCC(std::integral_constant<int,3>,
                                                const V& centroid,
                                                std::vector<int>& id,
                                                const std::vector<V>& P)
{
    typedef typename std::vector<Empty>::size_type size_type;

    id.clear();

    // get size_type for all the vectors we are using
    V c,d1,d2,dr,dn,cross,d ;
    std::vector<typename V::value_type> ai ;

    d1 = P[1] - P[0] ;    // two reference vectors
    d2 = P[2] - P[0] ;

    cross[0] = d1[1]*d2[2] - d1[2]*d2[1] ;    // cross product
    cross[1] = d1[2]*d2[0] - d1[0]*d2[2] ;
    cross[2] = d1[0]*d2[1] - d1[1]*d2[0] ;

    if (((centroid - P[0])*cross)<0)   // good orientation ?
    {
        dr = d1 ;
        dr /= dr.two_norm()  ;       // 'x-axis' unit vector
        dn = dr ;
        dn *= -(d2*dr) ;
        dn += d2 ;
        dn /= dn.two_norm()  ;       // 'y-axis' unit vector
    }
    else
    {
        dr = d2 ;
        dr /= dr.two_norm()  ;       // 'y-axis' unit vector
        dn = dr ;
        dn *= -(d1*dr) ;
        dn += d1 ;
        dn /= dn.two_norm()  ;        // 'x-axis' unit vector
    }

    // definition of angles, using projection on the local reference, ie by scalarly multipliying by dr and dn resp.
    for ( size_type j=1 ; j < P.size() ; j++)
    {
        ai.push_back(atan2((P[j]-P[0])*dn,(P[j]-P[0])*dr)) ;
        id.push_back(j) ;
    }

    // sort according to increasing angles
    for ( size_type j=1; j < ai.size(); j++) {
        for ( size_type i=0; i < j; i++) {
            if (ai[j]<ai[i]) {
                std::swap<typename V::value_type>(ai[i],ai[j]) ;
                std::swap<int>(id[i],id[j]) ;
            }
        }
    }

    id.insert(id.begin(),0);
}

template<class CM>
void IntersectionComputation<CM>::orderPointsCC(std::integral_constant<int,2>,
                                                const V& centroid,
                                                std::vector<int>& id,
                                                const  std::vector<V>& P)
{
    typedef typename std::vector<Empty>::size_type size_type;

    // get size_type for all the vectors we are using
    typedef typename std::vector<Empty>::size_type size_type;

    std::vector<typename V::value_type> ai(P.size());
    id.resize(P.size());

    // definition of angles
    for ( size_type i=0; i < P.size(); i++) {
        ai[i] = atan2(P[i][1]-centroid[1],P[i][0]-centroid[0]);
        id[i] = i;
    }

    // sort according to increasing angles
    for ( size_type j=1; j < ai.size(); j++) {
        for ( size_type i=0; i < j; i++) if (ai[j]<ai[i]) {
            std::swap<typename V::value_type>(ai[i],ai[j]);
            std::swap<int>(id[i],id[j]);
        }
    }
}

} /* namespace Dune::GridGlue */
} /* namespace Dune */
