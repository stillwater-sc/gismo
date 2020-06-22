/** @file geometry_example.cpp

    @brief create D-patch coupling

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <iostream>
#include <gismo.h>
#include <typeinfo>

using namespace gismo;

int main(int argc, char* argv[])
{
    bool plot = false;
    index_t numRefine  = 0;
    index_t numElevate = 0;
    index_t geometry = 0;
    std::string input("surfaces/simple.xml");

    gsCmdLine cmd("Tutorial on gsGeometry class.");
    cmd.addPlainString("filename", "G+Smo input geometry file.", input);
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addInt( "g", "geometry", "which geometry",  geometry );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // ======================================================================
    // reading the geometry
    // ======================================================================

    gsMultiPatch<> mpBspline,mp4p;
    // gsReadFile<>(input, mpBspline);

    /*
        |-----------|-----------|
        |           |           |
        |     5     |     4     |
        |           |           |
        |-----------|-----------|
                    |  1  |  3  |
                    |-----------|
                    |  0  |  2  |
                    |-----------|
    */

    if (geometry==0)
    {
        mpBspline.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); //
        mpBspline = mpBspline.uniformSplit(); // patches 0 to 3

        mpBspline.embed(3);
        mpBspline.addAutoBoundaries();
    }
    else if (geometry==1)
    {
        mpBspline.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); //
        mpBspline = mpBspline.uniformSplit(); // patches 0 to 3
        mpBspline.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // patch 4
        mpBspline.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // patch 5
        mpBspline.embed(3);

        gsMatrix<> ones;
        ones = gsMatrix<>::Ones(mpBspline.patch(4).coefs().rows(),1); // patch 4

        // translate second patch to the right
        mpBspline.patch(4).coefs().col(1) += ones; // patch 4

        ones = gsMatrix<>::Ones(mpBspline.patch(5).coefs().rows(),1); // patch 5
        mpBspline.patch(5).coefs().col(0) -= ones; // patch 5
        mpBspline.patch(5).coefs().col(1) += ones; // patch 5


        mpBspline.addInterface(1,boundary::side::north,4,boundary::side::south);
        mpBspline.addInterface(3,boundary::side::north,4,boundary::side::south);
        mpBspline.addInterface(4,boundary::side::west,5,boundary::side::east);

        mpBspline.addAutoBoundaries();
    }
    else
        GISMO_ERROR("Geometry with index "<<geometry<<" unknown.");

    mpBspline.addAutoBoundaries();
    gsDebugVar(mpBspline.nInterfaces());
    gsDebugVar(mpBspline.nBoundary());
    gsDebugVar(mpBspline.nBoxes());

    gsInfo<<mpBspline<<"\n";
    gsInfo<<"Geometry has "<<mpBspline.nPatches()<<" patches\n";
    mpBspline.checkConsistency();
    gsInfo<<mpBspline.getMaxValence()<<"\n";


    // VERTICES WITH MU=2; These are by defintiion only located on the boundaries


    // ALL VERTICES WITH MU>=3

    std::vector< std::vector<patchCorner> > EVs, OVs;

    mpBspline.getEVs(EVs);
    mpBspline.getOVs(OVs);


    // // for each corner in EVs and OVs, get corner list
    // std::vector<patchCorner> surroundingCorners;
    // for (sakjhdasjkdhaksjdhsakjdh)
    //     mpBspline.getCornerList(it,surroundingCorners);

    std::vector<std::vector<patchCorner> > cornerLists;
    std::vector<patchCorner> cornerList;
    patchCorner c;
    patchSide s;
    for(index_t i = 0;i<mpBspline.nPatches();++i)
    {
        for(int j=1;j<=4;++j)
        {
            s=patchSide(i,j);
            if (mpBspline.isInterface(s))
                gsDebug<<"Patch "<<i<<", side "<<j<<" is an interface\n";
            else if (mpBspline.isBoundary(s))
                gsDebug<<"Patch "<<i<<", side "<<j<<" is a boundary\n";
            else
                gsDebug<<"Patch "<<i<<", side "<<j<<" is ???\n";
        }
    }

    for(index_t i = 0;i<mpBspline.nPatches();++i)
    {
        for(int j=1;j<=4;++j)
        {
            c=patchCorner(i,j);
            bool isCycle = mpBspline.getCornerList(c,cornerList);
                gsDebug<<"Patch "<<i<<", corner "<<j<<" has valence "<<cornerList.size()<<"\n";
            bool alreadyReached = false;
            for(size_t k = 0;k<cornerList.size();++k)
                if(cornerList[k].patch<i)
                    alreadyReached = true;
            if(!alreadyReached)
            {
                cornerLists.push_back(cornerList);
            }
        }
    }

    // ----------------------------------------------------------------------

    // p-refine
    if (numElevate!=0)
        mpBspline.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        mpBspline.uniformRefine();

    if (plot)
    {
        gsWriteParaview<>( mpBspline, "mp", 1000, true);
        gsWriteParaview<>( mpBspline.basis(0), "basis", 1000, true);
    }

    if (mpBspline.domainDim()!=3)
        mpBspline.embed(3);

    gsInfo<<"Basis: "<< mpBspline.basis(0)<<"\n";

    gsMultiBasis<> bases(mpBspline);

    gsDofMapper mapper(bases);
    mapper.finalize();
    gsDebugVar(mapper.allFree());
    gsDebugVar(mapper.size());


    gsMatrix<index_t> indices1, indices2, corners;
    index_t p1, p2;
    boxSide s1, s2;
    gsBasis<> * basis1;
    gsBasis<> * basis2;
    boundaryInterface iface;
    corners.resize(2,1);

    // Interior interfaces
    for(gsBoxTopology::const_iiterator iit = mpBspline.iBegin(); iit!= mpBspline.iEnd(); iit++)
    {
        // MATCH THE INTERIOR COEFFICIENTS
        p1 = iit->first().patch;
        p2 = iit->second().patch;

        s1 = iit->first().side();
        s2 = iit->second().side();

        gsInfo<<"Interface between patch "<<p1<<" (side "<<s1<<") and patch "<<p2<<" (side "<<s2<<")\n";

        basis1 = &mpBspline.basis(p1);
        basis2 = &mpBspline.basis(p2);

        mpBspline.getInterface(patchSide(p1,s1),iface);
        gsDebugVar(iface);
        basis1->matchWith(iface,*basis2,indices1,indices2);

        gsInfo<<"\tInterface indices (matching):\n";
        gsInfo<<"\t\tBoundary 1: "<<indices1.transpose()<<"\n";
        gsInfo<<"\t\tBoundary 2: "<<indices2.transpose()<<"\n";

        gsInfo<<"\tNext-to-interface indices (matching):\n";
        indices1 = basis1->boundaryOffset(s1,1);
        indices2 = basis2->boundaryOffset(s2,1);
        gsInfo<<"\t\tBoundary 1: "<<indices1.transpose()<<"\n";
        gsInfo<<"\t\tBoundary 2: "<<indices2.transpose()<<"\n";


        // const index_t b1 = ;
        // const index_t b2 = s2.direction();

        // Flip the index vector if the directions of the indices do not match
            // gsDebugVar(iit->dirOrientation());
            // gsDebugVar(iit->dirMap());
        gsVector<bool> dirOr = iit->dirOrientation();
        for (short_t k = 0; k<mpBspline.dim(); ++k )
        {
            if ( k == s1.direction() ) // skip ?
                continue;

            if ( ! dirOr[k] ) // flip ?
            {
                gsDebug<<"\t\tReversed direction\n";
                indices1.reverse();
            }
        }

        // COUPLE/MATCH COEFFICIENTS of INTERIOR ONLY
        //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            size_t n = indices1.size();
            // get the indices of the basis functions adjacent to the interface
            GISMO_ASSERT(indices1.cols()==1,"Columns of indices should be 1. Someting went wrong");
            indices1 = indices1.block(1,0,n-2,1);
            indices2 = indices2.block(1,0,n-2,1);
            // corners.topRows(1) = indices1.topRows(1);
            // corners.bottomRows(1) = indices1.bottomRows(1);
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


        // TREATMENT FOR THE BOUNDARY CONTROLPOINTS/BASISFUNCTIONS HERE!
        //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        GISMO_ASSERT(indices1.rows()==indices2.rows(),"Indices do not have the same size. Someting went wrong?");
        // for (index_t k=0; k!= indices1.rows(); k++)
        // {

        // }
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    }
    // Boundary interfaces
    for(gsBoxTopology::const_biterator bit = mpBspline.bBegin(); bit!= mpBspline.bEnd(); bit++)
    {
            p1 = bit->patch;
            s1 = bit->side();
            gsDebug<<"Boundary of patch "<<p1<<" (side "<<s1<<")\n";

            basis1 = &mpBspline.basis(p1);
            indices1 = basis1->boundary(s1);
            gsInfo<<"\tBoundary indices: "<<indices1.transpose()<<"\n";
            indices1 = basis1->boundaryOffset(s1,1);
            gsInfo<<"\tNext-to-boundary indices: "<<indices1.transpose()<<"\n";

            // GET THE INDICES OF THE NON-VERTEX BASIS FUNCTIONS
            // QUESTION: THIS ONLY WORKS IF THE CORNERS ARE THE ENDS OF THE INDEX VECTORS. IS THAT THE CASE????
            //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            GISMO_ASSERT(indices1.cols()==1,"Columns of indices should be 1. Someting went wrong");
            n = indices1.size();
            corners.topRows(1) = indices1.topRows(1);
            corners.bottomRows(1) = indices1.bottomRows(1);
            indices1 = indices1.block(1,0,n-2,1);
            //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            // TREATMENT FOR THE BOUNDARY CONTROLPOINTS/BASISFUNCTIONS HERE!
            //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            // for (index_t k=0; k!= indices1.rows(); k++)
            // {

            // }
            //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    }

    // Vertices
    std::vector<bool> marked;
    marked.resize(4*mpBspline.nPatches()); // [p1c1, p1c2,p1c3,p1c4,p2c1,...,p2c4,...,pNc1,...pNc4]
    std::fill(marked.begin(), marked.end(), false);

    index_t cIdx;
    boxCorner corner;
    for(index_t p=0; p!=mpBspline.nPatches(); p++) // for all patches
    {
        gsInfo<<"Patch "<<p<<"\n";
        for (index_t c=1; c!=5; c++) // for all corners [corners run from 1 up to 4]
        {
            if (marked[4*p+c-1])
                continue;

            corner = boxCorner(c);
            gsInfo<<"\tcorner "<<corner<<"\n";

            basis1 = &mpBspline.basis(p);
            cIdx = basis1->functionAtCorner(corner);
            gsInfo<<"\t\tfunction index "<<cIdx<<"\n";


            std::vector<patchCorner> cornerList;
            patchCorner pcorner = patchCorner(p,corner);
            bool isCycle = mpBspline.getCornerList(pcorner,cornerList);

            gsInfo<<"\t\tvalence "<<cornerList.size()<<"\n";
            gsInfo<<"\t\ttype:"<<( isCycle ? " interior vertex " : " boundary vertex " )<<"\n";


            // TREATMENT FOR THE SURROUNDING VERTICES HERE!
            //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            // for (std::vector<patchCorner>::iterator cit = cornerList.begin(); cit != cornerList.end(); cit++)
            // {
                // gsInfo<<" patch "<<*cit<<"; ";
            // }
            //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


            // we can use this to store marked corners later
            marked[4*p + c - 1] = true; //corners run from 1 up to 4
        }
    }

    // ----------------------------------------------------------------------

    // // Cast all patches of the mp object to THB splines
    // gsMultiPatch<> mp;
    // // gsTensorBSpline<2,real_t> *geo;
    // gsTHBSpline<2,real_t> thb;
    // for (index_t k=0; k!=mpBspline.nPatches(); ++k)
    // {
    //     gsTensorBSpline<2,real_t> *geo = dynamic_cast< gsTensorBSpline<2,real_t> * > (&mpBspline.patch(k));
    //     thb = gsTHBSpline<2,real_t>(*geo);
    //     mp.addPatch(thb);
    // }

    // gsMultiBasis<> dbasis(mp);

    // gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    // gsInfo << dbasis.basis(0)<<"\n";




    // gsDebugVar(mp.patch(0).coefs());
    return 0;
}


