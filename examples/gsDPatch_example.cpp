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
    std::string input("surfaces/simple.xml");

    gsCmdLine cmd("Tutorial on gsGeometry class.");
    cmd.addPlainString("filename", "G+Smo input geometry file.", input);
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
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
    size_t n;
    index_t p1, p2;
    boxSide s1, s2;
    corners.resize(2,1);
    for(gsBoxTopology::const_iiterator it = mpBspline.iBegin(); it!= mpBspline.iEnd(); it++)
    {
        // use getContainedCorners(...) in gsBoundary.h?



        // MATCH THE INTERIOR COEFFICNETS
            p1 = it->first().patch;
            p2 = it->second().patch;

            s1 = it->first().side();
            s2 = it->second().side();

            gsDebug<<"Interface between patch "<<p1<<" (side "<<s1<<") and patch "<<p2<<" (side "<<s2<<")\n";

            //>>>>>>
            // COUPLE/MATCH COEFFICIENTS of INTERIOR ONLY

            // get the indices of the basis functions adjacent to the interface


            // gsDebugVar(it->first().patch);
            // gsDebugVar();
            // gsDebugVar(it->second());

            // indices = bases.basis(i).boundary(it);
            // GISMO_ASSERT(indices.cols()==1,"Columns of indices should be 1. Someting went wrong");
            // n = indices.size();
            // corners.topRows(1) = indices.topRows(1);
            // corners.bottomRows(1) = indices.bottomRows(1);
            // indices = indices.block(1,0,n-2,1);

            //<<<<<<

        // FIND INTERIOR VERTICES
            std::vector<boxCorner> corners;
            std::vector<patchCorner> cornerList;
            patchCorner corner;

            s1.getContainedCorners(mpBspline.dim(),corners);
            for (index_t k = 0; k != corners.size(); k++)
            {
                corner = static_cast<patchCorner &>(corners[k]);
                bool isCycle = mpBspline.getCornerList(corner,cornerList);
                gsInfo<<corners[k]<<" (valence = "<<cornerList.size()<<")\n";
            }

            gsInfo<<"----\n";
    }

    for(gsBoxTopology::const_biterator it = mpBspline.bBegin(); it!= mpBspline.bEnd(); it++)
    {
            p1 = it->patch;
            s1 = it->side();
            gsDebug<<"Boundary of patch "<<p1<<" (side "<<s1<<")\n";

            std::vector<boxCorner> corners;
            std::vector<patchCorner> cornerList;
            patchCorner corner;
            s1.getContainedCorners(mpBspline.dim(),corners);
            for (index_t k = 0; k != corners.size(); k++)
            {
                corner = static_cast<patchCorner &>(corners[k]);
                bool isCycle = mpBspline.getCornerList(corner,cornerList);
                gsInfo<<corners[k]<<" (valence = "<<cornerList.size()<<") --> surroundingCorners";
                for (std::vector<patchCorner>::iterator it = cornerList.begin(); it != cornerList.end(); it++)
                {
                    gsInfo<<" patch "<<it->patch<<"; ";
                }
                gsInfo<<"\n";
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


