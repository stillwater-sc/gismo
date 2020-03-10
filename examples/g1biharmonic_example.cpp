/** @file biharmonic_multiPatch.cpp

    @brief A Biharmonic example for ONLY TWO-PATCHES

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/
# include <gismo.h>
# include <gsG1Basis/gsG1AuxiliaryEdgeMultiplePatches.h>
# include <gsG1Basis/gsG1AuxiliaryVertexMultiplePatches.h>
# include <gsAssembler/gsG1BiharmonicAssembler.h>
# include <gsG1Basis/gsG1System.h>
# include <gsG1Basis/gsG1Mapper_pascal.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    // Geometry data
    index_t geometry = 0; // Which geometry

    index_t numRefine = 4;
    index_t numDegree = 0;
    index_t regularity = 1;

    // For the spline space of the gluing data
    index_t p_tilde = 1;
    index_t r_tilde = 0;

    index_t threads = 1;

    bool plot = false;
    bool direct = false;
    bool local = false;
    bool loop = false;
    bool local_g1 = false;

    gsCmdLine cmd("Example for solving the biharmonic problem.");
    cmd.addInt("k", "refine", "Number of refinement steps", numRefine);
    cmd.addInt("p", "p_tilde", "Polynomial degree for tilde{p}", p_tilde);
    cmd.addInt("r", "r_tilde", "Regularity for tilde{r}", r_tilde);
    cmd.addSwitch( "plot", "Plot result in ParaView format", plot );
    cmd.addSwitch( "direct", "Construction of the G1 basis functions", direct );
    cmd.addSwitch( "loop", "If you want to solve several levels", loop );
    cmd.addSwitch( "local_g1", "If you want to solve several levels", local_g1 );
    cmd.addInt("g", "geometry", "Geometry", geometry);
    cmd.addInt("t", "threads", "Threads", threads);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // ======= Solution =========
    gsFunctionExpr<> source  ("256*pi*pi*pi*pi*(4*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
    gsFunctionExpr<> laplace ("-16*pi*pi*(2*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
    gsFunctionExpr<> solVal("x",2);
    gsFunctionExpr<>sol1der ("-4*pi*(cos(4*pi*y) - 1)*sin(4*pi*x)",
                             "-4*pi*(cos(4*pi*x) - 1)*sin(4*pi*y)",2);
    gsFunctionExpr<>sol2der ("-16*pi^2*(cos(4*pi*y) - 1)*cos(4*pi*x)",
                             "-16*pi^2*(cos(4*pi*x) - 1)*cos(4*pi*y)",
                             " 16*pi^2*sin(4*pi*x)*sin(4*pi*y)", 2);
    gsFunctionWithDerivatives<real_t> solution(solVal, sol1der, sol2der);

    // ======= Geometry =========
    std::string string_geo;
    switch(geometry)
    {
        case 0:
            string_geo = "planar/multiPatches/4_square_diagonal.xml";
            numDegree = 2; // 2 == degree 3
            break;
        case 1:
            string_geo = "planar/multiPatches/6_square_diagonal.xml";
            numDegree = 2; // 2 == degree 3
            break;
        case 2:
            string_geo = "planar/multiPatches/square_curved.xml";
            numDegree = 0; // 0 == degree 3
            break;
        case 3:
            string_geo = "planar/twoPatches/square_curved.xml";
            numDegree = 0; // 0 == degree 3
            break;
        case 4:
            string_geo = "planar/twoPatches/square_curved_deg_5.xml";
            numDegree = 0; // 0 == degree 3
            break;
        case 5:
            string_geo = "planar/twoPatches/square_curved_deg_7.xml";
            numDegree = 0; // 0 == degree 3
            break;
        case 6:
            string_geo = "planar/twoPatches/square_non_conform.xml";
            numDegree = 2; // 0 == degree 3
            break;
        case 7:
            string_geo = "planar/twoPatches/square_bent.xml";
            numDegree = 0; // 0 == degree 3
            break;
        case 8:
            string_geo = "planar/twoPatches/square_complex_bent.xml";
            numDegree = 0; // 0 == degree 3
            break;
        case 9:
            string_geo = "planar/twoPatches/square_diagonal.xml";
            numDegree = 2; // 2 == degree 3
            break;
        case 10:
            string_geo = "KirchhoffLoveGeo/geo_fivePatchDiffParam.xml";
            numDegree = 2; // 2 == degree 3
            break;
        default:
            gsInfo << "No geometry is used! \n";
            break;
    }

    gsFileData<> fd(string_geo);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    gsMultiPatch<> multiPatch;
    fd.getId(0, multiPatch); // id=0: Multipatch domain
    multiPatch.computeTopology();


    gsWriteParaview(multiPatch,"geometry",5000,true);

    gsOptionList optionList;
    optionList.addInt("p_tilde","Grad",p_tilde);
    optionList.addInt("r_tilde","Reg",r_tilde);
    optionList.addInt("regularity","Regularity of the initial geometry",regularity);
    optionList.addSwitch("local","Local projection for gluing data",local);
    optionList.addSwitch("direct","Local projection for gluing data",direct);
    optionList.addSwitch("plot","Plot in Paraview",plot);
    optionList.addInt("refine","Refinement",numRefine);
    optionList.addInt("degree","Degree",numDegree);

    //multiPatch.patch(1).degreeElevate(1,0);
    multiPatch.degreeElevate(optionList.getInt("degree"));

    multiPatch.uniformRefine_withSameRegularity(optionList.getInt("refine"), optionList.getInt("regularity"));
    gsMultiBasis<> mb(multiPatch);


    gsG1System<real_t> g1System(multiPatch, mb);

    // ########### EDGE FUNCTIONS ###########
    // Interface loop
    for (size_t numInt = 0; numInt < multiPatch.interfaces().size(); numInt++ )
    {
        const boundaryInterface & item = multiPatch.interfaces()[numInt];

        std::string fileName;
        std::string basename = "InterfaceBasisFunctions" + util::to_string(numInt);
        gsParaviewCollection collection(basename);

        /*
        gsG1AuxiliaryEdgeMultiplePatches first(multiPatch, item.first().patch);
        first.computeG1EdgeBasis(optionList,item.first().m_index,false);

        gsG1AuxiliaryEdgeMultiplePatches second(multiPatch, item.second().patch);
        second.computeG1EdgeBasis(optionList,item.second().m_index,false);
        */
        gsG1AuxiliaryEdgeMultiplePatches singleInt(multiPatch, item.first().patch, item.second().patch);
        singleInt.computeG1InterfaceBasis(optionList);

        for (size_t i = 0; i < singleInt.getSinglePatch(0).getG1Basis().nPatches(); i++)
        {
            gsMultiPatch<> edgeSingleBF;

            edgeSingleBF.addPatch(singleInt.getSinglePatch(0).getG1Basis().patch(i));
            edgeSingleBF.addPatch(singleInt.getSinglePatch(1).getG1Basis().patch(i));

            fileName = basename + "_0_" + util::to_string(i);
            gsField<> temp_field(multiPatch.patch(item.first().patch),edgeSingleBF.patch(0));
            gsWriteParaview(temp_field,fileName,5000);
            collection.addTimestep(fileName,i,"0.vts");
            fileName = basename + "_1_" + util::to_string(i);
            gsField<> temp_field_1(multiPatch.patch(item.second().patch),edgeSingleBF.patch(1));
            gsWriteParaview(temp_field_1,fileName,5000);
            collection.addTimestep(fileName,i,"0.vts");
        }
        collection.save();
    }
    // Boundaries loop
    for (size_t numBdy = 0; numBdy < multiPatch.boundaries().size(); numBdy++ )
    {
        const patchSide & bit = multiPatch.boundaries()[numBdy];

        std::string fileName;
        std::string basename = "BoundaryBasisFunctions" + util::to_string(numBdy);
        gsParaviewCollection collection(basename);

        gsInfo << "Patch: " << bit.patch << "\n";
        gsInfo << "m_index: " << bit.m_index << "\n";
        gsG1AuxiliaryEdgeMultiplePatches singleBdy(multiPatch, bit.patch);
        singleBdy.computeG1BoundaryBasis(optionList, bit.m_index);

        for (size_t i = 0; i < singleBdy.getSinglePatch(0).getG1Basis().nPatches(); i++)
        {
            gsMultiPatch<> edgeSingleBF;

            edgeSingleBF.addPatch(singleBdy.getSinglePatch(0).getG1Basis().patch(i));

            fileName = basename + "_0_" + util::to_string(i);
            gsField<> temp_field(multiPatch.patch(bit.patch),edgeSingleBF.patch(0));
            gsWriteParaview(temp_field,fileName,5000);
            collection.addTimestep(fileName,i,"0.vts");
        }
        collection.save();
    }

// NEW NEW NEW NEW NEW NEW NEW NEW NEW

    gsBoundaryConditions<> bcInfo, bcInfo2;
    for (gsMultiPatch<>::const_biterator bit = multiPatch.bBegin(); bit != multiPatch.bEnd(); ++bit)
    {
        bcInfo.addCondition( *bit, condition_type::dirichlet, &solVal ); // = 0
        bcInfo2.addCondition(*bit, condition_type::neumann, &laplace ); // = 0
    }

    // BiharmonicAssembler
    gsG1BiharmonicAssembler<real_t> g1BiharmonicAssembler(multiPatch, mb, bcInfo, bcInfo2, source);
    g1BiharmonicAssembler.assemble();

    //g1BiharmonicAssembler.computeDirichletDofsL2Proj(g1_edges, g1_vertices, g1Mapper );
    //g1BiharmonicAssembler.constructDirichletSolution(g1_edges, g1_vertices, g1Mapper );

    //g1System.setGlobalMapper(g1BiharmonicAssembler.get_mapper());
    //g1System.initialize();
    //g1System.assemble();














} // main