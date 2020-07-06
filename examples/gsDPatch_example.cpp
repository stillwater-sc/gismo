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


void writeToCSVfile(std::string name, gsMatrix<> matrix)
{
    std::ofstream file(name.c_str());
    for(int  i = 0; i < matrix.rows(); i++){
        for(int j = 0; j < matrix.cols(); j++){
           std::string str = std::to_string(matrix(i,j));
           if(j+1 == matrix.cols()){
               file<<str;
           }else{
               file<<str<<',';
           }
        }
        file<<'\n';
    }
  }

template <class T>
class gsDPatch
{
    protected:
        const gsMultiPatch<T> & m_patches;
        gsMultiBasis<T> m_bases;
        mutable std::vector<boxCorner> m_bcorners;
        mutable std::vector<patchCorner> m_pcorners;
        mutable boxCorner m_bcorner;
        mutable patchCorner m_pcorner;
        index_t m_dim;
        mutable std::vector< std::vector<gsMatrix<index_t> > > m_sides; // structure: m_sides[s][0] is the 0th line of indices of interface s, m_sides[s][k] is the kth line of indices
        mutable std::vector<bool> m_sideCheck;
        mutable std::vector<bool> m_vertCheck;
        mutable std::vector<bool> m_basisCheck;

        mutable gsDofMapper m_mapModified,m_mapOriginal;

        mutable gsSparseMatrix<T> m_matrix;

        mutable size_t m_size;

        mutable std::vector<std::pair<patchCorner,index_t>> m_bVertices;
        mutable std::vector<std::pair<patchCorner,index_t>> m_iVertices;
        // mutable std::vector<patchSide> m_boundaries;
        // mutable std::vector<patchSide> m_interfaces;

    public:
        gsDPatch(const gsMultiPatch<T> & mp) : m_patches(mp), m_dim(mp.dim())
        {
            m_patches.checkConsistency();
            m_dim = m_patches.dim();
            size_t nSides = 2*m_patches.nInterfaces() + m_patches.nBoundary();
            size_t nVerts = 4*m_patches.nPatches();

            m_sides = std::vector<std::vector<gsMatrix<index_t>> >(nSides, std::vector<gsMatrix<index_t>>(2));

            m_sideCheck.resize(nSides);
            std::fill(m_sideCheck.begin(), m_sideCheck.end(), false);
            m_vertCheck.resize(nVerts);
            std::fill(m_vertCheck.begin(), m_vertCheck.end(), false);

            m_bases = gsMultiBasis<T>(m_patches);

            m_mapModified = gsDofMapper(m_bases);
            m_mapOriginal = gsDofMapper(m_bases);

            m_size = -1;
        }

    GISMO_CLONE_FUNCTION(gsDPatch)

        void sideInfo() const
        {
            gsInfo<<"**D-Patch Side info**\n";
            for(index_t i = 0;i<m_patches.nPatches();++i)
            {
                for(int j=1;j<=4;++j)
                {
                    if (m_patches.isInterface(patchSide(i,j)))
                        gsInfo<<"Patch "<<i<<", side "<<j<<" is an interface\n";
                    else if (m_patches.isBoundary(patchSide(i,j)))
                        gsInfo<<"Patch "<<i<<", side "<<j<<" is a boundary\n";
                    else
                        gsInfo<<"Patch "<<i<<", side "<<j<<" is ???\n";
                }
            }
        }

        void cornerInfo() const
        {
            gsInfo<<"**D-Patch Corner info**\n";
            for(index_t i = 0;i<m_patches.nPatches();++i)
            {
                for(int j=1;j<=4;++j)
                {
                    bool isCycle = m_patches.getCornerList(patchCorner(i,j),m_pcorners);
                    gsInfo<<"Patch "<<i<<", corner "<<j<<" has valence "<<m_pcorners.size()<<"\n";
                }
            }
        }

        void checkMesh() const
        {
            /* this function will check the mesh for the following:
                - no T-junction extensions intersect,
                - no one-bay face extension subdivides a face in the 3-disk neighborhood of an extraordinary point, and,
                - no extraordinary point lies within the 3-disk neighborhood of another.
            */
        }

        std::vector<std::pair<patchCorner,index_t>> boundaryVertices() const
        {

        }

        std::vector<std::pair<patchCorner,index_t>> interiorVertices() const
        {

        }

        const index_t indexFromSides(index_t index1, const patchSide side1, index_t index2, const patchSide side2)
        {
            /*
                Finds the index index1 away from side1 and index2 from side2
                index 1 is the index parallel to side 1
                index 2 is the index parallel to side 2
            */
            GISMO_ASSERT(side1.patch==side2.patch,"Sides must be from the same patch");
            GISMO_ASSERT(side1.side().direction()!=side2.side().direction(),"Sides must have different direction");
            index_t index;

            gsBasis<T> * basis = &m_patches.basis(side1.patch);

            gsVector<index_t> indices1 = static_cast<gsVector<index_t>>(basis->boundaryOffset(side1.side(),index2));

            index_t n = indices1.rows();
            if (side1.side()==1) //west
                if (side2.side().parameter()==0) //south
                    index = indices1(index1);
                else
                    index = indices1(n-1-index1); //north
            else if (side1.side()==2) //east
                if (side2.side().parameter()==0) //south
                    index = indices1(index1);
                else
                    index = indices1(n-1-index1); //north
            else if (side1.side()==3) //south
                if (side2.side().parameter()==0) //west
                    index = indices1(index1);
                else
                    index = indices1(n-1-index1); //east
            else if (side1.side()==4) //north
                if (side2.side().parameter()==0) //west
                    index = indices1(index1);
                else
                    index = indices1(n-1-index1); //east

            return index;
        }

        const gsVector<index_t> indicesFromVert(index_t index, const patchCorner corner, const patchSide side, index_t offset = 0)
        {
            /*
                Finds indices i1,...,in in the direction of side away from the vertex

            */

            gsVector<index_t> result(index);

            gsBasis<T> * basis = &m_patches.basis(side.patch);

            gsVector<index_t> indices = static_cast<gsVector<index_t>>(basis->boundaryOffset(side.side(),offset));
            if (side.side()==1) //west
            {
                if (corner.corner()==1)//southwest
                    result = indices.head(index);
                else if (corner.corner()==3) //northwest
                    result = indices.tail(index);
                else
                    GISMO_ERROR(corner.corner() << "is not on side "<<side.side()<<"!");
            }
            else if (side.side()==2) //east
            {
                if (corner.corner()==2)//southeast
                    result = indices.head(index);
                else if (corner.corner()==4) //northeast
                    result = indices.tail(index);
                else
                    GISMO_ERROR(corner.corner() << "is not on side "<<side.side()<<"!");
            }
            else if (side.side()==3) //south
            {
                if (corner.corner()==1)//southwest
                    result = indices.head(index);
                else if (corner.corner()==2) //southeast
                    result = indices.tail(index);
                else
                    GISMO_ERROR(corner.corner() << "is not adjacent to side "<<side.side()<<"!");
            }
            else if (side.side()==4) //north
            {
                if (corner.corner()==3)//northwest
                    result = indices.head(index);
                else if (corner.corner()==4) //northeast
                    result = indices.tail(index);
                else
                    GISMO_ERROR(corner.corner() << "is not adjacent to side "<<side.side()<<"!");
            }
            return result;
        }

        const index_t indexFromVert(index_t index, const patchCorner corner, const patchSide side, index_t offset = 0)
        {
            /*
                Finds indices i in the direction of side away from the vertex
                if index = 0, the corner index is requested
            */

            index_t result;

            gsBasis<T> * basis = &m_patches.basis(corner.patch);

            if ((index==0) && (offset==0))
            {
                result = basis->functionAtCorner(corner.corner());
            }
            else
            {
                gsVector<index_t> indices = static_cast<gsVector<index_t>>(basis->boundaryOffset(side.side(),offset));
                index_t end = indices.rows()-1;
                if (side.side()==1) //west
                {
                    if (corner.corner()==1)//southwest
                        result = indices.at(index);
                    else if (corner.corner()==3) //northwest
                        result = indices.at(end-index);
                    else
                        GISMO_ERROR(corner.corner() << "is not adjacent to side "<<side.side()<<"!");
                }
                else if (side.side()==2) //east
                {
                    if (corner.corner()==2)//southeast
                        result = indices.at(index);
                    else if (corner.corner()==4) //northeast
                        result = indices.at(end-index);
                    else
                        GISMO_ERROR(corner.corner() << "is not adjacent to side "<<side.side()<<"!");
                }
                else if (side.side()==3) //south
                {
                    if (corner.corner()==1)//southwest
                        result = indices.at(index);
                    else if (corner.corner()==2) //southeast
                        result = indices.at(end-index);
                    else
                        GISMO_ERROR(corner.corner() << "is not adjacent to side "<<side.side()<<"!");
                }
                else if (side.side()==4) //north
                {
                    if (corner.corner()==3)//northwest
                        result = indices.at(index);
                    else if (corner.corner()==4) //northeast
                        result = indices.at(end-index);
                    else
                        GISMO_ERROR(corner.corner() << "is not adjacent to side "<<side.side()<<"!");
                }
            }
            return result;
        }
        const index_t indexInCorner(const patchCorner corner) {return indexFromVert(0,corner,patchSide(corner.patch,1)); };

    protected: //functions
        const index_t sideIndex( index_t patch,  boxSide bside)     const   { return 4*patch + bside - 1; }
        const index_t sideIndex( patchSide pside )                  const   { return this->sideIndex(pside.patch, pside.side()); }

        const index_t vertIndex( index_t patch,  boxCorner corner)  const   { return 4*patch + corner -1; }

        const index_t index(index_t patch, index_t i, index_t j) const
        {
            index_t final = 0;
            return final;
        }

        const std::pair<index_t,bool> vertexData(const patchCorner corner) const
        {
            std::vector<patchCorner> corners;
            std::pair<index_t,bool> output;
            output.second = m_patches.getCornerList(corner,corners); // bool is true if interior vertex
            output.first = corners.size();
            return output;
        }
        const index_t getValence( patchCorner corner) const { return this->vertexInfo(corner).first(); }
        const bool isInteriorVertex( patchCorner corner) const { return this->vertexInfo(corner).first(); }

        const void vertexInfo(patchCorner corner) const
        {
            std::pair<index_t,bool> data = this->vertexData(corner);
            gsInfo<<"Patch "<<corner.patch<<", corner "<<corner<<" has valence "<<data.first<<" and is "<<(data.second ? "an interior vertex" : "a boundary vertex")<<"\n";

        }
        const void sideInfo(patchSide side) const
        {
            gsInfo<<"Patch "<<side.patch<<", side "<<side<<" is "<<(m_patches.isBoundary(side) ? "a boundary side" : "an interface")<<"\n";
        }


    public:
        const void matrix_into(gsSparseMatrix<T> & matrix) const { matrix = m_matrix; }
        const gsSparseMatrix<T> matrix() const { gsSparseMatrix<T> result; matrix_into(result); return result; }

        void initialize() // also initialize the mappers!
        {
            size_t nPatches = m_patches.nPatches();
            size_t nInterfaces = m_patches.nInterfaces();
            size_t nBoundaries = m_patches.nBoundary();

            size_t tmp;
            m_size = tmp = 0;
            // number of interior basis functions
            for (index_t k=0; k!=nPatches; k++)
                tmp += (m_patches.basis(k).component(0).size()-4)*(m_patches.basis(k).component(1).size()-4);
            gsDebug<<"Number of interior DoFs: "<<tmp<<"\n";
            m_size += tmp;

            // interfaces
            gsBasis<T> * basis1;
            gsBasis<T> * basis2;
            gsVector<index_t> indices1,indices2;
            tmp = 0;
            for(gsBoxTopology::const_iiterator iit = m_patches.iBegin(); iit!= m_patches.iEnd(); iit++)
            {
                basis1 = &m_patches.basis(iit->first().patch);
                basis2 = &m_patches.basis(iit->second().patch);
                tmp += basis1->boundary(iit->first().side()).size() - 4;
                tmp += basis2->boundary(iit->second().side()).size() - 4;
            }
            gsDebug<<"Number of interface DoFs: "<<tmp<<"\n";
            m_size += tmp;

            // boundaries
            tmp = 0;
            for(gsBoxTopology::const_biterator bit = m_patches.bBegin(); bit!= m_patches.bEnd(); bit++)
            {
                basis1 = &m_patches.basis(bit->patch);
                tmp += 2*(basis1->boundary(bit->side()).size() - 4);
            }
            gsDebug<<"Number of boundary DoFs: "<<tmp<<"\n";
            m_size += tmp;

            // vertices
            tmp = 0;
            std::vector<bool> passed(m_patches.nPatches()*4);
            std::fill(passed.begin(), passed.end(), false);

            std::vector<patchCorner> corners;
            index_t corn = 0;
            for (index_t p=0; p!=m_patches.nPatches(); p++)
                for (index_t c=1; c<=4; c++)
                {
                    index_t idx = vertIndex(p,c);
                    if (!passed.at(idx))
                    {
                        m_patches.getCornerList(patchCorner(p,c),corners);

                        for (index_t k=0; k!=corners.size(); k++)
                            passed.at(vertIndex(corners[k].patch,corners[k])) = true;

                        std::pair<index_t,bool> vdata = vertexData(patchCorner(p,c)); // corner c
                        if (!vdata.second) // boundary vertex
                            tmp += 4;
                        else
                            tmp += vdata.first; // valence;

                        corn +=1;
                    }
                }
            // gsDebug<<"Number of unique corners: "<<corn<<"\n";

            gsDebug<<"Number of vertex DoFs: "<<tmp<<"\n";

            m_size += tmp;

            m_matrix.resize(m_size,m_bases.totalSize());
        }

        void initializeMapper() // also initialize the mappers!
        {
            // interfaces
            std::fill(m_sideCheck.begin(), m_sideCheck.end(), false);
            std::fill(m_vertCheck.begin(), m_vertCheck.end(), false);

            gsBasis<T> * basis;

            std::vector<index_t> patches(2);
            std::vector<patchSide> psides(2);
            gsVector<index_t> indices;
            std::vector<patchCorner> pcorners;
            patchCorner pcorner;
            size_t end;
            index_t cidx, sidx;
            std::pair<index_t,bool> vdata1, vdata2;
            index_t jmin, jmax;

            for(gsBoxTopology::const_iiterator iit = m_patches.iBegin(); iit!= m_patches.iEnd(); iit++)
            {
                sidx = sideIndex( iit->second().patch,iit->second().side());
                if (m_sideCheck.at(sidx))
                    continue;

                patches[0] = iit->first().patch;
                patches[1] = iit->second().patch;
                psides[0] = patchSide(iit->first().patch,iit->first().side()); // the interface on the first patch
                psides[1] = patchSide(iit->second().patch,iit->second().side()); // the interface on the second patch

                for (index_t p = 0; p != 2; p++)
                {
                    sidx = sideIndex( patches[p] ,psides[p] );
                    if (m_sideCheck.at(sidx))
                        continue;

                    /*
                        o o o @ X | X @ o o o               X: eliminated DoFs by interface/boundary rule
                        o o o @ X | X @ o o o               @: modified DoFs by interface rule
                        o o o ? ? | ? ? o o o               o: preserved DoFs (interior)
                        @ @ ? x x | x x ? @ @               ?: Depends on the vertex (X and @ if not (interior vertex & valence = 4))
                        X X ? x x | x x ? X X               x: handled in vertex rule
                        ----------|----------
                                  | x x ? X X
                                  | x x ? @ @
                                  | ? ? o o o
                                  | X @ o o o
                                  | X @ o o o

                    */
                    basis = &m_patches.basis(patches[p]);
                    indices = static_cast<gsVector<index_t>>( basis->boundary(psides[p]) );
                    end = indices.size()-1;

                    patchSide(iit->first().patch,iit->first().side()).getContainedCorners(m_dim,pcorners);
                    vdata1 = this->vertexData(pcorners[0]);
                    vdata2 = this->vertexData(pcorners[1]);

                    jmin = 2;
                    jmax = end-2;
                    if (vdata1.first!=4 && vdata1.second) // corner 1 is interior vertex with valence unequal to 4
                    {
                        cidx = basis->functionAtCorner(pcorners[0]);
                        if (cidx==indices.at(0))
                            jmin = 3;
                        else if (cidx==indices.at(indices.size()-1))
                            jmax = end-3;
                        else
                            GISMO_ERROR("Place unknown?");
                    }
                    else if (vdata2.first!=4 && vdata2.second) // corner 2 is interior vertex with valence unequal to 4
                    {
                        cidx = basis->functionAtCorner(pcorners[0]);
                        if (cidx==indices.at(0))
                            jmin = 3;
                        else if (cidx==indices.at(indices.size()-1))
                            jmax = end-3;
                        else
                            GISMO_ERROR("Place unknown?");
                    }
                    indices = indices.segment(jmin,jmax-jmin+1); // start at jmin and have length jmax-jmin+1
                    m_mapModified.markBoundary(patches[p], indices);

                    gsDebug<<"Eliminated "<<indices.transpose()<<" of basis "<<patches[p]<<"\n";
                    m_sideCheck.at(sidx) = true;
                }
            }

            // boundaries
            for(gsBoxTopology::const_biterator bit = m_patches.bBegin(); bit!= m_patches.bEnd(); bit++)
            {
                sidx = sideIndex(bit->patch,bit->side());
                m_sideCheck.at(sidx) = true;
            }

            // vertices
            for (index_t p=0; p!=m_patches.nPatches(); p++)
            {
                for (index_t c=1; c<=4; c++)
                {
                    cidx = vertIndex(p,c);
                    if (m_vertCheck.at(cidx))
                        continue;

                    pcorner = patchCorner(p,c);
                    this->vertexInfo(pcorner);

                    // get valence and vertex info of corner
                    std::pair<index_t,bool> vdata = vertexData(pcorner); // corner c
                    if (!vdata.second) // boundary vertex
                    {
                        if (vdata.first==2) //valence = 2
                        {
                            /*
                                o o o @ X |i| X @ o o o                 x: eliminated DoFs by vertex rule
                                o o o @ X |n| X @ o o o                 X: eliminated DoFs by interface/boundary rule
                                o o o @ X |t| X @ o o o                 o: preserved DoFs (interior)
                                o o o * x |e| x * o o o                 @: modified DoFs by interface rule
                                o o o * x |r| x * o o o                 *: modified DoFs by vertex rule (unique DoFs)
                                -----------------------                 %: modified DoFs by vertex rule (matched DoFs)
                                -boundary-| | -boundary-
                                -----------------------
                            */
                            // we mark the nodes belonging to the interface
                            pcorner.getContainingSides(m_dim,psides);
                            for (index_t p=0; p!=psides.size(); p++)
                            {
                                if (m_patches.isInterface(psides[p]))
                                {
                                    m_mapModified.markBoundary(pcorner.patch,indicesFromVert(2,pcorner,psides[p]));
                                    gsDebug<<"Eliminated "<<indicesFromVert(2,pcorner,psides[p]).transpose()<<" of basis "<<pcorner.patch<<"\n";
                                }

                            }
                        }
                        else if (vdata.first==3) //valence = 3
                        {
                            /*
                                o o o @ X |i| X @ o o o                 x: eliminated DoFs by vertex rule
                                o o o @ X |n| X @ o o o                 X: eliminated DoFs by interface/boundary rule
                                o o o @ X |t| X @ o o o                 o: preserved DoFs (interior)
                                @ @ @ * x |e| x * @ @ @                 @: modified DoFs by interface rule
                                X X X x % |r| % x X X X                 *: modified DoFs by vertex rule (unique DoFs)
                                -----------------------                 %: modified DoFs by vertex rule (matched DoFs)
                                -boundary-| |-interface
                                -----------------------
                                         b| | % x X X X
                                         o| | x * @ @ @
                                         u| | X @ o o o
                                         n| | X @ o o o
                                         d| | X @ o o o

                            */
                            // We handle all corners associated to pcorner
                            m_patches.getCornerList(pcorner,pcorners);
                            for (index_t c=0; c!=pcorners.size(); c++)
                            {
                                // Eliminate their 0,1 and 1,0 vertices
                                pcorners[c].getContainingSides(m_dim,psides);
                                m_mapModified.eliminateDof(indexFromVert(1,pcorners[c],psides[0]),pcorners[c].patch);
                                m_mapModified.eliminateDof(indexFromVert(1,pcorners[c],psides[1]),pcorners[c].patch);
                                // m_mapModified.markTagged(indexFromSides(1,psides[0],1,psides[1]),pcorners[c].patch);
                                gsDebug<<"Eliminated "<<indexFromSides(1,psides[0],1,psides[1])<<" of basis "<<pcorners[c].patch<<"\n";

                                // And match the 0,0 vertex (i.e. the corner) to the corner that is first in the list pcorners.
                                if (c!=0)
                                {
                                    patchSide pseudo = patchSide(pcorner.patch,1); // this side does not contribute since we use index = 0 in indexFromVert
                                    // m_mapModified.matchDof(pcorners[0].patch,indexFromVert(0,pcorners[0],pseudo),pcorners[c].patch,indexFromVert(0,pcorners[c],pseudo));
                                    gsDebug<<"Matched "<<indexFromVert(0,pcorners[0],pseudo)<<" of basis "<<pcorners[0].patch<<" with "<<indexFromVert(0,pcorners[c],pseudo)<<" of basis "<<pcorners[c].patch<<"\n";

                                }
                                // mark the vertex as passed
                                m_vertCheck[ vertIndex(pcorners[c].patch, pcorners[c].corner()) ] = true;
                            }
                        }
                        else if (vdata.first==1) //valence = 1
                        {  } // do nothing
                        else
                            GISMO_ERROR("Boundary vertex with valence = "<<vdata.first<<" has no implementation");
                    }
                    else // interior vertex
                    {
                        if (vdata.first==4)
                        {
                            /*
                                o o o @ X |i| X @ o o o                 x: eliminated DoFs by vertex rule
                                o o o @ X |n| X @ o o o                 X: eliminated DoFs by interface/boundary rule
                                o o o @ X |t| X @ o o o                 o: preserved DoFs (interior)
                                @ @ @ * x |e| x * @ @ @                 @: modified DoFs by interface rule
                                X X X x x |r| x x X X X                 *: modified DoFs by vertex rule (unique DoFs)
                                -----------------------
                                -boundary-| |-interface
                                -----------------------
                                X X X x x |i| x x X X X
                                @ @ @ * x |n| x * @ @ @
                                o o o @ X |t| X @ o o o
                                o o o @ X |e| X @ o o o
                                o o o @ X |r| X @ o o o

                            */
                            // we mark the nodes belonging to the interfaces (both sides bordering the vertex)
                            pcorner.getContainingSides(m_dim,psides);
                            for (index_t p=0; p!=psides.size(); p++)
                            {
                                m_mapModified.markBoundary(pcorner.patch,indicesFromVert(2,pcorner,psides[p]));
                                gsDebug<<"Eliminated "<<indicesFromVert(2,pcorner,psides[p]).transpose()<<" of basis "<<pcorner.patch<<"\n";
                            }
                        }
                        else
                            GISMO_ERROR("Extraordinary vertex!! -- to be implemented\n");
                    }

                    // label vertex as processed
                    m_vertCheck[ cidx ] = true;
                }
            }
            m_mapModified.finalize();
            m_mapOriginal.finalize();
            gsDebugVar(m_mapModified.freeSize());
            // gsDebugVar(m_mapModified.coupledSize());
            // gsDebugVar(m_mapModified.boundarySize());

            gsDebugVar(m_mapModified.size());
            gsDebugVar(m_bases.totalSize());
        }


        void computeMatrix()
        {
            // TO DO; IMPLEMENT THE CHECKS
            m_basisCheck.resize(m_size);
            std::fill(m_basisCheck.begin(), m_basisCheck.end(), false);

            std::fill(m_sideCheck.begin(), m_sideCheck.end(), false);
            std::fill(m_vertCheck.begin(), m_vertCheck.end(), false);

            //

            std::vector<std::vector<patchCorner> > cornerLists;
            m_patches.getEVs(cornerLists);
            if (cornerLists.size()!=0)
                GISMO_ERROR("Refinement procedure needed since mesh contains an extraordinary vertex! This is currently not implemented....");

            m_matrix.resize( m_mapModified.freeSize(), m_bases.totalSize() );

            std::vector<gsBasis<T> *> basis(2);
            gsBasis<T> * basis1;
            gsBasis<T> * basis2;
            std::vector<gsMatrix<index_t>> indices(2); // interface indices
            std::vector<gsMatrix<index_t>> oindices(2); // interface indices (with offset 1)
            gsMatrix<index_t> tmpIndices;
            size_t n, m;
            gsVector<bool> dirOr;

            boundaryInterface iface;
            std::vector<patchCorner> pcorners;
            patchCorner pcorner;
            std::vector<boxSide> bsides;
            std::vector<patchSide> psides(2);
            std::vector<index_t> patches(2);
            std::vector<index_t> minJ(2);
            std::vector<index_t> maxJ(2);
            index_t cidx, sidx;
            index_t colIdx,rowIdx1,rowIdx2;

            for(gsBoxTopology::const_iiterator iit = m_patches.iBegin(); iit!= m_patches.iEnd(); iit++)
            {
                sideInfo(iit->first());
                sideInfo(iit->second());

                // get the bases belonging to both patches
                basis[0] = &m_patches.basis(iit->first().patch);
                basis[1] = &m_patches.basis(iit->second().patch);

                // get a boundaryInterface object for the interface and obtain the matched indices along this interface
                m_patches.getInterface(patchSide(iit->first().patch,iit->first().side()),iface);
                basis[0]->matchWith(iface,*basis[1],indices[0],indices[1]);

                // now we treat the offset of the interface
                oindices[0] = basis[0]->boundaryOffset(iit->first().side(),1);
                oindices[1] = basis[1]->boundaryOffset(iit->second().side(),1);

                // Flip the index vector if the directions of the indices do not match
                dirOr = iit->dirOrientation();
                for (short_t k = 0; k<m_patches.dim(); ++k )
                {
                    if ( k == iit->first().side().direction() ) // skip ?
                        continue;

                    if ( ! dirOr[k] ) // flip ?
                    {
                        gsDebug<<"\t\tReversed direction\n";
                        oindices[0].reverseInPlace();
                    }
                }

                minJ[0] = minJ[1] = 2;
                maxJ[0] = indices[0].size()-3;
                maxJ[1] = indices[1].size()-3;

                // store relevant data in vectors with size 2 for the loop over the associated corners
                patches[0] = iit->first().patch;
                patches[1] = iit->second().patch;
                psides[0] = patchSide(iit->first().patch,iit->first().side()); // the interface on the first patch
                psides[1] = patchSide(iit->second().patch,iit->second().side()); // the interface on the second patch

                psides[0].getContainedCorners(m_dim,pcorners); //corners on first patch
                for (index_t c=0; c!=pcorners.size(); c++)
                {
                    index_t idx1 = vertIndex(patches[0],pcorners[c]) ;
                    index_t idx2 = vertIndex(patches[0],pcorners[c]) ;
                    // label vertex as processed
                    if(m_vertCheck[ idx1 ] && m_vertCheck[ idx2] )
                        continue;

                    this->vertexInfo(patchCorner(patches[0],pcorners[c]));

                    // get valence and vertex info of corner
                    std::pair<index_t,bool> vdata = vertexData(pcorners[c]); // corner c
                    if (!vdata.second) // boundary vertex
                    {
                        if (vdata.first==1) //valence = 1
                            gsInfo<<"this case should not exist... (interface handling)\n";
                        else if (vdata.first==2) //valence = 2
                        {
                            /*
                                o o o @ X |i| X @ o o o                 x: eliminated DoFs by vertex rule
                                o o o @ X |n| X @ o o o                 X: eliminated DoFs by interface/boundary rule
                                o o o @ X |t| X @ o o o                 o: preserved DoFs (interior)
                                o o o * x |e| x * o o o                 @: modified DoFs by interface rule
                                o o o * x |r| x * o o o                 *: modified DoFs by vertex rule (unique DoFs)
                                -----------------------                 %: modified DoFs by vertex rule (matched DoFs)
                                -boundary-| | -boundary-
                                -----------------------
                            */
                            std::vector<patchCorner> otherCorners;
                            m_patches.getCornerList(pcorners[c],otherCorners);
                            patchCorner otherCorner = otherCorners[1]; // because the first entry is pcorner[c]
                            for (index_t i=0; i!=2; i++)
                            {
                                // get indices
                                index_t bj0_p1 = indexFromVert(i,pcorners[c],psides[0],0);
                                index_t bj1_p1 = indexFromVert(i,pcorners[c],psides[0],1);
                                index_t bj0_p2 = indexFromVert(i,otherCorner,psides[1],0);
                                index_t bj1_p2 = indexFromVert(i,otherCorner,psides[1],1);

                                rowIdx1 = m_mapModified.index(bj1_p1,patches[0]);
                                rowIdx2 = m_mapModified.index(bj1_p2,patches[1]);

                                gsDebugVar(rowIdx1);
                                gsDebugVar(rowIdx2);

                                colIdx = m_mapOriginal.index(bj0_p1,patches[0]);
                                m_matrix(rowIdx1,colIdx) = m_matrix(rowIdx2,colIdx) = 0.5;
                                colIdx = m_mapOriginal.index(bj0_p2,patches[1]);
                                m_matrix(rowIdx2,colIdx) = m_matrix(rowIdx2,colIdx) = 0.5;
                            }
                        }
                        else if (vdata.first==3) //valence = 3
                        {
                            /*
                                o o o @ X | X @ o o o               x: eliminated DoFs by vertex rule
                                o o o @ X | X @ o o o               X: eliminated DoFs by interface/boundary rule
                                o o o @ X | X @ o o o               o: preserved DoFs (interior)
                                @ @ @ * x | x * @ @ @               @: modified DoFs by interface rule
                                X X X x % | % x X X X               *: modified DoFs by vertex rule (unique DoFs)
                                ----------|----------               %: modified DoFs by vertex rule (matched DoFs)
                                          | % x X X X
                                          | x * @ @ @
                                          | X @ o o o
                                          | X @ o o o
                                          | X @ o o o

                            */
                            // the (0,0) basis functions (relative to the corner) of all patches have weight 1/4
                            gsBasis<T> * tmpBasis;
                            std::vector<patchCorner> otherCorners;
                            m_patches.getCornerList(pcorners[c],otherCorners);
                            // find the patch corner which shares the interface
                            patchCorner otherCorner;
                            if (otherCorners[1].patch == patches[1])
                                otherCorner = otherCorners[1];
                            else if (otherCorners[2].patch == patches[1])
                                otherCorner = otherCorners[2];
                            else
                                GISMO_ERROR("HUH?");

                            index_t b10_p1 = indexFromVert(1,pcorners[c],psides[0],0); // index from vertex pcorners[c] along side psides[0] with offset 0.
                            index_t b11_p1 = indexFromVert(1,pcorners[c],psides[0],1); // point 1,1
                            index_t b10_p2 = indexFromVert(1,otherCorner,psides[1],0); // point 0,1
                            index_t b11_p2 = indexFromVert(1,otherCorner,psides[1],1);

                            gsDebug<<"DoF (1,1) of patch "<<patches[0]<<" has index "<<b11_p1<<"\n";
                            gsDebug<<"DoF (1,1) of patch "<<patches[1]<<" has index "<<b11_p2<<"\n";

                            rowIdx1 = m_mapModified.index(b11_p1,patches[0]);
                            rowIdx2 = m_mapModified.index(b11_p2,patches[1]);
                            // 1. give the basis function a weight 1 from itself
                            colIdx = m_mapOriginal.index(b11_p1,patches[0]);
                            m_matrix(rowIdx1,colIdx) = 1;
                            colIdx = m_mapOriginal.index(b11_p2,patches[1]);
                            m_matrix(rowIdx2,colIdx) = 1;

                            // 2. give the basis function a weight 1/4 from the other (0,0) basis functions
                            index_t index;
                            for (index_t k =0; k!=otherCorners.size(); k++)
                            {
                                tmpBasis = &m_patches.basis(otherCorners[k].patch);
                                index = tmpBasis->functionAtCorner(otherCorners[k].corner());
                                colIdx = m_mapOriginal.index(index,otherCorners[k].patch);
                                if (!m_basisCheck[rowIdx1])
                                    m_matrix(rowIdx1,colIdx) = 0.25;
                                if (!m_basisCheck[rowIdx2])
                                    m_matrix(rowIdx2,colIdx) = 0.25;
                            }

                            // 3. add weight 1/2 from the interface bcs.
                            colIdx = m_mapOriginal.index(b10_p1,patches[0]);
                            m_matrix(rowIdx1,colIdx) = 0.5;
                            m_matrix(rowIdx2,colIdx) = 0.5;

                            colIdx = m_mapOriginal.index(b10_p2,patches[1]);
                            m_matrix(rowIdx1,colIdx) = 0.5;
                            m_matrix(rowIdx2,colIdx) = 0.5;

                            gsInfo<<"to be implemented\n";
                        }
                        else
                            GISMO_ERROR("Boundary vertex with valence = "<<vdata.first<<" has no implementation");
                    }
                    else // interior vertex
                    {
                        if (vdata.first==4)
                        {
                            gsInfo<<"Ordinary vertex!! (interior vertex with valence = 4) -- to be implemented\n";
                        }
                        else
                        {
                            // patchSide(iit->first().patch,iit->first().side()).getContainedCorners(m_dim,pcorners);
                            // vdata = this->vertexData(pcorners[0]);

                            // if (vdata.first!=4 && vdata.second) // corner 1 is interior vertex with valence unequal to 4
                            // {
                            //     cidx = basis->functionAtCorner(pcorners[0]);
                            //     if (cidx==indices.at(0))
                            //         jmin = 3;
                            //     else if (cidx==indices.at(indices.size()-1))
                            //         jmax = end-3;
                            //     else
                            //         GISMO_ERROR("Place unknown?");
                            // }
                            // else if (vdata2.first!=4 && vdata2.second) // corner 2 is interior vertex with valence unequal to 4
                            // {
                            //     cidx = basis->functionAtCorner(pcorners[0]);
                            //     if (cidx==indices.at(0))
                            //         jmin = 3;
                            //     else if (cidx==indices.at(indices.size()-1))
                            //         jmax = end-3;
                            //     else
                            //         GISMO_ERROR("Place unknown?");
                            // }
                            // cidx = basis[p]->functionAtCorner(pcorners[c]);
                            // if (cidx==indices[p].at(0))
                            //     minJ[p] = 3;
                            // else if (cidx==indices[p].at(indices[p].size()-1))
                            //     maxJ[p] = indices[p].size()-4;
                            // else
                            //     GISMO_ERROR("Place unknown?");

                            gsInfo<<"Extraordinary vertex!! -- to be implemented\n";
                        }
                    }
                    // label vertex as processed
                    m_vertCheck[ idx1 ] = true;
                    m_vertCheck[ idx2 ] = true;
                }
                /*
                    Now we treat the interior degrees of freedom. For all interface DoFs, a linear combination of basis functions is constructed and the DoFs located on the interface are removed from the mapper.
                */
                // indices[p] = m_sides.at( sideIndex(patches[p],psides[p]) ).at(0) = indices[p].block(minJ[p],0,maxJ[p]-minJ[p]+1,1);



                // TO DO: fill H-matrix entries.

                // label side as processed
                m_sideCheck[ sideIndex(patches[0], psides[0]) ] = true;
                m_sideCheck[ sideIndex(patches[1], psides[1]) ] = true;

                m_basisCheck[rowIdx1] = true;
                m_basisCheck[rowIdx2] = true;
            }



            // boundaries
            for(gsBoxTopology::const_biterator bit = m_patches.bBegin(); bit!= m_patches.bEnd(); bit++)
            {
                // Mark sides as done
                m_sideCheck[ sideIndex(bit->patch, bit->side()) ] = true;
            }

            patchCorner corner;
            gsInfo<<"Last corners\n";
            // iterate over the other vertices
            for (index_t p=0; p!=m_patches.nPatches(); p++)
                for (index_t c=1; c<=4; c++)
                {
                    if (m_vertCheck[ vertIndex(p,c) ]) // vertex already covered?
                        continue;

                    this->vertexInfo(patchCorner(p,c));


                    std::pair<index_t,bool> vdata = vertexData(patchCorner(p,c)); // corner c
                    if (!vdata.second) // boundary vertex
                    {
                        if (vdata.first!=1) //valence = 1
                            gsInfo<<"this case should not exist...\n";
                    }
                    else
                        gsInfo<<"this case should not exist...\n";

                    // label vertex as processed
                    m_vertCheck[ vertIndex(p,c) ] = true;
                }

            // gsDebugVar(m_matrix.toDense());

            bool checkSides = std::all_of(m_sideCheck.begin(), m_sideCheck.end(), [](bool m_sideCheck) { return m_sideCheck; });
            GISMO_ASSERT(checkSides,"Not all sides are checked");
            bool checkVerts = std::all_of(m_vertCheck.begin(), m_vertCheck.end(), [](bool m_vertCheck) { return m_vertCheck; });
            GISMO_ASSERT(checkVerts,"Not all sides are checked");

            gsDebugVar(m_matrix.toDense());
            gsDebugVar(m_mapModified.index(0,1));

        }


        void getObjects()
        {
            patchSide side;
            patchCorner corner;
            std::vector<patchCorner> pcorners;
            std::pair<patchCorner,index_t> vertex;
            bool interior;
            for (index_t p=0; p!=m_patches.nPatches(); p++)
                for (index_t s=1; s<=4; s++)
                {
                    //     side = patchSide(p,s);
                    //     if (m_patches.isBoundary(side))
                    //         m_boundaries.push_back(side);
                    //     else if (m_patches.isInterface(side))
                    //         m_interfaces.push_back(side);

                    corner = patchCorner(p,s);
                    interior = m_patches.getCornerList(corner,pcorners);


                    vertex.first() = corner;
                    vertex.second() = pcorners.size(); // valence
                    if (interior)
                        m_iVertices.push_back(vertex);
                    else
                        m_bVertices.push_back(vertex);
                }
        }

        void constructMap()
        {

                index_t p1, p2;
                gsBasis<T> * basis1;
                gsBasis<T> * basis2;
                gsMatrix<index_t> indices1,indices2, oindices1, oindices2;
                size_t n, m;
                gsVector<bool> dirOr;

                boundaryInterface iface;
                std::vector<patchCorner> corners;
                std::vector<boxSide> bsides;
                patchSide side;

                index_t jmin,jmax;
                // interfaces
                for(gsBoxTopology::const_iiterator iit = m_patches.iBegin(); iit!= m_patches.iEnd(); iit++)
                {
                    basis1 = &m_patches.basis(iit->first().patch);
                    basis2 = &m_patches.basis(iit->second().patch);
                    side = patchSide(iit->first().patch,iit->first().side());
                    side.getContainedCorners(m_dim,corners);

                    m_patches.getInterface(patchSide(iit->first().patch,iit->first().side()),iface);
                    basis1->matchWith(iface,*basis2,indices1,indices2);

                    // remove boundary vertices and write to m_sides
                    n = indices1.size();
                    m = indices2.size();

                    m_sides.at(p1).at(0) = indices1.block(3,0,n-4,1);
                    m_sides.at(p2).at(0) = indices2.block(3,0,n-4,1);

                    // eliminate the DoFs corresponding to the interface from the mapper.
                    m_mapModified.markBoundary(iit->first().patch,m_sides.at(p1).at(0));
                    m_mapModified.markBoundary(iit->second().patch,m_sides.at(p2).at(0));

                    // now we treat the offset of the interface
                    oindices1 = basis1->boundaryOffset(iit->first().side(),1);
                    oindices2 = basis2->boundaryOffset(iit->second().side(),1);

                    // Flip the index vector if the directions of the indices do not match
                    dirOr = iit->dirOrientation();
                    for (short_t k = 0; k<m_patches.dim(); ++k )
                    {
                        if ( k == iit->first().side().direction() ) // skip ?
                            continue;

                        if ( ! dirOr[k] ) // flip ?
                        {
                            gsDebug<<"\t\tReversed direction\n";
                            oindices1.reverseInPlace();
                        }
                    }

                    n = oindices1.size();
                    m_sides.at(p1).at(1) = oindices1.block(3,0,n-4,1);
                    m = oindices2.size();
                    m_sides.at(p2).at(1) = oindices2.block(3,0,n-4,1);

                    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    // Mark sides as done
                    m_sideCheck.at(p1) = m_sideCheck.at(p2) = true;
                }

                // boundaries
                for(gsBoxTopology::const_biterator bit = m_patches.bBegin(); bit!= m_patches.bEnd(); bit++)
                {
                    p1 = sideIndex(bit->patch, bit->side());

                    basis1 = &m_patches.basis(bit->patch);
                    indices1 = basis1->boundary(bit->side());

                    // GET THE INDICES OF THE NON-VERTEX BASIS FUNCTIONS
                    // QUESTION: THIS ONLY WORKS IF THE CORNERS ARE THE ENDS OF THE INDEX VECTORS. IS THAT THE CASE????
                    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    n = indices1.size();
                    // get the indices of the basis functions adjacent to the interface
                    m_sides.at(p1).at(0) = indices1.block(1,0,n-2,1);
                    n = oindices1.size();
                    oindices1 = basis1->boundaryOffset(bit->side(),1);
                    m_sides.at(p1).at(1) = oindices1.block(1,0,n-2,1);

                    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    // Mark sides as done
                    m_sideCheck.at(p1) = true;
                }

                // boundary vertices
                for (std::vector<std::pair<patchCorner,index_t>>::iterator it = m_bVertices.begin(); it!=m_bVertices.end(); ++it)
                {
                    basis1 = &m_patches.basis(it->first().patch);
                    if(it->second()==1)
                    {
                        gsInfo<<"Keep the dofs\n";
                    }
                    else if (it->second()==2)
                    {
                        gsInfo<<"Eliminate the dofs on the interface\n";
                        // eliminate the vertex itself
                        m_mapModified.eliminateDof(basis1.functionAtCorner(it->first()));

                        it->first().getContainingSides(m_dim,bsides); // is always two
                        // if (m_patches.isInterface(bsides[0])) // then this is the interface
                        // else// the other side is the interface
                    }
                    else if (it->second()==3)
                    {

                    }
                    else
                        GISMO_ERROR("No rule for boundary vertex of valence " <<it->second()<<". Something went wrong?");
                }
                // interface vertices
                for (std::vector<std::pair<patchCorner,index_t>>::iterator it = m_iVertices.begin(); it!=m_iVertices.end(); ++it)
                {

                }



                bool check = std::all_of(m_sideCheck.begin(), m_sideCheck.end(), [](bool m_sideCheck) { return m_sideCheck; });
                GISMO_ASSERT(check,"Not all sides are checked");
            }

        void computeSides()
        {
            std::fill(m_sideCheck.begin(), m_sideCheck.end(), false);

            index_t p1, p2;
            gsBasis<T> * basis1;
            gsBasis<T> * basis2;
            gsMatrix<index_t> indices1,indices2, oindices1, oindices2;
            size_t n;
            gsVector<bool> dirOr;

            boundaryInterface iface;

            // interfaces
            for(gsBoxTopology::const_iiterator iit = m_patches.iBegin(); iit!= m_patches.iEnd(); iit++)
            {
                // MATCH THE INTERIOR COEFFICIENTS
                p1 = sideIndex(iit->first().patch, iit->first().side());
                p2 = sideIndex(iit->second().patch, iit->second().side());

                gsDebugVar(m_sides.at(0).size());

                // gsInfo<<"Interface between patch "<<p1<<" (side "<<s1<<") and patch "<<p2<<" (side "<<s2<<")\n";

                basis1 = &m_patches.basis(iit->first().patch);
                basis2 = &m_patches.basis(iit->second().patch);

                m_patches.getInterface(patchSide(iit->first().patch,iit->first().side()),iface);
                basis1->matchWith(iface,*basis2,indices1,indices2);

                // remove boundary vertices and write to m_sides
                n = indices1.size();
                m_sides.at(p1).at(0) = indices1.block(1,0,n-2,1);
                n = indices2.size();
                m_sides.at(p2).at(0) = indices2.block(1,0,n-2,1);

                // now we treat the offset of the interface
                oindices1 = basis1->boundaryOffset(iit->first().side(),1);
                oindices2 = basis2->boundaryOffset(iit->second().side(),1);

                // Flip the index vector if the directions of the indices do not match
                dirOr = iit->dirOrientation();
                for (short_t k = 0; k<m_patches.dim(); ++k )
                {
                    if ( k == iit->first().side().direction() ) // skip ?
                        continue;

                    if ( ! dirOr[k] ) // flip ?
                    {
                        gsDebug<<"\t\tReversed direction\n";
                        oindices1.reverseInPlace();
                    }
                }

                n = oindices1.size();
                m_sides.at(p1).at(1) = oindices1.block(1,0,n-2,1);
                n = oindices2.size();
                m_sides.at(p2).at(1) = oindices2.block(1,0,n-2,1);

                // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                // Mark sides as done
                m_sideCheck.at(p1) = m_sideCheck.at(p2) = true;
            }

            // boundaries
            for(gsBoxTopology::const_biterator bit = m_patches.bBegin(); bit!= m_patches.bEnd(); bit++)
            {
                p1 = sideIndex(bit->patch, bit->side());

                basis1 = &m_patches.basis(bit->patch);
                indices1 = basis1->boundary(bit->side());

                // GET THE INDICES OF THE NON-VERTEX BASIS FUNCTIONS
                // QUESTION: THIS ONLY WORKS IF THE CORNERS ARE THE ENDS OF THE INDEX VECTORS. IS THAT THE CASE????
                //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                n = indices1.size();
                // get the indices of the basis functions adjacent to the interface
                m_sides.at(p1).at(0) = indices1.block(1,0,n-2,1);
                n = oindices1.size();
                oindices1 = basis1->boundaryOffset(bit->side(),1);
                m_sides.at(p1).at(1) = oindices1.block(1,0,n-2,1);



                // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                // Mark sides as done
                m_sideCheck.at(p1) = true;
            }

            bool check = std::all_of(m_sideCheck.begin(), m_sideCheck.end(), [](bool m_sideCheck) { return m_sideCheck; });
            GISMO_ASSERT(check,"Not all sides are checked");
        }
};

int main(int argc, char* argv[])
{
    bool plot = false;
    bool write = false;
    index_t numRefine  = 0;
    index_t numElevate = 0;
    index_t geometry = 0;
    std::string input;

    gsCmdLine cmd("Tutorial on gsGeometry class.");
    cmd.addPlainString("filename", "G+Smo input geometry file.", input);
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addInt( "g", "geometry", "which geometry",  geometry );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("write", "Write point data", write);
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
    if (input.empty())
    {
        if (geometry==0)
        {
            mpBspline.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); //
            mpBspline = mpBspline.uniformSplit(); // patches 0 to 3

            mpBspline.embed(3);
            mpBspline.addAutoBoundaries();
        }
        else if (geometry==1)
        {
            mpBspline.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // patch 0
            mpBspline.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // patch 1
            mpBspline.embed(3);

            gsMatrix<> ones;
            ones = gsMatrix<>::Ones(mpBspline.patch(1).coefs().rows(),1); // patch 1

            // translate second patch to the right
            mpBspline.patch(1).coefs().col(0) += ones; // patch 1



            mpBspline.addInterface(0,boundary::side::east,1,boundary::side::west);

            mpBspline.addAutoBoundaries();
        }
        else if (geometry==2)
        {
            mpBspline.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // patch 0
            mpBspline.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // patch 1
            mpBspline.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // patch 2
            mpBspline.embed(3);

            gsMatrix<> ones;
            ones = gsMatrix<>::Ones(mpBspline.patch(1).coefs().rows(),1); // patch 1

            // translate second patch up
            mpBspline.patch(1).coefs().col(1) += ones; // patch 1

            // translate third patch up and left
            ones = gsMatrix<>::Ones(mpBspline.patch(2).coefs().rows(),1); // patch 2
            mpBspline.patch(2).coefs().col(0) -= ones; // patch 2
            mpBspline.patch(2).coefs().col(1) += ones; // patch 2


            mpBspline.addInterface(0,boundary::side::north,1,boundary::side::south);
            mpBspline.addInterface(1,boundary::side::west,2,boundary::side::east);

            mpBspline.addAutoBoundaries();
        }
        else if (geometry==3)
        {
            mpBspline.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); //
            mpBspline = mpBspline.uniformSplit(); // patches 0 to 3
            mpBspline.addAutoBoundaries();
        }
        else if (geometry==4)
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
    }
    else
        gsReadFile<>(input, mpBspline);


    mpBspline.addAutoBoundaries();
    gsDebugVar(mpBspline.nInterfaces());
    gsDebugVar(mpBspline.nBoundary());
    gsDebugVar(mpBspline.nBoxes());

    gsInfo<<mpBspline<<"\n";
    gsInfo<<"Geometry has "<<mpBspline.nPatches()<<" patches\n";
    mpBspline.checkConsistency();
    gsInfo<<mpBspline.getMaxValence()<<"\n";


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

    gsDPatch<real_t> dpatch(mpBspline);
    // dpatch.sideInfo();
    // dpatch.cornerInfo();
    dpatch.initialize();
    dpatch.initializeMapper();
    dpatch.computeMatrix();

    gsSparseMatrix<> matrix;
    dpatch.matrix_into(matrix);

    if (write)
        writeToCSVfile("matrix.csv",matrix.toDense());

    // fileIO
    std::ofstream boundaries, oboundaries,interfaces, ointerfaces,vertices, overtices;
    if (write)
    {
        ointerfaces.open("interfacesOffset.csv",std::ofstream::out);
        interfaces.open("interfaces.csv",std::ofstream::out);
        oboundaries.open("boundariesOffset.csv",std::ofstream::out);
        boundaries.open("boundaries.csv",std::ofstream::out);
        overtices.open("verticesOffset.csv",std::ofstream::out);
        vertices.open("vertices.csv",std::ofstream::out);
    }

    gsMatrix<index_t> indices1, indices2, oindices1, oindices2, corners;
    gsMatrix<> coefs1,coefs2;
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

        coefs1 = mpBspline.patch(p1).coefs();
        coefs2 = mpBspline.patch(p2).coefs();

        gsInfo<<"Interface between patch "<<p1<<" (side "<<s1<<") and patch "<<p2<<" (side "<<s2<<")\n";

        basis1 = &mpBspline.basis(p1);
        basis2 = &mpBspline.basis(p2);

        mpBspline.getInterface(patchSide(p1,s1),iface);
        basis1->matchWith(iface,*basis2,indices1,indices2);

        // mapper.matchDofs(p1,indices1,p2,indices2);
        mapper.markBoundary(p1,indices1);
        mapper.markBoundary(p2,indices2);

        gsInfo<<"\tInterface indices (matching):\n";
        gsInfo<<"\t\tBoundary 1: "<<indices1.transpose()<<"\n";
        gsInfo<<"\t\tBoundary 2: "<<indices2.transpose()<<"\n";

        gsInfo<<"\tNext-to-interface indices (matching):\n";
        oindices1 = basis1->boundaryOffset(s1,1);
        oindices2 = basis2->boundaryOffset(s2,1);

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
                indices1.reverseInPlace();
            }
        }

        gsInfo<<"\t\tBoundary 1: "<<oindices1.transpose()<<"\n";
        gsInfo<<"\t\tBoundary 2: "<<oindices2.transpose()<<"\n";


        // const index_t b1 = ;
        // const index_t b2 = s2.direction();

        // COUPLE/MATCH COEFFICIENTS of INTERIOR ONLY
        //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            size_t n = indices1.size();
            // get the indices of the basis functions adjacent to the interface
            GISMO_ASSERT(indices1.cols()==1,"Columns of indices should be 1. Someting went wrong");
            indices1 = indices1.block(1,0,n-2,1);
            indices2 = indices2.block(1,0,n-2,1);
            // corners.topRows(1) = indices1.topRows(1);
            // corners.bottomRows(1) = indices1.bottomRows(1);

            // gsInfo<<"\tCorresponding coefficients (matching):\n";
            // gsInfo<<"\t\tBoundary 1: \n";
            for (index_t k=0; k!= indices1.size(); k++)
                interfaces<<std::to_string(coefs1(indices1.at(k),0))<<","<<std::to_string(coefs1(indices1.at(k),1))<<","<<std::to_string(coefs1(indices1.at(k),2))<<",\n";
                // gsInfo<<coefs1.row(indices1.at(k))<<"\n";
            // gsInfo<<"\t\tBoundary 2: \n";
            for (index_t k=0; k!= indices2.size(); k++)
                interfaces<<std::to_string(coefs2(indices2.at(k),0))<<","<<std::to_string(coefs2(indices2.at(k),1))<<","<<std::to_string(coefs2(indices2.at(k),2))<<",\n";
                // gsInfo<<coefs2.row(indices2.at(k))<<"\n";

            n = oindices1.size();
            // get the indices of the basis functions adjacent to the interface
            GISMO_ASSERT(oindices1.cols()==1,"Columns of indices should be 1. Someting went wrong");
            oindices1 = oindices1.block(1,0,n-2,1);
            oindices2 = oindices2.block(1,0,n-2,1);
            // corners.topRows(1) = indices1.topRows(1);
            // corners.bottomRows(1) = indices1.bottomRows(1);
            if (write)
            {
            // gsInfo<<"\tCorresponding coefficients (matching):\n";
            // gsInfo<<"\t\tBoundary 1: \n";
                for (index_t k=0; k!= oindices1.size(); k++)
                    ointerfaces<<std::to_string(coefs1(oindices1.at(k),0))<<","<<std::to_string(coefs1(oindices1.at(k),1))<<","<<std::to_string(coefs1(oindices1.at(k),2))<<",\n";
                    // gsInfo<<coefs1.row(indices1.at(k))<<"\n";
                // gsInfo<<"\t\tBoundary 2: \n";

                for (index_t k=0; k!= oindices2.size(); k++)
                    ointerfaces<<std::to_string(coefs2(oindices2.at(k),0))<<","<<std::to_string(coefs2(oindices2.at(k),1))<<","<<std::to_string(coefs2(oindices2.at(k),2))<<",\n";
                    // gsInfo<<coefs2.row(indices2.at(k))<<"\n";
            }
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
            coefs1 = mpBspline.patch(p1).coefs();
            s1 = bit->side();
            gsDebug<<"Boundary of patch "<<p1<<" (side "<<s1<<")\n";

            basis1 = &mpBspline.basis(p1);
            indices1 = basis1->boundary(s1);
            gsInfo<<"\tBoundary indices: "<<indices1.transpose()<<"\n";
            oindices1 = basis1->boundaryOffset(s1,1);
            gsInfo<<"\tNext-to-boundary indices: "<<oindices1.transpose()<<"\n";

            // GET THE INDICES OF THE NON-VERTEX BASIS FUNCTIONS
            // QUESTION: THIS ONLY WORKS IF THE CORNERS ARE THE ENDS OF THE INDEX VECTORS. IS THAT THE CASE????
            //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            size_t n = indices1.size();
            // get the indices of the basis functions adjacent to the interface
            GISMO_ASSERT(indices1.cols()==1,"Columns of indices should be 1. Someting went wrong");
            indices1 = indices1.block(1,0,n-2,1);

            if (write)
            {
                for (index_t k=0; k!= indices1.size(); k++)
                    boundaries<<std::to_string(coefs1(indices1.at(k),0))<<","<<std::to_string(coefs1(indices1.at(k),1))<<","<<std::to_string(coefs1(indices1.at(k),2))<<",\n";
            }

            GISMO_ASSERT(oindices1.cols()==1,"Columns of indices should be 1. Someting went wrong");
            oindices1 = oindices1.block(1,0,n-2,1);

            if (write)
            {
                for (index_t k=0; k!= oindices1.size(); k++)
                    oboundaries<<std::to_string(coefs1(oindices1.at(k),0))<<","<<std::to_string(coefs1(oindices1.at(k),1))<<","<<std::to_string(coefs1(oindices1.at(k),2))<<",\n";
            }

            // corners.topRows(1) = indices1.topRows(1);
            // corners.bottomRows(1) = indices1.bottomRows(1);
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
    std::vector<boxSide> sides;
    gsVector<index_t> overtices1, overtices2;
    for(index_t p=0; p!=mpBspline.nPatches(); p++) // for all patches
    {
        gsInfo<<"Patch "<<p<<"\n";
        for (index_t c=1; c!=5; c++) // for all corners [corners run from 1 up to 4]
        {
            if (marked[4*p+c-1])
                continue;

            coefs1 = mpBspline.patch(p).coefs();

            corner = boxCorner(c);
            gsInfo<<"\tcorner "<<corner<<"\n";

            basis1 = &mpBspline.basis(p);
            cIdx = basis1->functionAtCorner(corner);
            gsInfo<<"\t\tfunction index "<<cIdx<<"\n";

            // SURROUNDING FUNCTIONS
            corner.getContainingSides(mpBspline.dim(),sides);

            GISMO_ASSERT(sides.size()==2,"Length must be 2");
            indices1 = basis1->boundaryOffset(sides[0],0);
            indices2 = basis1->boundaryOffset(sides[1],0);
            // gsDebugVar(cIdx);
            // gsDebugVar(indices1.transpose());
            // gsDebugVar(indices2.transpose());

            // xIdx and yIdx tell at which index of
            index_t xIdx,yIdx;
            index_t N = indices1.size()-1;
            if (cIdx==indices1(0,0))
                xIdx = 0;
            else if (cIdx==indices1(N,0))
                xIdx = 1;
            else
                GISMO_ERROR("index is not on the end points?");

            if (cIdx==indices2(0,0))
                yIdx = 0;
            else if (cIdx==indices2(N,0))
                yIdx = 1;
            else
                GISMO_ERROR("index is not on the end points?");

            // gsDebugVar(xIdx);
            // gsDebugVar(yIdx);



            overtices1 = basis1->boundaryOffset(sides[0],1);
            overtices2 = basis1->boundaryOffset(sides[1],1);


            overtices1 = xIdx==1 ? overtices1.tail(1) : overtices1.head(1);
            overtices2 = yIdx==1 ? overtices2.tail(2) : overtices2.head(2);

            // gsDebugVar(overtices1);
            // gsDebugVar(overtices2);

            overtices1 = basis1->boundaryOffset(sides[0],2);
            overtices2 = basis1->boundaryOffset(sides[1],2);

            overtices1 = xIdx==1 ? overtices1.tail(2) : overtices1.head(2);
            overtices2 = yIdx==1 ? overtices2.tail(3) : overtices2.head(3);

            // gsDebugVar(overtices1);
            // gsDebugVar(overtices2);


            // bool minX = (min(oindices1.array()) < cIdx).any(); // does there exist an index in oindices1 that is smaller than cIdx?
            // bool maxX = (min(oindices1.array()) < cIdx).any(); // does there exist an index in oindices1 that is smaller than cIdx?
            // bool minY = (min(oindices2.array()) > cIdx).any(); // does there exist an index in oindices2 that is smaller than cIdx?
            // bool maxY = (min(oindices2.array()) > cIdx).any(); // does there exist an index in oindices2 that is smaller than cIdx?
            // gsDebugVar(minX);
            // gsDebugVar(maxX);
            // gsDebugVar(minY);
            // gsDebugVar(maxY);

            // overtices1.resize(3,1);
            // if (minX && minY) // south west
            // {
            //     overtices1.topRows(1)=oindices1.topRows(1);
            //     overtices1.bottomRows(2)=oindices2.topRows(2);
            // }
            // else if (minX && !minY) // south east
            // {
            //     overtices1.topRows(1)=oindices1.topRows(1);
            //     overtices1.bottomRows(2)=oindices2.bottomRows(2);
            // }
            // else if (X && !Y)
            // {
            //     overtices1.topRows(1)=oindices1.bottomRows(1);
            //     overtices1.bottomRows(2)=oindices2.topRows(2);
            // }
            // else if (X && Y)
            // {
            //     overtices1.topRows(1)=oindices1.bottomRows(1);
            //     overtices1.bottomRows(2)=oindices2.bottomRows(2);
            // }

            // gsDebugVar(overtices1);

            std::vector<patchCorner> cornerList;
            patchCorner pcorner = patchCorner(p,corner);
            bool isCycle = mpBspline.getCornerList(pcorner,cornerList);

            gsInfo<<"\t\tvalence "<<cornerList.size()<<"\n";
            gsInfo<<"\t\ttype:"<<( isCycle ? " interior vertex " : " boundary vertex " )<<"\n";

            if (write)
                vertices<<std::to_string(coefs1(cIdx,0))<<","<<std::to_string(coefs1(cIdx,1))<<","<<std::to_string(coefs1(cIdx,2))<<",\n";


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

    // gsDebugVar(mpBspline.basis(0).component(0));
    // gsDebugVar(mpBspline.basis(0).component(0).size());
    // gsWriteParaview(mpBspline.basis(0),"basis");

    mapper.finalize();
    // gsDebugVar(mapper.allFree());
    // gsDebugVar(mapper.size());
    // gsDebugVar(mapper.freeSize());
    // gsDebugVar(mapper.coupledSize());
    // gsDebugVar(mapper.index(0,1));
    // gsDebugVar(mapper.index(1,1));
    // gsDebugVar(mapper.index(2,1));
    // gsDebugVar(mapper.index(3,1));
    // gsDebugVar(mapper.index(4,1));
    // gsDebugVar(mapper.index(5,1));



    if (write)
    {

        interfaces.close();
        ointerfaces.close();
        boundaries.close();
        oboundaries.close();
        vertices.close();
        overtices.close();
    }

    gsDebugVar(dpatch.indexFromSides(1,patchSide(0,1),1,patchSide(0,4)));
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


    // std::vector<index_t> boxes(10);
    // boxes[0] = 1;
    // boxes[1] = boxes[2] = 0;
    // boxes[3] = boxes[4] = 2;

    // boxes[5] = 1;
    // boxes[6] = boxes[7] = 2;
    // boxes[8] = boxes[9] = 4;
    // mp.patch(0).refineElements(boxes);

    // gsDebugVar(mp.patch(0).coefs());

    // gsWriteParaview(mp.patch(0),"mp2",1000,true);
    // gsWriteParaview(mp.basis(0),"basis2",1000);

    // gsDebugVar(mp.basis(0).boundary(boxSide(2)));



    return 0;
}


