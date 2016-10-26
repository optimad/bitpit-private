/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

# include <cassert>
# include <memory>

# include "bitpit_common.hpp"
# include "bitpit_operators.hpp"
# include "bitpit_surfunstructured.hpp"

# include "levelSetKernel.hpp"
# include "levelSetCartesian.hpp"
# include "levelSetOctree.hpp"

# include "levelSetObject.hpp"
# include "levelSetCachedObject.hpp"
# include "levelSetSegmentation.hpp"
# include "levelSetMask.hpp"

namespace bitpit {

/*!
	@ingroup levelset
	@class  LevelSetMask
	@brief Implements the levelset around a set of cells or interfaces of the kernel     
*/

/*!
 * Destructor
 */
LevelSetMask::~LevelSetMask() {
};

/*!
 * Constructor
 * @param[in] id identifier of object
 * @param[in] mask the list of inner cells
 * @param[in] mesh the mesh hosting the cells
 */
LevelSetMask::LevelSetMask(int id, const std::unordered_set<long> &mask, const VolumeKernel &mesh ) :LevelSetSegmentation(id) {

    std::unordered_map<long,long> envelopeToMesh ;
    std::unordered_set<long> flipList ;

    SurfUnstructured *envelope = new SurfUnstructured(extractCellEnvelope(mask,mesh,&envelopeToMesh)) ;

    orientByMesh(envelopeToMesh, mesh, *envelope, flipList); 
    orientByMask(envelopeToMesh, mesh, mask, flipList); 
    flip(*envelope, flipList);

    std::unique_ptr<SurfUnstructured> segmentation(envelope) ;
    setSegmentation(std::move(segmentation));
};

/*!
 * Constructor
 * @param[in] id identifier of object
 * @param[in] list the list of interfaces
 * @param[in] mesh the mesh hosting the cells
 */
LevelSetMask::LevelSetMask(int id, const std::vector<long> &list, const VolumeKernel &mesh ) :LevelSetSegmentation(id) {

    std::unordered_map<long,long> envelopeToMesh ;
    std::unordered_set<long> flipList ;

    SurfUnstructured *envelope = new SurfUnstructured(extractFaceEnvelope(list,mesh,&envelopeToMesh)) ;

    orientByMesh(envelopeToMesh, mesh, *envelope, flipList); 
    flip(*envelope,flipList);

    std::unique_ptr<SurfUnstructured> segmentation(envelope) ;
    setSegmentation(std::move(segmentation));
};

/*!
 * Constructor
 * @param[in] id identifier of object
 * @param[in] list the list of interfaces
 * @param[in] mark list of interfaces to be flipped
 * @param[in] mesh the mesh hosting the cells
 */
LevelSetMask::LevelSetMask(int id, const std::vector<long> &list, const std::unordered_set<long> &mark, const VolumeKernel &mesh ) :LevelSetSegmentation(id) {

    std::unordered_map<long,long> envelopeToMesh ;
    std::unordered_set<long> flipList ;

    SurfUnstructured *envelope = new SurfUnstructured(extractFaceEnvelope(list,mesh,&envelopeToMesh)) ;

    orientByMesh(envelopeToMesh, mesh, *envelope, flipList); 
    orientByMark(envelopeToMesh, mark, flipList); 
    flip(*envelope,flipList);

    std::unique_ptr<SurfUnstructured> segmentation(envelope) ;
    setSegmentation(std::move(segmentation));
};

/*!
 * Extracts the external envelope and create a new patch from it
 * The external envelope is composed by all the outer faces of the masked cells
 * @param[in] mask the list of inner cells
 * @param[in] mesh the mesh hosting the cells
 * @param[in] envelopeToMesh pointer to map which hosts the ndex association between the cells of the envelope and the faces of the mesh.
 * If the nullptr is passed, a local map will be used.
 * @return surface mesh
*/
SurfUnstructured LevelSetMask::extractCellEnvelope(const std::unordered_set<long> &mask, const VolumeKernel &mesh, std::unordered_map<long,long> *envelopeToMesh){

    std::vector<long> list;

    for( const long &cellIndex : mask){
        const auto &cell = mesh.getCell(cellIndex);
        const long *adjacencies = cell.getAdjacencies();
        const long *interfaceIndex = cell.getInterfaces();

        int interfaceCount= cell.getInterfaceCount() ;

        for(int i=0; i<interfaceCount; i++){
            const long &neigh = adjacencies[i];

            if(mask.count(neigh)==0){
                list.push_back(interfaceIndex[i]);
            }
        }
    }

    return extractFaceEnvelope(list,mesh,envelopeToMesh);
}

/*!
 * Extracts the external envelope and create a new patch from it.
 * The external envelope is composed by all interfaces in the list.
 * @param[in] list the list of interfaces
 * @param[in] mesh the mesh hosting the cells
 * @param[in] envelopeToMesh pointer to map which hosts the index association between the cells of the envelope and the faces of the mesh.
 * If the nullptr is passed, a local map will be used.
 * @return surface mesh
*/
SurfUnstructured LevelSetMask::extractFaceEnvelope(const std::vector<long> &list, const VolumeKernel &mesh, std::unordered_map<long,long> *envelopeToMesh){

    std::unordered_map<long,long> localEnvelopeToMesh;

    if(envelopeToMesh==nullptr){
        envelopeToMesh = &localEnvelopeToMesh;
    }

    SurfUnstructured envelope(0,mesh.getDimension()-1,mesh.getDimension());

	// ====================================================================== //
	// RESIZE DATA STRUCTURES                                                 //
	// ====================================================================== //
    long nVertices(0);
    long nCells(list.size());

    for( const long &faceIndex : list){
        auto const &interface = mesh.getInterface(faceIndex);
        nVertices += interface.getVertexCount() ;
    }

	envelope.reserveVertices(nVertices);
	envelope.reserveCells(nCells);

	// ====================================================================== //
	// LOOP OVER CELLS                                                        //
	// ====================================================================== //
	std::unordered_map<long,long> vertexMap;

    for( const long &faceIndex : list){
        auto const &interface = mesh.getInterface(faceIndex);
        const long *faceConnect = interface.getConnect();
        int nFaceVertices = interface.getVertexCount();

        // Add face vertices to the envelope and get face
        // connectivity in the envelope
        std::unique_ptr<long[]> faceEnvelopeConnect = std::unique_ptr<long[]>(new long[nFaceVertices]);
        for (int j = 0; j < nFaceVertices; ++j) {
        	long vertexId = faceConnect[j];
        
        	// If the vertex is not yet in the envelope
        	// add it.
        	if (vertexMap.count(vertexId) == 0) {
        		const Vertex &vertex = mesh.getVertex(vertexId);
        		auto envelopeVertex = envelope.addVertex(vertex);
        		vertexMap[vertexId] = envelopeVertex->getId();
        	}
        
        	// Update face connectivity in the envelope
        	faceEnvelopeConnect[j] = vertexMap.at(vertexId);
        }

        // Add face to envelope
        ElementInfo::Type faceType = interface.getType();
        PatchKernel::CellIterator cellItr = envelope.addCell(faceType, true, std::move(faceEnvelopeConnect));
        envelopeToMesh->insert({{cellItr->getId(),faceIndex}});
	}

    envelope.squeeze();
    envelope.buildAdjacencies();

    envelope.getVTK().setName("geometry_002") ;
    envelope.write() ;

    return envelope;

}

/*!
 * Orients the face normals according to a mask of interior cells.
 * If a face owner is not within this list the normal will be flipped.
 * @param[in] envelope surface mesh
 * @param[in] mesh the mesh hosting the cells
 * @param[in] envelopeToMesh index association between the cells of the envelope and the faces of the mesh.
 * @param[in,out] flipList the list of faces to flip
*/
void LevelSetMask::orientByMesh(const std::unordered_map<long,long> &envelopeToMesh, const VolumeKernel &mesh, SurfUnstructured &envelope, std::unordered_set<long> &flipList){

    for( auto element : envelopeToMesh){
        long enveIndex = element.first ;
        long faceIndex = element.second ;

        std::array<double,3> facetNormal = envelope.evalFacetNormal(enveIndex);
        std::array<double,3> interfaceNormal = mesh.evalInterfaceNormal(faceIndex);

        bool needFlip = (dotProduct(facetNormal,interfaceNormal)<0) ;
        updateFlipList(flipList,enveIndex,needFlip) ;

    }
}

/*!
 * Orients the face normals according to a mask of interior cells.
 * If a face owner is not within this list the normal will be flipped.
 * @param[in] mask the list of interior cells
 * @param[in] mesh the mesh hosting the cells
 * @param[in] envelopeToMesh pointer to map which hosts the index association between the cells of the envelope and the faces of the mesh.
 * @param[in,out] flipList the list of faces to flip
*/
void LevelSetMask::orientByMask(const std::unordered_map<long,long> &envelopeToMesh, const VolumeKernel &mesh, const std::unordered_set<long> &mask, std::unordered_set<long> &flipList){

    for( auto element : envelopeToMesh){
        long enveIndex = element.first ;
        long faceIndex = element.second ;

        auto const &interface = mesh.getInterface(faceIndex);
        long ownerId=interface.getOwner();

        bool needFlip = (mask.count(ownerId)==0) ;
        updateFlipList( flipList, enveIndex, needFlip);

    }
}

/*!
 * Orients the face normals according to a marked faces.
 * If a face is within this list the normal will be flipped.
 * @param[in] envelopeToMesh pointer to map which hosts the index association between the cells of the envelope and the faces of the mesh.
 * @param[in] mark the list of faces to flip
 * @param[in,out] flipList the list of faces to flip
*/
void LevelSetMask::orientByMark( const std::unordered_map<long,long> &envelopeToMesh, const std::unordered_set<long> &mark, std::unordered_set<long> &flipList){

    for( auto element : envelopeToMesh){
        long enveIndex = element.first ;
        long faceIndex = element.second ;

        bool needFlip = (mark.count(faceIndex)!=0);
        updateFlipList( flipList, enveIndex, needFlip);

    }
}

/*!
 * Flips the normal of list of cells
 * @param[in] envelope surface mesh
 * @param[in] flipList the list of faces to flip
*/
void LevelSetMask::flip(SurfUnstructured &envelope, const std::unordered_set<long> &flipList){

    for( long index : flipList){

        auto &face = envelope.getCell(index); 
        const long *faceConnect = face.getConnect();
        int nFaceVertices = face.getVertexCount();
        std::unique_ptr<long[]> faceEnvelopeConnect = std::unique_ptr<long[]>(new long[nFaceVertices]);

        for (int j = 0; j < nFaceVertices; ++j) {
            faceEnvelopeConnect[j] = faceConnect[j] ;
        }

        for (int j = 0; j < (int) std::ceil(nFaceVertices)/2; ++j) {
            int top = j;
            int end = nFaceVertices-1-j;
            std::swap( faceEnvelopeConnect[top], faceEnvelopeConnect[end]);
        }

        face.setConnect(std::move(faceEnvelopeConnect)) ;
    }
}

/*!
 * Updates the list of faces which need to be flipped.
 * If a face is already within the set to be flipped and it should be flipped again, it will be removed from the set.
 * On the contary, if it was not present in the set, it will be added.
 * The set remains unaltered if the face does not need to be flipped.
 *
 * @param[in,out] flipList the list of faces which need to be flipped
 * @param[in] index the index of the face to be updated
 * @param[in] flip indicator if it needs to be flipped
 */
void LevelSetMask::updateFlipList( std::unordered_set<long> &flipList, const long &index, const bool &flip){

    if(flip){
        auto flipItr = flipList.find(index) ;
        if( flipItr==flipList.end()){
            flipList.insert(index);
        } else {
            flipList.erase(flipItr) ;
        }
    }
}

}

