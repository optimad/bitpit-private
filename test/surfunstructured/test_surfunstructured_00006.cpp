// ========================================================================== //
//           ** BitPit mesh ** Test 006 for class SurfUnstructured **         //
//                                                                            //
// ID renumbering of a  class SurfUnstructured.                               //
// ========================================================================== //
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

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //

// Standard Template Library
# include <array>
# include <vector>
# include <iostream>

// BitPit
# include "bitpit_common.hpp"                                                 // Utilities and common definitions
# include "bitpit_IO.hpp"                                                     // Input/output
# include "bitpit_operators.hpp"                                              // STL containers operators
# include "bitpit_patchkernel.hpp"                                            // BitPit base patch
# include "bitpit_surfunstructured.hpp"                                           // BitPit surftri patch

// ========================================================================== //
// NAMESPACES                                                                 //
// ========================================================================== //
using namespace std;
using namespace bitpit;

// ========================================================================== //
// GENERATE A TEST NON-MANIFOLD SURFACE TRIANGULATION FOR TESTS.              //
// ========================================================================== //
void generateTestTriangulation(
    SurfUnstructured                &mesh
) {

// ========================================================================== //
// void generateTestTriangulation(                                            //
//     SurfUnstructured            &mesh)                                     //
//                                                                            //
// Generate a non-manifold surface triangulation for tests.                   //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - mesh    : SurfUnstructured, surface mesh patch                           //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
long                    nV = 27;
long                    nS = 33;

// Counters
int                     i;

// ========================================================================== //
// INITIALIZE TRIANGULATION                                                   //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //

    // Reserve memory for vertex & cell storage ----------------------------- //
    mesh.reserveVertices(nV);
    mesh.reserveCells(nS);
}

// ========================================================================== //
// GENERATE VERTEX LIST                                                       //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    double                      off = -0.5;
    array<double, 3>            vertex;
    
    // 0-row ---------------------------------------------------------------- //
    vertex.fill(0.0);
    for (i = 0; i < 8; ++i) {
        vertex[0] = double(i);
        mesh.addVertex(vertex);
    } //next i

    // 1-row ---------------------------------------------------------------- //
    vertex.fill(0.0);
    vertex[1] = 1.0;
    for (i = 0; i < 9; ++i) {
        vertex[0] = double(i) + off;
        mesh.addVertex(vertex);
    } //next i

    // 2-row ---------------------------------------------------------------- //
    vertex.fill(0.0);
    vertex[1] = 2.0;
    for (i = 0; i < 8; ++i) {
        vertex[0] = double(i);
        mesh.addVertex(vertex);
    } //next i

    // Orthogonal element(s) ------------------------------------------------ //
    vertex[0] = 0.5*(mesh.getVertex(3)[0] + mesh.getVertex(12)[0]);
    vertex[1] = 0.5*(mesh.getVertex(3)[1] + mesh.getVertex(12)[1]);
    vertex[2] = 0.5 * sqrt(3.0);
    mesh.addVertex(vertex);
    vertex[0] = 0.5*(mesh.getVertex(12)[0] + mesh.getVertex(21)[0]);
    vertex[1] = 0.5*(mesh.getVertex(12)[1] + mesh.getVertex(21)[1]);
    vertex[2] = 0.5 * sqrt(3.0);
    mesh.addVertex(vertex);
}

// ========================================================================== //
// GENERATE CONNECTIVITY                                                      //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    long int                    off;
    vector<long>                connectivity(3);

    // 0-row ---------------------------------------------------------------- //
    off = 8;
    for (i = 0; i < 7; ++i) {
        connectivity[0] = i;
        connectivity[1] = i + 1 + off;
        connectivity[2] = i + off;
        mesh.addCell(ElementInfo::TRIANGLE, true, connectivity);
        connectivity[0] = i;
        connectivity[1] = i + 1;
        connectivity[2] = i + 1 + off;
        mesh.addCell(ElementInfo::TRIANGLE, true, connectivity);
    } //next i
    connectivity[0] = i;
    connectivity[1] = i + 1 + off;
    connectivity[2] = i + off;
    mesh.addCell(ElementInfo::TRIANGLE, true, connectivity);

    // 1-row ---------------------------------------------------------------- //
    off = 9;
    for (i = 8; i < 15; ++i) {
        connectivity[0] = i;
        connectivity[1] = i + 1;
        connectivity[2] = i + off;
        mesh.addCell(ElementInfo::TRIANGLE, true, connectivity);
        connectivity[0] = i + 1;
        connectivity[1] = i + 1 + off;
        connectivity[2] = i + off;
        mesh.addCell(ElementInfo::TRIANGLE, true, connectivity);
    } //next i
    connectivity[0] = i;
    connectivity[1] = i + 1;
    connectivity[2] = i + off;
    mesh.addCell(ElementInfo::TRIANGLE, true, connectivity);

    // Orthogonal element --------------------------------------------------- //
    connectivity[0] = 3;
    connectivity[1] = 12;
    connectivity[2] = 25;
    mesh.addCell(ElementInfo::TRIANGLE, true, connectivity);
    connectivity[0] = 12;
    connectivity[1] = 26;
    connectivity[2] = 25;
    mesh.addCell(ElementInfo::TRIANGLE, true, connectivity);
    connectivity[0] = 12;
    connectivity[1] = 21;
    connectivity[2] = 26;
    mesh.addCell(ElementInfo::TRIANGLE, true, connectivity);
    
}

return;
}

// ========================================================================== //
// SUBTEST #001 Creating a mesh a renumber its ids.		                      //
//	returning 0 if successfull, 1 if failed.	
// ========================================================================== //
int subtest_001() 
{

// Local variables
	SurfUnstructured                        mesh(0);
	mesh.setExpert(true);

//	generate a Dummy Triangulation
    generateTestTriangulation(mesh);
    mesh.buildAdjacencies();
	mesh.buildInterfaces();
	
	int nVert, nCell, nInterf;
	nVert	= mesh.getVertexCount();
	nCell	= mesh.getCellCount();
	nInterf	= mesh.getInterfaceCount();
	
	long offV = 13;
	long offC = 21;
	long offI = 1012;
	long counter;
	//creating checkVectors
	std::vector<long>	checkV(nVert), checkC(nCell), checkI(nInterf);
	
	counter = offV;
	for(auto & val : checkV){
		val = counter;
		++counter;
	}
	
	counter = offC;
	for(auto & val : checkC){
		val = counter;
		++counter;
	}
	
	counter = offI;
	for(auto & val : checkI){
		val = counter;
		++counter;
	}
	
	//renumbering mesh
	mesh.renumberPatch(offV,offC,offI);
	
 	std::vector<long> newVertexID = mesh.getVertices().getIds(true);
	std::vector<long> newCellID = mesh.getCells().getIds(true);
	std::vector<long> newInterfaceID = mesh.getInterfaces().getIds(true);
	
	std::cout<<newVertexID<<std::endl;
	std::cout<<newCellID<<std::endl;
	std::cout<<newInterfaceID<<std::endl;
	
	bool check;
	
	//check renumbering vertices
	check = (newVertexID.size() == nVert);
	counter=0;
	for(auto & id:newVertexID){
		check = check && (id == checkV[counter]);
		++counter;
	}
	if(check)	return 1;
	
	//check renumbering cells
	check = (newCellID.size() == nCell);
	counter=0;
	for(auto & id:newCellID){
		check = check && (id == checkC[counter]);
		++counter;
	}
	if(check)	return 1;
	
	//check renumbering interfaces
	check = (newInterfaceID.size() == nInterf);
	counter=0;
	for(auto & id:newInterfaceID){
		check = check && (id == checkI[counter]);
		++counter;
	}
	if(check)	return 1;
	
	
	mesh.write("renumberedMesh");
	return 0;
}

// ========================================================================== //
// MAIN FOR TEST #00006                                                       //
// ========================================================================== //
int main(
    void
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variabels
int                             err = 0;

// ========================================================================== //
// RUN SUB-TEST #001                                                          //
// ========================================================================== //
err = subtest_001();
if (err > 0) return(10 + err);

return err;

}
