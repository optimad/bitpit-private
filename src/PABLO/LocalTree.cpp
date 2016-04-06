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

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "bitpit_common.hpp"

#include "LocalTree.hpp"
#include <map>
#include <unordered_map>

namespace bitpit {

// =================================================================================== //
// NAME SPACES                                                                         //
// =================================================================================== //
using namespace std;

// =================================================================================== //
// CLASS IMPLEMENTATION                                                                    //
// =================================================================================== //

// =================================================================================== //
// CONSTRUCTORS AND OPERATORS
// =================================================================================== //

/*!Dimensional and default constructor.
 * \param[in] maxlevel Maximum refinement level of the octree.
 * \param[in] dim Space dimension of octree.
 */
LocalTree::LocalTree(int8_t maxlevel, uint8_t dim):m_firstDesc(dim, maxlevel),m_lastDesc(dim, maxlevel){
	m_dim = dim;
	m_global.setGlobal(maxlevel, m_dim);
	Octant oct0(m_dim, m_global.m_maxLevel);
	Octant octf(m_dim,m_global.m_maxLevel,0,0,0, m_global.m_maxLevel);
	Octant octl(m_dim,m_global.m_maxLevel,m_global.m_maxLength-1,m_global.m_maxLength-1,(m_dim-2)*(m_global.m_maxLength-1), m_global.m_maxLevel);
	m_octants.clear();
	m_octants.push_back(oct0);
	m_firstDesc = octf;
	m_lastDesc = octl;
	m_ghosts.clear();
	m_sizeGhosts = 0;
	m_localMaxDepth = 0;
	m_balanceCodim = 1;
	m_periodic.resize(m_dim*2);
};

/*!Default destructor.
 */
LocalTree::~LocalTree(){};

// =================================================================================== //
// METHODS
// =================================================================================== //

// =================================================================================== //
// BASIC GET/SET METHODS
// =================================================================================== //

/*!Get the first descentant octant of the octree.
 * \return Constant reference to the first descendant of the octree.
 */
const Octant&
LocalTree::getFirstDesc() const{
	return m_firstDesc;
};

/*!Get the last descentant octant of the octree.
 * \return Constant reference to the last descendant of the octree.
 */
const Octant&
LocalTree::getLastDesc() const{
	return m_lastDesc;
};

/*! Get the number of the ghosts for the local partition of the tree.
 * \return Number of ghosts.
 */
uint32_t
LocalTree::getSizeGhost() const{
	return m_sizeGhosts;
};

/*! Get the number of the octants in the local tree.
 * \return Number of local octants.
 */
uint32_t
LocalTree::getNumOctants() const{
	return m_octants.size();
};

/*! Get max depth reached in local tree
 * \return Max depth in local partition of the octree.
 */
uint8_t
LocalTree::getLocalMaxDepth() const{
	return m_localMaxDepth;
};

/** Get refinement/coarsening marker for idx-th octant
 * \param[in] idx Local index of the target octant.
 * \return Marker of the octant.
 */
int8_t
LocalTree::getMarker(int32_t idx){
	return m_octants[idx].getMarker();
};

/** Get refinement/coarsening marker for idx-th octant
 * \param[in] idx Local index of the target octant.
 * \return Level of the octant.
 */
uint8_t
LocalTree::getLevel(int32_t idx){
	return m_octants[idx].getLevel();
};

/** Compute the Morton index of the idx-th octant (without level).
 * \param[in] idx Local index of the target octant.
 * \return Morton index of the octant.
 */
uint64_t
LocalTree::computeMorton(int32_t idx){
	return m_octants[idx].computeMorton();
};

/** Get refinement/coarsening marker for idx-th ghost octant
 * \param[in] idx Local index of the target ghost octant.
 * \return Level of the ghost octant.
 */
uint8_t
LocalTree::getGhostLevel(int32_t idx){
	return m_ghosts[idx].getLevel();
};

/** Compute the Morton index of the idx-th ghost octant (without level).
 * \param[in] idx Local index of the target octant.
 * \return Morton index of the octant.
 */
uint64_t
LocalTree::computeGhostMorton(int32_t idx){
	return m_ghosts[idx].computeMorton();
};

/** Get if balancing-blocked idx-th octant
 * \param[in] idx Local index of the target octant.
 * \return Has the octant to be balanced?
 */
bool
LocalTree::getBalance(int32_t idx){
	return m_octants[idx].getBalance();
};

/*! Get the codimension for 2:1 balancing
 * \return Maximum codimension of the entity through which the 2:1 balance is performed.
 */
uint8_t
LocalTree::getBalanceCodim() const{
	return m_balanceCodim;
};

/** Set refinement/coarsening marker for idx-th octant
 * \param[in] idx Local index of the target octant.
 * \param[in] marker Refinement marker for the target octant.
 */
void
LocalTree::setMarker(int32_t idx, int8_t marker){
	m_octants[idx].setMarker(marker);
};

/** Set if balancing-blocked idx-th octant
 * \param[in] idx Local index of the target octant.
 * \param[in] balance Has the octant to be balanced?
 */
void
LocalTree::setBalance(int32_t idx, bool balance){
	m_octants[idx].setBalance(balance);
};

/*! Set the codimension for 2:1 balancing
 * \param[in] Maximum codimension of the entity through which the 2:1 balance is performed.
 */
void
LocalTree::setBalanceCodim(uint8_t b21codim){
	m_balanceCodim = b21codim;
};

/*!Set the first descentant octant of the octree.
 */
void
LocalTree::setFirstDesc(){
	octvector::const_iterator firstOctant = m_octants.begin();
	m_firstDesc = Octant(m_dim, m_global.m_maxLevel, firstOctant->m_x, firstOctant->m_y, firstOctant->m_z, m_global.m_maxLevel);
};

/*!Set the last descentant octant of the octree.
 */
void
LocalTree::setLastDesc(){
	octvector::const_iterator lastOctant = m_octants.end() - 1;
	uint32_t x,y,z,delta;
	delta = (uint32_t)(1<<((uint8_t)m_global.m_maxLevel - lastOctant->m_level)) - 1;
	x = lastOctant->m_x + delta;
	y = lastOctant->m_y + delta;
	z = lastOctant->m_z + (m_dim-2)*delta;
	m_lastDesc = Octant(m_dim, m_global.m_maxLevel,x,y,z, m_global.m_maxLevel);
};

/*! Set the periodic condition of the boundaries.
 * \param[in] periodic Vector with the periodic conditions (true/false) of each boundary.
 */
void
LocalTree::setPeriodic(bvector & periodic){
	m_periodic = periodic;
};

// =================================================================================== //
// OTHER GET/SET METHODS
// =================================================================================== //

// =================================================================================== //
// OTHER METHODS
// =================================================================================== //

/*!Extract an octant of the octree.
 * \param[in] idx Local index of the target octant.
 * \return Reference to the idx-th octant of the octree.
 */
Octant&
LocalTree::extractOctant(uint32_t idx){
	return m_octants[idx];
};

/*!Extract an octant of the octree.
 * \param[in] idx Local index of the target octant.
 * \return Constant reference to the idx-th octant of the octree.
 */
const Octant&
LocalTree::extractOctant(uint32_t idx) const{
	return m_octants[idx];
};

/*!Extract a ghost octant of the octree.
 * \param[in] idx Local index of the target ghost octant.
 * \return Reference to the idx-th ghost octant of the octree.
 */
Octant&
LocalTree::extractGhostOctant(uint32_t idx) {
	return m_ghosts[idx];
};

/*!Extract a ghost octant of the octree.
 * \param[in] idx Local index of the target ghost octant.
 * \return Constant reference to the idx-th ghost octant of the octree.
 */
const Octant&
LocalTree::extractGhostOctant(uint32_t idx) const{
	return m_ghosts[idx];
};

// =================================================================================== //

/*! Refine local tree: refine one time octants with marker >0
 * \param[out] mapidx mapidx[i] = index in old octants vector of the new i-th octant (index of father if octant is new after refinement)
 * \return	true if refinement done
 */
bool
LocalTree::refine(u32vector & mapidx){

	u32vector		last_child_index;
	octvector 		children(0,Octant(m_dim, m_global.m_maxLevel));
	uint32_t 		idx, nocts, ilastch;
	uint32_t 		offset = 0, blockidx;
	uint32_t		mapsize = mapidx.size();
	uint8_t 		nchm1 = m_global.m_nchildren-1, ich;
	bool 			dorefine = false;

	nocts = m_octants.size();
	for (idx=0; idx<nocts; idx++){
		if(m_octants[idx].getMarker() > 0 && m_octants[idx].getLevel() < m_global.m_maxLevel){
			last_child_index.push_back(idx+nchm1+offset);
			offset += nchm1;
		}
		else{
			if (m_octants[idx].m_marker > 0){
				m_octants[idx].m_marker = 0;
				m_octants[idx].m_info[15] = false;
			}
		}
	}
	if (offset > 0){
		if(mapsize > 0){
			mapidx.resize(m_octants.size()+offset);
		}
		m_octants.resize(m_octants.size()+offset, Octant(m_dim, m_global.m_maxLevel));
		blockidx = last_child_index[0]-nchm1;
		idx = m_octants.size();
		ilastch = last_child_index.size()-1;

		while (idx>blockidx){
			idx--;
			if(idx == last_child_index[ilastch]){
				children = m_octants[idx-offset].buildChildren();
				for (ich=0; ich<m_global.m_nchildren; ich++){
					m_octants[idx-ich] = (children[nchm1-ich]);
					if(mapsize>0) mapidx[idx-ich]  = mapidx[idx-offset];
				}
				offset -= nchm1;
				idx -= nchm1;
				//Update local max depth
				if (children[0].getLevel() > m_localMaxDepth){
					m_localMaxDepth = children[0].getLevel();
				}
				if (children[0].getMarker() > 0 && mapsize == 0){
					//More Refinement to do
					dorefine = true;
				}
				//delete []children;
				if (ilastch != 0){
					ilastch--;
				}
			}
			else {
				m_octants[idx] = m_octants[idx-offset];
				if(mapsize>0) mapidx[idx]  = mapidx[idx-offset];
			}
		}
	}

	octvector(m_octants).swap(m_octants);
	nocts = m_octants.size();

	setFirstDesc();
	setLastDesc();

	return dorefine;

};

// =================================================================================== //
/*! Coarse local tree: coarse one time family of octants with marker <0
 * (if at least one octant of family has marker>=0 set marker=0 for the entire family)
 * \param[out] mapidx mpaidx[i] = index in old octants vector of the new i-th octant (index of first child if octant is new after coarsening)
 * \return	true is coarsening done
 */
bool
LocalTree::coarse(u32vector & mapidx){

	u32vector		first_child_index;
	Octant			father(m_dim, m_global.m_maxLevel);
	uint32_t 		nocts, nocts0;
	uint32_t 		idx, idx2;
	uint32_t 		offset;
	uint32_t 		idx2_gh;
	uint32_t 		nidx;
	uint32_t		mapsize = mapidx.size();
	int8_t 			markerfather, marker;
	uint8_t 		nbro, nend;
	uint8_t 		nchm1 = m_global.m_nchildren-1;
	bool 			docoarse = false;
	bool 			wstop = false;

	//------------------------------------------ //
	// Initialization

	nbro = nend = 0;
	nidx = offset = 0;

	idx2_gh = 0;

	nocts = nocts0 = m_octants.size();
	m_sizeGhosts = m_ghosts.size();


	// Init first and last desc (even if already calculated)
	setFirstDesc();
	setLastDesc();

	//------------------------------------------ //

	// Set index for start and end check for ghosts
	if (m_ghosts.size()){
		while(idx2_gh < m_sizeGhosts && m_ghosts[idx2_gh].computeMorton() < m_lastDesc.computeMorton()){
			idx2_gh++;
		}
		idx2_gh = min((m_sizeGhosts-1), idx2_gh);
	}

	// Check and coarse internal octants
	for (idx=0; idx<nocts; idx++){
		if(m_octants[idx].getMarker() < 0 && m_octants[idx].getLevel() > 0){
			nbro = 0;
			father = m_octants[idx].buildFather();
			// Check if family is to be refined
			for (idx2=idx; idx2<idx+m_global.m_nchildren; idx2++){
				if (idx2<nocts){
					if(m_octants[idx2].getMarker() < 0 && m_octants[idx2].buildFather() == father){
						nbro++;
					}
				}
			}
			if (nbro == m_global.m_nchildren){
				nidx++;
				first_child_index.push_back(idx);
				idx = idx2-1;
			}
		}
	}
	uint32_t nblock = nocts;
	uint32_t nfchild = first_child_index.size();
	if (nidx!=0){
		nblock = nocts - nidx*nchm1;
		nidx = 0;
		for (idx=0; idx<nblock; idx++){
			if (idx+offset < nocts){
				if (nidx < nfchild){
					if (idx+offset == first_child_index[nidx]){
						markerfather = -m_global.m_maxLevel;
						father = m_octants[idx+offset].buildFather();
						for (uint32_t iii=0; iii<17; iii++){
							father.m_info[iii] = false;
						}
						for(idx2=0; idx2<m_global.m_nchildren; idx2++){
							if (idx2 < nocts){
								if (markerfather < m_octants[idx+offset+idx2].getMarker()+1){
									markerfather = m_octants[idx+offset+idx2].getMarker()+1;
								}
								for (uint32_t iii=0; iii<17; iii++){
									father.m_info[iii] = father.m_info[iii] || m_octants[idx+offset+idx2].m_info[iii];
								}
							}
						}
						father.m_info[13] = true;
						father.m_info[15] = true;
						father.setMarker(markerfather);
						if (markerfather < 0 && mapsize == 0){
							docoarse = true;
						}
						m_octants[idx] = father;
						if(mapsize > 0) mapidx[idx] = mapidx[idx+offset];
						offset += nchm1;
						nidx++;
					}
					else{
						m_octants[idx] = m_octants[idx+offset];
						if(mapsize > 0) mapidx[idx] = mapidx[idx+offset];
					}
				}
				else{
					m_octants[idx] = m_octants[idx+offset];
					if(mapsize > 0) mapidx[idx] = mapidx[idx+offset];
				}
			}
		}
	}
	m_octants.resize(nblock, Octant(m_dim, m_global.m_maxLevel));
	octvector(m_octants).swap(m_octants);
	nocts = m_octants.size();
	if(mapsize > 0){
		mapidx.resize(nocts);
	}

	// End on ghosts
	if (m_ghosts.size() && nocts > 0){
		if (m_ghosts[idx2_gh].buildFather() == m_octants[nocts-1].buildFather()){
			father = m_ghosts[idx2_gh].buildFather();
			for (uint32_t iii=0; iii<17; iii++){
				father.m_info[iii] = false;
			}
			markerfather = m_ghosts[idx2_gh].getMarker()+1;
			nbro = 0;
			idx = idx2_gh;
			marker = m_ghosts[idx].getMarker();
			while(marker < 0 && m_ghosts[idx].buildFather() == father){
				nbro++;
				if (markerfather < m_ghosts[idx].getMarker()+1){
					markerfather = m_ghosts[idx].getMarker()+1;
				}
				for (uint32_t iii=0; iii<m_global.m_nfaces; iii++){
					father.m_info[iii] = father.m_info[iii] || m_ghosts[idx].m_info[iii];
				}
				father.m_info[14] = father.m_info[14] || m_ghosts[idx].m_info[14];
				idx++;
				if(idx == m_sizeGhosts){
					break;
				}
				marker = m_ghosts[idx].getMarker();
			}
			nend = 0;
			idx = nocts-1;
			marker = m_octants[idx].getMarker();
			while(marker < 0 && m_octants[idx].buildFather() == father){
				nbro++;
				nend++;
				if (markerfather < m_octants[idx].getMarker()+1){
					markerfather = m_octants[idx].getMarker()+1;
				}
				idx--;
				if (wstop){
					break;
				}
				if (idx==0){
					wstop = true;
				}
				marker = m_octants[idx].getMarker();
			}
			if (nbro == m_global.m_nchildren){
				offset = nend;
			}
			else{
				nend = 0;
			}
		}
		if (nend != 0){
			for (idx=0; idx < nend; idx++){
				for (uint32_t iii=0; iii<16; iii++){
					father.m_info[iii] = father.m_info[iii] || m_octants[nocts-idx-1].m_info[iii];
				}
			}
			father.m_info[13] = true;
			father.m_info[15] = true;
			if (markerfather < 0 && mapsize == 0){
				docoarse = true;
			}
			father.setMarker(markerfather);
			m_octants.resize(nocts-offset, Octant(m_dim, m_global.m_maxLevel));
			m_octants.push_back(father);
			octvector(m_octants).swap(m_octants);
			nocts = m_octants.size();
			if(mapsize > 0){
				mapidx.resize(nocts);
			}
		}

	}

	// Set final first and last desc
	if(nocts>0){
		setFirstDesc();
		setLastDesc();
	}
	return docoarse;

};

// =================================================================================== //

/*! Refine local tree: refine one time all the octants
 * \param[out] mapidx mpaidx[i] = index in old octants vector of the new i-th octant (index of father if octant is new after refinement)
 * \return	true if refinement done
 */
bool
LocalTree::globalRefine(u32vector & mapidx){

	uint32_t 	idx, nocts;
	bool 		dorefine = false;

	nocts = m_octants.size();
	for (idx=0; idx<nocts; idx++){
		m_octants[idx].setMarker(1);
	}

	dorefine = refine(mapidx);

	return dorefine;

};

// =================================================================================== //

/*! Refine local tree: corse one time all the octants
 * \param[out] mapidx mpaidx[i] = index in old octants vector of the new i-th octant (index of first child if octant is new after coarsening)
 * \return	true if coarsening done
 */
bool
LocalTree::globalCoarse(u32vector & mapidx){

	uint32_t 	idx, nocts;
	bool 		dorefine = false;

	nocts = m_octants.size();
	for (idx=0; idx<nocts; idx++){
		m_octants[idx].setMarker(-1);
	}
	for (idx=0; idx<m_sizeGhosts; idx++){
		m_ghosts[idx].setMarker(-1);
	}

	dorefine = coarse(mapidx);

	return dorefine;

};

// =================================================================================== //
/*! Delete overlapping octants after coarse local tree. Check first and last descendants
 * of process before and after the local process
 * \param[in] lastDescPre Morton of last descendant of the previous process (rank-1).
 * \param[in] firstDescPost Morton of first descendant of the next process (rank+1).
 * \param[out] mapidx mpaidx[i] = index in old octants vector of the new i-th octant (index of first child if octant is new after coarsening)
 */
void
LocalTree::checkCoarse(uint64_t lastDescPre,
		uint64_t firstDescPost,
		u32vector & mapidx){

	BITPIT_UNUSED(firstDescPost);

	uint32_t		idx;
	uint32_t 		nocts;
	uint32_t 		mapsize = mapidx.size();
	uint64_t 		Morton;
	uint8_t 		toDelete = 0;

	nocts = getNumOctants();
	idx = 0;
	Morton = m_octants[idx].computeMorton();

	while(Morton <= lastDescPre && idx < nocts-1 && Morton != 0){
		// To delete, the father is in proc before me
		toDelete++;
		idx++;
		Morton = m_octants[idx].computeMorton();
	}
	if (nocts>toDelete){
		for(idx=0; idx<nocts-toDelete; idx++){
			m_octants[idx] = m_octants[idx+toDelete];
			if (mapsize>0) mapidx[idx] = mapidx[idx+toDelete];
		}
		m_octants.resize(nocts-toDelete, Octant(m_dim, m_global.m_maxLevel));
		if (mapsize>0){
			mapidx.resize(nocts-toDelete);
		}
	}
	else{
		m_octants.clear();
		mapidx.clear();
	}
	nocts = getNumOctants();

	setFirstDesc();
	setLastDesc();

};

// =================================================================================== //
/*! Update max depth reached in local tree
 */
void
LocalTree::updateLocalMaxDepth(){

	uint32_t noctants = getNumOctants();
	uint32_t i;

	m_localMaxDepth = 0;
	for(i = 0; i < noctants; i++){
		if(m_octants[i].getLevel() > m_localMaxDepth){
			m_localMaxDepth = m_octants[i].getLevel();
		}
	}
};

// =================================================================================== //

/*! Finds neighbours of idx-th octant through iface in vector m_octants or m_ghosts.
 * Returns a vector (empty if iface is a bound face) with the index of neighbours
 * in their structure (octants or ghosts) and sets isghost[i] = true if the
 * i-th neighbour is ghost in the local tree.
 * \param[in] idx Local index of the target local octant.
 * \param[in] iface local index of the face.
 * \param[out] neighbours Vector with the local indices of the neighbours (size = 0 if boundary face).
 * \param[out] isghost Vector with the information about the identity of each neighbour (is the i-th neighours a ghost octant?).
 */
void
LocalTree::findNeighbours(uint32_t idx, uint8_t iface, u32vector & neighbours, vector<bool> & isghost){

	uint64_t  		Morton, Mortontry;
	uint32_t 		noctants = getNumOctants();
	uint32_t 		idxtry;
	Octant* 		oct = &m_octants[idx];
	uint32_t 		size = oct->getSize();

	int8_t 			cxyz[3] = {0,0,0};
	for (int idim=0; idim<m_dim; idim++){
		cxyz[idim] = m_global.m_normals[iface][idim];
	}

	isghost.clear();
	neighbours.clear();

	// Default if iface is nface<iface<0
	if (iface > m_global.m_nfaces){
		return;
	}

	// Check if octants face is a process boundary
	if (oct->m_info[6+iface] == false){

		// Check if octants face is a boundary
		if (oct->m_info[iface] == false){

			//Build Morton number of virtual neigh of same size
			Octant samesizeoct(m_dim, oct->m_level, int32_t(oct->m_x)+int32_t(cxyz[0]*size), int32_t(oct->m_y)+int32_t(cxyz[1]*size), int32_t(oct->m_z)+int32_t(cxyz[2]*size), m_global.m_maxLevel);
			Morton = samesizeoct.computeMorton();
			// Search morton in octants
			// If a even face morton is lower than morton of oct, if odd higher
			// ---> can i search only before or after idx in octants
			Mortontry = oct->computeMorton();
			int32_t jump = (Mortontry > Morton) ? int32_t(idx/2+1) : int32_t((noctants -idx)/2+1);
			idxtry = uint32_t(idx + ((Mortontry < Morton) - (Mortontry > Morton))*jump);
			if (idxtry > noctants-1) idxtry = noctants-1;
			while(abs(jump) > 0){
				Mortontry = m_octants[idxtry].computeMorton();
				jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
				idxtry += jump;
				if (idxtry > noctants-1){
					if (jump > 0){
						idxtry = noctants - 1;
						jump = 0;
					}
					else if (jump < 0){
						idxtry = 0;
						jump = 0;
					}
				}
			}
			Mortontry = m_octants[idxtry].computeMorton();
			if(Mortontry == Morton && m_octants[idxtry].m_level == oct->m_level){
				//Found neighbour of same size
				isghost.push_back(false);
				neighbours.push_back(idxtry);
				return;
			}
			else{
				// Step until the mortontry lower than morton (one idx of distance)
				{
					while(Mortontry < Morton){
						idxtry++;
						if(idxtry > noctants-1){
							idxtry = noctants-1;
							Mortontry = m_octants[idxtry].computeMorton();
							break;
						}
						Mortontry = m_octants[idxtry].computeMorton();
					}
					while(Mortontry > Morton){
						idxtry--;
						if(idxtry > noctants-1){
							idxtry = 0;
							Mortontry = m_octants[idxtry].computeMorton();
							break;
						}
						Mortontry = m_octants[idxtry].computeMorton();
					}
				}
				if(Mortontry == Morton && m_octants[idxtry].m_level == oct->m_level){
					//Found neighbour of same size
					isghost.push_back(false);
					neighbours.push_back(idxtry);
					return;
				}
				// Compute Last discendent of virtual octant of same size
				Octant last_desc = samesizeoct.buildLastDesc();
				uint64_t Mortonlast = last_desc.computeMorton();
				Mortontry = m_octants[idxtry].computeMorton();
				int32_t Dx[3] = {0,0,0};
				int32_t Dxstar[3] = {0,0,0};
				u32array3 coord = oct->getCoord();
				u32array3 coordtry = m_octants[idxtry].getCoord();
				u32array3 coord1 = { {1,1,1} };
				u32array3 coordtry1 = { {1,1,1} };
				uint8_t level = oct->m_level;
				uint8_t leveltry = m_octants[idxtry].getLevel();
				while(Mortontry <= Mortonlast && idxtry < noctants){
					for (int idim=0; idim<m_dim; idim++){
						Dx[idim] 		= int32_t(int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]));
						Dxstar[idim]	= int32_t((cxyz[idim]-1)/2)*(m_octants[idxtry].getSize()) + int32_t((cxyz[idim]+1)/2)*size;
						coord1[idim] 	= coord[idim] + size;
						coordtry1[idim] = coordtry[idim] + m_octants[idxtry].getSize();
					}
					leveltry = m_octants[idxtry].getLevel();

					if (Dx[0] == Dxstar[0] && Dx[1] == Dxstar[1] && Dx[m_dim-1] == Dxstar[m_dim-1]){
						if (leveltry > level){
							if((abs(cxyz[0])*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[1])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[2])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1])))){
								neighbours.push_back(idxtry);
								isghost.push_back(false);
							}
						}
						if (leveltry < level){
							if((abs(cxyz[0])*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[1])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[2])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1])))){
								neighbours.push_back(idxtry);
								isghost.push_back(false);
							}
						}
					}

					idxtry++;
					if(idxtry>noctants-1){
						break;
					}
					Mortontry = m_octants[idxtry].computeMorton();
					coordtry = m_octants[idxtry].getCoord();
				}
				return;
			}
		}
		else if(m_periodic[iface]){
			// Check if octants face is a periodic boundary
			findPeriodicNeighbours(oct, iface, neighbours, isghost);
			return;
		}
		else{
		// Boundary Face
		return;
		}


	}
	else{
		// Check if octants face is a boundary
		if (oct->m_info[iface] == false){

			// IF OCTANT FACE IS A PROCESS BOUNDARY SEARCH ALSO IN GHOSTS
			if (m_ghosts.size()>0){

				// Search in ghosts
				uint32_t idxghost = uint32_t(m_sizeGhosts/2);
				Octant* octghost = &m_ghosts[idxghost];

				//Build Morton number of virtual neigh of same size
				Octant samesizeoct(m_dim, oct->m_level, int32_t(oct->m_x)+int32_t(cxyz[0]*size), int32_t(oct->m_y)+int32_t(cxyz[1]*size), int32_t(oct->m_z)+int32_t(cxyz[2]*size), m_global.m_maxLevel);
				Morton = samesizeoct.computeMorton(); //mortonEncode_magicbits(oct->m_x-size,oct->m_y,oct->m_z);
				// Search morton in octants
				// If a even face morton is lower than morton of oct, if odd higher
				// ---> can i search only before or after idx in octants
				Mortontry = octghost->computeMorton();
				int32_t jump = (Mortontry > Morton) ? int32_t(idxghost/2+1) : int32_t((m_sizeGhosts -idxghost)/2+1);
				idxtry = uint32_t(idxghost +((Mortontry<Morton)-(Mortontry>Morton))*jump);
				if (idxtry > m_ghosts.size()-1) idxtry = m_ghosts.size()-1;
				while(abs(jump) > 0){
					Mortontry = m_ghosts[idxtry].computeMorton();
					jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
					idxtry += jump;
					if (idxtry > m_ghosts.size()-1){
						if (jump > 0){
							idxtry = m_ghosts.size() - 1;
							Mortontry = m_ghosts[idxtry].computeMorton();
							jump = 0;
						}
						else if (jump < 0){
							idxtry = 0;
							Mortontry = m_ghosts[idxtry].computeMorton();
							jump = 0;
						}
					}
				}
				Mortontry = m_ghosts[idxtry].computeMorton();
				if(Mortontry == Morton && m_ghosts[idxtry].m_level == oct->m_level){
					//Found neighbour of same size
					isghost.push_back(true);
					neighbours.push_back(idxtry);
					return;
				}
				else{
					// Step until the mortontry lower than morton (one idx of distance)
					{
						while(Mortontry < Morton){
							idxtry++;
							if(idxtry > m_ghosts.size()-1){
								idxtry = m_ghosts.size()-1;
								Mortontry = m_ghosts[idxtry].computeMorton();
								break;
							}
							Mortontry = m_ghosts[idxtry].computeMorton();
						}
						while(m_ghosts[idxtry].computeMorton() > Morton){
							idxtry--;
							if(idxtry > m_ghosts.size()-1){
								idxtry = 0;
								Mortontry = m_ghosts[idxtry].computeMorton();
								break;
							}
							Mortontry = m_ghosts[idxtry].computeMorton();
						}
					}
					if(idxtry < m_sizeGhosts){
						if(m_ghosts[idxtry].computeMorton() == Morton && m_ghosts[idxtry].m_level == oct->m_level){
							//Found neighbour of same size
							isghost.push_back(true);
							neighbours.push_back(idxtry);
							return;
						}
						// Compute Last discendent of virtual octant of same size
						Octant last_desc = samesizeoct.buildLastDesc();
						uint64_t Mortonlast = last_desc.computeMorton();
						Mortontry = m_ghosts[idxtry].computeMorton();
						int32_t Dx[3] = {0,0,0};
						int32_t Dxstar[3] = {0,0,0};
						u32array3 coord = oct->getCoord();
						u32array3 coordtry = m_ghosts[idxtry].getCoord();
						u32array3 coord1 = { {1,1,1} };
						u32array3 coordtry1 = { {1,1,1} };
						uint8_t level = oct->m_level;
						uint8_t leveltry = m_octants[idxtry].getLevel();
						while(Mortontry <= Mortonlast && idxtry < m_sizeGhosts){
							for (int idim=0; idim<m_dim; idim++){
								Dx[idim] 		= int32_t(int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]));
								Dxstar[idim]	= int32_t((cxyz[idim]-1)/2)*(m_ghosts[idxtry].getSize()) + int32_t((cxyz[idim]+1)/2)*size;
								coord1[idim] 	= coord[idim] + size;
								coordtry1[idim] = coordtry[idim] + m_ghosts[idxtry].getSize();
							}
							leveltry = m_ghosts[idxtry].getLevel();

							if (Dx[0] == Dxstar[0] && Dx[1] == Dxstar[1] && Dx[m_dim-1] == Dxstar[m_dim-1]){
								if (leveltry > level){
									if((abs(cxyz[0])*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[1])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[2])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1])))){
										neighbours.push_back(idxtry);
										isghost.push_back(true);
									}
								}
								if (leveltry < level){
									if((abs(cxyz[0])*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[1])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[2])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1])))){
										neighbours.push_back(idxtry);
										isghost.push_back(true);
									}
								}
							}

							idxtry++;
							if(idxtry>m_sizeGhosts-1){
								break;
							}
							Mortontry = m_ghosts[idxtry].computeMorton();
							coordtry = m_ghosts[idxtry].getCoord();
						}
					}
				}

				uint32_t lengthneigh = 0;
				uint32_t sizeneigh = neighbours.size();
				for (idxtry=0; idxtry<sizeneigh; idxtry++){
					lengthneigh += m_ghosts[neighbours[idxtry]].getArea();
				}
				if (lengthneigh < oct->getArea()){
					// Search in octants

					// Check if octants face is a boundary
					if (oct->m_info[iface] == false){

						//Build Morton number of virtual neigh of same size
						Octant samesizeoct(m_dim, oct->m_level, int32_t(oct->m_x)+int32_t(cxyz[0]*size), int32_t(oct->m_y)+int32_t(cxyz[1]*size), int32_t(oct->m_z)+int32_t(cxyz[2]*size), m_global.m_maxLevel);
						Morton = samesizeoct.computeMorton();
						// Search morton in octants
						// If a even face morton is lower than morton of oct, if odd higher
						// ---> can i search only before or after idx in octants
						Mortontry = oct->computeMorton();
						int32_t jump = (Mortontry > Morton) ? int32_t(idx/2+1) : int32_t((noctants -idx)/2+1);
						idxtry = uint32_t(idx + ((Mortontry < Morton) - (Mortontry > Morton))*jump);
						if (idxtry > noctants-1) idxtry = noctants-1;
						while(abs(jump) > 0){
							Mortontry = m_octants[idxtry].computeMorton();
							jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
							idxtry += jump;
							if (idxtry > noctants-1){
								if (jump > 0){
									idxtry = noctants - 1;
									jump = 0;
								}
								else if (jump < 0){
									idxtry = 0;
									jump = 0;
								}
							}
						}
						Mortontry = m_octants[idxtry].computeMorton();
						if(Mortontry == Morton && m_octants[idxtry].m_level == oct->m_level){
							//Found neighbour of same size
							isghost.push_back(false);
							neighbours.push_back(idxtry);
							return;
						}
						else{
							// Step until the mortontry lower than morton (one idx of distance)
							{
								while(Mortontry < Morton){
									idxtry++;
									if(idxtry > noctants-1){
										idxtry = noctants-1;
										Mortontry = m_octants[idxtry].computeMorton();
										break;
									}
									Mortontry = m_octants[idxtry].computeMorton();
								}
								while(Mortontry > Morton){
									idxtry--;
									if(idxtry > noctants-1){
										idxtry = 0;
										Mortontry = m_octants[idxtry].computeMorton();
										break;
									}
									Mortontry = m_octants[idxtry].computeMorton();
								}
							}
							if(Mortontry == Morton && m_octants[idxtry].m_level == oct->m_level){
								//Found neighbour of same size
								isghost.push_back(false);
								neighbours.push_back(idxtry);
								return;
							}
							// Compute Last discendent of virtual octant of same size
							Octant last_desc = samesizeoct.buildLastDesc();
							uint64_t Mortonlast = last_desc.computeMorton();
							Mortontry = m_octants[idxtry].computeMorton();
							int32_t Dx[3] = {0,0,0};
							int32_t Dxstar[3] = {0,0,0};
							u32array3 coord = oct->getCoord();
							u32array3 coordtry = m_octants[idxtry].getCoord();
							u32array3 coord1 = { {1,1,1} };
							u32array3 coordtry1 = { {1,1,1} };
							uint8_t level = oct->m_level;
							uint8_t leveltry = m_octants[idxtry].getLevel();
							while(Mortontry <= Mortonlast && idxtry < noctants){
								for (int idim=0; idim<m_dim; idim++){
									Dx[idim] 		= int32_t(int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]));
									Dxstar[idim]	= int32_t((cxyz[idim]-1)/2)*(m_octants[idxtry].getSize()) + int32_t((cxyz[idim]+1)/2)*size;
									coord1[idim] 	= coord[idim] + size;
									coordtry1[idim] = coordtry[idim] + m_octants[idxtry].getSize();
								}
								leveltry = m_octants[idxtry].getLevel();

								if (Dx[0] == Dxstar[0] && Dx[1] == Dxstar[1] && Dx[m_dim-1] == Dxstar[m_dim-1]){
									if (leveltry > level){
										if((abs(cxyz[0])*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[1])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[2])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1])))){
											neighbours.push_back(idxtry);
											isghost.push_back(false);
										}
									}
									if (leveltry < level){
										if((abs(cxyz[0])*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[1])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[2])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1])))){
											neighbours.push_back(idxtry);
											isghost.push_back(false);
										}
									}
								}

								idxtry++;
								if(idxtry>noctants-1){
									break;
								}
								Mortontry = m_octants[idxtry].computeMorton();
								coordtry = m_octants[idxtry].getCoord();
							}
							return;
						}
					}
				}
				return;
			}
		}
		else if(m_periodic[iface]){
			// Check if octants face is a periodic boundary
			findPeriodicNeighbours(oct, iface, neighbours, isghost);
			return;
		}
		else{
		// Boundary Face
		return;
		}

	}
};

// =================================================================================== //

/*! Finds neighbours of octant through iface in vector m_octants or m_ghosts.
 * Returns a vector (empty if iface is a bound face) with the index of neighbours
 * in their structure (octants or ghosts) and sets isghost[i] = true if the
 * i-th neighbour is ghost in the local tree.
 * \param[in] oct Pointer to the target local octant.
 * \param[in] iface local index of the face.
 * \param[out] neighbours Vector with the local indices of the neighbours (size = 0 if boundary face).
 * \param[out] isghost Vector with the information about the identity of each neighbour (is the i-th neighours a ghost octant?).
 */
void
LocalTree::findNeighbours(Octant* oct, uint8_t iface, u32vector & neighbours, vector<bool> & isghost){

	uint64_t  Morton, Mortontry;
	uint32_t  noctants = getNumOctants();
	uint32_t idxtry;
	uint32_t size = oct->getSize();

	int8_t 			cxyz[3] = {0,0,0};
	for (int idim=0; idim<m_dim; idim++){
		cxyz[idim] = m_global.m_normals[iface][idim];
	}

	isghost.clear();
	neighbours.clear();

	// Default if iface is nface<iface<0
	if (iface > m_global.m_nfaces){
		return;
	}

	// Check if octants face is a process boundary
	if (oct->m_info[6+iface] == false){

		// Check if octants face is a boundary
		if (oct->m_info[iface] == false){

			//Build Morton number of virtual neigh of same size
			Octant samesizeoct(m_dim, oct->m_level, int32_t(oct->m_x)+int32_t(cxyz[0]*size), int32_t(oct->m_y)+int32_t(cxyz[1]*size), int32_t(oct->m_z)+int32_t(cxyz[2]*size), m_global.m_maxLevel);
			Morton = samesizeoct.computeMorton();
			// Search morton in octants
			// If a even face morton is lower than morton of oct, if odd higher
			// ---> can i search only before or after idx in octants
			int32_t jump = int32_t((noctants)/2+1);
			idxtry = uint32_t(jump);
			Mortontry = oct->computeMorton();
			while(abs(jump) > 0){
				Mortontry = m_octants[idxtry].computeMorton();
				jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
				idxtry += jump;
				if (idxtry > noctants-1){
					if (jump > 0){
						idxtry = noctants - 1;
						jump = 0;
					}
					else if (jump < 0){
						idxtry = 0;
						jump = 0;
					}
				}
			}
			Mortontry = m_octants[idxtry].computeMorton();
			if(Mortontry == Morton && m_octants[idxtry].m_level == oct->m_level){
				//Found neighbour of same size
				isghost.push_back(false);
				neighbours.push_back(idxtry);
				return;
			}
			else{
				// Step until the mortontry lower than morton (one idx of distance)
				{
					while(Mortontry < Morton){
						idxtry++;
						if(idxtry > noctants-1){
							idxtry = noctants-1;
							Mortontry = m_octants[idxtry].computeMorton();
							break;
						}
						Mortontry = m_octants[idxtry].computeMorton();
					}
					while(Mortontry > Morton){
						idxtry--;
						if(idxtry > noctants-1){
							idxtry = 0;
							Mortontry = m_octants[idxtry].computeMorton();
							break;
						}
						Mortontry = m_octants[idxtry].computeMorton();
					}
				}
				if(Mortontry == Morton && m_octants[idxtry].m_level == oct->m_level){
					//Found neighbour of same size
					isghost.push_back(false);
					neighbours.push_back(idxtry);
					return;
				}
				// Compute Last discendent of virtual octant of same size
				Octant last_desc = samesizeoct.buildLastDesc();
				uint64_t Mortonlast = last_desc.computeMorton();
				Mortontry = m_octants[idxtry].computeMorton();
				int32_t Dx[3] = {0,0,0};
				int32_t Dxstar[3] = {0,0,0};
				u32array3 coord = oct->getCoord();
				u32array3 coordtry = m_octants[idxtry].getCoord();
				u32array3 coord1 = { {1,1,1} };
				u32array3 coordtry1 = { {1,1,1} };
				uint8_t level = oct->m_level;
				uint8_t leveltry = m_octants[idxtry].getLevel();
				while(Mortontry <= Mortonlast && idxtry < noctants){
					for (int idim=0; idim<m_dim; idim++){
						Dx[idim] 		= int32_t(int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]));
						Dxstar[idim]	= int32_t((cxyz[idim]-1)/2)*(m_octants[idxtry].getSize()) + int32_t((cxyz[idim]+1)/2)*size;
						coord1[idim] 	= coord[idim] + size;
						coordtry1[idim] = coordtry[idim] + m_octants[idxtry].getSize();
					}

					leveltry = m_octants[idxtry].getLevel();

					if (Dx[0] == Dxstar[0] && Dx[1] == Dxstar[1] && Dx[m_dim-1] == Dxstar[m_dim-1]){
						if (leveltry > level){
							if((abs(cxyz[0])*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[1])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[2])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1])))){
								neighbours.push_back(idxtry);
								isghost.push_back(false);
							}
						}
						if (leveltry < level){
							if((abs(cxyz[0])*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[1])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[2])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1])))){
								neighbours.push_back(idxtry);
								isghost.push_back(false);
							}
						}
					}

					idxtry++;
					if(idxtry>noctants-1){
						break;
					}
					Mortontry = m_octants[idxtry].computeMorton();
					coordtry = m_octants[idxtry].getCoord();
				}
				return;
			}
		}
		else if(m_periodic[iface]){
			// Check if octants face is a periodic boundary
			findPeriodicNeighbours(oct, iface, neighbours, isghost);
			return;
		}
		else{
		// Boundary Face
		return;
		}

	}
	else{
		// Check if octants face is a boundary
		if (oct->m_info[iface] == false){

			// IF OCTANT FACE IS A PROCESS BOUNDARY SEARCH ALSO IN GHOSTS
			if (m_ghosts.size()>0){

				// Search in ghosts
				uint32_t idxghost = uint32_t(m_sizeGhosts/2);
				Octant* octghost = &m_ghosts[idxghost];

				//Build Morton number of virtual neigh of same size
				Octant samesizeoct(m_dim, oct->m_level, int32_t(oct->m_x)+int32_t(cxyz[0]*size), int32_t(oct->m_y)+int32_t(cxyz[1]*size), int32_t(oct->m_z)+int32_t(cxyz[2]*size), m_global.m_maxLevel);
				Morton = samesizeoct.computeMorton(); //mortonEncode_magicbits(oct->m_x-size,oct->m_y,oct->m_z);
				// Search morton in octants
				// If a even face morton is lower than morton of oct, if odd higher
				// ---> can i search only before or after idx in octants
				Mortontry = octghost->computeMorton();
				int32_t jump = (Mortontry > Morton) ? int32_t(idxghost/2+1) : int32_t((m_sizeGhosts -idxghost)/2+1);
				idxtry = uint32_t(idxghost +((Mortontry<Morton)-(Mortontry>Morton))*jump);
				if (idxtry > m_ghosts.size()-1) idxtry = m_ghosts.size()-1;
				while(abs(jump) > 0){
					Mortontry = m_ghosts[idxtry].computeMorton();
					jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
					idxtry += jump;
					if (idxtry > m_ghosts.size()-1){
						if (jump > 0){
							idxtry = m_ghosts.size() - 1;
							Mortontry = m_ghosts[idxtry].computeMorton();
							jump = 0;
						}
						else if (jump < 0){
							idxtry = 0;
							Mortontry = m_ghosts[idxtry].computeMorton();
							jump = 0;
						}
					}
				}
				Mortontry = m_ghosts[idxtry].computeMorton();
				if(Mortontry == Morton && m_ghosts[idxtry].m_level == oct->m_level){
					//Found neighbour of same size
					isghost.push_back(true);
					neighbours.push_back(idxtry);
					return;
				}
				else{
					// Step until the mortontry lower than morton (one idx of distance)
					{
						while(Mortontry < Morton){
							idxtry++;
							if(idxtry > m_ghosts.size()-1){
								idxtry = m_ghosts.size()-1;
								Mortontry = m_ghosts[idxtry].computeMorton();
								break;
							}
							Mortontry = m_ghosts[idxtry].computeMorton();
						}
						while(m_ghosts[idxtry].computeMorton() > Morton){
							idxtry--;
							if(idxtry > m_ghosts.size()-1){
								idxtry = 0;
								Mortontry = m_ghosts[idxtry].computeMorton();
								break;
							}
							Mortontry = m_ghosts[idxtry].computeMorton();
						}
					}
					if(idxtry < m_sizeGhosts){
						if(m_ghosts[idxtry].computeMorton() == Morton && m_ghosts[idxtry].m_level == oct->m_level){
							//Found neighbour of same size
							isghost.push_back(true);
							neighbours.push_back(idxtry);
							return;
						}
						// Compute Last discendent of virtual octant of same size
						Octant last_desc = samesizeoct.buildLastDesc();
						uint64_t Mortonlast = last_desc.computeMorton();
						Mortontry = m_ghosts[idxtry].computeMorton();
						int32_t Dx[3] = {0,0,0};
						int32_t Dxstar[3] = {0,0,0};
						u32array3 coord = oct->getCoord();
						u32array3 coordtry = m_ghosts[idxtry].getCoord();
						u32array3 coord1 = { {1,1,1} };
						u32array3 coordtry1 = { {1,1,1} };
						uint8_t level = oct->m_level;
						uint8_t leveltry = m_octants[idxtry].getLevel();
						while(Mortontry <= Mortonlast && idxtry < m_sizeGhosts){
							for (int idim=0; idim<m_dim; idim++){
								Dx[idim] 		= int32_t(int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]));
								Dxstar[idim]	= int32_t((cxyz[idim]-1)/2)*(m_ghosts[idxtry].getSize()) + int32_t((cxyz[idim]+1)/2)*size;
								coord1[idim] 	= coord[idim] + size;
								coordtry1[idim] = coordtry[idim] + m_ghosts[idxtry].getSize();
							}
							leveltry = m_ghosts[idxtry].getLevel();

							if (Dx[0] == Dxstar[0] && Dx[1] == Dxstar[1] && Dx[m_dim-1] == Dxstar[m_dim-1]){
								if (leveltry > level){
									if((abs(cxyz[0])*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[1])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[2])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1])))){
										neighbours.push_back(idxtry);
										isghost.push_back(true);
									}
								}
								if (leveltry < level){
									if((abs(cxyz[0])*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[1])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[2])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1])))){
										neighbours.push_back(idxtry);
										isghost.push_back(true);
									}
								}
							}

							idxtry++;
							if(idxtry>m_sizeGhosts-1){
								break;
							}
							Mortontry = m_ghosts[idxtry].computeMorton();
							coordtry = m_ghosts[idxtry].getCoord();
						}
					}
				}

				uint32_t lengthneigh = 0;
				uint32_t sizeneigh = neighbours.size();
				for (idxtry=0; idxtry<sizeneigh; idxtry++){
					lengthneigh += m_ghosts[neighbours[idxtry]].getArea();
				}
				if (lengthneigh < oct->getArea()){
					// Search in octants

					// Check if octants face is a boundary
					if (oct->m_info[iface] == false){

						//Build Morton number of virtual neigh of same size
						Octant samesizeoct(m_dim, oct->m_level, int32_t(oct->m_x)+int32_t(cxyz[0]*size), int32_t(oct->m_y)+int32_t(cxyz[1]*size), int32_t(oct->m_z)+int32_t(cxyz[2]*size), m_global.m_maxLevel);
						Morton = samesizeoct.computeMorton();
						// Search morton in octants
						// If a even face morton is lower than morton of oct, if odd higher
						// ---> can i search only before or after idx in octants
						int32_t jump = (int32_t((noctants)/2+1));
						idxtry = uint32_t(jump);
						while(abs(jump) > 0){
							Mortontry = m_octants[idxtry].computeMorton();
							jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
							idxtry += jump;
							if (idxtry > noctants-1){
								if (jump > 0){
									idxtry = noctants - 1;
									jump = 0;
								}
								else if (jump < 0){
									idxtry = 0;
									jump = 0;
								}
							}
						}
						Mortontry = m_octants[idxtry].computeMorton();
						if(Mortontry == Morton && m_octants[idxtry].m_level == oct->m_level){
							//Found neighbour of same size
							isghost.push_back(false);
							neighbours.push_back(idxtry);
							return;
						}
						else{
							// Step until the mortontry lower than morton (one idx of distance)
							{
								while(Mortontry < Morton){
									idxtry++;
									if(idxtry > noctants-1){
										idxtry = noctants-1;
										Mortontry = m_octants[idxtry].computeMorton();
										break;
									}
									Mortontry = m_octants[idxtry].computeMorton();
								}
								while(Mortontry > Morton){
									idxtry--;
									if(idxtry > noctants-1){
										idxtry = 0;
										Mortontry = m_octants[idxtry].computeMorton();
										break;
									}
									Mortontry = m_octants[idxtry].computeMorton();
								}
							}
							if(Mortontry == Morton && m_octants[idxtry].m_level == oct->m_level){
								//Found neighbour of same size
								isghost.push_back(false);
								neighbours.push_back(idxtry);
								return;
							}
							// Compute Last discendent of virtual octant of same size
							Octant last_desc = samesizeoct.buildLastDesc();
							uint64_t Mortonlast = last_desc.computeMorton();
							Mortontry = m_octants[idxtry].computeMorton();
							int32_t Dx[3] = {0,0,0};
							int32_t Dxstar[3] = {0,0,0};
							u32array3 coord = oct->getCoord();
							u32array3 coordtry = m_octants[idxtry].getCoord();
							u32array3 coord1 = { {1,1,1} };
							u32array3 coordtry1 = { {1,1,1} };
							uint8_t level = oct->m_level;
							uint8_t leveltry = m_octants[idxtry].getLevel();
							while(Mortontry <= Mortonlast && idxtry < noctants){
								for (int idim=0; idim<m_dim; idim++){
									Dx[idim] 		= int32_t(int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]));
									Dxstar[idim]	= int32_t((cxyz[idim]-1)/2)*(m_octants[idxtry].getSize()) + int32_t((cxyz[idim]+1)/2)*size;
									coord1[idim] 	= coord[idim] + size;
									coordtry1[idim] = coordtry[idim] + m_octants[idxtry].getSize();
								}
								leveltry = m_octants[idxtry].getLevel();

								if (Dx[0] == Dxstar[0] && Dx[1] == Dxstar[1] && Dx[m_dim-1] == Dxstar[m_dim-1]){
									if (leveltry > level){
										if((abs(cxyz[0])*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[1])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[2])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1])))){
											neighbours.push_back(idxtry);
											isghost.push_back(false);
										}
									}
									if (leveltry < level){
										if((abs(cxyz[0])*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[1])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[2])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1])))){
											neighbours.push_back(idxtry);
											isghost.push_back(false);
										}
									}
								}

								idxtry++;
								if(idxtry>noctants-1){
									break;
								}
								Mortontry = m_octants[idxtry].computeMorton();
								coordtry = m_octants[idxtry].getCoord();
							}
							return;
						}
					}
				}
				return;
			}
		}
		else if(m_periodic[iface]){
			// Check if octants face is a periodic boundary
			findPeriodicNeighbours(oct, iface, neighbours, isghost);
			return;
		}
		else{
		// Boundary Face
		return;
		}
	}
};

// =================================================================================== //

/*! Finds neighbours of idx-th ghost octant through iface in vector m_octants.
 * Returns a vector (empty if iface is not the pbound face for ghost)
 * with the index of neighbours in the structure octants.
 * \param[in] oct Local index of the target ghost octant.
 * \param[in] iface local index of the face.
 * \param[out] neighbours Vector with the local indices of the local octant neighbours (size = 0 if boundary face).
 */
void
LocalTree::findGhostNeighbours(uint32_t const idx, uint8_t iface, u32vector & neighbours){

	uint64_t  Morton, Mortontry;
	uint32_t  noctants = getNumOctants();
	uint32_t idxtry;
	Octant* oct = &m_ghosts[idx];
	uint32_t size = oct->getSize();

	int8_t 			cxyz[3] = {0,0,0};
	for (int idim=0; idim<m_dim; idim++){
		cxyz[idim] = m_global.m_normals[iface][idim];
	}

	neighbours.clear();

	// Default if iface is nface<iface<0
	if (iface > m_global.m_nfaces){
		return;
	}

	// Check if octants face is a process boundary
	if (oct->m_info[6+iface] == true){

		// Check if octants face is a boundary
		if (oct->m_info[iface] == false){

			//Build Morton number of virtual neigh of same size
			Octant samesizeoct(m_dim, oct->m_level, int32_t(oct->m_x)+int32_t(cxyz[0]*size), int32_t(oct->m_y)+int32_t(cxyz[1]*size), int32_t(oct->m_z)+int32_t(cxyz[2]*size), m_global.m_maxLevel);
			Morton = samesizeoct.computeMorton();
			// Search morton in octants
			// If a even face morton is lower than morton of oct, if odd higher
			// ---> can i search only before or after idx in octants
			int32_t jump = getNumOctants()/2;
			idxtry = uint32_t(getNumOctants()/2);
			Mortontry = m_octants[idxtry].computeMorton();
//			jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
			while(abs(jump) > 0){
				Mortontry = m_octants[idxtry].computeMorton();
				jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
				idxtry += jump;
				if (idxtry > noctants-1){
					if (jump > 0){
						idxtry = noctants - 1;
						Mortontry = m_octants[idxtry].computeMorton();
						jump = 0;
					}
					else if (jump < 0){
						idxtry = 0;
						Mortontry = m_octants[idxtry].computeMorton();
						jump = 0;
					}
				}
			}
			Mortontry = m_octants[idxtry].computeMorton();
			if(Mortontry == Morton && m_octants[idxtry].m_level == oct->m_level){
				//Found neighbour of same size
				neighbours.push_back(idxtry);
				return;
			}
			else{
				// Step until the mortontry lower than morton (one idx of distance)
				{
					while(Mortontry < Morton){
						idxtry++;
						if(idxtry > noctants-1){
							idxtry = noctants-1;
							Mortontry = m_octants[idxtry].computeMorton();
							break;
						}
						Mortontry = m_octants[idxtry].computeMorton();
					}
					while(Mortontry > Morton){
						idxtry--;
						if(idxtry > noctants-1){
							idxtry = 0;
							Mortontry = m_octants[idxtry].computeMorton();
							break;
						}
						Mortontry = m_octants[idxtry].computeMorton();
					}
				}
				if(Mortontry == Morton && m_octants[idxtry].m_level == oct->m_level){
					//Found neighbour of same size
					neighbours.push_back(idxtry);
					return;
				}
				// Compute Last discendent of virtual octant of same size
				Octant last_desc = samesizeoct.buildLastDesc();
				uint64_t Mortonlast = last_desc.computeMorton();
				Mortontry = m_octants[idxtry].computeMorton();
				int32_t Dx[3] = {0,0,0};
				int32_t Dxstar[3] = {0,0,0};
				u32array3 coord = oct->getCoord();
				u32array3 coordtry = m_octants[idxtry].getCoord();
				u32array3 coord1 = { {1,1,1} };
				u32array3 coordtry1 = { {1,1,1} };
				uint8_t level = oct->m_level;
				while(Mortontry <= Mortonlast && idxtry < noctants){
					for (int idim=0; idim<m_dim; idim++){
						Dx[idim] 		= int32_t(int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]));
						Dxstar[idim]	= int32_t((cxyz[idim]-1)/2)*(m_octants[idxtry].getSize()) + int32_t((cxyz[idim]+1)/2)*size;
						coord1[idim] 	= coord[idim] + size;
						coordtry1[idim] = coordtry[idim] + m_octants[idxtry].getSize();
					}
					uint8_t leveltry = m_octants[idxtry].getLevel();

					if (Dx[0] == Dxstar[0] && Dx[1] == Dxstar[1] && Dx[m_dim-1] == Dxstar[m_dim-1]){
						if (leveltry > level){
							if((abs(cxyz[0])*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[1])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[2])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1])))){
								neighbours.push_back(idxtry);
							}
						}
						if (leveltry < level){
							if((abs(cxyz[0])*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[1])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[2])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1])))){
								neighbours.push_back(idxtry);
							}
						}
					}

					idxtry++;
					if(idxtry>noctants-1){
						break;
					}
					Mortontry = m_octants[idxtry].computeMorton();
					coordtry = m_octants[idxtry].getCoord();
				}
				return;
			}
		}
		else if(m_periodic[iface]){
			// Check if octants face is a periodic boundary
			findGhostPeriodicNeighbours(oct, iface, neighbours);
			return;
		}
		else{
			// Boundary Face
			return;
		}
	}

};

// =================================================================================== //

/*! Finds neighbours of octant through a periodic iface in vector m_octants or m_ghosts.
 * Returns a vector with the index of neighbours in their structure (octants or ghosts)
 * and sets isghost[i] = true if the i-th neighbour is ghost in the local tree.
 * \param[in] oct Pointer to the target local octant.
 * \param[in] iface local index of the face.
 * \param[out] neighbours Vector with the local indices of the neighbours
 * \param[out] isghost Vector with the information about the identity of each neighbour (is the i-th neighours a ghost octant?).
 */
void
LocalTree::findPeriodicNeighbours(Octant* oct, uint8_t iface, u32vector & neighbours, vector<bool> & isghost){

	uint64_t  Morton, Mortontry;
	uint32_t  noctants = getNumOctants();
	uint32_t idxtry;
	uint32_t size = oct->getSize();

	int8_t 			cxyz[3] = {0,0,0};
	for (int idim=0; idim<m_dim; idim++){
		cxyz[idim] = m_global.m_normals[iface][idim];
	}

	isghost.clear();
	neighbours.clear();

	// Default if iface is nface<iface<0
	if (iface > m_global.m_nfaces){
		return;
	}

	// Check if octants face is a process boundary
	if (oct->m_info[6+iface] == false){

		//Build Morton number of virtual neigh of same size
		Octant samesizeoct = oct->computePeriodicOctant(iface);
		Morton = samesizeoct.computeMorton();
		// Search morton in octants
		// If a even face morton is lower than morton of oct, if odd higher
		// ---> can i search only before or after idx in octants
		int32_t jump = int32_t((noctants)/2+1);
		idxtry = uint32_t(jump);
		Mortontry = m_octants[idxtry].computeMorton();
		while(abs(jump) > 0){
			Mortontry = m_octants[idxtry].computeMorton();
			jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
			idxtry += jump;
			if (idxtry > noctants-1){
				if (jump > 0){
					idxtry = noctants - 1;
					jump = 0;
				}
				else if (jump < 0){
					idxtry = 0;
					jump = 0;
				}
			}
		}
		Mortontry = m_octants[idxtry].computeMorton();
		if(Mortontry == Morton && m_octants[idxtry].m_level == oct->m_level){
			//Found neighbour of same size
			isghost.push_back(false);
			neighbours.push_back(idxtry);
			return;
		}
		else{
			// Step until the mortontry lower than morton (one idx of distance)
			{
				while(Mortontry < Morton){
					idxtry++;
					if(idxtry > noctants-1){
						idxtry = noctants-1;
						Mortontry = m_octants[idxtry].computeMorton();
						break;
					}
					Mortontry = m_octants[idxtry].computeMorton();
				}
				while(Mortontry > Morton){
					idxtry--;
					if(idxtry > noctants-1){
						idxtry = 0;
						Mortontry = m_octants[idxtry].computeMorton();
						break;
					}
					Mortontry = m_octants[idxtry].computeMorton();
				}
			}

			if(Mortontry == Morton && m_octants[idxtry].m_level == oct->m_level){
				//Found neighbour of same size
				isghost.push_back(false);
				neighbours.push_back(idxtry);
				return;
			}
			// Compute Last discendent of virtual octant of same size
			Octant last_desc = samesizeoct.buildLastDesc();
			uint64_t Mortonlast = last_desc.computeMorton();
			Mortontry = m_octants[idxtry].computeMorton();
			int64_t Dx[3] = {0,0,0};
			int64_t Dxstar[3] = {0,0,0};
			array<int64_t,3> coord = oct->getPeriodicCoord(iface);
			u32array3 coordtry = m_octants[idxtry].getCoord();
			array<int64_t,3> coord1 = { {1,1,1} };
			u32array3 coordtry1 = { {1,1,1} };
			uint8_t level = oct->m_level;
			uint8_t leveltry = m_octants[idxtry].getLevel();
			while(Mortontry <= Mortonlast && idxtry < noctants){
				for (int idim=0; idim<m_dim; idim++){
					Dx[idim] 		= int32_t(int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]));
					Dxstar[idim]	= int32_t(int32_t((cxyz[idim]-1)/2)*(m_octants[idxtry].getSize()) + int32_t((cxyz[idim]+1)/2)*size);
					coord1[idim] 	= coord[idim] + size;
					coordtry1[idim] = coordtry[idim] + m_octants[idxtry].getSize();
				}

				leveltry = m_octants[idxtry].getLevel();

				if (Dx[0] == Dxstar[0] && Dx[1] == Dxstar[1] && Dx[m_dim-1] == Dxstar[m_dim-1]){
					if (leveltry > level){
						if((abs(cxyz[0])*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[1])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[2])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1])))){
							neighbours.push_back(idxtry);
							isghost.push_back(false);
						}
					}
					if (leveltry < level){
						if((abs(cxyz[0])*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[1])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[2])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1])))){
							neighbours.push_back(idxtry);
							isghost.push_back(false);
						}
					}
				}

				idxtry++;
				if(idxtry>noctants-1){
					break;
				}
				Mortontry = m_octants[idxtry].computeMorton();
				coordtry = m_octants[idxtry].getCoord();
			}
			return;
		}
	}
	else{

			// IF OCTANT FACE IS A PROCESS BOUNDARY SEARCH ALSO IN GHOSTS
			if (m_ghosts.size()>0){

				// Search in ghosts
				uint32_t idxghost = uint32_t(m_sizeGhosts/2);
				Octant* octghost = &m_ghosts[idxghost];

				//Build Morton number of virtual neigh of same size
				Octant samesizeoct = oct->computePeriodicOctant(iface);
				Morton = samesizeoct.computeMorton();
				// Search morton in octants
				// If a even face morton is lower than morton of oct, if odd higher
				// ---> can i search only before or after idx in octants
				Mortontry = octghost->computeMorton();
				int32_t jump = (Mortontry > Morton) ? int32_t(idxghost/2+1) : int32_t((m_sizeGhosts -idxghost)/2+1);
				idxtry = uint32_t(idxghost +((Mortontry<Morton)-(Mortontry>Morton))*jump);
				if (idxtry > m_ghosts.size()-1) idxtry = m_ghosts.size()-1;
				while(abs(jump) > 0){
					Mortontry = m_ghosts[idxtry].computeMorton();
					jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
					idxtry += jump;
					if (idxtry > m_ghosts.size()-1){
						if (jump > 0){
							idxtry = m_ghosts.size() - 1;
							Mortontry = m_ghosts[idxtry].computeMorton();
							jump = 0;
						}
						else if (jump < 0){
							idxtry = 0;
							Mortontry = m_ghosts[idxtry].computeMorton();
							jump = 0;
						}
					}
				}
				Mortontry = m_ghosts[idxtry].computeMorton();
				if(Mortontry == Morton && m_ghosts[idxtry].m_level == oct->m_level){
					//Found neighbour of same size
					isghost.push_back(true);
					neighbours.push_back(idxtry);
					return;
				}
				else{
					// Step until the mortontry lower than morton (one idx of distance)
					{
						while(Mortontry < Morton){
							idxtry++;
							if(idxtry > m_ghosts.size()-1){
								idxtry = m_ghosts.size()-1;
								Mortontry = m_ghosts[idxtry].computeMorton();
								break;
							}
							Mortontry = m_ghosts[idxtry].computeMorton();
						}
						while(m_ghosts[idxtry].computeMorton() > Morton){
							idxtry--;
							if(idxtry > m_ghosts.size()-1){
								idxtry = 0;
								Mortontry = m_ghosts[idxtry].computeMorton();
								break;
							}
							Mortontry = m_ghosts[idxtry].computeMorton();
						}
					}
					if(idxtry < m_sizeGhosts){
						if(m_ghosts[idxtry].computeMorton() == Morton && m_ghosts[idxtry].m_level == oct->m_level){
							//Found neighbour of same size
							isghost.push_back(true);
							neighbours.push_back(idxtry);
							return;
						}
						// Compute Last discendent of virtual octant of same size
						Octant last_desc = samesizeoct.buildLastDesc();
						uint64_t Mortonlast = last_desc.computeMorton();
						Mortontry = m_ghosts[idxtry].computeMorton();
						int32_t Dx[3] = {0,0,0};
						int32_t Dxstar[3] = {0,0,0};
						array<int64_t,3> coord = oct->getPeriodicCoord(iface);
						u32array3 coordtry = m_ghosts[idxtry].getCoord();
						array<int64_t,3> coord1 = { {1,1,1} };
						u32array3 coordtry1 = { {1,1,1} };
						uint8_t level = oct->m_level;
						uint8_t leveltry = m_octants[idxtry].getLevel();
						while(Mortontry <= Mortonlast && idxtry < m_sizeGhosts){
							for (int idim=0; idim<m_dim; idim++){
								Dx[idim] 		= int32_t(int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]));
								Dxstar[idim]	= int32_t((cxyz[idim]-1)/2)*(m_ghosts[idxtry].getSize()) + int32_t((cxyz[idim]+1)/2)*size;
								coord1[idim] 	= coord[idim] + size;
								coordtry1[idim] = coordtry[idim] + m_ghosts[idxtry].getSize();
							}
							leveltry = m_ghosts[idxtry].getLevel();

							if (Dx[0] == Dxstar[0] && Dx[1] == Dxstar[1] && Dx[m_dim-1] == Dxstar[m_dim-1]){
								if (leveltry > level){
									if((abs(cxyz[0])*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[1])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[2])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1])))){
										neighbours.push_back(idxtry);
										isghost.push_back(true);
									}
								}
								if (leveltry < level){
									if((abs(cxyz[0])*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[1])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[2])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1])))){
										neighbours.push_back(idxtry);
										isghost.push_back(true);
									}
								}
							}

							idxtry++;
							if(idxtry>m_sizeGhosts-1){
								break;
							}
							Mortontry = m_ghosts[idxtry].computeMorton();
							coordtry = m_ghosts[idxtry].getCoord();
						}
					}
				}

				uint32_t lengthneigh = 0;
				uint32_t sizeneigh = neighbours.size();
				for (idxtry=0; idxtry<sizeneigh; idxtry++){
					lengthneigh += m_ghosts[neighbours[idxtry]].getArea();
				}
				if (lengthneigh < oct->getArea()){
					// Search in octants

						//Build Morton number of virtual neigh of same size
						Octant samesizeoct = oct->computePeriodicOctant(iface);
						Morton = samesizeoct.computeMorton();
						// Search morton in octants
						// If a even face morton is lower than morton of oct, if odd higher
						// ---> can i search only before or after idx in octants
						int32_t jump = (int32_t((noctants)/2+1));
						idxtry = uint32_t(jump);
						while(abs(jump) > 0){
							Mortontry = m_octants[idxtry].computeMorton();
							jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
							idxtry += jump;
							if (idxtry > noctants-1){
								if (jump > 0){
									idxtry = noctants - 1;
									jump = 0;
								}
								else if (jump < 0){
									idxtry = 0;
									jump = 0;
								}
							}
						}
						Mortontry = m_octants[idxtry].computeMorton();
						if(Mortontry == Morton && m_octants[idxtry].m_level == oct->m_level){
							//Found neighbour of same size
							isghost.push_back(false);
							neighbours.push_back(idxtry);
							return;
						}
						else{
							// Step until the mortontry lower than morton (one idx of distance)
							{
								while(Mortontry < Morton){
									idxtry++;
									if(idxtry > noctants-1){
										idxtry = noctants-1;
										Mortontry = m_octants[idxtry].computeMorton();
										break;
									}
									Mortontry = m_octants[idxtry].computeMorton();
								}
								while(Mortontry > Morton){
									idxtry--;
									if(idxtry > noctants-1){
										idxtry = 0;
										Mortontry = m_octants[idxtry].computeMorton();
										break;
									}
									Mortontry = m_octants[idxtry].computeMorton();
								}
							}
							if(Mortontry == Morton && m_octants[idxtry].m_level == oct->m_level){
								//Found neighbour of same size
								isghost.push_back(false);
								neighbours.push_back(idxtry);
								return;
							}
							// Compute Last discendent of virtual octant of same size
							Octant last_desc = samesizeoct.buildLastDesc();
							uint64_t Mortonlast = last_desc.computeMorton();
							Mortontry = m_octants[idxtry].computeMorton();
							int32_t Dx[3] = {0,0,0};
							int32_t Dxstar[3] = {0,0,0};
							array<int64_t,3> coord = oct->getPeriodicCoord(iface);
							u32array3 coordtry = m_octants[idxtry].getCoord();
							array<int64_t,3> coord1 = { {1,1,1} };
							u32array3 coordtry1 = { {1,1,1} };
							uint8_t level = oct->m_level;
							uint8_t leveltry = m_octants[idxtry].getLevel();
							while(Mortontry <= Mortonlast && idxtry < noctants){
								for (int idim=0; idim<m_dim; idim++){
									Dx[idim] 		= int32_t(int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]));
									Dxstar[idim]	= int32_t((cxyz[idim]-1)/2)*(m_octants[idxtry].getSize()) + int32_t((cxyz[idim]+1)/2)*size;
									coord1[idim] 	= coord[idim] + size;
									coordtry1[idim] = coordtry[idim] + m_octants[idxtry].getSize();
								}
								leveltry = m_octants[idxtry].getLevel();

								if (Dx[0] == Dxstar[0] && Dx[1] == Dxstar[1] && Dx[m_dim-1] == Dxstar[m_dim-1]){
									if (leveltry > level){
										if((abs(cxyz[0])*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[1])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[2])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1])))){
											neighbours.push_back(idxtry);
											isghost.push_back(false);
										}
									}
									if (leveltry < level){
										if((abs(cxyz[0])*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[1])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[2])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1])))){
											neighbours.push_back(idxtry);
											isghost.push_back(false);
										}
									}
								}

								idxtry++;
								if(idxtry>noctants-1){
									break;
								}
								Mortontry = m_octants[idxtry].computeMorton();
								coordtry = m_octants[idxtry].getCoord();
							}
							return;
						}
				}
				return;
			}
	}
};

// =================================================================================== //

/*! Finds neighbours of idx-th ghost octant through a periodic iface in vector m_octants.
 * Returns a vector with the index of neighbours in the structure octants.
 * \param[in] oct Pointer to the target ghost octant.
 * \param[in] iface local index of the face.
 * \param[out] neighbours Vector with the local indices of the local octant neighbours.
 */
void
LocalTree::findGhostPeriodicNeighbours(Octant* oct, uint8_t iface, u32vector & neighbours){

	uint64_t  Morton, Mortontry;
	uint32_t  noctants = getNumOctants();
	uint32_t idxtry;
	uint32_t size = oct->getSize();

	int8_t 			cxyz[3] = {0,0,0};
	for (int idim=0; idim<m_dim; idim++){
		cxyz[idim] = m_global.m_normals[iface][idim];
	}

	neighbours.clear();

	// Default if iface is nface<iface<0
	if (iface > m_global.m_nfaces){
		return;
	}

	// Check if octants face is a process boundary
	if (oct->m_info[6+iface] == true){

		//Build Morton number of virtual neigh of same size
		Octant samesizeoct = oct->computePeriodicOctant(iface);
		Morton = samesizeoct.computeMorton();
		// Search morton in octants
		// If a even face morton is lower than morton of oct, if odd higher
		// ---> can i search only before or after idx in octants
		int32_t jump = getNumOctants()/2;
		idxtry = uint32_t(getNumOctants()/2);
		Mortontry = m_octants[idxtry].computeMorton();
//		jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
		while(abs(jump) > 0){

			Mortontry = m_octants[idxtry].computeMorton();
			jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
			idxtry += jump;
			if (idxtry > noctants-1){
				if (jump > 0){
					idxtry = noctants - 1;
					Mortontry = m_octants[idxtry].computeMorton();
					jump = 0;
				}
				else if (jump < 0){
					idxtry = 0;
					Mortontry = m_octants[idxtry].computeMorton();
					jump = 0;
				}
			}
		}
		Mortontry = m_octants[idxtry].computeMorton();
		if(Mortontry == Morton && m_octants[idxtry].m_level == oct->m_level){
			//Found neighbour of same size
			neighbours.push_back(idxtry);
			return;
		}
		else{
			// Step until the mortontry lower than morton (one idx of distance)
			{
				while(Mortontry < Morton){
					idxtry++;
					if(idxtry > noctants-1){
						idxtry = noctants-1;
						Mortontry = m_octants[idxtry].computeMorton();
						break;
					}
					Mortontry = m_octants[idxtry].computeMorton();
				}
				while(Mortontry > Morton){
					idxtry--;
					if(idxtry > noctants-1){
						idxtry = 0;
						Mortontry = m_octants[idxtry].computeMorton();
						break;
					}
					Mortontry = m_octants[idxtry].computeMorton();
				}
			}
			if(Mortontry == Morton && m_octants[idxtry].m_level == oct->m_level){
				//Found neighbour of same size
				neighbours.push_back(idxtry);
				return;
			}
			// Compute Last discendent of virtual octant of same size
			Octant last_desc = samesizeoct.buildLastDesc();
			uint64_t Mortonlast = last_desc.computeMorton();
			Mortontry = m_octants[idxtry].computeMorton();
			int32_t Dx[3] = {0,0,0};
			int32_t Dxstar[3] = {0,0,0};
			array<int64_t,3> coord = oct->getPeriodicCoord(iface);

			u32array3 coordtry = m_octants[idxtry].getCoord();
			array<int64_t,3> coord1 = { {1,1,1} };
			u32array3 coordtry1 = { {1,1,1} };
			uint8_t level = oct->m_level;
			while(Mortontry <= Mortonlast && idxtry < noctants){
				for (int idim=0; idim<m_dim; idim++){
					Dx[idim] 		= int32_t(int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]));
					Dxstar[idim]	= int32_t((cxyz[idim]-1)/2)*(m_octants[idxtry].getSize()) + int32_t((cxyz[idim]+1)/2)*size;
					coord1[idim] 	= coord[idim] + size;
					coordtry1[idim] = coordtry[idim] + m_octants[idxtry].getSize();
				}
				uint8_t leveltry = m_octants[idxtry].getLevel();

				if (Dx[0] == Dxstar[0] && Dx[1] == Dxstar[1] && Dx[m_dim-1] == Dxstar[m_dim-1]){
					if (leveltry > level){
						if((abs(cxyz[0])*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[1])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[2])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1])))){
							neighbours.push_back(idxtry);
						}
					}
					if (leveltry < level){
						if((abs(cxyz[0])*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[1])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[2])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1])))){
							neighbours.push_back(idxtry);
						}
					}
				}

				idxtry++;
				if(idxtry>noctants-1){
					break;
				}
				Mortontry = m_octants[idxtry].computeMorton();
				coordtry = m_octants[idxtry].getCoord();
			}
			return;
		}
	}

};

// =================================================================================== //

/*! Pre-processing for 2:1 balancing of local tree. Check if there are broken families over processes.
 * \param[in] internal Set to true if the interior octants have to be checked.
 */
void
 LocalTree::preBalance21(bool internal){

	Octant 			father(m_dim, m_global.m_maxLevel), lastdesc(m_dim, m_global.m_maxLevel);
	uint64_t 		mortonld;
	uint32_t 		nocts;
	uint32_t 		idx, idx2, idx0, last_idx;
	uint32_t 		idx1_gh, idx2_gh;
	int8_t 			marker;
	uint8_t 		nbro;

	//------------------------------------------ //
	// Initialization

	nbro = 0;
	idx2_gh = idx0 = 0;
	idx1_gh=0;

	nocts   = m_octants.size();
	m_sizeGhosts = m_ghosts.size();
	last_idx=nocts-1;

	//Clean index of ghost brothers in case of coarsening a broken family
	m_lastGhostBros.clear();

	// Set index for start and end check for ghosts
	if (m_ghosts.size()){
		while(m_ghosts[idx2_gh].computeMorton() <= m_lastDesc.computeMorton()){
			idx2_gh++;
			if (idx2_gh > m_sizeGhosts-1) break;
		}
		idx2_gh = min((m_sizeGhosts-1), idx2_gh);

		while(m_ghosts[idx1_gh].computeMorton() <= m_octants[0].computeMorton()){
			idx1_gh++;
			if (idx1_gh > m_sizeGhosts-1) break;
		}
		idx1_gh-=1;
		if (idx1_gh > m_sizeGhosts-1) idx1_gh=0;
	}

	// End on ghosts
	if (m_ghosts.size() && nocts > 0){
		if (m_ghosts[idx1_gh].buildFather()==m_octants[0].buildFather()){
			father = m_ghosts[idx1_gh].buildFather();
			nbro = 0;
			idx = idx1_gh;
			marker = m_ghosts[idx].getMarker();
			while(marker < 0 && m_ghosts[idx].buildFather() == father){
				nbro++;
				if (idx==0)
					break;
				idx--;
				marker = m_ghosts[idx].getMarker();
			}
			idx = 0;
			while(idx<nocts && m_octants[idx].buildFather() == father){
				if(m_octants[idx].getMarker()<0)
					nbro++;
				idx++;
				if(idx==nocts)
					break;
			}
			if (nbro != m_global.m_nchildren && idx!=nocts-1){
				for(uint32_t ii=0; ii<idx; ii++){
					if (m_octants[ii].getMarker()<0){
						m_octants[ii].setMarker(0);
						m_octants[ii].m_info[15]=true;
					}
				}
			}
		}

		if (m_ghosts[idx2_gh].buildFather()==m_octants[nocts-1].buildFather()){
			father = m_ghosts[idx2_gh].buildFather();
			nbro = 0;
			idx = idx2_gh;
			marker = m_ghosts[idx].getMarker();
			while(marker < 0 && m_ghosts[idx].buildFather() == father){

				//Add ghost index to structure for mapper in case of coarsening a broken family
				m_lastGhostBros.push_back(idx);

				nbro++;
				idx++;
				if(idx == m_sizeGhosts){
					break;
				}
				marker = m_ghosts[idx].getMarker();
			}
			idx = nocts-1;
			while(m_octants[idx].buildFather() == father ){
				if (m_octants[idx].getMarker()<0)
					nbro++;
				if (idx==0)
					break;
				idx--;
			}
			last_idx=idx;
			if (nbro != m_global.m_nchildren && idx!=nocts-1){
				for(uint32_t ii=idx+1; ii<nocts; ii++){
					if (m_octants[ii].getMarker()<0){
						m_octants[ii].setMarker(0);
						m_octants[ii].m_info[15]=true;
					}
				}
				//Clean ghost index to structure for mapper in case of coarsening a broken family
				m_lastGhostBros.clear();
			}
		}
	}

	// Check first internal octants
	if (internal){
		father = m_octants[0].buildFather();
		lastdesc = father.buildLastDesc();
		mortonld = lastdesc.computeMorton();
		nbro = 0;
		for (idx=0; idx<m_global.m_nchildren; idx++){
			if (idx<nocts){
				// Check if family is complete or to be checked in the internal loop (some brother refined)
				if (m_octants[idx].computeMorton() <= mortonld){
					nbro++;
				}
			}
		}
		if (nbro != m_global.m_nchildren)
			idx0 = nbro;

		// Check and coarse internal octants
		for (idx=idx0; idx<nocts; idx++){
			if(m_octants[idx].getMarker() < 0 && m_octants[idx].getLevel() > 0){
				nbro = 0;
				father = m_octants[idx].buildFather();
				// Check if family is to be coarsened
				for (idx2=idx; idx2<idx+m_global.m_nchildren; idx2++){
					if (idx2<nocts){
						if(m_octants[idx2].getMarker() < 0 && m_octants[idx2].buildFather() == father){
							nbro++;
						}
					}
				}
				if (nbro == m_global.m_nchildren){
					idx = idx2-1;
				}
				else{
					if (idx<=last_idx){
						m_octants[idx].setMarker(0);
						m_octants[idx].m_info[15]=true;
					}
				}
			}
		}
	}
};

// =================================================================================== //

/*! Pre-processing for 2:1 balancing of local tree. Check if there are broken families over processes.
 * \param[out] newmodified Vector of indices of interior octants checked and whose marker is modified.
 */
void
LocalTree::preBalance21(u32vector& newmodified){

	Octant 				father(m_dim, m_global.m_maxLevel), lastdesc(m_dim, m_global.m_maxLevel);
	uint64_t 			mortonld;
	uint32_t 			nocts;
	uint32_t 			idx, idx2, idx0, last_idx;
	uint32_t 			idx1_gh, idx2_gh;
	int8_t 			marker;
	uint8_t 			nbro;

	//------------------------------------------ //
	// Initialization

	nbro = 0;
	idx2_gh = idx0 = 0;
	idx1_gh=0;

	nocts   = m_octants.size();
	m_sizeGhosts = m_ghosts.size();
	last_idx=nocts-1;

	//Clean index of ghost brothers in case of coarsening a broken family
	m_lastGhostBros.clear();

	// Set index for start and end check for ghosts
	if (m_ghosts.size()){
		while(m_ghosts[idx2_gh].computeMorton() <= m_lastDesc.computeMorton()){
			idx2_gh++;
			if (idx2_gh > m_sizeGhosts-1) break;
		}
		idx2_gh = min((m_sizeGhosts-1), idx2_gh);

		while(m_ghosts[idx1_gh].computeMorton() <= m_octants[0].computeMorton()){
			idx1_gh++;
			if (idx1_gh > m_sizeGhosts-1) break;
		}
		idx1_gh-=1;
		if (idx1_gh > m_sizeGhosts-1) idx1_gh = 0;
	}

	// End on ghosts
	if (m_ghosts.size() && nocts > 0){
		if (m_ghosts[idx1_gh].buildFather()==m_octants[0].buildFather()){
			father = m_ghosts[idx1_gh].buildFather();
			nbro = 0;
			idx = idx1_gh;
			marker = m_ghosts[idx].getMarker();
			while(marker < 0 && m_ghosts[idx].buildFather() == father){

				//Add ghost index to structure for mapper in case of coarsening a broken family
				m_lastGhostBros.push_back(idx);

				nbro++;
				if (idx==0)
					break;
				idx--;
				marker = m_ghosts[idx].getMarker();
			}
			idx = 0;
			while(idx<nocts && m_octants[idx].buildFather() == father){
				if (m_octants[idx].getMarker()<0)
					nbro++;
				idx++;
				if(idx==nocts)
					break;
			}
			if (nbro != m_global.m_nchildren && idx!=nocts-1){
				for(uint32_t ii=0; ii<idx; ii++){
					if (m_octants[ii].getMarker()<0){
						m_octants[ii].setMarker(0);
						m_octants[ii].m_info[15]=true;
						newmodified.push_back(ii);
					}
				}
				//Clean index of ghost brothers in case of coarsening a broken family
				m_lastGhostBros.clear();
			}
		}

		if (m_ghosts[idx2_gh].buildFather()==m_octants[nocts-1].buildFather()){
			father = m_ghosts[idx2_gh].buildFather();
			nbro = 0;
			idx = idx2_gh;
			marker = m_ghosts[idx].getMarker();
			while(marker < 0 && m_ghosts[idx].buildFather() == father){
				nbro++;
				idx++;
				if(idx == m_sizeGhosts){
					break;
				}
				marker = m_ghosts[idx].getMarker();
			}
			idx = nocts-1;
			while(m_octants[idx].buildFather() == father){
				if (m_octants[idx].getMarker()<0)
					nbro++;
				idx--;
				if (idx==0)
					break;
			}
			last_idx=idx;
			if (nbro != m_global.m_nchildren && idx!=nocts-1){
				for(uint32_t ii=idx+1; ii<nocts; ii++){
					if (m_octants[ii].getMarker()<0){
						m_octants[ii].setMarker(0);
						m_octants[ii].m_info[15]=true;
						newmodified.push_back(ii);
					}
				}
			}
		}
	}

	// Check first internal octants
	father = m_octants[0].buildFather();
	lastdesc = father.buildLastDesc();
	mortonld = lastdesc.computeMorton();
	nbro = 0;
	for (idx=0; idx<m_global.m_nchildren; idx++){
		// Check if family is complete or to be checked in the internal loop (some brother refined)
		if (idx<nocts){
			if (m_octants[idx].computeMorton() <= mortonld){
				nbro++;
			}
		}
	}
	if (nbro != m_global.m_nchildren)
		idx0 = nbro;

	// Check and coarse internal octants
	for (idx=idx0; idx<nocts; idx++){
		if(m_octants[idx].getMarker() < 0 && m_octants[idx].getLevel() > 0){
			nbro = 0;
			father = m_octants[idx].buildFather();
			// Check if family is to be coarsened
			for (idx2=idx; idx2<idx+m_global.m_nchildren; idx2++){
				if (idx2<nocts){
					if(m_octants[idx2].getMarker() < 0 && m_octants[idx2].buildFather() == father){
						nbro++;
					}
				}
			}
			if (nbro == m_global.m_nchildren){
				idx = idx2-1;
			}
			else{
				if (idx<=last_idx){
					m_octants[idx].setMarker(0);
					m_octants[idx].m_info[15]=true;
					newmodified.push_back(idx);
				}
			}
		}
	}
};

// =================================================================================== //

/*! 2:1 balancing on level a local tree already adapted (balance only the octants with info[15] = false) (refinement wins!)
 * \param[in] doInterior Set to false if the interior octants are already balanced.
 * \return True if balanced done with some markers modification.
 */
bool
LocalTree::localBalance(bool doInterior){

	uint32_t			sizeneigh, modsize;
	u32vector		 	neigh;
	u32vector		 	modified, newmodified;
	uint32_t 			i, idx;
	uint8_t				iface, iedge, inode;
	int8_t				targetmarker;
	vector<bool> 		isghost;
	bool				Bdone = false;
	bool				Bedge = ((m_balanceCodim>1) && (m_dim==3));
	bool				Bnode = (m_balanceCodim==m_dim);

	octvector::iterator 	obegin, oend, it;
	u32vector::iterator 	ibegin, iend, iit;

	//If interior octants have to be balanced
	if(doInterior){
		// First loop on the octants
		obegin = m_octants.begin();
		oend = m_octants.end();
		idx = 0;
		for (it=obegin; it!=oend; it++){
			if (it->getBalance() && (it->getMarker() != 0 || it->m_info[15]) ){
				targetmarker = min(m_global.m_maxLevel, int8_t(m_octants[idx].getLevel() + m_octants[idx].getMarker()));

				//Balance through faces
				for (iface=0; iface<m_global.m_nfaces; iface++){
//					if(!it->getBound(iface)){
						findNeighbours(idx, iface, neigh, isghost);
						sizeneigh = neigh.size();
						for(i=0; i<sizeneigh; i++){
							if (!isghost[i]){
								{
									if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) > (targetmarker + 1) ){
										m_octants[idx].setMarker(m_octants[neigh[i]].getLevel()+m_octants[neigh[i]].getMarker()-1-m_octants[idx].getLevel());
										m_octants[idx].m_info[15] = true;
										modified.push_back(idx);
										Bdone = true;
									}
									else if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
										m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
										m_octants[neigh[i]].m_info[15] = true;
										modified.push_back(neigh[i]);
										Bdone = true;
									}
								};
							}
							else{
								{
									if((m_ghosts[neigh[i]].getLevel() + m_ghosts[neigh[i]].getMarker()) > (targetmarker + 1) ){
										m_octants[idx].setMarker(m_ghosts[neigh[i]].getLevel()+m_ghosts[neigh[i]].getMarker()-1-m_octants[idx].getLevel());
										m_octants[idx].m_info[15] = true;
										modified.push_back(idx);
										Bdone = true;
									}
								};
							}
						}
//					}
						targetmarker = min(m_global.m_maxLevel, int8_t(m_octants[idx].getLevel() + m_octants[idx].getMarker()));
				}

				if (Bedge){
					//Balance through edges
					for (iedge=0; iedge<m_global.m_nedges; iedge++){
						//if(!it->getBound(m_global.m_edgeFace[iedge][0]) && !it->getBound(m_global.m_edgeFace[iedge][1])){
							findEdgeNeighbours(idx, iedge, neigh, isghost);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if (!isghost[i]){
									{
										if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) > (targetmarker + 1) ){
											m_octants[idx].setMarker(m_octants[neigh[i]].getLevel()+m_octants[neigh[i]].getMarker()-1-m_octants[idx].getLevel());
											m_octants[idx].m_info[15] = true;
											modified.push_back(idx);
											Bdone = true;
										}
										else if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
											m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
											m_octants[neigh[i]].m_info[15] = true;
											modified.push_back(neigh[i]);
											Bdone = true;
										}
									};
								}
								else{
									if((m_ghosts[neigh[i]].getLevel() + m_ghosts[neigh[i]].getMarker()) > (targetmarker + 1) ){
										m_octants[idx].setMarker(m_ghosts[neigh[i]].getLevel()+m_ghosts[neigh[i]].getMarker()-1-m_octants[idx].getLevel());
										m_octants[idx].m_info[15] = true;
										modified.push_back(idx);
										Bdone = true;
									}
								}
							}
						//}
							targetmarker = min(m_global.m_maxLevel, int8_t(m_octants[idx].getLevel() + m_octants[idx].getMarker()));
					}
				}

				if (Bnode){
					//Balance through nodes
					for (inode=0; inode<m_global.m_nnodes; inode++){
						//if(!it->getBound(m_global.m_nodeFace[inode][0]) && !it->getBound(m_global.m_nodeFace[inode][1]) && !it->getBound(m_global.m_nodeFace[inode][m_dim-1])){
							findNodeNeighbours(idx, inode, neigh, isghost);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if (!isghost[i]){
									{
										if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) > (targetmarker + 1) ){
											m_octants[idx].setMarker(m_octants[neigh[i]].getLevel()+m_octants[neigh[i]].getMarker()-1-m_octants[idx].getLevel());
											m_octants[idx].m_info[15] = true;
											modified.push_back(idx);
											Bdone = true;
										}
										else if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
											m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
											m_octants[neigh[i]].m_info[15] = true;
											modified.push_back(neigh[i]);
											Bdone = true;
										}
									};
								}
								else{
									if((m_ghosts[neigh[i]].getLevel() + m_ghosts[neigh[i]].getMarker()) > (targetmarker + 1) ){
										m_octants[idx].setMarker(m_ghosts[neigh[i]].getLevel()+m_ghosts[neigh[i]].getMarker()-1-m_octants[idx].getLevel());
										m_octants[idx].m_info[15] = true;
										modified.push_back(idx);
										Bdone = true;
									}
								}
							}
						//}
							targetmarker = min(m_global.m_maxLevel, int8_t(m_octants[idx].getLevel() + m_octants[idx].getMarker()));
					}
				}

			}
			idx++;
		}
		// Loop on ghost octants (influence over interior borders)
		obegin = m_ghosts.begin();
		oend = m_ghosts.end();
		idx = 0;
		for (it=obegin; it!=oend; it++){
			if (!it->getNotBalance() && (it->getMarker() != 0 || it->m_info[15]) ){
				targetmarker = min(m_global.m_maxLevel, int8_t(it->getLevel()+it->getMarker()));

				//Balance through faces
				for (iface=0; iface<m_global.m_nfaces; iface++){
					if(it->getPbound(iface) == true){
						neigh.clear();
						findGhostNeighbours(idx, iface, neigh);
						sizeneigh = neigh.size();
						for(i=0; i<sizeneigh; i++){
							if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
								m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
								m_octants[neigh[i]].m_info[15] = true;
								modified.push_back(neigh[i]);
								Bdone = true;
							}
						}
					}
					targetmarker = min(m_global.m_maxLevel, int8_t(it->getLevel()+it->getMarker()));
				}

				if (Bedge){
					//Balance through edges
					for (iedge=0; iedge<m_global.m_nedges; iedge++){
						//if(it->getPbound(m_global.m_edgeFace[iedge][0]) == true || it->getPbound(m_global.m_edgeFace[iedge][1]) == true){
							neigh.clear();
							findGhostEdgeNeighbours(idx, iedge, neigh);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
									m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
									m_octants[neigh[i]].m_info[15] = true;
									modified.push_back(neigh[i]);
									Bdone = true;
								}
							}
						//}
							targetmarker = min(m_global.m_maxLevel, int8_t(it->getLevel()+it->getMarker()));
					}
				}

				if (Bnode){
					//Balance through nodes
					for (inode=0; inode<m_global.m_nnodes; inode++){
						//if(it->getPbound(m_global.m_nodeFace[inode][0]) == true || it->getPbound(m_global.m_nodeFace[inode][1]) == true || it->getPbound(m_global.m_nodeFace[inode][m_dim-1]) == true){
							neigh.clear();
							findGhostNodeNeighbours(idx, inode, neigh);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
									m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
									m_octants[neigh[i]].m_info[15] = true;
									modified.push_back(neigh[i]);
									Bdone = true;
								}
							}
						//}
							targetmarker = min(m_global.m_maxLevel, int8_t(it->getLevel()+it->getMarker()));
					}
				}

			}
			idx++;
		}

		// While loop for iterative balancing
		u32vector().swap(newmodified);
		modsize = modified.size();
		while(modsize!=0){
			ibegin = modified.begin();
			iend = modified.end();
			for (iit=ibegin; iit!=iend; iit++){
				idx = *iit;
				if (!m_octants[idx].getNotBalance()){
					targetmarker = min(m_global.m_maxLevel, int8_t(m_octants[idx].getLevel()+m_octants[idx].getMarker()));

					//Balance through faces
					for (iface=0; iface<m_global.m_nfaces; iface++){
						if(!m_octants[idx].getPbound(iface)){
							findNeighbours(idx, iface, neigh, isghost);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if (!isghost[i]){
									{
										if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
											m_octants[idx].setMarker(m_octants[neigh[i]].getLevel()+m_octants[neigh[i]].getMarker()-m_octants[idx].getLevel()-1);
											m_octants[idx].m_info[15] = true;
											newmodified.push_back(idx);
											Bdone = true;
										}
										else if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
											m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
											m_octants[neigh[i]].m_info[15] = true;
											newmodified.push_back(neigh[i]);
											Bdone = true;
										}
									};
								}
							}
						}
						targetmarker = min(m_global.m_maxLevel, int8_t(m_octants[idx].getLevel()+m_octants[idx].getMarker()));
					}

					if (Bedge){
						//Balance through edges
						for (iedge=0; iedge<m_global.m_nedges; iedge++){
							//if(!m_octants[idx].getPbound(m_global.m_edgeFace[iedge][0]) || !m_octants[idx].getPbound(m_global.m_edgeFace[iedge][1])){
								findEdgeNeighbours(idx, iedge, neigh, isghost);
								sizeneigh = neigh.size();
								for(i=0; i<sizeneigh; i++){
									if (!isghost[i]){
										{
											if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
												m_octants[idx].setMarker(m_octants[neigh[i]].getLevel()+m_octants[neigh[i]].getMarker()-m_octants[idx].getLevel()-1);
												m_octants[idx].m_info[15] = true;
												newmodified.push_back(idx);
												Bdone = true;
											}
											else if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
												m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
												m_octants[neigh[i]].m_info[15] = true;
												newmodified.push_back(neigh[i]);
												Bdone = true;
											}
										};
									}
								}
							//}
								targetmarker = min(m_global.m_maxLevel, int8_t(m_octants[idx].getLevel()+m_octants[idx].getMarker()));
						}
					}

					if (Bnode){
						//Balance through nodes
						for (inode=0; inode<m_global.m_nnodes; inode++){
							//if(!m_octants[idx].getPbound(m_global.m_nodeFace[inode][0]) || !m_octants[idx].getPbound(m_global.m_nodeFace[inode][1]) || !m_octants[idx].getPbound(m_global.m_nodeFace[inode][m_dim-1])){
								findNodeNeighbours(idx, inode, neigh, isghost);
								sizeneigh = neigh.size();
								for(i=0; i<sizeneigh; i++){
									if (!isghost[i]){
										{
											if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
												m_octants[idx].setMarker(m_octants[neigh[i]].getLevel()+m_octants[neigh[i]].getMarker()-m_octants[idx].getLevel()-1);
												m_octants[idx].m_info[15] = true;
												newmodified.push_back(idx);
												Bdone = true;
											}
											else if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
												m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
												m_octants[neigh[i]].m_info[15] = true;
												newmodified.push_back(neigh[i]);
												Bdone = true;
											}
										};
									}
								}
							//}
								targetmarker = min(m_global.m_maxLevel, int8_t(m_octants[idx].getLevel()+m_octants[idx].getMarker()));
						}
					}

				}
			}
			preBalance21(newmodified);
			u32vector().swap(modified);
			swap(modified,newmodified);
			modsize = modified.size();
			u32vector().swap(newmodified);
		}// end while

	}
	else{

		// Loop on ghost octants (influence over interior borders)
		obegin = m_ghosts.begin();
		oend = m_ghosts.end();
		idx = 0;
		for (it=obegin; it!=oend; it++){
			if (!it->getNotBalance() && it->m_info[15]){
				targetmarker = min(m_global.m_maxLevel, int8_t(it->getLevel()+it->getMarker()));

				//Balance through faces
				for (iface=0; iface<m_global.m_nfaces; iface++){
					if(it->getPbound(iface) == true){
						neigh.clear();
						findGhostNeighbours(idx, iface, neigh);
						sizeneigh = neigh.size();
						for(i=0; i<sizeneigh; i++){
							if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
								m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
								m_octants[neigh[i]].m_info[15] = true;
								modified.push_back(neigh[i]);
								Bdone = true;
							}
						}
					}
					targetmarker = min(m_global.m_maxLevel, int8_t(it->getLevel()+it->getMarker()));
				}

				if (Bedge){
					//Balance through edges
					for (iedge=0; iedge<m_global.m_nedges; iedge++){
						//if(it->getPbound(m_global.m_edgeFace[iedge][0]) == true || it->getPbound(m_global.m_edgeFace[iedge][1]) == true){
							neigh.clear();
							findGhostEdgeNeighbours(idx, iedge, neigh);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
									m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
									m_octants[neigh[i]].m_info[15] = true;
									modified.push_back(neigh[i]);
									Bdone = true;
								}
							}
						//}
							targetmarker = min(m_global.m_maxLevel, int8_t(it->getLevel()+it->getMarker()));
					}
				}

				if (Bnode){
					//Balance through nodes
					for (inode=0; inode<m_global.m_nnodes; inode++){
						//if(it->getPbound(m_global.m_nodeFace[inode][0]) == true || it->getPbound(m_global.m_nodeFace[inode][1]) == true || it->getPbound(m_global.m_nodeFace[inode][m_dim-1]) == true){
							neigh.clear();
							findGhostNodeNeighbours(idx, inode, neigh);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
									m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
									m_octants[neigh[i]].m_info[15] = true;
									modified.push_back(neigh[i]);
									Bdone = true;
								}
							}
						//}
							targetmarker = min(m_global.m_maxLevel, int8_t(it->getLevel()+it->getMarker()));
					}
				}

			}
			idx++;
		}

		// While loop for iterative balancing
		u32vector().swap(newmodified);
		modsize = modified.size();
		while(modsize!=0){
			ibegin = modified.begin();
			iend = modified.end();
			for (iit=ibegin; iit!=iend; iit++){
				idx = *iit;
				if (!m_octants[idx].getNotBalance()){
					targetmarker = min(m_global.m_maxLevel, int8_t(m_octants[idx].getLevel()+m_octants[idx].getMarker()));

					//Balance through faces
					for (iface=0; iface<m_global.m_nfaces; iface++){
						if(!m_octants[idx].getPbound(iface)){
							findNeighbours(idx, iface, neigh, isghost);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if (!isghost[i]){
									{
										if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
											m_octants[idx].setMarker(m_octants[neigh[i]].getLevel()+m_octants[neigh[i]].getMarker()-m_octants[idx].getLevel()-1);
											m_octants[idx].m_info[15] = true;
											newmodified.push_back(idx);
											Bdone = true;
										}
										else if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
											m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
											m_octants[neigh[i]].m_info[15] = true;
											newmodified.push_back(neigh[i]);
											Bdone = true;
										}
									};
								}
							}
						}
						targetmarker = min(m_global.m_maxLevel, int8_t(m_octants[idx].getLevel()+m_octants[idx].getMarker()));
					}

					if (Bedge){
						//Balance through edges
						for (iedge=0; iedge<m_global.m_nedges; iedge++){
							//if(!m_octants[idx].getPbound(m_global.m_edgeFace[iedge][0]) || !m_octants[idx].getPbound(m_global.m_edgeFace[iedge][1])){
								findEdgeNeighbours(idx, iedge, neigh, isghost);
								sizeneigh = neigh.size();
								for(i=0; i<sizeneigh; i++){
									if (!isghost[i]){
										{
											if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
												m_octants[idx].setMarker(m_octants[neigh[i]].getLevel()+m_octants[neigh[i]].getMarker()-m_octants[idx].getLevel()-1);
												m_octants[idx].m_info[15] = true;
												newmodified.push_back(idx);
												Bdone = true;
											}
											else if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
												m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
												m_octants[neigh[i]].m_info[15] = true;
												newmodified.push_back(neigh[i]);
												Bdone = true;
											}
										};
									}
								}
							//}
								targetmarker = min(m_global.m_maxLevel, int8_t(m_octants[idx].getLevel()+m_octants[idx].getMarker()));
						}
					}

					if (Bnode){
						//Balance through nodes
						for (inode=0; inode<m_global.m_nnodes; inode++){
							//if(!m_octants[idx].getPbound(m_global.m_nodeFace[inode][0]) || !m_octants[idx].getPbound(m_global.m_nodeFace[inode][1]) || !m_octants[idx].getPbound(m_global.m_nodeFace[inode][m_dim-1])){
								findNodeNeighbours(idx, inode, neigh, isghost);
								sizeneigh = neigh.size();
								for(i=0; i<sizeneigh; i++){
									if (!isghost[i]){
										{
											if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
												m_octants[idx].setMarker(m_octants[neigh[i]].getLevel()+m_octants[neigh[i]].getMarker()-m_octants[idx].getLevel()-1);
												m_octants[idx].m_info[15] = true;
												newmodified.push_back(idx);
												Bdone = true;
											}
											else if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
												m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
												m_octants[neigh[i]].m_info[15] = true;
												newmodified.push_back(neigh[i]);
												Bdone = true;
											}
										};
									}
								}
							//}
								targetmarker = min(m_global.m_maxLevel, int8_t(m_octants[idx].getLevel()+m_octants[idx].getMarker()));
						}
					}

				}
			}
			preBalance21(newmodified);
			u32vector().swap(modified);
			swap(modified,newmodified);
			modsize = modified.size();
			u32vector().swap(newmodified);
		}// end while
		obegin = oend = m_octants.end();
		ibegin = iend = modified.end();
	}
	return Bdone;
	// Pay attention : info[15] may be true after local balance for some octants
};

// =================================================================================== //

/*! 2:1 balancing on level a local tree already adapted (balance all the local octants and new or modified ghost octants) (refinement wins!)
 * \param[in] doInterior Set to false if the interior octants are already balanced.
 * \return True if balanced done with some markers modification.
 */
bool
LocalTree::localBalanceAll(bool doInterior){
	// Local variables
	uint32_t			sizeneigh, modsize;
	u32vector		 	neigh;
	u32vector		 	modified, newmodified;
	uint32_t 			i, idx;
	uint8_t				iface, iedge, inode;
	int8_t				targetmarker;
	vector<bool> 		isghost;
	bool				Bdone = false;
	bool				Bedge = ((m_balanceCodim>1) && (m_dim==3));
	bool				Bnode = (m_balanceCodim==m_dim);

	octvector::iterator 	obegin, oend, it;
	u32vector::iterator 	ibegin, iend, iit;


	//If interior octants have to be balanced
	if(doInterior){
		// First loop on the octants
		obegin = m_octants.begin();
		oend = m_octants.end();
		idx = 0;
		for (it=obegin; it!=oend; it++){
			if ((!it->getNotBalance()) && ((it->m_info[15]) || (it->getMarker()!=0) || ((it->getIsNewC()) || (it->getIsNewR())))){
				targetmarker = min(m_global.m_maxLevel, int8_t(m_octants[idx].getLevel() + m_octants[idx].getMarker()));

				//Balance through faces
				for (iface=0; iface<m_global.m_nfaces; iface++){
//					if(!it->getBound(iface)){
						findNeighbours(idx, iface, neigh, isghost);
						sizeneigh = neigh.size();
						for(i=0; i<sizeneigh; i++){
							if (!isghost[i]){
								{
									if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) > (targetmarker + 1) ){
										m_octants[idx].setMarker(m_octants[neigh[i]].getLevel()+m_octants[neigh[i]].getMarker()-1-m_octants[idx].getLevel());
										m_octants[idx].m_info[15] = true;
										modified.push_back(idx);
										Bdone = true;
									}
									else if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
										m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
										m_octants[neigh[i]].m_info[15] = true;
										modified.push_back(neigh[i]);
										Bdone = true;
									}
								};
							}
							else{
								{
									if((m_ghosts[neigh[i]].getLevel() + m_ghosts[neigh[i]].getMarker()) > (targetmarker + 1) ){
										m_octants[idx].setMarker(m_ghosts[neigh[i]].getLevel()+m_ghosts[neigh[i]].getMarker()-1-m_octants[idx].getLevel());
										m_octants[idx].m_info[15] = true;
										modified.push_back(idx);
										Bdone = true;
									}
								};

							}
						}
//					}
						targetmarker = min(m_global.m_maxLevel, int8_t(m_octants[idx].getLevel() + m_octants[idx].getMarker()));
				}

				if (Bedge){
					//Balance through edges
					for (iedge=0; iedge<m_global.m_nedges; iedge++){
						//if(!it->getBound(m_global.m_edgeFace[iedge][0]) && !it->getBound(m_global.m_edgeFace[iedge][1])){
							findEdgeNeighbours(idx, iedge, neigh, isghost);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if (!isghost[i]){
									{
										if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) > (targetmarker + 1) ){
											m_octants[idx].setMarker(m_octants[neigh[i]].getLevel()+m_octants[neigh[i]].getMarker()-1-m_octants[idx].getLevel());
											m_octants[idx].m_info[15] = true;
											modified.push_back(idx);
											Bdone = true;
										}
										else if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
											m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
											m_octants[neigh[i]].m_info[15] = true;
											modified.push_back(neigh[i]);
											Bdone = true;
										}
									};
								}
								else{
									if((m_ghosts[neigh[i]].getLevel() + m_ghosts[neigh[i]].getMarker()) > (targetmarker + 1) ){
										m_octants[idx].setMarker(m_ghosts[neigh[i]].getLevel()+m_ghosts[neigh[i]].getMarker()-1-m_octants[idx].getLevel());
										m_octants[idx].m_info[15] = true;
										modified.push_back(idx);
										Bdone = true;
									}
								}
							}
						//}
							targetmarker = min(m_global.m_maxLevel, int8_t(m_octants[idx].getLevel() + m_octants[idx].getMarker()));
					}
				}

				if (Bnode){
					//Balance through nodes
					for (inode=0; inode<m_global.m_nnodes; inode++){
						//if(!it->getBound(m_global.m_nodeFace[inode][0]) && !it->getBound(m_global.m_nodeFace[inode][1]) && !it->getBound(m_global.m_nodeFace[inode][m_dim-1])){
							findNodeNeighbours(idx, inode, neigh, isghost);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if (!isghost[i]){
									{
										if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) > (targetmarker + 1) ){
											m_octants[idx].setMarker(m_octants[neigh[i]].getLevel()+m_octants[neigh[i]].getMarker()-1-m_octants[idx].getLevel());
											m_octants[idx].m_info[15] = true;
											modified.push_back(idx);
											Bdone = true;
										}
										else if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
											m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
											m_octants[neigh[i]].m_info[15] = true;
											modified.push_back(neigh[i]);
											Bdone = true;
										}
									};
								}
								else{
									if((m_ghosts[neigh[i]].getLevel() + m_ghosts[neigh[i]].getMarker()) > (targetmarker + 1) ){
										m_octants[idx].setMarker(m_ghosts[neigh[i]].getLevel()+m_ghosts[neigh[i]].getMarker()-1-m_octants[idx].getLevel());
										m_octants[idx].m_info[15] = true;
										modified.push_back(idx);
										Bdone = true;
									}
								}
							}
						//}
							targetmarker = min(m_global.m_maxLevel, int8_t(m_octants[idx].getLevel() + m_octants[idx].getMarker()));
					}
				}

			}
			idx++;
		}
		// Loop on ghost octants (influence over interior borders)
		obegin = m_ghosts.begin();
		oend = m_ghosts.end();
		idx = 0;
		for (it=obegin; it!=oend; it++){
			if (!it->getNotBalance() && (it->m_info[15] || (it->getIsNewC() || it->getIsNewR()))){
				targetmarker = min(m_global.m_maxLevel, int8_t(it->getLevel()+it->getMarker()));

				//Balance through faces
				for (iface=0; iface<m_global.m_nfaces; iface++){
					if(it->getPbound(iface) == true){
						neigh.clear();
						findGhostNeighbours(idx, iface, neigh);
						sizeneigh = neigh.size();
						for(i=0; i<sizeneigh; i++){
							if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
								m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
								m_octants[neigh[i]].m_info[15] = true;
								modified.push_back(neigh[i]);
								Bdone = true;
							}
						}
					}
					targetmarker = min(m_global.m_maxLevel, int8_t(it->getLevel()+it->getMarker()));
				}

				if (Bedge){
					//Balance through edges
					for (iedge=0; iedge<m_global.m_nedges; iedge++){
						//if(it->getPbound(m_global.m_edgeFace[iedge][0]) == true || it->getPbound(m_global.m_edgeFace[iedge][1]) == true){
							neigh.clear();
							findGhostEdgeNeighbours(idx, iedge, neigh);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
									m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
									m_octants[neigh[i]].m_info[15] = true;
									modified.push_back(neigh[i]);
									Bdone = true;
								}
							}
						//}
							targetmarker = min(m_global.m_maxLevel, int8_t(it->getLevel()+it->getMarker()));
					}
				}

				if (Bnode){
					//Balance through nodes
					for (inode=0; inode<m_global.m_nnodes; inode++){
						//if(it->getPbound(m_global.m_nodeFace[inode][0]) == true || it->getPbound(m_global.m_nodeFace[inode][1]) == true || it->getPbound(m_global.m_nodeFace[inode][m_dim-1]) == true){
							neigh.clear();
							findGhostNodeNeighbours(idx, inode, neigh);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
									m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
									m_octants[neigh[i]].m_info[15] = true;
									modified.push_back(neigh[i]);
									Bdone = true;
								}
							}
						//}
							targetmarker = min(m_global.m_maxLevel, int8_t(it->getLevel()+it->getMarker()));
					}
				}

			}
			idx++;
		}

		// While loop for iterative balancing
		u32vector().swap(newmodified);
		modsize = modified.size();
		while(modsize!=0){
			ibegin = modified.begin();
			iend = modified.end();
			for (iit=ibegin; iit!=iend; iit++){
				idx = *iit;
				if (!m_octants[idx].getNotBalance()){
					targetmarker = min(m_global.m_maxLevel, int8_t(m_octants[idx].getLevel()+m_octants[idx].getMarker()));

					//Balance through faces
					for (iface=0; iface<m_global.m_nfaces; iface++){
						if(!m_octants[idx].getPbound(iface)){
							findNeighbours(idx, iface, neigh, isghost);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if (!isghost[i]){
									{
										if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
											m_octants[idx].setMarker(m_octants[neigh[i]].getLevel()+m_octants[neigh[i]].getMarker()-m_octants[idx].getLevel()-1);
											m_octants[idx].m_info[15] = true;
											newmodified.push_back(idx);
											Bdone = true;
										}
										else if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
											m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
											m_octants[neigh[i]].m_info[15] = true;
											newmodified.push_back(neigh[i]);
											Bdone = true;
										}
									};
								}
							}
						}
						targetmarker = min(m_global.m_maxLevel, int8_t(m_octants[idx].getLevel()+m_octants[idx].getMarker()));
					}

					if (Bedge){
						//Balance through edges
						for (iedge=0; iedge<m_global.m_nedges; iedge++){
							//if(!m_octants[idx].getPbound(m_global.m_edgeFace[iedge][0]) || !m_octants[idx].getPbound(m_global.m_edgeFace[iedge][1])){
								findEdgeNeighbours(idx, iedge, neigh, isghost);
								sizeneigh = neigh.size();
								for(i=0; i<sizeneigh; i++){
									if (!isghost[i]){
										{
											if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
												m_octants[idx].setMarker(m_octants[neigh[i]].getLevel()+m_octants[neigh[i]].getMarker()-m_octants[idx].getLevel()-1);
												m_octants[idx].m_info[15] = true;
												newmodified.push_back(idx);
												Bdone = true;
											}
											else if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
												m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
												m_octants[neigh[i]].m_info[15] = true;
												newmodified.push_back(neigh[i]);
												Bdone = true;
											}
										};
									}
								}
							//}
								targetmarker = min(m_global.m_maxLevel, int8_t(m_octants[idx].getLevel()+m_octants[idx].getMarker()));
						}
					}

					if (Bnode){
						//Balance through nodes
						for (inode=0; inode<m_global.m_nnodes; inode++){
							//if(!m_octants[idx].getPbound(m_global.m_nodeFace[inode][0]) || !m_octants[idx].getPbound(m_global.m_nodeFace[inode][1]) || !m_octants[idx].getPbound(m_global.m_nodeFace[inode][m_dim-1])){
								findNodeNeighbours(idx, inode, neigh, isghost);
								sizeneigh = neigh.size();
								for(i=0; i<sizeneigh; i++){
									if (!isghost[i]){
										{
											if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
												m_octants[idx].setMarker(m_octants[neigh[i]].getLevel()+m_octants[neigh[i]].getMarker()-m_octants[idx].getLevel()-1);
												m_octants[idx].m_info[15] = true;
												newmodified.push_back(idx);
												Bdone = true;
											}
											else if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
												m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
												m_octants[neigh[i]].m_info[15] = true;
												newmodified.push_back(neigh[i]);
												Bdone = true;
											}
										};
									}
								}
							//}
								targetmarker = min(m_global.m_maxLevel, int8_t(m_octants[idx].getLevel()+m_octants[idx].getMarker()));
						}
					}

				}
			}
			preBalance21(newmodified);
			u32vector().swap(modified);
			swap(modified,newmodified);
			modsize = modified.size();
			u32vector().swap(newmodified);
		}// end while

	}
	else{

		// Loop on ghost octants (influence over interior borders)
		obegin = m_ghosts.begin();
		oend = m_ghosts.end();
		idx = 0;
		for (it=obegin; it!=oend; it++){
			if (!it->getNotBalance() && (it->m_info[15] || (it->getIsNewC() || it->getIsNewR()))){
				targetmarker = min(m_global.m_maxLevel, int8_t(it->getLevel()+it->getMarker()));

				//Balance through faces
				for (iface=0; iface<m_global.m_nfaces; iface++){
					if(it->getPbound(iface) == true){
						neigh.clear();
						findGhostNeighbours(idx, iface, neigh);
						sizeneigh = neigh.size();
						for(i=0; i<sizeneigh; i++){
							if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
								m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
								m_octants[neigh[i]].m_info[15] = true;
								modified.push_back(neigh[i]);
								Bdone = true;
							}
						}
					}
					targetmarker = min(m_global.m_maxLevel, int8_t(it->getLevel()+it->getMarker()));
				}

				if (Bedge){
					//Balance through edges
					for (iedge=0; iedge<m_global.m_nedges; iedge++){
						//if(it->getPbound(m_global.m_edgeFace[iedge][0]) == true || it->getPbound(m_global.m_edgeFace[iedge][1]) == true){
							neigh.clear();
							findGhostEdgeNeighbours(idx, iedge, neigh);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
									m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
									m_octants[neigh[i]].m_info[15] = true;
									modified.push_back(neigh[i]);
									Bdone = true;
								}
							}
						//}
							targetmarker = min(m_global.m_maxLevel, int8_t(it->getLevel()+it->getMarker()));
					}
				}

				if (Bnode){
					//Balance through nodes
					for (inode=0; inode<m_global.m_nnodes; inode++){
						//if(it->getPbound(m_global.m_nodeFace[inode][0]) == true || it->getPbound(m_global.m_nodeFace[inode][1]) == true || it->getPbound(m_global.m_nodeFace[inode][m_dim-1]) == true){
							neigh.clear();
							findGhostNodeNeighbours(idx, inode, neigh);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
									m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
									m_octants[neigh[i]].m_info[15] = true;
									modified.push_back(neigh[i]);
									Bdone = true;
								}
							}
						//}
							targetmarker = min(m_global.m_maxLevel, int8_t(it->getLevel()+it->getMarker()));
					}
				}

			}
			idx++;
		}

		// While loop for iterative balancing
		u32vector().swap(newmodified);
		modsize = modified.size();
		while(modsize!=0){
			ibegin = modified.begin();
			iend = modified.end();
			for (iit=ibegin; iit!=iend; iit++){
				idx = *iit;
				if (!m_octants[idx].getNotBalance()){
					targetmarker = min(m_global.m_maxLevel, int8_t(m_octants[idx].getLevel()+m_octants[idx].getMarker()));

					//Balance through faces
					for (iface=0; iface<m_global.m_nfaces; iface++){
						if(!m_octants[idx].getPbound(iface)){
							findNeighbours(idx, iface, neigh, isghost);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if (!isghost[i]){
									{
										if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
											m_octants[idx].setMarker(m_octants[neigh[i]].getLevel()+m_octants[neigh[i]].getMarker()-m_octants[idx].getLevel()-1);
											m_octants[idx].m_info[15] = true;
											newmodified.push_back(idx);
											Bdone = true;
										}
										else if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
											m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
											m_octants[neigh[i]].m_info[15] = true;
											newmodified.push_back(neigh[i]);
											Bdone = true;
										}
									};
								}
							}
						}
						targetmarker = min(m_global.m_maxLevel, int8_t(m_octants[idx].getLevel()+m_octants[idx].getMarker()));
					}

					if (Bedge){
						//Balance through edges
						for (iedge=0; iedge<m_global.m_nedges; iedge++){
							//if(!m_octants[idx].getPbound(m_global.m_edgeFace[iedge][0]) || !m_octants[idx].getPbound(m_global.m_edgeFace[iedge][1])){
								findEdgeNeighbours(idx, iedge, neigh, isghost);
								sizeneigh = neigh.size();
								for(i=0; i<sizeneigh; i++){
									if (!isghost[i]){
										{
											if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
												m_octants[idx].setMarker(m_octants[neigh[i]].getLevel()+m_octants[neigh[i]].getMarker()-m_octants[idx].getLevel()-1);
												m_octants[idx].m_info[15] = true;
												newmodified.push_back(idx);
												Bdone = true;
											}
											else if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
												m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
												m_octants[neigh[i]].m_info[15] = true;
												newmodified.push_back(neigh[i]);
												Bdone = true;
											}
										};
									}
								}
							//}
								targetmarker = min(m_global.m_maxLevel, int8_t(m_octants[idx].getLevel()+m_octants[idx].getMarker()));
						}
					}

					if (Bnode){
						//Balance through nodes
						for (inode=0; inode<m_global.m_nnodes; inode++){
							//if(!m_octants[idx].getPbound(m_global.m_nodeFace[inode][0]) || !m_octants[idx].getPbound(m_global.m_nodeFace[inode][1]) || !m_octants[idx].getPbound(m_global.m_nodeFace[inode][m_dim-1])){
								findNodeNeighbours(idx, inode, neigh, isghost);
								sizeneigh = neigh.size();
								for(i=0; i<sizeneigh; i++){
									if (!isghost[i]){
										{
											if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
												m_octants[idx].setMarker(m_octants[neigh[i]].getLevel()+m_octants[neigh[i]].getMarker()-m_octants[idx].getLevel()-1);
												m_octants[idx].m_info[15] = true;
												newmodified.push_back(idx);
												Bdone = true;
											}
											else if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
												m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
												m_octants[neigh[i]].m_info[15] = true;
												newmodified.push_back(neigh[i]);
												Bdone = true;
											}
										};
									}
								}
							//}
								targetmarker = min(m_global.m_maxLevel, int8_t(m_octants[idx].getLevel()+m_octants[idx].getMarker()));
						}
					}

				}
			}
			preBalance21(newmodified);
			u32vector().swap(modified);
			swap(modified,newmodified);
			modsize = modified.size();
			u32vector().swap(newmodified);
		}// end while
		obegin = oend = m_octants.end();
		ibegin = iend = modified.end();
	}
	return Bdone;
	// Pay attention : info[15] may be true after local balance for some octants
};

// =================================================================================== //

/*! Finds neighbours of idx-th octant through iedge in vector m_octants.
 * Returns a vector (empty if iedge is a bound edge) with the index of neighbours
 * in their structure (octants or ghosts) and sets isghost[i] = true if the
 * i-th neighbour is ghost in the local tree.
 * \param[in] idx Local index of the target local octant.
 * \param[in] iedge local index of the edge.
 * \param[out] neighbours Vector with the local indices of the neighbours (size = 0 if boundary edge).
 * \param[out] isghost Vector with the information about the identity of each neighbour (is the i-th neighours a ghost octant?).
 */
void
LocalTree::findEdgeNeighbours(uint32_t idx, uint8_t iedge, u32vector & neighbours, vector<bool> & isghost){

	uint64_t  		Morton, Mortontry;
	uint32_t  		noctants = getNumOctants();
	uint32_t 		idxtry;
	Octant* 		oct = &m_octants[idx];
	uint32_t 		size = oct->getSize();
	uint8_t 		iface1, iface2;
	int32_t 		Dx, Dy, Dz;
	int32_t 		Dxstar,Dystar,Dzstar;

	//Alternative to switch case
	int8_t cx = m_global.m_edgeCoeffs[iedge][0];
	int8_t cy = m_global.m_edgeCoeffs[iedge][1];
	int8_t cz = m_global.m_edgeCoeffs[iedge][2];

	isghost.clear();
	neighbours.clear();

	// Default if iedge is nface<iedge<0
	if (iedge > m_global.m_nfaces*2){
		return;
	}

	// Check if octants edge is a process boundary
	iface1 = m_global.m_edgeFace[iedge][0];
	iface2 = m_global.m_edgeFace[iedge][1];

	// Check if octants edge is a boundary
	if (oct->m_info[iface1] == false && oct->m_info[iface2] == false){

		//Build Morton number of virtual neigh of same size
		Octant samesizeoct(m_dim, oct->m_level, oct->m_x+cx*size, oct->m_y+cy*size, oct->m_z+cz*size, m_global.m_maxLevel);
		Morton = samesizeoct.computeMorton(); //mortonEncode_magicbits(oct->m_x-size,oct->m_y,oct->m_z);

		//SEARCH IN GHOSTS

		if (m_ghosts.size()>0){
			// Search in ghosts
			uint32_t idxghost = uint32_t(m_sizeGhosts/2);
			Octant* octghost = &m_ghosts[idxghost];

			// Search morton in octants
			// If a even face morton is lower than morton of oct, if odd higher
			// ---> can i search only before or after idx in octants
			int32_t jump = int32_t(idxghost/2);
			idxtry = uint32_t(((octghost->computeMorton()<Morton)-(octghost->computeMorton()>Morton))*jump);
			if (idxtry > m_ghosts.size()-1) idxtry = m_ghosts.size()-1;
			while(abs(jump) > 0){
				Mortontry = m_ghosts[idxtry].computeMorton();
				jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
				idxtry += jump;
				if (idxtry > m_ghosts.size()-1){
					if (jump > 0){
						idxtry = m_ghosts.size() - 1;
						jump = 0;
					}
					else if (jump < 0){
						idxtry = 0;
						jump = 0;
					}
				}
			}
			if(m_ghosts[idxtry].computeMorton() == Morton && m_ghosts[idxtry].m_level == oct->m_level){
				//Found neighbour of same size
				isghost.push_back(true);
				neighbours.push_back(idxtry);
				return;
			}
			else{
				// Step until the mortontry lower than morton (one idx of distance)
				{
					while(m_ghosts[idxtry].computeMorton() < Morton){
						idxtry++;
						if(idxtry > m_ghosts.size()-1){
							idxtry = m_ghosts.size()-1;
							break;
						}
					}
					while(m_ghosts[idxtry].computeMorton() > Morton){
						idxtry--;
						if(idxtry > m_ghosts.size()-1){
							idxtry = 0;
							break;
						}
					}
				}
				if(idxtry < m_sizeGhosts){
					if(m_ghosts[idxtry].computeMorton() == Morton && m_ghosts[idxtry].m_level == oct->m_level){
						//Found neighbour of same size
						isghost.push_back(true);
						neighbours.push_back(idxtry);
						return;
					}
					// Compute Last discendent of virtual octant of same size
					Octant last_desc = samesizeoct.buildLastDesc();
					uint64_t Mortonlast = last_desc.computeMorton();
					Mortontry = m_ghosts[idxtry].computeMorton();
					while(Mortontry <= Mortonlast && idxtry < m_ghosts.size()){
						Dx = int32_t(abs(cx))*(-int32_t(oct->m_x) + int32_t(m_ghosts[idxtry].m_x));
						Dy = int32_t(abs(cy))*(-int32_t(oct->m_y) + int32_t(m_ghosts[idxtry].m_y));
						Dz = int32_t(abs(cz))*(-int32_t(oct->m_z) + int32_t(m_ghosts[idxtry].m_z));
						Dxstar = int32_t((cx-1)/2)*(m_ghosts[idxtry].getSize()) + int32_t((cx+1)/2)*size;
						Dystar = int32_t((cy-1)/2)*(m_ghosts[idxtry].getSize()) + int32_t((cy+1)/2)*size;
						Dzstar = int32_t((cz-1)/2)*(m_ghosts[idxtry].getSize()) + int32_t((cz+1)/2)*size;

						uint32_t x0 = oct->m_x;
						uint32_t x1 = x0 + size;
						uint32_t y0 = oct->m_y;
						uint32_t y1 = y0 + size;
						uint32_t z0 = oct->m_z;
						uint32_t z1 = z0 + size;
						uint32_t x0try = m_ghosts[idxtry].m_x;
						uint32_t x1try = x0try + m_ghosts[idxtry].getSize();
						uint32_t y0try = m_ghosts[idxtry].m_y;
						uint32_t y1try = y0try + m_ghosts[idxtry].getSize();
						uint32_t z0try = m_ghosts[idxtry].m_z;
						uint32_t z1try = z0try + m_ghosts[idxtry].getSize();
						uint8_t level = oct->m_level;
						uint8_t leveltry = m_ghosts[idxtry].getLevel();

						if (Dx == Dxstar && Dy == Dystar && Dz == Dzstar){
							if (leveltry > level){
								if((abs(cx)*abs(cz)*((y0try>=y0)*(y0try<y1))) + (abs(cy)*abs(cz)*((x0try>=x0)*(x0try<x1))) + (abs(cx)*abs(cy)*((z0try>=z0)*(z0try<z1)))){
									neighbours.push_back(idxtry);
									isghost.push_back(true);
								}
							}
							if (leveltry < level){
								if((abs(cx)*abs(cz)*((y0>=y0try)*(y0<y1try))) + (abs(cy)*abs(cz)*((x0>=x0try)*(x0<x1try))) + (abs(cx)*abs(cy)*((z0>=z0try)*(z0<z1try)))){
									neighbours.push_back(idxtry);
									isghost.push_back(true);
								}
							}
						}

						idxtry++;
						if(idxtry>m_sizeGhosts-1){
							break;
						}
						Mortontry = m_ghosts[idxtry].computeMorton();
					}
				}
			}
		}
		// Search in octants
		//Build Morton number of virtual neigh of same size
		// Search morton in octants
		// If a even face morton is lower than morton of oct, if odd higher
		// ---> can i search only before or after idx in octants
		Mortontry = oct->computeMorton();
		int32_t jump = (Mortontry > Morton) ? int32_t(idx/2+1) : int32_t((noctants -idx)/2+1);
		idxtry = uint32_t(idx +((Mortontry<Morton)-(Mortontry>Morton))*jump);
		if (idxtry > noctants-1)
			idxtry = noctants-1;
		while(abs(jump) > 0){
			Mortontry = m_octants[idxtry].computeMorton();
			jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
			idxtry += jump;
			if (idxtry > m_octants.size()-1){
				if (jump > 0){
					idxtry = m_octants.size() - 1;
					jump = 0;
				}
				else if (jump < 0){
					idxtry = 0;
					jump = 0;
				}
			}
		}
		if(m_octants[idxtry].computeMorton() == Morton && m_octants[idxtry].m_level == oct->m_level){
			//Found neighbour of same size
			isghost.push_back(false);
			neighbours.push_back(idxtry);
			return;
		}
		else{
			// Step until the mortontry lower than morton (one idx of distance)
			{
				while(m_octants[idxtry].computeMorton() < Morton){
					idxtry++;
					if(idxtry > noctants-1){
						idxtry = noctants-1;
						break;
					}
				}
				while(m_octants[idxtry].computeMorton() > Morton){
					idxtry--;
					if(idxtry > noctants-1){
						idxtry = 0;
						break;
					}
				}
			}
			if (idxtry < noctants){
				if(m_octants[idxtry].computeMorton() == Morton && m_octants[idxtry].m_level == oct->m_level){
					//Found neighbour of same size
					isghost.push_back(false);
					neighbours.push_back(idxtry);
					return;
				}
				// Compute Last discendent of virtual octant of same size
				Octant last_desc = samesizeoct.buildLastDesc();
				uint64_t Mortonlast = last_desc.computeMorton();
				Mortontry = m_octants[idxtry].computeMorton();
				while(Mortontry <= Mortonlast && idxtry <= noctants-1){
					Dx = int32_t(abs(cx))*(-int32_t(oct->m_x) + int32_t(m_octants[idxtry].m_x));
					Dy = int32_t(abs(cy))*(-int32_t(oct->m_y) + int32_t(m_octants[idxtry].m_y));
					Dz = int32_t(abs(cz))*(-int32_t(oct->m_z) + int32_t(m_octants[idxtry].m_z));
					Dxstar = int32_t((cx-1)/2)*(m_octants[idxtry].getSize()) + int32_t((cx+1)/2)*size;
					Dystar = int32_t((cy-1)/2)*(m_octants[idxtry].getSize()) + int32_t((cy+1)/2)*size;
					Dzstar = int32_t((cz-1)/2)*(m_octants[idxtry].getSize()) + int32_t((cz+1)/2)*size;

					uint32_t x0 = oct->m_x;
					uint32_t x1 = x0 + size;
					uint32_t y0 = oct->m_y;
					uint32_t y1 = y0 + size;
					uint32_t z0 = oct->m_z;
					uint32_t z1 = z0 + size;
					uint32_t x0try = m_octants[idxtry].m_x;
					uint32_t x1try = x0try + m_octants[idxtry].getSize();
					uint32_t y0try = m_octants[idxtry].m_y;
					uint32_t y1try = y0try + m_octants[idxtry].getSize();
					uint32_t z0try = m_octants[idxtry].m_z;
					uint32_t z1try = z0try + m_octants[idxtry].getSize();
					uint8_t level = oct->m_level;
					uint8_t leveltry = m_octants[idxtry].getLevel();

					if (Dx == Dxstar && Dy == Dystar && Dz == Dzstar){
						if (leveltry > level){
							if((abs(cx)*abs(cz)*((y0try>=y0)*(y0try<y1))) + (abs(cy)*abs(cz)*((x0try>=x0)*(x0try<x1))) + (abs(cx)*abs(cy)*((z0try>=z0)*(z0try<z1)))){
								neighbours.push_back(idxtry);
								isghost.push_back(false);
							}
						}
						if (leveltry < level){
							if((abs(cx)*abs(cz)*((y0>=y0try)*(y0<y1try))) + (abs(cy)*abs(cz)*((x0>=x0try)*(x0<x1try))) + (abs(cx)*abs(cy)*((z0>=z0try)*(z0<z1try)))){
								neighbours.push_back(idxtry);
								isghost.push_back(false);
							}
						}
					}

					idxtry++;
					if(idxtry>noctants-1){
						break;
					}
					Mortontry = m_octants[idxtry].computeMorton();
				}
			}
		}
		return;
	}
	else{
		// Boundary Face
		return;
	}
};

// =================================================================================== //

/*! Finds neighbours of idx-th octant through iedge in vector m_octants.
 * Returns a vector (empty if iedge is a bound edge) with the index of neighbours
 * in their structure (octants or ghosts) and sets isghost[i] = true if the
 * i-th neighbour is ghost in the local tree.
 * \param[in] oct Pointer to the target local octant.
 * \param[in] iedge local index of the edge.
 * \param[out] neighbours Vector with the local indices of the neighbours (size = 0 if boundary edge).
 * \param[out] isghost Vector with the information about the identity of each neighbour (is the i-th neighours a ghost octant?).
 */
void
LocalTree::findEdgeNeighbours(Octant* oct, uint8_t iedge, u32vector & neighbours, vector<bool> & isghost){

	uint64_t  		Morton, Mortontry;
	uint32_t  		noctants = getNumOctants();
	uint32_t 		idxtry;
	uint32_t 		size = oct->getSize();
	uint8_t 		iface1, iface2;
	int32_t 		Dx, Dy, Dz;
	int32_t 		Dxstar,Dystar,Dzstar;

	//Alternative to switch case
	int8_t cx = m_global.m_edgeCoeffs[iedge][0];
	int8_t cy = m_global.m_edgeCoeffs[iedge][1];
	int8_t cz = m_global.m_edgeCoeffs[iedge][2];

	isghost.clear();
	neighbours.clear();

	// Default if iedge is nface<iedge<0
	if (iedge > m_global.m_nfaces*2){
		return;
	}

	// Check if octants edge is a process boundary
	iface1 = m_global.m_edgeFace[iedge][0];
	iface2 = m_global.m_edgeFace[iedge][1];

	// Check if octants edge is a boundary
	if (oct->m_info[iface1] == false && oct->m_info[iface2] == false){

		//Build Morton number of virtual neigh of same size
		Octant samesizeoct(m_dim, oct->m_level, oct->m_x+cx*size, oct->m_y+cy*size, oct->m_z+cz*size, m_global.m_maxLevel);
		Morton = samesizeoct.computeMorton(); //mortonEncode_magicbits(oct->m_x-size,oct->m_y,oct->m_z);

		//SEARCH IN GHOSTS
		if (m_ghosts.size()>0){
			// Search in ghosts
			uint32_t idxghost = uint32_t(m_sizeGhosts/2);
			Octant* octghost = &m_ghosts[idxghost];

			// Search morton in octants
			// If a even face morton is lower than morton of oct, if odd higher
			// ---> can i search only before or after idx in octants
			Mortontry = octghost->computeMorton();
			int32_t jump = int32_t(idxghost/2+1);
			idxtry = uint32_t(idxghost +((Mortontry<Morton)-(Mortontry>Morton))*jump);
			while(abs(jump) > 0){
				Mortontry = m_ghosts[idxtry].computeMorton();
				jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
				idxtry += jump;
				if (idxtry > m_ghosts.size()-1){
					if (jump > 0){
						idxtry = m_ghosts.size() - 1;
						jump = 0;
					}
					else if (jump < 0){
						idxtry = 0;
						jump = 0;
					}
				}
			}
			if(m_ghosts[idxtry].computeMorton() == Morton && m_ghosts[idxtry].m_level == oct->m_level){
				//Found neighbour of same size
				isghost.push_back(true);
				neighbours.push_back(idxtry);
				return;
			}
			else{
				// Step until the mortontry lower than morton (one idx of distance)
				{
					while(m_ghosts[idxtry].computeMorton() < Morton){
						idxtry++;
						if(idxtry > m_ghosts.size()-1){
							idxtry = m_ghosts.size()-1;
							break;
						}
					}
					while(m_ghosts[idxtry].computeMorton() > Morton){
						idxtry--;
						if(idxtry > m_ghosts.size()-1){
							idxtry = 0;
							break;
						}
					}
				}
				if(idxtry < m_sizeGhosts){
					if(m_ghosts[idxtry].computeMorton() == Morton && m_ghosts[idxtry].m_level == oct->m_level){
						//Found neighbour of same size
						isghost.push_back(true);
						neighbours.push_back(idxtry);
						return;
					}
					// Compute Last discendent of virtual octant of same size
					Octant last_desc = samesizeoct.buildLastDesc();
					uint64_t Mortonlast = last_desc.computeMorton();
					Mortontry = m_ghosts[idxtry].computeMorton();
					while(Mortontry <= Mortonlast && idxtry < m_ghosts.size()){
						Dx = int32_t(abs(cx))*(-int32_t(oct->m_x) + int32_t(m_ghosts[idxtry].m_x));
						Dy = int32_t(abs(cy))*(-int32_t(oct->m_y) + int32_t(m_ghosts[idxtry].m_y));
						Dz = int32_t(abs(cz))*(-int32_t(oct->m_z) + int32_t(m_ghosts[idxtry].m_z));
						Dxstar = int32_t((cx-1)/2)*(m_ghosts[idxtry].getSize()) + int32_t((cx+1)/2)*size;
						Dystar = int32_t((cy-1)/2)*(m_ghosts[idxtry].getSize()) + int32_t((cy+1)/2)*size;
						Dzstar = int32_t((cz-1)/2)*(m_ghosts[idxtry].getSize()) + int32_t((cz+1)/2)*size;

						uint32_t x0 = oct->m_x;
						uint32_t x1 = x0 + size;
						uint32_t y0 = oct->m_y;
						uint32_t y1 = y0 + size;
						uint32_t z0 = oct->m_z;
						uint32_t z1 = z0 + size;
						uint32_t x0try = m_ghosts[idxtry].m_x;
						uint32_t x1try = x0try + m_ghosts[idxtry].getSize();
						uint32_t y0try = m_ghosts[idxtry].m_y;
						uint32_t y1try = y0try + m_ghosts[idxtry].getSize();
						uint32_t z0try = m_ghosts[idxtry].m_z;
						uint32_t z1try = z0try + m_ghosts[idxtry].getSize();
						uint8_t level = oct->m_level;
						uint8_t leveltry = m_ghosts[idxtry].getLevel();

						if (Dx == Dxstar && Dy == Dystar && Dz == Dzstar){
							if (leveltry > level){
								if((abs(cx)*abs(cz)*((y0try>=y0)*(y0try<y1))) + (abs(cy)*abs(cz)*((x0try>=x0)*(x0try<x1))) + (abs(cx)*abs(cy)*((z0try>=z0)*(z0try<z1)))){
									neighbours.push_back(idxtry);
									isghost.push_back(true);
								}
							}
							if (leveltry < level){
								if((abs(cx)*abs(cz)*((y0>=y0try)*(y0<y1try))) + (abs(cy)*abs(cz)*((x0>=x0try)*(x0<x1try))) + (abs(cx)*abs(cy)*((z0>=z0try)*(z0<z1try)))){
									neighbours.push_back(idxtry);
									isghost.push_back(true);
								}
							}
						}

						idxtry++;
						if(idxtry>m_sizeGhosts-1){
							break;
						}
						Mortontry = m_ghosts[idxtry].computeMorton();
					}
				}
			}
		}

		// Search in octants
		//Build Morton number of virtual neigh of same size
		// Search morton in octants
		// If a even face morton is lower than morton of oct, if odd higher
		// ---> can i search only before or after idx in octants
		int32_t jump = int32_t((noctants)/2+1);
		idxtry = uint32_t(jump);
		if (idxtry > noctants-1)
			idxtry = noctants-1;
		while(abs(jump) > 0){
			Mortontry = m_octants[idxtry].computeMorton();
			jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
			idxtry += jump;
			if (idxtry > m_octants.size()-1){
				if (jump > 0){
					idxtry = m_octants.size() - 1;
					jump = 0;
				}
				else if (jump < 0){
					idxtry = 0;
					jump = 0;
				}
			}
		}
		if(m_octants[idxtry].computeMorton() == Morton && m_octants[idxtry].m_level == oct->m_level){
			//Found neighbour of same size
			isghost.push_back(false);
			neighbours.push_back(idxtry);
			return;
		}
		else{
			// Step until the mortontry lower than morton (one idx of distance)
			{
				while(m_octants[idxtry].computeMorton() < Morton){
					idxtry++;
					if(idxtry > noctants-1){
						idxtry = noctants-1;
						break;
					}
				}
				while(m_octants[idxtry].computeMorton() > Morton){
					idxtry--;
					if(idxtry > noctants-1){
						idxtry = 0;
						break;
					}
				}
			}
			if (idxtry < noctants){
				if(m_octants[idxtry].computeMorton() == Morton && m_octants[idxtry].m_level == oct->m_level){
					//Found neighbour of same size
					isghost.push_back(false);
					neighbours.push_back(idxtry);
					return;
				}
				// Compute Last discendent of virtual octant of same size
				Octant last_desc = samesizeoct.buildLastDesc();
				uint64_t Mortonlast = last_desc.computeMorton();
				Mortontry = m_octants[idxtry].computeMorton();
				while(Mortontry <= Mortonlast && idxtry <= noctants-1){
					Dx = int32_t(abs(cx))*(-int32_t(oct->m_x) + int32_t(m_octants[idxtry].m_x));
					Dy = int32_t(abs(cy))*(-int32_t(oct->m_y) + int32_t(m_octants[idxtry].m_y));
					Dz = int32_t(abs(cz))*(-int32_t(oct->m_z) + int32_t(m_octants[idxtry].m_z));
					Dxstar = int32_t((cx-1)/2)*(m_octants[idxtry].getSize()) + int32_t((cx+1)/2)*size;
					Dystar = int32_t((cy-1)/2)*(m_octants[idxtry].getSize()) + int32_t((cy+1)/2)*size;
					Dzstar = int32_t((cz-1)/2)*(m_octants[idxtry].getSize()) + int32_t((cz+1)/2)*size;

					uint32_t x0 = oct->m_x;
					uint32_t x1 = x0 + size;
					uint32_t y0 = oct->m_y;
					uint32_t y1 = y0 + size;
					uint32_t z0 = oct->m_z;
					uint32_t z1 = z0 + size;
					uint32_t x0try = m_octants[idxtry].m_x;
					uint32_t x1try = x0try + m_octants[idxtry].getSize();
					uint32_t y0try = m_octants[idxtry].m_y;
					uint32_t y1try = y0try + m_octants[idxtry].getSize();
					uint32_t z0try = m_octants[idxtry].m_z;
					uint32_t z1try = z0try + m_octants[idxtry].getSize();
					uint8_t level = oct->m_level;
					uint8_t leveltry = m_octants[idxtry].getLevel();

					if (Dx == Dxstar && Dy == Dystar && Dz == Dzstar){
						if (leveltry > level){
							if((abs(cx)*abs(cz)*((y0try>=y0)*(y0try<y1))) + (abs(cy)*abs(cz)*((x0try>=x0)*(x0try<x1))) + (abs(cx)*abs(cy)*((z0try>=z0)*(z0try<z1)))){
								neighbours.push_back(idxtry);
								isghost.push_back(false);
							}
						}
						if (leveltry < level){
							if((abs(cx)*abs(cz)*((y0>=y0try)*(y0<y1try))) + (abs(cy)*abs(cz)*((x0>=x0try)*(x0<x1try))) + (abs(cx)*abs(cy)*((z0>=z0try)*(z0<z1try)))){
								neighbours.push_back(idxtry);
								isghost.push_back(false);
							}
						}
					}

					idxtry++;
					if(idxtry>noctants-1){
						break;
					}
					Mortontry = m_octants[idxtry].computeMorton();
				}
			}
		}
		return;
	}
	else{
		// Boundary Face
		return;
	}

};

// =================================================================================== //

/*! Finds neighbours of idx-th ghost through iedge in vector m_octants.
 * Returns a vector (empty if iedge is not a pbound edge) with the index of neighbours
 * in the structure m_octants.
 * \param[in] idx Local index of the target ghost octant.
 * \param[in] iedge local index of the edge.
 * \param[out] neighbours Vector with the local indices of the local octants neighbours (size = 0 if boundary edge).
 */
void
LocalTree::findGhostEdgeNeighbours(uint32_t idx, uint8_t iedge, u32vector & neighbours){

	uint64_t  		Morton, Mortontry;
	uint32_t  		noctants = getNumOctants();
	uint32_t 		idxtry;
	Octant* 		oct = &m_ghosts[idx];
	uint32_t 		size = oct->getSize();
	uint8_t 		iface1, iface2;
	int32_t 		Dx, Dy, Dz;
	int32_t 		Dxstar,Dystar,Dzstar;

	//Alternative to switch case
	int8_t cx = m_global.m_edgeCoeffs[iedge][0];
	int8_t cy = m_global.m_edgeCoeffs[iedge][1];
	int8_t cz = m_global.m_edgeCoeffs[iedge][2];

	neighbours.clear();

	// Default if iedge is nface<iedge<0
	if (iedge > m_global.m_nfaces*2){
		return;
	}

	// Check if octants edge is a process boundary
	iface1 = m_global.m_edgeFace[iedge][0];
	iface2 = m_global.m_edgeFace[iedge][1];

	// Check if octants edge is a pboundary edge
	if (oct->m_info[iface1+6] == true || oct->m_info[iface2+6] == true){

		//Build Morton number of virtual neigh of same size
		Octant samesizeoct(m_dim, oct->m_level, oct->m_x+cx*size, oct->m_y+cy*size, oct->m_z+cz*size, m_global.m_maxLevel);
		Morton = samesizeoct.computeMorton();

		//Build Morton number of virtual neigh of same size
		// Search morton in octants
		// If a even face morton is lower than morton of oct, if odd higher
		// ---> can i search only before or after idx in octants
		int32_t jump = getNumOctants()/2;
		idxtry = uint32_t(getNumOctants()/2);
		while(abs(jump) > 0){
			Mortontry = m_octants[idxtry].computeMorton();
			jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
			idxtry += jump;
			if (idxtry > m_octants.size()-1){
				if (jump > 0){
					idxtry = m_octants.size() - 1;
					jump = 0;
				}
				else if (jump < 0){
					idxtry = 0;
					jump = 0;
				}
			}
		}
		if(m_octants[idxtry].computeMorton() == Morton && m_octants[idxtry].m_level == oct->m_level){
			//Found neighbour of same size
			neighbours.push_back(idxtry);
			return;
		}
		else{
			// Step until the mortontry lower than morton (one idx of distance)
			{
				while(m_octants[idxtry].computeMorton() < Morton){
					idxtry++;
					if(idxtry > noctants-1){
						idxtry = noctants-1;
						break;
					}
				}
				while(m_octants[idxtry].computeMorton() > Morton){
					idxtry--;
					if(idxtry > noctants-1){
						idxtry = 0;
						break;
					}
				}
			}
			if (idxtry < noctants){
				if(m_octants[idxtry].computeMorton() == Morton && m_octants[idxtry].m_level == oct->m_level){
					//Found neighbour of same size
					neighbours.push_back(idxtry);
					return;
				}
				// Compute Last discendent of virtual octant of same size
				Octant last_desc = samesizeoct.buildLastDesc();
				uint64_t Mortonlast = last_desc.computeMorton();
				Mortontry = m_octants[idxtry].computeMorton();
				while(Mortontry <= Mortonlast && idxtry <= noctants-1){
					Dx = int32_t(abs(cx))*(-int32_t(oct->m_x) + int32_t(m_octants[idxtry].m_x));
					Dy = int32_t(abs(cy))*(-int32_t(oct->m_y) + int32_t(m_octants[idxtry].m_y));
					Dz = int32_t(abs(cz))*(-int32_t(oct->m_z) + int32_t(m_octants[idxtry].m_z));
					Dxstar = int32_t((cx-1)/2)*(m_octants[idxtry].getSize()) + int32_t((cx+1)/2)*size;
					Dystar = int32_t((cy-1)/2)*(m_octants[idxtry].getSize()) + int32_t((cy+1)/2)*size;
					Dzstar = int32_t((cz-1)/2)*(m_octants[idxtry].getSize()) + int32_t((cz+1)/2)*size;

					uint32_t x0 = oct->m_x;
					uint32_t x1 = x0 + size;
					uint32_t y0 = oct->m_y;
					uint32_t y1 = y0 + size;
					uint32_t z0 = oct->m_z;
					uint32_t z1 = z0 + size;
					uint32_t x0try = m_octants[idxtry].m_x;
					uint32_t x1try = x0try + m_octants[idxtry].getSize();
					uint32_t y0try = m_octants[idxtry].m_y;
					uint32_t y1try = y0try + m_octants[idxtry].getSize();
					uint32_t z0try = m_octants[idxtry].m_z;
					uint32_t z1try = z0try + m_octants[idxtry].getSize();
					uint8_t level = oct->m_level;
					uint8_t leveltry = m_octants[idxtry].getLevel();

					if (Dx == Dxstar && Dy == Dystar && Dz == Dzstar){
						if (leveltry > level){
							if((abs(cx)*abs(cz)*((y0try>=y0)*(y0try<y1))) + (abs(cy)*abs(cz)*((x0try>=x0)*(x0try<x1))) + (abs(cx)*abs(cy)*((z0try>=z0)*(z0try<z1)))){
								neighbours.push_back(idxtry);
							}
						}
						if (leveltry < level){
							if((abs(cx)*abs(cz)*((y0>=y0try)*(y0<y1try))) + (abs(cy)*abs(cz)*((x0>=x0try)*(x0<x1try))) + (abs(cx)*abs(cy)*((z0>=z0try)*(z0<z1try)))){
								neighbours.push_back(idxtry);
							}
						}
					}
					idxtry++;
					if(idxtry>noctants-1){
						break;
					}
					Mortontry = m_octants[idxtry].computeMorton();
				}
			}
		}
		return;
	}
	else{
		// Boundary Face
		return;
	}
};

// =================================================================================== //

/*! Finds neighbours of idx-th octant through inode in vector m_octants.
 * Returns a vector (empty if inode is a bound node) with the index of neighbours
 * in their structure (octants or ghosts) and sets isghost[i] = true if the
 * i-th neighbour is ghost in the local tree.
 * \param[in] oct Pointer to the target local octant.
 * \param[in] inode local index of the node.
 * \param[out] neighbours Vector with the local indices of the neighbours (size = 0 if boundary node).
 * \param[out] isghost Vector with the information about the identity of each neighbour (is the i-th neighours a ghost octant?).
 */
void
LocalTree::findNodeNeighbours(Octant* oct, uint8_t inode, u32vector & neighbours, vector<bool> & isghost){

	uint64_t  	Morton, Mortontry;
	uint32_t  	noctants = getNumOctants();
	uint32_t 	idxtry;
	uint32_t 	size = oct->getSize();
	uint8_t 	iface1, iface2, iface3;

	int8_t 			cxyz[3] = {0,0,0};
	for (int idim=0; idim<m_dim; idim++){
		cxyz[idim] = m_global.m_nodeCoeffs[inode][idim];
	}

	isghost.clear();
	neighbours.clear();

	// Default if inode is nnodes<inode<0
	if (inode > m_global.m_nnodes){
		return;
	}

	// Check if octants node is a boundary
	iface1 = m_global.m_nodeFace[inode][0];
	iface2 = m_global.m_nodeFace[inode][1];
	iface3 = m_global.m_nodeFace[inode][m_dim-1];

	// Check if octants node is a boundary
	if (oct->m_info[iface1] == false && oct->m_info[iface2] == false && oct->m_info[iface3] == false){

		//Build Morton number of virtual neigh of same size
		Octant samesizeoct(m_dim, oct->m_level, oct->m_x+cxyz[0]*size, oct->m_y+cxyz[1]*size, oct->m_z+cxyz[2]*size, m_global.m_maxLevel);
		Morton = samesizeoct.computeMorton(); //mortonEncode_magicbits(oct->m_x-size,oct->m_y,oct->m_z);

		//SEARCH IN GHOSTS

		if (m_ghosts.size()>0){
			// Search in ghosts
			uint32_t idxghost = uint32_t(m_sizeGhosts/2);
			Octant* octghost = &m_ghosts[idxghost];

			// Search morton in octants
			// If a even face morton is lower than morton of oct, if odd higher
			// ---> can i search only before or after idx in octants
			Mortontry = octghost->computeMorton();
			int32_t jump = (Mortontry > Morton) ? int32_t(idxghost/2+1) : int32_t((m_sizeGhosts -idxghost)/2+1);
			idxtry = uint32_t(idxghost +((Mortontry<Morton)-(Mortontry>Morton))*jump);
			if (idxtry > m_sizeGhosts-1)
				idxtry = m_sizeGhosts-1;
			while(abs(jump) > 0){
				Mortontry = m_ghosts[idxtry].computeMorton();
				jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
				idxtry += jump;
				if (idxtry > m_ghosts.size()-1){
					if (jump > 0){
						idxtry = m_ghosts.size() - 1;
						jump = 0;
					}
					else if (jump < 0){
						idxtry = 0;
						jump = 0;
					}
				}
			}
			if(m_ghosts[idxtry].computeMorton() == Morton && m_ghosts[idxtry].m_level == oct->m_level){
				//Found neighbour of same size
				isghost.push_back(true);
				neighbours.push_back(idxtry);
				return;
			}
			else{
				// Step until the mortontry lower than morton (one idx of distance)
				{
					while(m_ghosts[idxtry].computeMorton() < Morton){
						idxtry++;
						if(idxtry > m_ghosts.size()-1){
							idxtry = m_ghosts.size()-1;
							break;
						}
					}
					while(m_ghosts[idxtry].computeMorton() > Morton){
						idxtry--;
						if(idxtry > m_ghosts.size()-1){
							idxtry = 0;
							break;
						}
					}
				}
				if(idxtry < m_sizeGhosts){
					if(m_ghosts[idxtry].computeMorton() == Morton && m_ghosts[idxtry].m_level == oct->m_level){
						//Found neighbour of same size
						isghost.push_back(true);
						neighbours.push_back(idxtry);
						return;
					}
					// Compute Last discendent of virtual octant of same size
					Octant last_desc = samesizeoct.buildLastDesc();
					uint64_t Mortonlast = last_desc.computeMorton();
					Mortontry = m_ghosts[idxtry].computeMorton();
					int32_t Dx[3] = {0,0,0};
					int32_t Dxstar[3] = {0,0,0};
					u32array3 coord = oct->getCoord();
					u32array3 coordtry = m_ghosts[idxtry].getCoord();
					u32array3 coord1 = { {1,1,1} };
					u32array3 coordtry1 = { {1,1,1} };
					while(Mortontry <= Mortonlast && idxtry < m_sizeGhosts){
						for (int idim=0; idim<m_dim; idim++){
						  //Dx[idim] 		= abs(int((abs(cxyz[idim]))*(-coord[idim] + coordtry[idim])));
							Dx[idim] 		= int32_t(int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]));
							Dxstar[idim]	= int32_t((cxyz[idim]-1)/2)*(m_ghosts[idxtry].getSize()) + int32_t((cxyz[idim]+1)/2)*size;
							coord1[idim] 	= coord[idim] + size;
							coordtry1[idim] = coordtry[idim] + m_ghosts[idxtry].getSize();
						}
						if (Dx[0] == Dxstar[0] && Dx[1] == Dxstar[1] && Dx[m_dim-1] == Dxstar[m_dim-1]){
							neighbours.push_back(idxtry);
							isghost.push_back(true);
							return;
						}
						idxtry++;
						if(idxtry>m_sizeGhosts-1){
							break;
						}
						Mortontry = m_ghosts[idxtry].computeMorton();
					}
				}
			}
		}

		// Search in octants
		//Build Morton number of virtual neigh of same size
		// Search morton in octants
		// If a even face morton is lower than morton of oct, if odd higher
		// ---> can i search only before or after idx in octants
		int32_t jump = (int32_t((noctants)/2+1));
		idxtry = uint32_t(jump);
		if (idxtry > noctants-1)
			idxtry = noctants-1;
		while(abs(jump) > 0){
			Mortontry = m_octants[idxtry].computeMorton();
			jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
			idxtry += jump;
			if (idxtry > m_octants.size()-1){
				if (jump > 0){
					idxtry = m_octants.size() - 1;
					jump = 0;
				}
				else if (jump < 0){
					idxtry = 0;
					jump = 0;
				}
			}
		}
		if(m_octants[idxtry].computeMorton() == Morton && m_octants[idxtry].m_level == oct->m_level){
			//Found neighbour of same size
			isghost.push_back(false);
			neighbours.push_back(idxtry);
			return;
		}
		else{
			// Step until the mortontry lower than morton (one idx of distance)
			{
				while(m_octants[idxtry].computeMorton() < Morton){
					idxtry++;
					if(idxtry > noctants-1){
						idxtry = noctants-1;
						break;
					}
				}
				while(m_octants[idxtry].computeMorton() > Morton){
					idxtry--;
					if(idxtry > noctants-1){
						idxtry = 0;
						break;
					}
				}
			}
			if (idxtry < noctants){
				if(m_octants[idxtry].computeMorton() == Morton && m_octants[idxtry].m_level == oct->m_level){
					//Found neighbour of same size
					isghost.push_back(false);
					neighbours.push_back(idxtry);
					return;
				}
				// Compute Last discendent of virtual octant of same size
				Octant last_desc = samesizeoct.buildLastDesc();
				uint64_t Mortonlast = last_desc.computeMorton();
				Mortontry = m_octants[idxtry].computeMorton();
				int32_t Dx[3] = {0,0,0};
				int32_t Dxstar[3] = {0,0,0};
				u32array3 coord = oct->getCoord();
				u32array3 coordtry = m_octants[idxtry].getCoord();
				u32array3 coord1 = { {1,1,1} };
				u32array3 coordtry1 = { {1,1,1} };
				while(Mortontry <= Mortonlast && idxtry <= noctants-1){
					for (int idim=0; idim<m_dim; idim++){
					  //Dx[idim] 		= abs(int((abs(cxyz[idim]))*(-coord[idim] + coordtry[idim])));
					  Dx[idim] 		= int32_t(int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]));
						Dxstar[idim]	= int32_t((cxyz[idim]-1)/2)*(m_octants[idxtry].getSize()) + int32_t((cxyz[idim]+1)/2)*size;
						coord1[idim] 	= coord[idim] + size;
						coordtry1[idim] = coordtry[idim] + m_octants[idxtry].getSize();
					}
					if (Dx[0] == Dxstar[0] && Dx[1] == Dxstar[1] && Dx[m_dim-1] == Dxstar[m_dim-1]){
						neighbours.push_back(idxtry);
						isghost.push_back(false);
						return;
					}
					idxtry++;
					if(idxtry>noctants-1){
						break;
					}
					Mortontry = m_octants[idxtry].computeMorton();
				}
			}
		}
		return;
	}
	else{
		// Boundary Node
		return;
	}
};

// =================================================================================== //

/*! Finds neighbours of idx-th octant through inode in vector m_octants.
 * Returns a vector (empty if inode is a bound node) with the index of neighbours
 * in their structure (octants or ghosts) and sets isghost[i] = true if the
 * i-th neighbour is ghost in the local tree.
 * \param[in] idx Local index of the target local octant.
 * \param[in] inode local index of the node.
 * \param[out] neighbours Vector with the local indices of the neighbours (size = 0 if boundary node).
 * \param[out] isghost Vector with the information about the identity of each neighbour (is the i-th neighours a ghost octant?).
 */
void
LocalTree::findNodeNeighbours(uint32_t idx, uint8_t inode, u32vector & neighbours, vector<bool> & isghost){

	uint64_t  		Morton, Mortontry;
	uint32_t  		noctants = getNumOctants();
	uint32_t 		idxtry;
	Octant* 		oct = &m_octants[idx];
	uint32_t 		size = oct->getSize();
	uint8_t 		iface1, iface2, iface3;

	int8_t 			cxyz[3] = {0,0,0};
	for (int idim=0; idim<m_dim; idim++){
		cxyz[idim] = m_global.m_nodeCoeffs[inode][idim];
	}

	isghost.clear();
	neighbours.clear();

	// Default if inode is nnodes<inode<0
	if (inode > m_global.m_nnodes){
		return;
	}

	// Check if octants node is a boundary
	iface1 = m_global.m_nodeFace[inode][0];
	iface2 = m_global.m_nodeFace[inode][1];
	iface3 = m_global.m_nodeFace[inode][m_dim-1];

	// Check if octants node is a boundary
	if (oct->m_info[iface1] == false && oct->m_info[iface2] == false && oct->m_info[iface3] == false){

		//Build Morton number of virtual neigh of same size
		Octant samesizeoct(m_dim, oct->m_level, oct->m_x+cxyz[0]*size, oct->m_y+cxyz[1]*size, oct->m_z+cxyz[2]*size, m_global.m_maxLevel);
		Morton = samesizeoct.computeMorton(); //mortonEncode_magicbits(oct->m_x-size,oct->m_y,oct->m_z);

		//SEARCH IN GHOSTS
		if (m_ghosts.size()>0){
			// Search in ghosts
			uint32_t idxghost = uint32_t(m_sizeGhosts/2);
			Octant* octghost = &m_ghosts[idxghost];

			// Search morton in octants
			// If a even face morton is lower than morton of oct, if odd higher
			// ---> can i search only before or after idx in octants
			Mortontry = octghost->computeMorton();
			int32_t jump = (Mortontry > Morton) ? int32_t(idxghost/2+1) : int32_t((m_sizeGhosts -idxghost)/2+1);
			idxtry = uint32_t(idxghost +((Mortontry<Morton)-(Mortontry>Morton))*jump);
			if (idxtry > m_sizeGhosts-1)
				idxtry = m_sizeGhosts-1;
			while(abs(jump) > 0){
				Mortontry = m_ghosts[idxtry].computeMorton();
				jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
				idxtry += jump;
				if (idxtry > m_ghosts.size()-1){
					if (jump > 0){
						idxtry = m_ghosts.size() - 1;
						jump = 0;
					}
					else if (jump < 0){
						idxtry = 0;
						jump = 0;
					}
				}
			}
			if(m_ghosts[idxtry].computeMorton() == Morton && m_ghosts[idxtry].m_level == oct->m_level){
				//Found neighbour of same size
				isghost.push_back(true);
				neighbours.push_back(idxtry);
				return;
			}
			else{
				// Step until the mortontry lower than morton (one idx of distance)
				{
					while(m_ghosts[idxtry].computeMorton() < Morton){
						idxtry++;
						if(idxtry > m_ghosts.size()-1){
							idxtry = m_ghosts.size()-1;
							break;
						}
					}
					while(m_ghosts[idxtry].computeMorton() > Morton){
						idxtry--;
						if(idxtry > m_ghosts.size()-1){
							idxtry = 0;
							break;
						}
					}
				}
				if(idxtry < m_sizeGhosts){
					if(m_ghosts[idxtry].computeMorton() == Morton && m_ghosts[idxtry].m_level == oct->m_level){
						//Found neighbour of same size
						isghost.push_back(true);
						neighbours.push_back(idxtry);
						return;
					}
					// Compute Last discendent of virtual octant of same size
					Octant last_desc = samesizeoct.buildLastDesc();
					uint64_t Mortonlast = last_desc.computeMorton();
					Mortontry = m_ghosts[idxtry].computeMorton();
					int32_t Dx[3] = {0,0,0};
					int32_t Dxstar[3] = {0,0,0};
					u32array3 coord = oct->getCoord();
					u32array3 coordtry = m_ghosts[idxtry].getCoord();
					u32array3 coord1 = { {1,1,1} };
					u32array3 coordtry1 = { {1,1,1} };
					while(Mortontry <= Mortonlast && idxtry < m_sizeGhosts){
						for (int idim=0; idim<m_dim; idim++){
						  //							Dx[idim] 		= abs(int((abs(cxyz[idim]))*(-coord[idim] + coordtry[idim])));
						  Dx[idim] 		= int32_t(int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]));
							Dxstar[idim]	= int32_t((cxyz[idim]-1)/2)*(m_ghosts[idxtry].getSize()) + int32_t((cxyz[idim]+1)/2)*size;
							coord1[idim] 	= coord[idim] + size;
							coordtry1[idim] = coordtry[idim] + m_ghosts[idxtry].getSize();
						}
						if (Dx[0] == Dxstar[0] && Dx[1] == Dxstar[1] && Dx[m_dim-1] == Dxstar[m_dim-1]){
							neighbours.push_back(idxtry);
							isghost.push_back(true);
							return;
						}
						idxtry++;
						if(idxtry>m_sizeGhosts-1){
							break;
						}
						Mortontry = m_ghosts[idxtry].computeMorton();
					}
				}
			}
		}

		// Search in octants
		//Build Morton number of virtual neigh of same size
		// Search morton in octants
		// If a even face morton is lower than morton of oct, if odd higher
		// ---> can i search only before or after idx in octants
		int32_t jump = (int32_t((noctants)/2+1));
		idxtry = uint32_t(jump);
		if (idxtry > noctants-1)
			idxtry = noctants-1;
		while(abs(jump) > 0){
			Mortontry = m_octants[idxtry].computeMorton();
			jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
			idxtry += jump;
			if (idxtry > m_octants.size()-1){
				if (jump > 0){
					idxtry = m_octants.size() - 1;
					jump = 0;
				}
				else if (jump < 0){
					idxtry = 0;
					jump = 0;
				}
			}
		}
		if(m_octants[idxtry].computeMorton() == Morton && m_octants[idxtry].m_level == oct->m_level){
			//Found neighbour of same size
			isghost.push_back(false);
			neighbours.push_back(idxtry);
			return;
		}
		else{
			// Step until the mortontry lower than morton (one idx of distance)
			{
				while(m_octants[idxtry].computeMorton() < Morton){
					idxtry++;
					if(idxtry > noctants-1){
						idxtry = noctants-1;
						break;
					}
				}
				while(m_octants[idxtry].computeMorton() > Morton){
					idxtry--;
					if(idxtry > noctants-1){
						idxtry = 0;
						break;
					}
				}
			}
			if (idxtry < noctants){
				if(m_octants[idxtry].computeMorton() == Morton && m_octants[idxtry].m_level == oct->m_level){
					//Found neighbour of same size
					isghost.push_back(false);
					neighbours.push_back(idxtry);
					return;
				}
				// Compute Last discendent of virtual octant of same size
				Octant last_desc = samesizeoct.buildLastDesc();
				uint64_t Mortonlast = last_desc.computeMorton();
				Mortontry = m_octants[idxtry].computeMorton();
				int32_t Dx[3] = {0,0,0};
				int32_t Dxstar[3] = {0,0,0};
				u32array3 coord = oct->getCoord();
				u32array3 coordtry = m_octants[idxtry].getCoord();
				u32array3 coord1 = { {1,1,1} };
				u32array3 coordtry1 = { {1,1,1} };
				while(Mortontry <= Mortonlast && idxtry <= noctants-1){
					for (int idim=0; idim<m_dim; idim++){
					  //Dx[idim] 		= abs(int((abs(cxyz[idim]))*(-coord[idim] + coordtry[idim])));
						Dx[idim] 		= int32_t(int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]));
						Dxstar[idim]	= int32_t((cxyz[idim]-1)/2)*(m_octants[idxtry].getSize()) + int32_t((cxyz[idim]+1)/2)*size;
						coord1[idim] 	= coord[idim] + size;
						coordtry1[idim] = coordtry[idim] + m_octants[idxtry].getSize();
					}
					if (Dx[0] == Dxstar[0] && Dx[1] == Dxstar[1] && Dx[m_dim-1] == Dxstar[m_dim-1]){
						neighbours.push_back(idxtry);
						isghost.push_back(false);
						return;
					}
					idxtry++;
					if(idxtry>noctants-1){
						break;
					}
					Mortontry = m_octants[idxtry].computeMorton();
				}
			}
		}
		return;
	}
	else{
		// Boundary Node
		return;
	}
};

// =================================================================================== //

/*! Finds neighbours of idx-th ghost through inode in vector m_octants.
 * Returns a vector (empty if inode is not a pbound node) with the index of neighbours
 * in the structure m_octants.
 * \param[in] idx Local index of the target ghost octant.
 * \param[in] inode local index of the node.
 * \param[out] neighbours Vector with the local indices of the local octants neighbours (size = 0 if boundary node).
 */
void
LocalTree::findGhostNodeNeighbours(uint32_t idx, uint8_t inode, u32vector & neighbours){

	uint64_t  		Morton, Mortontry;
	uint32_t  		noctants = getNumOctants();
	uint32_t 		idxtry;
	Octant* 		oct = &m_ghosts[idx];
	uint32_t 		size = oct->getSize();
	uint8_t 		iface1, iface2, iface3;

	int8_t 			cxyz[3] = {0,0,0};
	for (int idim=0; idim<m_dim; idim++){
		cxyz[idim] = m_global.m_nodeCoeffs[inode][idim];
	}

	neighbours.clear();

	// Default if inode is nnodes<inode<0
	if (inode > m_global.m_nnodes){
		return;
	}

	// Check if octants node is a boundary
	iface1 = m_global.m_nodeFace[inode][0];
	iface2 = m_global.m_nodeFace[inode][1];
	iface3 = m_global.m_nodeFace[inode][m_dim-1];

	// Check if octants node is a pboundary node
	if (oct->m_info[iface1+6] == true || oct->m_info[iface2+6] == true || oct->m_info[iface3+6] == true){

		//Build Morton number of virtual neigh of same size
		Octant samesizeoct(m_dim, oct->m_level, oct->m_x+cxyz[0]*size, oct->m_y+cxyz[1]*size, oct->m_z+cxyz[2]*size, m_global.m_maxLevel);
		Morton = samesizeoct.computeMorton(); //mortonEncode_magicbits(oct->m_x-size,oct->m_y,oct->m_z);
		int32_t jump = noctants/2;
		idxtry = jump;
		while(abs(jump) > 0){
			Mortontry = m_octants[idxtry].computeMorton();
			jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
			idxtry += jump;
		}
		if(m_octants[idxtry].computeMorton() == Morton && m_octants[idxtry].m_level == oct->m_level){
			//Found neighbour of same size
			neighbours.push_back(idxtry);
			return;
		}
		else{
			// Step until the mortontry lower than morton (one idx of distance)
			{
				while(m_octants[idxtry].computeMorton() < Morton){
					idxtry++;
					if(idxtry > noctants-1){
						idxtry = noctants-1;
						break;
					}
				}
				while(m_octants[idxtry].computeMorton() > Morton){
					idxtry--;
					if(idxtry > noctants-1){
						idxtry = 0;
						break;
					}
				}
			}
			if (idxtry < noctants){
				if(m_octants[idxtry].computeMorton() == Morton && m_octants[idxtry].m_level == oct->m_level){
					//Found neighbour of same size
					neighbours.push_back(idxtry);
					return;
				}
				// Compute Last discendent of virtual octant of same size
				Octant last_desc = samesizeoct.buildLastDesc();
				uint64_t Mortonlast = last_desc.computeMorton();
				Mortontry = m_octants[idxtry].computeMorton();
				int32_t Dx[3] = {0,0,0};
				int32_t Dxstar[3] = {0,0,0};
				u32array3 coord = oct->getCoord();
				u32array3 coordtry = m_octants[idxtry].getCoord();
				u32array3 coord1 = { {1,1,1} };
				u32array3 coordtry1 = { {1,1,1} };
				while(Mortontry <= Mortonlast && idxtry <= noctants-1){
					for (int idim=0; idim<m_dim; idim++){
					  //Dx[idim] 		= abs(int((abs(cxyz[idim]))*(-coord[idim] + coordtry[idim])));
						Dx[idim] 		= int32_t(int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]));
						Dxstar[idim]	= int32_t((cxyz[idim]-1)/2)*(m_octants[idxtry].getSize()) + int32_t((cxyz[idim]+1)/2)*size;
						coord1[idim] 	= coord[idim] + size;
						coordtry1[idim] = coordtry[idim] + m_octants[idxtry].getSize();
					}
					if (Dx[0] == Dxstar[0] && Dx[1] == Dxstar[1] && Dx[m_dim-1] == Dxstar[m_dim-1]){
						neighbours.push_back(idxtry);
					}
					idxtry++;
					Mortontry = m_octants[idxtry].computeMorton();
				}
			}
		}
		return;
	}
	else{
		// Boundary Node
		return;
	}
};

// =================================================================================== //

/*! Compute and store in m_intersections the intersections of the local tree.
 */
void
LocalTree::computeIntersections() {

		octvector::iterator 	it, obegin, oend;
		Intersection 			intersection;
		u32vector 				neighbours;
		vector<bool>			isghost;
		uint32_t 				counter, idx;
		uint32_t 				i, nsize;
		uint8_t 				iface, iface2;

		m_intersections.clear();
		m_intersections.reserve(2*3*m_octants.size());

		counter = idx = 0;

		// Loop on ghosts
		obegin = m_ghosts.begin();
		oend = m_ghosts.end();
		for (it = obegin; it != oend; it++){
			for (iface = 0; iface < m_dim; iface++){
				iface2 = iface*2;
				findGhostNeighbours(idx, iface2, neighbours);
				nsize = neighbours.size();
				if (!(it->m_info[iface2])){
					//Internal intersection
					for (i = 0; i < nsize; i++){
						intersection.m_dim = m_dim;
						intersection.m_finer = getGhostLevel(idx) >= getLevel((int)neighbours[i]);
						intersection.m_out = intersection.m_finer;
						intersection.m_outisghost = intersection.m_finer;
						intersection.m_owners[0]  = neighbours[i];
						intersection.m_owners[1] = idx;
						intersection.m_iface = m_global.m_oppFace[iface2] - (getGhostLevel(idx) >= getLevel((int)neighbours[i]));
						intersection.m_isnew = false;
						intersection.m_isghost = true;
						intersection.m_bound = false;
						intersection.m_pbound = true;
						m_intersections.push_back(intersection);
						counter++;
					}
				}
				else{
					//Periodic intersection
					for (i = 0; i < nsize; i++){
						intersection.m_dim = m_dim;
						intersection.m_finer = getGhostLevel(idx) >= getLevel((int)neighbours[i]);
						intersection.m_out = intersection.m_finer;
						intersection.m_outisghost = intersection.m_finer;
						intersection.m_owners[0]  = neighbours[i];
						intersection.m_owners[1] = idx;
						intersection.m_iface = m_global.m_oppFace[iface2] - (getGhostLevel(idx) >= getLevel((int)neighbours[i]));
						intersection.m_isnew = false;
						intersection.m_isghost = true;
						intersection.m_bound = true;
						intersection.m_pbound = true;
						m_intersections.push_back(intersection);
						counter++;
					}
				}
			}
			idx++;
		}

		// Loop on octants
		idx=0;
		obegin = m_octants.begin();
		oend = m_octants.end();
		for (it = obegin; it != oend; it++){
			for (iface = 0; iface < m_dim; iface++){
				iface2 = iface*2;
				findNeighbours(idx, iface2, neighbours, isghost);
				nsize = neighbours.size();
				if (nsize) {
					if (!(it->m_info[iface2])){
						//Internal intersection
						for (i = 0; i < nsize; i++){
							if (isghost[i]){
								intersection.m_dim = m_dim;
								intersection.m_owners[0] = idx;
								intersection.m_owners[1] = neighbours[i];
								intersection.m_finer = (nsize>1);
								intersection.m_out = (nsize>1);
								intersection.m_outisghost = (nsize>1);
								intersection.m_iface = iface2 + (nsize>1);
								intersection.m_isnew = false;
								intersection.m_isghost = true;
								intersection.m_bound = false;
								intersection.m_pbound = true;
								m_intersections.push_back(intersection);
								counter++;
							}
							else{
								intersection.m_dim = m_dim;
								intersection.m_owners[0] = idx;
								intersection.m_owners[1] = neighbours[i];
								intersection.m_finer = (nsize>1);
								intersection.m_out = (nsize>1);
								intersection.m_outisghost = false;
								intersection.m_iface = iface2 + (nsize>1);
								intersection.m_isnew = false;
								intersection.m_isghost = false;
								intersection.m_bound = false;
								intersection.m_pbound = false;
								m_intersections.push_back(intersection);
								counter++;
							}
						}
					}
					else{
						//Periodic intersection
						for (i = 0; i < nsize; i++){
							if (isghost[i]){
								intersection.m_dim = m_dim;
								intersection.m_owners[0] = idx;
								intersection.m_owners[1] = neighbours[i];
								intersection.m_finer = (nsize>1);
								intersection.m_out = intersection.m_finer;
								intersection.m_outisghost = intersection.m_finer;
								intersection.m_iface = iface2 + (nsize>1);
								intersection.m_isnew = false;
								intersection.m_isghost = true;
								intersection.m_bound = true;
								intersection.m_pbound = true;
								m_intersections.push_back(intersection);
								counter++;
							}
							else{
								intersection.m_dim = m_dim;
								intersection.m_owners[0] = idx;
								intersection.m_owners[1] = neighbours[i];
								intersection.m_finer = (nsize>1);
								intersection.m_out = intersection.m_finer;
								intersection.m_outisghost = false;
								intersection.m_iface = iface2 + (nsize>1);
								intersection.m_isnew = false;
								intersection.m_isghost = false;
								intersection.m_bound = true;
								intersection.m_pbound = false;
								m_intersections.push_back(intersection);
								counter++;
							}
						}
					}
				}
				else{
					//Boundary intersection
					intersection.m_dim = m_dim;
					intersection.m_owners[0] = idx;
					intersection.m_owners[1] = idx;
					intersection.m_finer = 0;
					intersection.m_out = 0;
					intersection.m_outisghost = false;
					intersection.m_iface = iface2;
					intersection.m_isnew = false;
					intersection.m_isghost = false;
					intersection.m_bound = true;
					intersection.m_pbound = false;
					m_intersections.push_back(intersection);
					counter++;
				}
				if (it->m_info[iface2+1]){
					if (!(m_periodic[iface2+1])){
						//Boundary intersection
						intersection.m_dim = m_dim;
						intersection.m_owners[0] = idx;
						intersection.m_owners[1] = idx;
						intersection.m_finer = 0;
						intersection.m_out = 0;
						intersection.m_outisghost = false;
						intersection.m_iface = iface2+1;
						intersection.m_isnew = false;
						intersection.m_isghost = false;
						intersection.m_bound = true;
						intersection.m_pbound = false;
						m_intersections.push_back(intersection);
						counter++;
					}
					else{
						//Periodic intersection
						findNeighbours(idx, iface2+1, neighbours, isghost);
						nsize = neighbours.size();
						for (i = 0; i < nsize; i++){
							if (isghost[i]){
								intersection.m_dim = m_dim;
								intersection.m_owners[0] = idx;
								intersection.m_owners[1] = neighbours[i];
								intersection.m_finer = (nsize>1);
								intersection.m_out = intersection.m_finer;
								intersection.m_outisghost = intersection.m_finer;
								intersection.m_iface = iface2 + (nsize>1);
								intersection.m_isnew = false;
								intersection.m_isghost = true;
								intersection.m_bound = true;
								intersection.m_pbound = true;
								m_intersections.push_back(intersection);
								counter++;
							}
							else{
								intersection.m_dim = m_dim;
								intersection.m_owners[0] = idx;
								intersection.m_owners[1] = neighbours[i];
								intersection.m_finer = (nsize>1);
								intersection.m_out = intersection.m_finer;
								intersection.m_outisghost = false;
								intersection.m_iface = iface2 + (nsize>1);
								intersection.m_isnew = false;
								intersection.m_isghost = false;
								intersection.m_bound = true;
								intersection.m_pbound = false;
								m_intersections.push_back(intersection);
								counter++;
							}
						}
					}
				}
			}
			idx++;
		}
		intervector(m_intersections).swap(m_intersections);
	}

// =================================================================================== //
/*! Find an input Morton in octants and return the local idx
 * \param[in] Morton Morton index to be found.
 * \return Local index of the target octant (=nocts if target Morton not found).
*/
uint32_t
LocalTree::findMorton(uint64_t Morton){

	uint32_t 		nocts = m_octants.size();
	uint32_t 		idx = nocts/2;
	uint64_t 		Mortontry = m_octants[idx].computeMorton();
	int32_t 		jump = nocts/2;

	while(abs(jump)>0){
		if (Mortontry == Morton){
			return idx;
		}
		jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
		idx += jump;
		if (idx > nocts){
			return nocts-1;
		}
		Mortontry = m_octants[idx].computeMorton();
	}
	if (Mortontry<Morton){
		for (uint32_t idx2=idx; idx2<nocts; idx2++){
			Mortontry = m_octants[idx2].computeMorton();
			if (Mortontry == Morton){
				return idx2;
			}
		}
	}
	else{
		for(uint32_t idx2=0; idx2<idx+1; idx2++){
			Mortontry = m_octants[idx2].computeMorton();
			if (Mortontry == Morton){
				return idx2;
			}
		}
	}
	return nocts;
};

// =================================================================================== //
/*! Find an input Morton in ghosts and return the local idx
 * \param[in] Morton Morton index to be found.
 * \return Index of the target ghost octant (=nghosts if target Morton not found).
*/
uint32_t
LocalTree::findGhostMorton(uint64_t Morton){
	uint32_t 		nocts = m_ghosts.size();
	uint32_t 		idx = nocts/2;
	uint64_t 		Mortontry = m_ghosts[idx].computeMorton();
	int32_t 		jump = nocts/2;

	while(abs(jump)>0){
		if (Mortontry == Morton){
			return idx;
		}
		jump = ((Mortontry<Morton)-(Mortontry>Morton))*abs(jump)/2;
		idx += jump;
		if (idx > nocts){
			return nocts;
		}
		Mortontry = m_ghosts[idx].computeMorton();
	}
	if (Mortontry<Morton){
		for (uint32_t idx2=idx; idx2<nocts; idx2++){
			Mortontry = m_ghosts[idx2].computeMorton();
			if (Mortontry == Morton){
				return idx2;
			}
		}
	}
	else{
		for(uint32_t idx2=0; idx2<idx; idx2++){
			Mortontry = m_ghosts[idx2].computeMorton();
			if (Mortontry == Morton){
				return idx2;
			}
		}
	}
	return nocts;
};

// =================================================================================== //

/** Compute the connectivity of octants and store the coordinates of nodes.
 */
void
LocalTree::computeConnectivity(){
	u32arr3vector                                octnodes;
	vector<uint64_t>                             mortonList;
	unordered_map<uint64_t, array<uint32_t, 3> > nodeCoords;
	unordered_map<uint64_t, vector<uint64_t> >   nodeOctants;
	uint32_t                                     noctants = getNumOctants();
	uint32_t                                     nghosts  = m_sizeGhosts;



	// Gather node information
	octnodes.reserve(m_global.m_nnodes);

	mortonList.reserve(noctants);
	nodeCoords.reserve(noctants);
	nodeOctants.reserve(noctants);

	for (uint64_t n = 0; n < (noctants + nghosts); n++){
		if (n < noctants) {
			uint32_t octantId = n;
			m_octants[octantId].getNodes(octnodes);
		} else {
			uint32_t octantId = n - noctants;
			m_ghosts[octantId].getNodes(octnodes);
		}

		for (auto &node : octnodes){
			uint64_t morton = keyXYZ(node[0], node[1], node[2], m_global.m_maxLevel);
			if (nodeCoords.count(morton) == 0) {
				mortonList.push_back(morton);
				nodeCoords.insert({{morton, std::move(node)}});
				nodeOctants[morton].reserve(8);
			}

			nodeOctants[morton].push_back(n);
		}
	}
	std::sort(mortonList.begin(), mortonList.end());

	// Build node list and connectivity
	m_nodes.reserve(mortonList.size());
	m_connectivity.resize(noctants);
	m_ghostsConnectivity.resize(nghosts);

	uint32_t nodeId = 0;
	for (auto &morton : mortonList) {
		m_nodes.emplace_back(std::move(nodeCoords.at(morton)));
		for (const auto &n : nodeOctants.at(morton)) {
			std::vector<uint32_t> *octantConnect;
			if (n < noctants) {
				uint32_t octantId = n;
				octantConnect = &(m_connectivity[octantId]);
			} else {
				uint32_t octantId = n - noctants;
				octantConnect = &(m_ghostsConnectivity[octantId]);
			}

			if (octantConnect->size() == 0) {
				octantConnect->reserve(m_global.m_nnodes);
			}
			octantConnect->push_back(nodeId);
		}
		nodeId++;
	}
};

/*! Clear nodes vector and connectivity of octants of local tree
*/
void
LocalTree::clearConnectivity(){
	u32arr3vector().swap(m_nodes);
	u32vector2D().swap(m_connectivity);
	u32vector2D().swap(m_ghostsConnectivity);
};

/*! Updates nodes vector and connectivity of octants of local tree
*/
void
LocalTree::updateConnectivity(){
	clearConnectivity();
	computeConnectivity();
};

// =================================================================================== //

}
