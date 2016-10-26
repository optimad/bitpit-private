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

# ifndef __BITPIT_LEVELSET_MASK_HPP__
# define __BITPIT_LEVELSET_MASK_HPP__

// Standard Template Library
# include <vector>
# include <set>

namespace bitpit{

class PatchKernel;
class SurfUnstructured;
class LevelSetSegmentation;

class LevelSetMask : public LevelSetSegmentation {

    private:
    SurfUnstructured extractCellEnvelope(const std::unordered_set<long> &, const VolumeKernel &, std::unordered_map<long,long> *envelopeToMesh=nullptr);
    SurfUnstructured extractFaceEnvelope(const std::vector<long> &, const VolumeKernel &, std::unordered_map<long,long> *envelopeToMesh=nullptr);

    void orientByMesh(const std::unordered_map<long,long> &, const VolumeKernel &, SurfUnstructured &, std::unordered_set<long> &);
    void orientByMask(const std::unordered_map<long,long> &, const VolumeKernel &, const std::unordered_set<long> &, std::unordered_set<long> &);
    void orientByMark(const std::unordered_map<long,long> &, const std::unordered_set<long> &,std::unordered_set<long> &);

    void flip(SurfUnstructured &, const std::unordered_set<long> &);
    void updateFlipList( std::unordered_set<long> &, const long &, const bool &);

    public:
    ~LevelSetMask();
    LevelSetMask(int, const std::unordered_set<long> &, const VolumeKernel &);
    LevelSetMask(int, const std::vector<long> &, const VolumeKernel &);
    LevelSetMask(int, const std::vector<long> &, const std::unordered_set<long> &, const VolumeKernel &);
};

}

#endif
