/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2019 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
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

#ifndef __BITPIT_POINTCLOUD_HPP__
#define __BITPIT_POINTCLOUD_HPP__

#include <array>
#include <vector>

#include "bitpit_IO.hpp"
#include "bitpit_patchkernel.hpp"

namespace bitpit {

class PointCloud : public PointKernel {

public:
    // Constructors
    PointCloud();
    PointCloud(int dimension);
    PointCloud(int id, int dimension);
    PointCloud(std::istream &stream);

    // Clone
    std::unique_ptr<PatchKernel> clone() const override;

    // Setters
    void setExpert(bool expert);

    // Search algorithms
    long locatePoint(const std::array<double, 3> &point) const override;

protected:
    PointCloud(const PointCloud &other) = default;

    int _getDumpVersion() const override;
    void _dump(std::ostream &stream) const override;
    void _restore(std::istream &stream) override;

};

}

#endif
