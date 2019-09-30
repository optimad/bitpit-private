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

#ifndef __BITPIT_SURFACE_KERNEL_HPP__
#define __BITPIT_SURFACE_KERNEL_HPP__

#include "patch_kernel.hpp"

namespace bitpit {

class SurfaceKernel : public PatchKernel {

public:
    // Static constant
    static const unsigned short SELECT_TRIANGLE;
    static const unsigned short SELECT_QUAD;
    static const unsigned short SELECT_ALL;
    static const std::map<ElementType, unsigned short>     m_selectionTypes;

    // Types definitions
    typedef double (SurfaceKernel::*eval_f_)(long, int &) const;

    void setSpaceDimension(int dimension);
    int getSpaceDimension(void) const;

    virtual ~SurfaceKernel();
    virtual double evalCellArea(long) const;
    virtual double evalEdgeLength(long, int) const;
    virtual double evalMinEdgeLength(long, int &) const;
    virtual double evalMaxEdgeLength(long, int &) const;
    virtual double evalAngleAtVertex(long, int) const;
    virtual double evalMinAngleAtVertex(long, int &) const;
    virtual double evalMaxAngleAtVertex(long, int &) const;
    virtual double evalAspectRatio(long, int &) const;
    virtual std::array<double, 3> evalFacetNormal(long) const;
    std::array<double, 3> evalEdgeNormal(long, int) const;
    std::array<double, 3> evalVertexNormal(long, int) const;
    virtual std::array<double, 3> evalLimitedVertexNormal(long, int, double ) const;
    double evalCellSize(long id) const override;

    bool adjustCellOrientation();
    bool adjustCellOrientation(long id, bool invert = false);
    void flipCellOrientation(long id);

    void displayQualityStats(std::ostream&, unsigned int padding = 0) const;
    std::vector<double> computeHistogram(eval_f_ funct_, std::vector<double> &bins, long &count, int n_intervals = 8, unsigned short mask = SELECT_ALL) const;

private:
    void initialize();

    bool compareSelectedTypes(unsigned short, ElementType) const;
    void displayHistogram(long, const std::vector<double>&, const std::vector<double>&, const std::string&, std::ostream&, unsigned int padding = 0) const;

    bool sameOrientationAtInterface(long id);

protected:
    int m_spaceDim;

    SurfaceKernel(bool expert);
    SurfaceKernel(int patch_dim, int space_dim, bool expert);
    SurfaceKernel(int id, int patch_dim, int space_dim, bool expert);

};

}

#endif
