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
#include "line_kernel.hpp"

namespace bitpit {

class SurfaceKernel : public PatchKernel {

public:
    // Static constant
    static const unsigned short SELECT_TRIANGLE;
    static const unsigned short SELECT_QUAD;
    static const unsigned short SELECT_ALL;

    // Types definitions
    typedef double (SurfaceKernel::*eval_f_)(long, int &) const;

    void setSpaceDimension(int dimension);
    int getSpaceDimension(void) const;

    int getVolumeCodimension() const override;
    int getSurfaceCodimension() const override;
    int getLineCodimension() const override;
    int getPointCodimension() const override;

    void extractEnvelope(LineKernel &envelope) const;

    virtual ~SurfaceKernel();
    virtual double evalCellArea(long) const;
    virtual double evalEdgeLength(long, int) const;
    virtual double evalMinEdgeLength(long, int &) const;
    virtual double evalMaxEdgeLength(long, int &) const;
    virtual double evalAngleAtVertex(long, int) const;
    virtual double evalMinAngleAtVertex(long, int &) const;
    virtual double evalMaxAngleAtVertex(long, int &) const;
    virtual double evalAspectRatio(long, int &) const;
    virtual std::array<double, 3> evalFacetNormal(long, const std::array<double, 3> &orientation = {{0., 0., 1.}}) const;
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
    static const std::map<ElementType, unsigned short> m_selectionTypes;

    void initialize();

    bool compareSelectedTypes(unsigned short, ElementType) const;
    void displayHistogram(long, const std::vector<double>&, const std::vector<double>&, const std::string&, std::ostream&, unsigned int padding = 0) const;

    bool sameOrientationAtInterface(long id);

protected:
    SurfaceKernel(bool expert);
    SurfaceKernel(int dimension, bool expert);
    SurfaceKernel(int id, int dimension, bool expert);

};

}

#endif
