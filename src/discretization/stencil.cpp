/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
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

#include <cassert>

#include "bitpit_operators.hpp"

#include "stencil.hpp"

/*!
    Output stream operator from class StencilScalar to communication buffer.

    \param[in] buffer is the output memory stream
    \param[in] stencil is the stencil to be streamed
    \result Returns the same output stream received in input.
*/
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buffer, const bitpit::StencilScalar &stencil)
{
    return operator<<(buffer, static_cast<const bitpit::BaseStencil<double> &>(stencil));
}

/*!
    Input stream operator from class StencilScalar to communication buffer.

    \param[in] buffer is the input memory stream
    \param[in] stencil is the stencil to be streamed
    \result Returns the same input stream received in input.
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, bitpit::StencilScalar &stencil)
{
    return operator>>(buffer, static_cast<bitpit::BaseStencil<double> &>(stencil));
}

/*!
    Output stream operator from class StencilVector to communication buffer.

    \param[in] buffer is the output memory stream
    \param[in] stencil is the stencil to be streamed
    \result Returns the same output stream received in input.
*/
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buffer, const bitpit::StencilVector &stencil)
{
    return operator<<(buffer, static_cast<const bitpit::BaseStencil<std::array<double, 3>> &>(stencil));
}

/*!
    Input stream operator from class StencilVector to communication buffer.

    \param[in] buffer is the input memory stream
    \param[in] stencil is the stencil to be streamed
    \result Returns the same input stream received in input.
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, bitpit::StencilVector &stencil)
{
    return operator>>(buffer, static_cast<bitpit::BaseStencil<std::array<double, 3>> &>(stencil));
}

namespace bitpit {

StencilScalar::StencilScalar()
    : BaseStencil(0.)
{
}

StencilScalar::StencilScalar(const BaseStencil &other)
    : BaseStencil(other)
{
}

StencilScalar::StencilScalar(BaseStencil &&other)
    : BaseStencil(other)
{
}

StencilVector::StencilVector()
    : BaseStencil(std::array<double, 3>{{0., 0., 0.}})
{
}

StencilVector::StencilVector(const BaseStencil &other)
    : BaseStencil(other)
{
}

StencilVector::StencilVector(BaseStencil &&other)
    : BaseStencil(other)
{
}

StencilScalar dotProduct(const StencilVector &stencil_A, const std::array<double, 3> &vector)
{
    StencilScalar stencil_B;

    const FlatVector2D<std::array<double, 3>> &weights_A = stencil_A.getWeights();

    const FlatVector2D<long> &pattern = stencil_A.getPattern();
    for (int i = 0; i < pattern.size(); ++i) {
        int nBucketItems = pattern.getItemCount(i);
        for (int j = 0; j < nBucketItems; ++j) {
            const long id = pattern.getItem(i, j);
            const std::array<double, 3> &weight_A = weights_A.getItem(i, j);
            stencil_B.addWeight(id, ::dotProduct(weight_A, vector));
        }
    }

    stencil_B.setConstant(::dotProduct(stencil_A.getConstant(), vector));

    return stencil_B;
}

StencilVector operator*(const StencilScalar &stencil_A, const std::array<double, 3> &vector)
{
    StencilVector stencil_B;

    const FlatVector2D<double> &weights_A = stencil_A.getWeights();

    const FlatVector2D<long> &pattern = stencil_A.getPattern();
    for (int i = 0; i < pattern.size(); ++i) {
        int nBucketItems = pattern.getItemCount(i);
        for (int j = 0; j < nBucketItems; ++j) {
            const long id = pattern.getItem(i, j);
            const double &weight_A = weights_A.getItem(i, j);
            stencil_B.addWeight(id, ::operator*(weight_A, vector));
        }
    }

    stencil_B.setConstant(::operator*(stencil_A.getConstant(), vector));

    return stencil_B;
}

}
