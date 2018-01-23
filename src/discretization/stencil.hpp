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

#ifndef __BTPIT_STENCIL_HPP__
#define __BTPIT_STENCIL_HPP__

#include <array>
#include <unordered_map>

#include "bitpit_common.hpp"
#include "bitpit_containers.hpp"

namespace bitpit {

template<typename weight_t>
class BaseStencil;

}

template<typename weight_t>
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buffer, const bitpit::BaseStencil<weight_t> &stencil);

template<typename weight_t>
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, bitpit::BaseStencil<weight_t> &stencil);

namespace bitpit {

template <typename weight_t>
class BaseStencil {

template<typename U>
friend bitpit::OBinaryStream& (::operator<<) (bitpit::OBinaryStream &buffer, const BaseStencil<U> &stencil);
template<typename U>
friend bitpit::IBinaryStream& (::operator>>) (bitpit::IBinaryStream &buffer, BaseStencil<U> &stencil);

public:
    BaseStencil(const weight_t &zero);
    BaseStencil(int nBuckets, const weight_t &zero);

    std::size_t size() const;
    std::size_t size(int bucket) const;

    int getBucketCount() const;

    const weight_t & getWeight(long id) const;
    const weight_t & getWeight(int bucket, long id) const;
    void addWeight(long id, const weight_t &weight);
    void addWeight(int bucket, long id, const weight_t &weight);
    void setWeight(long id, const weight_t &weight);
    void setWeight(int bucket, long id, const weight_t &weight);

    const weight_t & getConstant() const;
    void setConstant(const weight_t &value);
    void sumConstant(const weight_t &value);

    const FlatVector2D<long> & getPattern() const;
    const FlatVector2D<weight_t> & getWeights() const;

    void clear();
    void flatten();
    void optimize(double tolerance = 1.e-12);
    void renumber(const std::unordered_map<long, long> &map);
    void addComplementToZero(const long id);

    void display(std::ostream &out, double factor = 1.) const;

    size_t getBinarySize() const;

    BaseStencil<weight_t> & operator*=(double factor);
    BaseStencil<weight_t> & operator/=(double factor);
    BaseStencil<weight_t> & operator+=(const BaseStencil<weight_t> &other);
    BaseStencil<weight_t> & operator-=(const BaseStencil<weight_t> &other);

private:
    weight_t m_zero;
    FlatVector2D<long> m_pattern;
    FlatVector2D<weight_t> m_weights;
    weight_t m_constant;

    weight_t * find(int bucket, long id) const;

    template<typename U = weight_t, typename std::enable_if<std::is_fundamental<U>::value>::type* = nullptr>
    bool isWeightNeglibile(int bucket, size_t k, double tolerance = 1.e-12);

    template<typename U = weight_t, typename std::enable_if<!std::is_fundamental<U>::value>::type* = nullptr>
    bool isWeightNeglibile(int bucket, size_t k, double tolerance = 1.e-12);

    void clear(int nBuckets);

};

}

// Template implementation
#include "stencil.tpp"

// Spcializations
namespace bitpit {

class StencilScalar;

bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buffer, const StencilScalar &stencil);
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, StencilScalar &stencil);

class StencilScalar : public BaseStencil<double> {

friend bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buffer, const StencilScalar &stencil);
friend bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, StencilScalar &stencil);

public:
    StencilScalar();
    StencilScalar(const BaseStencil &);
    StencilScalar(BaseStencil &&);

};

class StencilVector;

bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buffer, const StencilVector &stencil);
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, StencilVector &stencil);

class StencilVector : public BaseStencil<std::array<double, 3>> {

friend bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buffer, const StencilVector &stencil);
friend bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, StencilVector &stencil);

public:
    StencilVector();
    StencilVector(const BaseStencil &other);
    StencilVector(BaseStencil &&other);

};

}

// Operators
template <typename weight_t>
bitpit::BaseStencil<weight_t> operator*(const bitpit::BaseStencil<weight_t> &stencil, double factor);

template <typename weight_t>
bitpit::BaseStencil<weight_t> operator*(double factor, const bitpit::BaseStencil<weight_t> &stencil);

template <typename weight_t>
bitpit::BaseStencil<weight_t> operator/(const bitpit::BaseStencil<weight_t> &stencil, double factor);

template <typename weight_t>
bitpit::BaseStencil<weight_t> operator+(const bitpit::BaseStencil<weight_t> &, const bitpit::BaseStencil<weight_t> &);

template <typename weight_t>
bitpit::BaseStencil<weight_t> operator-(const bitpit::BaseStencil<weight_t> &, const bitpit::BaseStencil<weight_t> &);

bitpit::StencilScalar dotProduct(const bitpit::StencilVector &, const std::array<double,3> &);

bitpit::StencilVector operator*(const bitpit::StencilScalar &, const std::array<double,3> &);


#endif
