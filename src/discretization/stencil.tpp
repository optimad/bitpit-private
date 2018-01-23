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

#ifndef __BITPIT_STENCIL_TPP__
#define __BITPIT_STENCIL_TPP__

/*!
    Output stream operator from class BaseStencil to communication buffer.

    \param[in] buffer is the output memory stream
    \param[in] stencil is the stencil to be streamed
    \result Returns the same output stream received in input.
*/
template<typename weight_t>
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buffer, const bitpit::BaseStencil<weight_t> &stencil)
{
    buffer << stencil.m_zero;
    buffer << stencil.m_pattern;
    buffer << stencil.m_weights;
    buffer << stencil.m_constant;

    return buffer;
}

/*!
    Input stream operator from class BaseStencil to communication buffer.

    \param[in] buffer is the input memory stream
    \param[in] stencil is the stencil to be streamed
    \result Returns the same input stream received in input.
*/
template<typename weight_t>
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, bitpit::BaseStencil<weight_t> &stencil)
{
    buffer >> stencil.m_zero;
    buffer >> stencil.m_pattern;
    buffer >> stencil.m_weights;
    buffer >> stencil.m_constant;

    return buffer;
}

namespace bitpit {

/**
* Constructor
*
* \param zero is the value to be used as zero
*/
template<typename weight_t>
BaseStencil<weight_t>::BaseStencil(const weight_t &zero)
    : m_zero(zero)
{
    clear();
}

/**
* Constructor
*
* \param nBuckets are the number of buckets the stencil is divided in
* \param zero is the value to be used as zero
*/
template<typename weight_t>
BaseStencil<weight_t>::BaseStencil(int nBuckets, const weight_t &zero)
    : m_zero(zero)
{
    clear(nBuckets);
}

/**
* Size of the stencil, expressed in number of weight
*
* \param nBuckets are the number of buckets the stencil is divided in
* \param zero is the value to be used as zero
*/
template<typename weight_t>
std::size_t BaseStencil<weight_t>::size() const
{
    return m_pattern.getItemCount();
}

template<typename weight_t>
std::size_t BaseStencil<weight_t>::size(int bucket) const
{
    return m_pattern.getItemCount(bucket);
}

template<typename weight_t>
int BaseStencil<weight_t>::getBucketCount() const
{
    return m_pattern.size();
}

template<typename weight_t>
const FlatVector2D<long> & BaseStencil<weight_t>::getPattern() const
{
    return m_pattern;
}

template<typename weight_t>
const FlatVector2D<weight_t> & BaseStencil<weight_t>::getWeights() const
{
    return m_weights;
}

template<typename weight_t>
void BaseStencil<weight_t>::addWeight(long id, const weight_t &weight)
{
    addWeight(0, id, weight);
}

template<typename weight_t>
void BaseStencil<weight_t>::addWeight(int bucket, long id, const weight_t &weight)
{
    m_pattern.pushBackItem(bucket, id);
    m_weights.pushBackItem(bucket, weight);
}

template<typename weight_t>
const weight_t & BaseStencil<weight_t>::getConstant() const
{
    return m_constant;
}

template<typename weight_t>
void BaseStencil<weight_t>::setConstant(const weight_t &value)
{
    m_constant = value;
}

template<typename weight_t>
void BaseStencil<weight_t>::sumConstant(const weight_t &value)
{
    m_constant += value;
}

template<typename weight_t>
void BaseStencil<weight_t>::clear()
{
    int nBuckets = std::max(getBucketCount(), 1);

    clear(nBuckets);
}

template<typename weight_t>
void BaseStencil<weight_t>::clear(int nBuckets)
{
    m_pattern.initialize(nBuckets, 0, -1);
    m_weights.initialize(nBuckets, 0, m_zero);
    m_constant = m_zero;
}

template<typename weight_t>
void BaseStencil<weight_t>::flatten()
{
    m_pattern.merge();
    m_weights.merge();
}

template<typename weight_t>
void BaseStencil<weight_t>::optimize(double tolerance)
{
    int nBuckets = getBucketCount();
    for (int i = 0; i < nBuckets; ++i) {
        int nBucketItems = size(i);
        for (int j = 0; j < nBucketItems; ++j) {
            if (isWeightNeglibile(i, j, tolerance)) {
                m_pattern.eraseItem(i, j);
                m_weights.eraseItem(i, j);
                --j;
                --nBucketItems;
            }
        }
    }
}

template<typename weight_t>
void BaseStencil<weight_t>::renumber(const std::unordered_map<long, long> &map)
{
    int nWeights = size();
    for (int k = 0; k < nWeights; ++k) {
        long id = m_pattern.rawGetItem(k);
        id = map.at(id);
    }
}

template<typename weight_t>
void BaseStencil<weight_t>::addComplementToZero(const long id)
{
    int nWeights = size();
    if (nWeights == 0) {
        return;
    }

    weight_t complement = -1. * m_weights.rawGetItem(0);
    for (int n = 1; n < nWeights; ++n) {
        complement -= m_weights.rawGetItem(n);
    }

    addWeight(id, complement);
}

template<typename weight_t>
template<typename U, typename std::enable_if<std::is_fundamental<U>::value>::type*>
bool BaseStencil<weight_t>::isWeightNeglibile(int bucket, size_t k, double tolerance)
{
    return (std::abs(m_weights.getItem(bucket, k)) <= tolerance);
}

template<typename weight_t>
template<typename U, typename std::enable_if<!std::is_fundamental<U>::value>::type*>
bool BaseStencil<weight_t>::isWeightNeglibile(int bucket, size_t k, double tolerance)
{
    return (norm2(m_weights.getItem(bucket, k)) <= tolerance);
}

template<typename weight_t>
weight_t * BaseStencil<weight_t>::find(int bucket, long id) const
{
    int nBucketItems = size(bucket);
    for (int j = 0; j < nBucketItems; ++j) {
        long guessId = m_pattern.get(bucket) + j;
        if (guessId == id) {
            weight_t *weight = m_weights.get(bucket) + j;
            return weight;
        }
    }

    return nullptr;
}

template<typename weight_t>
void BaseStencil<weight_t>::display(std::ostream &out, double factor) const
{
    int nBuckets = getBucketCount();

    weight_t sum = m_zero;
    for (int i = 0; i < nBuckets; ++i) {
        int nBucketItems = size(i);
        out << " bucket : " << i << " n. bucket items : " << nBucketItems <<  std::endl;
        for (int j = 0; j < nBucketItems; ++j) {
            long id = m_pattern.getItem(i, j);
            weight_t weight = m_weights.getItem(i, j);
            out << "   id: " << id << " weight: " << (factor * weight) << std::endl;
            sum += factor * weight;
        }
    }

    out << " constant : " << factor *m_constant << std::endl;
    out << " sum      : " << sum << std::endl;
}

/*!
    Returns the buffer size (in bytes) required to store the stencil.

    \result The buffer size (in bytes) required to store the stencil.
*/
template<typename weight_t>
size_t BaseStencil<weight_t>::getBinarySize() const
{
    return (sizeof(m_zero) + m_pattern.getBinarySize() + m_weights.getBinarySize() + sizeof(m_constant));
}

template<typename weight_t>
BaseStencil<weight_t>& BaseStencil<weight_t>::operator*=(double factor)
{
    for (int k = 0; k < m_weights.getItemCount(); ++k) {
        weight_t &weight = m_weights.rawGetItem(k);
        weight *= factor;
    }
    m_constant *= factor;

    return *this;
}

template<typename weight_t>
BaseStencil<weight_t>& BaseStencil<weight_t>::operator/=(double factor)
{
    for (int k = 0; k < m_weights.getItemCount(); ++k) {
        weight_t &weight = m_weights.rawGetItem(k);
        weight /= factor;
    }
    m_constant /= factor;

    return *this;
}

template<typename weight_t>
BaseStencil<weight_t>& BaseStencil<weight_t>::operator+=(const BaseStencil<weight_t> &other)
{
    int nBuckets = getBucketCount();
    int other_nBuckets = other.getBucketCount();
    if (nBuckets != other_nBuckets) {
        throw std::runtime_error("Stencil must have the same number of buckets.");
    }

    for (int i = 0; i < other_nBuckets; ++i) {
        const int other_nBucketItems = other.m_weights.getItemCount(i);
        for (int j = 0; j < other_nBucketItems; ++j) {
            long id = other.m_pattern.get(i) + j;

            weight_t *weight = find(i, id);
            weight_t *other_weight = other.m_weights.get(i) + j;
            if (weight) {
                weight += other_weight;
            } else {
                addWeight(i, id, other_weight);
            }
        }
    }

    m_constant += other.m_constant;

    return *this;
}

template<typename weight_t>
BaseStencil<weight_t>& BaseStencil<weight_t>::operator-=(const BaseStencil<weight_t> &other)
{
    int nBuckets = getBucketCount();
    int other_nBuckets = other.getBucketCount();
    if (nBuckets != other_nBuckets) {
        throw std::runtime_error("Stencil must have the same number of buckets.");
    }

    for (int i = 0; i < other_nBuckets; ++i) {
        const int other_nBucketItems = other.m_weights.getItemCount(i);
        for (int j = 0; j < other_nBucketItems; ++j) {
            long id = other.m_pattern.get(i) + j;

            weight_t *weight = find(i, id);
            weight_t *other_weight = other.m_weights.get(i) + j;
            if (weight) {
                weight -= other_weight;
            } else {
                addWeight(i, id, -1. * other_weight);
            }
        }
    }

    m_constant -= other.m_constant;

    return *this;
}

}

template<typename weight_t>
bitpit::BaseStencil<weight_t> operator*(const bitpit::BaseStencil<weight_t> &stencil_A, double factor)
{
    return (factor * stencil_A);
}

template<typename weight_t>
bitpit::BaseStencil<weight_t> operator*(double factor, const bitpit::BaseStencil<weight_t> &stencil_A)
{
    bitpit::BaseStencil<weight_t> stencil_B(stencil_A);
    stencil_B *= factor;

    return stencil_B;
}

template<typename weight_t>
bitpit::BaseStencil<weight_t> operator/(const bitpit::BaseStencil<weight_t> &stencil_A, double factor)
{
    bitpit::BaseStencil<weight_t> stencil_B(stencil_A);
    stencil_B /= factor;

    return stencil_B;
}

template<typename weight_t>
bitpit::BaseStencil<weight_t> operator+(const bitpit::BaseStencil<weight_t> &stencil_A, const bitpit::BaseStencil<weight_t> &value)
{
    bitpit::BaseStencil<weight_t> stencil_B(stencil_A);
    stencil_B += value;

    return stencil_B;
}

template<typename weight_t>
bitpit::BaseStencil<weight_t> operator-(const bitpit::BaseStencil<weight_t> &stencil_A, const bitpit::BaseStencil<weight_t> &value)
{
    bitpit::BaseStencil<weight_t> stencil_B(stencil_A);
    stencil_B -= value;

    return stencil_B;
}
#endif
