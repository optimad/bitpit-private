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

#include "volume_kernel.hpp"

namespace bitpit {

/*!
	\class VolumeKernel
	\ingroup volumepatches

	\brief The VolumeKernel class provides an interface for defining
	volume patches.

	VolumeKernel is the base class for defining voulme patches.
*/

/*!
	Creates a new patch.

	\param expert if true, the expert mode will be enabled
*/
VolumeKernel::VolumeKernel(bool expert)
	: PatchKernel(expert)
{
}

/*!
	Creates a new patch.

	\param dimension is the dimension of the patch
	\param expert if true, the expert mode will be enabled
*/
VolumeKernel::VolumeKernel(int dimension, bool expert)
	: PatchKernel(dimension, expert)
{
}

/*!
	Creates a new patch.

	\param id is the id that will be assigned to the patch
	\param dimension is the dimension of the patch
	\param expert if true, the expert mode will be enabled
*/
VolumeKernel::VolumeKernel(int id, int dimension, bool expert)
	: PatchKernel(id, dimension, expert)
{
}

/*!
	Destroys the patch.
*/
VolumeKernel::~VolumeKernel()
{

}

/*!
	Get the codimension of the patch in the volume space.

	\result The codimension of the patch in the volume space.
*/
int VolumeKernel::getVolumeCodimension() const
{
	return 0;
}

/*!
	Get the codimension of the patch in the surface space.

	\result The codimension of the patch in the surface space.
*/
int VolumeKernel::getSurfaceCodimension() const
{
	return -1;
}

/*!
	Get the codimension of the patch in the line space.

	\result The codimension of the patch in the line space.
*/
int VolumeKernel::getLineCodimension() const
{
	return -2;
}

/*!
	Get the codimension of the patch in the point space.

	\result The codimension of the patch in the point space.
*/
int VolumeKernel::getPointCodimension() const
{
	return -3;
}

/*!
	Checks if the specified point is inside the patch.

	\param[in] x is the x coordinate of the point
	\param[in] y is the y coordinate of the point
	\param[in] z is the z coordinate of the point
	\result Returns true if the point is inside the patch, false otherwise.
 */
bool VolumeKernel::isPointInside(double x, double y, double z) const
{
	return isPointInside({{x, y, z}});
}

/*!
	Checks if the specified point is inside a cell.

	\param[in] id is the index of the cells
	\param[in] x is the x coordinate of the point
	\param[in] y is the y coordinate of the point
	\param[in] z is the z coordinate of the point
	\result Returns true if the point is inside the cell, false otherwise.
 */
bool VolumeKernel::isPointInside(long id, double x, double y, double z) const
{
	return isPointInside(id, {{x, y, z}});
}

}
