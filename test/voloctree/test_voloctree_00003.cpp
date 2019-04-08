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

#include <array>
#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#endif

#include "bitpit_common.hpp"
#include "bitpit_IO.hpp"
#include "bitpit_voloctree.hpp"

using namespace bitpit;

/*!
* Subtest 001
*
* Testing point localization in a 2D patch.
*/
int subtest_001()
{
	std::array<double, 3> origin = {{0., 0., 0.}};
	double length = 20;
	double dh = 0.5;

	std::array<double, 3> point;
	std::vector<std::array<double, 3>> pointList;

	pointList.push_back(origin);
	for (int i = 0; i < 3; ++i) {
		point[i] = origin[i] + length / 4;
		pointList.push_back(point);
	}
	for (int i = 0; i < 3; ++i) {
		point[i] = origin[i] + length / 2;
		pointList.push_back(point);
	}
	for (int i = 0; i < 3; ++i) {
		point[i] = origin[i] + length / 2 + 0.01;
		pointList.push_back(point);
	}
	for (int i = 0; i < 3; ++i) {
		point[i] = origin[i] - length / 2;
		pointList.push_back(point);
	}
	for (int i = 0; i < 3; ++i) {
		point[i] = origin[i] - length / 4;
		pointList.push_back(point);
	}

	log::cout() << "  >> 2D octree patch" << "\n";

	VolOctree *patch_2D = new VolOctree(2, origin, length, dh);
	patch_2D->getVTK().setName("octree_uniform_patch_2D");
	patch_2D->buildInterfaces();
	patch_2D->update();

	log::cout() << "\n  >> 2D location test" << std::endl;
	log::cout() << std::endl;

	for (auto testPoint : pointList) {
		log::cout() << "Point [" << testPoint[0] << ", " << testPoint[1] << ", " << testPoint[2] << "] ";
		if (patch_2D->isPointInside(testPoint)) {
			log::cout() << " is inside the patch" << std::endl;
		} else {
			log::cout() << " is outside the patch" << std::endl;
		}

		log::cout() << "Point [" << testPoint[0] << ", " << testPoint[1] << ", " << testPoint[2] << "] ";
		log::cout() << " is inside the element " << patch_2D->locatePoint(testPoint) << std::endl;

		log::cout() << "Point [" << testPoint[0] << ", " << testPoint[1] << ", " << testPoint[2] << "] ";
		if (patch_2D->isPointInside(testPoint[0], testPoint[1], testPoint[2])) {
			log::cout() << " is inside the patch" << std::endl;
		} else {
			log::cout() << " is outside the patch" << std::endl;
		}

		log::cout() << "Point [" << testPoint[0] << ", " << testPoint[1] << ", " << testPoint[2] << "] ";
		log::cout() << " is inside the element " << patch_2D->locatePoint(testPoint[0], testPoint[1], testPoint[2]) << std::endl;
	}

	log::cout() << std::endl;

	delete patch_2D;

	return 0;
}

/*!
* Subtest 002
*
* Testing point localization in a 3D patch.
*/
int subtest_002()
{
	std::array<double, 3> origin = {{0., 0., 0.}};
	double length = 20;
	double dh = 0.5;

	std::array<double, 3> point;
	std::vector<std::array<double, 3>> pointList;

	pointList.push_back(origin);
	for (int i = 0; i < 3; ++i) {
		point[i] = origin[i] + length / 4;
		pointList.push_back(point);
	}
	for (int i = 0; i < 3; ++i) {
		point[i] = origin[i] + length / 2;
		pointList.push_back(point);
	}
	for (int i = 0; i < 3; ++i) {
		point[i] = origin[i] + length / 2 + 0.01;
		pointList.push_back(point);
	}
	for (int i = 0; i < 3; ++i) {
		point[i] = origin[i] - length / 2;
		pointList.push_back(point);
	}
	for (int i = 0; i < 3; ++i) {
		point[i] = origin[i] - length / 4;
		pointList.push_back(point);
	}

	log::cout() << "  >> 3D octree patch" << "\n";

	VolOctree *patch_3D = new VolOctree(3, origin, length, dh);
	patch_3D->getVTK().setName("octree_uniform_patch_3D");
	patch_3D->buildInterfaces();
	patch_3D->update();

	log::cout() << "\n  >> 3D location test" << std::endl;
	log::cout() << std::endl;

	for (auto testPoint : pointList) {
		log::cout() << "Point [" << testPoint[0] << ", " << testPoint[1] << ", " << testPoint[2] << "] ";
		if (patch_3D->isPointInside(testPoint)) {
			log::cout() << " is inside the patch" << std::endl;
		} else {
			log::cout() << " is outside the patch" << std::endl;
		}

		log::cout() << "Point [" << testPoint[0] << ", " << testPoint[1] << ", " << testPoint[2] << "] ";
		log::cout() << " is inside the element " << patch_3D->locatePoint(testPoint) << std::endl;

		log::cout() << "Point [" << testPoint[0] << ", " << testPoint[1] << ", " << testPoint[2] << "] ";
		if (patch_3D->isPointInside(testPoint[0], testPoint[1], testPoint[2])) {
			log::cout() << " is inside the patch" << std::endl;
		} else {
			log::cout() << " is outside the patch" << std::endl;
		}

		log::cout() << "Point [" << testPoint[0] << ", " << testPoint[1] << ", " << testPoint[2] << "] ";
		log::cout() << " is inside the element " << patch_3D->locatePoint(testPoint[0], testPoint[1], testPoint[2]) << std::endl;
	}

	log::cout() << std::endl;

	delete patch_3D;

	return 0;
}

/*!
* Main program.
*/
int main(int argc, char *argv[])
{
#if BITPIT_ENABLE_MPI==1
	MPI_Init(&argc,&argv);
#else
	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);
#endif

	// Initialize the logger
	log::manager().initialize(log::COMBINED);

	// Seed the random function
	std::srand(1);

	// Run the subtests
	log::cout() << "Testing point localization in octree patches" << std::endl;

	int status;
	try {
		status = subtest_001();
		if (status != 0) {
			return status;
		}

		status = subtest_002();
		if (status != 0) {
			return status;
		}
	} catch (const std::exception &exception) {
		log::cout() << exception.what();
		exit(1);
	}

#if BITPIT_ENABLE_MPI==1
	MPI_Finalize();
#endif
}
