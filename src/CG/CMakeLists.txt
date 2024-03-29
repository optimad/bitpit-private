#---------------------------------------------------------------------------
#
#  bitpit
#
#  Copyright (C) 2015-2019 OPTIMAD engineering Srl
#
#  -------------------------------------------------------------------------
#  License
#  This file is part of bitpit.
#
#  bitpit is free software: you can redistribute it and/or modify it
#  under the terms of the GNU Lesser General Public License v3 (LGPL)
#  as published by the Free Software Foundation.
#
#  bitpit is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
#
#---------------------------------------------------------------------------*/

# CMake settings
cmake_minimum_required(VERSION 2.8)
initialize_lto_property()

# Name of the current module
get_filename_component(MODULE_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)

# Add library to targets
addModuleIncludeDirectories(${MODULE_NAME})

file(GLOB SOURCE_FILES "*.cpp")
set(CG_SOURCES "${SOURCE_FILES}" CACHE INTERNAL "Sources of CG module" FORCE)
unset(SOURCE_FILES)

file(GLOB HEADER_FILES "*.hpp" "*.tpp")
set(CG_HEADERS "${HEADER_FILES}" CACHE INTERNAL "Headers of CG module" FORCE)
unset(HEADER_FILES)

if (NOT "${CG_SOURCES}" STREQUAL "")
    set(CG_TARGET_OBJECT "CG_TARGET_OBJECT")
    add_library(${CG_TARGET_OBJECT} OBJECT ${CG_SOURCES})
    set_lto_property(${CG_TARGET_OBJECT})
endif ()

install(FILES ${CG_HEADERS} DESTINATION include/${PROJECT_NAME})
