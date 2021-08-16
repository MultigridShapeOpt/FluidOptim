# Author: Jose Pinzon
# Source: https://github.com/MultigridShapeOpt
# *
# * This file is a part of the FluidOptim UG4 plugin under development at 
# * the Research Group Approximation and Optimization, Hamburg University
# * and as part of the project SENSUS (LFF-GK11).
# *
# * This library is free software; you can redistribute it and/or
# * modify it under the terms of the GNU Lesser General Public
# * License as published by the Free Software Foundation; either
# * version 2.1 of the License, or (at your option) any later version.
# *
# * This library is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# * Lesser General Public License for more details.



# Name of your plugin and sources.
set(pluginName	FluidOptim)
set(SOURCES		fluid_optim_plugin.cpp
				DeformationStateEquation.cpp
				TransformedNavierStokes.cpp
				BoundaryControl.cpp
				AdjointSystem.cpp
				DeformationAdjointSystem.cpp
				DesignEquation.cpp
				SurfaceDesignEquation.cpp
				SurfaceRHS.cpp
				ExtensionEquation.cpp
				SurfaceBoundaryControl.cpp
	)


################################################################################
# The code below doesn't have to be changed (usually)
################################################################################
cmake_minimum_required(VERSION 2.6)

project(UG_PLUGIN_${pluginName})

# include the definitions and dependencies for ug-plugins.
include(${UG_ROOT_CMAKE_PATH}/ug_plugin_includes.cmake)

if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(${CMAKE_CURRENT_SOURCE_DIR} ${SOURCES})
else(buildEmbeddedPlugins)
	# create a shared library from the sources and link it against ug4.
	add_library(${pluginName} SHARED ${SOURCES})
	target_link_libraries (${pluginName} ug4)
endif(buildEmbeddedPlugins)
