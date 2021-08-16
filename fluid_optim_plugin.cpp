/* Author: Jose Pinzon
   Source: https://github.com/MultigridShapeOpt
  *
  * This file is a part of the FluidOptim UG4 plugin under development at 
  * the Research Group Approximation and Optimization, Hamburg University
  * and as part of the project SENSUS (LFF-GK11).
  *
  * This library is free software; you can redistribute it and/or
  * modify it under the terms of the GNU Lesser General Public
  * License as published by the Free Software Foundation; either
  * version 2.1 of the License, or (at your option) any later version.
  *
  * This library is distributed in the hope that it will be useful,
  * but WITHOUT ANY WARRANTY; without even the implied warranty of
  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  * Lesser General Public License for more details.
*/

#include "bridge/util.h"
#include "bridge/util_domain_dependent.h"
#include "bridge/util_domain_algebra_dependent.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_algebra/operator/debug_writer.h"

#include "fluid_optim_tools.h"
#include "gf_tools.h"
#include "DeformationStateEquation.h"
#include "TransformedNavierStokes.h"
#include "BoundaryControl.h"
#include "AdjointSystem.h"
#include "DeformationAdjointSystem.h"
#include "DesignEquation.h"
#include "SurfaceDesignEquation.h"
#include "SurfaceRHS.h"
#include "SurfaceBoundaryControl.h"
#include "ExtensionEquation.h"


using namespace std;
using namespace ug::bridge;

namespace ug{
namespace FluidOptim{

struct Functionality
{

/**
 * Function called for the registration of Domain dependent parts
 * of the plugin. All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
	static const int dim = TDomain::dim;
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();

//	TransformedNavierStokes Class with debugging functions, Non-linear
	{
		typedef TransformedNavierStokes<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("TransformedNavierStokes").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Functions#Subset(s)")
			.add_method("set_nonlinear", &T::set_nonlinear, "", "Non-linear problem")
			.add_method("set_picard", &T::set_picard,"","Picard Iteration")
			.add_method("no_deformation", &T::no_deformation, "", "Identity DeformationGradient")
			.add_method("use_terms_on_current_domain", &T::use_terms_on_current_domain,"","Use current domain, with deformation")
			.add_method("normal_diffusion", &T::normal_diffusion, "", "Test")
			.add_method("normal_pressure", &T::normal_pressure, "", "Test")
			.add_method("normal_continuity", &T::normal_continuity, "", "Test")
			.add_method("normal_convection", &T::normal_convection, "", "Test")
			.add_method("normal_stabilization", &T::normal_stabilization, "", "Test")
			.add_method("normal_jacobian_diffusion", &T::normal_jacobian_diffusion, "", "Test")
			.add_method("normal_jacobian_pressure", &T::normal_jacobian_pressure, "", "Test")
			.add_method("normal_jacobian_continuity", &T::normal_jacobian_continuity, "", "Test")
			.add_method("normal_jacobian_stabilization", &T::normal_jacobian_stabilization, "", "Test")
			.add_method("normal_jacobian_vectorconvection", &T::normal_jacobian_vectorconvection, "", "Test")
			.add_method("normal_jacobian_newtonderivative", &T::normal_jacobian_newtonderivative, "", "Test")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(number)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
	#ifdef UG_FOR_LUA
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(const char*)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
	#endif
			.add_method("set_deformation_vector_d1", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_deformation_vector_d1),"","DeformationVectord1")
			.add_method("set_deformation_vector_d1", static_cast<void (T::*)(number)>(&T::set_deformation_vector_d1),"","DeformationVectord1")
	#ifdef UG_FOR_LUA
			.add_method("set_deformation_vector_d1", static_cast<void (T::*)(const char*)>(&T::set_deformation_vector_d1),"","DeformationVectord1")
			.add_method("set_deformation_vector_d1", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_deformation_vector_d1),"","DeformationVectord1")
	#endif
			.add_method("set_deformation_vector_d2", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_deformation_vector_d2),"","DeformationVectord2")
			.add_method("set_deformation_vector_d2", static_cast<void (T::*)(number)>(&T::set_deformation_vector_d2),"","DeformationVectord2")
	#ifdef UG_FOR_LUA
			.add_method("set_deformation_vector_d2", static_cast<void (T::*)(const char*)>(&T::set_deformation_vector_d2),"","DeformationVectord2")
			.add_method("set_deformation_vector_d2", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_deformation_vector_d2),"","DeformationVectord2")
	#endif
			.add_method("set_deformation_vector_d3", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_deformation_vector_d3),"","DeformationVectord3")
			.add_method("set_deformation_vector_d3", static_cast<void (T::*)(number)>(&T::set_deformation_vector_d3),"","DeformationVectord3")
			.add_method("set_compare_value", &T::set_compare_value)
			.add_method("compare_deformation", &T::compare_deformation)
			.add_method("set_quad_order", &T::set_quad_order)
			.add_method("set_stabilization", &T::set_stabilization)
			.add_method("set_stabilization_type", &T::set_stabilization_type)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name,"TransformedNavierStokes",tag);
	}	
//	DeformationStateEquation Class, Non-linear
	{
		typedef DeformationStateEquation<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("DeformationStateEquation").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.add_method("set_quad_order", &T::set_quad_order, "", "quad order")
			.add_method("set_nonlinear", &T::set_nonlinear, "", "Non-linear problem")
			.add_method("set_picard", &T::set_picard,"","Picard Iteration")
			.add_method("set_symmetric", &T::set_symmetric,"","Symmetric Gradient")
			.add_method("set_extension_factor", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_extension_factor),"","Extension Factor")
			.add_method("set_extension_factor", static_cast<void (T::*)(number)>(&T::set_extension_factor),"","Extension Factor")
			//.add_method("set_nu_extension", &T::set_nu_extension,"","Extension Factor")
			//.add_method("set_deformation_vector_d1", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_deformation_vector_d1),"","DeformationVectord1")
			//.add_method("set_deformation_vector_d1", static_cast<void (T::*)(number)>(&T::set_deformation_vector_d1),"","DeformationVectord1")
			//.add_method("set_deformation_vector_d2", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_deformation_vector_d2),"","DeformationVectord2")
			//.add_method("set_deformation_vector_d2", static_cast<void (T::*)(number)>(&T::set_deformation_vector_d2),"","DeformationVectord2")
			//.add_method("vector", &T::vector)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "DeformationStateEquation", tag);
	}
//	Boundary Control rhs provider
	{
		typedef BoundaryControl<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("BoundaryControl").append(suffix);
		reg.add_class_<T,TBase>(name,grp)
			.template add_constructor<void (*)(const char*)>("Function(s)")
			.add_method("set_fixed_boundary_control", static_cast<void (T::*)(number, const char*, const char*)>(&T::set_fixed_boundary_control))
			.add_method("set_boundary_inner_subsets", &T::set_boundary_inner_subsets,"","Boundary and Inner subset")
			.add_method("set_variable_boundary_control", &T::set_variable_boundary_control, "", "Lua BoundaryControl")
			.add_method("set_boundary_control", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_boundary_control), "", "Lua BoundaryControl")
			.add_method("set_boundary_control", static_cast<void (T::*)(number)>(&T::set_boundary_control), "", "Lua BoundaryControl")
#ifdef UG_FOR_LUA
			.add_method("set_boundary_control", static_cast<void (T::*)(const char*)>(&T::set_boundary_control), "", "Lua BoundaryControl")
			.add_method("set_boundary_control", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_boundary_control), "", "Lua BoundaryControl")
#endif	
			.add_method("check_boundary_control", &T::check_boundary_control,"","Check that import is set")
			//.add_method("check_value_boundary_control", &T::check_value_boundary_control,"","Check Value of Import")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name,"BoundaryControl",tag);
	}
//	AdjointSystem Class, Linear
	{
		typedef AdjointSystem<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("AdjointSystem").append(suffix);
		reg.add_class_<T,TBase>(name,grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Functions#Subset(s)")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(number)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
	#ifdef UG_FOR_LUA
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(const char*)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
	#endif
			.add_method("set_deformation_vector_d1", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_deformation_vector_d1),"","DeformationVectord1")
			.add_method("set_deformation_vector_d1", static_cast<void (T::*)(number)>(&T::set_deformation_vector_d1),"","DeformationVectord1")
	#ifdef UG_FOR_LUA
			.add_method("set_deformation_vector_d1", static_cast<void (T::*)(const char*)>(&T::set_deformation_vector_d1),"","DeformationVectord1")
			.add_method("set_deformation_vector_d1", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_deformation_vector_d1),"","DeformationVectord1")
	#endif
			.add_method("set_deformation_vector_d2", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_deformation_vector_d2),"","DeformationVectord2")
			.add_method("set_deformation_vector_d2", static_cast<void (T::*)(number)>(&T::set_deformation_vector_d2),"","DeformationVectord2")
	#ifdef UG_FOR_LUA
			.add_method("set_deformation_vector_d2", static_cast<void (T::*)(const char*)>(&T::set_deformation_vector_d2),"","DeformationVectord2")
			.add_method("set_deformation_vector_d2", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_deformation_vector_d2),"","DeformationVectord2")
	#endif
			.add_method("set_deformation_vector_d3", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_deformation_vector_d3),"","DeformationVectord3")
			.add_method("set_deformation_vector_d3", static_cast<void (T::*)(number)>(&T::set_deformation_vector_d3),"","DeformationVectord3")
			.add_method("set_velocity_vector_d1", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_velocity_vector_d1),"","VelocityGradientd1")
			.add_method("set_velocity_vector_d1", static_cast<void (T::*)(number)>(&T::set_velocity_vector_d1),"","VelocityGradientd1")
			.add_method("set_velocity_vector_d2", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_velocity_vector_d2),"","VelocityGradientd2")
			.add_method("set_velocity_vector_d2", static_cast<void (T::*)(number)>(&T::set_velocity_vector_d2),"","VelocityGradientd2")
			.add_method("set_velocity_vector_d3", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_velocity_vector_d3),"","VelocityGradientd3")
			.add_method("set_velocity_vector_d3", static_cast<void (T::*)(number)>(&T::set_velocity_vector_d3),"","VelocityGradientd3")
			.add_method("set_velocity_d1", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_velocity_d1),"","VelocityValued1")
			.add_method("set_velocity_d1", static_cast<void (T::*)(number)>(&T::set_velocity_d1),"","VelocityValued1")
			.add_method("set_velocity_d2", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_velocity_d2),"","VelocityValued2")
			.add_method("set_velocity_d2", static_cast<void (T::*)(number)>(&T::set_velocity_d2),"","VelocityValued2")
			.add_method("set_velocity_d3", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_velocity_d3),"","VelocityValued3")
			.add_method("set_velocity_d3", static_cast<void (T::*)(number)>(&T::set_velocity_d3),"","VelocityValued3")
			.add_method("set_stabilization", &T::set_stabilization)
			.add_method("set_stabilization_type", &T::set_stabilization_type)
			.add_method("set_quad_order", &T::set_quad_order, "", "quad order")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name,"AdjointSystem",tag);
		
	}
	//	DeformationAdjointSystem Class, Linear
	{
		typedef DeformationAdjointSystem<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("DeformationAdjointSystem").append(suffix);
		reg.add_class_<T,TBase>(name,grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Functions#Subset(s)")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(number)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
	#ifdef UG_FOR_LUA
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(const char*)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
	#endif
			.add_method("set_deformation_vector_d1", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_deformation_vector_d1),"","DeformationVectord1")
			.add_method("set_deformation_vector_d1", static_cast<void (T::*)(number)>(&T::set_deformation_vector_d1),"","DeformationVectord1")
			.add_method("set_deformation_vector_d2", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_deformation_vector_d2),"","DeformationVectord2")
			.add_method("set_deformation_vector_d2", static_cast<void (T::*)(number)>(&T::set_deformation_vector_d2),"","DeformationVectord2")
			.add_method("set_deformation_vector_d3", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_deformation_vector_d3),"","DeformationVectord3")
			.add_method("set_deformation_vector_d3", static_cast<void (T::*)(number)>(&T::set_deformation_vector_d3),"","DeformationVectord3")
			.add_method("set_deformation_d1", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_deformation_d1),"","DeformationValued1")
			.add_method("set_deformation_d1", static_cast<void (T::*)(number)>(&T::set_deformation_d1),"","DeformationValued1")
			.add_method("set_deformation_d2", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_deformation_d2),"","DeformationValued2")
			.add_method("set_deformation_d2", static_cast<void (T::*)(number)>(&T::set_deformation_d2),"","DeformationValued2")
			.add_method("set_deformation_d3", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_deformation_d3),"","DeformationValued3")
			.add_method("set_deformation_d3", static_cast<void (T::*)(number)>(&T::set_deformation_d3),"","DeformationValued3")
			.add_method("set_velocity_vector_d1", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_velocity_vector_d1),"","VelocityGradientd1")
			.add_method("set_velocity_vector_d1", static_cast<void (T::*)(number)>(&T::set_velocity_vector_d1),"","VelocityGradientd1")
			.add_method("set_velocity_vector_d2", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_velocity_vector_d2),"","VelocityGradientd2")
			.add_method("set_velocity_vector_d2", static_cast<void (T::*)(number)>(&T::set_velocity_vector_d2),"","VelocityGradientd2")
			.add_method("set_velocity_vector_d3", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_velocity_vector_d3),"","VelocityGradientd3")
			.add_method("set_velocity_vector_d3", static_cast<void (T::*)(number)>(&T::set_velocity_vector_d3),"","VelocityGradientd3")
			.add_method("set_velocity_d1", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_velocity_d1),"","VelocityValued1")
			.add_method("set_velocity_d1", static_cast<void (T::*)(number)>(&T::set_velocity_d1),"","VelocityValued1")
			.add_method("set_velocity_d2", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_velocity_d2),"","VelocityValued2")
			.add_method("set_velocity_d2", static_cast<void (T::*)(number)>(&T::set_velocity_d2),"","VelocityValued2")
			.add_method("set_velocity_d3", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_velocity_d3),"","VelocityValued3")
			.add_method("set_velocity_d3", static_cast<void (T::*)(number)>(&T::set_velocity_d3),"","VelocityValued3")
			.add_method("set_adjoint_velocity_vector_d1", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_adjoint_velocity_vector_d1),"","AdjointVelocityGradientd1")
			.add_method("set_adjoint_velocity_vector_d1", static_cast<void (T::*)(number)>(&T::set_adjoint_velocity_vector_d1),"","AdjointVelocityGradientd1")
			.add_method("set_adjoint_velocity_vector_d2", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_adjoint_velocity_vector_d2),"","AdjointVelocityGradientd2")
			.add_method("set_adjoint_velocity_vector_d2", static_cast<void (T::*)(number)>(&T::set_adjoint_velocity_vector_d2),"","AdjointVelocityGradientd2")
			.add_method("set_adjoint_velocity_vector_d3", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_adjoint_velocity_vector_d3),"","AdjointVelocityGradientd3")
			.add_method("set_adjoint_velocity_vector_d3", static_cast<void (T::*)(number)>(&T::set_adjoint_velocity_vector_d3),"","AdjointVelocityGradientd3")
			.add_method("set_pressure", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_pressure),"","PressureScalar")
			.add_method("set_pressure", static_cast<void (T::*)(number)>(&T::set_pressure),"","PressureScalar")
			.add_method("set_adjoint_pressure", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_adjoint_pressure),"","PressureScalar")
			.add_method("set_adjoint_pressure", static_cast<void (T::*)(number)>(&T::set_adjoint_pressure),"","AdjointPressureScalar")
			.add_method("set_adjoint_velocity_d1", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_adjoint_velocity_d1),"","VelocityValued1")
			.add_method("set_adjoint_velocity_d1", static_cast<void (T::*)(number)>(&T::set_adjoint_velocity_d1),"","VelocityValued1")
			.add_method("set_adjoint_velocity_d2", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_adjoint_velocity_d2),"","VelocityValued2")
			.add_method("set_adjoint_velocity_d2", static_cast<void (T::*)(number)>(&T::set_adjoint_velocity_d2),"","VelocityValued2")
			.add_method("set_adjoint_velocity_d3", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_adjoint_velocity_d3),"","VelocityValued3")
			.add_method("set_adjoint_velocity_d3", static_cast<void (T::*)(number)>(&T::set_adjoint_velocity_d3),"","VelocityValued3")
			.add_method("set_extension_factor", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_extension_factor),"","Extension Factor")
			.add_method("set_extension_factor", static_cast<void (T::*)(number)>(&T::set_extension_factor),"","Extension Factor")
			.add_method("set_determinant_threshold", &T::set_determinant_threshold)
			.add_method("set_beta", &T::set_beta)
			.add_method("set_lambda_vol", &T::set_lambda_vol)
			.add_method("set_lambda_barycenter", &T::set_lambda_barycenter)
			.add_method("set_nu_penalty", &T::set_nu_penalty,"","Geometrical Penalty")
			.add_method("set_volume_defect", &T::set_volume_defect,"","Geometrical Penalty")
			.add_method("set_barycenter_defect", &T::set_barycenter_defect,"","Geometrical Penalty")
			.add_method("set_quad_order", &T::set_quad_order, "", "Quad order")
			.add_method("set_stabilization_scale", &T::set_stabilization_scale, "", "Scaling of stabilization sensitivity")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name,"DeformationAdjointSystem",tag);
		
	}	
	//	Design Equation: Gradient of Boundary Control variable, linear problem
	{
		typedef DesignEquation<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("DesignEquation").append(suffix);
		reg.add_class_<T,TBase>(name,grp)
			.template add_constructor<void (*)(const char*)>("Function(s)")
			.add_method("set_boundary_inner_subsets", &T::set_boundary_inner_subsets,"","Boundary and Inner subset")
			.add_method("set_boundary_control", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_boundary_control), "", "Lua BoundaryControl")
			.add_method("set_boundary_control", static_cast<void (T::*)(number)>(&T::set_boundary_control), "", "Lua BoundaryControl")
			.add_method("set_adjoint_deformation_d1", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_adjoint_deformation_d1),"","DeformationVectord1")
			.add_method("set_adjoint_deformation_d1", static_cast<void (T::*)(number)>(&T::set_adjoint_deformation_d1),"","DeformationVectord1")
			.add_method("set_adjoint_deformation_d2", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_adjoint_deformation_d2),"","DeformationVectord2")
			.add_method("set_adjoint_deformation_d2", static_cast<void (T::*)(number)>(&T::set_adjoint_deformation_d2),"","DeformationVectord2")
			.add_method("set_adjoint_deformation_d3", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_adjoint_deformation_d3),"","DeformationVectord3")
			.add_method("set_adjoint_deformation_d3", static_cast<void (T::*)(number)>(&T::set_adjoint_deformation_d3),"","DeformationVectord3")
			.add_method("set_quad_order", &T::set_quad_order, "", "quad order")
			.add_method("set_alpha", &T::set_alpha, "", "regularization parameter")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name,"DesignEquation",tag);
	}
	//	Surface Design Equation: provides only a mass-matrix, ideally on a surface subset
	{
		typedef SurfaceDesignEquation<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("SurfaceDesignEquation").append(suffix);
		reg.add_class_<T,TBase>(name,grp)
			.template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset")
			.add_method("set_boundary_control", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_boundary_control), "", "Lua BoundaryControl")
			.add_method("set_boundary_control", static_cast<void (T::*)(number)>(&T::set_boundary_control), "", "Lua BoundaryControl")
			.add_method("set_adjoint_deformation_d1", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_adjoint_deformation_d1),"","DeformationVectord1")
			.add_method("set_adjoint_deformation_d1", static_cast<void (T::*)(number)>(&T::set_adjoint_deformation_d1),"","DeformationVectord1")
			.add_method("set_adjoint_deformation_d2", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_adjoint_deformation_d2),"","DeformationVectord2")
			.add_method("set_adjoint_deformation_d2", static_cast<void (T::*)(number)>(&T::set_adjoint_deformation_d2),"","DeformationVectord2")
			.add_method("set_adjoint_deformation_d3", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_adjoint_deformation_d3),"","DeformationVectord3")
			.add_method("set_adjoint_deformation_d3", static_cast<void (T::*)(number)>(&T::set_adjoint_deformation_d3),"","DeformationVectord3")
			.add_method("set_quad_order", &T::set_quad_order, "", "quad order")
			.add_method("set_alpha", &T::set_alpha, "", "regularization parameter")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name,"SurfaceDesignEquation",tag);
	}
	//	Surface RHS:provides rhs on surface of obstacle
	{
		typedef SurfaceRHS<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("SurfaceRHS").append(suffix);
		reg.add_class_<T,TBase>(name,grp)
			.template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset")
			.add_method("set_boundary_control", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_boundary_control), "", "Lua BoundaryControl")
			.add_method("set_boundary_control", static_cast<void (T::*)(number)>(&T::set_boundary_control), "", "Lua BoundaryControl")
			.add_method("set_adjoint_deformation_d1", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_adjoint_deformation_d1),"","DeformationVectord1")
			.add_method("set_adjoint_deformation_d1", static_cast<void (T::*)(number)>(&T::set_adjoint_deformation_d1),"","DeformationVectord1")
			.add_method("set_adjoint_deformation_d2", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_adjoint_deformation_d2),"","DeformationVectord2")
			.add_method("set_adjoint_deformation_d2", static_cast<void (T::*)(number)>(&T::set_adjoint_deformation_d2),"","DeformationVectord2")
			.add_method("set_adjoint_deformation_d3", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_adjoint_deformation_d3),"","DeformationVectord3")
			.add_method("set_adjoint_deformation_d3", static_cast<void (T::*)(number)>(&T::set_adjoint_deformation_d3),"","DeformationVectord3")
			.add_method("set_quad_order", &T::set_quad_order, "", "quad order")
			.add_method("set_alpha", &T::set_alpha, "", "regularization parameter")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name,"SurfaceRHS",tag);
	}
		//	SurfaceBoundaryControl :provides rhs on surface of obstacle for the deformation equation
	{
		typedef SurfaceBoundaryControl<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("SurfaceBoundaryControl").append(suffix);
		reg.add_class_<T,TBase>(name,grp)
			.template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset")
			.add_method("set_boundary_control", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_boundary_control), "", "Lua BoundaryControl")
			.add_method("set_boundary_control", static_cast<void (T::*)(number)>(&T::set_boundary_control), "", "Lua BoundaryControl")
			.add_method("set_quad_order", &T::set_quad_order, "", "quad order")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name,"SurfaceBoundaryControl",tag);
	}
	//	Extension Equation
	{
		typedef ExtensionEquation<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("ExtensionEquation").append(suffix);
		reg.add_class_<T,TBase>(name,grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset")
			.add_method("set_deformation_vector_d1", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_deformation_vector_d1),"","DeformationVectord1")
			.add_method("set_deformation_vector_d1", static_cast<void (T::*)(number)>(&T::set_deformation_vector_d1),"","DeformationVectord1")
			.add_method("set_deformation_vector_d2", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_deformation_vector_d2),"","DeformationVectord2")
			.add_method("set_deformation_vector_d2", static_cast<void (T::*)(number)>(&T::set_deformation_vector_d2),"","DeformationVectord2")
			.add_method("set_deformation_vector_d3", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_deformation_vector_d3),"","DeformationVectord3")
			.add_method("set_deformation_vector_d3", static_cast<void (T::*)(number)>(&T::set_deformation_vector_d3),"","DeformationVectord3")
			.add_method("set_deformation_d1", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_deformation_d1),"","DeformationValued1")
			.add_method("set_deformation_d1", static_cast<void (T::*)(number)>(&T::set_deformation_d1),"","DeformationValued1")
			.add_method("set_deformation_d2", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_deformation_d2),"","DeformationValued2")
			.add_method("set_deformation_d2", static_cast<void (T::*)(number)>(&T::set_deformation_d2),"","DeformationValued2")
			.add_method("set_deformation_d3", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_deformation_d3),"","DeformationValued3")
			.add_method("set_deformation_d3", static_cast<void (T::*)(number)>(&T::set_deformation_d3),"","DeformationValued3")
			.add_method("set_adjoint_deformation_d1", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_adjoint_deformation_d1),"","DeformationVectord1")
			.add_method("set_adjoint_deformation_d1", static_cast<void (T::*)(number)>(&T::set_adjoint_deformation_d1),"","DeformationVectord1")
			.add_method("set_adjoint_deformation_d2", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_adjoint_deformation_d2),"","DeformationVectord2")
			.add_method("set_adjoint_deformation_d2", static_cast<void (T::*)(number)>(&T::set_adjoint_deformation_d2),"","DeformationVectord2")
			.add_method("set_adjoint_deformation_d3", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_adjoint_deformation_d3),"","DeformationVectord3")
			.add_method("set_adjoint_deformation_d3", static_cast<void (T::*)(number)>(&T::set_adjoint_deformation_d3),"","DeformationVectord3")
			.add_method("set_extension_factor", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_extension_factor),"","Extension Factor")
			.add_method("set_extension_factor", static_cast<void (T::*)(number)>(&T::set_extension_factor),"","Extension Factor")
			.add_method("set_quad_order", &T::set_quad_order, "", "quad order")
			.add_method("set_extension_upper", &T::set_extension_upper, "", "upper extension factor limit")			
			.add_method("set_extension_lower", &T::set_extension_lower, "", "lower extension factor limit")
			.add_method("set_chi", &T::set_chi, "", "chi value")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name,"ExtensionEquation",tag);
	}
}

template <int dim>
static void Dimension(Registry& reg, string grp)
{
	string dimSuffix = GetDimensionSuffix<dim>();
	string dimTag = GetDimensionTag<dim>();

	
}

template <typename TDomain, typename TAlgebra>
static void DomainAlgebra(Registry& reg, string grp)
{
	static const int dim = TDomain::dim;
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

	typedef ApproximationSpace<TDomain> approximation_space_type;
	typedef GridFunction<TDomain, TAlgebra> function_type;
	
	//volume defect
	{
		reg.add_function("VolumeDefect", &VolumeDefect<function_type>, grp);
	}
	//barycenter defect
	{
		reg.add_function("BarycenterDefect", &BarycenterDefect<function_type>, grp);
	}
	//Drag of objective function
	{
		reg.add_function("Drag", &Drag<function_type>, grp);
	}
	//Boundary control of objective function
	{
		reg.add_function("Control", &Control<function_type>, grp);
	}
	//Threshold of objective function
	{
		reg.add_function("Threshold", &Threshold<function_type>, grp);
	}
	//Threshold of objective function
	{
		reg.add_function("Extension", &Extension<function_type>, grp);
	}
	{
		reg.add_function("SurfaceIntegral", &SurfaceIntegral<function_type>, grp);
	}
	{
		reg.add_function("VolumeIntegral", &VolumeIntegral<function_type>, grp);
	}
 	//Limits of objective function
	{
		reg.add_function("LimitGF", &LimitGF<function_type>, grp);
	}
	{
		reg.add_function("L2VecProd", &L2VecProd<function_type>, grp);
	}
	{
		reg.add_function("ControlGradient", &ControlGradient<function_type>, grp);
	}
}

}; // end Functionality

/**
 * Class exporting the functionality of the plugin restricted to 2 and 3 spatial
 * dimensions. All functionality that is to be used in scripts or visualization
 * only in 2d and 3d must be registered here.
 */

} 


/**
 * This function is called when the plugin is loaded.
 */
extern "C" void
InitUGPlugin_FluidOptim(Registry* reg, string grp)
{
	grp.append("/SpatialDisc/ElemDisc");
	typedef FluidOptim::Functionality Functionality;

	try{
		
		RegisterDomain2d3dDependent<Functionality>(*reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(*reg,grp);
		RegisterDimensionDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// namespace ug
