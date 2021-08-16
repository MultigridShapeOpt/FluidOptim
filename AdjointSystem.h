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
 
#ifndef ADJOINT_SYSTEM
#define ADJOINT_SYSTEM

#include <stdio.h>
#include <string>

#include "bridge/bridge.h"

#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/disc_util/fe_geom.h"
#include "lib_disc/local_finite_element/lagrange/lagrange.h"
#include "lib_disc/quadrature/gauss/gauss_quad.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"


namespace ug{
	namespace FluidOptim{
		
		template<typename TDomain>
		class AdjointSystem: public IElemDisc<TDomain>
		{
			//CONSTRUCTOR's
			public:
				AdjointSystem(const char* functions, const char* subsets);
				AdjointSystem(const std::vector<std::string>& vFct, const std::vector<std::string>& vSubset);
				
			protected:
			///	Base class type, why do we have to add it?
			typedef IElemDisc<TDomain> base_type;
				
			///Also self-type is defined
			typedef AdjointSystem<TDomain> this_type;
			//PUBLIC MEMBERS:
			public:
				///	quadrature order
				bool m_bQuadOrderUserDef;
				int m_quadOrder;

				///	current shape function set, these are the identifiers for the Local Finite Elements used
				LFEID m_vLFEID;
				LFEID m_pLFEID;

				///	current element
				GridObject* m_pElem;
				
				///	World dimension
				static const int dim = base_type::dim;

				//stabilization types
				number m_stab_type;
			public:

				void set_quad_order(size_t order);
				///	type of trial space for each function used
				virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);
				///Functions for imports
				void set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > user);
				void set_kinematic_viscosity(number val);
				#ifdef UG_FOR_LUA
				void set_kinematic_viscosity(const char* fctName);
				void set_kinematic_viscosity(LuaFunctionHandle fct);
				#endif
				
				void set_stabilization(number stab)
				{
					m_stabParam=stab;
				}
				void set_stabilization_type(number stab)
				{
					m_stab_type=stab;
				}
			protected:
				///	register utils
			///	\{
				void register_all_funcs(const LFEID& vLfeID, const LFEID& pLfeID, const int quadOrder);
				template<typename TElem, typename VGeom, typename PGeom> void register_func();
			/// \}
			///	Data import for kinematic viscosity
				DataImport<number, dim> m_imKinViscosity;

				/// abbreviation for pressure
				static const size_t _P_ = dim;
				
				/// Stabilization parameter
				number m_stabParam;
				
			//PUBLIC ASSEMBLY FUNCTIONS AND PREP FUNCTIONS, these are copy/pasted from navier_stokes_fe.h
			public:
			
				template<typename TElem, typename VGeom, typename PGeom>
				void prep_elem_loop(const ReferenceObjectID roid, const int si);

				template<typename TElem, typename VGeom, typename PGeom>
				void prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[]);

			///	finishes the loop over all elements
				template<typename TElem, typename VGeom, typename PGeom>
				void fsh_elem_loop();
				
			///	assembles the local stiffness matrix using a finite volume scheme
				template<typename TElem, typename VGeom, typename PGeom>
				void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

			///	assembles the local mass matrix using a finite volume scheme
				template<typename TElem, typename VGeom, typename PGeom>
				void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

			///	assembles the stiffness part of the local defect
				template<typename TElem, typename VGeom, typename PGeom>
				void add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

			///	assembles the mass part of the local defect
				template<typename TElem, typename VGeom, typename PGeom>
				void add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

			///	assembles the local right hand side
				template<typename TElem, typename VGeom, typename PGeom>
				void add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[]); 
				
			protected:
				void IdentityDeformationGradient(MathMatrix<dim, dim>& GradU, const size_t ip);
				void VelocityGradient(MathMatrix<dim, dim>& GradV, const size_t ip);
				void VelocityVector(MathVector<dim>& V, const size_t ip);
				DataImport<MathVector<dim>,dim> m_imDeformationVectord1;//row1 Deformation Gradient
				DataImport<MathVector<dim>,dim> m_imDeformationVectord2;//row2 Deformation Gradient
				DataImport<MathVector<dim>,dim> m_imDeformationVectord3;//row3 Deformation Gradient
				
				DataImport<MathVector<dim>,dim> m_imVelocityGradientd1;//row1 Velocity Gradient
				DataImport<MathVector<dim>,dim> m_imVelocityGradientd2;//row2 Velocity Gradient
				DataImport<MathVector<dim>,dim> m_imVelocityGradientd3;//row3 Velocity Gradient
				
				DataImport<number, dim> m_imVelocityd1;//flow on x-direction
				DataImport<number, dim> m_imVelocityd2;//flow on y-direction
				DataImport<number, dim> m_imVelocityd3;//flow on z-direction
				
			public:
				void set_velocity_d1(SmartPtr<CplUserData<number, dim> > user);
				void set_velocity_d1(number val);
				void set_velocity_d2(SmartPtr<CplUserData<number, dim> > user);
				void set_velocity_d2(number val);
				void set_velocity_d3(SmartPtr<CplUserData<number, dim> > user);
				void set_velocity_d3(number val);
				
				void set_velocity_vector_d1(SmartPtr<CplUserData<MathVector<dim>,dim>> user);
				void set_velocity_vector_d1(number val);
				void set_velocity_vector_d2(SmartPtr<CplUserData<MathVector<dim>,dim>> user);
				void set_velocity_vector_d2(number val);
				void set_velocity_vector_d3(SmartPtr<CplUserData<MathVector<dim>,dim>> user);
				void set_velocity_vector_d3(number val);
				
				void set_deformation_vector_d1(SmartPtr<CplUserData<MathVector<dim>, dim> > user);
				void set_deformation_vector_d1(number val);
				#ifdef UG_FOR_LUA
				void set_deformation_vector_d1(const char* fctName);
				void set_deformation_vector_d1(LuaFunctionHandle fct);
				#endif
				//TODO:set default values, maybe in constructor?
				void set_deformation_vector_d2(SmartPtr<CplUserData<MathVector<dim>, dim> > user);
				void set_deformation_vector_d2(number val);
				#ifdef UG_FOR_LUA
				void set_deformation_vector_d2(const char* fctName);
				void set_deformation_vector_d2(LuaFunctionHandle fct);
				#endif
				void set_deformation_vector_d3(SmartPtr<CplUserData<MathVector<dim>, dim> > user);
				void set_deformation_vector_d3(number val);
				
			//For debugging purposes
			protected:
				bool m_bNoDeformation;//DeformationGradient is equal to Identity
			
				
		};//end AdjointSystem
	}//end namespace FluidOpt
}//end namespace ug
#endif /*ADJOINT_SYSTEM*/