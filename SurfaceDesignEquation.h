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
 
 
#ifndef SURAFCE_DESIGN_EQUATION
#define SURAFCE_DESIGN_EQUATION

//These includes were taken from restricted_deformation_elasticity.h

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
		
		/*!\brief Deformation from the obstacle surface to the rest of volume.
		*
		*The model is composed of a linear term, nonlinear term (similar to convection), and
		*a boundary term. The vector \f$ \vec{w}\f$ is the variable
		*/
		
		template<typename TDomain>
		class SurfaceDesignEquation: public IElemDisc<TDomain>
		{
			protected:
				///	Base class type
				typedef IElemDisc<TDomain> base_type;
				
				///Also self-type is defined
				typedef SurfaceDesignEquation<TDomain> this_type;
			
			//PUBLIC MEMBERS:
			public:
				///	integration quadrature order
				bool m_bQuadOrderUserDef;
				int m_quadOrder;

				///	current shape function set, these are the identifiers for the Local Finite Elements used
				LFEID m_lfeID;
				
				///	current inner subset
				int m_si;
				
				///current order of disc scheme
				int m_order;

				///	current element
				GridObject* m_pElem;
				
				///	World dimension
				static const int dim = base_type::dim;
				
				///Domain type
				///Position type
			public:
				
				/*
					Functions for setting of the Variable Boundary Control through Lua Functions/Constant/
					After, functions for checking the set and value of the Import
				*/
				//Boundary Control D
				void set_boundary_control(SmartPtr<CplUserData<number, dim>> user_data);
				void set_boundary_control(number control);
				
				// Deformation values to form deformation vector
				void set_adjoint_deformation_d1(SmartPtr<CplUserData<number, dim> > user);
				void set_adjoint_deformation_d1(number val);
				void set_adjoint_deformation_d2(SmartPtr<CplUserData<number, dim> > user);
				void set_adjoint_deformation_d2(number val);
				void set_adjoint_deformation_d3(SmartPtr<CplUserData<number, dim> > user);
				void set_adjoint_deformation_d3(number val);

				
			//BORDER CONTROL VARIABLE
			protected:
				//Set vector L (dim sized)
				void AdjointDeformationVector(MathVector<dim>& L, const size_t ip);
				//Set vector W (dim sized)

				//Adjoint Deformation

				//border control
				DataImport<number,dim> m_imD;
				//Adjoint Deformation
				DataImport<number,dim> m_imAdjointDeformationd1;
				DataImport<number,dim> m_imAdjointDeformationd2;
				DataImport<number,dim> m_imAdjointDeformationd3;
				//regularization parameter
				number m_alpha;
								
				//TODO: this is used to directly access the components of extension factor dof K, below for boundary control gradient G
				static const size_t _G_ = 0;

			//CONSTRUCTOR
			public:
				///constructor
				SurfaceDesignEquation(const char* functions, const char* subsets);
				///destructor
				//virtual ~SurfaceDesignEquation();
				void set_alpha(number al){
					this->m_alpha=al;
				}
				//sets and gets and other virtual shit that has to be implemented
				/// sets the quad order
				void set_quad_order(const size_t order)
				{
					m_quadOrder=order;
					m_bQuadOrderUserDef=true;
				}
				int get_quad_order()
				{
					return m_quadOrder;
				}
				virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);
				
				void register_all_funcs(int order, int qorder);
				
				template<typename TElem, typename TFEGeom>
				void register_func();// here we register the implemented IElemDisc methods 
				
			
			//PUBLIC ASSEMBLY FUNCTIONS AND PREP FUNCTIONS, these are copy/pasted from navier_stokes_fe.h
			public:
			
				template<typename TElem, typename TFEGeom>
				void prep_elem_loop(const ReferenceObjectID roid, const int si);

				template<typename TElem, typename TFEGeom>
				void prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[]);

				template<typename TElem, typename TFEGeom>
				void fsh_elem_loop();
			
				template<typename TElem, typename TFEGeom>
				void add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[]); 
				
				template<typename TElem, typename TFEGeom>
				void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);
			
			/* 
				Useless add methods, only rhs is assembled here on the obstacle's surface
			*/
			protected:
				template<typename TElem, typename TFEGeom>
				void add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]){}
				template<typename TElem, typename TFEGeom>
				void add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]){}
				template<typename TElem, typename TFEGeom>
				void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]){}
				
				
		
		};//end SurfaceDesignEquation
		
	}//end namespace FluidOpt
}//end namespace ug
#endif /*DESIGN_EQUATION*/