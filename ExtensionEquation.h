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
 
#ifndef EXTENSION_EQUATION
#define EXTENSION_EQUATION

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
		
		template<typename TDomain>
		class ExtensionEquation: public IElemDisc<TDomain>
		{
			protected:
				///	Base class type
				typedef IElemDisc<TDomain> base_type;
				
				///Also self-type is defined
				typedef ExtensionEquation<TDomain> this_type;
			
			//PUBLIC MEMBERS:
			public:
				///	integration quadrature order
				bool m_bQuadOrderUserDef;
				int m_quadOrder;

				///	current shape function set, these are the identifiers for the Local Finite Elements used
				LFEID m_lfeID;
				
				
				///current order of disc scheme
				int m_order;

				///	current element
				GridObject* m_pElem;
				
				///	World dimension
				static const int dim = base_type::dim;
				
				///Domain type
				///Position type
			public:
		
				// Deformation values to form deformation vector
				void set_adjoint_deformation_d1(SmartPtr<CplUserData<number, dim> > user);
				void set_adjoint_deformation_d1(number val);
				void set_adjoint_deformation_d2(SmartPtr<CplUserData<number, dim> > user);
				void set_adjoint_deformation_d2(number val);
				void set_adjoint_deformation_d3(SmartPtr<CplUserData<number, dim> > user);
				void set_adjoint_deformation_d3(number val);
				// Deformation values to form deformation vector
				void set_deformation_d1(SmartPtr<CplUserData<number, dim> > user);
				void set_deformation_d1(number val);
				void set_deformation_d2(SmartPtr<CplUserData<number, dim> > user);
				void set_deformation_d2(number val);
				void set_deformation_d3(SmartPtr<CplUserData<number, dim> > user);
				void set_deformation_d3(number val);
				// Deformation vectors for DeformationGradient and Deformation 
				void set_deformation_vector_d1(SmartPtr<CplUserData<MathVector<dim>, dim> > user);
				void set_deformation_vector_d1(number val);
				void set_deformation_vector_d2(SmartPtr<CplUserData<MathVector<dim>, dim> > user);
				void set_deformation_vector_d2(number val);
				void set_deformation_vector_d3(SmartPtr<CplUserData<MathVector<dim>, dim> > user);
				void set_deformation_vector_d3(number val);
				
				//Extension factor
				void set_extension_factor(SmartPtr<CplUserData<number, dim> > user);
				void set_extension_factor(number val);
			//BORDER CONTROL VARIABLE
			protected:
				//Set adjoint deformation vector L (dim sized)
				void AdjointDeformationVector(MathVector<dim>& L, const size_t ip);
				//Set deformation vector W (dim sized)
				void DeformationVector(MathVector<dim>& W, const size_t ip);
				//Set deformation gradient Grad_W
				void DeformationGradient(MathMatrix<dim, dim>& GradW, const size_t ip);

				//Adjoint Deformation
				DataImport<number,dim> m_imAdjointDeformationd1;
				DataImport<number,dim> m_imAdjointDeformationd2;
				DataImport<number,dim> m_imAdjointDeformationd3;
				
				//Deformation
				DataImport<number,dim> m_imDeformationd1;
				DataImport<number,dim> m_imDeformationd2;
				DataImport<number,dim> m_imDeformationd3;
				
				//Deformation Gradient
				DataImport<MathVector<dim>,dim> m_imDeformationVectord1;//row1 Deformation Gradient
				DataImport<MathVector<dim>,dim> m_imDeformationVectord2;//row2 Deformation Gradient
				DataImport<MathVector<dim>,dim> m_imDeformationVectord3;//row3 Deformation Gradient
				//ExtensionFactor
				DataImport<number, dim> m_imNuExtensionValue;
				
				number m_nu_ext_upper;
				number m_nu_ext_lower;
				number m_chi;
				
				
				//TODO: this is used to directly access the components of extension factor dof K, below for boundary control gradient G

				static const size_t _K_ = 0;
			//CONSTRUCTOR
			public:
				///constructor
				ExtensionEquation(const char* functions, const char* subsets);

				void set_extension_upper(number up){
					this->m_nu_ext_upper=up;
				}
				void set_extension_lower(number low){
					this->m_nu_ext_lower=low;
				}
				void set_chi(number val){
					this->m_chi = val;
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
				
				
		
		};//end DesignEquation
		
	}//end namespace FluidOpt
}//end namespace ug
#endif /*EXTENSION_EQUATION*/