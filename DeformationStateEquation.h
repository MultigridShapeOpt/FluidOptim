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

#ifndef DEFORMATION_STATE_EQUATION
#define DEFORMATION_STATE_EQUATION

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


template <typename TDomain>
class DeformationStateEquation
	: public IElemDisc<TDomain>
{
	private:
	///	Base class type
		typedef IElemDisc<TDomain> base_type;

	///	own type
		typedef DeformationStateEquation<TDomain> this_type;

	///	base element type of associated domain
		typedef typename domain_traits<TDomain::dim>::grid_base_object TBaseElem;

	public:
	///	Domain type
		typedef typename base_type::domain_type domain_type;

	///	World dimension
		static const int dim = base_type::dim;

	///	Position type
		typedef typename base_type::position_type position_type;

	public:
	/// constructor
		DeformationStateEquation(const char* functions, const char* subsets);

	///	destructor
		virtual	~DeformationStateEquation();


	///	sets the quad order
		void set_quad_order(const size_t order) {m_quadOrder = order; m_bQuadOrderUserDef = true;}
	///	gets the quad order
		int get_quad_order() {return m_quadOrder;}
		
		void set_nonlinear(bool nl);
		void set_picard(bool pic);
		void set_symmetric(bool sym);
		//void set_nu_extension(number n);
		// This function actually sets the value for the data through DataImport.set_data()
		void set_extension_factor(SmartPtr<CplUserData<number, dim>> user_data);
		// Functions below will create a SmartPtr with the provided parameter
		void set_extension_factor(number extension);
	///	type of trial space for each function used
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
		{
		//	check number
			if(vLfeID.size() != (size_t)dim)
				UG_THROW("SmallStrainMechanics: Needs exactly "<<dim<<" functions.");

			//	check & set order
			for(size_t i = 0; i < vLfeID.size(); ++i)
				if(vLfeID[i].order() < 1)
					UG_THROW("SmallStrainMechanics: Adaptive order not implemented.");

		//	remember lfeID;
			m_lfeID = vLfeID[0];

		//	set order
			m_order = vLfeID[0].order();

		//	update assemble functions
			set_assemble_funcs();
		}


	private:
		size_t num_fct() const {return dim;}

		///	sets the requested assembling routines
		void set_assemble_funcs();

		void register_all_fe_funcs(int order, int quadOrder);
		template <typename TElem, typename TFEGeom>
		void register_fe_func();

		///	assemble methods
		template<typename TElem, typename TFEGeom>
		void prep_timestep_elem(const number time, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		template<typename TElem, typename TFEGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);

		template<typename TElem, typename TFEGeom>
		void prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[]);

		template<typename TElem, typename TFEGeom>
		void fsh_elem_loop();

		template<typename TElem, typename TFEGeom>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		template<typename TElem, typename TFEGeom>
		void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		template<typename TElem, typename TFEGeom>
		void add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		template<typename TElem, typename TFEGeom>
		void add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		template<typename TElem, typename TFEGeom>
		void add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		template<typename TElem, typename TFEGeom>
		void fsh_timestep_elem(const number time, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);


	///	updates the geometry for a given element
		void update_geo_elem(TBaseElem* elem, DimFEGeometry<dim>& geo);


	protected:
		number m_variable;
		template<typename TFEGeom>
		void DisplacementGradient(MathMatrix<dim, dim>& GradU, const size_t ip,
									const TFEGeom& geo, const LocalVector& u)
		{
			//	loop shape-functions at one integration point ip in order
			//	to compute local_grad(ip,i): \frac{\partial N_i}{\eps_ip}
			//	and global_grad(ip,i): \frac{\partial N_i}{\X_ip}

			for (size_t i = 0; i < (size_t) dim; ++i)
				for (size_t J = 0; J < (size_t) dim; ++J)
				{
					GradU[i][J] = 0.0;

					//	compute GradU: displacementGradient
					for (size_t a = 0; a < geo.num_sh(); ++a)
						GradU[i][J] += geo.global_grad(ip, a)[J] * u(i, a);
				}
		}

	private:
	

	///	current order of disc scheme
		int m_order;

	///	current shape function set
		LFEID m_lfeID;

	///	current integration order
		bool m_bQuadOrderUserDef;
		int m_quadOrder;
		
		bool m_bNonLinear;
		bool m_bPicardIteration;
		bool m_bSymmetric;
		//number m_nu_extension;
	protected:
	//Imports for m_nu_extension
	DataImport<number,dim> m_imExtensionFactor;//K
/* 	
	public:
	void set_extension_factor(SmartPtr<CplUserData<number, dim> > user);
	void set_extension_factor(number val); */
		
};//end class

}
}

#endif /* DEFORMATION_STATE_EQUATION */
