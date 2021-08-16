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
 
#include "BoundaryControl.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/disc_util/fe_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

 namespace ug{
	 namespace FluidOptim{
		 		 
		////////////////////////////////////////////////////////////////////////////////
		//	Functions that need registration
		////////////////////////////////////////////////////////////////////////////////
		
		template<typename TDomain>
		void BoundaryControl<TDomain>::
		set_variable_boundary_control(bool control)
		{
			this->m_bVariableBoundaryControl=control;
		}
		
		template<typename TDomain>
		void BoundaryControl<TDomain>::
		set_fixed_boundary_control(number n, const char* BndSubsets, const char* InnerSubsets)
		{
			this->m_boundary_control=n;
			this->set_boundary_inner_subsets(BndSubsets,InnerSubsets);
			/*
			// remember bnd name
			m_BndSubsetNames = std::string(BndSubsets);

			// schedule for assembling on all inner subsets
			std::vector<std::string> vSubsets = this->symb_subsets();
			std::vector<std::string> vNew = TokenizeTrimString(InnerSubsets);
			for(size_t i = 0; i < vNew.size(); ++i)
				if(std::find(vSubsets.begin(), vSubsets.end(), vNew[i]) == vSubsets.end())
					vSubsets.push_back(vNew[i]);
			this->set_subsets(vSubsets); 
			*/
		}
		template<typename TDomain>
		void BoundaryControl<TDomain>::
		set_boundary_inner_subsets(const char* BndSubsets, const char* InnerSubsets)
		{
			m_BndSubsetNames = std::string(BndSubsets);

			// schedule for assembling on all inner subsets
			std::vector<std::string> vSubsets = this->symb_subsets();
			std::vector<std::string> vNew = TokenizeTrimString(InnerSubsets);
			for(size_t i = 0; i < vNew.size(); ++i)
				if(std::find(vSubsets.begin(), vSubsets.end(), vNew[i]) == vSubsets.end())
					vSubsets.push_back(vNew[i]);
			this->set_subsets(vSubsets);
		}
		
		// /Boundary Control with LUA
		template<typename TDomain>
		void BoundaryControl<TDomain>::
		set_boundary_control(SmartPtr<CplUserData<number,dim> > user_data)
		{
			this->m_imBoundaryControl.set_data(user_data);
		}
		
		template<typename TDomain>
		void BoundaryControl<TDomain>::
		set_boundary_control(number control)
		{
			this->set_boundary_control(make_sp(new ConstUserNumber<dim>(control)));
		}
		
		#ifdef UG_FOR_LUA
		template<typename TDomain>
		void BoundaryControl<TDomain>::
		set_boundary_control(const char* fctName)
		{
			set_boundary_control(LuaUserDataFactory<number,dim>::create(fctName));
		}
		template<typename TDomain>
		void BoundaryControl<TDomain>::
		set_boundary_control(LuaFunctionHandle fct)
		{
			set_boundary_control(make_sp(new LuaUserData<number,dim>(fct)));
		}
		#endif
		template<typename TDomain>
		bool BoundaryControl<TDomain>::
		check_boundary_control()
		{
			return m_imBoundaryControl.data_given();
		}
		/* template<typename TDomain>
		bool BoundaryControl<TDomain>::
		check_value_boundary_control(number val)
		{
			if(*(this->m_imBoundaryControl.user_data())==val)
				return true;
			else
				return false;		
		} */
		
		template<typename TDomain>
		void BoundaryControl<TDomain>::
		prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
		{
		//	check number
			if(vLfeID.size() != TDomain::dim)
				UG_THROW("BoundaryControl: needs exactly dim function.");

		//	check that not ADAPTIVE
			if(vLfeID[0].order() < 1)
				UG_THROW("BoundaryControl: Adaptive order not implemented.");

		//	set order
			m_lfeID = vLfeID[0];
			m_order = vLfeID[0].order();
			m_quadOrder = 2*m_order+1;

			register_all_funcs(m_order);
	
		}
		
		////////////////////////////////////////////////////////////////////////////////
		//	assembling functions
		////////////////////////////////////////////////////////////////////////////////
		
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void BoundaryControl<TDomain>::
		prep_elem_loop(const ReferenceObjectID roid, const int si)
		{
			this->m_si = si;

		//	register subsetIndex at Geometry
			TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

			try{
				geo.update_local(roid, m_lfeID, m_quadOrder);
			}
			UG_CATCH_THROW("ShapeDiffVolumePerimeter::prep_elem_loop:"
								" Cannot update Finite Element Geometry.");

		//	get subset ids for boundary
			m_vBndSSGrp = this->approx_space()->subset_grp_by_name(m_BndSubsetNames.c_str());

		//	request subset indices as boundary subset. This will force the
		//	creation of boundary subsets when calling geo.upd ate
			for(size_t s = 0; s < m_vBndSSGrp.size(); ++s){
				const int si = m_vBndSSGrp[s];
				geo.add_boundary_subset(si);
			}
			
			static const int refDim = TElem::dim;
		//  set local positions
			m_imBoundaryControl.template set_local_ips(geo.local_ips(),geo.num_ip(),false);
		}
		
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void BoundaryControl<TDomain>::
		prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
		{
		//  update Geometry for this element
			TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);
			try{
				geo.update(elem, vCornerCoords, m_lfeID, m_quadOrder);
				geo.update_boundary_faces(elem, vCornerCoords,
						   m_quadOrder,
						   &(this->subset_handler()));
			}
			UG_CATCH_THROW("BoundaryControl::prep_elem: "
								"Cannot update Finite Element Geometry.");
								
		//	set global positions
			m_imBoundaryControl.set_global_ips(geo.global_ips(),geo.num_ip());
		}
		
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void BoundaryControl<TDomain>::
		fsh_elem_loop()
		{
		//	remove subsetIndex from Geometry
			TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID,m_order);

		//	unrequest subset indices as boundary subset. This will force the
		//	creation of boundary subsets when calling geo.update
			for(size_t s = 0; s < m_vBndSSGrp.size(); ++s){
				const int si = m_vBndSSGrp[s];
				geo.remove_boundary_subset(si);
			}
		}
		
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void BoundaryControl<TDomain>::
		add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{
			TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);
			typedef typename TFEGeom::BF BF;
			
			if(m_boundary_control != 0.0 || m_bVariableBoundaryControl){
				
				for(size_t s = 0; s < m_vBndSSGrp.size(); ++s){
					
					const int si = m_vBndSSGrp[s];
					const std::vector<BF>& vBF = geo.bf(si);

					for(size_t b = 0; b < vBF.size(); ++b){
						for(size_t ip = 0; ip < vBF[b].num_ip(); ++ip){
							const number D = m_imBoundaryControl[ip];
							for(size_t sh = 0; sh < vBF[b].num_sh(); ++sh){
								for(int d1 = 0; d1 < dim; ++d1){
									
									MathVector<dim> normal = vBF[b].normal();
									const number vol = VecTwoNorm(normal);
									VecScale(normal, normal, 1.0/vol);
									if(m_bVariableBoundaryControl)
									{
										d(d1,sh) += D * vBF[b].shape(ip,sh)* normal[d1]*vBF[b].weight(ip);
									}
									if(!m_bVariableBoundaryControl)
									{
										d(d1,sh) += m_boundary_control * vBF[b].shape(ip,sh)* normal[d1]*vBF[b].weight(ip);
									}

									//const MathVector<dim>& grad = vBF[b].global_grad(ip, sh);
									
								}

							}
						}
					}
				}
			}
		}
		////////////////////////////////////////////////////////////////////////////////
		//	Constructor
		////////////////////////////////////////////////////////////////////////////////

		template<typename TDomain>
		BoundaryControl<TDomain>::BoundaryControl(const char* function)
		: IElemDisc<TDomain>(function, ""),m_order(1),
		m_lfeID(LFEID::LAGRANGE, TDomain::dim, m_order)
		{
			this->clear_add_fct();
			this->m_bVariableBoundaryControl=false;
			m_boundary_control = 0.0;//no rhs
			
			//register imports
			this->register_import(m_imBoundaryControl);
			
			m_imBoundaryControl.set_rhs_part();
			
		}
		////////////////////////////////////////////////////////////////////////////////
		//	register assemble functions
		////////////////////////////////////////////////////////////////////////////////
		#ifdef UG_DIM_1
		template<>
		void BoundaryControl<Domain1d>::register_all_funcs(int order)
		{
			//	RegularEdge
			register_func<RegularEdge, DimFEGeometry<dim> > ();
		}
		#endif
		#ifdef UG_DIM_2
		template<>
		void BoundaryControl<Domain2d>::register_all_funcs(int order)
		{
			register_func<Triangle, DimFEGeometry<dim, 2> >();
			register_func<Quadrilateral, DimFEGeometry<dim, 2> >();
		}
		#endif
		#ifdef UG_DIM_3
		template<>
		void BoundaryControl<Domain3d>::register_all_funcs(int order)
		{
			register_func<Tetrahedron, DimFEGeometry<dim, 3> >();
			register_func<Prism, DimFEGeometry<dim, 3> >();
			register_func<Pyramid, DimFEGeometry<dim, 3> >();
			register_func<Hexahedron, DimFEGeometry<dim, 3> >();
			register_func<Octahedron, DimFEGeometry<dim, 3> >();
		}
		#endif
		
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void BoundaryControl<TDomain>::register_func()
		{
			ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
			typedef this_type T;

			this->clear_add_fct(id);

			this->set_prep_elem_loop_fct(id, &T::template prep_elem_loop<TElem, TFEGeom>);
			this->set_prep_elem_fct(	 id, &T::template prep_elem<TElem, TFEGeom>);
			this->set_add_rhs_elem_fct(	 id, &T::template add_rhs_elem<TElem, TFEGeom>);
			this->set_fsh_elem_loop_fct( id, &T::template fsh_elem_loop<TElem, TFEGeom>);

			this->set_add_jac_A_elem_fct(	 id, &T::template add_jac_A_elem<TElem, TFEGeom>);
			this->set_add_jac_M_elem_fct(	 id, &T::template add_jac_M_elem<TElem, TFEGeom>);
			this->set_add_def_A_elem_fct(	 id, &T::template add_def_A_elem<TElem, TFEGeom>);
			this->set_add_def_M_elem_fct(	 id, &T::template add_def_M_elem<TElem, TFEGeom>);
		}
		
		
		////////////////////////////////////////////////////////////////////////////////
		//	explicit template instantiations
		////////////////////////////////////////////////////////////////////////////////
		#ifdef UG_DIM_1
		template class BoundaryControl<Domain1d>;
		#endif
		#ifdef UG_DIM_2
		template class BoundaryControl<Domain2d>;
		#endif
		#ifdef UG_DIM_3
		template class BoundaryControl<Domain3d>;
		#endif

	 }//end namespace FluidOptim
 }//end namespace ug4