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
 
#include "DesignEquation.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/disc_util/fe_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

 namespace ug{
	 namespace FluidOptim{
		
		//Calculates the velocities at the ip and puts them into a vector 
		template<typename TDomain>
		void DesignEquation<TDomain>::
		AdjointDeformationVector(MathVector<dim>& L, const size_t ip)
		{
			number w1= m_imAdjointDeformationd1[ip];
			number w2= m_imAdjointDeformationd2[ip];
			
			L[0]=w1;
			L[1]=w2;
			if(this->dim == 3){
				number w3 = m_imAdjointDeformationd3[ip];
				L[2] = w3;
			} 
		}

		//***************IMPORT: ADJOINT DEFORMATION VECTOR
		//For the d1
		template<typename TDomain>
		void DesignEquation<TDomain>::
		set_adjoint_deformation_d1(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imAdjointDeformationd1.set_data(user);
		}
		template<typename TDomain>
		void DesignEquation<TDomain>::
		set_adjoint_deformation_d1(number val)
		{
			set_adjoint_deformation_d1(make_sp(new ConstUserNumber<dim>(val)));
		}
				
		//For the d2 
		template<typename TDomain>
		void DesignEquation<TDomain>::
		set_adjoint_deformation_d2(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imAdjointDeformationd2.set_data(user);
		}
		template<typename TDomain>
		void DesignEquation<TDomain>::
		set_adjoint_deformation_d2(number val)
		{
			set_adjoint_deformation_d2(make_sp(new ConstUserNumber<dim>(val)));
		}
		
		//For the d3
		template<typename TDomain>
		void DesignEquation<TDomain>::
		set_adjoint_deformation_d3(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imAdjointDeformationd3.set_data(user);
		}
		template<typename TDomain>
		void DesignEquation<TDomain>::
		set_adjoint_deformation_d3(number val)
		{
			set_adjoint_deformation_d3(make_sp(new ConstUserNumber<dim>(val)));
		}
		
		////////////////////////////////////////////////////////////////////////////////
		//	Functions that need registration
		////////////////////////////////////////////////////////////////////////////////				
		template<typename TDomain>
		void DesignEquation<TDomain>::
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
		
		template<typename TDomain>
		void DesignEquation<TDomain>::
		set_boundary_control(SmartPtr<CplUserData<number,dim> > user_data)
		{
			this->m_imD.set_data(user_data);
		}
		
		template<typename TDomain>
		void DesignEquation<TDomain>::
		set_boundary_control(number control)
		{
			this->set_boundary_control(make_sp(new ConstUserNumber<dim>(control)));
		}
		
		template<typename TDomain>
		void DesignEquation<TDomain>::
		prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
		{
	/* 	//	check number
			if(vLfeID.size() != TDomain::dim)
				UG_THROW("DesignEquation: needs exactly dim function.");
 */
		//	check that not ADAPTIVE
			if(vLfeID[0].order() < 1)
				UG_THROW("DesignEquation: Adaptive order not implemented.");

		//	set order
			m_lfeID = vLfeID[0];
			m_order = vLfeID[0].order();
			//m_quadOrder = 2*m_order+1;
		//	set default quadrature order if not set by user
			if (!m_bQuadOrderUserDef) {
				m_quadOrder = 2 * m_order + 1;
			}
			//	set all non-set orders
			else
			{
				if (m_quadOrder < 0){
					m_quadOrder = 2 * m_order + 1;
				}
			}
			register_all_funcs(m_order);
	
		}
		
		////////////////////////////////////////////////////////////////////////////////
		//	assembling functions
		////////////////////////////////////////////////////////////////////////////////
		
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void DesignEquation<TDomain>::
		prep_elem_loop(const ReferenceObjectID roid, const int si)
		{
			this->m_si = si;

		//	register subsetIndex at Geometry
			TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

			try{
				geo.update_local(roid, m_lfeID, m_quadOrder);
			}
			UG_CATCH_THROW("DesignEquation::prep_elem_loop:"
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
			//BoundaryControl
			m_imD.template set_local_ips(geo.local_ips(),geo.num_ip(),false);

			//AdjointDeformation Vector
			m_imAdjointDeformationd1.template set_local_ips(geo.local_ips(),geo.num_ip(),false);
			m_imAdjointDeformationd2.template set_local_ips(geo.local_ips(),geo.num_ip(),false);
			m_imAdjointDeformationd3.template set_local_ips(geo.local_ips(),geo.num_ip(),false);

		}
		
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void DesignEquation<TDomain>::
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
			UG_CATCH_THROW("DesignEquation::prep_elem: "
								"Cannot update Finite Element Geometry.");
								
		//	set global positions
			//BoundaryControl
			m_imD.set_global_ips(geo.global_ips(),geo.num_ip());

			//AdjointDeformation vector
			m_imAdjointDeformationd1.set_global_ips(geo.global_ips(),geo.num_ip());
			m_imAdjointDeformationd2.set_global_ips(geo.global_ips(),geo.num_ip());
			m_imAdjointDeformationd3.set_global_ips(geo.global_ips(),geo.num_ip());

		}
		
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void DesignEquation<TDomain>::
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
		void DesignEquation<TDomain>::
		add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{
			TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);
			typedef typename TFEGeom::BF BF;
			
			size_t fcts= this-> num_fct();
				
			for(size_t s = 0; s < m_vBndSSGrp.size(); ++s){
					
				const int si = m_vBndSSGrp[s];
				const std::vector<BF>& vBF = geo.bf(si);
				//Mass matrix-like
				for(size_t b = 0; b < vBF.size(); ++b){
					for(size_t ip = 0; ip < vBF[b].num_ip(); ++ip){
						for(size_t sh1 = 0; sh1 < vBF[b].num_sh(); ++sh1){
							//for(size_t index = 0; index < fcts; ++index){
								for(size_t sh2 = 0; sh2 < vBF[b].num_sh(); ++sh2){
									J(_G_,sh1,_G_,sh2) += vBF[b].shape(ip,sh1)
															 *vBF[b].shape(ip,sh2) 
															 *vBF[b].weight(ip);
								}
							//}
						}
					}					
				}				
			}			
		}//end jac func
		
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void DesignEquation<TDomain>::
		add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{
			TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);
			typedef typename TFEGeom::BF BF;
			
			size_t fcts= this-> num_fct();
				
			for(size_t s = 0; s < m_vBndSSGrp.size(); ++s){
					
				const int si = m_vBndSSGrp[s];
				const std::vector<BF>& vBF = geo.bf(si);

				for(size_t b = 0; b < vBF.size(); ++b){
					for(size_t ip = 0; ip < vBF[b].num_ip(); ++ip){
						const number D = m_imD[ip];
						
						for(size_t sh = 0; sh < vBF[b].num_sh(); ++sh){
							//for(size_t index = 0; index < fcts; ++index){
								MathVector<dim> vAdjDef;VecSet(vAdjDef,0.0);
								AdjointDeformationVector(vAdjDef,ip);
								MathVector<dim> normal = vBF[b].normal();
								const number vol = VecTwoNorm(normal);
								VecScale(normal, normal, 1.0/vol);
								
								d(_G_,sh) -= m_alpha * D * vBF[b].shape(ip,sh) * vBF[b].weight(ip) 
											   + VecDot(normal,vAdjDef) * vBF[b].shape(ip,sh) * vBF[b].weight(ip);
												
							//}

						}
					}
				}
			}
		}
		
		////////////////////////////////////////////////////////////////////////////////
		//	Constructor
		////////////////////////////////////////////////////////////////////////////////

		template<typename TDomain>
		DesignEquation<TDomain>::DesignEquation(const char* function)
		: IElemDisc<TDomain>(function, ""),m_order(1),
		m_lfeID(LFEID::LAGRANGE, TDomain::dim, m_order)
		{
			this->clear_add_fct();
			
			//	set defaults
			m_order = 1;
			m_bQuadOrderUserDef = false;
			m_quadOrder = -1;
			
			//register imports
			//BoundaryControl
			this->register_import(m_imD);

			//AdjointDeformationVector
			this->register_import(m_imAdjointDeformationd1);
			this->register_import(m_imAdjointDeformationd2);
			this->register_import(m_imAdjointDeformationd3);


			
		}
		////////////////////////////////////////////////////////////////////////////////
		//	register assemble functions
		////////////////////////////////////////////////////////////////////////////////
		#ifdef UG_DIM_1
		template<>
		void DesignEquation<Domain1d>::register_all_funcs(int order)
		{
			//	RegularEdge
			register_func<RegularEdge, DimFEGeometry<dim> > ();
		}
		#endif
		#ifdef UG_DIM_2
		template<>
		void DesignEquation<Domain2d>::register_all_funcs(int order)
		{
			
			register_func<Triangle, DimFEGeometry<dim, 2> >();
			register_func<Quadrilateral, DimFEGeometry<dim, 2> >();
		}
		#endif
		#ifdef UG_DIM_3
		template<>
		void DesignEquation<Domain3d>::register_all_funcs(int order)
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
		void DesignEquation<TDomain>::register_func()
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
		template class DesignEquation<Domain1d>;
		#endif
		#ifdef UG_DIM_2
		template class DesignEquation<Domain2d>;
		#endif
		#ifdef UG_DIM_3
		template class DesignEquation<Domain3d>;
		#endif
	 }//end namespace FluidOptim
 }//end namespace ug4