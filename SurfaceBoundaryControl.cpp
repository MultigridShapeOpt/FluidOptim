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
 
#include "SurfaceBoundaryControl.h"
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
		void SurfaceBoundaryControl<TDomain>::
		set_boundary_control(SmartPtr<CplUserData<number,dim> > user_data)
		{
			this->m_imD.set_data(user_data);
		}
		
		template<typename TDomain>
		void SurfaceBoundaryControl<TDomain>::
		set_boundary_control(number control)
		{
			this->set_boundary_control(make_sp(new ConstUserNumber<dim>(control)));
		}
		
		template<typename TDomain>
		void SurfaceBoundaryControl<TDomain>::
		prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
		{
		//	check that not ADAPTIVE
			if(vLfeID[0].order() < 1)
				UG_THROW("SurfaceBoundaryControl: Adaptive order not implemented.");

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
			register_all_funcs(m_order, m_quadOrder);
	
		}
		
		////////////////////////////////////////////////////////////////////////////////
		//	assembling functions
		////////////////////////////////////////////////////////////////////////////////
		
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void SurfaceBoundaryControl<TDomain>::
		prep_elem_loop(const ReferenceObjectID roid, const int si)
		{
		//	register subsetIndex at Geometry
			TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

			try{
				geo.update_local(roid, m_lfeID, m_quadOrder);
			}
			UG_CATCH_THROW("SurfaceBoundaryControl::prep_elem_loop:"
								" Cannot update Finite Element Geometry.");

			
			static const int refDim = TElem::dim;
		//  set local positions
			//BoundaryControl
			m_imD.template set_local_ips<refDim>(geo.local_ips(),geo.num_ip(),false);
		}
		
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void SurfaceBoundaryControl<TDomain>::
		prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
		{
		//  update Geometry for this element
			TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);
			try{
				geo.update(elem, vCornerCoords, m_lfeID, m_quadOrder);

			}
			UG_CATCH_THROW("SurfaceBoundaryControl::prep_elem: "
								"Cannot update Finite Element Geometry.");
								
		//	set global positions
			//BoundaryControl
			m_imD.set_global_ips(geo.global_ips(),geo.num_ip());

		}
		
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void SurfaceBoundaryControl<TDomain>::
		fsh_elem_loop()
		{	}
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void SurfaceBoundaryControl<TDomain>::
		add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{	}//end jac func
		
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void SurfaceBoundaryControl<TDomain>::
		add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{
			TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);
			
			ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
			static const int refDim = TElem::dim;
			
			MathVector<dim> vNormal;VecSet(vNormal,0.0);
			MathVector<dim> vUnitNormal;VecSet(vUnitNormal,0.0);
			ElementNormal<dim>(id,vNormal,vCornerCoords);
			const number normalNorm = VecTwoNorm(vNormal);
			VecScale(vUnitNormal, vNormal, 1.0/normalNorm);
			
			//testing because we are suspicious people
			number ksi=0.05;//small number for tolerance
			
			if(VecTwoNorm(vUnitNormal) < (1-ksi) || VecTwoNorm(vUnitNormal) > (1+ksi))
			{
				UG_THROW("SurfaceBoundaryControl::add_rhs_elem: UnitNormal is not unit");
			}
				
			
			for (size_t ip = 0; ip < geo.num_ip(); ++ip){
				const number D = m_imD[ip];
				
				for(size_t sh = 0; sh < geo.num_sh(); ++sh){		
					for(int d1 = 0; d1 < dim; ++d1){
						
						d(d1,sh) += D * geo.shape(ip,sh) * vUnitNormal[d1]*geo.weight(ip);
					
					}
				}//end sh for
			}//end ip for
		}
		
		////////////////////////////////////////////////////////////////////////////////
		//	Constructor
		////////////////////////////////////////////////////////////////////////////////

		template<typename TDomain>
		SurfaceBoundaryControl<TDomain>::SurfaceBoundaryControl(const char* function, const char* subsets)
		: IElemDisc<TDomain>(function, subsets),m_order(1),
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
			
		}
		////////////////////////////////////////////////////////////////////////////////
		//	register assemble functions
		////////////////////////////////////////////////////////////////////////////////
		#ifdef UG_DIM_1
		template<>
		void SurfaceBoundaryControl<Domain1d>::register_all_funcs(int order, int quadOrder)
		{
			//	RegularEdge
			register_func<RegularEdge, DimFEGeometry<dim> > ();
		}
		#endif
		#ifdef UG_DIM_2
		template<>
		void SurfaceBoundaryControl<Domain2d>::register_all_funcs(int order, int quadOrder)
		{
			register_func<RegularEdge, DimFEGeometry<dim, 1> >();
			register_func<Triangle, DimFEGeometry<dim, 2> >();
			register_func<Quadrilateral, DimFEGeometry<dim, 2> >();
			if (quadOrder != 2 * order + 1) {
				register_func<Triangle, DimFEGeometry<dim> > ();
				register_func<Quadrilateral, DimFEGeometry<dim> > ();
			}

		}
		#endif
		#ifdef UG_DIM_3
		template<>
		void SurfaceBoundaryControl<Domain3d>::register_all_funcs(int order, int quadOrder)
		{
			register_func<Triangle, DimFEGeometry<dim, 2> >();
			register_func<Quadrilateral, DimFEGeometry<dim, 2> >();
			
			register_func<Tetrahedron, DimFEGeometry<dim, 3> >();
			register_func<Prism, DimFEGeometry<dim, 3> >();
			register_func<Pyramid, DimFEGeometry<dim, 3> >();
			register_func<Hexahedron, DimFEGeometry<dim, 3> >();
			register_func<Octahedron, DimFEGeometry<dim, 3> >();
		}
		#endif
		
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void SurfaceBoundaryControl<TDomain>::register_func()
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
		template class SurfaceBoundaryControl<Domain1d>;
		#endif
		#ifdef UG_DIM_2
		template class SurfaceBoundaryControl<Domain2d>;
		#endif
		#ifdef UG_DIM_3
		template class SurfaceBoundaryControl<Domain3d>;
		#endif
	 }//end namespace FluidOptim
 }//end namespace ug4