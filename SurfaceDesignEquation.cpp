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
 
#include "SurfaceDesignEquation.h"
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
		void SurfaceDesignEquation<TDomain>::
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
		void SurfaceDesignEquation<TDomain>::
		set_adjoint_deformation_d1(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imAdjointDeformationd1.set_data(user);
		}
		template<typename TDomain>
		void SurfaceDesignEquation<TDomain>::
		set_adjoint_deformation_d1(number val)
		{
			set_adjoint_deformation_d1(make_sp(new ConstUserNumber<dim>(val)));
		}
				
		//For the d2 
		template<typename TDomain>
		void SurfaceDesignEquation<TDomain>::
		set_adjoint_deformation_d2(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imAdjointDeformationd2.set_data(user);
		}
		template<typename TDomain>
		void SurfaceDesignEquation<TDomain>::
		set_adjoint_deformation_d2(number val)
		{
			set_adjoint_deformation_d2(make_sp(new ConstUserNumber<dim>(val)));
		}
		
		//For the d3
		template<typename TDomain>
		void SurfaceDesignEquation<TDomain>::
		set_adjoint_deformation_d3(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imAdjointDeformationd3.set_data(user);
		}
		template<typename TDomain>
		void SurfaceDesignEquation<TDomain>::
		set_adjoint_deformation_d3(number val)
		{
			set_adjoint_deformation_d3(make_sp(new ConstUserNumber<dim>(val)));
		}
		
		////////////////////////////////////////////////////////////////////////////////
		//	Functions that need registration
		////////////////////////////////////////////////////////////////////////////////				

		
		template<typename TDomain>
		void SurfaceDesignEquation<TDomain>::
		set_boundary_control(SmartPtr<CplUserData<number,dim> > user_data)
		{
			this->m_imD.set_data(user_data);
		}
		
		template<typename TDomain>
		void SurfaceDesignEquation<TDomain>::
		set_boundary_control(number control)
		{
			this->set_boundary_control(make_sp(new ConstUserNumber<dim>(control)));
		}
		
		template<typename TDomain>
		void SurfaceDesignEquation<TDomain>::
		prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
		{
	/* 	//	check number
			if(vLfeID.size() != TDomain::dim)
				UG_THROW("SurfaceDesignEquation: needs exactly dim function.");
 */
		//	check that not ADAPTIVE
			if(vLfeID[0].order() < 1)
				UG_THROW("SurfaceDesignEquation: Adaptive order not implemented.");

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
		void SurfaceDesignEquation<TDomain>::
		prep_elem_loop(const ReferenceObjectID roid, const int si)
		{
			this->m_si = si;

		//	register subsetIndex at Geometry
			TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

			try{
				geo.update_local(roid, m_lfeID, m_quadOrder);
			}
			UG_CATCH_THROW("SurfaceDesignEquation::prep_elem_loop:"
								" Cannot update Finite Element Geometry.");


			
			static const int refDim = TElem::dim;
		//  set local positions
			//BoundaryControl
			m_imD.template set_local_ips<refDim>(geo.local_ips(),geo.num_ip(),false);

			//AdjointDeformation Vector
			m_imAdjointDeformationd1.template set_local_ips<refDim>(geo.local_ips(),geo.num_ip(),false);
			m_imAdjointDeformationd2.template set_local_ips<refDim>(geo.local_ips(),geo.num_ip(),false);
			m_imAdjointDeformationd3.template set_local_ips<refDim>(geo.local_ips(),geo.num_ip(),false);

		}
		
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void SurfaceDesignEquation<TDomain>::
		prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
		{
		//  update Geometry for this element
			TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);
			try{
				geo.update(elem, vCornerCoords, m_lfeID, m_quadOrder);

			}
			UG_CATCH_THROW("SurfaceDesignEquation::prep_elem: "
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
		void SurfaceDesignEquation<TDomain>::
		fsh_elem_loop()
		{	}
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void SurfaceDesignEquation<TDomain>::
		add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{
			TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);
							
			//Mass-Matrix like for _K_
			for (size_t ip = 0; ip < geo.num_ip(); ++ip){			
				////////////////////////////////////////////////////
				// Stiffness/Diffusion Matrix Terms
				////////////////////////////////////////////////////
				for(size_t sh1=0;sh1<geo.num_sh();++sh1){
					for(size_t sh2=0;sh2<geo.num_sh();++sh2){
						J(_G_,sh1,_G_,sh2) += geo.shape(ip,sh1)
											*geo.shape(ip,sh2) 
											*geo.weight(ip);
					}//end sh2
				}//end sh1
			}//end ip 			
		}//end jac func
		
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void SurfaceDesignEquation<TDomain>::
		add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{	}
		
		////////////////////////////////////////////////////////////////////////////////
		//	Constructor
		////////////////////////////////////////////////////////////////////////////////

		template<typename TDomain>
		SurfaceDesignEquation<TDomain>::SurfaceDesignEquation(const char* function, const char* subsets)
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
		void SurfaceDesignEquation<Domain1d>::register_all_funcs(int order, int quadOrder)
		{
			//	RegularEdge
			register_func<RegularEdge, DimFEGeometry<dim> > ();
		}
		#endif
		#ifdef UG_DIM_2
		template<>
		void SurfaceDesignEquation<Domain2d>::register_all_funcs(int order, int quadOrder)
		{
			/* register_func<RegularEdge, DimFEGeometry<dim, order> > ();
			register_func<RegularEdge, DimFEGeometry<dim, 2> > ();
			register_func<Triangle, DimFEGeometry<dim, 2> >();
			register_func<Quadrilateral, DimFEGeometry<dim, 2> >(); */
			register_func<RegularEdge, DimFEGeometry<dim, 1> >();
			
			if (quadOrder != 2 * order + 1) {
				register_func<Triangle, DimFEGeometry<dim> > ();
				register_func<Quadrilateral, DimFEGeometry<dim> > ();
			}

			//	special compiled cases

			//	Triangle
			switch (order)
			{
				case 1:{
					typedef FEGeometry<Triangle, dim, LagrangeLSFS<ReferenceTriangle, 1> ,
							GaussQuadrature<ReferenceTriangle, 3> > FEGeom;
					register_func<Triangle, FEGeom > (); break;
				}
				case 2:{
					typedef FEGeometry<Triangle, dim, LagrangeLSFS<ReferenceTriangle, 2> ,
							GaussQuadrature<ReferenceTriangle, 5> > FEGeom;
					register_func<Triangle, FEGeom > (); break;
				}
				case 3:{
					typedef FEGeometry<Triangle, dim, LagrangeLSFS<ReferenceTriangle, 3> ,
							GaussQuadrature<ReferenceTriangle, 7> > FEGeom;
					register_func<Triangle, FEGeom > (); break;
				}
				default: register_func<Triangle, DimFEGeometry<dim> > (); break;
			}
		}
		#endif
		#ifdef UG_DIM_3
		template<>
		void SurfaceDesignEquation<Domain3d>::register_all_funcs(int order, int quadOrder)
		{
			register_func<Triangle, DimFEGeometry<dim, 2> > ();
			register_func<Quadrilateral, DimFEGeometry<dim, 2> > ();
			
			register_func<Tetrahedron, DimFEGeometry<dim, 3> >();
			register_func<Prism, DimFEGeometry<dim, 3> >();
			register_func<Pyramid, DimFEGeometry<dim, 3> >();
			register_func<Hexahedron, DimFEGeometry<dim, 3> >();
			register_func<Octahedron, DimFEGeometry<dim, 3> >();
		}
		#endif
		
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void SurfaceDesignEquation<TDomain>::register_func()
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
		template class SurfaceDesignEquation<Domain1d>;
		#endif
		#ifdef UG_DIM_2
		template class SurfaceDesignEquation<Domain2d>;
		#endif
		#ifdef UG_DIM_3
		template class SurfaceDesignEquation<Domain3d>;
		#endif
	 }//end namespace FluidOptim
 }//end namespace ug4