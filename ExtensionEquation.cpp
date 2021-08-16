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
 
#include "ExtensionEquation.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/disc_util/fe_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

 namespace ug{
	 namespace FluidOptim{
		////////////////////////////////////////////////////////////////////////////////
		//	PRIVATE USE FUNCTIONS
		////////////////////////////////////////////////////////////////////////////////
		//Calculates the velocities at the ip and puts them into a vector 
		template<typename TDomain>
		void ExtensionEquation<TDomain>::
		AdjointDeformationVector(MathVector<dim>& L, const size_t ip)
		{
			number w1= m_imAdjointDeformationd1[ip];
			number w2= m_imAdjointDeformationd2[ip];
			
			L[0]=w1;
			L[1]=w2;
			//TODO:better way to do this...cleaner...but this should work
			if(this->dim == 3){
				number w3= m_imAdjointDeformationd3[ip];
				L[2]=w3;
			}
		}
		//Calculation of the Identity Deformation Gradient
		template<typename TDomain>
		void ExtensionEquation<TDomain>::
		DeformationGradient(MathMatrix<dim, dim>& GradW, const size_t ip)
		{
			//1.- calculate gradients of each deformation component u1,u2
			MathVector<dim> vdeformation_d1=m_imDeformationVectord1[ip];
			MathVector<dim> vdeformation_d2=m_imDeformationVectord2[ip];
			//2.- assign the gradients to each row using assign(vector, row)
			GradW.assign(vdeformation_d1,0);
			GradW.assign(vdeformation_d2,1);
			
			if(this->dim == 3){
				MathVector<dim> vdeformation_d3=m_imDeformationVectord3[ip];
				GradW.assign(vdeformation_d3,2);
			}
			
		}
		//Calculates the velocities at the ip and puts them into a vector 
		template<typename TDomain>
		void ExtensionEquation<TDomain>::
		DeformationVector(MathVector<dim>& w, const size_t ip)
		{
			number w1= m_imDeformationd1[ip];
			number w2= m_imDeformationd2[ip];
			
			w[0]=w1;
			w[1]=w2;
			
			if(this->dim == 3){
				number w3= m_imDeformationd3[ip];
				w[2]=w3;
			}
		}
		////////////////////////////////////////////////////////////////////////////////
		//	IMPORT SETTING FUNCTIONS
		////////////////////////////////////////////////////////////////////////////////
		//***************IMPORT: DEFORMATION GRADIENT ROWS
		//For the d1
		template<typename TDomain>
		void ExtensionEquation<TDomain>::
		set_adjoint_deformation_d1(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imAdjointDeformationd1.set_data(user);
		}
		template<typename TDomain>
		void ExtensionEquation<TDomain>::
		set_adjoint_deformation_d1(number val)
		{
			set_adjoint_deformation_d1(make_sp(new ConstUserNumber<dim>(val)));
		}
				
		//For the d2 
		template<typename TDomain>
		void ExtensionEquation<TDomain>::
		set_adjoint_deformation_d2(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imAdjointDeformationd2.set_data(user);
		}
		template<typename TDomain>
		void ExtensionEquation<TDomain>::
		set_adjoint_deformation_d2(number val)
		{
			set_adjoint_deformation_d2(make_sp(new ConstUserNumber<dim>(val)));
		}
		//For the d3 
		template<typename TDomain>
		void ExtensionEquation<TDomain>::
		set_adjoint_deformation_d3(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imAdjointDeformationd3.set_data(user);
		}
		template<typename TDomain>
		void ExtensionEquation<TDomain>::
		set_adjoint_deformation_d3(number val)
		{
			set_adjoint_deformation_d3(make_sp(new ConstUserNumber<dim>(val)));
		}
		
		//***************IMPORT: DEFORMATION VECTOR COMPONENTS
		//d1
		template<typename TDomain>
		void ExtensionEquation<TDomain>::
		set_deformation_d1(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imDeformationd1.set_data(user);
		}
		template<typename TDomain>
		void ExtensionEquation<TDomain>::
		set_deformation_d1(number val)
		{
			set_deformation_d1(make_sp(new ConstUserNumber<dim>(val)));
		}		
		//d2
		template<typename TDomain>
		void ExtensionEquation<TDomain>::
		set_deformation_d2(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imDeformationd2.set_data(user);
		}
		template<typename TDomain>
		void ExtensionEquation<TDomain>::
		set_deformation_d2(number val)
		{
			set_deformation_d2(make_sp(new ConstUserNumber<dim>(val)));
		}
		//d3
		template<typename TDomain>
		void ExtensionEquation<TDomain>::
		set_deformation_d3(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imDeformationd3.set_data(user);
		}
		template<typename TDomain>
		void ExtensionEquation<TDomain>::
		set_deformation_d3(number val)
		{
			set_deformation_d3(make_sp(new ConstUserNumber<dim>(val)));
		}
		//***************IMPORT: DEFORMATION GRADIENT ROWS, until 3d
		//For the d1
		template<typename TDomain>
		void ExtensionEquation<TDomain>::
		set_deformation_vector_d1(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imDeformationVectord1.set_data(user);
		}
		template<typename TDomain>
		void ExtensionEquation<TDomain>::
		set_deformation_vector_d1(number val)
		{
			set_deformation_vector_d1(make_sp(new ConstUserVector<dim>(val)));
		}
		
		//For the d2 
		template<typename TDomain>
		void ExtensionEquation<TDomain>::
		set_deformation_vector_d2(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imDeformationVectord2.set_data(user);
		}
		template<typename TDomain>
		void ExtensionEquation<TDomain>::
		set_deformation_vector_d2(number val)
		{
			set_deformation_vector_d2(make_sp(new ConstUserVector<dim>(val)));
		}
		
		//For the d3 
		template<typename TDomain>
		void ExtensionEquation<TDomain>::
		set_deformation_vector_d3(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imDeformationVectord3.set_data(user);
		}
		template<typename TDomain>
		void ExtensionEquation<TDomain>::
		set_deformation_vector_d3(number val)
		{
			set_deformation_vector_d3(make_sp(new ConstUserVector<dim>(val)));
		}
		//***************IMPORT: Extension factor
		template<typename TDomain>
		void ExtensionEquation<TDomain>::
		set_extension_factor(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imNuExtensionValue.set_data(user);
		}
		template<typename TDomain>
		void ExtensionEquation<TDomain>::
		set_extension_factor(number val)
		{
			set_extension_factor(make_sp(new ConstUserNumber<dim>(val)));
		}	

		
		template<typename TDomain>
		void ExtensionEquation<TDomain>::
		prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
		{

		//	check that not ADAPTIVE
			if(vLfeID[0].order() < 1)
				UG_THROW("ExtensionEquation: Adaptive order not implemented.");

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
		void ExtensionEquation<TDomain>::
		prep_elem_loop(const ReferenceObjectID roid, const int si)
		{
	
		//	register subsetIndex at Geometry
			TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

			try{
				geo.update_local(roid, m_lfeID, m_quadOrder);
			}
			UG_CATCH_THROW("ExtensionEquation::prep_elem_loop:"
								" Cannot update Finite Element Geometry.");

		//  set local positions

			static const int refDim = TElem::dim;
			//ExtensionFactor
			m_imNuExtensionValue.template set_local_ips<refDim>(geo.local_ips(),geo.num_ip(),false);
			//AdjointDeformation Vector
			m_imAdjointDeformationd1.template set_local_ips<refDim>(geo.local_ips(),geo.num_ip(),false);
			m_imAdjointDeformationd2.template set_local_ips<refDim>(geo.local_ips(),geo.num_ip(),false);
			m_imAdjointDeformationd3.template set_local_ips<refDim>(geo.local_ips(),geo.num_ip(),false);
			//Deformation Vector
			m_imDeformationd1.template set_local_ips<refDim>(geo.local_ips(),geo.num_ip(),false);
			m_imDeformationd2.template set_local_ips<refDim>(geo.local_ips(),geo.num_ip(),false);
			m_imDeformationd3.template set_local_ips<refDim>(geo.local_ips(),geo.num_ip(),false);
			//DeformationGradient Matrix
			m_imDeformationVectord1.template set_local_ips<refDim>(geo.local_ips(),geo.num_ip(),false);
			m_imDeformationVectord2.template set_local_ips<refDim>(geo.local_ips(),geo.num_ip(),false);
			m_imDeformationVectord3.template set_local_ips<refDim>(geo.local_ips(),geo.num_ip(),false);
		}
		
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void ExtensionEquation<TDomain>::
		prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
		{
		//  update Geometry for this element
			TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);
			try{
				geo.update(elem, vCornerCoords, m_lfeID, m_quadOrder);
			}
			UG_CATCH_THROW("ExtensionEquation::prep_elem: "
								"Cannot update Finite Element Geometry.");
								
		//	set global positions

			//ExtensionFactor
			m_imNuExtensionValue.set_global_ips(geo.global_ips(),geo.num_ip());
			//AdjointDeformation vector
			m_imAdjointDeformationd1.set_global_ips(geo.global_ips(),geo.num_ip());
			m_imAdjointDeformationd2.set_global_ips(geo.global_ips(),geo.num_ip());
			m_imAdjointDeformationd3.set_global_ips(geo.global_ips(),geo.num_ip());
			//Deformation Vector
			m_imDeformationd1.set_global_ips(geo.global_ips(),geo.num_ip());
			m_imDeformationd2.set_global_ips(geo.global_ips(),geo.num_ip());
			m_imDeformationd3.set_global_ips(geo.global_ips(),geo.num_ip());
			//DeformationGradient Matrix
			m_imDeformationVectord1.set_global_ips(geo.global_ips(),geo.num_ip());
			m_imDeformationVectord2.set_global_ips(geo.global_ips(),geo.num_ip());
			m_imDeformationVectord3.set_global_ips(geo.global_ips(),geo.num_ip());
		}
		
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void ExtensionEquation<TDomain>::
		fsh_elem_loop()
		{	}
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void ExtensionEquation<TDomain>::
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
						J(_K_,sh1,_K_,sh2) += geo.shape(ip,sh1)
											*geo.shape(ip,sh2) 
											*geo.weight(ip);
					}//end sh2
				}//end sh1
			}//end ip 
			
		}//end jac func
		
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void ExtensionEquation<TDomain>::
		add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{
			TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

			for (size_t ip = 0; ip < geo.num_ip(); ++ip){
				const number ext = m_imNuExtensionValue[ip];
				//DEFORMATION IMPORT VECTOR AND GRADIENT
				//Declaration
				MathVector<dim> Deformation;VecSet(Deformation,0.0);	
				MathMatrix<dim, dim> DefGradient;MatSet(DefGradient,0.0);
				//Calculation
				DeformationVector(Deformation,ip);
				DeformationGradient(DefGradient, ip);
				//Adjoint deformation L
				MathVector<dim> vAdjDefVector;VecSet(vAdjDefVector,0.0);
				MathVector<dim> vTemp;VecSet(vTemp,0.0);
				AdjointDeformationVector(vAdjDefVector,ip);
				//Operations for dot product
				MatVecMult(vTemp,DefGradient,Deformation);
				number const scalar_product = VecDot(vTemp,vAdjDefVector);
				number const scale = m_chi*(ext-0.5*(m_nu_ext_upper+m_nu_ext_lower));
				
				for(size_t sh = 0; sh < geo.num_sh(); ++sh){
					d(_K_,sh) -= scale*geo.shape(ip,sh)*geo.weight(ip);
					d(_K_,sh) += scalar_product*geo.shape(ip,sh)*geo.weight(ip);
				}
			}//end ip 
		}
		
		////////////////////////////////////////////////////////////////////////////////
		//	Constructor
		////////////////////////////////////////////////////////////////////////////////

		template<typename TDomain>
		ExtensionEquation<TDomain>::ExtensionEquation(const char* function, const char* subsets)
		: IElemDisc<TDomain>(function, subsets),m_order(1),
		m_lfeID(LFEID::LAGRANGE, TDomain::dim, m_order)
		{
			this->clear_add_fct();
			
			//	set defaults
			m_order = 1;
			m_bQuadOrderUserDef = false;
			m_quadOrder = -1;
			
			m_nu_ext_lower=0.0;
			m_nu_ext_upper=0.0;
			m_chi=0.0;
			
			//register imports

			//ExtensionFactor
			this->register_import(m_imNuExtensionValue);
			//AdjointDeformationVector
			this->register_import(m_imAdjointDeformationd1);
			this->register_import(m_imAdjointDeformationd2);
			this->register_import(m_imAdjointDeformationd3);
			//DeformationVector
			this->register_import(m_imDeformationd1);
			this->register_import(m_imDeformationd2);
			this->register_import(m_imDeformationd3);
			//DeformationGradient
			this->register_import(m_imDeformationVectord1);
			this->register_import(m_imDeformationVectord2);
			this->register_import(m_imDeformationVectord3);

			
		}
		////////////////////////////////////////////////////////////////////////////////
		//	register assemble functions
		////////////////////////////////////////////////////////////////////////////////
		#ifdef UG_DIM_1
		template<>
		void ExtensionEquation<Domain1d>::register_all_funcs(int order, int quadOrder)
		{
			//	RegularEdge
			register_func<RegularEdge, DimFEGeometry<dim> > ();
		}
		#endif
		#ifdef UG_DIM_2
		template<>
		void ExtensionEquation<Domain2d>::register_all_funcs(int order, int quadOrder)
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
		void ExtensionEquation<Domain3d>::register_all_funcs(int order, int quadOrder)
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
		void ExtensionEquation<TDomain>::register_func()
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
		template class ExtensionEquation<Domain1d>;
		#endif
		#ifdef UG_DIM_2
		template class ExtensionEquation<Domain2d>;
		#endif
		#ifdef UG_DIM_3
		template class ExtensionEquation<Domain3d>;
		#endif
	 }//end namespace FluidOptim
 }//end namespace ug4