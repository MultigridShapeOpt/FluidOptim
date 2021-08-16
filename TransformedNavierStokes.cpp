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
 
 // for various user data
 //also taken from restricted_deformation_elasticity.cpp
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

#include "lib_disc/spatial_disc/disc_util/fe_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"
#include "lib_disc/local_finite_element/lagrange/lagrange.h"
#include "lib_disc/local_finite_element/lagrange/lagrangep1.h"
#include "lib_disc/quadrature/gauss/gauss_quad.h"

#include "TransformedNavierStokes.h"
 
 
 namespace ug{
	 namespace FluidOptim{

		//Calculation of the Identity Deformation Gradient
		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		IdentityDeformationGradient(MathMatrix<dim, dim>& GradU, const size_t ip)
		{
			//1.- calculate gradients of each deformation component u1,u2
			MathVector<dim> vdeformation_d1=m_imDeformationVectord1[ip];
			MathVector<dim> vdeformation_d2=m_imDeformationVectord2[ip];
			//2.- assign the gradients to each row using assign(vector, row)
			GradU.assign(vdeformation_d1,0);
			GradU.assign(vdeformation_d2,1);
			
			if(this->dim == 3){
				MathVector<dim> vdeformation_d3=m_imDeformationVectord3[ip];
				GradU.assign(vdeformation_d3,2);
			}
			
			//3.- add identity on the diagonal
			for(int i=0; i < dim; ++i)
			{
				GradU[i][i] += 1.0;
			}
		}
		///deformation gradient
		//For the d1
		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		set_deformation_vector_d1(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imDeformationVectord1.set_data(user);
		}
		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		set_deformation_vector_d1(number val)
		{
			set_deformation_vector_d1(make_sp(new ConstUserVector<dim>(val)));
		}
		#ifdef UG_FOR_LUA
		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		set_deformation_vector_d1(const char* fctName)
		{
			set_deformation_vector_d1(LuaUserDataFactory<MathVector<dim>, dim>::create(fctName));
		}

		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		set_deformation_vector_d1(LuaFunctionHandle fct)
		{
			set_deformation_vector_d1(make_sp(new LuaUserData<MathVector<dim>, dim>(fct)));
		}
		#endif
		
		//For the d2 
		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		set_deformation_vector_d2(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imDeformationVectord2.set_data(user);
		}
		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		set_deformation_vector_d2(number val)
		{
			set_deformation_vector_d2(make_sp(new ConstUserVector<dim>(val)));
		}
		#ifdef UG_FOR_LUA
		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		set_deformation_vector_d2(const char* fctName)
		{
			set_deformation_vector_d2(LuaUserDataFactory<MathVector<dim>, dim>::create(fctName));
		}

		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		set_deformation_vector_d2(LuaFunctionHandle fct)
		{
			set_deformation_vector_d2(make_sp(new LuaUserData<MathVector<dim>, dim>(fct)));
		}
		#endif
		//For the d3
		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		set_deformation_vector_d3(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imDeformationVectord3.set_data(user);
		}
		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		set_deformation_vector_d3(number val)
		{
			set_deformation_vector_d3(make_sp(new ConstUserVector<dim>(val)));
		}
		////////////////////////////////////////////////////////////////////////////////
		//	general
		////////////////////////////////////////////////////////////////////////////////
		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::set_quad_order(size_t order)
		{
			m_quadOrder = order;
			m_bQuadOrderUserDef = true;
		}
		
		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
		{
			if(bNonRegularGrid)
				UG_THROW("TransformedNavierStokes: only implemented for regular grids.");

			//	check number
			if(vLfeID.size() != dim+1)
				UG_THROW("TransformedNavierStokes: Needs exactly "<<dim+1<<" functions.");

			for(int d = 1; d < dim; ++d)
				if(vLfeID[0] != vLfeID[d])
					UG_THROW("TransformedNavierStokes: trial spaces for velocity expected to be"
					" identical for all velocity components.");

			//	remember lfeID;
			m_vLFEID = vLfeID[0];
			m_pLFEID = vLfeID[dim];

			if(!m_bQuadOrderUserDef) m_quadOrder = 2*m_vLFEID.order()+1;

			//	update assemble functions
			register_all_funcs(m_vLFEID, m_pLFEID, m_quadOrder);
		}
		// FUNCTIONS THAT NEED REGISTERING
		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		use_terms_on_current_domain()
		{
			this->normal_continuity(true);
			this->normal_convection(true);
			this->normal_diffusion(true);
			this->normal_pressure(true);
			this->normal_stabilization(true);
			//on jacobian
			this->normal_jacobian_continuity(true);
			this->normal_jacobian_diffusion(true);
			this->normal_jacobian_newtonderivative(true);
			this->normal_jacobian_pressure(true);
			this->normal_jacobian_stabilization(true);
			this->normal_jacobian_vectorconvection(true);
		}
		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		set_nonlinear(bool option)
		{
			this->m_bNonLinear=option;
		}
		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		set_picard(bool option)
		{
			this->m_bPicard=option;
		}
		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		no_deformation(bool def)
		{
			this->m_bNoDeformation=def;
		}

		
		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		normal_diffusion(bool dif)
		{
			this->m_bNormalDiffusion=dif;
		}
		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		normal_continuity(bool cont)
		{
			this->m_bNormalContinuity=cont;
		}
		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		normal_pressure(bool pres)
		{
			this->m_bNormalPressure=pres;
		}
		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		normal_stabilization(bool stab)
		{
			this->m_bNormalStabilization=stab;
		}
		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		normal_convection(bool conv)
		{
			this->m_bNormalConvection=conv;
		}
		///JACOBIAN NORMAL WASTE OF TIME
		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		normal_jacobian_diffusion(bool dif)
		{
			this->m_bNormalJacobianDiffusion=dif;
		}
		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		normal_jacobian_pressure(bool press)
		{
			this->m_bNormalJacobianPressure=press;
		}
		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		normal_jacobian_continuity(bool cont)
		{
			this->m_bNormalJacobianContinuity=cont;
		}
		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		normal_jacobian_stabilization(bool stab)
		{
			this->m_bNormalJacobianStabilization=stab;
		}
		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		normal_jacobian_vectorconvection(bool vec)
		{
			this->m_bNormalJacobianVectorConvection=vec;
		}
		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		normal_jacobian_newtonderivative(bool nton)
		{
			this->m_bNormalJacobianNewtonDerivative=nton;
		}
		////IMPORTS
		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imKinViscosity.set_data(user);
		}
		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		set_kinematic_viscosity(number val)
		{
			if(val == 0.0) set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> >());
			else set_kinematic_viscosity(make_sp(new ConstUserNumber<dim>(val)));
		}
		#ifdef UG_FOR_LUA
		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		set_kinematic_viscosity(const char* fctName)
		{
			set_kinematic_viscosity(LuaUserDataFactory<number,dim>::create(fctName));
		}

		template<typename TDomain>
		void TransformedNavierStokes<TDomain>::
		set_kinematic_viscosity(LuaFunctionHandle fct)
		{
			set_kinematic_viscosity(make_sp(new LuaUserData<number,dim>(fct)));
		}
		#endif
		
		////////////////////////////////////////////////////////////////////////////////
		// Assembling functions
		////////////////////////////////////////////////////////////////////////////////
		//Element looping functions
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void TransformedNavierStokes<TDomain>::prep_elem_loop(const ReferenceObjectID roid, const int si)
		{
			DimFEGeometry<dim>& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);
			DimFEGeometry<dim>& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_quadOrder);
			try
			{
				vgeo.update_local(roid, m_vLFEID, m_quadOrder);
				pgeo.update_local(roid, m_pLFEID, m_quadOrder);
			}
			UG_CATCH_THROW("FluidOptim::TransformedNavierStokes: "
			"Cannot update Finite Element Geometry.");
			m_imKinViscosity.set_local_ips(vgeo.local_ips(), vgeo.num_ip(),false);
			
			m_imDeformationVectord1.set_local_ips(vgeo.local_ips(), vgeo.num_ip(),false);
			m_imDeformationVectord2.set_local_ips(vgeo.local_ips(), vgeo.num_ip(),false);
			m_imDeformationVectord3.set_local_ips(vgeo.local_ips(), vgeo.num_ip(),false);
		}
		
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void TransformedNavierStokes<TDomain>::prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
		{
			m_pElem = elem;

			// 	Update Geometry for this element
			DimFEGeometry<dim>& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);
			DimFEGeometry<dim>& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_quadOrder);
			try{
				vgeo.update(elem, vCornerCoords, m_vLFEID, m_quadOrder);
				pgeo.update(elem, vCornerCoords, m_pLFEID, m_quadOrder);
			}
			UG_CATCH_THROW("TransformedNavierStokes: Cannot update "
							"Finite Element Geometry.");
			static const int refDim = TElem::dim;
			m_imKinViscosity.set_global_ips(vgeo.global_ips(), vgeo.num_ip());

			m_imDeformationVectord1.set_global_ips(vgeo.global_ips(), vgeo.num_ip());
			m_imDeformationVectord2.set_global_ips(vgeo.global_ips(), vgeo.num_ip());
			m_imDeformationVectord3.set_global_ips(vgeo.global_ips(), vgeo.num_ip());
		}
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void TransformedNavierStokes<TDomain>::fsh_elem_loop()
		{
		}
		///	assembles the local stiffness matrix 
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void TransformedNavierStokes<TDomain>::
		add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{
		 	//	request geometry
					const DimFEGeometry<dim>& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);
					const DimFEGeometry<dim>& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_quadOrder);
					
										
					for (size_t ip = 0; ip < vgeo.num_ip(); ++ip){

						//TODO:verify that approach is correct, jacobian
						MathMatrix<dim, dim> UnitDeformationGradient;MatSet(UnitDeformationGradient,0.0);
						MathMatrix<dim, dim> InvUnitDeformationGradient;MatSet(InvUnitDeformationGradient,0.0);
						number dDF=0.0;
						
						IdentityDeformationGradient(UnitDeformationGradient, ip);
						Inverse(InvUnitDeformationGradient,UnitDeformationGradient);
						dDF=Determinant(UnitDeformationGradient);
						
												
						////////////////////////////////////////////////////
						// Diffusive Term (Momentum Equation)
						////////////////////////////////////////////////////

						const number scale = m_imKinViscosity[ip]* vgeo.weight(ip);
						if(m_bNormalJacobianDiffusion == true)
						{
							for (int vdim = 0; vdim < dim; ++vdim){
								for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){

									for (size_t ush = 0; ush < vgeo.num_sh(); ++ush){
										for (int udim = 0; udim < dim; ++udim) {

											J(vdim, vsh, vdim, ush) -=  scale
																	   * vgeo.global_grad(ip, ush)[udim]
																	   * vgeo.global_grad(ip, vsh)[udim];


										}
									}
								}
							}
						}
						if(m_bNormalJacobianDiffusion == false)
						{
							
							/* for (int vdim = 0; vdim < dim; ++vdim){
								for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){

									for (size_t ush = 0; ush < vgeo.num_sh(); ++ush){
										
										MathVector<dim> v_TransformedGradient1; VecSet(v_TransformedGradient1,0.0);								
										MathVector<dim> v_TransformedGradient2; VecSet(v_TransformedGradient2,0.0);

										TransposedMatVecMult(v_TransformedGradient1,InvUnitDeformationGradient,vgeo.global_grad(ip,ush));
										TransposedMatVecMult(v_TransformedGradient2,InvUnitDeformationGradient,vgeo.global_grad(ip,vsh));
										for (int udim = 0; udim < dim; ++udim) {

											J(vdim, vsh, vdim, ush) +=  scale
																	   * v_TransformedGradient1[udim]
																	   * v_TransformedGradient2[udim]*dDF;


										}
									}
								}
							} */
							for(int d1=0; d1 <  dim; ++d1){
								for(size_t sh1=0;sh1<vgeo.num_sh();++sh1){
									for(size_t sh2=0;sh2<vgeo.num_sh();++sh2){
										for(int d2=0 ; d2< dim ;++d2){

											MathMatrix<dim,dim> M1;MatSet(M1,0.0);
											MathMatrix<dim,dim> M2;MatSet(M2,0.0);
											M1.assign(vgeo.global_grad(ip,sh1),d1);
											M2.assign(vgeo.global_grad(ip,sh2),d2);
											
											MathMatrix<dim,dim> TransM1;MatSet(TransM1,0.0);
											MathMatrix<dim,dim> TransM2;MatSet(TransM2,0.0);
											MatMultiply(TransM1,M1,InvUnitDeformationGradient);
											MatMultiply(TransM2,M2,InvUnitDeformationGradient);
											
											number contraction = MatContraction(TransM1,TransM2);

											J(d1,sh1,d2,sh2) -= contraction*scale*dDF;
										}//d2
									}//end sh2
								}//end sh1
							}//d1
						}
						if(m_bNonLinear == true)
						{
							////////////////////////////////////////////////////
							// Linearization Terms
							////////////////////////////////////////////////////
							
							////////////////////////////////////////////////////
							//Vector-Convection Matrix (N)
							////////////////////////////////////////////////////
							//interpolate velocity at ip, assign quasi-laplacian terms in d1,d1
							
							if(m_bNormalJacobianVectorConvection == true)
							{
								// 	Interpolate Velocity at ip
								MathVector<dim> Vel;
								for(int d = 0; d < dim; ++d){
									Vel[d] = 0.0;
									for(size_t sh = 0; sh < vgeo.num_sh(); ++sh)
										Vel[d] += u(d, sh) * vgeo.shape(ip, sh);
								}

								// linearization of u \nabla u w.r.t second u (i.e. keeping first as old iterate)
								for (int vdim = 0; vdim < dim; ++vdim){
									for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){

										for (size_t ush = 0; ush < vgeo.num_sh(); ++ush){
											for (int udim = 0; udim < dim; ++udim) {

												J(vdim, vsh, vdim, ush) -= Vel[udim]
																		   * vgeo.global_grad(ip, ush)[udim]
																		   * vgeo.shape(ip, vsh)
																		   * vgeo.weight(ip);
											}
										}
									}
								}
							}
							if(m_bNormalJacobianVectorConvection == false)
							{
								// 	Interpolate Velocity at ip
								MathVector<dim> Vel;
								for(int d = 0; d < dim; ++d)
								{
									Vel[d] = 0.0;
									for(size_t sh = 0; sh < vgeo.num_sh(); ++sh)
										Vel[d] += u(d, sh) * vgeo.shape(ip, sh);
								}
								for(int vdim=0; vdim < dim ; vdim++){
									for(size_t vsh=0; vsh < vgeo.num_sh() ;++vsh){
										for(size_t ush=0; ush < vgeo.num_sh(); ush++){

											MathVector<dim> v_TransformedGradient; VecSet(v_TransformedGradient,0.0);								
											TransposedMatVecMult(v_TransformedGradient,InvUnitDeformationGradient,vgeo.global_grad(ip,ush));
											
											for(int udim=0; udim < dim; udim++){
												J(vdim,vsh,vdim,ush) -= Vel[udim]*v_TransformedGradient[udim]*vgeo.shape(ip,vsh)*vgeo.weight(ip)*dDF;
											}
										}
									}

								}//end sh1
							}
							
							if(m_bPicard == false)
							{
								////////////////////////////////////////////////////
								//Newton-Derivative Matrix 
								////////////////////////////////////////////////////
								
								if(m_bNormalJacobianNewtonDerivative == true)
								{
									// 	Interpolate Functional Matrix of velocity at ip
									MathMatrix<dim, dim> gradVel;
									for(int d1 = 0; d1 < dim; ++d1){
										for(int d2 = 0; d2 <dim; ++d2){
											gradVel(d1, d2) = 0.0;
											for(size_t sh = 0; sh < vgeo.num_sh(); ++sh)
												gradVel(d1, d2) += u(d1, sh) * vgeo.global_grad(ip, sh)[d2];
										}
									}

									// linearization of u \nabla u w.r.t first u (i.e. keeping second as old iterate)
									for (int vdim = 0; vdim < dim; ++vdim){
										for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){

											for (size_t ush = 0; ush < vgeo.num_sh(); ++ush){
												for (int udim = 0; udim < dim; ++udim) {

													J(vdim, vsh, udim, ush) -= gradVel(vdim, udim)
																			   * vgeo.shape(ip, ush)
																			   * vgeo.shape(ip, vsh)
																			   * vgeo.weight(ip);
												}
											}
										}
									}
								}
								if(m_bNormalJacobianNewtonDerivative == false)
								{
									// 	Interpolate Functional Matrix of velocity at ip
									MathMatrix<dim, dim> gradVel;
									for(int d1 = 0; d1 < dim; ++d1){
										for(int d2 = 0; d2 <dim; ++d2){
											gradVel(d1, d2) = 0.0;
											for(size_t sh = 0; sh < vgeo.num_sh(); ++sh)
												gradVel(d1, d2) += u(d1, sh) * vgeo.global_grad(ip, sh)[d2];
										}
									}
																		
									MathMatrix<dim,dim> TransformedGradVel; MatSet(TransformedGradVel, 0.0);
									MatMultiply(TransformedGradVel,gradVel,InvUnitDeformationGradient);
									for(int vdim=0 ; vdim < dim; vdim++){
										for(size_t vsh=0 ; vsh < vgeo.num_sh() ; vsh++){
											for(size_t ush=0 ; ush<vgeo.num_sh() ; ush++){
												for(int udim=0 ; udim < dim ; udim++){
													J(vdim,vsh,udim,ush) -= TransformedGradVel(vdim,udim)
																			*vgeo.shape(ip,ush)
																			*vgeo.shape(ip,vsh)
																			*vgeo.weight(ip)
																			*dDF;
												}
											}
										}

									}//end sh1
								}
							}
							
						}

						////////////////////////////////////////////////////
						// Pressure Term (Momentum Equation)
						////////////////////////////////////////////////////
						
						if(m_bNormalJacobianPressure == true)
						{
							for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){
								for (int vdim = 0; vdim < dim; ++vdim){
									for (size_t psh = 0; psh < pgeo.num_sh(); ++psh)
										J(vdim, vsh, _P_, psh) +=     pgeo.shape(ip, psh)
																	* vgeo.global_grad(ip, vsh)[vdim]
																	* vgeo.weight(ip);
								}
							}
						}
						if(m_bNormalJacobianPressure == false)
						{
							for(size_t vsh=0;vsh<vgeo.num_sh();vsh++){
								for(int vdim=0 ; vdim < dim ; ++vdim){
									for(size_t psh=0 ; psh < pgeo.num_sh() ; ++psh){
										
										MathMatrix<dim, dim> GradTestV;MatSet(GradTestV,0.0);
										MathMatrix<dim, dim> TransformedGradTestV;MatSet(TransformedGradTestV,0.0);
										GradTestV.assign(vgeo.global_grad(ip,vsh),vdim);
										MatMultiply(TransformedGradTestV,GradTestV,InvUnitDeformationGradient);
										number TrTGradTestV = Trace(TransformedGradTestV);
										
										number trace = 0.0;//gonna store here the dot product between the deformations in Xi direction and the gradient of the shape func
										
										for(size_t i=0 ; i < (size_t) dim; ++i){
											trace += vgeo.global_grad(ip,vsh)[i]*InvUnitDeformationGradient[i][vdim];
										}
										//J(vdim,vsh,_P_,psh) -= pgeo.shape(ip,psh)*trace*vgeo.weight(ip)*dDF;
										J(vdim,vsh,_P_,psh) += pgeo.shape(ip,psh)*TrTGradTestV*vgeo.weight(ip)*dDF;
									}//psh
								}//vdim		
							}//vsh
						}
						

						

						////////////////////////////////////////////////////
						// Continuity Equation (conservation of mass)
						////////////////////////////////////////////////////
						if(m_bNormalJacobianContinuity == true)
						{
							for (size_t psh = 0; psh < pgeo.num_sh(); ++psh){
								for (size_t ush = 0; ush < vgeo.num_sh(); ++ush){
									for (int udim = 0; udim < dim; ++udim) {
											J(_P_, psh, udim, ush) +=
															 vgeo.global_grad(ip, ush)[udim]
														   * pgeo.shape(ip, psh)
														   * vgeo.weight(ip);
									}
								}
							}
						}
						if(m_bNormalJacobianContinuity == false)
						{
							//calculate the trace before assigning value
							for(size_t psh=0; psh < pgeo.num_sh() ; ++psh)
							{
								for(size_t vsh=0; vsh < vgeo.num_sh() ; ++vsh)
								{
									for(size_t vdim=0; vdim < (size_t) dim; ++vdim)
									{
										MathMatrix<dim, dim> GradTestV;MatSet(GradTestV,0.0);
										MathMatrix<dim, dim> TransformedGradTestV;MatSet(TransformedGradTestV,0.0);
										GradTestV.assign(vgeo.global_grad(ip,vsh),vdim);
										MatMultiply(TransformedGradTestV,GradTestV,InvUnitDeformationGradient);
										number TrTGradTestV = Trace(TransformedGradTestV);
										
										number trace=0.0;
										for(size_t i=0;i<(size_t) dim;++i){//physically meaningless index
											trace += vgeo.global_grad(ip,vsh)[i]*InvUnitDeformationGradient[i][vdim];
										}
											//J(_P_,psh,vdim,vsh) -= trace*pgeo.shape(ip,psh)*vgeo.weight(ip)*dDF;
											J(_P_,psh,vdim,vsh) += TrTGradTestV*pgeo.shape(ip,psh)*vgeo.weight(ip)*dDF;

									}
								}
							}
						}
						
					}
					if(m_stab_type ==  0.0 && m_stabParam != 0.0)
					{
							const number scale = m_stabParam * ElementDiameterSq<GridObject, TDomain>(*m_pElem, *this->domain());
							if(m_bNormalJacobianStabilization == true)
							{
								for (size_t ip = 0; ip < pgeo.num_ip(); ++ip){
									for (size_t psh = 0; psh < pgeo.num_sh(); ++psh){
										for (size_t psh2 = 0; psh2 < pgeo.num_sh(); ++psh2){
											J(_P_, psh,_P_, psh2) += scale
														   * VecDot(pgeo.global_grad(ip, psh), pgeo.global_grad(ip, psh2))
														   * pgeo.weight(ip);
										}
									}
								}
							}
						
							if(m_bNormalJacobianStabilization == false)
							{
								for(size_t ip=0; ip<pgeo.num_ip();++ip)
								{								
									//TODO:verify that approach is correct, jacobian
									MathMatrix<dim, dim> UnitDeformationGradient;
									MathMatrix<dim, dim> InvUnitDeformationGradient;
									number dDF=0.0;
									
									IdentityDeformationGradient(UnitDeformationGradient, ip);
									Inverse(InvUnitDeformationGradient,UnitDeformationGradient);
									dDF=Determinant(UnitDeformationGradient);
									
									MathVector<dim> v_tempVec1;VecSet(v_tempVec1,0.0);
									MathVector<dim> v_tempVec2;VecSet(v_tempVec2,0.0);
									for(size_t psh1=0;psh1<pgeo.num_sh();++psh1)
									{
										TransposedMatVecMult(v_tempVec1,InvUnitDeformationGradient,pgeo.global_grad(ip,psh1));
										
										for(size_t psh2=0;psh2<pgeo.num_sh();psh2++){
											TransposedMatVecMult(v_tempVec2,InvUnitDeformationGradient,pgeo.global_grad(ip,psh2));
											J(_P_,psh1,_P_,psh2)+= scale*VecDot(v_tempVec1,v_tempVec2)*pgeo.weight(ip)*dDF;
										}
									}
								}
							}											
					}//case pressure_gradient
					else if(m_stab_type != 0.0)
					{
						for(size_t ip=0; ip<pgeo.num_ip();++ip)
						{
							const number scale = m_stab_type/m_imKinViscosity[ip];
							//TODO:verify that approach is correct, jacobian
							MathMatrix<dim, dim> UnitDeformationGradient;
							MathMatrix<dim, dim> InvUnitDeformationGradient;
							number dDF=0.0;
								
							IdentityDeformationGradient(UnitDeformationGradient, ip);
							Inverse(InvUnitDeformationGradient,UnitDeformationGradient);
							dDF=Determinant(UnitDeformationGradient);
								
							for(size_t psh1=0;psh1<pgeo.num_sh();++psh1)
							{	
								for(size_t psh2=0;psh2<pgeo.num_sh();psh2++){
									J(_P_,psh1,_P_,psh2)+= scale*(pgeo.shape(ip,psh1)*pgeo.shape(ip,psh2)-1.0/(pgeo.num_sh()*pgeo.num_sh()))*pgeo.weight(ip)*dDF;
								}
							}
						}						
					}//case pressure_projection

		}
		
		
		///	assembles the stiffness part of the local defect, for nonlinear stuff
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void TransformedNavierStokes<TDomain>::
		add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{
			//	request geometry
				const DimFEGeometry<dim>& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);
				const DimFEGeometry<dim>& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_quadOrder);

				// loop integration points, note: pgeo and vgeo have same ip
				for (size_t ip = 0; ip < vgeo.num_ip(); ++ip){
					
					//TODO:verify that approach is correct
					MathMatrix<dim, dim> UnitDeformationGradient;MatSet(UnitDeformationGradient,0.0);
					MathMatrix<dim, dim> InvUnitDeformationGradient;MatSet(InvUnitDeformationGradient,0.0);
					number dDF=0.0;
					
					IdentityDeformationGradient(UnitDeformationGradient, ip);
					Inverse(InvUnitDeformationGradient,UnitDeformationGradient);
					dDF=Determinant(UnitDeformationGradient);
					
					// 	Interpolate Functional Matrix of velocity at ip
					MathMatrix<dim, dim> gradVel;MatSet(gradVel,0.0);
					for(int d1 = 0; d1 < dim; ++d1){
						for(int d2 = 0; d2 <dim; ++d2){
							gradVel(d1, d2) = 0.0;
							for(size_t sh = 0; sh < vgeo.num_sh(); ++sh)
								gradVel(d1, d2) += u(d1, sh) * vgeo.global_grad(ip, sh)[d2];
						}
					}

					number divu = 0.0;
					for (size_t ush = 0; ush < vgeo.num_sh(); ++ush){
						for (int udim = 0; udim < dim; ++udim) {
							divu += u(udim, ush) * vgeo.global_grad(ip, ush)[udim];
						}
					}
					
					//using the new approach
					MathMatrix<dim,dim> TransformedGradVel;MatSet(TransformedGradVel, 0.0);
					MatMultiply(TransformedGradVel,gradVel,InvUnitDeformationGradient);
					MathMatrix<dim,dim> Difference;//new
					number trace_div = Trace(TransformedGradVel);
					
					/*number trace_div = 0.0;
					
					 for (size_t ush = 0; ush < vgeo.num_sh(); ++ush){
						
						MathVector<dim> v_TransformedGradient; VecSet(v_TransformedGradient,0.0);
						MathMatrix<dim, dim> TranposedInverse; Transpose( TranposedInverse , InvUnitDeformationGradient );
						MatVecMult(v_TransformedGradient,TranposedInverse,vgeo.global_grad(ip,ush));
						
						for (int udim = 0; udim < dim; ++udim) {
							trace_div += u(udim, ush) * v_TransformedGradient[udim];
						}
					} */
					
					//TODO:remove this, just for debugging
					if(m_bCompareDeformation == true)
					{
						MatSubtract(Difference,TransformedGradVel,gradVel);
						if(Trace(Difference) != 0.0)
							UG_THROW("TransformedNavierStokes: The trace of the TransformedGradVel-gradVel with new approach != 0.0");
					}
					
					////////////////////////////////////////////////////
					// Diffusive Term (Momentum Equation)
					////////////////////////////////////////////////////
					if(m_bNormalDiffusion == true){
						const number scale = m_imKinViscosity[ip] * vgeo.weight(ip);
						for (int vdim = 0; vdim < dim; ++vdim){
							for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){
								for (int udim = 0; udim < dim; ++udim) {

									d(vdim, vsh) -=  scale * gradVel(vdim, udim)
												   * vgeo.global_grad(ip, vsh)[udim];
								}
							}
						}
					}
					if(m_bNormalDiffusion == false){
						const number scale = m_imKinViscosity[ip] * vgeo.weight(ip);
						for (int vdim = 0; vdim < dim; ++vdim){
							for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){
								
								//Transformation of vector
								MathVector<dim> v_TransformedGradient; VecSet(v_TransformedGradient,0.0);
								MathMatrix<dim, dim> TranposedInverse; Transpose( TranposedInverse , InvUnitDeformationGradient );
								MatVecMult(v_TransformedGradient,TranposedInverse,vgeo.global_grad(ip,vsh));
								
								MathMatrix<dim, dim> GradTestV;MatSet(GradTestV,0.0);
								MathMatrix<dim, dim> TransformedGradTestV;MatSet(TransformedGradTestV,0.0);
								GradTestV.assign(vgeo.global_grad(ip,vsh),vdim);
								MatMultiply(TransformedGradTestV,GradTestV,InvUnitDeformationGradient);
								number contract = MatContraction(TransformedGradVel,TransformedGradTestV);
								d(vdim, vsh) -= scale*contract*dDF;
	
							}
						}
						
					}
					if(m_bNonLinear)
					{
						if(m_bNormalConvection == true)
						{
						// 	Interpolate Velocity at ip
							MathVector<dim> Vel;VecSet(Vel,0.0);
							for(int d1 = 0; d1 < dim; ++d1){
								Vel[d1] = 0.0;
								for(size_t sh = 0; sh < vgeo.num_sh(); ++sh)
									Vel[d1] += u(d1, sh) * vgeo.shape(ip, sh);
							}

							MathVector<dim> convFlux;VecSet(convFlux,0.0);
							MatVecMult(convFlux, gradVel, Vel);

							for (int vdim = 0; vdim < dim; ++vdim){
								for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){
									d(vdim, vsh) -= convFlux[vdim]
												   * vgeo.shape(ip, vsh) * vgeo.weight(ip);
								}
							}
						}
						if(m_bNormalConvection == false)
						{
							////////////////////////////////////////////////////
							// Convective Term (Momentum Equation)
							////////////////////////////////////////////////////
							// 	Interpolation of the Velocity at ip
							MathVector<dim> Vel;VecSet(Vel,0.0);
							for(int d1 = 0; d1 < dim; ++d1){
								Vel[d1] = 0.0;
								for(size_t sh = 0; sh < vgeo.num_sh(); ++sh)
									Vel[d1] += u(d1, sh) * vgeo.shape(ip, sh);
							}

							MathVector<dim> TConvFlux;VecSet(TConvFlux,0.0);
							MatVecMult(TConvFlux, TransformedGradVel, Vel);	
							for (int vdim = 0; vdim < dim; ++vdim){
								for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){
									d(vdim, vsh) -=  TConvFlux[vdim]* vgeo.shape(ip, vsh) * vgeo.weight(ip) * dDF;
								}
							}
						}
						
					}
				
					////////////////////////////////////////////////////
					// Pressure Term (Momentum Equation)
					////////////////////////////////////////////////////
					number pressure = 0.0;
					for (size_t psh = 0; psh < pgeo.num_sh(); ++psh)
						pressure += u(_P_, psh) * pgeo.shape(ip, psh);
					
					if(m_bNormalPressure == true)
					{
						for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){
							for (int vdim = 0; vdim < dim; ++vdim){
								d(vdim, vsh) += pressure
											   * vgeo.global_grad(ip, vsh)[vdim]
											   * vgeo.weight(ip);
							}
						}
					}
					if(m_bNormalPressure == false)
					{
						for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){
							for (int vdim = 0; vdim < dim; ++vdim){
								MathMatrix<dim, dim> GradTestV;MatSet(GradTestV,0.0);
								MathMatrix<dim, dim> TransformedGradTestV;MatSet(TransformedGradTestV,0.0);
								GradTestV.assign(vgeo.global_grad(ip,vsh),vdim);
								MatMultiply(TransformedGradTestV,GradTestV,InvUnitDeformationGradient);
								
								number TrTGradTestV = Trace(TransformedGradTestV);
								
								/*number trace=0.0;//gonna store here the dot product between the deformations in Xi direction and the gradient of the shape func
								 for(size_t i=0;i< (size_t)dim;++i){
									trace += vgeo.global_grad(ip,vsh)[i]*InvUnitDeformationGradient[i][vdim];
								} */
								d(vdim, vsh) += pressure * TrTGradTestV *vgeo.weight(ip)*dDF;
								//d(vdim, vsh) -= pressure*trace*vgeo.weight(ip)*dDF;
							}
						}
					}
					

					
					

					////////////////////////////////////////////////////
					// Continuity Equation (conservation of mass)
					////////////////////////////////////////////////////
					if(m_bNormalContinuity == true)
					{
						for (size_t psh = 0; psh < pgeo.num_sh(); ++psh){
							d(_P_, psh) += divu * pgeo.shape(ip, psh)* vgeo.weight(ip);
						}
					}
					if(m_bNormalContinuity == false)
					{
							for (size_t psh = 0; psh < pgeo.num_sh(); ++psh){
								d(_P_, psh) += pgeo.shape(ip, psh)*trace_div*vgeo.weight(ip)*dDF;
							}
					}
					
				}

				// stabilization
				if(m_stab_type == 0.0 && m_stabParam != 0.0)
				{
					MathMatrix<dim, dim> DF_inv;//for storing the gradient of the deformation and inverse

					const number scale = m_stabParam * ElementDiameterSq<GridObject, TDomain>(*m_pElem, *this->domain());
					if(m_bNormalStabilization == true)
					{
						for (size_t ip = 0; ip < pgeo.num_ip(); ++ip){

							MathVector<dim> pressGrad; VecSet(pressGrad, 0.0);
							for (size_t psh = 0; psh < pgeo.num_sh(); ++psh)
								VecScaleAppend(pressGrad, u(_P_, psh), pgeo.global_grad(ip, psh));

							for (size_t psh = 0; psh < pgeo.num_sh(); ++psh){
								d(_P_, psh) += scale
											   * VecDot(pgeo.global_grad(ip, psh), pressGrad)
											   * pgeo.weight(ip);
							}
						}
					}//end m_bNormalStabilization
					
					if(m_bNormalStabilization ==  false)
					{
						for (size_t ip = 0; ip < pgeo.num_ip(); ++ip){
								//TODO:verify that approach is correct
								MathMatrix<dim, dim> UnitDeformationGradient;MatSet(UnitDeformationGradient,0.0);
								MathMatrix<dim, dim> InvUnitDeformationGradient;MatSet(InvUnitDeformationGradient,0.0);
								number dDF=0.0;
					
								IdentityDeformationGradient(UnitDeformationGradient, ip);
								Inverse(InvUnitDeformationGradient,UnitDeformationGradient);
								dDF=Determinant(UnitDeformationGradient);
			
								//Calculate the pressure gradient
								MathVector<dim> vPressureGradient; VecSet(vPressureGradient, 0.0);
								MathVector<dim> vTransformedPressureGrad;VecSet(vTransformedPressureGrad,0.0);
								for (size_t psh = 0; psh < pgeo.num_sh(); ++psh)
									VecScaleAppend(vPressureGradient, u(_P_, psh), pgeo.global_grad(ip, psh));
								//Transform the gradient
								
								TransposedMatVecMult(vTransformedPressureGrad,InvUnitDeformationGradient,vPressureGradient);
								
								//We have to do the same for the shape functions
								
								
								for (size_t psh = 0; psh < pgeo.num_sh(); ++psh){
									MathVector<dim> vTransformedGrad;VecSet(vTransformedGrad,0.0);
									TransposedMatVecMult(vTransformedGrad,InvUnitDeformationGradient,pgeo.global_grad(ip, psh));
									d(_P_, psh) += scale* VecDot(vTransformedGrad, vTransformedPressureGrad)* pgeo.weight(ip)*dDF;
								}								
						}// ps ip
					
					}//end transformed grad stab
				}//case pressure_gradient
				else if (m_stab_type != 0.0)
				{
					number pressure_average = 0.0;
					for (size_t psh = 0; psh < pgeo.num_sh(); ++psh) {
						pressure_average +=  1.0/pgeo.num_sh()*u(_P_, psh);
					}
					for (size_t ip = 0; ip < pgeo.num_ip(); ++ip){
						
						const number scale = m_stab_type/m_imKinViscosity[ip];
						//TODO:verify that approach is correct
						MathMatrix<dim, dim> UnitDeformationGradient;MatSet(UnitDeformationGradient,0.0);
						MathMatrix<dim, dim> InvUnitDeformationGradient;MatSet(InvUnitDeformationGradient,0.0);
						number dDF=0.0;
					
						IdentityDeformationGradient(UnitDeformationGradient, ip);
						Inverse(InvUnitDeformationGradient,UnitDeformationGradient);
						dDF=Determinant(UnitDeformationGradient);
								
						number pressure = 0.0;
						for (size_t psh = 0; psh < pgeo.num_sh(); ++psh)
							pressure += u(_P_, psh) * pgeo.shape(ip, psh);
							
						for (size_t psh = 0; psh < pgeo.num_sh(); ++psh){
							//d(_P_, psh) += m_imKinViscosity[ip]*pgeo.shape(ip,psh)*pressure*(1.0 - 1.0/(pgeo.num_sh()*pgeo.num_sh())) * pgeo.weight(ip)*dDF;
							//d(_P_, psh) += pressure*(pgeo.shape(ip,psh)*pgeo.shape(ip,psh) - (1.0/(pgeo.num_sh()*pgeo.num_sh()))) * pgeo.weight(ip)*dDF;
							d(_P_, psh) += scale*(pressure * pgeo.shape(ip,psh) - pressure_average/pgeo.num_sh()) * pgeo.weight(ip)*dDF;
						}
					}//ip end
				}//case pressure_projection

		}
		///	assembles the local right hand side
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void TransformedNavierStokes<TDomain>::add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{
			
		}
		//Mass Matrix methods will be assembled if there is time dependency
		///	assembles the local mass matrix using a finite volume scheme
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void TransformedNavierStokes<TDomain>::add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{
		}
		///	assembles the mass part of the local defect
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void TransformedNavierStokes<TDomain>::add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{
		}
		
		////////////////////////////////////////////////////////////////////////////////
		//	Constructor
		////////////////////////////////////////////////////////////////////////////////

		template<typename TDomain>
		TransformedNavierStokes<TDomain>::
		TransformedNavierStokes(const char* functions, const char* subsets)
		 : IElemDisc<TDomain>(functions,subsets)
		{
		//	set defaults
			m_bQuadOrderUserDef = false;
			m_quadOrder = -1;
		//	check number of functions
			if(this->num_fct() != dim+1)
				UG_THROW("Wrong number of functions: The ElemDisc 'TransformedNavierStokes'"
					   " needs exactly "<<dim+1<<" symbolic function.");
		//Default value assigned	   
		m_stabParam=0.0;// no stab
		m_stab_type=0.0;
		m_bNonLinear=false;
		m_bPicard=false;
		m_bNormalContinuity=false;
		m_bNormalDiffusion=false;
		m_bNormalPressure=false;
		m_bNormalStabilization=false;
		m_bNormalConvection=false;
		//jacobian defaults
		m_bNormalJacobianDiffusion=false;
		m_bNormalJacobianPressure=false;
		m_bNormalJacobianContinuity=false;
		m_bNormalJacobianStabilization=false;
		m_bNormalJacobianNewtonDerivative=false;
		m_bNormalJacobianVectorConvection=false;
		
		//deformation vectors debugging
		m_bCompareDeformation=false;
		m_deformationValue=0.0;
		
		//	register imports
			this->register_import(m_imKinViscosity);
			
			//For the deformation vectors d1 d2 d3
			this->register_import(m_imDeformationVectord1);
			this->register_import(m_imDeformationVectord2);
			this->register_import(m_imDeformationVectord3);
			
			this->clear_add_fct();
		}
		
		template<typename TDomain>
		TransformedNavierStokes<TDomain>::
		TransformedNavierStokes(const std::vector<std::string>& vFct,
                                    const std::vector<std::string>& vSubset)
		 : IElemDisc<TDomain>(vFct,vSubset)
		{
		//	set defaults
			m_bQuadOrderUserDef = false;
			m_quadOrder = -1;
			
		//	check number of functions
			if(this->num_fct() != dim+1)
				UG_THROW("Wrong number of functions: The ElemDisc 'TransformedNavierStokes'"
					   " needs exactly "<<dim+1<<" symbolic function.");
			//Default value assigned	   
			m_stabParam=0.0;
			m_stab_type=0.0;
			m_bNonLinear=false;
			m_bPicard=false;
			m_bNormalContinuity=false;
			m_bNormalDiffusion=false;
			m_bNormalPressure=false;
			m_bNormalStabilization=false;
			m_bNormalConvection=false;
			//jacobian defaults
			m_bNormalJacobianDiffusion=false;
			m_bNormalJacobianPressure=false;
			m_bNormalJacobianContinuity=false;
			m_bNormalJacobianStabilization=false;
			m_bNormalJacobianNewtonDerivative=false;
			m_bNormalJacobianVectorConvection=false;
			//deformation vectors debugging
			m_bCompareDeformation=false;
			m_deformationValue=0.0;
		//	register imports
			this->register_import(m_imKinViscosity);

			//For the deformation vectors d1 d2 d3
			this->register_import(m_imDeformationVectord1);
			this->register_import(m_imDeformationVectord2);
			this->register_import(m_imDeformationVectord3);
			
			this->clear_add_fct();
		}

		////////////////////////////////////////////////////////////////////////////////
		//	register assemble functions
		////////////////////////////////////////////////////////////////////////////////
		#ifdef UG_DIM_1
		template<>
		void TransformedNavierStokes<Domain1d>::
		register_all_funcs(const LFEID& vLfeID, const LFEID& pLfeID, const int quadOrder)
		{
			UG_THROW("Not implemented.");
		}
		#endif

		#ifdef UG_DIM_2
		template<>
		void TransformedNavierStokes<Domain2d>::
		register_all_funcs(const LFEID& vLfeID, const LFEID& pLfeID, const int quadOrder)
		{
			typedef DimFEGeometry<dim> FVGeom;
			register_func<Triangle, FVGeom, FVGeom >();
			register_func<Quadrilateral, FVGeom, FVGeom >();
		}
		#endif

		#ifdef UG_DIM_3
		template<>
		void TransformedNavierStokes<Domain3d>::
		register_all_funcs(const LFEID& vLfeID, const LFEID& pLfeID, const int quadOrder)
		{
			typedef DimFEGeometry<dim> FVGeom;
			register_func<Tetrahedron, FVGeom, FVGeom >();
			register_func<Prism, FVGeom, FVGeom >();
			register_func<Hexahedron, FVGeom, FVGeom >();
		}
		#endif
		
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void TransformedNavierStokes<TDomain>::register_func()
		{
			ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
			typedef this_type T;

			this->clear_add_fct(id);
			this->set_prep_elem_loop_fct(	id, &T::template prep_elem_loop<TElem, VGeom, PGeom>);
			this->set_prep_elem_fct(	 	id, &T::template prep_elem<TElem, VGeom, PGeom>);
			this->set_fsh_elem_loop_fct( 	id, &T::template fsh_elem_loop<TElem, VGeom, PGeom>);
			this->set_add_jac_A_elem_fct(	id, &T::template add_jac_A_elem<TElem, VGeom, PGeom>);
			this->set_add_jac_M_elem_fct(	id, &T::template add_jac_M_elem<TElem, VGeom, PGeom>);
			this->set_add_def_A_elem_fct(	id, &T::template add_def_A_elem<TElem, VGeom, PGeom>);
			this->set_add_def_M_elem_fct(	id, &T::template add_def_M_elem<TElem, VGeom, PGeom>);
			this->set_add_rhs_elem_fct(	id, &T::template add_rhs_elem<TElem, VGeom, PGeom>);
		}
		
		////////////////////////////////////////////////////////////////////////////////
		//	explicit template instantiations
		////////////////////////////////////////////////////////////////////////////////

		#ifdef UG_DIM_2
		template class TransformedNavierStokes<Domain2d>;
		#endif
		#ifdef UG_DIM_3
		template class TransformedNavierStokes<Domain3d>;
		#endif
	 }//namespace TransformedNavierStokes
 }//namespace ug4