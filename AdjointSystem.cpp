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

#include "AdjointSystem.h"

namespace ug{
	namespace FluidOptim{
		////////////////////////////////////////////////////////////////////////////////
		//NON-REGISTERED FUNCTIONS: THEY USE IMPORTS TO CALCULATE, NOT PUBLIC
		////////////////////////////////////////////////////////////////////////////////
		
		//Calculate a matrix using the imports for velocity gradient
		template<typename TDomain>
		void AdjointSystem<TDomain>::
		VelocityGradient(MathMatrix<dim, dim>& GradV, const size_t ip)
		{
			MathVector<dim> velocity_d1=m_imVelocityGradientd1[ip];
			MathVector<dim> velocity_d2=m_imVelocityGradientd2[ip];
			
			GradV.assign(velocity_d1,0);
			GradV.assign(velocity_d2,1);
			if(this->dim == 3){
				MathVector<dim> velocity_d3=m_imVelocityGradientd3[ip];
				GradV.assign(velocity_d3,2);
			}
		}
		//Calculates the velocities at the ip and puts them into a vector 
		template<typename TDomain>
		void AdjointSystem<TDomain>::
		VelocityVector(MathVector<dim>& V, const size_t ip)
		{
			number v1= m_imVelocityd1[ip];
			number v2= m_imVelocityd2[ip];
			
			V[0]=v1;
			V[1]=v2;
			
			if(this->dim == 3){
				number v3= m_imVelocityd3[ip];			
				V[2]=v3;
			}
		}
		//Calculation of the Identity Deformation Gradient
		template<typename TDomain>
		void AdjointSystem<TDomain>::
		IdentityDeformationGradient(MathMatrix<dim, dim>& GradU, const size_t ip)
		{
			//1.- calculate gradients of each deformation component u1,u2
			MathVector<dim> vdeformation_d1=m_imDeformationVectord1[ip];
			MathVector<dim> vdeformation_d2=m_imDeformationVectord2[ip];
			//2.- assign the gradients to each row using assign(vector, row)
			GradU.assign(vdeformation_d1,0);
			GradU.assign(vdeformation_d2,1);
			
			if(this->dim == 3){
				MathVector<dim> vdeformation_d3 = m_imDeformationVectord3[ip];
				GradU.assign(vdeformation_d3,2);
			}
			//3.- add identity on the diagonal
			for(int i=0; i < dim; ++i)
			{
				GradU[i][i] += 1.0;
			}
		}
		
		////////////////////////////////////////////////////////////////////////////////
		//	IMPORT SETTING FUNCTIONS
		////////////////////////////////////////////////////////////////////////////////
		
		//Let's not use the LUA handlers for these imports, since we pass a CplUserData derived class directly :)
		
		//***************IMPORT: VELOCITY VECTOR COMPONENTS
		//d1
		template<typename TDomain>
		void AdjointSystem<TDomain>::
		set_velocity_d1(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imVelocityd1.set_data(user);
		}
		template<typename TDomain>
		void AdjointSystem<TDomain>::
		set_velocity_d1(number val)
		{
			set_velocity_d1(make_sp(new ConstUserNumber<dim>(val)));
		}		
		//d2
		template<typename TDomain>
		void AdjointSystem<TDomain>::
		set_velocity_d2(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imVelocityd2.set_data(user);
		}
		template<typename TDomain>
		void AdjointSystem<TDomain>::
		set_velocity_d2(number val)
		{
			set_velocity_d2(make_sp(new ConstUserNumber<dim>(val)));
		}
		//d3
		template<typename TDomain>
		void AdjointSystem<TDomain>::
		set_velocity_d3(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imVelocityd3.set_data(user);
		}
		template<typename TDomain>
		void AdjointSystem<TDomain>::
		set_velocity_d3(number val)
		{
			set_velocity_d3(make_sp(new ConstUserNumber<dim>(val)));
		}
		//***************IMPORT: VELOCITY GRADIENT ROWS
		//For the d1
		template<typename TDomain>
		void AdjointSystem<TDomain>::
		set_velocity_vector_d1(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imVelocityGradientd1.set_data(user);
		}
		template<typename TDomain>
		void AdjointSystem<TDomain>::
		set_velocity_vector_d1(number val)
		{
			set_velocity_vector_d1(make_sp(new ConstUserVector<dim>(val)));
		}
		//For the d2
		template<typename TDomain>
		void AdjointSystem<TDomain>::
		set_velocity_vector_d2(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imVelocityGradientd2.set_data(user);
		}
		template<typename TDomain>
		void AdjointSystem<TDomain>::
		set_velocity_vector_d2(number val)
		{
			set_velocity_vector_d2(make_sp(new ConstUserVector<dim>(val)));
		}
		//For the d3
		template<typename TDomain>
		void AdjointSystem<TDomain>::
		set_velocity_vector_d3(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imVelocityGradientd3.set_data(user);
		}
		template<typename TDomain>
		void AdjointSystem<TDomain>::
		set_velocity_vector_d3(number val)
		{
			set_velocity_vector_d3(make_sp(new ConstUserVector<dim>(val)));
		}
		//***************IMPORT: DEFORMATION GRADIENT ROWS
		//For the d1
		template<typename TDomain>
		void AdjointSystem<TDomain>::
		set_deformation_vector_d1(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imDeformationVectord1.set_data(user);
		}
		template<typename TDomain>
		void AdjointSystem<TDomain>::
		set_deformation_vector_d1(number val)
		{
			set_deformation_vector_d1(make_sp(new ConstUserVector<dim>(val)));
		}
		#ifdef UG_FOR_LUA
		template<typename TDomain>
		void AdjointSystem<TDomain>::
		set_deformation_vector_d1(const char* fctName)
		{
			set_deformation_vector_d1(LuaUserDataFactory<MathVector<dim>, dim>::create(fctName));
		}

		template<typename TDomain>
		void AdjointSystem<TDomain>::
		set_deformation_vector_d1(LuaFunctionHandle fct)
		{
			set_deformation_vector_d1(make_sp(new LuaUserData<MathVector<dim>, dim>(fct)));
		}
		#endif
		
		//For the d2 
		template<typename TDomain>
		void AdjointSystem<TDomain>::
		set_deformation_vector_d2(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imDeformationVectord2.set_data(user);
		}
		template<typename TDomain>
		void AdjointSystem<TDomain>::
		set_deformation_vector_d2(number val)
		{
			set_deformation_vector_d2(make_sp(new ConstUserVector<dim>(val)));
		}
		#ifdef UG_FOR_LUA
		template<typename TDomain>
		void AdjointSystem<TDomain>::
		set_deformation_vector_d2(const char* fctName)
		{
			set_deformation_vector_d2(LuaUserDataFactory<MathVector<dim>, dim>::create(fctName));
		}

		template<typename TDomain>
		void AdjointSystem<TDomain>::
		set_deformation_vector_d2(LuaFunctionHandle fct)
		{
			set_deformation_vector_d2(make_sp(new LuaUserData<MathVector<dim>, dim>(fct)));
		}
		#endif
		//For the d3
		template<typename TDomain>
		void AdjointSystem<TDomain>::
		set_deformation_vector_d3(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imDeformationVectord3.set_data(user);
		}
		template<typename TDomain>
		void AdjointSystem<TDomain>::
		set_deformation_vector_d3(number val)
		{
			set_deformation_vector_d3(make_sp(new ConstUserVector<dim>(val)));
		}
		////////////////////////////////////////////////////////////////////////////////
		//	general
		////////////////////////////////////////////////////////////////////////////////
		template<typename TDomain>
		void AdjointSystem<TDomain>::
		set_quad_order(size_t order)
		{
			m_quadOrder = order;
			m_bQuadOrderUserDef = true;
		}
		
		template<typename TDomain>
		void AdjointSystem<TDomain>::
		prepare_setting(const std::
		vector<LFEID>& vLfeID, bool bNonRegularGrid)
		{
			if(bNonRegularGrid)
				UG_THROW("AdjointSystem: only implemented for regular grids.");

			//	check number
			if(vLfeID.size() != dim+1)
				UG_THROW("AdjointSystem: Needs exactly "<<dim+1<<" functions.");

			for(int d = 1; d < dim; ++d)
				if(vLfeID[0] != vLfeID[d])
					UG_THROW("AdjointSystem: trial spaces for velocity expected to be"
					" identical for all velocity components.");

			//	remember lfeID;
			m_vLFEID = vLfeID[0];
			m_pLFEID = vLfeID[dim];

			if(!m_bQuadOrderUserDef) m_quadOrder = 2*m_vLFEID.order()+1;

			//	update assemble functions
			register_all_funcs(m_vLFEID, m_pLFEID, m_quadOrder);
		}
		
		template<typename TDomain>
		void AdjointSystem<TDomain>::
		set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imKinViscosity.set_data(user);
		}
		template<typename TDomain>
		void AdjointSystem<TDomain>::
		set_kinematic_viscosity(number val)
		{
			if(val == 0.0) set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> >());
			else set_kinematic_viscosity(make_sp(new ConstUserNumber<dim>(val)));
		}
		#ifdef UG_FOR_LUA
		template<typename TDomain>
		void AdjointSystem<TDomain>::
		set_kinematic_viscosity(const char* fctName)
		{
			set_kinematic_viscosity(LuaUserDataFactory<number,dim>::create(fctName));
		}

		template<typename TDomain>
		void AdjointSystem<TDomain>::
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
		void AdjointSystem<TDomain>::prep_elem_loop(const ReferenceObjectID roid, const int si)
		{
			DimFEGeometry<dim>& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);
			DimFEGeometry<dim>& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_quadOrder);
			try
			{
				vgeo.update_local(roid, m_vLFEID, m_quadOrder);
				pgeo.update_local(roid, m_pLFEID, m_quadOrder);
			}
			UG_CATCH_THROW("FluidOptim::AdjointSystem: "
			"Cannot update Finite Element Geometry.");
			m_imKinViscosity.set_local_ips(vgeo.local_ips(), vgeo.num_ip(),false);
			
			m_imDeformationVectord1.set_local_ips(vgeo.local_ips(), vgeo.num_ip(),false);
			m_imDeformationVectord2.set_local_ips(vgeo.local_ips(), vgeo.num_ip(),false);
			m_imDeformationVectord3.set_local_ips(vgeo.local_ips(), vgeo.num_ip(),false);
			
			m_imVelocityGradientd1.set_local_ips(vgeo.local_ips(), vgeo.num_ip(),false);
			m_imVelocityGradientd2.set_local_ips(vgeo.local_ips(), vgeo.num_ip(),false);
			m_imVelocityGradientd3.set_local_ips(vgeo.local_ips(), vgeo.num_ip(),false);
			
			m_imVelocityd1.set_local_ips(vgeo.local_ips(), vgeo.num_ip(),false);
			m_imVelocityd2.set_local_ips(vgeo.local_ips(), vgeo.num_ip(),false);
			m_imVelocityd3.set_local_ips(vgeo.local_ips(), vgeo.num_ip(),false);
		}
		
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void AdjointSystem<TDomain>::prep_elem(const LocalVector& u, GridObject* elem, 
											   const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
		{
			m_pElem = elem;

			// 	Update Geometry for this element
			DimFEGeometry<dim>& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);
			DimFEGeometry<dim>& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_quadOrder);
			try{
				vgeo.update(elem, vCornerCoords, m_vLFEID, m_quadOrder);
				pgeo.update(elem, vCornerCoords, m_pLFEID, m_quadOrder);
			}
			UG_CATCH_THROW("AdjointSystem: Cannot update "
							"Finite Element Geometry.");

			m_imKinViscosity.set_global_ips(vgeo.global_ips(), vgeo.num_ip());
			
			m_imDeformationVectord1.set_global_ips(vgeo.global_ips(), vgeo.num_ip());
			m_imDeformationVectord2.set_global_ips(vgeo.global_ips(), vgeo.num_ip());
			m_imDeformationVectord3.set_global_ips(vgeo.global_ips(), vgeo.num_ip());
			
			m_imVelocityGradientd1.set_global_ips(vgeo.global_ips(), vgeo.num_ip());
			m_imVelocityGradientd2.set_global_ips(vgeo.global_ips(), vgeo.num_ip());
			m_imVelocityGradientd3.set_global_ips(vgeo.global_ips(), vgeo.num_ip());
			
			m_imVelocityd1.set_global_ips(vgeo.global_ips(), vgeo.num_ip());
			m_imVelocityd2.set_global_ips(vgeo.global_ips(), vgeo.num_ip());
			m_imVelocityd3.set_global_ips(vgeo.global_ips(), vgeo.num_ip());
			
		}
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void AdjointSystem<TDomain>::fsh_elem_loop()
		{
		}
		///	assembles the local stiffness matrix 
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void AdjointSystem<TDomain>::
		add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{
			//	request geometry
			const DimFEGeometry<dim>& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);
			const DimFEGeometry<dim>& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_quadOrder);
													
			for (size_t ip = 0; ip < vgeo.num_ip(); ++ip){
				//TODO: think about declaring here or outside loop? What is best for parallelism?
				MathMatrix<dim, dim> UnitDeformationGradient;MatSet(UnitDeformationGradient,0.0);
				MathMatrix<dim, dim> InvUnitDeformationGradient;MatSet(InvUnitDeformationGradient,0.0);
				number dDF=0.0;
				
				IdentityDeformationGradient(UnitDeformationGradient, ip);
				Inverse(InvUnitDeformationGradient,UnitDeformationGradient);
				dDF=Determinant(UnitDeformationGradient);
				
				////////////////////////////////////////////////////
				// Stiffness/Diffusion Matrix Terms
				////////////////////////////////////////////////////
				const number scale = m_imKinViscosity[ip]* vgeo.weight(ip);
				
				for (int vdim = 0; vdim < dim; ++vdim){
					for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){
						for (size_t ush = 0; ush < vgeo.num_sh(); ++ush){
										
							MathVector<dim> v_TransformedGradient1; VecSet(v_TransformedGradient1,0.0);								
							MathVector<dim> v_TransformedGradient2; VecSet(v_TransformedGradient2,0.0);

							TransposedMatVecMult(v_TransformedGradient1,InvUnitDeformationGradient,vgeo.global_grad(ip,ush));
							TransposedMatVecMult(v_TransformedGradient2,InvUnitDeformationGradient,vgeo.global_grad(ip,vsh));
							for (int udim = 0; udim < dim; ++udim) {
								J(vdim, vsh, vdim, ush) -=  scale
															* v_TransformedGradient1[udim]
															* v_TransformedGradient2[udim]*dDF;

							}
						}
					}
				}
				
				/////////////////////////////////////////////////////////////
				// Newton-Derivative Matrix, uses Gradient of Velocity at ip 
				/////////////////////////////////////////////////////////////
				// 	Interpolate Functional Matrix of velocity at ip
				MathMatrix<dim, dim> GradVel;MatSet(GradVel,0.0);
				VelocityGradient(GradVel,ip);
																	
				MathMatrix<dim,dim> TransformedGradVel; MatSet(TransformedGradVel, 0.0);
				MatMultiply(TransformedGradVel,GradVel,InvUnitDeformationGradient);
				
				for(int vdim=0 ; vdim < dim; vdim++){
					for(size_t vsh=0 ; vsh < vgeo.num_sh() ; vsh++){
						for(size_t ush=0 ; ush<vgeo.num_sh() ; ush++){
							for(int udim=0 ; udim < dim ; udim++){
								J(vdim,vsh,udim,ush) -= TransformedGradVel(udim,vdim)
											 			*vgeo.shape(ip,ush)
														*vgeo.shape(ip,vsh)
														*vgeo.weight(ip)
														*dDF; 
							}
						}
					}

				}//end sh1
				////////////////////////////////////////////////////
				// Vector-Convection Matrix, uses Velocity value at ip
				////////////////////////////////////////////////////
				
				// 	Interpolate Velocity at ip
				MathVector<dim> Vel;VecSet(Vel,0.0);
				VelocityVector(Vel,ip);
				for(int vdim=0; vdim < dim ; vdim++){
					for(size_t vsh=0; vsh < vgeo.num_sh() ;++vsh){
						for(size_t ush=0; ush < vgeo.num_sh(); ush++){
							
							
							MathMatrix<dim,dim> GradientQ;MatSet(GradientQ,0.0);
							GradientQ.assign(vgeo.global_grad(ip,vsh),vdim);
							MathMatrix<dim,dim> TransformedGradientQ;MatSet(TransformedGradientQ,0.0);
							MatMultiply(TransformedGradientQ,GradientQ,InvUnitDeformationGradient);
							//Perform GradQInvDF*V
							MathVector<dim> vTemp1;VecSet(vTemp1,0.0);//MatVec multiplication
							MatVecMult(vTemp1,TransformedGradientQ,Vel);
							//Perform dot product inside loop, 
							for(int udim=0; udim < dim; udim++){
								MathVector<dim> vShapes;VecSet(vShapes,0.0);vShapes[udim]=vgeo.shape(ip,ush);
								number vec_prod = VecProd(vShapes,vTemp1);
								
								J(vdim,vsh,udim,ush) -= vec_prod*vgeo.weight(ip)*dDF;
							} 
 							
							

						}
					}
				}//end sh1
				
				
				////////////////////////////////////////////////////
				// Adjoint Pressure Term H
				////////////////////////////////////////////////////
				
				for(int vdim=0 ; vdim < dim ; ++vdim){
					for(size_t vsh=0;vsh<vgeo.num_sh();vsh++){
						for(size_t psh=0 ; psh < pgeo.num_sh() ; ++psh){
										
							
							//Trace(GradientQ*InvDF)
							MathMatrix<dim, dim> GradientAdjoint;MatSet(GradientAdjoint,0.0);
							GradientAdjoint.assign(vgeo.global_grad(ip,vsh),vdim);
							MathMatrix<dim, dim> TransformedAdjointVelGradient;MatSet(TransformedAdjointVelGradient,0.0);
							MatMultiply(TransformedAdjointVelGradient,GradientAdjoint,InvUnitDeformationGradient);
							//trace calculation
							number TrTransAdjVel = Trace(TransformedAdjointVelGradient);
							//assign to matrix
							J(vdim,vsh,_P_,psh) += TrTransAdjVel*pgeo.shape(ip,psh)*vgeo.weight(ip)*dDF;

							
						}//psh
					}//vdim		
				}//vsh
				
				////////////////////////////////////////////////////
				// Continuity Equation
				////////////////////////////////////////////////////
				//calculate the trace before assigning value
				for(size_t psh=0; psh < pgeo.num_sh() ; ++psh){
					for(size_t vsh=0; vsh < vgeo.num_sh() ; ++vsh){
						for(size_t vdim=0; vdim < (size_t) dim; ++vdim){
								
							//Trace(GradientQ*InvDF)
							MathMatrix<dim, dim> GradientAdjoint;MatSet(GradientAdjoint,0.0);
							GradientAdjoint.assign(vgeo.global_grad(ip,vsh),vdim);
							MathMatrix<dim, dim> TransformedAdjointVelGradient;MatSet(TransformedAdjointVelGradient,0.0);
							MatMultiply(TransformedAdjointVelGradient,GradientAdjoint,InvUnitDeformationGradient);
							number TrTransAdjVel = Trace(TransformedAdjointVelGradient);
							
							//assign value to matrix
							J(_P_,psh,vdim,vsh) += TrTransAdjVel*pgeo.shape(ip,psh)*vgeo.weight(ip)*dDF;

						}
					}
				}
			}//end ip
			
			////////////////////////////////////////////////////
			// Stabilization Term
			////////////////////////////////////////////////////
			if(m_stab_type == 0.0 && m_stabParam != 0.0)
			{
				const number scale = m_stabParam * ElementDiameterSq<GridObject, TDomain>(*m_pElem, *this->domain());
				for(size_t ip=0; ip<pgeo.num_ip();++ip){
					//TODO:verify that approach is correct, jacobian
					MathMatrix<dim, dim> UnitDeformationGradient;
					MathMatrix<dim, dim> InvUnitDeformationGradient;
					number dDF=0.0;
									
					IdentityDeformationGradient(UnitDeformationGradient, ip);
					Inverse(InvUnitDeformationGradient,UnitDeformationGradient);
					dDF=Determinant(UnitDeformationGradient);
									
					MathVector<dim> v_tempVec1;VecSet(v_tempVec1,0.0);
					MathVector<dim> v_tempVec2;VecSet(v_tempVec2,0.0);
					for(size_t psh1=0;psh1<pgeo.num_sh();++psh1){
						TransposedMatVecMult(v_tempVec1,InvUnitDeformationGradient,pgeo.global_grad(ip,psh1));
										
						for(size_t psh2=0;psh2<pgeo.num_sh();psh2++){
							TransposedMatVecMult(v_tempVec2,InvUnitDeformationGradient,pgeo.global_grad(ip,psh2));
							J(_P_,psh1,_P_,psh2)+= scale*VecDot(v_tempVec1,v_tempVec2)*pgeo.weight(ip)*dDF;
						}
					}
				}//end pressure ip		
			}//end stab param !=0 
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
							J(_P_,psh1,_P_,psh2)+= scale*(pgeo.shape(ip,psh1)*pgeo.shape(ip,psh2)-1.0/(pgeo.num_sh()*pgeo.num_sh()))*pgeo.weight(ip)*dDF;						}
					}
				}						
			}//case pressure_projection
		}
		///	assembles the local stiffness matrix 
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void AdjointSystem<TDomain>::
		add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{
			const DimFEGeometry<dim>& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);
			const DimFEGeometry<dim>& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_quadOrder);
													
			for (size_t ip = 0; ip < vgeo.num_ip(); ++ip){
				//1.-Need Gradient of the velocity field (matrix) and the Deformation Gradient
				MathMatrix<dim, dim> VelGrad;MatSet(VelGrad,0.0);
				MathMatrix<dim, dim> UnitDeformationGradient;MatSet(UnitDeformationGradient,0.0);
				MathMatrix<dim, dim> InvUnitDeformationGradient;MatSet(InvUnitDeformationGradient,0.0);
				number dDF=0.0;
				
				VelocityGradient(VelGrad,ip);
				IdentityDeformationGradient(UnitDeformationGradient, ip);
				Inverse(InvUnitDeformationGradient,UnitDeformationGradient);
				dDF=Determinant(UnitDeformationGradient);
				
				const number scale= m_imKinViscosity[ip]*vgeo.weight(ip)*dDF;
				
				//2.- We need to transform the VelGrad
				MathMatrix<dim, dim> TransformedGradVel;MatSet(TransformedGradVel,0.0);				
				MatMultiply(TransformedGradVel,VelGrad,InvUnitDeformationGradient);
				
				//3.- Assign value to d
				for (int vdim = 0; vdim < dim; ++vdim){
					for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){
						
						//4.- Transformation of shape function gradient at shape function
						/* MathVector<dim> v_TransformedGradient; VecSet(v_TransformedGradient,0.0);
						MathMatrix<dim, dim> TranposedInverse; Transpose( TranposedInverse , InvUnitDeformationGradient );
						MatVecMult(v_TransformedGradient,TranposedInverse,vgeo.global_grad(ip,vsh)); */
						
						//different approach now
						MathMatrix<dim, dim> GradientAdjoint;MatSet(GradientAdjoint,0.0);
						GradientAdjoint.assign(vgeo.global_grad(ip,vsh),vdim);
						MathMatrix<dim, dim> TransformedAdjointVelGradient;MatSet(TransformedAdjointVelGradient,0.0);
						MatMultiply(TransformedAdjointVelGradient,GradientAdjoint,InvUnitDeformationGradient);
						number contraction_VelAdjoint = MatContraction(TransformedGradVel,TransformedAdjointVelGradient);
						//assign to rhs
						d(vdim, vsh) -= scale*contraction_VelAdjoint;
						
						/* for (int udim = 0; udim < dim; ++udim){
							d(vdim, vsh) +=  scale*TransformedGradVel(vdim, udim)*TransformedGradient[udim];
						} */
					}
				}
				
			}
		}
		///	assembles the local right hand side
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void AdjointSystem<TDomain>::
		add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{	}
		//Mass Matrix methods will be assembled if there is time dependency
		///	assembles the local mass matrix using a finite volume scheme
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void AdjointSystem<TDomain>::
		add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{	}
		///	assembles the mass part of the local defect
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void AdjointSystem<TDomain>::
		add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{	}
		
		////////////////////////////////////////////////////////////////////////////////
		//	Constructor
		////////////////////////////////////////////////////////////////////////////////

		template<typename TDomain>
		AdjointSystem<TDomain>::
		AdjointSystem(const char* functions, const char* subsets)
		:IElemDisc<TDomain>(functions,subsets)
		{
		//	set defaults
			m_bQuadOrderUserDef = false;
			m_quadOrder = -1;
			m_stab_type=0.0;
		//	check number of functions
			if(this->num_fct() != dim+1)
				UG_THROW("Wrong number of functions: The ElemDisc 'AdjointSystem'"
					   " needs exactly "<<dim+1<<" symbolic function.");
			//Default value assigned	   
			m_stabParam=0.0;// no stab
		//	register imports
			this->register_import(m_imKinViscosity);			
			//For the deformation vectors d1 d2 d3
			this->register_import(m_imDeformationVectord1);
			this->register_import(m_imDeformationVectord2);
			this->register_import(m_imDeformationVectord3);
			
			this->register_import(m_imVelocityGradientd1);
			this->register_import(m_imVelocityGradientd2);
			this->register_import(m_imVelocityGradientd3);
			
			this->register_import(m_imVelocityd1);
			this->register_import(m_imVelocityd2);
			this->register_import(m_imVelocityd3);
			
			this->clear_add_fct();
		}
		
		template<typename TDomain>
		AdjointSystem<TDomain>::
		AdjointSystem(const std::vector<std::string>& vFct,
                      const std::vector<std::string>& vSubset)
		:IElemDisc<TDomain>(vFct,vSubset)
		{
		//	set defaults		
			m_bQuadOrderUserDef = false;
			m_quadOrder = -1;
		//	check number of functions
			if(this->num_fct() != dim+1)
				UG_THROW("Wrong number of functions: The ElemDisc 'AdjointSystem'"
					   " needs exactly "<<dim+1<<" symbolic function.");
			//Default value assigned	   
			m_stabParam=0.0;// no stab
			m_stab_type=0.0;
		//	register imports
			this->register_import(m_imKinViscosity);
			
			//For the deformation vectors d1 d2 d3
			this->register_import(m_imDeformationVectord1);
			this->register_import(m_imDeformationVectord2);
			this->register_import(m_imDeformationVectord3);
			
			this->register_import(m_imVelocityGradientd1);
			this->register_import(m_imVelocityGradientd2);
			this->register_import(m_imVelocityGradientd3);
			
			this->register_import(m_imVelocityd1);
			this->register_import(m_imVelocityd2);
			this->register_import(m_imVelocityd3);
			
			this->clear_add_fct();
		}
		
		////////////////////////////////////////////////////////////////////////////////
		//	register assemble functions
		////////////////////////////////////////////////////////////////////////////////
		#ifdef UG_DIM_1
		template<>
		void AdjointSystem<Domain1d>::
		register_all_funcs(const LFEID& vLfeID, const LFEID& pLfeID, const int quadOrder)
		{
			UG_THROW("Not implemented.");
		}
		#endif

		#ifdef UG_DIM_2
		template<>
		void AdjointSystem<Domain2d>::
		register_all_funcs(const LFEID& vLfeID, const LFEID& pLfeID, const int quadOrder)
		{
			typedef DimFEGeometry<dim> FVGeom;
			register_func<Triangle, FVGeom, FVGeom >();
			register_func<Quadrilateral, FVGeom, FVGeom >();
		}
		#endif

		#ifdef UG_DIM_3
		template<>
		void AdjointSystem<Domain3d>::
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
		void AdjointSystem<TDomain>::register_func()
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
		template class AdjointSystem<Domain2d>;
		#endif
		#ifdef UG_DIM_3
		template class AdjointSystem<Domain3d>;
		#endif
	}//end namespace FluidOpt
}//end namespace ug