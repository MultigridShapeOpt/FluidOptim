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

#include "DeformationStateEquation.h"

// for various user data
#include "bindings/lua/lua_user_data.h"
#include "lib_disc/spatial_disc/user_data/user_data.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"
#include "lib_disc/spatial_disc/disc_util/fe_geom.h"

using namespace std;
namespace ug{
namespace FluidOptim{
////////////////////////////////////////////////////////////////////////////////
//NON-REGISTERED FUNCTIONS: THEY USE IMPORTS TO CALCULATE, NOT PUBLIC
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//	General
////////////////////////////////////////////////////////////////////////////////
template<typename TDomain>
void DeformationStateEquation<TDomain>::
set_nonlinear(bool nl)
{
	this->m_bNonLinear=nl;
}

template<typename TDomain>
void DeformationStateEquation<TDomain>::
set_picard(bool pic)
{
	this->m_bPicardIteration=pic;
}

template<typename TDomain>
void DeformationStateEquation<TDomain>::
set_symmetric(bool sym)
{
	this->m_bSymmetric=sym;
}

/* template<typename TDomain>
void DeformationStateEquation<TDomain>::
set_nu_extension(number n)
{
	this->m_nu_extension=n;
} */
template<typename TDomain>
void DeformationStateEquation<TDomain>::
set_extension_factor(SmartPtr<CplUserData<number,dim> > user_data)
{
	this->m_imExtensionFactor.set_data(user_data);
}
		
template<typename TDomain>
void DeformationStateEquation<TDomain>::
set_extension_factor(number control)
{
	this->set_extension_factor(make_sp(new ConstUserNumber<dim>(control)));
}
template<typename TDomain>
void DeformationStateEquation<TDomain>::
update_geo_elem(TBaseElem* elem, DimFEGeometry<dim>& geo)
{
	SmartPtr<TDomain> dom = this->domain();

	typedef typename IElemDisc<TDomain>::domain_type::position_accessor_type
			position_accessor_type;
	const position_accessor_type& aaPos = dom->position_accessor();

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];

	//	get vertices and extract corner coordinates
	const size_t numVertices = elem->num_vertices();
	for (size_t i = 0; i < numVertices; ++i) {
		coCoord[i] = aaPos[elem->vertex(i)];
	};

	//	prepare geometry for type and order
   	try{
		geo.update(elem, &(coCoord[0]), m_lfeID, m_quadOrder);
	}
	UG_CATCH_THROW("DeformationStateEquation::update_geo_elem:"
					" Cannot update Finite Element Geometry.");
}
	
template<typename TDomain>
template<typename TElem, typename TFEGeom>
void DeformationStateEquation<TDomain>::
prep_timestep_elem(const number time, const LocalVector& u,
		GridObject* elem, const MathVector<dim> vCornerCoords[])
{
}

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void DeformationStateEquation<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{

	//	request geometry
	TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

	//	prepare geometry for type and order
	try{
		geo.update_local(roid, m_lfeID, m_quadOrder);
	}UG_CATCH_THROW("DeformationStateEquation::prep_elem_loop:"
					" Cannot update Finite Element Geometry.");
					
	//ExtensionFactor
	m_imExtensionFactor.template set_local_ips(geo.local_ips(),geo.num_ip(),false);
}

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void DeformationStateEquation<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
{
	//	request geometry
	TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

	try{
		geo.update(elem, vCornerCoords, m_lfeID, m_quadOrder);
	}
	UG_CATCH_THROW("DeformationStateEquation::prep_elem:"
					" Cannot update Finite Element Geometry.");
	//ExtensionFactor
	m_imExtensionFactor.set_global_ips(geo.global_ips(),geo.num_ip());

}

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void DeformationStateEquation<TDomain>::fsh_elem_loop()
{}

//  assemble stiffness jacobian
template<typename TDomain>
template<typename TElem, typename TFEGeom>

void DeformationStateEquation<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u,
		GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);


	for (size_t ip = 0; ip < geo.num_ip(); ++ip){
		const number K = m_imExtensionFactor[ip];
		for(int d1=0; d1 <  dim; ++d1){
			for(size_t sh1=0;sh1<geo.num_sh();++sh1){
				for(size_t sh2=0;sh2<geo.num_sh();++sh2){
					for(int d2=0 ; d2< dim ;++d2){
					/* 	//Laplacian Matrix
						J(d1,sh1,d1,sh2) += geo.global_grad(ip,sh2)[d2]*geo.global_grad(ip,sh1)[d2]*geo.weight(ip);
						
						if(m_bSymmetric){							
							//J(d2,sh2,d1,sh1) += geo.global_grad(ip,sh1)[d1]*geo.global_grad(ip,sh2)[d2]*geo.weight(ip);
							J(d1,sh2,d2,sh1) += geo.global_grad(ip,sh1)[d1]*geo.global_grad(ip,sh2)[d2]*geo.weight(ip);							
						}  */
						
						MathMatrix<dim,dim> M1;MatSet(M1,0.0);
						MathMatrix<dim,dim> M2;MatSet(M2,0.0);//Transpose M1
						MathMatrix<dim,dim> M3;MatSet(M3,0.0);
						MathMatrix<dim,dim> MatSum;MatSet(MatSum,0.0);
						M1.assign(geo.global_grad(ip,sh1),d1);
						Transpose(M2,M1);
						MatAdd(MatSum,M1,M2);
						M3.assign(geo.global_grad(ip,sh2),d2);
						
						if(m_bSymmetric)
							J(d1,sh1,d2,sh2) += MatContraction(MatSum,M3)*geo.weight(ip);
						else
							J(d1,sh1,d2,sh2) += MatContraction(M1,M3)*geo.weight(ip);
						
					}//d2
				}//end sh2
			}//end sh1
		}//d1
		
		
		if(m_bNonLinear)
		{
			MathVector<dim> Def;
			for(int d1 = 0; d1 < dim; ++d1){
				Def[d1] = 0.0;
				for(size_t sh = 0; sh < geo.num_sh(); ++sh){
					Def[d1] += u(d1, sh)*geo.shape(ip, sh);
				}	
			}
			for(size_t d1=0;d1<(size_t) dim;++d1){
				for(size_t sh1=0;sh1<geo.num_sh();++sh1){
					for(size_t sh2=0;sh2<geo.num_sh();++sh2){
						for(size_t d2=0;d2<(size_t) dim;++d2){
							J(d1,sh1,d1,sh2) += K
												*Def[d2]
												*geo.global_grad(ip,sh2)[d2]
												*geo.shape(ip,sh1)
												*geo.weight(ip);
						}//d2
					}//end sh2
				}//end sh1
			}//d1
			if(!m_bPicardIteration)
			{
				MathMatrix<dim, dim> GradU;
				DisplacementGradient<TFEGeom>(GradU, ip, geo, u);
				
				for(size_t d1=0;d1<(size_t) dim;++d1){
					for(size_t sh1=0;sh1<geo.num_sh();++sh1){
						for(size_t sh2=0;sh2<geo.num_sh();++sh2){
							for(size_t d2=0;d2<(size_t) dim;++d2){
								J(d1,sh1,d2,sh2) += K
													*GradU(d1,d2)
													*geo.shape(ip,sh2)
													*geo.shape(ip,sh1)
													*geo.weight(ip);
							}//d2
						}//end sh2
					}//end sh1
				}//d1
			}//end Picard
		}//end if nonlinear
		

	} //end(ip)

}

//  assemble stiffness defect
template<typename TDomain>
template<typename TElem, typename TFEGeom>
void DeformationStateEquation<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u,
		GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

	MathMatrix<dim, dim> GradU;
	MathMatrix<dim,dim> TransposeGradU;
	MathMatrix<dim,dim> SymGradU;
	for (size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		
		const number K = m_imExtensionFactor[ip];
	 	DisplacementGradient<TFEGeom>(GradU, ip, geo, u);
		if(m_bSymmetric)
		{
			Transpose(TransposeGradU,GradU);
			MatAdd(SymGradU,GradU,TransposeGradU);
			if(SymGradU(0,1) != SymGradU(1,0))
				UG_THROW("DeformationStateEquation: the Deformation Gradient was not Correctly Symmetrized");
		}
		for(size_t d1 = 0; d1 < (size_t) TDomain::dim; ++d1){
			for (size_t sh = 0; sh < geo.num_sh(); ++sh){ // loop shape functions
				for (size_t d2 = 0; d2 < (size_t) TDomain::dim; ++d2){ // loop components
						if(!m_bSymmetric)
							d(d1, sh) +=   GradU(d1,d2)* geo.global_grad(ip, sh)[d2]*geo.weight(ip);
						if(m_bSymmetric)
							d(d1, sh) +=   SymGradU(d1,d2)* geo.global_grad(ip, sh)[d2]*geo.weight(ip);
				} 
			}
		}
		
		
		if(m_bNonLinear)
		{	
			////////////////////////////////////////////////////
			// Convective Term 
			////////////////////////////////////////////////////
			//Interpolate variable at ip
			MathVector<dim> Def;
			for(int d1 = 0; d1 < dim; ++d1){
				Def[d1] = 0.0;
				for(size_t sh = 0; sh < geo.num_sh(); ++sh)
					Def[d1] += u(d1, sh)*geo.shape(ip, sh);
			}
			
			MathVector<dim> ConvectiveDeformation;
			MatVecMult(ConvectiveDeformation, GradU, Def);
			
			for (int d1 = 0; d1 < dim; ++d1){
				for (size_t sh = 0; sh < geo.num_sh(); ++sh){
					//d(d1, sh) +=  m_nu_extension*ConvectiveDeformation[d1]*geo.shape(ip, sh)*geo.weight(ip);
					d(d1, sh) +=  K*ConvectiveDeformation[d1]*geo.shape(ip, sh)*geo.weight(ip);
				}
			}
		}//end if nonlinear
	}//end (ip)

}


//	assemble mass jacobian
template<typename TDomain>
template<typename TElem, typename TFEGeom>
void DeformationStateEquation<TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u,
		GridObject* elem, const MathVector<dim> vCornerCoords[])
{
}
//  assemble mass-defect
template<typename TDomain>
template<typename TElem, typename TFEGeom>
void DeformationStateEquation<TDomain>::
add_def_M_elem(LocalVector& d, const LocalVector& u,
		GridObject* elem, const MathVector<dim> vCornerCoords[])
{

}

//  assemble right-hand-side d(i,sh)
template<typename TDomain>
template<typename TElem, typename TFEGeom>
void DeformationStateEquation<TDomain>::
add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
{

}
template<typename TDomain>
template<typename TElem, typename TFEGeom>
void DeformationStateEquation<TDomain>::
fsh_timestep_elem(const number time, const LocalVector& u,
		GridObject* elem, const MathVector<dim> vCornerCoords[])
{

}



template <typename TDomain>
void
DeformationStateEquation<TDomain>::
set_assemble_funcs()
{
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

	register_all_fe_funcs(m_order, m_quadOrder);

}

////////////////////////////////////////////////////////////////////////////////
//	Constructor
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain>
DeformationStateEquation<TDomain>::
DeformationStateEquation(const char* functions, const char* subsets) :
			IElemDisc<TDomain> (functions, subsets)
			
			
{
	//	check number of functions
	if (this->num_fct() != (size_t) dim)
		UG_THROW("Wrong number of functions: The ElemDisc 'DeformationStateEquation'"
				" needs exactly "<<dim<<" symbolic function.");
	//	set defaults
	m_order = 1;
	m_bQuadOrderUserDef = false;
	m_quadOrder = -1;
	
	m_bNonLinear=false;
	m_bPicardIteration=false;
	m_bSymmetric=false;
	//m_nu_extension=1;
	//	update assemble functions
	set_assemble_funcs();
	//ExtensionFactor
	this->register_import(m_imExtensionFactor);
}


template <typename TDomain>
DeformationStateEquation<TDomain>::
~DeformationStateEquation()
{

}


////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////
#ifdef UG_DIM_1
template<>
void DeformationStateEquation<Domain1d>::register_all_fe_funcs(int order,
		int quadOrder)
{
	//	RegularEdge
	UG_THROW("Not implemented.");
}
#endif
#ifdef UG_DIM_2
template<>
void DeformationStateEquation<Domain2d>::register_all_fe_funcs(int order,
		int quadOrder)
{
	
	register_fe_func<Triangle, DimFEGeometry<dim, 2> >();
	register_fe_func<Quadrilateral, DimFEGeometry<dim, 2> >();
/* 	if (quadOrder != 2 * order + 1) {
		register_fe_func<Triangle, DimFEGeometry<dim> > ();
		register_fe_func<Quadrilateral, DimFEGeometry<dim> > ();
	} */

/* 	//	special compiled cases

	//	Triangle
	switch (order)
	{
		case 1:{
			typedef FEGeometry<Triangle, dim, LagrangeLSFS<ReferenceTriangle, 1> ,
					GaussQuadrature<ReferenceTriangle, 3> > FEGeom;
			register_fe_func<Triangle, FEGeom > (); break;
		}
		case 2:{
			typedef FEGeometry<Triangle, dim, LagrangeLSFS<ReferenceTriangle, 2> ,
					GaussQuadrature<ReferenceTriangle, 5> > FEGeom;
			register_fe_func<Triangle, FEGeom > (); break;
		}
		case 3:{
			typedef FEGeometry<Triangle, dim, LagrangeLSFS<ReferenceTriangle, 3> ,
					GaussQuadrature<ReferenceTriangle, 7> > FEGeom;
			register_fe_func<Triangle, FEGeom > (); break;
		}
		default: register_fe_func<Triangle, DimFEGeometry<dim> > (); break;
	}

	//	Quadrilateral
	switch (order)
	{
		case 1:{
			typedef FEGeometry<Quadrilateral, dim, LagrangeLSFS<
					ReferenceQuadrilateral, 1> , GaussQuadrature<
					ReferenceQuadrilateral, 3> > FEGeom;
			register_fe_func<Quadrilateral, FEGeom > (); break;
		}
		case 2:{
			typedef FEGeometry<Quadrilateral, dim, LagrangeLSFS<
					ReferenceQuadrilateral, 2> , GaussQuadrature<
					ReferenceQuadrilateral, 7> > FEGeom;
			register_fe_func<Quadrilateral, FEGeom > (); break;
		}
		case 3:{
			typedef FEGeometry<Quadrilateral, dim, LagrangeLSFS<
					ReferenceQuadrilateral, 3> , GaussQuadrature<
					ReferenceQuadrilateral, 11> > FEGeom;
			register_fe_func<Quadrilateral, FEGeom > (); break;
		}
		default: register_fe_func<Quadrilateral, DimFEGeometry<dim> > (); break;
	} */
}
#endif

#ifdef UG_DIM_3
template<>
void DeformationStateEquation<Domain3d>::register_all_fe_funcs(int order,
		int quadOrder)
{
	if (quadOrder != 2 * order + 1) {
		register_fe_func<Tetrahedron, DimFEGeometry<dim> > ();
		register_fe_func<Prism, DimFEGeometry<dim> > ();
		register_fe_func<Pyramid, DimFEGeometry<dim> > ();
		register_fe_func<Hexahedron, DimFEGeometry<dim> > ();
	}

	//	special compiled cases

	//	Tetrahedron
	switch (order) 
	{
		case 1:{
			typedef FEGeometry<Tetrahedron, dim, LagrangeLSFS<
				ReferenceTetrahedron, 1> , GaussQuadrature<
				ReferenceTetrahedron, 3> > FEGeom;
			register_fe_func<Tetrahedron, FEGeom > (); break;
		}
		case 2:{
			typedef FEGeometry<Tetrahedron, dim, LagrangeLSFS<
					ReferenceTetrahedron, 2> , GaussQuadrature<
					ReferenceTetrahedron, 5> > FEGeom;
			register_fe_func<Tetrahedron, FEGeom > (); break;
		}
		case 3:{
			typedef FEGeometry<Tetrahedron, dim, LagrangeLSFS<
					ReferenceTetrahedron, 3> , GaussQuadrature<
					ReferenceTetrahedron, 7> > FEGeom;
			register_fe_func<Tetrahedron, FEGeom > (); break;
		}
		default: register_fe_func<Tetrahedron, DimFEGeometry<dim> > (); break;

	}

	//	Prism
	switch (order)
	{
		case 1:{
			typedef FEGeometry<Prism, dim, LagrangeLSFS<ReferencePrism, 1> ,
					GaussQuadrature<ReferencePrism, 2> > FEGeom;
			register_fe_func<Prism, FEGeom > (); break;
		}
		default:
			register_fe_func<Prism, DimFEGeometry<dim> > (); break;
	}

	//	Pyramid
	switch (order)
	{
		default: register_fe_func<Pyramid, DimFEGeometry<dim> > (); break;
	}

	//	Hexahedron
	switch (order)
	{
		case 1:{
			if (quadOrder == 2){
				typedef FEGeometry<Hexahedron, dim, LagrangeLSFS<
				ReferenceHexahedron, 1> , GaussQuadrature<
				ReferenceHexahedron, 2> > FEGeom;
				register_fe_func<Hexahedron, FEGeom > (); break;
			}
			else{
				typedef FEGeometry<Hexahedron, dim, LagrangeLSFS<
				ReferenceHexahedron, 1> , GaussQuadrature<
				ReferenceHexahedron, 3> > FEGeom;
				register_fe_func<Hexahedron, FEGeom > (); break;
			}
		}
		case 2:{
			typedef FEGeometry<Hexahedron, dim, LagrangeLSFS<
					ReferenceHexahedron, 2> , GaussQuadrature<
					ReferenceHexahedron, 7> > FEGeom;
			register_fe_func<Hexahedron, FEGeom > (); break;
		}
		case 3:{
			typedef FEGeometry<Hexahedron, dim, LagrangeLSFS<
					ReferenceHexahedron, 3> , GaussQuadrature<
					ReferenceHexahedron, 11> > FEGeom;
			register_fe_func<Hexahedron, FEGeom > (); break;
		}
		default: register_fe_func<Hexahedron, DimFEGeometry<dim> > (); break;
	}
}
#endif
template<typename TDomain>
template<typename TElem, typename TFEGeom>
void DeformationStateEquation<TDomain>::register_fe_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;
	static const int refDim = reference_element_traits<TElem>::dim;
	
	this->clear_add_fct(id);

	this->set_prep_elem_loop_fct(id, &T::template prep_elem_loop<TElem, TFEGeom>);
	this->set_prep_elem_fct(id, &T::template prep_elem<TElem, TFEGeom>);
	this->set_fsh_elem_loop_fct(id, &T::template fsh_elem_loop<TElem, TFEGeom>);

	this->set_prep_timestep_elem_fct(id, &T::template prep_timestep_elem<TElem, TFEGeom>);
	this->set_fsh_timestep_elem_fct(id, &T::template fsh_timestep_elem<TElem, TFEGeom>);

	this->set_add_jac_A_elem_fct(id, &T::template add_jac_A_elem<TElem, TFEGeom>);
	this->set_add_jac_M_elem_fct(id, &T::template add_jac_M_elem<TElem, TFEGeom>);
	this->set_add_def_A_elem_fct(id, &T::template add_def_A_elem<TElem, TFEGeom>);
	this->set_add_def_M_elem_fct(id, &T::template add_def_M_elem<TElem, TFEGeom>);
	this->set_add_rhs_elem_fct(id, &T::template add_rhs_elem<TElem, TFEGeom>);
	
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////
#ifdef UG_DIM_1
template class DeformationStateEquation<Domain1d> ;
#endif
#ifdef UG_DIM_2
template class DeformationStateEquation<Domain2d> ;
#endif
#ifdef UG_DIM_3
template class DeformationStateEquation<Domain3d> ;
#endif 


}
}

