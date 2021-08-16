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
#ifndef FLUID_OPTIM_TOOLS
#define FLUID_OPTIM_TOOLS

#include "lib_disc/function_spaces/integrate.h"
#include <vector>
#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_grid/lg_base.h"
#include "lib_grid/tools/periodic_boundary_manager.h"
#include "common/profiler/profiler.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/fvcr_geom.h"
#include "lib_disc/quadrature/quadrature_provider.h"
#include "lib_disc/function_spaces/grid_function_global_user_data.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

namespace ug{
namespace FluidOptim{

template <typename TGridFunction>
number SurfaceIntegral(SmartPtr<TGridFunction> spGridFct, SmartPtr<TGridFunction> spGridFct2,
                              const char* vCmp,
                              const char* BndSubsets, const char* InnerSubsets,
                              int quadOrder){
	static const int dim = TGridFunction::dim;
	static const int WorldDim = TGridFunction::dim;
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object grid_base_object;
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator const_iterator;


//	read subsets
	SubsetGroup innerSSGrp(spGridFct->domain()->subset_handler());
	if(InnerSubsets != NULL){
		innerSSGrp.add(TokenizeString(InnerSubsets));
		if(!SameDimensionsInAllSubsets(innerSSGrp))
			UG_THROW("SurfaceIntegral: Subsets '"<<InnerSubsets<<"' do not have same dimension."
					 "Can not integrate on subsets of different dimensions.");
	}
	else{
		innerSSGrp.add_all();
		RemoveLowerDimSubsets(innerSSGrp);
	}

//	read subsets
	SubsetGroup bndSSGrp(spGridFct->domain()->subset_handler());
	if(BndSubsets != NULL)
		bndSSGrp.add(TokenizeString(BndSubsets));
	else
		UG_THROW("SurfaceIntegral: No boundary subsets passed.");

//	get function group
	const FunctionGroup vFctID = spGridFct->fct_grp_by_name(vCmp);
	std::vector<LFEID> vLFEID;
	for(size_t fct = 0; fct < vFctID.size(); ++fct){
		vLFEID.push_back(spGridFct->lfeid(vFctID[fct]));
	}
	//the Grid Function Group 2
	const FunctionGroup vFctID2 = spGridFct2->fct_grp_by_name(vCmp);
	std::vector<LFEID> vLFEID2;
	for(size_t fct = 0; fct < vFctID2.size(); ++fct){
		vLFEID2.push_back(spGridFct2->lfeid(vFctID2[fct]));
	}
	

//	reset the result
	number dot_product = 0;

//	loop subsets
	for(size_t i = 0; i < innerSSGrp.size(); ++i)
	{
	//	get subset index
		const int si = innerSSGrp[i];

	//	skip empty subset
		if(innerSSGrp.dim(i) == DIM_SUBSET_EMPTY_GRID) continue;

	//	check dimension
		if(innerSSGrp.dim(i) != dim)
			UG_THROW("SurfaceIntegral: Dimension of inner subset is "<<
					 innerSSGrp.dim(i)<<", but only World Dimension "<<dim<<
					 " subsets can be used for inner subsets.");

	//	note: this iterator is for the base elements, e.g. Face and not
	//			for the special type, e.g. Triangle, Quadrilateral
		const_iterator iterBegin = spGridFct->template begin<grid_base_object>(si);
		const_iterator iterEnd = spGridFct->template end<grid_base_object>(si);
		const_iterator iter = iterBegin;

		typename domain_traits<TGridFunction::dim>::position_accessor_type& aaPos
			= spGridFct->domain()->position_accessor();
		const ISubsetHandler* ish = spGridFct->domain()->subset_handler().get();
		Grid& grid = *spGridFct->domain()->grid();

	//	this is the base element type (e.g. Face). This is the type when the
	//	iterators above are dereferenciated.
		typedef typename domain_traits<dim>::element_type Element;
		typedef typename domain_traits<dim>::side_type Side;

	//	vector of corner coordinates of element corners (to be filled for each elem)
		std::vector<MathVector<WorldDim> > vCorner;
		std::vector<int> vSubsetIndex;

	// 	iterate over all elements
		for(; iter != iterEnd; ++iter)
		{
		//	get element
			Element* pElem = *iter;

		//	get all corner coordinates
			CollectCornerCoordinates(vCorner, *pElem, aaPos, true);

		//	get reference object id
			const ReferenceObjectID elemRoid = pElem->reference_object_id();

		//	get sides
			typename Grid::traits<Side>::secure_container vSide;
			grid.associated_elements_sorted(vSide, pElem);
			vSubsetIndex.resize(vSide.size());
			for(size_t i = 0; i < vSide.size(); ++i)
				vSubsetIndex[i] = ish->get_subset_index(vSide[i]);

			DimReferenceMapping<dim, WorldDim>& rMapping
				= ReferenceMappingProvider::get<dim, WorldDim>(elemRoid, vCorner);

			const DimReferenceElement<dim>& rRefElem
				= ReferenceElementProvider::get<dim>(elemRoid);

		//	get element values
			std::vector<DoFIndex> vInd;
			std::vector<std::vector<number> > vvValue(vFctID.size());
			for(size_t fct = 0; fct < vvValue.size(); ++fct){
				spGridFct->dof_indices(pElem, vFctID[fct], vInd);
				vvValue[fct].resize(vInd.size());
				for(size_t sh = 0; sh < vInd.size(); ++sh)
					vvValue[fct][sh] = DoFRef(*spGridFct, vInd[sh]);
			}
			//for Grid Function 2
			std::vector<DoFIndex> vInd2;
			std::vector<std::vector<number> > vvValue2(vFctID2.size());
			for(size_t fct = 0; fct < vvValue2.size(); ++fct){
				spGridFct2->dof_indices(pElem, vFctID2[fct], vInd2);
				vvValue2[fct].resize(vInd2.size());
				for(size_t sh = 0; sh < vInd2.size(); ++sh)
					vvValue2[fct][sh] = DoFRef(*spGridFct2, vInd2[sh]);
			}
			
			const static int _C_ = 0;

		//	loop sub elements
			for(size_t side = 0; side < vSide.size(); ++side)
			{
			//	check if side used
				if(!bndSSGrp.contains(vSubsetIndex[side])) continue;

			//	get side
				Side* pSide = vSide[side];

				std::vector<MathVector<WorldDim> > vSideCorner(rRefElem.num(dim-1, side, 0));
				std::vector<MathVector<dim> > vLocalSideCorner(rRefElem.num(dim-1, side, 0));
				for(size_t co = 0; co < vSideCorner.size(); ++co){
					vSideCorner[co] = vCorner[rRefElem.id(dim-1, side, 0, co)];
					vLocalSideCorner[co] = rRefElem.corner(rRefElem.id(dim-1, side, 0, co));
				}

			//	side quad rule
				const ReferenceObjectID sideRoid = pSide->reference_object_id();
				const QuadratureRule<dim-1>& rSideQuadRule
						= QuadratureRuleProvider<dim-1>::get(sideRoid, quadOrder);

			//	quadrature points
				const number* vWeight = rSideQuadRule.weights();
				const size_t nip = rSideQuadRule.size();
				std::vector<MathVector<dim> > vLocalIP(nip);
				std::vector<MathVector<dim> > vGlobalIP(nip);

				DimReferenceMapping<dim-1, dim>& map
					= ReferenceMappingProvider::get<dim-1, dim>(sideRoid, vLocalSideCorner);

				for(size_t ip = 0; ip < nip; ++ip)
					map.local_to_global(vLocalIP[ip], rSideQuadRule.point(ip));
				
				//for(size_t ip = 0; ip < nip; ++ip)
					//vLocalIP[ip]=rSideQuadRule.point(ip); 
				
				for(size_t ip = 0; ip < nip; ++ip)
					rMapping.local_to_global(vGlobalIP[ip], vLocalIP[ip]);

			//	compute transformation matrices
				DimReferenceMapping<dim-1, dim>& map2
					= ReferenceMappingProvider::get<dim-1, dim>(sideRoid, vSideCorner);
				std::vector<MathMatrix<dim-1, WorldDim> > vJT(nip);
				map2.jacobian_transposed(&(vJT[0]), rSideQuadRule.points(), nip);

				std::vector<MathMatrix<dim, WorldDim> > vElemJT(nip);
				rMapping.jacobian_transposed(&(vElemJT[0]), &vLocalIP[0], nip);

			//	loop integration points
				for(size_t ip = 0; ip < nip; ++ip)
				{

				//	1. Interpolate grid function 1
					const LocalShapeFunctionSet<dim>& rTrialSpaceP =
							LocalFiniteElementProvider::get<dim>(elemRoid, vLFEID[_C_]);
					std::vector<number> vShapeGF1;
					rTrialSpaceP.shapes(vShapeGF1, vLocalIP[ip]);

					number value1 = 0.0;
					for(size_t sh = 0; sh < vvValue[_C_].size(); ++sh)
						value1 += vShapeGF1[sh] * vvValue[_C_][sh];
					
				//	2. Interpolate grid function 2
					const LocalShapeFunctionSet<dim>& rTrialSpaceP2 =
							LocalFiniteElementProvider::get<dim>(elemRoid, vLFEID2[_C_]);
					std::vector<number> vShapeGF2;
					rTrialSpaceP2.shapes(vShapeGF2, vLocalIP[ip]);

					number value2 = 0.0;
					for(size_t sh = 0; sh < vvValue2[_C_].size(); ++sh)
						value2 += vShapeGF2[sh] * vvValue2[_C_][sh];

				//	get quadrature weight
					const number weightIP = vWeight[ip];

				//	get determinate of mapping
					const number det = SqrtGramDeterminant(vJT[ip]);

				//	add contribution of integration point
					dot_product +=  value1*value2*weightIP*det;

				}
			} // end bf
		} // end elem
	} // end subsets

#ifdef UG_PARALLEL
	// sum over processes
	if(pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		number local = dot_product;
		com.allreduce(&local, &dot_product, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
	}
#endif
	return dot_product;
}

template <typename TGridFunction>
number VolumeIntegral(SmartPtr<TGridFunction> spGridFct, SmartPtr<TGridFunction> spGridFct2, 
					  const char* vCmp, const char* IntegrationSubsets, int quadOrder){
						
	number scalar_product=0.0;
	const static int _C_ = 0;
	static const int dim = TGridFunction::dim;
	static const int WorldDim = TGridFunction::dim;
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object grid_base_object;
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator const_iterator;
	
	 std::vector<DoFIndex> multInd;
	//	TODO: think about, we have 2 GridFncs, which one do we reference for the geometry access?
	// The subsets of definition for the GridFncs must match (ideally as far as i understand), also the 
	// order of the interpolations, else we have to do duplicate some code here. 
	//	Create SubsetGroup from the given subsets for integration, if empty use all
	SubsetGroup intSSGrp(spGridFct->domain()->subset_handler());
	intSSGrp.add(TokenizeString(IntegrationSubsets));
	
	
	//	Create function group from GridFunction components
	const FunctionGroup vFctID = spGridFct->fct_grp_by_name(vCmp);//should be of size 1 (only 1 function)
	
	std::vector<LFEID> vLFEID;
	//Extract and store the LFE id from the FunctionPattern through the FunctionGroup
	for(size_t fct = 0; fct < vFctID.size(); ++fct){
		vLFEID.push_back(spGridFct->lfeid(vFctID[fct]));
	}
	
	//	Create function group from GridFunction components
	const FunctionGroup vFctID2 = spGridFct2->fct_grp_by_name(vCmp);//should be of size 1 (only 1 function)
	
	std::vector<LFEID> vLFEID2;
	//Extract and store the LFE id from the FunctionPattern through the FunctionGroup
	for(size_t fct = 0; fct < vFctID2.size(); ++fct){
		vLFEID2.push_back(spGridFct2->lfeid(vFctID2[fct]));
	}
	
	//Loop subsets assigned to the integration subset group
	for(size_t i = 0; i < intSSGrp.size(); ++i)
	{
		//Store subset index
		const int si = intSSGrp[i];
		//TODO:the element iterators are a reason why the integration sets must match...else we would only have the 
		//set intersections in common, I guess
		const_iterator iterBegin = spGridFct->template begin<grid_base_object>(si);
		const_iterator iterEnd = spGridFct->template end<grid_base_object>(si);
		const_iterator iter = iterBegin;
		typename domain_traits<TGridFunction::dim>::position_accessor_type& aaPos= spGridFct->domain()->
																				   position_accessor();
		//const ISubsetHandler* ish = spGridFct->domain()->subset_handler().get();
		//Grid& grid = *spGridFct->domain()->grid();
		
		typedef typename domain_traits<dim>::element_type Element;

		//vector of MathVectors sized dim, will store the corners of each face
		std::vector<MathVector<WorldDim> > vCorner;
		std::vector<int> vSubsetIndex;
		//For storing a Jacobian per ip for a given element
		//std::vector<MathMatrix<dim, WorldDim> > vJT;
		std::vector<MathMatrix<dim, WorldDim> > vElemJT;

		for(; iter != iterEnd; ++iter)
		{
			//Pointer to element obtained from dereferencing the iterator...see note above
			Element* pElem = *iter;
		
			CollectCornerCoordinates(vCorner, *pElem, aaPos, true);
			
			//This tells us which kind of geometric object it is, e.g. ROID_TRIANGLE 	
			const ReferenceObjectID elemRoid = pElem->reference_object_id();
			
			try{
				//Do we need the type or use best?
				const QuadratureRule<dim>& rQuadRule = QuadratureRuleProvider<dim>::get(elemRoid, quadOrder);

			//	get reference element mapping by reference object id
				DimReferenceMapping<dim, WorldDim>& rMapping = ReferenceMappingProvider::
																get<dim, WorldDim>(elemRoid,vCorner);

			//	number of integration points
				const number* vWeight = rQuadRule.weights();
				const size_t numIP = rQuadRule.size();
				std::vector<MathVector<dim> > vLocalIP(numIP);
				std::vector<MathVector<dim> > vGlobalIP(numIP);
				
				//	Get  Values in the Element
				//Stores integers with dof identifiers
				std::vector<DoFIndex> vInd;
				//vector of vectors sized #functions in function group
				std::vector<std::vector<number> > vvValue(vFctID.size());
				//We visit each vector inside the vvValue (1 for each symbolic function)
				for(size_t fct = 0; fct < vvValue.size(); ++fct){
					spGridFct->dof_indices(pElem, vFctID[fct], vInd);
					vvValue[fct].resize(vInd.size());
					for(size_t sh = 0; sh < vInd.size(); ++sh){
						//DoFRef will give the value of GridFunction at index provided (dof)
						vvValue[fct][sh] = DoFRef(*spGridFct, vInd[sh]);
					}//end sh value for
				}//end fct value for
				
				//fct2
				//	Get  Values in the Element
				//Stores integers with dof identifiers
				std::vector<DoFIndex> vInd2;
				//vector of vectors sized #functions in function group
				std::vector<std::vector<number> > vvValue2(vFctID2.size());
				//We visit each vector inside the vvValue2 (1 for each symbolic function)
				for(size_t fct = 0; fct < vvValue2.size(); ++fct){
					spGridFct2->dof_indices(pElem, vFctID2[fct], vInd2);
					vvValue2[fct].resize(vInd2.size());
					for(size_t sh = 0; sh < vInd2.size(); ++sh){
						//DoFRef will give the value of GridFunction at index provided (dof)
						vvValue2[fct][sh] = DoFRef(*spGridFct2, vInd2[sh]);
					}//end sh value for
				}//end fct value for
								
				for(size_t ip = 0; ip < numIP; ++ip)
					vLocalIP[ip]=rQuadRule.point(ip); 
				
				for(size_t ip = 0; ip < numIP; ++ip)
					rMapping.local_to_global(vGlobalIP[ip], vLocalIP[ip]);
				
				
				//Transformation matrices
				//std::vector<MathMatrix<dim, WorldDim> > vElemJT(numIP);
				vElemJT.resize(numIP);
				rMapping.jacobian_transposed(&(vElemJT[0]), &vLocalIP[0], numIP);
				
				//element integral value, reset to zero for each element
				number elemValue = 0;
				//loop integration points
				for(size_t ip = 0; ip < numIP; ++ip)
				{
				//	1. Interpolate pressure at ip
					const LocalShapeFunctionSet<dim>& rTrialSpaceP =
							LocalFiniteElementProvider::get<dim>(elemRoid, vLFEID[_C_]);
					std::vector<number> vShapeP;
					rTrialSpaceP.shapes(vShapeP, vLocalIP[ip]);

					number value1 = 0.0;
					for(size_t sh = 0; sh < vvValue[_C_].size(); ++sh)
						value1 += vShapeP[sh] * vvValue[_C_][sh];
					
				//	1. Interpolate value2
					const LocalShapeFunctionSet<dim>& rTrialSpaceP2 =
							LocalFiniteElementProvider::get<dim>(elemRoid, vLFEID2[_C_]);
					std::vector<number> vShapeP2;
					rTrialSpaceP2.shapes(vShapeP2, vLocalIP[ip]);

					number value2 = 0.0;
					for(size_t sh = 0; sh < vvValue2[_C_].size(); ++sh)
						value2 += vShapeP2[sh] * vvValue2[_C_][sh];
					
					//5. Do the integration
					//	get quadrature weight
					const number weightIP = vWeight[ip];
					//const number weightIP = rQuadRule.weight(ip);
					//	get determinate of mapping
					const number det = SqrtGramDeterminant(vElemJT[ip]);
					//0.5*chi*(ext-0.5*(ext_upper+ext_lower))*(ext-0.5*(ext_upper+ext_lower))
					//0.5*chi*(extension_factor-0.5*(ext_upper+ext_lower))*(extension_factor-0.5*(ext_upper+ext_lower));
					elemValue +=  value1 * value2 * weightIP * det ;					
				}//end ip for
				
				scalar_product += elemValue;	
				
			}UG_CATCH_THROW("SumValuesOnElems failed.");//end try
		}//end element iteratior for
	}//end intSSGrp for
	
#ifdef UG_PARALLEL
	// sum over processes
	if(pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		number local = scalar_product;
		com.allreduce(&local, &scalar_product, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
	}
#endif
	return scalar_product;

}//end volumeintegral

template <typename TGridFunction>
number Extension(SmartPtr<TGridFunction> spGridFct, number ext_upper, number ext_lower, 
					const char* IntegrationSubsets, const char* vCmp, int quadOrder){
						
	number extension_value=0.0;
	const static int _C_ = 0;
	static const int dim = TGridFunction::dim;
	static const int WorldDim = TGridFunction::dim;
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object grid_base_object;
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator const_iterator;
	
	 std::vector<DoFIndex> multInd;
	//	TODO: think about, we have 2 GridFncs, which one do we reference for the geometry access?
	// The subsets of definition for the GridFncs must match (ideally as far as i understand), also the 
	// order of the interpolations, else we have to do duplicate some code here. 
	//	Create SubsetGroup from the given subsets for integration, if empty use all
	SubsetGroup intSSGrp(spGridFct->domain()->subset_handler());
	intSSGrp.add(TokenizeString(IntegrationSubsets));
	
	
	//	Create function group from GridFunction components
	const FunctionGroup vFctID = spGridFct->fct_grp_by_name(vCmp);//should be of size 1 (only 1 function)
	
	std::vector<LFEID> vLFEID;
	//Extract and store the LFE id from the FunctionPattern through the FunctionGroup
	for(size_t fct = 0; fct < vFctID.size(); ++fct){
		vLFEID.push_back(spGridFct->lfeid(vFctID[fct]));
	}
	
	//Loop subsets assigned to the integration subset group
	for(size_t i = 0; i < intSSGrp.size(); ++i)
	{
		//Store subset index
		const int si = intSSGrp[i];
		//TODO:the element iterators are a reason why the integration sets must match...else we would only have the 
		//set intersections in common, I guess
		const_iterator iterBegin = spGridFct->template begin<grid_base_object>(si);
		const_iterator iterEnd = spGridFct->template end<grid_base_object>(si);
		const_iterator iter = iterBegin;
		typename domain_traits<TGridFunction::dim>::position_accessor_type& aaPos= spGridFct->domain()->
																				   position_accessor();
		//const ISubsetHandler* ish = spGridFct->domain()->subset_handler().get();
		//Grid& grid = *spGridFct->domain()->grid();
		
		typedef typename domain_traits<dim>::element_type Element;

		//vector of MathVectors sized dim, will store the corners of each face
		std::vector<MathVector<WorldDim> > vCorner;
		std::vector<int> vSubsetIndex;
		//For storing a Jacobian per ip for a given element
		//std::vector<MathMatrix<dim, WorldDim> > vJT;
		std::vector<MathMatrix<dim, WorldDim> > vElemJT;

		for(; iter != iterEnd; ++iter)
		{
			//Pointer to element obtained from dereferencing the iterator...see note above
			Element* pElem = *iter;
		
			CollectCornerCoordinates(vCorner, *pElem, aaPos, true);
			
			//This tells us which kind of geometric object it is, e.g. ROID_TRIANGLE 	
			const ReferenceObjectID elemRoid = pElem->reference_object_id();
			
			try{
				//Do we need the type or use best?
				const QuadratureRule<dim>& rQuadRule = QuadratureRuleProvider<dim>::get(elemRoid, quadOrder);

			//	get reference element mapping by reference object id
				DimReferenceMapping<dim, WorldDim>& rMapping = ReferenceMappingProvider::
																get<dim, WorldDim>(elemRoid,vCorner);

			//	number of integration points
				const number* vWeight = rQuadRule.weights();
				const size_t numIP = rQuadRule.size();
				std::vector<MathVector<dim> > vLocalIP(numIP);
				std::vector<MathVector<dim> > vGlobalIP(numIP);
				
				//	Get Velocity Values in the Element
				//Stores integers with dof identifiers
				std::vector<DoFIndex> vInd;
				//vector of vectors sized #functions in function group
				std::vector<std::vector<number> > vvValue(vFctID.size());
				//We visit each vector inside the vvValue (1 for each symbolic function)
				for(size_t fct = 0; fct < vvValue.size(); ++fct){
					spGridFct->dof_indices(pElem, vFctID[fct], vInd);
					vvValue[fct].resize(vInd.size());
					for(size_t sh = 0; sh < vInd.size(); ++sh){
						//DoFRef will give the value of GridFunction at index provided (dof)
						vvValue[fct][sh] = DoFRef(*spGridFct, vInd[sh]);
					}//end sh value for
				}//end fct value for
								
				for(size_t ip = 0; ip < numIP; ++ip)
					vLocalIP[ip]=rQuadRule.point(ip); 
				
				for(size_t ip = 0; ip < numIP; ++ip)
					rMapping.local_to_global(vGlobalIP[ip], vLocalIP[ip]);
				
				
				//Transformation matrices
				//std::vector<MathMatrix<dim, WorldDim> > vElemJT(numIP);
				vElemJT.resize(numIP);
				rMapping.jacobian_transposed(&(vElemJT[0]), &vLocalIP[0], numIP);
				
				//element integral value, reset to zero for each element
				number elemExtension = 0;
				//loop integration points
				for(size_t ip = 0; ip < numIP; ++ip)
				{
				//	1. Interpolate pressure at ip
					const LocalShapeFunctionSet<dim>& rTrialSpaceP =
							LocalFiniteElementProvider::get<dim>(elemRoid, vLFEID[_C_]);
					std::vector<number> vShapeP;
					rTrialSpaceP.shapes(vShapeP, vLocalIP[ip]);

					number extension_factor = 0.0;
					for(size_t sh = 0; sh < vvValue[_C_].size(); ++sh)
						extension_factor += vShapeP[sh] * vvValue[_C_][sh];
					
					//5. Do the integration
					//	get quadrature weight
					const number weightIP = vWeight[ip];
					//const number weightIP = rQuadRule.weight(ip);
					//	get determinate of mapping
					const number det = SqrtGramDeterminant(vElemJT[ip]);
					//0.5*chi*(ext-0.5*(ext_upper+ext_lower))*(ext-0.5*(ext_upper+ext_lower))
					//0.5*chi*(extension_factor-0.5*(ext_upper+ext_lower))*(extension_factor-0.5*(ext_upper+ext_lower));
					elemExtension += (extension_factor-0.5*(ext_upper+ext_lower))*(extension_factor-0.5*(ext_upper+ext_lower)) 
										* weightIP * det ;					
				}//end ip for
				
				extension_value += elemExtension;	
				
			}UG_CATCH_THROW("SumValuesOnElems failed.");//end try
		}//end element iteratior for
	}//end intSSGrp for
	
#ifdef UG_PARALLEL
	// sum over processes
	if(pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		number local = extension_value;
		com.allreduce(&local, &extension_value, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
	}
#endif
	return extension_value;

}//end threshold

template <typename TGridFunction>
number Threshold(SmartPtr<TGridFunction> spGridFct, number threshold_value, 
					const char* IntegrationSubsets, const char* vCmp, int quadOrder){
						
	number threshold=0.0;
	
	static const int dim = TGridFunction::dim;
	static const int WorldDim = TGridFunction::dim;
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object grid_base_object;
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator const_iterator;
	
	 std::vector<DoFIndex> multInd;

	
	//	Create SubsetGroup from the given subsets for integration, if empty use all
	SubsetGroup intSSGrp(spGridFct->domain()->subset_handler());
	intSSGrp.add(TokenizeString(IntegrationSubsets));
	
	
	//	Create function group from GridFunction components
	const FunctionGroup vFctID = spGridFct->fct_grp_by_name(vCmp);
	// If every variable has its own different Local Finite Element, then each can be accessed here individually
	// usually (for instance for a deformation field), the LFEs will be the same. 
	std::vector<LFEID> vLFEID;
	//Extract and store the LFE id from the FunctionPattern through the FunctionGroup
	for(size_t fct = 0; fct < vFctID.size(); ++fct){
		vLFEID.push_back(spGridFct->lfeid(vFctID[fct]));
	}
	
	//Loop subsets assigned to the integration subset group
	for(size_t i = 0; i < intSSGrp.size(); ++i)
	{
		//Store subset index
		const int si = intSSGrp[i];
		
	//	note: this iterator is for the base elements, e.g. Face and not
	//	for the special type, e.g. Triangle, Quadrilateral
	// With these we will access every element from the subset
		const_iterator iterBegin = spGridFct->template begin<grid_base_object>(si);
		const_iterator iterEnd = spGridFct->template end<grid_base_object>(si);
		const_iterator iter = iterBegin;
		typename domain_traits<TGridFunction::dim>::position_accessor_type& aaPos= spGridFct->domain()->position_accessor();
;
		
	//	this is the base element type (e.g. Face). This is the type when the
	//	iterators above are dereferenciated.
	//  What he said...
		typedef typename domain_traits<dim>::element_type Element;
		//typedef typename domain_traits<dim>::side_type Side;

	//	vector of corner coordinates of element corners (to be filled for each elem)
		//vector of MathVectors sized dim, will store the corners of each face
		std::vector<MathVector<WorldDim> > vCorner;
		std::vector<int> vSubsetIndex;
		//For storing a Jacobian per ip for a given element
		//std::vector<MathMatrix<dim, WorldDim> > vJT;
		std::vector<MathMatrix<dim, WorldDim> > vElemJT;
		
		
		//We start to access every element using our iterators (give a read to the topic online or on c++ literature)
		//...pointer shit...STL related things...
		for(; iter != iterEnd; ++iter)
		{
			//Pointer to element obtained from dereferencing the iterator...see note above
			Element* pElem = *iter;
		
			//	get all corner coordinatesg
			CollectCornerCoordinates(vCorner, *pElem, aaPos, true);
			
			//This tells us which kind of geometric object it is, e.g. ROID_TRIANGLE 	
			const ReferenceObjectID elemRoid = pElem->reference_object_id();
			
			try{
				//Do we need the type or use best?
				const QuadratureRule<dim>& rQuadRule = QuadratureRuleProvider<dim>::get(elemRoid, quadOrder);

			//	get reference element mapping by reference object id
				DimReferenceMapping<dim, WorldDim>& rMapping = ReferenceMappingProvider::get<dim, WorldDim>(elemRoid,vCorner);

			//	number of integration points
				const number* vWeight = rQuadRule.weights();
				const size_t numIP = rQuadRule.size();
				std::vector<MathVector<dim> > vLocalIP(numIP);
				std::vector<MathVector<dim> > vGlobalIP(numIP);
				
				//	get element values
				//Stores integers with dof identifiers
				std::vector<DoFIndex> vInd;
				//vector of vectors sized #functions in function group
				std::vector<std::vector<number> > vvValue(vFctID.size());
				//We visit each vector inside the vvValue (1 for each symbolic function)
				for(size_t fct = 0; fct < vvValue.size(); ++fct){
					
					spGridFct->dof_indices(pElem, vFctID[fct], vInd);
					
					//The dof vector (integers) now contains the indices for the corresponding function in vFctID
					//Then we resize the vector inside vvValue according to how many dofs we have in this function 
					vvValue[fct].resize(vInd.size());
					for(size_t sh = 0; sh < vInd.size(); ++sh){
						//DoFRef will give the value of GridFunction at index provided (dof)
						vvValue[fct][sh] = DoFRef(*spGridFct, vInd[sh]);
					}//end sh value for
				}//end fct value for
				
				
				for(size_t ip = 0; ip < numIP; ++ip)
					vLocalIP[ip]=rQuadRule.point(ip); 
				
				for(size_t ip = 0; ip < numIP; ++ip)
					rMapping.local_to_global(vGlobalIP[ip], vLocalIP[ip]);
				
				//Transformation matrices
				vElemJT.resize(numIP);
				rMapping.jacobian_transposed(&(vElemJT[0]), &vLocalIP[0], numIP);
				
				//element integral value
				number intValElem = 0;
				//loop integration points
				for(size_t ip = 0; ip < numIP; ++ip)
				{

				// 	1. Interpolate Functional Matrix of Identity Gradient at ip
					std::vector<MathVector<dim> > vvLocGradV[dim];
					std::vector<MathVector<dim> > vvGradV[dim];
					MathMatrix<dim, dim> JTInv;
					Inverse(JTInv, vElemJT[ip]);
					for(int d1 = 0; d1 < dim; ++d1){						
						const LocalShapeFunctionSet<dim>& rTrialSpaceP = LocalFiniteElementProvider::get<dim>(elemRoid, vLFEID[d1]);
						rTrialSpaceP.grads(vvLocGradV[d1], vLocalIP[ip]);
							vvGradV[d1].resize(vvLocGradV[d1].size());
							for(size_t sh = 0; sh < vvGradV[d1].size(); ++sh)
								MatVecMult(vvGradV[d1][sh], JTInv, vvLocGradV[d1][sh]);
					}

					MathMatrix<dim, dim> IdGrad;
					for(int d1 = 0; d1 < dim; ++d1){
						for(int d2 = 0; d2 < dim; ++d2){
							IdGrad(d1, d2) = 0.0;
							for(size_t sh = 0; sh < vvValue[d1].size(); ++sh){									
								IdGrad(d1, d2) += vvValue[d1][sh] * vvGradV[d1][sh][d2];
							}
						}
					}
					for(int d1=0;d1<dim; ++d1 ){
						IdGrad(d1,d1) += 1;
					}
					const number detDF= Determinant(IdGrad);
				
					//	get quadrature weight
					const number weightIP = vWeight[ip];

					//	get determinate of mapping
					const number det = SqrtGramDeterminant(vElemJT[ip]);
					number pos_part= ((threshold_value-detDF) < 0)? 0:(threshold_value-detDF);//(m_detThreshold-dDF);
					
					intValElem +=  pos_part*pos_part * weightIP * det ;				
				}//end ip for
				
				threshold += intValElem;	
				
			}UG_CATCH_THROW("SumValuesOnElems failed.");
			//end try
			
		}//end element iteratior for
		
		
	}//end intSSGrp for
	
#ifdef UG_PARALLEL
	// sum over processes
	if(pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		number local = threshold;
		com.allreduce(&local, &threshold, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
	}
#endif
	return threshold;

}//end threshold
	
template <typename TGridFunction>
number Control(SmartPtr<TGridFunction> spGridFct,
                              const char* vCmp,
                              const char* IntegrationSubsets,                           
                              int quadOrder){
	number control=0.0;
	const static int _C_ = 0;
	static const int dim = TGridFunction::dim;
	static const int WorldDim = TGridFunction::dim;
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object grid_base_object;
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator const_iterator;
	
	 std::vector<DoFIndex> multInd;
	//	TODO: think about, we have 2 GridFncs, which one do we reference for the geometry access?
	// The subsets of definition for the GridFncs must match (ideally as far as i understand), also the 
	// order of the interpolations, else we have to do duplicate some code here. 
	//	Create SubsetGroup from the given subsets for integration, if empty use all
	SubsetGroup intSSGrp(spGridFct->domain()->subset_handler());
	intSSGrp.add(TokenizeString(IntegrationSubsets));
	
	
	//	Create function group from GridFunction components
	const FunctionGroup vFctID = spGridFct->fct_grp_by_name(vCmp);//should be of size 1 (only 1 function)
	
	std::vector<LFEID> vLFEID;
	//Extract and store the LFE id from the FunctionPattern through the FunctionGroup
	for(size_t fct = 0; fct < vFctID.size(); ++fct){
		vLFEID.push_back(spGridFct->lfeid(vFctID[fct]));
	}
	
	//Loop subsets assigned to the integration subset group
	for(size_t i = 0; i < intSSGrp.size(); ++i)
	{
		//Store subset index
		const int si = intSSGrp[i];
		//TODO:the element iterators are a reason why the integration sets must match...else we would only have the 
		//set intersections in common, I guess
		const_iterator iterBegin = spGridFct->template begin<grid_base_object>(si);
		const_iterator iterEnd = spGridFct->template end<grid_base_object>(si);
		const_iterator iter = iterBegin;
		typename domain_traits<TGridFunction::dim>::position_accessor_type& aaPos= spGridFct->domain()->
																				   position_accessor();
		//const ISubsetHandler* ish = spGridFct->domain()->subset_handler().get();
		//Grid& grid = *spGridFct->domain()->grid();
		
		typedef typename domain_traits<dim>::element_type Element;

		//vector of MathVectors sized dim, will store the corners of each face
		std::vector<MathVector<WorldDim> > vCorner;
		std::vector<int> vSubsetIndex;
		//For storing a Jacobian per ip for a given element
		//std::vector<MathMatrix<dim, WorldDim> > vJT;
		std::vector<MathMatrix<dim, WorldDim> > vElemJT;

		for(; iter != iterEnd; ++iter)
		{
			//Pointer to element obtained from dereferencing the iterator...see note above
			Element* pElem = *iter;
		
			CollectCornerCoordinates(vCorner, *pElem, aaPos, true);
			
			//This tells us which kind of geometric object it is, e.g. ROID_TRIANGLE 	
			const ReferenceObjectID elemRoid = pElem->reference_object_id();
			
			try{
				//Do we need the type or use best?
				const QuadratureRule<dim>& rQuadRule = QuadratureRuleProvider<dim>::get(elemRoid, quadOrder);

			//	get reference element mapping by reference object id
				DimReferenceMapping<dim, WorldDim>& rMapping = ReferenceMappingProvider::
																get<dim, WorldDim>(elemRoid,vCorner);

			//	number of integration points
				const number* vWeight = rQuadRule.weights();
				const size_t numIP = rQuadRule.size();
				std::vector<MathVector<dim> > vLocalIP(numIP);
				std::vector<MathVector<dim> > vGlobalIP(numIP);
				
				//	Get Velocity Values in the Element
				//Stores integers with dof identifiers
				std::vector<DoFIndex> vInd;
				//vector of vectors sized #functions in function group
				std::vector<std::vector<number> > vvValue(vFctID.size());
				//We visit each vector inside the vvValue (1 for each symbolic function)
				for(size_t fct = 0; fct < vvValue.size(); ++fct){
					spGridFct->dof_indices(pElem, vFctID[fct], vInd);
					vvValue[fct].resize(vInd.size());
					for(size_t sh = 0; sh < vInd.size(); ++sh){
						//DoFRef will give the value of GridFunction at index provided (dof)
						vvValue[fct][sh] = DoFRef(*spGridFct, vInd[sh]);
					}//end sh value for
				}//end fct value for
				
				
				
				
				for(size_t ip = 0; ip < numIP; ++ip)
					vLocalIP[ip]=rQuadRule.point(ip); 
				
				for(size_t ip = 0; ip < numIP; ++ip)
					rMapping.local_to_global(vGlobalIP[ip], vLocalIP[ip]);
				
				
				//Transformation matrices
				//std::vector<MathMatrix<dim, WorldDim> > vElemJT(numIP);
				vElemJT.resize(numIP);
				rMapping.jacobian_transposed(&(vElemJT[0]), &vLocalIP[0], numIP);
				
				//element integral value, reset to zero for each element
				number elemControl = 0;
				//loop integration points
				for(size_t ip = 0; ip < numIP; ++ip)
				{
				//	1. Interpolate pressure at ip
					const LocalShapeFunctionSet<dim>& rTrialSpaceP =
							LocalFiniteElementProvider::get<dim>(elemRoid, vLFEID[_C_]);
					std::vector<number> vShapeP;
					rTrialSpaceP.shapes(vShapeP, vLocalIP[ip]);

					number boundary_control = 0.0;
					for(size_t sh = 0; sh < vvValue[_C_].size(); ++sh)
						boundary_control += vShapeP[sh] * vvValue[_C_][sh];
					
					//5. Do the integration
					//	get quadrature weight
					const number weightIP = vWeight[ip];
					//const number weightIP = rQuadRule.weight(ip);
					//	get determinate of mapping
					const number det = SqrtGramDeterminant(vElemJT[ip]);
					
					elemControl += boundary_control*boundary_control * weightIP * det ;					
				}//end ip for
				
				control += elemControl;	
				
			}UG_CATCH_THROW("SumValuesOnElems failed.");//end try
		}//end element iteratior for
	}//end intSSGrp for
	
#ifdef UG_PARALLEL
	// sum over processes
	if(pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		number local = control;
		com.allreduce(&local, &control, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
	}
#endif
	return control;
}//end control	
	
template <typename TGridFunction>
number Drag(SmartPtr<TGridFunction> spGridFctDef,
			SmartPtr<TGridFunction> spGridFctVel,
                              const char* vCmpVel,
							  const char* vCmpDef,
                              const char* IntegrationSubsets,                           
                              int quadOrder){
	number drag=0.0;
	static const int dim = TGridFunction::dim;
	static const int WorldDim = TGridFunction::dim;
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object grid_base_object;
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator const_iterator;
	
	 std::vector<DoFIndex> multInd;
	//	TODO: think about, we have 2 GridFncs, which one do we reference for the geometry access?
	// The subsets of definition for the GridFncs must match (ideally as far as i understand), also the 
	// order of the interpolations, else we have to do duplicate some code here. 
	//	Create SubsetGroup from the given subsets for integration, if empty use all
	SubsetGroup intSSGrp(spGridFctVel->domain()->subset_handler());
	intSSGrp.add(TokenizeString(IntegrationSubsets));
	
	
	//	Create function group from GridFunction components
	const FunctionGroup vFctID = spGridFctVel->fct_grp_by_name(vCmpVel);
	const FunctionGroup vFctID_Def = spGridFctDef->fct_grp_by_name(vCmpDef);
	
	std::vector<LFEID> vLFEID;
	//Extract and store the LFE id from the FunctionPattern through the FunctionGroup
	for(size_t fct = 0; fct < vFctID.size(); ++fct){
		vLFEID.push_back(spGridFctVel->lfeid(vFctID[fct]));
	}
	
	std::vector<LFEID> vLFEID_DEF;
	//Extract and store the LFE id from the FunctionPattern through the FunctionGroup
	for(size_t fct = 0; fct < vFctID_Def.size(); ++fct){
		vLFEID_DEF.push_back(spGridFctDef->lfeid(vFctID_Def[fct]));
	}
	
	//Loop subsets assigned to the integration subset group
	for(size_t i = 0; i < intSSGrp.size(); ++i)
	{
		//Store subset index
		const int si = intSSGrp[i];
		//TODO:the element iterators are a reason why the integration sets must match...else we would only have the 
		//set intersections in common, I guess
		const_iterator iterBegin = spGridFctVel->template begin<grid_base_object>(si);
		const_iterator iterEnd = spGridFctVel->template end<grid_base_object>(si);
		const_iterator iter = iterBegin;
		typename domain_traits<TGridFunction::dim>::position_accessor_type& aaPos= spGridFctVel->domain()->
																				   position_accessor();
		//const ISubsetHandler* ish = spGridFct->domain()->subset_handler().get();
		//Grid& grid = *spGridFct->domain()->grid();
		
		typedef typename domain_traits<dim>::element_type Element;

		//vector of MathVectors sized dim, will store the corners of each face
		std::vector<MathVector<WorldDim> > vCorner;
		std::vector<int> vSubsetIndex;
		//For storing a Jacobian per ip for a given element
		//std::vector<MathMatrix<dim, WorldDim> > vJT;
		std::vector<MathMatrix<dim, WorldDim> > vElemJT;

		for(; iter != iterEnd; ++iter)
		{
			//Pointer to element obtained from dereferencing the iterator...see note above
			Element* pElem = *iter;
		
			CollectCornerCoordinates(vCorner, *pElem, aaPos, true);
			
			//This tells us which kind of geometric object it is, e.g. ROID_TRIANGLE 	
			const ReferenceObjectID elemRoid = pElem->reference_object_id();
			
			try{
				//Do we need the type or use best?
				const QuadratureRule<dim>& rQuadRule = QuadratureRuleProvider<dim>::get(elemRoid, quadOrder);

			//	get reference element mapping by reference object id
				DimReferenceMapping<dim, WorldDim>& rMapping = ReferenceMappingProvider::get<dim, WorldDim>(elemRoid,vCorner);

			//	number of integration points
				const number* vWeight = rQuadRule.weights();
				const size_t numIP = rQuadRule.size();
				std::vector<MathVector<dim> > vLocalIP(numIP);
				std::vector<MathVector<dim> > vGlobalIP(numIP);
				
				//	Get Velocity Values in the Element
				//Stores integers with dof identifiers
				std::vector<DoFIndex> vInd;
				//vector of vectors sized #functions in function group
				std::vector<std::vector<number> > vvValueVel(vFctID.size());
				//We visit each vector inside the vvValue (1 for each symbolic function)
				for(size_t fct = 0; fct < vvValueVel.size(); ++fct){
					spGridFctVel->dof_indices(pElem, vFctID[fct], vInd);
					vvValueVel[fct].resize(vInd.size());
					for(size_t sh = 0; sh < vInd.size(); ++sh){
						//DoFRef will give the value of GridFunction at index provided (dof)
						vvValueVel[fct][sh] = DoFRef(*spGridFctVel, vInd[sh]);
					}//end sh value for
				}//end fct value for
				
				//	Get Deformation Values in the Element
				//Stores integers with dof identifiers
				std::vector<DoFIndex> vIndDef;
				//vector of vectors sized #functions in function group
				std::vector<std::vector<number> > vvValueDef(vFctID_Def.size());
				//We visit each vector inside the vvValue (1 for each symbolic function)
				for(size_t fct = 0; fct < vvValueDef.size(); ++fct){
					spGridFctDef->dof_indices(pElem, vFctID_Def[fct], vIndDef);
					vvValueDef[fct].resize(vIndDef.size());
					for(size_t sh = 0; sh < vIndDef.size(); ++sh){
						//DoFRef will give the value of GridFunction at index provided (dof)
						vvValueDef[fct][sh] = DoFRef(*spGridFctDef, vIndDef[sh]);
					}//end sh value for
				}//end fct value for
				
				
				
				for(size_t ip = 0; ip < numIP; ++ip)
					vLocalIP[ip]=rQuadRule.point(ip); 
				
				for(size_t ip = 0; ip < numIP; ++ip)
					rMapping.local_to_global(vGlobalIP[ip], vLocalIP[ip]);
				
				//Transformation matrices
				//std::vector<MathMatrix<dim, WorldDim> > vElemJT(numIP);
				vElemJT.resize(numIP);
				rMapping.jacobian_transposed(&(vElemJT[0]), &vLocalIP[0], numIP);
				
				//element integral value, reset to zero for each element
				number elemDrag = 0;
				//loop integration points
				for(size_t ip = 0; ip < numIP; ++ip)
				{
				// 	1. Interpolate Functional Matrix of Identity Gradient at ip
				
				//Stores the gradients (a vector size dim ) in a matrix with dim columns
					std::vector<MathVector<dim> > vvLocGradV[dim];
					std::vector<MathVector<dim> > vvGradV[dim];
					MathMatrix<dim, dim> JTInv;
					Inverse(JTInv, vElemJT[ip]);
					for(int d1 = 0; d1 < dim; ++d1){
						
						const LocalShapeFunctionSet<dim>& rTrialSpaceP = 
						LocalFiniteElementProvider::get<dim>(elemRoid, vLFEID[d1]);
						// grads returns returns all gradients evaluated at a several points
						// This function returns the gradients of all Shape Functions at several element-local evaluation point in an array.
						//Parameters
						//	[out]	vvGrad	Vector of gradients
						//	[in]	vLocPos	Vector of Position on reference element

						rTrialSpaceP.grads(vvLocGradV[d1], vLocalIP[ip]);

							vvGradV[d1].resize(vvLocGradV[d1].size());
							for(size_t sh = 0; sh < vvGradV[d1].size(); ++sh)
								MatVecMult(vvGradV[d1][sh], JTInv, vvLocGradV[d1][sh]);
					}
					
				// 	1.1 Interpolate Functional Matrix of Velocity Gradient
				
				//Stores the gradients (a vector size dim ) in a matrix with dim columns
					std::vector<MathVector<dim> > vvLocGradV_Def[dim];
					std::vector<MathVector<dim> > vvGradV_Def[dim];
					MathMatrix<dim, dim> JTInv_Def;
					Inverse(JTInv_Def, vElemJT[ip]);
					for(int d1 = 0; d1 < dim; ++d1){
						
						const LocalShapeFunctionSet<dim>& rTrialSpaceP_Def = 
						LocalFiniteElementProvider::get<dim>(elemRoid, vLFEID_DEF[d1]);
						// grads returns returns all gradients evaluated at a several points
						// This function returns the gradients of all Shape Functions at several element-local evaluation point in an array.
						//Parameters
						//	[out]	vvGrad	Vector of gradients
						//	[in]	vLocPos	Vector of Position on reference element

						rTrialSpaceP_Def.grads(vvLocGradV_Def[d1], vLocalIP[ip]);

							vvGradV_Def[d1].resize(vvLocGradV_Def[d1].size());
							for(size_t sh = 0; sh < vvGradV_Def[d1].size(); ++sh)
								MatVecMult(vvGradV_Def[d1][sh], JTInv, vvLocGradV_Def[d1][sh]);
					}

					MathMatrix<dim, dim> IdGrad;
					for(int d1 = 0; d1 < dim; ++d1){
						for(int d2 = 0; d2 < dim; ++d2){
							IdGrad(d1, d2) = 0.0;
							for(size_t sh = 0; sh < vvValueDef[d1].size(); ++sh){
								IdGrad(d1, d2) += vvValueDef[d1][sh] * vvGradV_Def[d1][sh][d2];
							}
						}
					}
					for(int d1=0;d1<dim; ++d1 ){
						IdGrad(d1,d1) += 1;
					}
					//2. Get det(DF) and Inverse(DF)
					const number detDF= Determinant(IdGrad);
					MathMatrix<dim, dim> InvDF; Inverse(InvDF,IdGrad);
					
					//3. Interpolate Velocity Gradient Matrix
					MathMatrix<dim, dim> VelGrad;
					for(int d1 = 0; d1 < dim; ++d1){
						for(int d2 = 0; d2 < dim; ++d2){
							VelGrad(d1, d2) = 0.0;
							for(size_t sh = 0; sh < vvValueVel[d1].size(); ++sh){
								VelGrad(d1, d2) += vvValueVel[d1][sh] * vvGradV[d1][sh][d2];
							}
						}
					}
					
					//4. Carry out gradient Transformations GradVel*InvDF
					MathMatrix<dim, dim> TransformedVelGrad;
					MatMultiply(TransformedVelGrad, VelGrad, InvDF);
					//5. Do the integration
					//	get quadrature weight
					const number weightIP = vWeight[ip];
					//const number weightIP = rQuadRule.weight(ip);
					//	get determinate of mapping
					const number det = SqrtGramDeterminant(vElemJT[ip]);
					
					elemDrag +=  MatContraction(TransformedVelGrad,TransformedVelGrad)* detDF * weightIP * det ;
					//elemDrag +=  MatContraction(VelGrad,VelGrad)* detDF * weightIP * det ;					
				}//end ip for
				
				drag += elemDrag;	
				
			}UG_CATCH_THROW("SumValuesOnElems failed.");//end try
		}//end element iteratior for
	}//end intSSGrp for
	
#ifdef UG_PARALLEL
	// sum over processes
	if(pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		number local = drag;
		com.allreduce(&local, &drag, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
	}
#endif
	return drag;
}//end drag


template <typename TGridFunction>
std::vector<number> BarycenterDefect(SmartPtr<TGridFunction> spGridFct,
                              const char* vCmp,
                              const char* IntegrationSubsets,                           
                              int quadOrder){
	number x_bary=0.0;
	number y_bary=0.0;
	number z_bary=0.0;
	

	static const int dim = TGridFunction::dim;
	static const int WorldDim = TGridFunction::dim;
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object grid_base_object;
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator const_iterator;
	
	std::vector<DoFIndex> multInd;
	 
	SubsetGroup intSSGrp(spGridFct->domain()->subset_handler());
	intSSGrp.add(TokenizeString(IntegrationSubsets));
	
	//	Create function group from GridFunction components
	const FunctionGroup vFctID = spGridFct->fct_grp_by_name(vCmp);
	// If every variable has its own different Local Finite Element, then each can be accessed here individually
	// usually (for instance for a deformation field), the LFEs will be the same. 
	std::vector<LFEID> vLFEID;
	//Extract and store the LFE id from the FunctionPattern through the FunctionGroup
	for(size_t fct = 0; fct < vFctID.size(); ++fct){
		vLFEID.push_back(spGridFct->lfeid(vFctID[fct]));
	}
	
	//Loop subsets assigned to the integration subset group
	for(size_t i = 0; i < intSSGrp.size(); ++i)
	{
		//Store subset index
		const int si = intSSGrp[i];
		
		const_iterator iterBegin = spGridFct->template begin<grid_base_object>(si);
		const_iterator iterEnd = spGridFct->template end<grid_base_object>(si);
		const_iterator iter = iterBegin;
		typename domain_traits<TGridFunction::dim>::position_accessor_type& aaPos= 
													spGridFct->domain()->position_accessor();
		
		typedef typename domain_traits<dim>::element_type Element;
		
		std::vector<MathVector<WorldDim> > vCorner;
		std::vector<int> vSubsetIndex;
		//For storing a Jacobian per ip for a given element
		//std::vector<MathMatrix<dim, WorldDim> > vJT;
		std::vector<MathMatrix<dim, WorldDim> > vElemJT;
		
		for(; iter != iterEnd; ++iter)
		{
			//Pointer to element obtained from dereferencing the iterator...see note above
			Element* pElem = *iter;
		
			CollectCornerCoordinates(vCorner, *pElem, aaPos, true);
			const ReferenceObjectID elemRoid = pElem->reference_object_id();
			
			try{
			//Do we need the type or use best?
				const QuadratureRule<dim>& rQuadRule = 
									QuadratureRuleProvider<dim>::get(elemRoid, quadOrder);

			//	get reference element mapping by reference object id
				DimReferenceMapping<dim, WorldDim>& rMapping = 
									ReferenceMappingProvider::get<dim, WorldDim>(elemRoid,vCorner);

			//	number of integration points
				const number* vWeight = rQuadRule.weights();
				const size_t numIP = rQuadRule.size();
				std::vector<MathVector<dim> > vLocalIP(numIP);
				std::vector<MathVector<dim> > vGlobalIP(numIP);
				
				
				
				for(size_t ip = 0; ip < numIP; ++ip)
					vLocalIP[ip]=rQuadRule.point(ip); 
				
				for(size_t ip = 0; ip < numIP; ++ip)
					rMapping.local_to_global(vGlobalIP[ip], vLocalIP[ip]);
				
				/* for(size_t ip = 0; ip < numIP; ++ip)
					rMapping.local_to_global(vGlobalIP[ip], *rQuadRule.points());  */
				
				//std::cout<<"LOCAL IPs"<<std::endl;
				//for(size_t ip = 0;ip< numIP; ++ip)
					//std::cout<<vLocalIP[ip][0]<<" "<<vLocalIP[ip][1]<<std::endl;
				
				
				/* for(size_t ip = 0; ip < numIP; ++ip)
					rMapping.local_to_global(vLocalIP[ip], rQuadRule.point(ip));
				
				for(size_t ip = 0; ip < numIP; ++ip)
					rMapping.local_to_global(vGlobalIP[ip], vLocalIP[ip]); */		
				
				
				
				
				//	get element values
				//Stores integers with dof identifiers
				std::vector<DoFIndex> vInd;
				//vector of vectors sized #functions in function group
				std::vector<std::vector<number> > vvValue(vFctID.size());
				for(size_t fct = 0; fct < vvValue.size(); ++fct){
					spGridFct->dof_indices(pElem, vFctID[fct], vInd);					
					//The dof vector (integers) now contains the indices for the corresponding function in vFctID
					//Then we resize the vector inside vvValue according to how many dofs we have in this function 
					vvValue[fct].resize(vInd.size());
					for(size_t sh = 0; sh < vInd.size(); ++sh){
						//DoFRef will give the value of GridFunction at index provided (dof)
						vvValue[fct][sh] = DoFRef(*spGridFct, vInd[sh]);
					}//end sh value for
				}//end fct value for
				
				vElemJT.resize(numIP);
				rMapping.jacobian_transposed(&(vElemJT[0]), &vLocalIP[0], numIP);
				
				//element integral value
				number integrandX = 0.0;
				number integrandY = 0.0;
				number integrandZ = 0.0;
				number detOld;
				for(size_t ip = 0; ip < numIP; ++ip)
				{
				//ip values
				number elemXBary = 0.0;
				number elemYBary = 0.0;
				number elemZBary = 0.0;
				//ip global pos
				MathVector<dim> vIP;VecSet(vIP,0.0);
				rMapping.local_to_global(vIP,vLocalIP[ip]);
					// 	1. Interpolate Functional Matrix of Identity Gradient at ip
				
				//Stores the gradients (a vector size dim ) in a matrix with dim columns
					std::vector<MathVector<dim> > vvLocGradV[dim];
					std::vector<MathVector<dim> > vvGradV[dim];
					MathMatrix<dim, dim> JTInv;
					Inverse(JTInv, vElemJT[ip]);
					for(int d1 = 0; d1 < dim; ++d1){						
						const LocalShapeFunctionSet<dim>& rTrialSpaceP = 
								LocalFiniteElementProvider::get<dim>(elemRoid, vLFEID[d1]);
						rTrialSpaceP.grads(vvLocGradV[d1], vLocalIP[ip]);
							vvGradV[d1].resize(vvLocGradV[d1].size());
							for(size_t sh = 0; sh < vvGradV[d1].size(); ++sh)
								MatVecMult(vvGradV[d1][sh], JTInv, vvLocGradV[d1][sh]);
					}
					MathMatrix<dim, dim> IdGrad;
					for(int d1 = 0; d1 < dim; ++d1){
						for(int d2 = 0; d2 < dim; ++d2){
							IdGrad(d1, d2) = 0.0;
							for(size_t sh = 0; sh < vvValue[d1].size(); ++sh){																	
								IdGrad(d1, d2) += vvValue[d1][sh] * vvGradV[d1][sh][d2];
							}
						}
					}
					for(int d1=0;d1<dim; ++d1 ){
						IdGrad(d1,d1) += 1;
					}
					const number detDF= Determinant(IdGrad);
					/* if(detDF != 1) 
						UG_THROW("VOLUME DEFECT: THE DETERMINANT OF THE DEFORMATION \n GRADIENT IS NOT 1 UNDER NO DEFORMATION CONDITION")
					*/
				//2.1- Interpolate coordinate on x, index 0
					const LocalShapeFunctionSet<dim>& rTrialSpaceP1 = 
									LocalFiniteElementProvider::get<dim>(elemRoid, vLFEID[0]);
					std::vector<number> vShapeP1;
					rTrialSpaceP1.shapes(vShapeP1, vLocalIP[ip]);

					for(size_t sh = 0; sh < vvValue[0].size(); ++sh){
						elemXBary += vShapeP1[sh] * vvValue[0][sh];//this gives W_x
					}
					elemXBary += vGlobalIP[ip][0];//w_x+ip_x
					//elemXBary += vIP[0];
				//2.2- Interpolate coordinate on y, index 2
					const LocalShapeFunctionSet<dim>& rTrialSpaceP2 = 
									LocalFiniteElementProvider::get<dim>(elemRoid, vLFEID[1]);
					std::vector<number> vShapeP2;
					rTrialSpaceP2.shapes(vShapeP2, vLocalIP[ip]);

					for(size_t sh = 0; sh < vvValue[1].size(); ++sh){
						elemYBary += vShapeP2[sh] * vvValue[1][sh];//this gives W_y
					}
					elemYBary += vGlobalIP[ip][1];
					//elemYBary += vIP[1];
				//2.3- Interpolate coordinate on z, index 2
					if(dim == 3){
						const LocalShapeFunctionSet<dim>& rTrialSpaceP3 = 
									LocalFiniteElementProvider::get<dim>(elemRoid, vLFEID[2]);
						std::vector<number> vShapeP3;
						rTrialSpaceP3.shapes(vShapeP3, vLocalIP[ip]);

						for(size_t sh = 0; sh < vvValue[2].size(); ++sh){
							elemZBary += vShapeP3[sh] * vvValue[2][sh];//this gives W_y
						}
						elemZBary += vGlobalIP[ip][2];
					}//end if dim 3
				//3.- Perform Integration (Summation)

				//	get quadrature weight
					const number weightIP = vWeight[ip];
					//std::cout<<vGlobalIP[ip][0]<<" "<<vGlobalIP[ip][1]<<" "<<weightIP<<std::endl;
				//	get determinate of mapping
					const number det = SqrtGramDeterminant(vElemJT[ip]);
					
						
					
					integrandX +=  elemXBary * detDF * weightIP * det ;
					integrandY +=  elemYBary * detDF * weightIP * det ;
					integrandZ +=  elemZBary * detDF * weightIP * det ;
					
					if(ip>0)
					{
					if(det != detOld) 
						UG_THROW("Unequal dets of JT")
					}
					
					detOld=det;
				}//end ip for
				x_bary += integrandX;
				y_bary += integrandY;
				z_bary += integrandZ;
			}UG_CATCH_THROW("SumValuesOnElems failed.");
		}//end elem for
		
	}//end intSSGrp for
#ifdef UG_PARALLEL
	// sum over processes
	if(pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		number local = x_bary;
		com.allreduce(&local, &x_bary, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
		local = y_bary;
		com.allreduce(&local, &y_bary, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
		local = z_bary;
		com.allreduce(&local, &z_bary, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
	}
#endif
	std::vector<number> vBarycenter(dim);VecSet(vBarycenter,0.0);
	vBarycenter[0] = x_bary;
	vBarycenter[1] = y_bary;
	
	if(dim == 3){
		vBarycenter[2] = z_bary;
	}
	return vBarycenter;
}//end barycenter calculation function

//template <typename TGridFunction>
//number VolumeDefect(SmartPtr<TGridFunction> spGridFct,const char* vCmp,const char* BndSubsets, const char* InnerSubsets,number kinVisco, number density,int quadOrder)
template <typename TGridFunction>
number VolumeDefect(SmartPtr<TGridFunction> spGridFct, number currentVolume, 
					const char* IntegrationSubsets, const char* vCmp, int quadOrder,
					bool noDeformation, number comp_val, bool compare){
						
	number newVolume=0.0;

	static const int dim = TGridFunction::dim;
	static const int WorldDim = TGridFunction::dim;
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object grid_base_object;
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator const_iterator;
	
	 std::vector<DoFIndex> multInd;

	//	coord and vertex array
	//MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	//Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];
	
	//	Create SubsetGroup from the given subsets for integration, if empty use all
	SubsetGroup intSSGrp(spGridFct->domain()->subset_handler());
	intSSGrp.add(TokenizeString(IntegrationSubsets));
	
/* 	if(IntegrationSubsets != NULL){
		intSSGrp.add(TokenizeString(IntegrationSubsets));
		if(!SameDimensionsInAllSubsets(intSSGrp))
			UG_THROW("VolumeDefect: Subsets '"<<IntegrationSubsets<<"' do not have same dimension."
					 "Can not integrate on subsets of different dimensions.");
	}
	else{
		intSSGrp.add_all();
		//RemoveLowerDimSubsets(intSSGrp);
	} */
	
	//	Create function group from GridFunction components
	const FunctionGroup vFctID = spGridFct->fct_grp_by_name(vCmp);
	// If every variable has its own different Local Finite Element, then each can be accessed here individually
	// usually (for instance for a deformation field), the LFEs will be the same. 
	std::vector<LFEID> vLFEID;
	//Extract and store the LFE id from the FunctionPattern through the FunctionGroup
	for(size_t fct = 0; fct < vFctID.size(); ++fct){
		vLFEID.push_back(spGridFct->lfeid(vFctID[fct]));
	}
	
	//Loop subsets assigned to the integration subset group
	for(size_t i = 0; i < intSSGrp.size(); ++i)
	{
		//Store subset index
		const int si = intSSGrp[i];
		
	//	note: this iterator is for the base elements, e.g. Face and not
	//	for the special type, e.g. Triangle, Quadrilateral
	// With these we will access every element from the subset
		const_iterator iterBegin = spGridFct->template begin<grid_base_object>(si);
		const_iterator iterEnd = spGridFct->template end<grid_base_object>(si);
		const_iterator iter = iterBegin;
		typename domain_traits<TGridFunction::dim>::position_accessor_type& aaPos= spGridFct->domain()->position_accessor();
		//const ISubsetHandler* ish = spGridFct->domain()->subset_handler().get();
		//Grid& grid = *spGridFct->domain()->grid();
		
	//	this is the base element type (e.g. Face). This is the type when the
	//	iterators above are dereferenciated.
	//  What he said...
		typedef typename domain_traits<dim>::element_type Element;
		//typedef typename domain_traits<dim>::side_type Side;

	//	vector of corner coordinates of element corners (to be filled for each elem)
		//vector of MathVectors sized dim, will store the corners of each face
		std::vector<MathVector<WorldDim> > vCorner;
		std::vector<int> vSubsetIndex;
		//For storing a Jacobian per ip for a given element
		//std::vector<MathMatrix<dim, WorldDim> > vJT;
		std::vector<MathMatrix<dim, WorldDim> > vElemJT;
		
		
		//We start to access every element using our iterators (give a read to the topic online or on c++ literature)
		//...pointer shit...STL related things...
		for(; iter != iterEnd; ++iter)
		{
			//Pointer to element obtained from dereferencing the iterator...see note above
			Element* pElem = *iter;
		
			//	get all corner coordinates
			//returns the corner coordinates of a geometric object
			//Returns the corner coordinated of a geometric object in a vector
			//This function collects the corner coordinates for a given geometric object in the order prescribed by the reference elements
			//[out]	vCornerCoordsOut	vector of corner coordinates
			//[in]	elem	Geometric Object
			//[in]	aaPos	AttachmentAccessor for Positions
			//[in]	clearContainer	empty container before filling
			CollectCornerCoordinates(vCorner, *pElem, aaPos, true);
			
			//This tells us which kind of geometric object it is, e.g. ROID_TRIANGLE 	
			const ReferenceObjectID elemRoid = pElem->reference_object_id();
			
			try{
				//Do we need the type or use best?
				const QuadratureRule<dim>& rQuadRule = QuadratureRuleProvider<dim>::get(elemRoid, quadOrder);

			//	get reference element mapping by reference object id
				DimReferenceMapping<dim, WorldDim>& rMapping = ReferenceMappingProvider::get<dim, WorldDim>(elemRoid,vCorner);

			//	number of integration points
				const number* vWeight = rQuadRule.weights();
				const size_t numIP = rQuadRule.size();
				std::vector<MathVector<dim> > vLocalIP(numIP);
				std::vector<MathVector<dim> > vGlobalIP(numIP);
			/* TODO:we are not using the sides as subelements, since we do not need normals nor tangents
			//TODO: do we have to get the sides as well? We do not need the normals...
			//pending if necessary, think about the following: it seems that this would give us
			//access to the edges and faces of the element (depending on 2d/3d)
			//	get sides
			//side could be edge or face or vertex
			typename Grid::traits<Side>::secure_container vSide;
			//Gets lower dimensional elements from the element
			grid.associated_elements_sorted(vSide, pElem);
			vSubsetIndex.resize(vSide.size());
			//TODO:think about what this does????
			for(size_t i = 0; i < vSide.size(); ++i)
				vSubsetIndex[i] = ish->get_subset_index(vSide[i]);
			*/
			
				
				//TODO:seems like this is only necessary for the sides
				//DimReferenceMapping<dim, WorldDim>& rMapping = ReferenceMappingProvider::get<dim, WorldDim>(elemRoid, vCorner);
				
				//TODO:seems like this is only necessary for the sides
				//const DimReferenceElement<dim>& rRefElem = ReferenceElementProvider::get<dim>(elemRoid);
				
				
				//	get element values
				//Stores integers with dof identifiers
				std::vector<DoFIndex> vInd;
				//vector of vectors sized #functions in function group
				std::vector<std::vector<number> > vvValue(vFctID.size());
				//We visit each vector inside the vvValue (1 for each symbolic function)
				for(size_t fct = 0; fct < vvValue.size(); ++fct){
					
					//extracts all multiindices for a function (sorted)
					//All Multi-Indices of a function living on the element (including the subelements) are extracted and stored in a std::vector.
					//The order of the indices is sorted, i.e. the dofs are provided as specified in the local dof set of the local finite element trial space. 
					//If bHang is set to true, also the DoFs on the Constrained Objects belonging to the constraining Subelements are extracted 
					//and added at the end of the indices. If bClear is set to true, the vector is cleared before insertion.
					//Parameters
					//[in]	elem	the element
					//[in]	fct	the function
					//[out]	ind	vector of multi indices
					//[in]	bHang	flag if extracting of constrained dofs required
					//[in]	bClear	flag if vector has to be clear before insertion

					spGridFct->dof_indices(pElem, vFctID[fct], vInd);
					
					//The dof vector (integers) now contains the indices for the corresponding function in vFctID
					//Then we resize the vector inside vvValue according to how many dofs we have in this function 
					vvValue[fct].resize(vInd.size());
					for(size_t sh = 0; sh < vInd.size(); ++sh){
						//DoFRef will give the value of GridFunction at index provided (dof)
						vvValue[fct][sh] = DoFRef(*spGridFct, vInd[sh]);
					}//end sh value for
				}//end fct value for
				
				
				for(size_t ip = 0; ip < numIP; ++ip)
					rMapping.local_to_global(vLocalIP[ip], rQuadRule.point(ip));
				
				//Obtain the positions of 
				for(size_t ip = 0; ip < numIP; ++ip)
					rMapping.local_to_global(vGlobalIP[ip], vLocalIP[ip]);
				
				//Transformation matrices
				//std::vector<MathMatrix<dim, WorldDim> > vElemJT(numIP);
				vElemJT.resize(numIP);
				rMapping.jacobian_transposed(&(vElemJT[0]), &vLocalIP[0], numIP);
				
				//element integral value
				number intValElem = 0;
				//loop integration points
				for(size_t ip = 0; ip < numIP; ++ip)
				{
						
				// 	1. Interpolate Functional Matrix of Identity Gradient at ip
				
				//Stores the gradients (a vector size dim ) in a matrix with dim columns
					std::vector<MathVector<dim> > vvLocGradV[dim];
					std::vector<MathVector<dim> > vvGradV[dim];
					MathMatrix<dim, dim> JTInv;
					Inverse(JTInv, vElemJT[ip]);
					for(int d1 = 0; d1 < dim; ++d1){
						
						const LocalShapeFunctionSet<dim>& rTrialSpaceP = LocalFiniteElementProvider::get<dim>(elemRoid, vLFEID[d1]);
						// grads returns returns all gradients evaluated at a several points
						// This function returns the gradients of all Shape Functions at several element-local evaluation point in an array.
						//Parameters
						//	[out]	vvGrad	Vector of gradients
						//	[in]	vLocPos	Vector of Position on reference element

						rTrialSpaceP.grads(vvLocGradV[d1], vLocalIP[ip]);

							vvGradV[d1].resize(vvLocGradV[d1].size());
							for(size_t sh = 0; sh < vvGradV[d1].size(); ++sh)
								MatVecMult(vvGradV[d1][sh], JTInv, vvLocGradV[d1][sh]);
					}

					MathMatrix<dim, dim> IdGrad;
					for(int d1 = 0; d1 < dim; ++d1){
						for(int d2 = 0; d2 < dim; ++d2){
							IdGrad(d1, d2) = 0.0;
							for(size_t sh = 0; sh < vvValue[d1].size(); ++sh){
								
								if(noDeformation ==  true){
						
									if(vvValue[d1][sh] != 0) UG_THROW("VOLUME DEFECT: vvValue != 0")
								}
								if(compare == true){
								
									if(vvValue[d1][sh] != comp_val) UG_THROW("VOLUME DEFECT: comparison failed")
								}
									
								IdGrad(d1, d2) += vvValue[d1][sh] * vvGradV[d1][sh][d2];
							}
						}
					}
					for(int d1=0;d1<dim; ++d1 ){
						IdGrad(d1,d1) += 1;
					}
					const number detDF= Determinant(IdGrad);
					
					if(noDeformation ==  true){
						
						if(detDF != 1) UG_THROW("VOLUME DEFECT: THE DETERMINANT OF THE DEFORMATION \n GRADIENT IS NOT 1 UNDER NO DEFORMATION CONDITION")
					}
					//	get quadrature weight
					const number weightIP = vWeight[ip];
					//const number weightIP = rQuadRule.weight(ip);
					//	get determinate of mapping
					const number det = SqrtGramDeterminant(vElemJT[ip]);
					
					intValElem +=  detDF * weightIP * det ;				
				}//end ip for
				
				newVolume += intValElem;	
				
			}UG_CATCH_THROW("SumValuesOnElems failed.");
			//end try
			
		}//end element iteratior for
		
		
	}//end intSSGrp for
	
#ifdef UG_PARALLEL
	// sum over processes
	if(pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		number local = newVolume;
		com.allreduce(&local, &newVolume, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
	}
#endif
	return newVolume-currentVolume;

}


}//end namespace FluidOptim
}//end namespace ug

#endif 

