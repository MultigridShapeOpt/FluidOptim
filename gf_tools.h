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

#ifndef GF_TOOLS
#define GF_TOOLS

#include <vector>
#include <string>
#include <cmath>  // for isinf, isnan
#include <boost/function.hpp>


#include "common/util/file_util.h"

#include "lib_algebra/cpu_algebra/sparsematrix_print.h"
#include "lib_algebra/operator/interface/matrix_operator.h"
#include "lib_algebra/operator/debug_writer.h"
#include "lib_algebra/operator/vector_writer.h"
#include "lib_algebra/common/matrixio/matrix_io_mtx.h"
#include "lib_algebra/common/connection_viewer_output.h"
#include "lib_algebra/common/csv_gnuplot_output.h"
#include "lib_disc/io/vtkoutput.h"
#include "lib_disc/spatial_disc/constraints/constraint_interface.h"
#include "lib_disc/dof_manager/dof_distribution.h"
#include "lib_disc/spatial_disc/user_data/user_data.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/common/geometry_util.h"
#include "lib_disc/function_spaces/integrate.h"
#include "lib_grid/algorithms/debug_util.h"  // for ElementDebugInfo
#include "lib_grid/tools/periodic_boundary_manager.h"

#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/dof_position_util.h"

#ifdef UG_PARALLEL
#include "pcl/pcl.h"
#endif

namespace ug {
//vCmp: symbolic for control variable
//vCmpAdj: symbolics for adjoint deformation L
//spGridFctAdj: n-dimensional gridfunc
//spGridFctCtrl: scalar grid function
//spGridFctOut: scalar grid function

template <typename TGridFunction>
void ControlGradient(SmartPtr<TGridFunction> spGridFctOut, SmartPtr<TGridFunction> spGridFctCtrl, SmartPtr<TGridFunction> spGridFctAdj,
                              const char* vCmp, const char* vCmpAdj, const number alpha,
                              const char* BndSubsets, const char* InnerSubsets,
                              int quadOrder){
	static const int dim = TGridFunction::dim;
	static const int WorldDim = TGridFunction::dim;
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object grid_base_object;
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator const_iterator;

//	We use the output spGridFctOut as pivot for all pointers of geomtry across the function
//	read subsets
	SubsetGroup innerSSGrp(spGridFctOut->domain()->subset_handler());
	if(InnerSubsets != NULL){
		innerSSGrp.add(TokenizeString(InnerSubsets));
		if(!SameDimensionsInAllSubsets(innerSSGrp))
			UG_THROW("ControlGradient: Subsets '"<<InnerSubsets<<"' do not have same dimension."
					 "Can not integrate on subsets of different dimensions.");
	}
	else{
		innerSSGrp.add_all();
		RemoveLowerDimSubsets(innerSSGrp);
	}

//	read subsets
	SubsetGroup bndSSGrp(spGridFctOut->domain()->subset_handler());
	if(BndSubsets != NULL)
		bndSSGrp.add(TokenizeString(BndSubsets));
	else
		UG_THROW("ControlGradient: No boundary subsets passed.");
	/*Group of functions represented by integers FunctionGroup is just a group of size_t, 
	representing some functions. The function group is based on a FunctionPattern and the
	integer represent the position of the function in the function pattern. Selection of 
	functions is best via usage of symbolic names of the functions.
	
	Each of them extracts the given functions provided in the string vCmp from their respective 
	GridFcts.
	The correspondig bFctID FunctionGroup class will have size equal to the number of elements in the string vCmp
	We need the Local Finite Element IDs for the quadrature rule eventually */
	
//	the Grid Function Group Outputs
	const FunctionGroup vFctIDOut = spGridFctOut->fct_grp_by_name(vCmp);
	std::vector<LFEID> vLFEID_out;
	for(size_t fct = 0; fct < vFctIDOut.size(); ++fct){
		vLFEID_out.push_back(spGridFctOut->lfeid(vFctIDOut[fct]));
	}
	//the Grid Function Group Control
	const FunctionGroup vFctIDCtrl = spGridFctCtrl->fct_grp_by_name(vCmp);
	std::vector<LFEID> vLFEID_ctrl;
	for(size_t fct = 0; fct < vFctIDCtrl.size(); ++fct){
		vLFEID_ctrl.push_back(spGridFctCtrl->lfeid(vFctIDCtrl[fct]));
	}
	//the Grid Function Group Adjoint
	const FunctionGroup vFctIDAdj = spGridFctAdj->fct_grp_by_name(vCmpAdj);
	std::vector<LFEID> vLFEID_adj;
	for(size_t fct = 0; fct < vFctIDAdj.size(); ++fct){
		vLFEID_adj.push_back(spGridFctAdj->lfeid(vFctIDAdj[fct]));
	}



//	loop subsets
	for(size_t i = 0; i < innerSSGrp.size(); ++i)
	{
	//	get subset index
		const int si = innerSSGrp[i];

	//	skip empty subset
		if(innerSSGrp.dim(i) == DIM_SUBSET_EMPTY_GRID) continue;

	//	check dimension
		if(innerSSGrp.dim(i) != dim)
			UG_THROW("ControlGradient: Dimension of inner subset is "<<
					 innerSSGrp.dim(i)<<", but only World Dimension "<<dim<<
					 " subsets can be used for inner subsets.");

	//	note: this iterator is for the base elements, e.g. Face and not
	//			for the special type, e.g. Triangle, Quadrilateral
		const_iterator iterBegin = spGridFctOut->template begin<grid_base_object>(si);
		const_iterator iterEnd = spGridFctOut->template end<grid_base_object>(si);
		const_iterator iter = iterBegin;

		typename domain_traits<TGridFunction::dim>::position_accessor_type& aaPos
			= spGridFctOut->domain()->position_accessor();
		const ISubsetHandler* ish = spGridFctOut->domain()->subset_handler().get();
		Grid& grid = *spGridFctOut->domain()->grid();

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
/*
			std::cout<<"Size of vFctAdj is: "<<vFctIDAdj.size()<<"\n";
				std::cout<<"Size of vFctAdj is: "<<vFctIDAdj[0]<<","<<vFctIDAdj[1]<<"\n";
				//Get values of adjoint deformation in vvValueAdj
				std::vector<DoFIndex> vIndAdj;
				std::vector<std::vector<number> > vvValueAdj(vFctIDAdj.size());
				for(size_t fct = 0; fct < vvValueAdj.size(); ++fct){
					spGridFctAdj->dof_indices(pElem, vFctIDAdj[fct], vIndAdj);
					std::cout<<"Size of vIndAdj is: "<<vIndAdj.size()<<"\n";
					vvValueAdj[fct].resize(vIndAdj.size());
						for(size_t sh = 0; sh < vIndAdj.size(); ++sh){
							std::cout<<"Value from dofref is: "<<DoFRef(*spGridFctAdj, vIndAdj[sh])<<"\n";
							vvValueAdj[fct][sh] = DoFRef(*spGridFctAdj, vIndAdj[sh]);
						}
				}
		
*/
		//	loop sides of element
			for(size_t side = 0; side < vSide.size(); ++side)
			{
			//	check if side used
				if(!bndSSGrp.contains(vSubsetIndex[side])) continue;

			//	get side
				Side* pSide = vSide[side];
				
				/*
				Store values of the GridFunc for the current geometric object (side/element). The vv stores them as follows:
					#rows: number of functions in FunctionGroup
					#columns: number of dofs of the given function on the geometric object
				The vector of DoFIndex vInd will have size equal to the number of dofs for the given function in the FunctionGroup. On each
				loop DoFIndex vector is overwritten. It stores the indeces as ints, which stems from MultiIndex class.
				Is as follows:
					-create Dof vector to be filled
					-create vv for storing the actual values
					-function dof_indices(geom_object, function, dof_vector) fills dof_vector accordingly
					-resize the respective row of the vv to dof_vector.size() entries (columns) which will store the values
					-for each entry of this row (column entries) correspondig to one function use DoFRef(GridFunc, dof_vector)
					 to fill the vv of values.
				The idx sh represents the node (where a shape function might be), when linear these are the corners of the element
				*/
				std::cout<<"Size of FunctionGroup vFctAdj is: "<<vFctIDAdj.size()<<"with elements ["<<vFctIDAdj[0]<<","<<vFctIDAdj[1]<<"]\n";
										
				std::vector<DoFIndex> vIndAdj;
				std::vector<std::vector<number> > vvValueAdj(vFctIDAdj.size());
				for(size_t fct = 0; fct < vvValueAdj.size(); ++fct){
					spGridFctAdj->dof_indices(pSide, vFctIDAdj[fct], vIndAdj);
					
					std::cout<<"Size of vIndAdj is: "<<vIndAdj.size()<<"\n";//verify that vIndAdj got assigned dof indices
					
					vvValueAdj[fct].resize(vIndAdj.size());
						for(size_t sh = 0; sh < vIndAdj.size(); ++sh){
							std::cout<<"Value from dofref is: "<<DoFRef(*spGridFctAdj, vIndAdj[sh])<<"\n";//see value gotten
							vvValueAdj[fct][sh] = DoFRef(*spGridFctAdj, vIndAdj[sh]);
						}
				}
				
				//for Grid Function Ctrl
				std::vector<DoFIndex> vIndCtrl;
				std::vector<std::vector<number> > vvValueCtrl(vFctIDCtrl.size());
				for(size_t fct = 0; fct < vvValueCtrl.size(); ++fct){
					spGridFctCtrl->dof_indices(pSide, vFctIDCtrl[fct], vIndCtrl);
					vvValueCtrl[fct].resize(vIndCtrl.size());
					for(size_t sh = 0; sh < vIndCtrl.size(); ++sh)
						vvValueCtrl[fct][sh] = DoFRef(*spGridFctCtrl, vIndCtrl[sh]);
				}
				
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
			
			// 	We have to decide which NORMAL to use, NOW: INNER normal
				MathVector<WorldDim> Normal;
				ElementNormal<WorldDim>(sideRoid, Normal, &vSideCorner[0]);
				VecNormalize(Normal, Normal);
				VecScale(Normal, Normal, 1); // inner normal
				
				std::cout<<"The normal Vector is: "<<Normal[0]<<","<<Normal[1]<<"\n";//check that it works
				
			//	quadrature points
				const number* vWeight = rSideQuadRule.weights();
				const size_t nip = rSideQuadRule.size();
				std::vector<MathVector<dim> > vLocalIP(nip);
				std::vector<MathVector<dim> > vGlobalIP(nip);

				DimReferenceMapping<dim-1, dim>& map
					= ReferenceMappingProvider::get<dim-1, dim>(sideRoid, vLocalSideCorner);

				for(size_t ip = 0; ip < nip; ++ip)
					map.local_to_global(vLocalIP[ip], rSideQuadRule.point(ip));
				
				for(size_t ip = 0; ip < nip; ++ip)
					rMapping.local_to_global(vGlobalIP[ip], vLocalIP[ip]);

			//	compute transformation matrices
				DimReferenceMapping<dim-1, dim>& map2
					= ReferenceMappingProvider::get<dim-1, dim>(sideRoid, vSideCorner);
				std::vector<MathMatrix<dim-1, WorldDim> > vJT(nip);
				map2.jacobian_transposed(&(vJT[0]), rSideQuadRule.points(), nip);

				std::vector<MathMatrix<dim, WorldDim> > vElemJT(nip);
				rMapping.jacobian_transposed(&(vElemJT[0]), &vLocalIP[0], nip);
				
				/*
				We do a numeric integratio and attempt to modify the values of the output grid_function
				*/
				for(size_t ip = 0; ip < nip; ++ip)
				{
					/*
					Interpolation of the values of the boundary control variable, we use the FEM ID for this. 
					The trialSpace object will fill the shape function vector with the corresponding values of the ShapeFunc at
					the integration point.
					The value then is taken from the multplication between the values of the grid_function component on the corners
					*/
					const LocalShapeFunctionSet<dim>& rTrialSpaceP2 = LocalFiniteElementProvider::get<dim>(elemRoid, vLFEID_ctrl[0]);
					std::vector<number> vShapeGF2;
					rTrialSpaceP2.shapes(vShapeGF2, vLocalIP[ip]);
					
					
					number value_ctrl = 0.0;
					for(size_t sh = 0; sh < vvValueCtrl[0].size(); ++sh)
						value_ctrl += vShapeGF2[sh] * vvValueCtrl[0][sh];
					
					//Interpolate Adjoints vector
					MathVector<dim> vLocalValues;VecSet(vLocalValues,0.0);
					for(int fct = 0; fct < dim; ++fct){
						
						const LocalShapeFunctionSet<dim>& rTrialSpaceAdj =
							LocalFiniteElementProvider::get<dim>(elemRoid, vLFEID_adj[fct]);
						std::vector<number> vShapeFunctions;
						rTrialSpaceAdj.shapes(vShapeFunctions, vLocalIP[ip]);

						for(size_t sh = 0; sh < vvValueAdj[fct].size(); ++sh)
							vLocalValues[fct] += vShapeFunctions[sh] * vvValueAdj[fct][sh];
						
					}//end for fct to fill vLocalValues
					
					//	get quadrature weight
					const number weightIP = vWeight[ip];
				//	get determinate of mapping
				
					const number det = SqrtGramDeterminant(vJT[ip]);
					std::cout<<"Currently at ip: "<<vLocalIP[ip][0]<<","<<vLocalIP[ip][1]<<"\n";
					std::cout<<"The local value at ip is: "<<vLocalValues[0]<<","<<vLocalValues[1]<<"\n";
					std::cout<<"The value of ctrl at ip is: "<<value_ctrl<<"\n";
					//perform the actual operation
					number dot_product= VecProd(vLocalValues, Normal) * det * weightIP + alpha * value_ctrl *det *weightIP;
					std::cout<<"For ip "<<ip<<" the dot prod is= "<<VecProd(vLocalValues, Normal)<<"\n";
					std::cout<<"For ip "<<ip<<" the value dot_prod is= "<<dot_product<<"\n";
					//assign values
					//get dofs indices from side
					std::vector<DoFIndex> vIndSide;
					std::cout<<"The size of vIndSide are: "<<vIndSide.size()<<"\n";
					spGridFctOut->dof_indices(pSide, 0, vIndSide);
					for(size_t dof = 0; dof < vIndSide.size(); ++dof)
					{
						number value_to_assign= DoFRef(*spGridFctOut, vIndSide[dof])-dot_product;
						std::cout<<"Value to substract from "<<DoFRef(*spGridFctOut, vIndSide[dof])<<"\n";
						std::cout<<"Value to assign "<<value_to_assign<<"\n";
						DoFRef(*spGridFctOut, vIndSide[dof]) = value_to_assign;
					}
				}
				
			//TODO: DEPRECATED here we copy one-to-one the dof values of d -> g, they should have same dof distribution, and only 1 function idx 0
				/*
					the dofs sometimes correspond to the sh which correspond to the corners...
					for the same element, the 
				*/
	/*			std::vector<DoFIndex> vIndCopy;
				spGridFctCtrl->dof_indices(pSide, 0, vIndCopy);//TODO:inner_dof_indices or dof_indices??????
				std::cout<<"Size of assigning index group: "<<vIndCopy.size()<<"\n";
				//std::cout<<"Size of assigning index group: "<<nInd<<"\n";
				for(size_t dof = 0; dof < vIndCopy.size(); ++dof)
				{
					std::cout<<"Value to modify "<<DoFRef(*spGridFctOut, vIndCopy[dof])<<" with "
								<< alpha * DoFRef(*spGridFctCtrl, vIndCopy[dof])<<"\n";	
					number value_to_assign = DoFRef(*spGridFctOut, vIndCopy[dof])-alpha*DoFRef(*spGridFctCtrl, vIndCopy[dof]);
					std::cout<<"Value to assign "<<value_to_assign<<"\n";
					DoFRef(*spGridFctOut, vIndCopy[dof]) = value_to_assign;
				}
				*/
			} // end bf
		} // end elem
	} // end subsets

}//end controlgradient



template <typename TBaseElem, typename TGridFunction>
static void LimitGFOnElements
(ConstSmartPtr<DoFDistribution> dd, SmartPtr<TGridFunction> vecOut, ConstSmartPtr<TGridFunction> vecIn,
 number upper_limit, number lower_limit){
	 
	typename DoFDistribution::traits<TBaseElem>::const_iterator iter, iterEnd;
	// loop the elements in the subset
	std::vector<DoFIndex> vInd;
	try
	{
		// iterate all elements (including SHADOW_RIM_COPY!)
		iter = dd->template begin<TBaseElem>(SurfaceView::ALL);
		iterEnd = dd->template end<TBaseElem>(SurfaceView::ALL);
		//element loop
		for (; iter != iterEnd; ++iter)
		{
			TBaseElem* elem = *iter;
			// loop indices at this element
			const size_t nInd = vecIn->inner_dof_indices(elem, 0, vInd);
			for (size_t dof = 0; dof < nInd; ++dof)
			{
				const number& val = DoFRef(*vecIn, vInd[dof]);
								
				if (val < lower_limit)
				{
					DoFRef(*vecOut, vInd[dof]) = lower_limit;
				}
				else if (val > upper_limit){
					DoFRef(*vecOut, vInd[dof]) = upper_limit;
				}else{
					DoFRef(*vecOut, vInd[dof]) = DoFRef(*vecIn, vInd[dof]);
				}
			}
			
		}
	}
	UG_CATCH_THROW("Error while limiting vector.")
}


template <typename TGridFunction>
void LimitGF
(SmartPtr<TGridFunction> limitedVecOut, ConstSmartPtr<TGridFunction> vecIn, 
 number upper_limit, number lower_limit)
{
	ConstSmartPtr<DoFDistribution> dd = vecIn->dof_distribution();

	if (dd->max_dofs(VERTEX))
		LimitGFOnElements<Vertex, TGridFunction>(dd, limitedVecOut, vecIn, upper_limit, lower_limit);
	if (dd->max_dofs(EDGE))
		LimitGFOnElements<Edge, TGridFunction>(dd, limitedVecOut, vecIn, upper_limit, lower_limit);
	if (dd->max_dofs(FACE))
		LimitGFOnElements<Face, TGridFunction>(dd, limitedVecOut, vecIn, upper_limit, lower_limit);
	if (dd->max_dofs(VOLUME))
		LimitGFOnElements<Volume, TGridFunction>(dd, limitedVecOut, vecIn, upper_limit, lower_limit);
}

template <typename TGridFunction>
number L2VecProd(TGridFunction* spGridFct1,
				TGridFunction* spGridFct2,
									const char* vCmp,
									const char* IntegrationSubsets,                           
									int quadOrder){
	number scalar_product=0.0;//final result

	const static int _C_ = 0;
	static const int dim = TGridFunction::dim;
	static const int WorldDim = TGridFunction::dim;
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object grid_base_object;
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator const_iterator;
	
	std::vector<DoFIndex> multInd;
	SubsetGroup intSSGrp(spGridFct1->domain()->subset_handler());
	intSSGrp.add(TokenizeString(IntegrationSubsets));
	
	//std::cout<<"L2VecProd: 1.- DEFINITION OF TYPES AND INITIALIZATION DONE\n";
	
	//	Create function group from GridFunction components
	const FunctionGroup vFctID_1 = spGridFct1->fct_grp_by_name(vCmp);//should be of size 1 (only 1 function)
	const FunctionGroup vFctID_2 = spGridFct2->fct_grp_by_name(vCmp);
	std::vector<LFEID> vLFEID1;
	//Extract and store the LFE id from the FunctionPattern through the FunctionGroup
	for(size_t fct = 0; fct < vFctID_1.size(); ++fct){
		vLFEID1.push_back(spGridFct1->lfeid(vFctID_1[fct]));
	}
	
	std::vector<LFEID> vLFEID2;
	//Extract and store the LFE id from the FunctionPattern through the FunctionGroup
	for(size_t fct = 0; fct < vFctID_2.size(); ++fct){
		vLFEID2.push_back(spGridFct2->lfeid(vFctID_2[fct]));
	}
	//std::cout<<"L2VecProd: 2.- CREATION OF FINITE ELEMENT IDS AND EXTRACTION OF COMPONENT DONE\n";
	//std::cout<<"L2VecProd: 3.- START OF LOOP FOR SUBSETS\n";
	//Loop subsets assigned to the integration subset group
	for(size_t i = 0; i < intSSGrp.size(); ++i)
	{
		//Store subset index
		const int si = intSSGrp[i];
		const_iterator iterBegin = spGridFct1->template begin<grid_base_object>(si);
		const_iterator iterEnd = spGridFct1->template end<grid_base_object>(si);
		const_iterator iter = iterBegin;
		typename domain_traits<TGridFunction::dim>::position_accessor_type& aaPos= spGridFct1->domain()->
																				   position_accessor();
		
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
				std::vector<DoFIndex> vInd1;
				//vector of vectors sized #functions in function group
				std::vector<std::vector<number> > vvValue1(vFctID_1.size());
				//We visit each vector inside the vvValue1 (1 for each symbolic function)
				for(size_t fct = 0; fct < vvValue1.size(); ++fct){
					spGridFct1->dof_indices(pElem, vFctID_1[fct], vInd1);
					vvValue1[fct].resize(vInd1.size());
					for(size_t sh = 0; sh < vInd1.size(); ++sh){
						//DoFRef will give the value of GridFunction at index provided (dof)
						vvValue1[fct][sh] = DoFRef(*spGridFct1, vInd1[sh]);
					}//end sh value for
				}//end fct value for
				
				//	Get Dummy GridFunction 2 values in the element
				//Stores integers with dof identifiers
				std::vector<DoFIndex> vInd2;
				//vector of vectors sized #functions in function group
				std::vector<std::vector<number> > vvValue2(vFctID_2.size());
				//We visit each vector inside the vvValue2 (1 for each symbolic function)
				for(size_t fct = 0; fct < vvValue2.size(); ++fct){
					spGridFct2->dof_indices(pElem, vFctID_2[fct], vInd2);
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
				//	1. Interpolate values for Dummy 1
					const LocalShapeFunctionSet<dim>& rTrialSpaceP1 =
							LocalFiniteElementProvider::get<dim>(elemRoid, vLFEID1[_C_]);
					std::vector<number> vShape1;
					rTrialSpaceP1.shapes(vShape1, vLocalIP[ip]);

					number value1 = 0.0;
					for(size_t sh = 0; sh < vvValue1[_C_].size(); ++sh)
						value1 +=  vShape1[sh]*vvValue1[_C_][sh];
					
				//	1. Interpolate values for Dummy 2
					const LocalShapeFunctionSet<dim>& rTrialSpaceP2 =
							LocalFiniteElementProvider::get<dim>(elemRoid, vLFEID2[_C_]);
					std::vector<number> vShape2;
					rTrialSpaceP2.shapes(vShape2, vLocalIP[ip]);

					number value2 = 0.0;
					for(size_t sh = 0; sh < vvValue2[_C_].size(); ++sh)
						value2 += vShape2[sh] * vvValue2[_C_][sh];
					
					//5. Do the integration
					//	get quadrature weight
					const number weightIP = vWeight[ip];
					//const number weightIP = rQuadRule.weight(ip);
					//	get determinate of mapping
					const number det = SqrtGramDeterminant(vElemJT[ip]);
					
					elemValue += value1*value2 * weightIP * det ;					
				}//end ip for
				
				scalar_product += elemValue;	
				
			}UG_CATCH_THROW("VecProd failed.");//end try
		}//end element iteratior for
	}//end intSSGrp for
	//std::cout<<"L2VecProd: 5.- END OF LOOP FOR SUBSETS AND RETURN\n";
#ifdef UG_PARALLEL
	// sum over processes
	if(pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		number local = scalar_product;
		com.allreduce(&local, &scalar_product, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
	}
#endif 
	//std::cout<<"L2VecProd: 6.- RESULT OF DOT PRODUCT: "<<scalar_product<<"\n";
	return scalar_product;
	
}//end L2VecProd




} // end namespace ug

#endif /* GF_TOOLS */
