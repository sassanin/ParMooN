/** ==========================================================================
#    This file is part of the finite element software ParMooN.
# 
#    ParMooN (cmg.cds.iisc.ac.in/parmoon) is a free finite element software  
#    developed by the research groups of Prof. Sashikumaar Ganesan (IISc, Bangalore),
#    Prof. Volker John (WIAS Berlin) and Prof. Gunar Matthies (TU-Dresden):
#
#    ParMooN is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as
#    published by the Free Software Foundation, either version 3 of the
#    License, or (at your option) any later version.
#
#    ParMooN is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with ParMooN.  If not, see <http://www.gnu.org/licenses/>.
#
#    If your company is selling a software using ParMooN, please consider 
#    the option to obtain a commercial license for a fee. Please send 
#    corresponding requests to sashi@iisc.ac.in

# =========================================================================*/ 
   
// =======================================================================
// @(#)Convolution.C        1.3 04/13/00
//
// Purpose:     convolute velocity and tensors
//
// Authors:     Volker John, Gunar Matthies
// =======================================================================

#include <Convolution.h>
#include <Database.h>
#ifdef __2D__
  #include <FEDatabase2D.h>
  #include <NodalFunctional2D.h>
#endif
#ifdef __3D__
  #include <FEDatabase3D.h>
  #include <NodalFunctional3D.h>
#endif
#include <Joint.h>
#include <LinAlg.h>
#include <MooNMD_Io.h>
#include <math.h>
#include <stdlib.h>

#include <string.h>

// compute value of Gaussian filter

// compute characteristic filter width
double CharacteristicFilterWidth(double h)
{
  double delta;
  double filter_constant = TDatabase::ParamDB->FILTER_WIDTH_CONSTANT;
  double filter_power = TDatabase::ParamDB->FILTER_WIDTH_POWER;

  if (filter_power==1)
    delta = filter_constant *h;
  else
    delta =  filter_constant*pow(h,filter_power);

  return(delta);  
}
#ifdef __2D__

double GaussianFilter(double delta, double dist_sq)
{
  double a;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;
  
  a = (gamma/Pi)*exp(-gamma * dist_sq/(delta*delta))/(delta*delta);
  return(a);
}

void  ConvoluteVelocity(TFEVectFunct2D *u, TFEVectFunct2D *uConv)
{
  int i,j,k,l,n,m, N_UsedElements, N_LocalUsedElements;
  int N_LocalUsedElementsNeigh;
  int N_Cells, N_Points, N_Parameters, N_, N_U, N_DOF, N_loc_dofConv;
  int Used[N_FEs2D], UsedNeigh[N_FEs2D], *N_BaseFunct;
  int neigh_i, N_Neigh, ll, lll, N_PointsNeigh, N_Edges, same_dof;
  FE2D LocalUsedElements[N_FEs2D], CurrentElement, CurrentElementConv;
  FE2D CurrentElementNeigh,LocalUsedElementsNeigh[N_FEs2D] ;
  TFE2D *Element;
  TNodalFunctional2D *nf;
  BaseFunct2D BaseFunct, *BaseFuncts, BaseFunctNeigh;
  TCollection *Coll;
  TBaseCell *cell, *neigh;
  TFE2D *ele;
  double *weights, *xi, *eta, *xi_ref, *eta_ref;
  double *weightsNeigh, *xiNeigh, *etaNeigh;
  double X_orig[MaxN_PointsForNodal2D], Y_orig[MaxN_PointsForNodal2D];
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D], AbsDetjkConv[MaxN_QuadPoints_2D];
  double XNeigh[MaxN_QuadPoints_2D], YNeigh[MaxN_QuadPoints_2D];
  double AbsDetjkNeigh[MaxN_QuadPoints_2D];
  RefTrans2D RefTrans;
  double *Param[MaxN_QuadPoints_2D], *aux;
  double FEValue_u1_QuadPoint[MaxN_QuadPoints_2D];
  double FEValue_u2_QuadPoint[MaxN_QuadPoints_2D];
  double FEValue_u1_QuadPointNeigh[MaxN_QuadPoints_2D];
  double FEValue_u2_QuadPointNeigh[MaxN_QuadPoints_2D];
  double *ExactVal[MaxN_QuadPoints_2D];
  int *DOF, ActiveBound, DirichletBound, end, last, *DOFConv, *DOFNeigh;
  double **OrigFEValues, *Orig, value, value1;
  double **OrigFEValuesNeigh;
  double FEFunctValues[MaxN_BaseFunctions2D];
  double FEFunctValues1[MaxN_BaseFunctions2D];
  double FEFunctValuesNeigh[MaxN_BaseFunctions2D];
  double FEFunctValues1Neigh[MaxN_BaseFunctions2D];
  int *GlobalNumbers, *BeginIndex, *GlobalNumbersConv, *BeginIndexConv;
  double LocError[4];
  double hK,delta;
  bool *SecondDer;
  int n_fespaces = 1, N_Derivatives = 1, index;
  int N_UConv, N_DOFConv;
  double *Values, *u_conv,distance_sq,h_fine,g, *u_values;
  double *x_conv, *y_conv;
  TJoint *joint;
  TFESpace2D *fespace, *fespaceConv; 

  // gives an array where the needed derivatives are described
  // defined in NavierStokes.h 
  // here we get as result {D00}
  MultiIndex2D NeededDerivatives[1] = { D00 };

  // get fe spaces of velocity u and of covoluted tensor
  fespace = u->GetFESpace2D();
  fespaceConv = uConv->GetFESpace2D();

  // gives a pointer to all available basis functions which
  // are described in the fedatabase
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  // gives the number of available basis functions from the
  // fedatabase
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  // say that we don't need second derivatives
  // first: allocate memory for an array of length n_fespaces 
  // second: set all entries to FALSE
  SecondDer = new bool[n_fespaces];
  for(i=0;i<n_fespaces;i++)
    SecondDer[i] = FALSE;

  // get information of the numbering of the degrees of freedom
  // first: get pointer to array where a global numbering is stored
  GlobalNumbers = fespace->GetGlobalNumbers();
  // second: get pointer to array where the start of numbering for
  // each cell is stored
  BeginIndex = fespace->GetBeginIndex();
 
  // get number of dof for one velocity component N_U
  // -> total number of dof is twice N_U
  N_U = u->GetLength();
  N_DOF = 2*N_U;
   
  // # dof in the new fe space
  N_UConv = fespaceConv->GetN_DegreesOfFreedom();
 
  N_DOFConv = 2*N_UConv;
  // allocate memory for values of convuluted function
  u_conv = uConv->GetValues(); 
  // initialize u_conv to 0
  memset(u_conv,0,N_DOFConv*SizeOfDouble);

  // get information of the numbering of the dof for the space of convolution
  // first: get pointer to array where a global numbering is stored
  GlobalNumbersConv = fespaceConv->GetGlobalNumbers();
  // second: get pointer to array where the start of numbering for
  // each cell is stored
  BeginIndexConv = fespaceConv->GetBeginIndex();

  // allocate memory for geometric position of nodes
  x_conv = new double[N_DOFConv]; 
  y_conv = new double[N_DOFConv]; 
   
// ########################################################################
// computing the convolution
// loop over all cells
// ########################################################################

  // get pointer to set of mesh cells which define the fe space
  Coll = fespace->GetCollection(); 
  // get number of mesh cells
  N_Cells = Coll->GetN_Cells();

  // direct a pointer to the values of the finite element vector
  // function
  // first part of the array : u1
  // second part of the array : u2
  Values = u->GetValues();           

  // store positition of each mesh cells in the collection of mesh 
  // cells: set ClipBoard to i
  for(i=0;i<N_Cells;i++)                       
  {
    cell = Coll->GetCell(i);
    cell->SetClipBoard(i);
  }

  // loop over all mesh cells to accumulate the values for the 
  // convolution 
  for(i=0;i<N_Cells;i++)
  {
    // get current mesh cell
    cell = Coll->GetCell(i);        
    // get cell diameter (longest distance between to vertices of the cell)
    hK = cell->GetDiameter();       
    delta = CharacteristicFilterWidth(hK);
 
    // ####################################################################
    // find local used finite elements on this cell
    // ####################################################################
    
    // have a look how much 1 are in Used ->   N_LocalUsedElements
    CurrentElement = fespace->GetFE2D(i, cell);
    LocalUsedElements[0] = CurrentElement;
    N_LocalUsedElements = 1;
    
    // ####################################################################
    // calculate values on original element
    // ####################################################################
    
    // From id of fe space the situation in the referenz mesh cell is 
    // known. Now, get information in the original mesh cell. 
    // Every finite element gets automatically a quadrature rule.
    // input: N_LocalUsedElements, LocalUsedElements, cell, SecondDer
    // output : N_Points - number of quadrature points
    //          xi - xi-values of quad. points in ref. cell
    //          eta - eta-values of quad. points in ref. cell
    //          weigths - weights fro quad. rule in quad. points 
    //          X - X-values of quad. points in original cell  
    //          Y - Y-values of quad. points in original cell  
    //          AbsDetjk - absolute value of Jacobian in quad. points
    //          xi, eta, weights, X, Y,  AbsDetjk are pointers to arrays

    TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                         Coll, cell, SecondDer,
                         N_Points, xi, eta, weights, X, Y, AbsDetjk);

    // compute geometric positions of the fe nodes
    // get id of finite element in current mesh cell
    CurrentElement = fespace->GetFE2D(i, cell);
    CurrentElementConv = fespaceConv->GetFE2D(i, cell);
    // get fe from its id
    Element = TFEDatabase2D::GetFE2D(CurrentElementConv);
    // get reference transformation 
    RefTrans = Element->GetRefTransID();
    // get pointer to the nodal functionals (fe nodes) of the fe 
    // (in ref. cell)
    nf = Element->GetNodalFunctional2D();
    // get number and coordinates of local dof in ref cell
    // xi_ref, eta_ref are pointers
    nf->GetPointsForAll(N_loc_dofConv, xi_ref, eta_ref);
      
    // get coordinates of fe nodes in original cell
    // input: RefTrans, N_loc_dof, xi_ref, eta_ref
    // output : X_orig, Y_orig - pointers to arrays with coordinates
    //          AbsDetjk - same as above    
    TFEDatabase2D::GetOrigFromRef(RefTrans,N_loc_dofConv, xi_ref, 
                                  eta_ref, X_orig, Y_orig, AbsDetjkConv);   
   
    // get pointer to basis functions of the current element
    BaseFunct = BaseFuncts[CurrentElement];  
    // # of basis functions, is the same as N_loc_dof
    N_ = N_BaseFunct[CurrentElement];        

    // find the part of the global index array where the information
    // for the current mesh cell are stored: DOF
    DOF = GlobalNumbers + BeginIndex[i];
    DOFConv = GlobalNumbersConv + BeginIndexConv[i];
    // copy the values of the finite element function from the global
    // array (Values) to local arrays 
    for(l=0;l<N_;l++)
    {
      FEFunctValues[l] = Values[DOF[l]];
      FEFunctValues1[l] = Values[DOF[l]+N_U];
    }

    // compute values of fe function in the quadrature points with
    // linear combinations of the local fe values and the basis functions
    for(k=0;k<N_Derivatives;k++)
    {
      // get fe values of the basis functions in the quad. points
      // in the current element (pointer to pointer)
      OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct,
                                            NeededDerivatives[k]);
      // accumulate the values of the fe function 
      // do it for all quad. points
      for(j=0;j<N_Points;j++)
      {
        // get pointer to array where the values for quad. point j are
        // stored
        Orig = OrigFEValues[j];
        // initialize 
        value = value1 = 0;
        // compute linear combination over all local dof
        for(l=0;l<N_;l++)
        {
          // value += u1[l] * phi[l]
          value += FEFunctValues[l] * Orig[l];
          value1 += FEFunctValues1[l] * Orig[l];
        } // endfor l
        // copy result to position in array 
        FEValue_u1_QuadPoint[j] = value;
        FEValue_u2_QuadPoint[j] = value1;
      } // endfor j
    } // endfor k
    
    // compute part of the convolution in all local fe nodes (dof)
    // which comes from current mesh cell 
    // loop over all local dof
    for(l=0;l<N_loc_dofConv;l++)         
    {
      // get global indes of local dof
      index = DOFConv[l];          
      // initialize
      value = value1 = 0;
      // loop over all quadrature points 
      for(j=0;j<N_Points;j++) 
      {
        // square of the distance between local dof and quad. point
        // in original cell 
        distance_sq = ((X_orig[l]-X[j])*(X_orig[l]-X[j])+
                       (Y_orig[l]-Y[j])*(Y_orig[l]-Y[j])); 
        // compute value of GaussianFilter
        // input: delta, distance_sq
        g = GaussianFilter(delta,distance_sq);  
        // apply quadrature rule in ref. cell 
        value += g*AbsDetjk[j]*weights[j]*FEValue_u1_QuadPoint[j];
        value1 += g*AbsDetjk[j]*weights[j]*FEValue_u2_QuadPoint[j];
      }
      // add local result to global array
      u_conv[index] += value;
      u_conv[index+N_UConv] += value1;
      // store global coordinates of the fe nodes (dof)
      x_conv[index] = X_orig[l];
      y_conv[index] = Y_orig[l];
      x_conv[index+N_UConv] = X_orig[l];
      y_conv[index+N_UConv] = Y_orig[l];
      
      // starting to compute the integral in the neighbour cells
      // number of edges
      N_Edges=cell->GetN_Edges();
      // loop over all edges
      for(ll=0;ll<N_Edges;ll++)                           
      {
        // get pointer to edge[l]
        joint=cell->GetJoint(ll); 
        
        // if boundary edge continue 
        if ((joint->GetType() == BoundaryEdge)||
            (joint->GetType() == IsoBoundEdge))
          continue;
        
        // get pointer to neighbour cell
        neigh=cell->GetJoint(ll)->GetNeighbour(cell);
        
        // get id of neighbour
        neigh_i = neigh->GetClipBoard();
        
        // get id of finite element in neigh mesh cell wrt the
        // finite element space of the convolution
        CurrentElementNeigh = fespaceConv->GetFE2D(neigh_i, neigh);
        
        // get pointer to basis functions of the neighbour element
        BaseFunctNeigh = BaseFuncts[CurrentElementNeigh]; 
 
        // # of basis functions, is the same as N_loc_dof
        N_Neigh = N_BaseFunct[CurrentElementNeigh]; 
     
        // find the part of the global index array where the information
        // for the neigh mesh cell are stored: DOF
        DOFNeigh = GlobalNumbersConv + BeginIndexConv[neigh_i];
        
        // check if computation on neigh cell necessary
        // for all dof
        same_dof = 0;
        for (j=0;j<N_Neigh;j++)
        {
          //cout << "N " << DOFNeigh[j] << " ind " << index << endl;
          if (DOFNeigh[j]==index)
          same_dof = 1;
        }
        // if the cells have the same dof, continue
        if (same_dof==1)
          continue;
        
        ////////////////////////////////////////////////////////////////////
        
        // if not, do something! 
        // i.e. the same thing as before
        
        // get id of finite element in neigh mesh cell wrt the finite
        // element space of the velocity
        CurrentElementNeigh = fespace->GetFE2D(neigh_i, neigh);
        
        // get pointer to basis functions of the neighbour element
        BaseFunctNeigh = BaseFuncts[CurrentElementNeigh]; 
 
        // # of basis functions, is the same as N_loc_dof
        N_Neigh = N_BaseFunct[CurrentElementNeigh]; 
     
        // find the part of the global index array where the information
        // for the neigh mesh cell are stored: DOF
        DOFNeigh = GlobalNumbers + BeginIndex[neigh_i];

        // initialize auxiliary array UsedNeigh
        memset(UsedNeigh, 0, N_FEs2D*SizeOfInt);
        
        // get id of finite element space in neigh cell 
        CurrentElementNeigh = fespace->GetFE2D(neigh_i, neigh);
        // set entry in auxiliary array to 1
        UsedNeigh[CurrentElementNeigh] = 1;
       
        // have a look how much 1 are in UsedNeigh -> N_LocalUsedElementsNeigh
        N_LocalUsedElementsNeigh = 0;
        memset(LocalUsedElementsNeigh, 0, SizeOfInt*N_FEs2D);
        j = 0;
        for(k=0;k<N_FEs2D;k++)
          if(UsedNeigh[k])
          {
            LocalUsedElementsNeigh[j] = (FE2D)k;
            j++;
          }
        N_LocalUsedElementsNeigh = j;
        
        // ####################################################################
        // calculate values on original element
        // ####################################################################
        
        // From id of fe space the situation in the referenz mesh cell is 
        // known. Now, get information in the original mesh cell. 
        // Every finite element gets automatically a quadrature rule.
        // input: N_LocalUsedElementsNeigh, LocalUsedElementsNeigh, 
        //        neigh, SecondDer
        // output : N_PointsNeigh - number of quadrature points
        //          xiNeigh - xi-values of quad. points in ref. cell
        //          etaNeigh - eta-values of quad. points in ref. cell
        //          weigthsNeigh - weights fro quad. rule in quad. points 
        //          XNeigh - X-values of quad. points in original cell  
        //          YNeigh - Y-values of quad. points in original cell  
        //          AbsDetjkNeigh - absolute value of Jacobian in quad. points
        //          xiNeigh, etaNeigh, weightsNeigh, XNeigh, YNeigh,  
        //          AbsDetjkNeigh are pointers to arrays

        TFEDatabase2D::GetOrig(N_LocalUsedElementsNeigh, LocalUsedElementsNeigh, 
                             Coll, neigh, SecondDer,
                             N_PointsNeigh, xiNeigh, etaNeigh, weightsNeigh,
                             XNeigh, YNeigh, AbsDetjkNeigh);
    
        // copy the values of the finite element function from the global
        // array (Values) to local arrays 
        for(lll=0;lll<N_Neigh;lll++)
        {
          FEFunctValuesNeigh[lll] = Values[DOFNeigh[lll]];
          FEFunctValues1Neigh[lll] = Values[DOFNeigh[lll]+N_U];
        }
        
        // compute values of fe function in the quadrature points with
        // linear combinations of the local fe values and the basis functions
        for(k=0;k<N_Derivatives;k++)
        {
          // get fe values of the basis functions in the quad. points
          // in the current element (pointer to pointer)
          OrigFEValuesNeigh = TFEDatabase2D::GetOrigElementValues(BaseFunctNeigh,
                                                       NeededDerivatives[k]);
          // accumulate the values of the fe function 
          // do it for all quad. points
          for(j=0;j<N_PointsNeigh;j++)
          {
            // get pointer to array where the values for quad. point j are
            // stored
            Orig = OrigFEValuesNeigh[j];
            // initialize 
            value = value1 = 0;
            // compute linear combination over all local dof
            for(lll=0;lll<N_Neigh;lll++)
            {
              // value += u1[l] * phi[l]
              value += FEFunctValuesNeigh[lll] * Orig[lll];
              value1 += FEFunctValues1Neigh[lll] * Orig[lll];
            } // endfor lll
            // copy result to position in array 
            FEValue_u1_QuadPointNeigh[j] = value;
            FEValue_u2_QuadPointNeigh[j] = value1;
          } // endfor j
        } // endfor k
        
        // compute part of the convolution in all local fe nodes (dof)
        // which comes from current mesh cell 
        // loop over all local dof
        
        // initialize
        value = value1 = 0;
        // loop over all quadrature points 
        for(j=0;j<N_PointsNeigh;j++) 
        { 
          // square of the distance between local dof and quad. point
          // in original cell 
          switch(TDatabase::ParamDB->INTERNAL_PERIODIC_IDENTITY)
          {
            case 0:
              distance_sq = ((X_orig[l]-XNeigh[j])*(X_orig[l]-XNeigh[j])+
                             (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j])); 
              break;
            case 1:
              exit (4711);
              break;
            case 2:
              // periodic square (1,1)^2 
                if (fabs(X_orig[l]-XNeigh[j])> 1.5)
                {
                  if (X_orig[l]-XNeigh[j]>1.5)
                    distance_sq = ((X_orig[l]-XNeigh[j]-2)*(X_orig[l]-XNeigh[j]-2)+
                                   (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j])); 
                  else
                    distance_sq = ((X_orig[l]-XNeigh[j]+2)*(X_orig[l]-XNeigh[j]+2)+
                                   (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j]));                                       
                }
                else
                  distance_sq = ((X_orig[l]-XNeigh[j])*(X_orig[l]-XNeigh[j])+
                                 (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j]));                   
              break;
          }
          // square of the distance between local dof and quad. point
          // in original cell 
          //distance_sq = ((X_orig[l]-XNeigh[j])*(X_orig[l]-XNeigh[j])+
          //                (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j])); 
          //cout << "X " << XNeigh[j] << " Y " << YNeigh[j] << endl;
          // compute value of GaussianFilter
          // input: delta, distance_sq
          g = GaussianFilter(delta,distance_sq);
          //cout << "dist " <<  sqrt(distance_sq) << endl;
          //
          // apply quadrature rule in ref. cell 
          value += g*AbsDetjkNeigh[j]*weightsNeigh[j]
                        *FEValue_u1_QuadPointNeigh[j];
          value1 += g*AbsDetjkNeigh[j]*weightsNeigh[j]
                        *FEValue_u2_QuadPointNeigh[j];
            
        }
        // add local result to global array
        u_conv[index] += value;
        u_conv[index+N_UConv] += value1; 
      } // end of ll loop (edges)
    } // end of l loop (local dof)
  } // endfor i (cells)
  
 
//  for (i=0;i<N_DOFConv; i++)
  //{
//    cout << i<< " x " <<   x_conv[i] << " y  " << y_conv[i] ;
  //  OutPut(" u_conv " << u_conv[i] << endl);
  //}

  // release memory
  delete SecondDer;
  delete x_conv;
  delete y_conv;

} // ConvoluteVelocity

void  ConvoluteVelocityFull(TFEVectFunct2D *u, TFEVectFunct2D *uConv)
{
  int i,j,k,l,n,m, N_UsedElements, N_LocalUsedElements;
  int N_LocalUsedElementsNeigh;
  int N_Cells, N_Points, N_Parameters, N_, N_U, N_DOF, N_loc_dofConv;
  int Used[N_FEs2D], UsedNeigh[N_FEs2D], *N_BaseFunct;
  int neigh_i, N_Neigh, ll, lll, N_PointsNeigh, N_Edges, same_dof;
  FE2D LocalUsedElements[N_FEs2D], CurrentElement, CurrentElementConv;
  FE2D CurrentElementNeigh,LocalUsedElementsNeigh[N_FEs2D] ;
  TFE2D *Element;
  TNodalFunctional2D *nf;
  BaseFunct2D BaseFunct, *BaseFuncts, BaseFunctNeigh;
  TCollection *Coll;
  TBaseCell *cell, *neigh;
  TFE2D *ele;
  double *weights, *xi, *eta, *xi_ref, *eta_ref;
  double *weightsNeigh, *xiNeigh, *etaNeigh;
  double X_orig[MaxN_PointsForNodal2D], Y_orig[MaxN_PointsForNodal2D];
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D], AbsDetjkConv[MaxN_QuadPoints_2D];
  double XNeigh[MaxN_QuadPoints_2D], YNeigh[MaxN_QuadPoints_2D];
  double AbsDetjkNeigh[MaxN_QuadPoints_2D];
  RefTrans2D RefTrans;
  double *Param[MaxN_QuadPoints_2D], *aux;
  double FEValue_u1_QuadPoint[MaxN_QuadPoints_2D];
  double FEValue_u2_QuadPoint[MaxN_QuadPoints_2D];
  double FEValue_u1_QuadPointNeigh[MaxN_QuadPoints_2D];
  double FEValue_u2_QuadPointNeigh[MaxN_QuadPoints_2D];
  double *ExactVal[MaxN_QuadPoints_2D];
  int *DOF, ActiveBound, DirichletBound, end, last, *DOFConv, *DOFNeigh;
  double **OrigFEValues, *Orig, value, value1;
  double **OrigFEValuesNeigh;
  double FEFunctValues[MaxN_BaseFunctions2D];
  double FEFunctValues1[MaxN_BaseFunctions2D];
  double FEFunctValuesNeigh[MaxN_BaseFunctions2D];
  double FEFunctValues1Neigh[MaxN_BaseFunctions2D];
  int *GlobalNumbers, *BeginIndex, *GlobalNumbersConv, *BeginIndexConv;
  double LocError[4];
  double hK, delta;
  bool *SecondDer;
  int n_fespaces = 1, N_Derivatives = 1, index;
  int N_UConv, N_DOFConv, *conv_comp;
  double *Values, *u_conv,distance_sq,h_fine,g;
  double *u_values, *x_conv, *y_conv;
  TJoint *joint;
  TVertex *vertex0;
  double x_vertex, y_vertex;
  TFESpace2D *fespace, *fespaceConv; 

  // gives an array where the needed derivatives are described
  // defined in NavierStokes.h 
  // here we get as result {D00}
  MultiIndex2D NeededDerivatives[1] = { D00 };

  // get fe spaces of velocity u and of covoluted tensor
  fespace = u->GetFESpace2D();
  fespaceConv = uConv->GetFESpace2D();

  // gives a pointer to all available basis functions which
  // are described in the fedatabase
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  // gives the number of available basis functions from the
  // fedatabase
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  // say that we don't need second derivatives
  // first: allocate memory for an array of length n_fespaces 
  // second: set all entries to FALSE
  SecondDer = new bool[n_fespaces];
  for(i=0;i<n_fespaces;i++)
    SecondDer[i] = FALSE;

  // get information of the numbering of the degrees of freedom
  // first: get pointer to array where a global numbering is stored
  GlobalNumbers = fespace->GetGlobalNumbers();
  // second: get pointer to array where the start of numbering for
  // each cell is stored
  BeginIndex = fespace->GetBeginIndex();
 
  // get number of dof for one velocity component N_U
  // -> total number of dof is twice N_U
  N_U = u->GetLength();
  N_DOF = 2*N_U;
  
  // get pointer to set of mesh cells which define the fe space
  Coll = fespace->GetCollection();
 
  // # dof in the new fe space
  N_UConv = fespaceConv->GetN_DegreesOfFreedom();
 
  N_DOFConv = 2*N_UConv;
  // allocate memory for values of convoluted function
  u_conv = uConv->GetValues(); 
  // initialize u_conv to 0
  memset(u_conv,0,N_DOFConv*SizeOfDouble);
  
  // allocate memory for values of convoluted function
  conv_comp = new int[N_UConv]; 
  // initialize conv_comp to 0
  memset(conv_comp,0,N_UConv*SizeOfInt);

  // get information of the numbering of the dof for the space of convolution
  // first: get pointer to array where a global numbering is stored
  GlobalNumbersConv = fespaceConv->GetGlobalNumbers();
  // second: get pointer to array where the start of numbering for
  // each cell is stored
  BeginIndexConv = fespaceConv->GetBeginIndex();

  // allocate memory for geometric position of nodes
  x_conv = new double[N_DOFConv]; 
  y_conv = new double[N_DOFConv]; 
   
// ########################################################################
// computing the convolution
// loop over all cells
// ########################################################################

  // get pointer to set of mesh cells which define the fe space
  Coll = fespace->GetCollection(); 
  // get number of mesh cells
  N_Cells = Coll->GetN_Cells();

  // direct a pointer to the values of the finite element vector
  // function
  // first part of the array : u1
  // second part of the array : u2
  Values = u->GetValues();           

  // store positition of each mesh cells in the collection of mesh 
  // cells: set ClipBoard to i
  for(i=0;i<N_Cells;i++)                       
  {
    cell = Coll->GetCell(i);
    cell->SetClipBoard(i);
  }

  // loop over all mesh cells to accumulate the values for the 
  // convolution 
  for(i=0;i<N_Cells;i++)
  {
    // get current mesh cell
    cell = Coll->GetCell(i);        
    // get cell diameter (longest distance between to vertices of the cell)
    hK = cell->GetDiameter();       
    delta =  CharacteristicFilterWidth(hK);
 
    // ####################################################################
    // find local used finite elements on this cell
    // ####################################################################
    
    // initialize auxiliary array Used
    memset(Used, 0, N_FEs2D*SizeOfInt);
 
    for(j=0;j<n_fespaces;j++)
    {
      // get id of finite element space in this cell 
      CurrentElement = fespace->GetFE2D(i, cell);
      // set entry in auxiliary array to 1
      Used[CurrentElement] = 1;
    }
    
    // have a look how much 1 are in Used ->   N_LocalUsedElements
    N_LocalUsedElements = 0;
    memset(LocalUsedElements, 0, SizeOfInt*N_FEs2D);
    j = 0;
    for(k=0;k<N_FEs2D;k++)
      if(Used[k])
      {
        LocalUsedElements[j] = (FE2D)k;
        j++;
      }
    N_LocalUsedElements = j;
    
    // ####################################################################
    // calculate values on original element
    // ####################################################################
    
    // From id of fe space the situation in the referenz mesh cell is 
    // known. Now, get information in the original mesh cell. 
    // Every finite element gets automatically a quadrature rule.
    // input: N_LocalUsedElements, LocalUsedElements, cell, SecondDer
    // output : N_Points - number of quadrature points
    //          xi - xi-values of quad. points in ref. cell
    //          eta - eta-values of quad. points in ref. cell
    //          weigths - weights fro quad. rule in quad. points 
    //          X - X-values of quad. points in original cell  
    //          Y - Y-values of quad. points in original cell  
    //          AbsDetjk - absolute value of Jacobian in quad. points
    //          xi, eta, weights, X, Y,  AbsDetjk are pointers to arrays

    TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                         Coll, cell, SecondDer,
                         N_Points, xi, eta, weights, X, Y, AbsDetjk);
    
    // compute geometric positions of the fe nodes
    // get id of finite element in current mesh cell
    CurrentElement = fespace->GetFE2D(i, cell);
    CurrentElementConv = fespaceConv->GetFE2D(i, cell);
    // get fe from its id
    Element = TFEDatabase2D::GetFE2D(CurrentElementConv);
    // get reference transformation 
    RefTrans = Element->GetRefTransID();
    // get pointer to the nodal functionals (fe nodes) of the fe 
    // (in ref. cell)
    nf = Element->GetNodalFunctional2D();
    // get number and coordinates of local dof in ref cell
    // xi_ref, eta_ref are pointers
    nf->GetPointsForAll(N_loc_dofConv, xi_ref, eta_ref);
      
    // get coordinates of fe nodes in original cell
    // input: RefTrans, N_loc_dof, xi_ref, eta_ref
    // output : X_orig, Y_orig - pointers to arrays with coordinates
    //          AbsDetjk - same as above    
    TFEDatabase2D::GetOrigFromRef(RefTrans,N_loc_dofConv, xi_ref, 
                                eta_ref, X_orig, Y_orig, AbsDetjkConv);   
   
    // get pointer to basis functions of the current element
    BaseFunct = BaseFuncts[CurrentElement];  
    // # of basis functions, is the same as N_loc_dof
    N_ = N_BaseFunct[CurrentElement];        

    // find the part of the global index array where the information
    // for the current mesh cell are stored: DOF
    DOF = GlobalNumbers + BeginIndex[i];
    DOFConv = GlobalNumbersConv + BeginIndexConv[i];
    // copy the values of the finite element function from the global
    // array (Values) to local arrays 
    for(l=0;l<N_;l++)
    {
      FEFunctValues[l] = Values[DOF[l]];
      FEFunctValues1[l] = Values[DOF[l]+N_U];
    }

    // compute values of fe function in the quadrature points with
    // linear combinations of the local fe values and the basis functions
    for(k=0;k<N_Derivatives;k++)
    {
      // get fe values of the basis functions in the quad. points
      // in the current element (pointer to pointer)
      OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct,
                                            NeededDerivatives[k]);
      // accumulate the values of the fe function 
      // do it for all quad. points
      for(j=0;j<N_Points;j++)
      {
        // get pointer to array where the values for quad. point j are
        // stored
        Orig = OrigFEValues[j];
        // initialize 
        value = value1 = 0;
        // compute linear combination over all local dof
        for(l=0;l<N_;l++)
        {
          // value += u1[l] * phi[l]
          value += FEFunctValues[l] * Orig[l];
          value1 += FEFunctValues1[l] * Orig[l];
        } // endfor l
        // copy result to position in array 
        FEValue_u1_QuadPoint[j] = value;
        FEValue_u2_QuadPoint[j] = value1;
      } // endfor j
    } // endfor k
    
    // compute part of the convolution in all local fe nodes (dof)
    // which comes from current mesh cell 
    // loop over all local dof
    for(l=0;l<N_loc_dofConv;l++)         
    {
      // get global index of local dof
      index = DOFConv[l];          

      // check if value is already computed
      // if yes : continue
      if (conv_comp[index])
        continue;

      // set conv_comp[index]
      conv_comp[index] = 1;
      // initialize
      value = value1 = 0;
      // loop over all quadrature points 
      for(j=0;j<N_Points;j++) 
      {
        // square of the distance between local dof and quad. point
        // in original cell 
        distance_sq = ((X_orig[l]-X[j])*(X_orig[l]-X[j])+
                       (Y_orig[l]-Y[j])*(Y_orig[l]-Y[j])); 
        // compute value of GaussianFilter
        // input: delta, distance_sq
        g = GaussianFilter(delta,distance_sq);  
        // apply quadrature rule in ref. cell 
        value += g*AbsDetjk[j]*weights[j]*FEValue_u1_QuadPoint[j];
        value1 += g*AbsDetjk[j]*weights[j]*FEValue_u2_QuadPoint[j];
      }
      // add local result to global array
      u_conv[index] += value;
      u_conv[index+N_UConv] += value1;
      // store global coordinates of the fe nodes (dof)
      x_conv[index] = X_orig[l];
      y_conv[index] = Y_orig[l];
      x_conv[index+N_UConv] = X_orig[l];
      y_conv[index+N_UConv] = Y_orig[l];
      
      // starting to compute the integral in ALL the cells
      
      // loop over all cells
      for(ll=0;ll<N_Cells;ll++)                           
      {
        // check if we are at the current mesh cell
        if (ll==i)
          continue;
        
        neigh = Coll->GetCell(ll); 
        // get one vertex of the neigh 
        vertex0 = neigh->GetVertex(0);
        // get coordinates of vertex
        x_vertex = vertex0->GetX();
        y_vertex = vertex0->GetY();
        // compute distance to current point where the convolution
        // is computed
        switch(TDatabase::ParamDB->INTERNAL_PERIODIC_IDENTITY)
        {
          case 0:
            distance_sq = ((X_orig[l]-x_vertex)*(X_orig[l]-x_vertex)+
                           (Y_orig[l]-y_vertex)*(Y_orig[l]-y_vertex));
            break;
          case 1:
            exit(4711);
            break;
          case 2:
            // periodic square (1,1)^2 
            if (fabs(X_orig[l]-x_vertex)> 1.5)
            {
              if (X_orig[l]-x_vertex>1.5)
                distance_sq = ((X_orig[l]-x_vertex-2)*(X_orig[l]-x_vertex-2)+
                               (Y_orig[l]-y_vertex)*(Y_orig[l]-y_vertex)); 
              else
                distance_sq = ((X_orig[l]-x_vertex+2)*(X_orig[l]-x_vertex+2)+
                               (Y_orig[l]-y_vertex)*(Y_orig[l]-y_vertex));                                       
            }
            else
              distance_sq = ((X_orig[l]-x_vertex)*(X_orig[l]-x_vertex)+
                             (Y_orig[l]-y_vertex)*(Y_orig[l]-y_vertex));                   
            break;
        }
        //        distance_sq = ((X_orig[l]-x_vertex)*(X_orig[l]-x_vertex)+
        //               (Y_orig[l]-y_vertex)*(Y_orig[l]-y_vertex));
        // check if neigh is close enough
        if (distance_sq > 2*delta *delta)
          continue;
        
        // get id of neighbour
        neigh_i = neigh->GetClipBoard();
        
        // get id of finite element in neigh mesh cell wrt the finite
        // element space of the velocity
        CurrentElementNeigh = fespace->GetFE2D(neigh_i, neigh);
        
        // get pointer to basis functions of the neighbour element
        BaseFunctNeigh = BaseFuncts[CurrentElementNeigh]; 
 
        // # of basis functions, is the same as N_loc_dof
        N_Neigh = N_BaseFunct[CurrentElementNeigh]; 
     
        // find the part of the global index array where the information
        // for the neigh mesh cell are stored: DOF
        DOFNeigh = GlobalNumbers + BeginIndex[neigh_i];

        // initialize auxiliary array UsedNeigh
        memset(UsedNeigh, 0, N_FEs2D*SizeOfInt);
        
        // set entry in auxiliary array to 1
        UsedNeigh[CurrentElementNeigh] = 1;
       
        // have a look how much 1 are in UsedNeigh -> N_LocalUsedElementsNeigh
        N_LocalUsedElementsNeigh = 0;
        memset(LocalUsedElementsNeigh, 0, SizeOfInt*N_FEs2D);
        j = 0;
        for(k=0;k<N_FEs2D;k++)
          if(UsedNeigh[k])
          {
            LocalUsedElementsNeigh[j] = (FE2D)k;
            j++;
          }
        N_LocalUsedElementsNeigh = j;
        
        // ####################################################################
        // calculate values on the original element
        // ####################################################################
        
        // From id of fe space the situation in the referenz mesh cell is 
        // known. Now, get information in the original mesh cell. 
        // Every finite element gets automatically a quadrature rule.
        // input: N_LocalUsedElementsNeigh, LocalUsedElementsNeigh, 
        //        neigh, SecondDer
        // output : N_PointsNeigh - number of quadrature points
        //          xiNeigh - xi-values of quad. points in ref. cell
        //          etaNeigh - eta-values of quad. points in ref. cell
        //          weigthsNeigh - weights fro quad. rule in quad. points 
        //          XNeigh - X-values of quad. points in original cell  
        //          YNeigh - Y-values of quad. points in original cell  
        //          AbsDetjkNeigh - absolute value of Jacobian in quad. points
        //          xiNeigh, etaNeigh, weightsNeigh, XNeigh, YNeigh,  
        //          AbsDetjkNeigh are pointers to arrays

        TFEDatabase2D::GetOrig(N_LocalUsedElementsNeigh, LocalUsedElementsNeigh, 
                             Coll, neigh, SecondDer,
                             N_PointsNeigh, xiNeigh, etaNeigh, weightsNeigh,
                             XNeigh, YNeigh, AbsDetjkNeigh);
    
        // copy the values of the finite element function from the global
        // array (Values) to local arrays 
        for(lll=0;lll<N_Neigh;lll++)
        {
          FEFunctValuesNeigh[lll] = Values[DOFNeigh[lll]];
          FEFunctValues1Neigh[lll] = Values[DOFNeigh[lll]+N_U];
        }
        
        // compute values of fe function in the quadrature points with
        // linear combinations of the local fe values and the basis functions
        for(k=0;k<N_Derivatives;k++)
        {
          // get fe values of the basis functions in the quad. points
          // in the current element (pointer to pointer)
          OrigFEValuesNeigh = TFEDatabase2D::GetOrigElementValues(BaseFunctNeigh,
                                                                NeededDerivatives[k]);
          // accumulate the values of the fe function 
          // do it for all quad. points
          for(j=0;j<N_PointsNeigh;j++)
          {
            // get pointer to array where the values for quad. point j are
            // stored
            Orig = OrigFEValuesNeigh[j];
            // initialize 
            value = value1 = 0;
            // compute linear combination over all local dof
            for(lll=0;lll<N_Neigh;lll++)
            {
              // value += u1[l] * phi[l]
              value += FEFunctValuesNeigh[lll] * Orig[lll];
              value1 += FEFunctValues1Neigh[lll] * Orig[lll];
            } // endfor lll
            // copy result to position in array 
            FEValue_u1_QuadPointNeigh[j] = value;
            FEValue_u2_QuadPointNeigh[j] = value1;
          } // endfor j
        } // endfor k
        
        // compute part of the convolution in all local fe nodes (dof)
        // which comes from current mesh cell 
        // loop over all local dof
        
        // initialize
        value = value1 = 0;
        // loop over all quadrature points 
        for(j=0;j<N_PointsNeigh;j++) 
        { 
          // square of the distance between local dof and quad. point
          // in original cell 
          switch(TDatabase::ParamDB->INTERNAL_PERIODIC_IDENTITY)
          {
            case 0:
              distance_sq = ((X_orig[l]-XNeigh[j])*(X_orig[l]-XNeigh[j])+
                             (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j])); 
              break;
            case 1:
              exit (4711);
              break;
            case 2:
              // periodic square (1,1)^2 
                if (fabs(X_orig[l]-XNeigh[j])> 1.5)
                {
                  if (X_orig[l]-XNeigh[j]>1.5)
                    distance_sq = ((X_orig[l]-XNeigh[j]-2)*(X_orig[l]-XNeigh[j]-2)+
                                   (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j])); 
                  else
                    distance_sq = ((X_orig[l]-XNeigh[j]+2)*(X_orig[l]-XNeigh[j]+2)+
                                   (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j]));                                       
                }
                else
                  distance_sq = ((X_orig[l]-XNeigh[j])*(X_orig[l]-XNeigh[j])+
                                 (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j]));                   
              break;
          }
          //distance_sq = ((X_orig[l]-XNeigh[j])*(X_orig[l]-XNeigh[j])+
          //               (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j])); 
          //cout << "X " << XNeigh[j] << " Y " << YNeigh[j] << endl;
          // compute value of GaussianFilter
          // input: delta, distance_sq
          g = GaussianFilter(delta,distance_sq);
          //cout << "dist " <<  sqrt(distance_sq) << endl;
          
          // apply quadrature rule in ref. cell 
          value += g*AbsDetjkNeigh[j]*weightsNeigh[j]*FEValue_u1_QuadPointNeigh[j];
          value1 += g*AbsDetjkNeigh[j]*weightsNeigh[j]*FEValue_u2_QuadPointNeigh[j];
       
        }
        // add local `result to global array
        u_conv[index] += value;
        u_conv[index+N_UConv] += value1; 
      } // end of ll loop (edges)
    } // end of l loop (local dof)
  } // endfor i (cells)
  
 
  // for (i=0;i<N_DOFConv; i++)
  //{
  //  cout << i<< " x " <<   x_conv[i] << " y  " << y_conv[i] ;
  //  cout << " u_conv " << u_conv[i] << endl;
  //}

  // release memory
  delete SecondDer;
  delete conv_comp;
  delete x_conv;
  delete y_conv;

} // ConvoluteVelocityFull

// ========================================================================
// convolute (grad w grad w^T)
// ========================================================================
void  ConvoluteDuTensor(TFEVectFunct2D *u, TFEVectFunct2D *duTensor)
{
  int i,j,k,l,n,m, N_UsedElements, N_LocalUsedElements;
  int N_LocalUsedElementsNeigh;
  int N_Cells, N_Points, N_Parameters, N_, N_U, N_DOF, N_loc_dofConv;
  int Used[N_FEs2D], UsedNeigh[N_FEs2D], *N_BaseFunct;
  int neigh_i, N_Neigh, ll, lll, N_PointsNeigh, N_Edges, same_dof;
  FE2D LocalUsedElements[N_FEs2D], CurrentElement, CurrentElementConv;
  FE2D CurrentElementNeigh,LocalUsedElementsNeigh[N_FEs2D] ;
  TFE2D *Element;
  TNodalFunctional2D *nf;
  BaseFunct2D BaseFunct, *BaseFuncts, BaseFunctNeigh;
  TCollection *Coll;
  TBaseCell *cell, *neigh;
  TFE2D *ele;
  double *weights, *xi, *eta, *xi_ref, *eta_ref;
  double *weightsNeigh, *xiNeigh, *etaNeigh;
  double X_orig[MaxN_PointsForNodal2D], Y_orig[MaxN_PointsForNodal2D];
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D], AbsDetjkConv[MaxN_QuadPoints_2D];
  double XNeigh[MaxN_QuadPoints_2D], YNeigh[MaxN_QuadPoints_2D];
  double AbsDetjkNeigh[MaxN_QuadPoints_2D];
  RefTrans2D RefTrans;
  double *Param[MaxN_QuadPoints_2D], *aux;
  double FEValue_D1u1_QuadPoint[MaxN_QuadPoints_2D];
  double FEValue_D1u2_QuadPoint[MaxN_QuadPoints_2D];
  double FEValue_D2u1_QuadPoint[MaxN_QuadPoints_2D];
  double FEValue_D2u2_QuadPoint[MaxN_QuadPoints_2D];
  double FEValue_D1u1_QuadPointNeigh[MaxN_QuadPoints_2D];
  double FEValue_D1u2_QuadPointNeigh[MaxN_QuadPoints_2D];
  double FEValue_D2u1_QuadPointNeigh[MaxN_QuadPoints_2D];
  double FEValue_D2u2_QuadPointNeigh[MaxN_QuadPoints_2D];
  double *ExactVal[MaxN_QuadPoints_2D];
  int *DOF, ActiveBound, DirichletBound, end, last, *DOFConv, *DOFNeigh;
  double **OrigFEValues, *Orig;
  double value, value1;
  double value11, value12, value21, value22;
  double **OrigFEValuesNeigh;
  double FEFunctValues[MaxN_BaseFunctions2D];
  double FEFunctValues1[MaxN_BaseFunctions2D];
  double FEFunctValuesNeigh[MaxN_BaseFunctions2D];
  double FEFunctValues1Neigh[MaxN_BaseFunctions2D];
  int *GlobalNumbers, *BeginIndex, *GlobalNumbersConv, *BeginIndexConv;
  double LocError[4];
  int N_UConv, N_DOFConv, index;
  double *Values, *u_conv,distance_sq,h_fine,g, *u_values;
  double *x_conv, *y_conv, delta, hK;
  TJoint *joint;
  TFESpace2D *fespace, *fespaceConv; 

  bool SecondDer[1] = { FALSE };
  int N_Derivatives = 2;
  MultiIndex2D NeededDerivatives[2] = { D10, D01 };

  // get fe spaces of velocity u and of covoluted tensor
  fespace = u->GetFESpace2D();
  fespaceConv = duTensor->GetFESpace2D();

  // gives a pointer to all available basis functions which
  // are described in the fedatabase
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  // gives the number of available basis functions from the
  // fedatabase
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  // get information of the numbering of the degrees of freedom
  // first: get pointer to array where a global numbering is stored
  GlobalNumbers = fespace->GetGlobalNumbers();
  // second: get pointer to array where the start of numbering for
  // each cell is stored
  BeginIndex = fespace->GetBeginIndex();
 
  // get number of dof for one velocity component N_U
  // -> total number of dof is twice N_U
  N_U = u->GetLength();
  N_DOF = 2*N_U;
  
  // get pointer to set of mesh cells which define the fe space
  Coll = fespace->GetCollection();
  // get number of mesh cells
  N_Cells = Coll->GetN_Cells();
 
  // # dof in the new fe space
  N_UConv = fespaceConv->GetN_DegreesOfFreedom();
 
  N_DOFConv = 4*N_UConv;
  // allocate memory for values of convoluted function
  u_conv = duTensor->GetValues(); 
  // initialize u_conv to 0
  memset(u_conv,0,N_DOFConv*SizeOfDouble);

  // get information of the numbering of the dof for the space of convolution
  // first: get pointer to array where a global numbering is stored
  GlobalNumbersConv = fespaceConv->GetGlobalNumbers();
  // second: get pointer to array where the start of numbering for
  // each cell is stored
  BeginIndexConv = fespaceConv->GetBeginIndex();

  // allocate memory for geometric position of nodes
  x_conv = new double[N_DOFConv]; 
  y_conv = new double[N_DOFConv]; 
   
// ########################################################################
// computing the convolution
// loop over all cells
// ########################################################################

  // direct a pointer to the values of the finite element vector
  // function
  // first part of the array : u1
  // second part of the array : u2
  Values = u->GetValues();           
  OutPut("U " << sqrt(Ddot(N_DOF,
                                 Values,
                                 Values)) << endl);


  // store positition of each mesh cells in the collection of mesh 
  // cells: set ClipBoard to i
  for(i=0;i<N_Cells;i++)                       
  {
    cell = Coll->GetCell(i);
    cell->SetClipBoard(i);
  }

  // loop over all mesh cells to accumulate the values for the 
  // convolution 
  for(i=0;i<N_Cells;i++)
  {
    // get current mesh cell
    cell = Coll->GetCell(i);        
    hK = cell->GetDiameter();
    delta =  CharacteristicFilterWidth(hK);

    // ####################################################################
    // find local used finite elements on this cell
    // ####################################################################
    CurrentElement = fespace->GetFE2D(i, cell);
    LocalUsedElements[0] = CurrentElement;
    N_LocalUsedElements = 1;
    
    // ####################################################################
    // calculate values on original element
    // ####################################################################
    
    // From id of fe space the situation in the referenz mesh cell is 
    // known. Now, get information in the original mesh cell. 
    // Every finite element gets automatically a quadrature rule.
    // input: N_LocalUsedElements, LocalUsedElements, cell, SecondDer
    // output : N_Points - number of quadrature points
    //          xi - xi-values of quad. points in ref. cell
    //          eta - eta-values of quad. points in ref. cell
    //          weigths - weights fro quad. rule in quad. points 
    //          X - X-values of quad. points in original cell  
    //          Y - Y-values of quad. points in original cell  
    //          AbsDetjk - absolute value of Jacobian in quad. points
    //          xi, eta, weights, X, Y,  AbsDetjk are pointers to arrays

    TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                         Coll, cell, SecondDer,
                         N_Points, xi, eta, weights, X, Y, AbsDetjk);
 
    // compute geometric positions of the fe nodes
    // get id of finite element in current mesh cell
    CurrentElement = fespace->GetFE2D(i, cell);
    CurrentElementConv = fespaceConv->GetFE2D(i, cell);
    // get fe from its id
    Element = TFEDatabase2D::GetFE2D(CurrentElementConv);
    // get reference transformation 
    RefTrans = Element->GetRefTransID();
    // get pointer to the nodal functionals (fe nodes) of the fe 
    // (in ref. cell)
    nf = Element->GetNodalFunctional2D();
    // get number and coordinates of local dof in ref cell
    // xi_ref, eta_ref are pointers
    nf->GetPointsForAll(N_loc_dofConv, xi_ref, eta_ref);
      
    // get coordinates of fe nodes in original cell
    // input: RefTrans, N_loc_dof, xi_ref, eta_ref
    // output : X_orig, Y_orig - pointers to arrays with coordinates
    //          AbsDetjk - same as above    
    TFEDatabase2D::GetOrigFromRef(RefTrans,N_loc_dofConv, xi_ref, 
                                eta_ref, X_orig, Y_orig, AbsDetjkConv);   
   
    // get pointer to basis functions of the current element
    BaseFunct = BaseFuncts[CurrentElement];  
    // # of basis functions, is the same as N_loc_dof
    N_ = N_BaseFunct[CurrentElement];        

    // find the part of the global index array where the information
    // for the current mesh cell are stored: DOF
    DOF = GlobalNumbers + BeginIndex[i];
    DOFConv = GlobalNumbersConv + BeginIndexConv[i];
    // copy the values of the finite element function from the global
    // array (Values) to local arrays 
    for(l=0;l<N_;l++)
    {
      FEFunctValues[l] = Values[DOF[l]];
      FEFunctValues1[l] = Values[DOF[l]+N_U];
    }

    // compute values of fe function in the quadrature points with
    // linear combinations of the local fe values and the basis functions
    for(k=0;k<N_Derivatives;k++)
    {
      // get fe values of the basis functions in the quad. points
      // in the current element (pointer to pointer)
      OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct,
                                            NeededDerivatives[k]);
      // accumulate the values of the fe function 
      // do it for all quad. points
      for(j=0;j<N_Points;j++)
      {
        // get pointer to array where the values for quad. point j are
        // stored
        Orig = OrigFEValues[j];
        // initialize 
        value = value1 = 0;
        // compute linear combination over all local dof
        for(l=0;l<N_;l++)
        {
          // value += u1[l] * phi[l]
          value += FEFunctValues[l] * Orig[l];
          value1 += FEFunctValues1[l] * Orig[l];
        } // endfor l
        // copy result to position in array 

        if(k==0)
        {
          FEValue_D1u1_QuadPoint[j] = value;
          FEValue_D1u2_QuadPoint[j] = value1;
        }
        else
        {
          FEValue_D2u1_QuadPoint[j] = value;
          FEValue_D2u2_QuadPoint[j] = value1;
        }
      } // endfor j
    } // endfor k
    
    // compute part of the convolution in all local fe nodes (dof)
    // which comes from current mesh cell 
    // loop over all local dof
    for(l=0;l<N_loc_dofConv;l++)         
    {
      // get global indes of local dof
      index = DOFConv[l];          
      // initialize
      value11 = value12 = value21 = value22 = 0;
      // loop over all quadrature points 
      for(j=0;j<N_Points;j++) 
      {
        // square of the distance between local dof and quad. point
        // in original cell 
        distance_sq = ((X_orig[l]-X[j])*(X_orig[l]-X[j])+
                       (Y_orig[l]-Y[j])*(Y_orig[l]-Y[j])); 
        // compute value of GaussianFilter
        // input: delta, distance_sq
        g = GaussianFilter(delta,distance_sq);  
        // apply quadrature rule in ref. cell 
        value11 += g*AbsDetjk[j]*weights[j]*
                   ( FEValue_D1u1_QuadPoint[j]*FEValue_D1u1_QuadPoint[j]
                    +FEValue_D2u1_QuadPoint[j]*FEValue_D2u1_QuadPoint[j]);
        value12 += g*AbsDetjk[j]*weights[j]*
                   ( FEValue_D1u1_QuadPoint[j]*FEValue_D1u2_QuadPoint[j]
                    +FEValue_D2u1_QuadPoint[j]*FEValue_D2u2_QuadPoint[j]);
        value21 += g*AbsDetjk[j]*weights[j]*
                   ( FEValue_D1u2_QuadPoint[j]*FEValue_D1u1_QuadPoint[j]
                    +FEValue_D2u2_QuadPoint[j]*FEValue_D2u1_QuadPoint[j]);
        value22 += g*AbsDetjk[j]*weights[j]*
                   ( FEValue_D1u2_QuadPoint[j]*FEValue_D1u2_QuadPoint[j]
                    +FEValue_D2u2_QuadPoint[j]*FEValue_D2u2_QuadPoint[j]);
      }
      // add local result to global array
      u_conv[index] += value11;
      u_conv[index+N_UConv] += value12;
      u_conv[index+2*N_UConv] += value21;
      u_conv[index+3*N_UConv] += value22;
      // store global coordinates of the fe nodes (dof)
      x_conv[index] = X_orig[l];
      y_conv[index] = Y_orig[l];
      x_conv[index+N_UConv] = X_orig[l];
      y_conv[index+N_UConv] = Y_orig[l];
      
      // starting to compute the integral in the neighbour cells
      // number of edges
      N_Edges=cell->GetN_Edges();
      // loop over all edges
      for(ll=0;ll<N_Edges;ll++)                           
      {
        // get pointer to edge[l]
        joint=cell->GetJoint(ll); 
        
        // if boundary edge continue 
        if ((joint->GetType() == BoundaryEdge)||
            (joint->GetType() == IsoBoundEdge))
          continue;
        
        // get pointer to neighbour cell
        neigh=cell->GetJoint(ll)->GetNeighbour(cell);
        
        // get id of neighbour
        neigh_i = neigh->GetClipBoard();
        
        // get id of finite element in neigh mesh cell wrt the
        // finite element space of the convolution
        CurrentElementNeigh = fespaceConv->GetFE2D(neigh_i, neigh);
        
        // get pointer to basis functions of the neighbour element
        BaseFunctNeigh = BaseFuncts[CurrentElementNeigh]; 
 
        // # of basis functions, is the same as N_loc_dof
        N_Neigh = N_BaseFunct[CurrentElementNeigh]; 
     
        // find the part of the global index array where the information
        // for the neigh mesh cell are stored: DOF
        DOFNeigh = GlobalNumbersConv + BeginIndexConv[neigh_i];
        
        // check if computation on neigh cell necessary
        // for all dof
        same_dof = 0;
        for (j=0;j<N_Neigh;j++)
        {
          //cout << "N " << DOFNeigh[j] << " ind " << index << endl;
          if (DOFNeigh[j]==index)
          same_dof = 1;
        }
        // if the cells have the same dof, continue
        if (same_dof==1)
          continue;
        
        ////////////////////////////////////////////////////////////////////
        
        // if not, do something! 
        // i.e. the same thing as before
        
        // get id of finite element in neigh mesh cell wrt the finite
        // element space of the velocity
        CurrentElementNeigh = fespace->GetFE2D(neigh_i, neigh);
        
        // get pointer to basis functions of the neighbour element
        BaseFunctNeigh = BaseFuncts[CurrentElementNeigh]; 
 
        // # of basis functions, is the same as N_loc_dof
        N_Neigh = N_BaseFunct[CurrentElementNeigh]; 
     
        // find the part of the global index array where the information
        // for the neigh mesh cell are stored: DOF
        DOFNeigh = GlobalNumbers + BeginIndex[neigh_i];

        // initialize auxiliary array UsedNeigh
        memset(UsedNeigh, 0, N_FEs2D*SizeOfInt);
        
        // get id of finite element space in neigh cell 
        CurrentElementNeigh = fespace->GetFE2D(neigh_i, neigh);
        // set entry in auxiliary array to 1
        UsedNeigh[CurrentElementNeigh] = 1;
       
        // have a look how much 1 are in UsedNeigh -> N_LocalUsedElementsNeigh
        N_LocalUsedElementsNeigh = 1;
        LocalUsedElementsNeigh[0] = CurrentElementNeigh;
        
        // ####################################################################
        // calculate values on original element
        // ####################################################################
        
        // From id of fe space the situation in the referenz mesh cell is 
        // known. Now, get information in the original mesh cell. 
        // Every finite element gets automatically a quadrature rule.
        // input: N_LocalUsedElementsNeigh, LocalUsedElementsNeigh, 
        //        neigh, SecondDer
        // output : N_PointsNeigh - number of quadrature points
        //          xiNeigh - xi-values of quad. points in ref. cell
        //          etaNeigh - eta-values of quad. points in ref. cell
        //          weigthsNeigh - weights fro quad. rule in quad. points 
        //          XNeigh - X-values of quad. points in original cell  
        //          YNeigh - Y-values of quad. points in original cell  
        //          AbsDetjkNeigh - absolute value of Jacobian in quad. points
        //          xiNeigh, etaNeigh, weightsNeigh, XNeigh, YNeigh,  
        //          AbsDetjkNeigh are pointers to arrays

        TFEDatabase2D::GetOrig(N_LocalUsedElementsNeigh, LocalUsedElementsNeigh, 
                             Coll, neigh, SecondDer,
                             N_PointsNeigh, xiNeigh, etaNeigh, weightsNeigh,
                             XNeigh, YNeigh, AbsDetjkNeigh);
    
        // copy the values of the finite element function from the global
        // array (Values) to local arrays 
        for(lll=0;lll<N_Neigh;lll++)
        {
          FEFunctValuesNeigh[lll] = Values[DOFNeigh[lll]];
          FEFunctValues1Neigh[lll] = Values[DOFNeigh[lll]+N_U];
        }
        
        // compute values of fe function in the quadrature points with
        // linear combinations of the local fe values and the basis functions
        for(k=0;k<N_Derivatives;k++)
        {
          // get fe values of the basis functions in the quad. points
          // in the current element (pointer to pointer)
          OrigFEValuesNeigh = TFEDatabase2D::GetOrigElementValues(BaseFunctNeigh,
                                                       NeededDerivatives[k]);
          // accumulate the values of the fe function 
          // do it for all quad. points
          for(j=0;j<N_PointsNeigh;j++)
          {
            // get pointer to array where the values for quad. point j are
            // stored
            Orig = OrigFEValuesNeigh[j];
            // initialize 
            value = value1 = 0;
            // compute linear combination over all local dof
            for(lll=0;lll<N_Neigh;lll++)
            {
              // value += u1[l] * phi[l]
              value += FEFunctValuesNeigh[lll] * Orig[lll];
              value1 += FEFunctValues1Neigh[lll] * Orig[lll];
            } // endfor lll
            // copy result to position in array 

            if(k==0)
            {
              FEValue_D1u1_QuadPointNeigh[j] = value;
              FEValue_D1u2_QuadPointNeigh[j] = value1;
            }
            else
            {
              FEValue_D2u1_QuadPointNeigh[j] = value;
              FEValue_D2u2_QuadPointNeigh[j] = value1;
            }
          } // endfor j
        } // endfor k
        
        // compute part of the convolution in all local fe nodes (dof)
        // which comes from current mesh cell 
        // loop over all local dof
        
        // initialize
        value11 = value12 = value21 = value22 = 0;
        // loop over all quadrature points 
        for(j=0;j<N_PointsNeigh;j++) 
        { 
          // square of the distance between local dof and quad. point
          // in original cell 
          distance_sq = ((X_orig[l]-XNeigh[j])*(X_orig[l]-XNeigh[j])+
                         (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j])); 
          //cout << "X " << XNeigh[j] << " Y " << YNeigh[j] << endl;
          // compute value of GaussianFilter
          // input: delta, distance_sq
          g = GaussianFilter(delta,distance_sq);
          //cout << "dist " <<  sqrt(distance_sq) << endl;
          //
          // apply quadrature rule in ref. cell 
          value11 += g*AbsDetjk[j]*weightsNeigh[j]*
                     ( FEValue_D1u1_QuadPointNeigh[j]*
                       FEValue_D1u1_QuadPointNeigh[j]
                      +FEValue_D2u1_QuadPointNeigh[j]*
                       FEValue_D2u1_QuadPointNeigh[j]);
          value12 += g*AbsDetjk[j]*weightsNeigh[j]*
                     ( FEValue_D1u1_QuadPointNeigh[j]*
                       FEValue_D1u2_QuadPointNeigh[j]
                      +FEValue_D2u1_QuadPointNeigh[j]*
                       FEValue_D2u2_QuadPointNeigh[j]);
          value21 += g*AbsDetjk[j]*weightsNeigh[j]*
                     ( FEValue_D1u2_QuadPointNeigh[j]*
                       FEValue_D1u1_QuadPointNeigh[j]
                      +FEValue_D2u2_QuadPointNeigh[j]*
                       FEValue_D2u1_QuadPointNeigh[j]);
          value22 += g*AbsDetjk[j]*weightsNeigh[j]*
                     ( FEValue_D1u2_QuadPointNeigh[j]*
                       FEValue_D1u2_QuadPointNeigh[j]
                      +FEValue_D2u2_QuadPointNeigh[j]*
                       FEValue_D2u2_QuadPointNeigh[j]);
        }
        // add local result to global array
        u_conv[index] += value11;
        u_conv[index+N_UConv] += value12; 
        u_conv[index+2*N_UConv] += value21; 
        u_conv[index+3*N_UConv] += value22; 
      } // end of ll loop (edges)
    } // end of l loop (local dof)
  } // endfor i (cells)
  
 
  OutPut("C0 " << sqrt(Ddot(N_UConv,
                                 u_conv,
                                 u_conv)) << endl);
  OutPut("C1 " << sqrt(Ddot(N_UConv,
                                 u_conv+N_UConv,
                                 u_conv+N_UConv)) << endl);
  OutPut("C2 " << sqrt(Ddot(N_UConv,
                                 u_conv+2*N_UConv,
                                 u_conv+2*N_UConv)) << endl);
  OutPut("C3 " << sqrt(Ddot(N_UConv,
                                 u_conv+3*N_UConv,
                                 u_conv+3*N_UConv)) << endl);

  // for (i=0;i<N_DOFConv; i++)
  // {
  //   cout << u_conv[i] << "  " << u_conv[i+N_UConv] << "  ";
  //   cout << u_conv[i+2*N_UConv] << "  " << u_conv[i+3*N_UConv] << endl;
  // }

  // release memory
  delete x_conv;
  delete y_conv;

} // ConvoluteDuTensor

// ========================================================================
// convolute D = (grad w grad w^T) by local integration
// use D_{21} = D_{12}
// input u
// output DuTensor
// ========================================================================
void  ConvoluteSymmetricTensor(TFEVectFunct2D *u, TFEVectFunct2D *duTensor)

{
  int i,j,k,l,n,m, N_UsedElements, N_LocalUsedElements;
  int N_LocalUsedElementsNeigh;
  int N_Cells, N_Points, N_Parameters, N_, N_U, N_DOF, N_loc_dofConv;
  int Used[N_FEs2D], UsedNeigh[N_FEs2D], *N_BaseFunct;
  int neigh_i, N_Neigh, ll, lll, N_PointsNeigh, N_Edges, same_dof;
  FE2D LocalUsedElements[N_FEs2D], CurrentElement, CurrentElementConv;
  FE2D CurrentElementNeigh,LocalUsedElementsNeigh[N_FEs2D] ;
  TFE2D *Element;
  TNodalFunctional2D *nf;
  BaseFunct2D BaseFunct, *BaseFuncts, BaseFunctNeigh;
  TCollection *Coll;
  TBaseCell *cell, *neigh;
  TFE2D *ele;
  double *weights, *xi, *eta, *xi_ref, *eta_ref;
  double *weightsNeigh, *xiNeigh, *etaNeigh;
  double X_orig[MaxN_PointsForNodal2D], Y_orig[MaxN_PointsForNodal2D];
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D], AbsDetjkConv[MaxN_QuadPoints_2D];
  double XNeigh[MaxN_QuadPoints_2D], YNeigh[MaxN_QuadPoints_2D];
  double AbsDetjkNeigh[MaxN_QuadPoints_2D];
  RefTrans2D RefTrans;
  double *Param[MaxN_QuadPoints_2D], *aux;
  double FEValue_D1u1_QuadPoint[MaxN_QuadPoints_2D];
  double FEValue_D1u2_QuadPoint[MaxN_QuadPoints_2D];
  double FEValue_D2u1_QuadPoint[MaxN_QuadPoints_2D];
  double FEValue_D2u2_QuadPoint[MaxN_QuadPoints_2D];
  double FEValue_D1u1_QuadPointNeigh[MaxN_QuadPoints_2D];
  double FEValue_D1u2_QuadPointNeigh[MaxN_QuadPoints_2D];
  double FEValue_D2u1_QuadPointNeigh[MaxN_QuadPoints_2D];
  double FEValue_D2u2_QuadPointNeigh[MaxN_QuadPoints_2D];
  double *ExactVal[MaxN_QuadPoints_2D];
  int *DOF, ActiveBound, DirichletBound, end, last, *DOFConv, *DOFNeigh;
  double **OrigFEValues, *Orig;
  double value, value1;
  double value11, value12, value22;
  double **OrigFEValuesNeigh;
  double FEFunctValues[MaxN_BaseFunctions2D];
  double FEFunctValues1[MaxN_BaseFunctions2D];
  double FEFunctValuesNeigh[MaxN_BaseFunctions2D];
  double FEFunctValues1Neigh[MaxN_BaseFunctions2D];
  int *GlobalNumbers, *BeginIndex, *GlobalNumbersConv, *BeginIndexConv;
  double LocError[4];
  int N_UConv, N_DOFConv, index;
  double *Values, *u_conv,distance_sq,h_fine,g, *u_values;
  double *x_conv, *y_conv, delta, hK;
  TJoint *joint;
  TFESpace2D *fespace, *fespaceConv; 

  bool SecondDer[1] = { FALSE };
  int N_Derivatives = 2;
  MultiIndex2D NeededDerivatives[2] = { D10, D01 };

  // get fe spaces of velocity u and of covoluted tensor
  fespace = u->GetFESpace2D();
  fespaceConv = duTensor->GetFESpace2D();

  // gives a pointer to all available basis functions which
  // are described in the fedatabase
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  // gives the number of available basis functions from the
  // fedatabase
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  // get information of the numbering of the degrees of freedom
  // first: get pointer to array where a global numbering is stored
  GlobalNumbers = fespace->GetGlobalNumbers();
  // second: get pointer to array where the start of numbering for
  // each cell is stored
  BeginIndex = fespace->GetBeginIndex();
 
  // get number of dof for one velocity component N_U
  // -> total number of dof is twice N_U
  N_U = u->GetLength();
  N_DOF = 2*N_U;
  
  // get pointer to set of mesh cells which define the fe space
  Coll = fespace->GetCollection();
  // get number of mesh cells
  N_Cells = Coll->GetN_Cells();
 
  // # dof in the new fe space
  N_UConv = fespaceConv->GetN_DegreesOfFreedom();
 
  N_DOFConv = 3*N_UConv;
  // allocate memory for values of convoluted function
  u_conv = duTensor->GetValues(); 
  // initialize u_conv to 0
  memset(u_conv,0,N_DOFConv*SizeOfDouble);

  // get information of the numbering of the dof for the space of convolution
  // first: get pointer to array where a global numbering is stored
  GlobalNumbersConv = fespaceConv->GetGlobalNumbers();
  // second: get pointer to array where the start of numbering for
  // each cell is stored
  BeginIndexConv = fespaceConv->GetBeginIndex();

  // allocate memory for geometric position of nodes
  x_conv = new double[N_DOFConv]; 
  y_conv = new double[N_DOFConv]; 
   
// ########################################################################
// computing the convolution
// loop over all cells
// ########################################################################

  // direct a pointer to the values of the finite element vector
  // function
  // first part of the array : u1
  // second part of the array : u2
  Values = u->GetValues();           

  // store positition of each mesh cells in the collection of mesh 
  // cells: set ClipBoard to i
  for(i=0;i<N_Cells;i++)                       
  {
    cell = Coll->GetCell(i);
    cell->SetClipBoard(i);
  }

  // loop over all mesh cells to accumulate the values for the 
  // convolution 
  for(i=0;i<N_Cells;i++)
  {
    // get current mesh cell
    cell = Coll->GetCell(i);        
    hK = cell->GetDiameter();
    delta =  CharacteristicFilterWidth(hK);

    // ####################################################################
    // find local used finite elements on this cell
    // ####################################################################
    CurrentElement = fespace->GetFE2D(i, cell);
    LocalUsedElements[0] = CurrentElement;
    N_LocalUsedElements = 1;
    
    // ####################################################################
    // calculate values on original element
    // ####################################################################
    
    // From id of fe space the situation in the referenz mesh cell is 
    // known. Now, get information in the original mesh cell. 
    // Every finite element gets automatically a quadrature rule.
    // input: N_LocalUsedElements, LocalUsedElements, cell, SecondDer
    // output : N_Points - number of quadrature points
    //          xi - xi-values of quad. points in ref. cell
    //          eta - eta-values of quad. points in ref. cell
    //          weigths - weights fro quad. rule in quad. points 
    //          X - X-values of quad. points in original cell  
    //          Y - Y-values of quad. points in original cell  
    //          AbsDetjk - absolute value of Jacobian in quad. points
    //          xi, eta, weights, X, Y,  AbsDetjk are pointers to arrays

    TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                         Coll, cell, SecondDer,
                         N_Points, xi, eta, weights, X, Y, AbsDetjk);
 
    // compute geometric positions of the fe nodes
    // get id of finite element in current mesh cell
    CurrentElement = fespace->GetFE2D(i, cell);
    CurrentElementConv = fespaceConv->GetFE2D(i, cell);
    // get fe from its id
    Element = TFEDatabase2D::GetFE2D(CurrentElementConv);
    // get reference transformation 
    RefTrans = Element->GetRefTransID();
    // get pointer to the nodal functionals (fe nodes) of the fe 
    // (in ref. cell)
    nf = Element->GetNodalFunctional2D();
    // get number and coordinates of local dof in ref cell
    // xi_ref, eta_ref are pointers
    nf->GetPointsForAll(N_loc_dofConv, xi_ref, eta_ref);
      
    // get coordinates of fe nodes in original cell
    // input: RefTrans, N_loc_dof, xi_ref, eta_ref
    // output : X_orig, Y_orig - pointers to arrays with coordinates
    //          AbsDetjk - same as above    
    TFEDatabase2D::GetOrigFromRef(RefTrans,N_loc_dofConv, xi_ref, 
                                eta_ref, X_orig, Y_orig, AbsDetjkConv);   
   
    // get pointer to basis functions of the current element
    BaseFunct = BaseFuncts[CurrentElement];  
    // # of basis functions, is the same as N_loc_dof
    N_ = N_BaseFunct[CurrentElement];        

    // find the part of the global index array where the information
    // for the current mesh cell are stored: DOF
    DOF = GlobalNumbers + BeginIndex[i];
    DOFConv = GlobalNumbersConv + BeginIndexConv[i];
    // copy the values of the finite element function from the global
    // array (Values) to local arrays 
    for(l=0;l<N_;l++)
    {
      FEFunctValues[l] = Values[DOF[l]];
      FEFunctValues1[l] = Values[DOF[l]+N_U];
    }

    // compute values of fe function in the quadrature points with
    // linear combinations of the local fe values and the basis functions
    for(k=0;k<N_Derivatives;k++)
    {
      // get fe values of the basis functions in the quad. points
      // in the current element (pointer to pointer)
      OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct,
                                            NeededDerivatives[k]);
      // accumulate the values of the fe function 
      // do it for all quad. points
      for(j=0;j<N_Points;j++)
      {
        // get pointer to array where the values for quad. point j are
        // stored
        Orig = OrigFEValues[j];
        // initialize 
        value = value1 = 0;
        // compute linear combination over all local dof
        for(l=0;l<N_;l++)
        {
          // value += u1[l] * phi[l]
          value += FEFunctValues[l] * Orig[l];
          value1 += FEFunctValues1[l] * Orig[l];
        } // endfor l
        // copy result to position in array 

        if(k==0)
        {
          FEValue_D1u1_QuadPoint[j] = value;
          FEValue_D1u2_QuadPoint[j] = value1;
        }
        else
        {
          FEValue_D2u1_QuadPoint[j] = value;
          FEValue_D2u2_QuadPoint[j] = value1;
        }
      } // endfor j
    } // endfor k
    
    // compute part of the convolution in all local fe nodes (dof)
    // which comes from current mesh cell 
    // loop over all local dof
    for(l=0;l<N_loc_dofConv;l++)         
    {
      // get global indes of local dof
      index = DOFConv[l];          
      // initialize
      value11 = value12 = value22 = 0;
      // loop over all quadrature points 
      for(j=0;j<N_Points;j++) 
      {
        // square of the distance between local dof and quad. point
        // in original cell 
        distance_sq = ((X_orig[l]-X[j])*(X_orig[l]-X[j])+
                       (Y_orig[l]-Y[j])*(Y_orig[l]-Y[j])); 
        // compute value of GaussianFilter
        // input: delta, distance_sq
        g = GaussianFilter(delta,distance_sq);  
        // apply quadrature rule in ref. cell 
        value11 += g*AbsDetjk[j]*weights[j]*
                   ( FEValue_D1u1_QuadPoint[j]*FEValue_D1u1_QuadPoint[j]
                    +FEValue_D2u1_QuadPoint[j]*FEValue_D2u1_QuadPoint[j]);
        value12 += g*AbsDetjk[j]*weights[j]*
                   ( FEValue_D1u1_QuadPoint[j]*FEValue_D1u2_QuadPoint[j]
                    +FEValue_D2u1_QuadPoint[j]*FEValue_D2u2_QuadPoint[j]);
        value22 += g*AbsDetjk[j]*weights[j]*
                   ( FEValue_D1u2_QuadPoint[j]*FEValue_D1u2_QuadPoint[j]
                    +FEValue_D2u2_QuadPoint[j]*FEValue_D2u2_QuadPoint[j]);
      }
      // add local result to global array
      u_conv[index] += value11;
      u_conv[index+N_UConv] += value12;
      u_conv[index+2*N_UConv] += value22;
      // store global coordinates of the fe nodes (dof)
      x_conv[index] = X_orig[l];
      y_conv[index] = Y_orig[l];
      x_conv[index+N_UConv] = X_orig[l];
      y_conv[index+N_UConv] = Y_orig[l];
      
      // starting to compute the integral in the neighbour cells
      // number of edges
      N_Edges=cell->GetN_Edges();
      // loop over all edges
      for(ll=0;ll<N_Edges;ll++)                           
      {
        // get pointer to edge[l]
        joint=cell->GetJoint(ll); 
        
        // if boundary edge continue 
        if ((joint->GetType() == BoundaryEdge)||
            (joint->GetType() == IsoBoundEdge))
          continue;
        
        // get pointer to neighbour cell
        neigh=cell->GetJoint(ll)->GetNeighbour(cell);
        
        // get id of neighbour
        neigh_i = neigh->GetClipBoard();
        
        // get id of finite element in neigh mesh cell wrt the
        // finite element space of the convolution
        CurrentElementNeigh = fespaceConv->GetFE2D(neigh_i, neigh);
        
        // get pointer to basis functions of the neighbour element
        BaseFunctNeigh = BaseFuncts[CurrentElementNeigh]; 
 
        // # of basis functions, is the same as N_loc_dof
        N_Neigh = N_BaseFunct[CurrentElementNeigh]; 
     
        // find the part of the global index array where the information
        // for the neigh mesh cell are stored: DOF
        DOFNeigh = GlobalNumbersConv + BeginIndexConv[neigh_i];
        
        // check if computation on neigh cell necessary
        // for all dof
        same_dof = 0;
        for (j=0;j<N_Neigh;j++)
        {
          //cout << "N " << DOFNeigh[j] << " ind " << index << endl;
          if (DOFNeigh[j]==index)
          same_dof = 1;
        }
        // if the cells have the same dof, continue
        if (same_dof==1)
          continue;
        
        ////////////////////////////////////////////////////////////////////
        
        // if not, do something! 
        // i.e. the same thing as before
        
        // get id of finite element in neigh mesh cell wrt the finite
        // element space of the velocity
        CurrentElementNeigh = fespace->GetFE2D(neigh_i, neigh);
        
        // get pointer to basis functions of the neighbour element
        BaseFunctNeigh = BaseFuncts[CurrentElementNeigh]; 
 
        // # of basis functions, is the same as N_loc_dof
        N_Neigh = N_BaseFunct[CurrentElementNeigh]; 
     
        // find the part of the global index array where the information
        // for the neigh mesh cell are stored: DOF
        DOFNeigh = GlobalNumbers + BeginIndex[neigh_i];

        // initialize auxiliary array UsedNeigh
        memset(UsedNeigh, 0, N_FEs2D*SizeOfInt);
        
        // get id of finite element space in neigh cell 
        CurrentElementNeigh = fespace->GetFE2D(neigh_i, neigh);
        // set entry in auxiliary array to 1
        UsedNeigh[CurrentElementNeigh] = 1;
       
        // have a look how much 1 are in UsedNeigh -> N_LocalUsedElementsNeigh
        N_LocalUsedElementsNeigh = 1;
        LocalUsedElementsNeigh[0] = CurrentElementNeigh;
        
        // ####################################################################
        // calculate values on original element
        // ####################################################################
        
        // From id of fe space the situation in the referenz mesh cell is 
        // known. Now, get information in the original mesh cell. 
        // Every finite element gets automatically a quadrature rule.
        // input: N_LocalUsedElementsNeigh, LocalUsedElementsNeigh, 
        //        neigh, SecondDer
        // output : N_PointsNeigh - number of quadrature points
        //          xiNeigh - xi-values of quad. points in ref. cell
        //          etaNeigh - eta-values of quad. points in ref. cell
        //          weigthsNeigh - weights fro quad. rule in quad. points 
        //          XNeigh - X-values of quad. points in original cell  
        //          YNeigh - Y-values of quad. points in original cell  
        //          AbsDetjkNeigh - absolute value of Jacobian in quad. points
        //          xiNeigh, etaNeigh, weightsNeigh, XNeigh, YNeigh,  
        //          AbsDetjkNeigh are pointers to arrays

        TFEDatabase2D::GetOrig(N_LocalUsedElementsNeigh, LocalUsedElementsNeigh, 
                             Coll, neigh, SecondDer,
                             N_PointsNeigh, xiNeigh, etaNeigh, weightsNeigh,
                             XNeigh, YNeigh, AbsDetjkNeigh);
    
        // copy the values of the finite element function from the global
        // array (Values) to local arrays 
        for(lll=0;lll<N_Neigh;lll++)
        {
          FEFunctValuesNeigh[lll] = Values[DOFNeigh[lll]];
          FEFunctValues1Neigh[lll] = Values[DOFNeigh[lll]+N_U];
        }
        
        // compute values of fe function in the quadrature points with
        // linear combinations of the local fe values and the basis functions
        for(k=0;k<N_Derivatives;k++)
        {
          // get fe values of the basis functions in the quad. points
          // in the current element (pointer to pointer)
          OrigFEValuesNeigh = TFEDatabase2D::GetOrigElementValues(BaseFunctNeigh,
                                                       NeededDerivatives[k]);
          // accumulate the values of the fe function 
          // do it for all quad. points
          for(j=0;j<N_PointsNeigh;j++)
          {
            // get pointer to array where the values for quad. point j are
            // stored
            Orig = OrigFEValuesNeigh[j];
            // initialize 
            value = value1 = 0;
            // compute linear combination over all local dof
            for(lll=0;lll<N_Neigh;lll++)
            {
              // value += u1[l] * phi[l]
              value += FEFunctValuesNeigh[lll] * Orig[lll];
              value1 += FEFunctValues1Neigh[lll] * Orig[lll];
            } // endfor lll
            // copy result to position in array 

            if(k==0)
            {
              FEValue_D1u1_QuadPointNeigh[j] = value;
              FEValue_D1u2_QuadPointNeigh[j] = value1;
            }
            else
            {
              FEValue_D2u1_QuadPointNeigh[j] = value;
              FEValue_D2u2_QuadPointNeigh[j] = value1;
            }
          } // endfor j
        } // endfor k
        
        // compute part of the convolution in all local fe nodes (dof)
        // which comes from current mesh cell 
        // loop over all local dof
        
        // initialize
        value11 = value12 = value22 = 0;
        // loop over all quadrature points 
        for(j=0;j<N_PointsNeigh;j++) 
        { 
          // square of the distance between local dof and quad. point
          // in original cell 
          switch(TDatabase::ParamDB->INTERNAL_PERIODIC_IDENTITY)
          {
            case 0:
              distance_sq = ((X_orig[l]-XNeigh[j])*(X_orig[l]-XNeigh[j])+
                             (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j])); 
              break;
            case 1:
              exit (4711);
              break;
            case 2:
              // periodic square (1,1)^2 
                if (fabs(X_orig[l]-XNeigh[j])> 1.5)
                {
                  if (X_orig[l]-XNeigh[j]>1.5)
                    distance_sq = ((X_orig[l]-XNeigh[j]-2)*(X_orig[l]-XNeigh[j]-2)+
                                   (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j])); 
                  else
                    distance_sq = ((X_orig[l]-XNeigh[j]+2)*(X_orig[l]-XNeigh[j]+2)+
                                   (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j]));                                       
                }
                else
                  distance_sq = ((X_orig[l]-XNeigh[j])*(X_orig[l]-XNeigh[j])+
                                 (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j]));                   
              break;
          }
            
          // compute value of GaussianFilter
          // input: delta, distance_sq
          g = GaussianFilter(delta,distance_sq);
          //cout << "dist " <<  sqrt(distance_sq) << endl;
          //
          // apply quadrature rule in ref. cell 
          value11 += g*AbsDetjk[j]*weightsNeigh[j]*
                     ( FEValue_D1u1_QuadPointNeigh[j]*
                       FEValue_D1u1_QuadPointNeigh[j]
                      +FEValue_D2u1_QuadPointNeigh[j]*
                       FEValue_D2u1_QuadPointNeigh[j]);
          value12 += g*AbsDetjk[j]*weightsNeigh[j]*
                     ( FEValue_D1u1_QuadPointNeigh[j]*
                       FEValue_D1u2_QuadPointNeigh[j]
                      +FEValue_D2u1_QuadPointNeigh[j]*
                       FEValue_D2u2_QuadPointNeigh[j]);
          value22 += g*AbsDetjk[j]*weightsNeigh[j]*
                     ( FEValue_D1u2_QuadPointNeigh[j]*
                       FEValue_D1u2_QuadPointNeigh[j]
                      +FEValue_D2u2_QuadPointNeigh[j]*
                       FEValue_D2u2_QuadPointNeigh[j]);
        }
        // add local result to global array
        u_conv[index] += value11;
        u_conv[index+N_UConv] += value12; 
        u_conv[index+2*N_UConv] += value22; 
      } // end of ll loop (edges)
    } // end of l loop (local dof)
  } // endfor i (cells)
  
  // for (i=0;i<N_DOFConv; i++)
  // {
  //   cout << u_conv[i] << "  " << u_conv[i+N_UConv] << "  ";
  //   cout << u_conv[i+2*N_UConv] << "  " << u_conv[i+3*N_UConv] << endl;
  // }

  // release memory
  delete x_conv;
  delete y_conv;

} // ConvoluteSymmetricTensor

// ========================================================================
// convolute D = (grad w grad w^T) by full integration
// use D_{21} = D_{12}
// input u
// output DuTensor
// ========================================================================

void  ConvoluteSymmetricTensorFull(TFEVectFunct2D *u, TFEVectFunct2D *duTensor)

{
  int i,j,k,l,n,m, N_UsedElements, N_LocalUsedElements;
  int N_LocalUsedElementsNeigh;
  int N_Cells, N_Points, N_Parameters, N_, N_U, N_DOF, N_loc_dofConv;
  int Used[N_FEs2D], UsedNeigh[N_FEs2D], *N_BaseFunct;
  int neigh_i, N_Neigh, ll, lll, N_PointsNeigh, N_Edges, same_dof;
  FE2D LocalUsedElements[N_FEs2D], CurrentElement, CurrentElementConv;
  FE2D CurrentElementNeigh,LocalUsedElementsNeigh[N_FEs2D] ;
  TFE2D *Element;
  TNodalFunctional2D *nf;
  BaseFunct2D BaseFunct, *BaseFuncts, BaseFunctNeigh;
  TCollection *Coll;
  TBaseCell *cell, *neigh;
  TFE2D *ele;
  double *weights, *xi, *eta, *xi_ref, *eta_ref;
  double *weightsNeigh, *xiNeigh, *etaNeigh;
  double X_orig[MaxN_PointsForNodal2D], Y_orig[MaxN_PointsForNodal2D];
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D], AbsDetjkConv[MaxN_QuadPoints_2D];
  double XNeigh[MaxN_QuadPoints_2D], YNeigh[MaxN_QuadPoints_2D];
  double AbsDetjkNeigh[MaxN_QuadPoints_2D];
  RefTrans2D RefTrans;
  double *Param[MaxN_QuadPoints_2D], *aux;
  double FEValue_D1u1_QuadPoint[MaxN_QuadPoints_2D];
  double FEValue_D1u2_QuadPoint[MaxN_QuadPoints_2D];
  double FEValue_D2u1_QuadPoint[MaxN_QuadPoints_2D];
  double FEValue_D2u2_QuadPoint[MaxN_QuadPoints_2D];
  double FEValue_D1u1_QuadPointNeigh[MaxN_QuadPoints_2D];
  double FEValue_D1u2_QuadPointNeigh[MaxN_QuadPoints_2D];
  double FEValue_D2u1_QuadPointNeigh[MaxN_QuadPoints_2D];
  double FEValue_D2u2_QuadPointNeigh[MaxN_QuadPoints_2D];
  double *ExactVal[MaxN_QuadPoints_2D];
  int *DOF, ActiveBound, DirichletBound, end, last, *DOFConv, *DOFNeigh;
  double **OrigFEValues, *Orig;
  double value, value1;
  double value11, value12, value22;
  double **OrigFEValuesNeigh;
  double FEFunctValues[MaxN_BaseFunctions2D];
  double FEFunctValues1[MaxN_BaseFunctions2D];
  double FEFunctValuesNeigh[MaxN_BaseFunctions2D];
  double FEFunctValues1Neigh[MaxN_BaseFunctions2D];
  int *GlobalNumbers, *BeginIndex, *GlobalNumbersConv, *BeginIndexConv;
  double LocError[4];
  int N_UConv, N_DOFConv, index, *conv_comp;
  double *Values, *u_conv,distance_sq,h_fine,g, *u_values;
  double *x_conv, *y_conv, delta, hK;
  TJoint *joint;
  TFESpace2D *fespace, *fespaceConv; 
  TVertex *vertex0;
  double x_vertex, y_vertex;

  bool SecondDer[1] = { FALSE };
  int N_Derivatives = 2;
  MultiIndex2D NeededDerivatives[2] = { D10, D01 };

  // get fe spaces of velocity u and of covoluted tensor
  fespace = u->GetFESpace2D();
  fespaceConv = duTensor->GetFESpace2D();

  // gives a pointer to all available basis functions which
  // are described in the fedatabase
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  // gives the number of available basis functions from the
  // fedatabase
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  // get information of the numbering of the degrees of freedom
  // first: get pointer to array where a global numbering is stored
  GlobalNumbers = fespace->GetGlobalNumbers();
  // second: get pointer to array where the start of numbering for
  // each cell is stored
  BeginIndex = fespace->GetBeginIndex();
 
  // get number of dof for one velocity component N_U
  // -> total number of dof is twice N_U
  N_U = u->GetLength();
  N_DOF = 2*N_U;
  
  // get pointer to set of mesh cells which define the fe space
  Coll = fespace->GetCollection();
  // get number of mesh cells
  N_Cells = Coll->GetN_Cells();
 
  // # dof in the new fe space
  N_UConv = fespaceConv->GetN_DegreesOfFreedom();
 
  N_DOFConv = 3*N_UConv;
  // allocate memory for values of convoluted function
  u_conv = duTensor->GetValues(); 
  // initialize u_conv to 0
  memset(u_conv,0,N_DOFConv*SizeOfDouble);

  // allocate memory for values of convoluted function
  conv_comp = new int[N_UConv]; 
  // initialize conv_comp to 0
  memset(conv_comp,0,N_UConv*SizeOfInt);

  // get information of the numbering of the dof for the space of convolution
  // first: get pointer to array where a global numbering is stored
  GlobalNumbersConv = fespaceConv->GetGlobalNumbers();
  // second: get pointer to array where the start of numbering for
  // each cell is stored
  BeginIndexConv = fespaceConv->GetBeginIndex();

  // allocate memory for geometric position of nodes
  x_conv = new double[N_DOFConv]; 
  y_conv = new double[N_DOFConv]; 
   
// ########################################################################
// computing the convolution
// loop over all cells
// ########################################################################

  // direct a pointer to the values of the finite element vector
  // function
  // first part of the array : u1
  // second part of the array : u2
  Values = u->GetValues();           

  // store positition of each mesh cells in the collection of mesh 
  // cells: set ClipBoard to i
  for(i=0;i<N_Cells;i++)                       
  {
    cell = Coll->GetCell(i);
    cell->SetClipBoard(i);
  }

  // loop over all mesh cells to accumulate the values for the 
  // convolution 
  for(i=0;i<N_Cells;i++)
  {
    // get current mesh cell
    cell = Coll->GetCell(i);        
    hK = cell->GetDiameter();
    delta =  CharacteristicFilterWidth(hK);

    // ####################################################################
    // find local used finite elements on this cell
    // ####################################################################
    CurrentElement = fespace->GetFE2D(i, cell);
    LocalUsedElements[0] = CurrentElement;
    N_LocalUsedElements = 1;
    
    // ####################################################################
    // calculate values on original element
    // ####################################################################
    
    // From id of fe space the situation in the referenz mesh cell is 
    // known. Now, get information in the original mesh cell. 
    // Every finite element gets automatically a quadrature rule.
    // input: N_LocalUsedElements, LocalUsedElements, cell, SecondDer
    // output : N_Points - number of quadrature points
    //          xi - xi-values of quad. points in ref. cell
    //          eta - eta-values of quad. points in ref. cell
    //          weigths - weights fro quad. rule in quad. points 
    //          X - X-values of quad. points in original cell  
    //          Y - Y-values of quad. points in original cell  
    //          AbsDetjk - absolute value of Jacobian in quad. points
    //          xi, eta, weights, X, Y,  AbsDetjk are pointers to arrays

    TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                         Coll, cell, SecondDer,
                         N_Points, xi, eta, weights, X, Y, AbsDetjk);
 
    // compute geometric positions of the fe nodes
    // get id of finite element in current mesh cell
    CurrentElement = fespace->GetFE2D(i, cell);
    CurrentElementConv = fespaceConv->GetFE2D(i, cell);
    // get fe from its id
    Element = TFEDatabase2D::GetFE2D(CurrentElementConv);
    // get reference transformation 
    RefTrans = Element->GetRefTransID();
    // get pointer to the nodal functionals (fe nodes) of the fe 
    // (in ref. cell)
    nf = Element->GetNodalFunctional2D();
    // get number and coordinates of local dof in ref cell
    // xi_ref, eta_ref are pointers
    nf->GetPointsForAll(N_loc_dofConv, xi_ref, eta_ref);
      
    // get coordinates of fe nodes in original cell
    // input: RefTrans, N_loc_dof, xi_ref, eta_ref
    // output : X_orig, Y_orig - pointers to arrays with coordinates
    //          AbsDetjk - same as above    
    TFEDatabase2D::GetOrigFromRef(RefTrans,N_loc_dofConv, xi_ref, 
                                eta_ref, X_orig, Y_orig, AbsDetjkConv);   
   
    // get pointer to basis functions of the current element
    BaseFunct = BaseFuncts[CurrentElement];  
    // # of basis functions, is the same as N_loc_dof
    N_ = N_BaseFunct[CurrentElement];        

    // find the part of the global index array where the information
    // for the current mesh cell are stored: DOF
    DOF = GlobalNumbers + BeginIndex[i];
    DOFConv = GlobalNumbersConv + BeginIndexConv[i];
    // copy the values of the finite element function from the global
    // array (Values) to local arrays 
    for(l=0;l<N_;l++)
    {
      FEFunctValues[l] = Values[DOF[l]];
      FEFunctValues1[l] = Values[DOF[l]+N_U];
    }

    // compute values of fe function in the quadrature points with
    // linear combinations of the local fe values and the basis functions
    for(k=0;k<N_Derivatives;k++)
    {
      // get fe values of the basis functions in the quad. points
      // in the current element (pointer to pointer)
      OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct,
                                            NeededDerivatives[k]);
      // accumulate the values of the fe function 
      // do it for all quad. points
      for(j=0;j<N_Points;j++)
      {
        // get pointer to array where the values for quad. point j are
        // stored
        Orig = OrigFEValues[j];
        // initialize 
        value = value1 = 0;
        // compute linear combination over all local dof
        for(l=0;l<N_;l++)
        {
          // value += u1[l] * phi[l]
          value += FEFunctValues[l] * Orig[l];
          value1 += FEFunctValues1[l] * Orig[l];
        } // endfor l
        // copy result to position in array 

        if(k==0)
        {
          FEValue_D1u1_QuadPoint[j] = value;
          FEValue_D1u2_QuadPoint[j] = value1;
        }
        else
        {
          FEValue_D2u1_QuadPoint[j] = value;
          FEValue_D2u2_QuadPoint[j] = value1;
        }
      } // endfor j
    } // endfor k
    
    // compute part of the convolution in all local fe nodes (dof)
    // which comes from current mesh cell 
    // loop over all local dof
    for(l=0;l<N_loc_dofConv;l++)         
    {
      // get global indes of local dof
      index = DOFConv[l];          
      // check if value is already computed
      // if yes : continue
      if (conv_comp[index])
        continue;

      // set conv_comp[index]
      conv_comp[index] = 1;

      // initialize
      value11 = value12 = value22 = 0;
      // loop over all quadrature points 
      for(j=0;j<N_Points;j++) 
      {
        // square of the distance between local dof and quad. point
        // in original cell 
        distance_sq = ((X_orig[l]-X[j])*(X_orig[l]-X[j])+
                       (Y_orig[l]-Y[j])*(Y_orig[l]-Y[j])); 
        // compute value of GaussianFilter
        // input: delta, distance_sq
        g = GaussianFilter(delta,distance_sq);  
        // apply quadrature rule in ref. cell 
        value11 += g*AbsDetjk[j]*weights[j]*
                   ( FEValue_D1u1_QuadPoint[j]*FEValue_D1u1_QuadPoint[j]
                    +FEValue_D2u1_QuadPoint[j]*FEValue_D2u1_QuadPoint[j]);
        value12 += g*AbsDetjk[j]*weights[j]*
                   ( FEValue_D1u1_QuadPoint[j]*FEValue_D1u2_QuadPoint[j]
                    +FEValue_D2u1_QuadPoint[j]*FEValue_D2u2_QuadPoint[j]);
        value22 += g*AbsDetjk[j]*weights[j]*
                   ( FEValue_D1u2_QuadPoint[j]*FEValue_D1u2_QuadPoint[j]
                    +FEValue_D2u2_QuadPoint[j]*FEValue_D2u2_QuadPoint[j]);
      }
      // add local result to global array
      u_conv[index] += value11;
      u_conv[index+N_UConv] += value12;
      u_conv[index+2*N_UConv] += value22;
      // store global coordinates of the fe nodes (dof)
      x_conv[index] = X_orig[l];
      y_conv[index] = Y_orig[l];
      x_conv[index+N_UConv] = X_orig[l];
      y_conv[index+N_UConv] = Y_orig[l];
      
      // starting to compute the integral in ALL the cells
      
      // loop over all cells
      for(ll=0;ll<N_Cells;ll++)                           
      {
        // check if we are at the current mesh cell
        if (ll==i)
          continue;

        neigh = Coll->GetCell(ll); 
        // get one vertex of the neigh 
        vertex0 = neigh->GetVertex(0);
        // get coordinates of vertex
        x_vertex = vertex0->GetX();
        y_vertex = vertex0->GetY();
        // compute distance to current point where the convolution
        // is computed
        switch(TDatabase::ParamDB->INTERNAL_PERIODIC_IDENTITY)
        {
          case 0:
            distance_sq = ((X_orig[l]-x_vertex)*(X_orig[l]-x_vertex)+
                           (Y_orig[l]-y_vertex)*(Y_orig[l]-y_vertex));
            break;
          case 1:
            exit(4711);
            break;
          case 2:
            // periodic square (1,1)^2 
            if (fabs(X_orig[l]-x_vertex)> 1.5)
            {
              if (X_orig[l]-x_vertex>1.5)
                distance_sq = ((X_orig[l]-x_vertex-2)*(X_orig[l]-x_vertex-2)+
                               (Y_orig[l]-y_vertex)*(Y_orig[l]-y_vertex)); 
              else
                distance_sq = ((X_orig[l]-x_vertex+2)*(X_orig[l]-x_vertex+2)+
                               (Y_orig[l]-y_vertex)*(Y_orig[l]-y_vertex));                                       
            }
            else
              distance_sq = ((X_orig[l]-x_vertex)*(X_orig[l]-x_vertex)+
                             (Y_orig[l]-y_vertex)*(Y_orig[l]-y_vertex));                   
            break;
        }
           
        // check if neigh is close enough
        if (distance_sq > 2*delta *delta)
          continue;
        // get id of neighbour
        neigh_i = neigh->GetClipBoard();

        // get id of finite element in neigh mesh cell wrt the finite
        // element space of the velocity
        CurrentElementNeigh = fespace->GetFE2D(neigh_i, neigh);
        
        // get pointer to basis functions of the neighbour element
        BaseFunctNeigh = BaseFuncts[CurrentElementNeigh]; 
 
        // # of basis functions, is the same as N_loc_dof
        N_Neigh = N_BaseFunct[CurrentElementNeigh]; 
     
        // find the part of the global index array where the information
        // for the neigh mesh cell are stored: DOF
        DOFNeigh = GlobalNumbers + BeginIndex[neigh_i];

        // initialize auxiliary array UsedNeigh
        memset(UsedNeigh, 0, N_FEs2D*SizeOfInt);
        
        // set entry in auxiliary array to 1
        UsedNeigh[CurrentElementNeigh] = 1;
       
        // have a look how much 1 are in UsedNeigh -> N_LocalUsedElementsNeigh
        N_LocalUsedElementsNeigh = 0;
        memset(LocalUsedElementsNeigh, 0, SizeOfInt*N_FEs2D);
        j = 0;
        for(k=0;k<N_FEs2D;k++)
          if(UsedNeigh[k])
          {
            LocalUsedElementsNeigh[j] = (FE2D)k;
            j++;
          }
        N_LocalUsedElementsNeigh = j;
      
        // ####################################################################
        // calculate values on original element
        // ####################################################################
        
        // From id of fe space the situation in the referenz mesh cell is 
        // known. Now, get information in the original mesh cell. 
        // Every finite element gets automatically a quadrature rule.
        // input: N_LocalUsedElementsNeigh, LocalUsedElementsNeigh, 
        //        neigh, SecondDer
        // output : N_PointsNeigh - number of quadrature points
        //          xiNeigh - xi-values of quad. points in ref. cell
        //          etaNeigh - eta-values of quad. points in ref. cell
        //          weigthsNeigh - weights fro quad. rule in quad. points 
        //          XNeigh - X-values of quad. points in original cell  
        //          YNeigh - Y-values of quad. points in original cell  
        //          AbsDetjkNeigh - absolute value of Jacobian in quad. points
        //          xiNeigh, etaNeigh, weightsNeigh, XNeigh, YNeigh,  
        //          AbsDetjkNeigh are pointers to arrays

        TFEDatabase2D::GetOrig(N_LocalUsedElementsNeigh, LocalUsedElementsNeigh, 
                             Coll, neigh, SecondDer,
                             N_PointsNeigh, xiNeigh, etaNeigh, weightsNeigh,
                             XNeigh, YNeigh, AbsDetjkNeigh);
    
        // copy the values of the finite element function from the global
        // array (Values) to local arrays 
        for(lll=0;lll<N_Neigh;lll++)
        {
          FEFunctValuesNeigh[lll] = Values[DOFNeigh[lll]];
          FEFunctValues1Neigh[lll] = Values[DOFNeigh[lll]+N_U];
        }
        
        // compute values of fe function in the quadrature points with
        // linear combinations of the local fe values and the basis functions
        for(k=0;k<N_Derivatives;k++)
        {
          // get fe values of the basis functions in the quad. points
          // in the current element (pointer to pointer)
          OrigFEValuesNeigh = TFEDatabase2D::GetOrigElementValues(BaseFunctNeigh,
                                                       NeededDerivatives[k]);
          // accumulate the values of the fe function 
          // do it for all quad. points
          for(j=0;j<N_PointsNeigh;j++)
          {
            // get pointer to array where the values for quad. point j are
            // stored
            Orig = OrigFEValuesNeigh[j];
            // initialize 
            value = value1 = 0;
            // compute linear combination over all local dof
            for(lll=0;lll<N_Neigh;lll++)
            {
              // value += u1[l] * phi[l]
              value += FEFunctValuesNeigh[lll] * Orig[lll];
              value1 += FEFunctValues1Neigh[lll] * Orig[lll];
            } // endfor lll
            // copy result to position in array 

            if(k==0)
            {
              FEValue_D1u1_QuadPointNeigh[j] = value;
              FEValue_D1u2_QuadPointNeigh[j] = value1;
            }
            else
            {
              FEValue_D2u1_QuadPointNeigh[j] = value;
              FEValue_D2u2_QuadPointNeigh[j] = value1;
            }
          } // endfor j
        } // endfor k
        
        // compute part of the convolution in all local fe nodes (dof)
        // which comes from current mesh cell 
        // loop over all local dof
        
        // initialize
        value11 = value12 = value22 = 0;
        // loop over all quadrature points 
        for(j=0;j<N_PointsNeigh;j++) 
        { 
          // square of the distance between local dof and quad. point
          // in original cell 
          switch(TDatabase::ParamDB->INTERNAL_PERIODIC_IDENTITY)
          {
            case 0:
              distance_sq = ((X_orig[l]-XNeigh[j])*(X_orig[l]-XNeigh[j])+
                             (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j])); 
              break;
            case 1:
              exit (4711);
              break;
            case 2:
              // periodic square (1,1)^2 
                if (fabs(X_orig[l]-XNeigh[j])> 1.5)
                {
                  if (X_orig[l]-XNeigh[j]>1.5)
                    distance_sq = ((X_orig[l]-XNeigh[j]-2)*(X_orig[l]-XNeigh[j]-2)+
                                   (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j])); 
                  else
                    distance_sq = ((X_orig[l]-XNeigh[j]+2)*(X_orig[l]-XNeigh[j]+2)+
                                   (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j]));                                       
                }
                else
                  distance_sq = ((X_orig[l]-XNeigh[j])*(X_orig[l]-XNeigh[j])+
                                 (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j]));                   
              break;
          }
          // square of the distance between local dof and quad. point
          // in original cell 
          //distance_sq = ((X_orig[l]-XNeigh[j])*(X_orig[l]-XNeigh[j])+
          //               (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j])); 
          //cout << "X " << XNeigh[j] << " Y " << YNeigh[j] << endl;
          // compute value of GaussianFilter
          // input: delta, distance_sq
          g = GaussianFilter(delta,distance_sq);
          //cout << "dist " <<  sqrt(distance_sq) << endl;
          //
          // apply quadrature rule in ref. cell 
          value11 += g*AbsDetjk[j]*weightsNeigh[j]*
                     ( FEValue_D1u1_QuadPointNeigh[j]*
                       FEValue_D1u1_QuadPointNeigh[j]
                      +FEValue_D2u1_QuadPointNeigh[j]*
                       FEValue_D2u1_QuadPointNeigh[j]);
          value12 += g*AbsDetjk[j]*weightsNeigh[j]*
                     ( FEValue_D1u1_QuadPointNeigh[j]*
                       FEValue_D1u2_QuadPointNeigh[j]
                      +FEValue_D2u1_QuadPointNeigh[j]*
                       FEValue_D2u2_QuadPointNeigh[j]);
          value22 += g*AbsDetjk[j]*weightsNeigh[j]*
                     ( FEValue_D1u2_QuadPointNeigh[j]*
                       FEValue_D1u2_QuadPointNeigh[j]
                      +FEValue_D2u2_QuadPointNeigh[j]*
                       FEValue_D2u2_QuadPointNeigh[j]);
        }
        // add local result to global array
        u_conv[index] += value11;
        u_conv[index+N_UConv] += value12; 
        u_conv[index+2*N_UConv] += value22; 
      } // end of ll loop (edges)
    } // end of l loop (local dof)
  } // endfor i (cells)
  
  /* for (i=0;i<N_UConv; i++)
    {
      cout << x_conv[i] << " " << y_conv[i];
      cout << "    " << Values[i] << " "<< Values[i+N_U] << "     ";
      cout << u_conv[i] << "  " << u_conv[i+N_UConv] << "  ";
      cout << u_conv[i+2*N_UConv] << endl;
      }*/

  // release memory
  delete x_conv;
  delete y_conv;
  delete conv_comp;

} // ConvoluteSymmetricTensorFull
// ========================================================================
// convolute D = (grad w grad w^T) by full integration
// use D_{21} = D_{12}
// input u
// output DuTensor for all dof at bounaries 
// ========================================================================
/*
void  ConvoluteSymmetricTensorFullBoundary(TFEVectFunct2D *u, TFEVectFunct2D *duTensor)

{
  int i,j,k,l,n,m, N_UsedElements, N_LocalUsedElements;
  int N_LocalUsedElementsNeigh;
  int N_Cells, N_Points, N_Parameters, N_, N_U, N_DOF, N_loc_dofConv;
  int Used[N_FEs2D], UsedNeigh[N_FEs2D], *N_BaseFunct;
  int neigh_i, N_Neigh, ll, lll, N_PointsNeigh, N_Edges, same_dof;
  FE2D LocalUsedElements[N_FEs2D], CurrentElement, CurrentElementConv;
  FE2D CurrentElementNeigh,LocalUsedElementsNeigh[N_FEs2D] ;
  TFE2D *Element;
  TNodalFunctional2D *nf;
  BaseFunct2D BaseFunct, *BaseFuncts, BaseFunctNeigh;
  TCollection *Coll;
  TBaseCell *cell, *neigh;
  TFE2D *ele;
  double *weights, *xi, *eta, *xi_ref, *eta_ref;
  double *weightsNeigh, *xiNeigh, *etaNeigh;
  double X_orig[MaxN_PointsForNodal2D], Y_orig[MaxN_PointsForNodal2D];
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D], AbsDetjkConv[MaxN_QuadPoints_2D];
  double XNeigh[MaxN_QuadPoints_2D], YNeigh[MaxN_QuadPoints_2D];
  double AbsDetjkNeigh[MaxN_QuadPoints_2D];
  RefTrans2D RefTrans;
  double *Param[MaxN_QuadPoints_2D], *aux;
  double FEValue_D1u1_QuadPoint[MaxN_QuadPoints_2D];
  double FEValue_D1u2_QuadPoint[MaxN_QuadPoints_2D];
  double FEValue_D2u1_QuadPoint[MaxN_QuadPoints_2D];
  double FEValue_D2u2_QuadPoint[MaxN_QuadPoints_2D];
  double FEValue_D1u1_QuadPointNeigh[MaxN_QuadPoints_2D];
  double FEValue_D1u2_QuadPointNeigh[MaxN_QuadPoints_2D];
  double FEValue_D2u1_QuadPointNeigh[MaxN_QuadPoints_2D];
  double FEValue_D2u2_QuadPointNeigh[MaxN_QuadPoints_2D];
  double *ExactVal[MaxN_QuadPoints_2D];
  int *DOF, ActiveBound, DirichletBound, end, last, *DOFConv, *DOFNeigh;
  double **OrigFEValues, *Orig;
  double value, value1;
  double value11, value12, value22;
  double **OrigFEValuesNeigh;
  double FEFunctValues[MaxN_BaseFunctions2D];
  double FEFunctValues1[MaxN_BaseFunctions2D];
  double FEFunctValuesNeigh[MaxN_BaseFunctions2D];
  double FEFunctValues1Neigh[MaxN_BaseFunctions2D];
  int *GlobalNumbers, *BeginIndex, *GlobalNumbersConv, *BeginIndexConv;
  double LocError[4];
  int N_UConv, N_DOFConv, index, *conv_comp;
  double *Values, *u_conv,distance_sq,h_fine,g, *u_values;
  double *x_conv, *y_conv, delta, hK;
  TJoint *joint;
  TFESpace2D *fespace, *fespaceConv; 
  TVertex *vertex0;
  double x_vertex, y_vertex;

  bool SecondDer[1] = { FALSE };
  int N_Derivatives = 2;
  MultiIndex2D NeededDerivatives[2] = { D10, D01 };

  // get fe spaces of velocity u and of covoluted tensor
  fespace = u->GetFESpace2D();
  fespaceConv = duTensor->GetFESpace2D();

  // gives a pointer to all available basis functions which
  // are described in the fedatabase
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  // gives the number of available basis functions from the
  // fedatabase
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  // get information of the numbering of the degrees of freedom
  // first: get pointer to array where a global numbering is stored
  GlobalNumbers = fespace->GetGlobalNumbers();
  // second: get pointer to array where the start of numbering for
  // each cell is stored
  BeginIndex = fespace->GetBeginIndex();
 
  // get number of dof for one velocity component N_U
  // -> total number of dof is twice N_U
  N_U = u->GetLength();
  N_DOF = 2*N_U;
  
  // get pointer to set of mesh cells which define the fe space
  Coll = fespace->GetCollection();
  // get number of mesh cells
  N_Cells = Coll->GetN_Cells();
 
  // # dof in the new fe space
  N_UConv = fespaceConv->GetN_DegreesOfFreedom();
 
  N_DOFConv = 3*N_UConv;
  // allocate memory for values of convoluted function
  u_conv = duTensor->GetValues(); 
  // initialize u_conv to 0
  memset(u_conv,0,N_DOFConv*SizeOfDouble);

  // allocate memory for values of convoluted function
  conv_comp = new int[N_UConv]; 
  // initialize conv_comp to 0
  memset(conv_comp,0,N_UConv*SizeOfInt);

  // get information of the numbering of the dof for the space of convolution
  // first: get pointer to array where a global numbering is stored
  GlobalNumbersConv = fespaceConv->GetGlobalNumbers();
  // second: get pointer to array where the start of numbering for
  // each cell is stored
  BeginIndexConv = fespaceConv->GetBeginIndex();

  // allocate memory for geometric position of nodes
  x_conv = new double[N_DOFConv]; 
  y_conv = new double[N_DOFConv]; 
   
// ########################################################################
// computing the convolution
// loop over all cells
// ########################################################################

  // direct a pointer to the values of the finite element vector
  // function
  // first part of the array : u1
  // second part of the array : u2
  Values = u->GetValues();           

  // store positition of each mesh cells in the collection of mesh 
  // cells: set ClipBoard to i
  for(i=0;i<N_Cells;i++)                       
  {
    cell = Coll->GetCell(i);
    cell->SetClipBoard(i);
  }

  // loop over all mesh cells to accumulate the values for the 
  // convolution 
  for(i=0;i<N_Cells;i++)
  {
    // get current mesh cell
    cell = Coll->GetCell(i);        
    N_Joints = cell->GetN_Edges();
    for(m=0;m<N_Joints;m++)
    {
      joint = cell->GetJoint(m);
      has_boundedge = 0;
      if( joint->GetType() == BoundaryEdge ||
          joint->GetType() == IsoBoundEdge)
      {
         has_boundedge++;+
         break;
      }
    }
    // there is no edge on the boundary 
    if (!has_boundedge)
      continue;
    hK = cell->GetDiameter();
    delta =  CharacteristicFilterWidth(hK);

    // ####################################################################
    // find local used finite elements on this cell
    // ####################################################################
    CurrentElement = fespace->GetFE2D(i, cell);
    LocalUsedElements[0] = CurrentElement;
    N_LocalUsedElements = 1;
    
    // ####################################################################
    // calculate values on original element
    // ####################################################################
    
    // From id of fe space the situation in the referenz mesh cell is 
    // known. Now, get information in the original mesh cell. 
    // Every finite element gets automatically a quadrature rule.
    // input: N_LocalUsedElements, LocalUsedElements, cell, SecondDer
    // output : N_Points - number of quadrature points
    //          xi - xi-values of quad. points in ref. cell
    //          eta - eta-values of quad. points in ref. cell
    //          weigths - weights fro quad. rule in quad. points 
    //          X - X-values of quad. points in original cell  
    //          Y - Y-values of quad. points in original cell  
    //          AbsDetjk - absolute value of Jacobian in quad. points
    //          xi, eta, weights, X, Y,  AbsDetjk are pointers to arrays

    TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                         cell, SecondDer,
                         N_Points, xi, eta, weights, X, Y, AbsDetjk);
 
    // compute geometric positions of the fe nodes
    // get id of finite element in current mesh cell
    CurrentElement = fespace->GetFE2D(i, cell);
    CurrentElementConv = fespaceConv->GetFE2D(i, cell);
    // get fe from its id
    Element = TFEDatabase2D::GetFE2D(CurrentElementConv);
    // get reference transformation 
    RefTrans = Element->GetRefTransID();
    // get pointer to the nodal functionals (fe nodes) of the fe 
    // (in ref. cell)
    nf = Element->GetNodalFunctional2D();
    // get number and coordinates of local dof in ref cell
    // xi_ref, eta_ref are pointers
    nf->GetPointsForAll(N_loc_dofConv, xi_ref, eta_ref);
      
    // get coordinates of fe nodes in original cell
    // input: RefTrans, N_loc_dof, xi_ref, eta_ref
    // output : X_orig, Y_orig - pointers to arrays with coordinates
    //          AbsDetjk - same as above    
    TFEDatabase2D::GetOrigFromRef(RefTrans,N_loc_dofConv, xi_ref, 
                                eta_ref, X_orig, Y_orig, AbsDetjkConv);   
   
    // get pointer to basis functions of the current element
    BaseFunct = BaseFuncts[CurrentElement];  
    // # of basis functions, is the same as N_loc_dof
    N_ = N_BaseFunct[CurrentElement];        

    // find the part of the global index array where the information
    // for the current mesh cell are stored: DOF
    DOF = GlobalNumbers + BeginIndex[i];
    DOFConv = GlobalNumbersConv + BeginIndexConv[i];
    // copy the values of the finite element function from the global
    // array (Values) to local arrays 
    for(l=0;l<N_;l++)
    {
      FEFunctValues[l] = Values[DOF[l]];
      FEFunctValues1[l] = Values[DOF[l]+N_U];
    }

    // compute values of fe function in the quadrature points with
    // linear combinations of the local fe values and the basis functions
    for(k=0;k<N_Derivatives;k++)
    {
      // get fe values of the basis functions in the quad. points
      // in the current element (pointer to pointer)
      OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct,
                                            NeededDerivatives[k]);
      // accumulate the values of the fe function 
      // do it for all quad. points
      for(j=0;j<N_Points;j++)
      {
        // get pointer to array where the values for quad. point j are
        // stored
        Orig = OrigFEValues[j];
        // initialize 
        value = value1 = 0;
        // compute linear combination over all local dof
        for(l=0;l<N_;l++)
        {
          // value += u1[l] * phi[l]
          value += FEFunctValues[l] * Orig[l];
          value1 += FEFunctValues1[l] * Orig[l];
        } // endfor l
        // copy result to position in array 

        if(k==0)
        {
          FEValue_D1u1_QuadPoint[j] = value;
          FEValue_D1u2_QuadPoint[j] = value1;
        }
        else
        {
          FEValue_D2u1_QuadPoint[j] = value;
          FEValue_D2u2_QuadPoint[j] = value1;
        }
      } // endfor j
    } // endfor k
    
    // compute part of the convolution in all local fe nodes (dof)
    // which comes from current mesh cell 
    // loop over all local dof
    for(l=0;l<N_loc_dofConv;l++)         
    {
      // get global indes of local dof
      index = DOFConv[l];          
      // check if value is already computed
      // if yes : continue
      if (conv_comp[index])
        continue;

      // set conv_comp[index]
      conv_comp[index] = 1;

      // initialize
      value11 = value12 = value22 = 0;
      // loop over all quadrature points 
      for(j=0;j<N_Points;j++) 
      {
        // square of the distance between local dof and quad. point
        // in original cell 
        distance_sq = ((X_orig[l]-X[j])*(X_orig[l]-X[j])+
                       (Y_orig[l]-Y[j])*(Y_orig[l]-Y[j])); 
        // compute value of GaussianFilter
        // input: delta, distance_sq
        g = GaussianFilter(delta,distance_sq);  
        // apply quadrature rule in ref. cell 
        value11 += g*AbsDetjk[j]*weights[j]*
                   ( FEValue_D1u1_QuadPoint[j]*FEValue_D1u1_QuadPoint[j]
                    +FEValue_D2u1_QuadPoint[j]*FEValue_D2u1_QuadPoint[j]);
        value12 += g*AbsDetjk[j]*weights[j]*
                   ( FEValue_D1u1_QuadPoint[j]*FEValue_D1u2_QuadPoint[j]
                    +FEValue_D2u1_QuadPoint[j]*FEValue_D2u2_QuadPoint[j]);
        value22 += g*AbsDetjk[j]*weights[j]*
                   ( FEValue_D1u2_QuadPoint[j]*FEValue_D1u2_QuadPoint[j]
                    +FEValue_D2u2_QuadPoint[j]*FEValue_D2u2_QuadPoint[j]);
      }
      // add local result to global array
      u_conv[index] += value11;
      u_conv[index+N_UConv] += value12;
      u_conv[index+2*N_UConv] += value22;
      // store global coordinates of the fe nodes (dof)
      x_conv[index] = X_orig[l];
      y_conv[index] = Y_orig[l];
      x_conv[index+N_UConv] = X_orig[l];
      y_conv[index+N_UConv] = Y_orig[l];
      
      // starting to compute the integral in ALL the cells
      
      // loop over all cells
      for(ll=0;ll<N_Cells;ll++)                           
      {
        // check if we are at the current mesh cell
        if (ll==i)
          continue;

        neigh = Coll->GetCell(ll); 
        // get one vertex of the neigh 
        vertex0 = neigh->GetVertex(0);
        // get coordinates of vertex
        x_vertex = vertex0->GetX();
        y_vertex = vertex0->GetY();
        // compute distance to current point where the convolution
        // is computed
        switch(TDatabase::ParamDB->INTERNAL_PERIODIC_IDENTITY)
        {
          case 0:
            distance_sq = ((X_orig[l]-x_vertex)*(X_orig[l]-x_vertex)+
                           (Y_orig[l]-y_vertex)*(Y_orig[l]-y_vertex));
            break;
          case 1:
            exit(4711);
            break;
          case 2:
            // periodic square (1,1)^2 
            if (fabs(X_orig[l]-x_vertex)> 1.5)
            {
              if (X_orig[l]-x_vertex>1.5)
                distance_sq = ((X_orig[l]-x_vertex-2)*(X_orig[l]-x_vertex-2)+
                               (Y_orig[l]-y_vertex)*(Y_orig[l]-y_vertex)); 
              else
                distance_sq = ((X_orig[l]-x_vertex+2)*(X_orig[l]-x_vertex+2)+
                               (Y_orig[l]-y_vertex)*(Y_orig[l]-y_vertex));                                       
            }
            else
              distance_sq = ((X_orig[l]-x_vertex)*(X_orig[l]-x_vertex)+
                             (Y_orig[l]-y_vertex)*(Y_orig[l]-y_vertex));                   
            break;
        }
           
        // check if neigh is close enough
        if (dstance_sq > 2*delta *delta)
          continue;
        // get id of neighbour
        neigh_i = neigh->GetClipBoard();

        // get id of finite element in neigh mesh cell wrt the finite
        // element space of the velocity
        CurrentElementNeigh = fespace->GetFE2D(neigh_i, neigh);
        
        // get pointer to basis functions of the neighbour element
        BaseFunctNeigh = BaseFuncts[CurrentElementNeigh]; 
 
        // # of basis functions, is the same as N_loc_dof
        N_Neigh = N_BaseFunct[CurrentElementNeigh]; 
     
        // find the part of the global index array where the information
        // for the neigh mesh cell are stored: DOF
        DOFNeigh = GlobalNumbers + BeginIndex[neigh_i];

        // initialize auxiliary array UsedNeigh
        memset(UsedNeigh, 0, N_FEs2D*SizeOfInt);
        
        // set entry in auxiliary array to 1
        UsedNeigh[CurrentElementNeigh] = 1;
       
        // have a look how much 1 are in UsedNeigh -> N_LocalUsedElementsNeigh
        N_LocalUsedElementsNeigh = 0;
        memset(LocalUsedElementsNeigh, 0, SizeOfInt*N_FEs2D);
        j = 0;
        for(k=0;k<N_FEs2D;k++)
          if(UsedNeigh[k])
          {
            LocalUsedElementsNeigh[j] = (FE2D)k;
            j++;
          }
        N_LocalUsedElementsNeigh = j;
      
        // ####################################################################
        // calculate values on original element
        // ####################################################################
        
        // From id of fe space the situation in the referenz mesh cell is 
        // known. Now, get information in the original mesh cell. 
        // Every finite element gets automatically a quadrature rule.
        // input: N_LocalUsedElementsNeigh, LocalUsedElementsNeigh, 
        //        neigh, SecondDer
        // output : N_PointsNeigh - number of quadrature points
        //          xiNeigh - xi-values of quad. points in ref. cell
        //          etaNeigh - eta-values of quad. points in ref. cell
        //          weigthsNeigh - weights fro quad. rule in quad. points 
        //          XNeigh - X-values of quad. points in original cell  
        //          YNeigh - Y-values of quad. points in original cell  
        //          AbsDetjkNeigh - absolute value of Jacobian in quad. points
        //          xiNeigh, etaNeigh, weightsNeigh, XNeigh, YNeigh,  
        //          AbsDetjkNeigh are pointers to arrays

        TFEDatabase2D::GetOrig(N_LocalUsedElementsNeigh, LocalUsedElementsNeigh, 
                             neigh, SecondDer,
                             N_PointsNeigh, xiNeigh, etaNeigh, weightsNeigh,
                             XNeigh, YNeigh, AbsDetjkNeigh);
    
        // copy the values of the finite element function from the global
        // array (Values) to local arrays 
        for(lll=0;lll<N_Neigh;lll++)
        {
          FEFunctValuesNeigh[lll] = Values[DOFNeigh[lll]];
          FEFunctValues1Neigh[lll] = Values[DOFNeigh[lll]+N_U];
        }
        
        // compute values of fe function in the quadrature points with
        // linear combinations of the local fe values and the basis functions
        for(k=0;k<N_Derivatives;k++)
        {
          // get fe values of the basis functions in the quad. points
          // in the current element (pointer to pointer)
          OrigFEValuesNeigh = TFEDatabase2D::GetOrigElementValues(BaseFunctNeigh,
                                                       NeededDerivatives[k]);
          // accumulate the values of the fe function 
          // do it for all quad. points
          for(j=0;j<N_PointsNeigh;j++)
          {
            // get pointer to array where the values for quad. point j are
            // stored
            Orig = OrigFEValuesNeigh[j];
            // initialize 
            value = value1 = 0;
            // compute linear combination over all local dof
            for(lll=0;lll<N_Neigh;lll++)
            {
              // value += u1[l] * phi[l]
              value += FEFunctValuesNeigh[lll] * Orig[lll];
              value1 += FEFunctValues1Neigh[lll] * Orig[lll];
            } // endfor lll
            // copy result to position in array 

            if(k==0)
            {
              FEValue_D1u1_QuadPointNeigh[j] = value;
              FEValue_D1u2_QuadPointNeigh[j] = value1;
            }
            else
            {
              FEValue_D2u1_QuadPointNeigh[j] = value;
              FEValue_D2u2_QuadPointNeigh[j] = value1;
            }
          } // endfor j
        } // endfor k
        
        // compute part of the convolution in all local fe nodes (dof)
        // which comes from current mesh cell 
        // loop over all local dof
        
        // initialize
        value11 = value12 = value22 = 0;
        // loop over all quadrature points 
        for(j=0;j<N_PointsNeigh;j++) 
        { 
          // square of the distance between local dof and quad. point
          // in original cell 
          switch(TDatabase::ParamDB->INTERNAL_PERIODIC_IDENTITY)
          {
            case 0:
              distance_sq = ((X_orig[l]-XNeigh[j])*(X_orig[l]-XNeigh[j])+
                             (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j])); 
              break;
            case 1:
              exit (4711);
              break;
            case 2:
              // periodic square (1,1)^2 
                if (fabs(X_orig[l]-XNeigh[j])> 1.5)
                {
                  if (X_orig[l]-XNeigh[j]>1.5)
                    distance_sq = ((X_orig[l]-XNeigh[j]-2)*(X_orig[l]-XNeigh[j]-2)+
                                   (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j])); 
                  else
                    distance_sq = ((X_orig[l]-XNeigh[j]+2)*(X_orig[l]-XNeigh[j]+2)+
                                   (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j]));                                       
                }
                else
                  distance_sq = ((X_orig[l]-XNeigh[j])*(X_orig[l]-XNeigh[j])+
                                 (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j]));                   
              break;
          }
          // square of the distance between local dof and quad. point
          // in original cell 
          //distance_sq = ((X_orig[l]-XNeigh[j])*(X_orig[l]-XNeigh[j])+
          //               (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j])); 
          //cout << "X " << XNeigh[j] << " Y " << YNeigh[j] << endl;
          // compute value of GaussianFilter
          // input: delta, distance_sq
          g = GaussianFilter(delta,distance_sq);
          //cout << "dist " <<  sqrt(distance_sq) << endl;
          //
          // apply quadrature rule in ref. cell 
          value11 += g*AbsDetjk[j]*weightsNeigh[j]*
                     ( FEValue_D1u1_QuadPointNeigh[j]*
                       FEValue_D1u1_QuadPointNeigh[j]
                      +FEValue_D2u1_QuadPointNeigh[j]*
                       FEValue_D2u1_QuadPointNeigh[j]);
          value12 += g*AbsDetjk[j]*weightsNeigh[j]*
                     ( FEValue_D1u1_QuadPointNeigh[j]*
                       FEValue_D1u2_QuadPointNeigh[j]
                      +FEValue_D2u1_QuadPointNeigh[j]*
                       FEValue_D2u2_QuadPointNeigh[j]);
          value22 += g*AbsDetjk[j]*weightsNeigh[j]*
                     ( FEValue_D1u2_QuadPointNeigh[j]*
                       FEValue_D1u2_QuadPointNeigh[j]
                      +FEValue_D2u2_QuadPointNeigh[j]*
                       FEValue_D2u2_QuadPointNeigh[j]);
        }
        // add local result to global array
        u_conv[index] += value11;
        u_conv[index+N_UConv] += value12; 
        u_conv[index+2*N_UConv] += value22; 
      } // end of ll loop (edges)
    } // end of l loop (local dof)
  } // endfor i (cells)
  
  // for (i=0;i<N_UConv; i++)
  //  {
  //    cout << x_conv[i] << " " << y_conv[i];
  //    cout << "    " << Values[i] << " "<< Values[i+N_U] << "     ";
  //    cout << u_conv[i] << "  " << u_conv[i+N_UConv] << "  ";
  //    cout << u_conv[i+2*N_UConv] << endl;
  //    }

  // release memory
  delete x_conv;
  delete y_conv;
  delete conv_comp;

} // ConvoluteSymmetricTensorFullBoundary
*/
#endif

#ifdef __3D__

double GaussianFilter3D(double delta, double dist_sq)
{
  double a,b;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;
  
  b = gamma/(Pi*delta*delta);
  a = b*sqrt(b)*exp(-gamma * dist_sq/(delta*delta));
  return(a);
}

void  ConvoluteVelocity3D(TFEVectFunct3D *u, TFEVectFunct3D *uConv)
{
  int i,j,k,l,n,m, N_UsedElements, N_LocalUsedElements;
  int N_LocalUsedElementsNeigh;
  int N_Cells, N_Points, N_Parameters, N_, N_U, N_DOF, N_loc_dofConv;
  int Used[N_FEs3D], UsedNeigh[N_FEs3D], *N_BaseFunct;
  int neigh_i, N_Neigh, ll, lll, N_PointsNeigh, N_Faces, same_dof;
  FE3D LocalUsedElements[N_FEs3D], CurrentElement, CurrentElementConv;
  FE3D CurrentElementNeigh,LocalUsedElementsNeigh[N_FEs3D] ;
  TFE3D *Element;
  TNodalFunctional3D *nf;
  BaseFunct3D BaseFunct, *BaseFuncts, BaseFunctNeigh;
  TCollection *Coll;
  TBaseCell *cell, *neigh;
  TFE3D *ele;
  double *weights, *xi, *eta, *zeta, *xi_ref, *eta_ref, *zeta_ref;
  double *weightsNeigh, *xiNeigh, *etaNeigh, *zetaNeigh;
  double X_orig[MaxN_PointsForNodal3D], Y_orig[MaxN_PointsForNodal3D];
  double Z_orig[MaxN_PointsForNodal3D];
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D], AbsDetjkConv[MaxN_QuadPoints_3D];
  double XNeigh[MaxN_QuadPoints_3D], YNeigh[MaxN_QuadPoints_3D];
  double ZNeigh[MaxN_QuadPoints_3D];
  double AbsDetjkNeigh[MaxN_QuadPoints_3D];
  RefTrans3D RefTrans;
  double *Param[MaxN_QuadPoints_3D], *aux;
  double FEValue_u1_QuadPoint[MaxN_QuadPoints_3D];
  double FEValue_u2_QuadPoint[MaxN_QuadPoints_3D];
  double FEValue_u3_QuadPoint[MaxN_QuadPoints_3D];
  double FEValue_u1_QuadPointNeigh[MaxN_QuadPoints_3D];
  double FEValue_u2_QuadPointNeigh[MaxN_QuadPoints_3D];
  double FEValue_u3_QuadPointNeigh[MaxN_QuadPoints_3D];
  double *ExactVal[MaxN_QuadPoints_3D];
  int *DOF, ActiveBound, DirichletBound, end, last, *DOFConv, *DOFNeigh;
  double **OrigFEValues, *Orig, value, value1, value2;
  double **OrigFEValuesNeigh;
  double FEFunctValues[MaxN_BaseFunctions3D];
  double FEFunctValues1[MaxN_BaseFunctions3D];
  double FEFunctValues2[MaxN_BaseFunctions3D];
  double FEFunctValuesNeigh[MaxN_BaseFunctions3D];
  double FEFunctValues1Neigh[MaxN_BaseFunctions3D];
  double FEFunctValues2Neigh[MaxN_BaseFunctions3D];
  int *GlobalNumbers, *BeginIndex, *GlobalNumbersConv, *BeginIndexConv;
  double LocError[4];
  double hK,delta;
  bool *SecondDer;
  int n_fespaces = 1, N_Derivatives = 1, index;
  int N_UConv, N_DOFConv;
  double *Values, *u_conv,distance_sq,h_fine,g, *u_values;
  double *x_conv, *y_conv, *z_conv;
  TJoint *joint;
  TFESpace3D *fespace, *fespaceConv; 

  TDatabase::ParamDB->INTERNAL_QUAD_RULE=1;

  // gives an array where the needed derivatives are described
  // defined in NavierStokes.h 
  // here we get as result {D00}
  MultiIndex3D NeededDerivatives[1] = { D000 };

  // get fe spaces of velocity u and of covoluted tensor
  fespace = u->GetFESpace3D();
  fespaceConv = uConv->GetFESpace3D();

  // gives a pointer to all available basis functions which
  // are described in the fedatabase
  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  // gives the number of available basis functions from the
  // fedatabase
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  // say that we don't need second derivatives
  // first: allocate memory for an array of length n_fespaces 
  // second: set all entries to FALSE
  SecondDer = new bool[n_fespaces];
  for(i=0;i<n_fespaces;i++)
    SecondDer[i] = FALSE;

  // get information of the numbering of the degrees of freedom
  // first: get pointer to array where a global numbering is stored
  GlobalNumbers = fespace->GetGlobalNumbers();
  // second: get pointer to array where the start of numbering for
  // each cell is stored
  BeginIndex = fespace->GetBeginIndex();
 
  // get number of dof for one velocity component N_U
  // -> total number of dof is twice N_U
  N_U = u->GetLength();
  N_DOF = 3*N_U;
  
  // get pointer to set of mesh cells which define the fe space
  Coll = fespace->GetCollection();
 
  // # dof in the new fe space
  N_UConv = fespaceConv->GetN_DegreesOfFreedom();
 
  N_DOFConv = 3*N_UConv;
  // allocate memory for values of convuluted function
  u_conv = uConv->GetValues(); 
  // initialize u_conv to 0
  memset(u_conv,0,N_DOFConv*SizeOfDouble);

  // get information of the numbering of the dof for the space of convolution
  // first: get pointer to array where a global numbering is stored
  GlobalNumbersConv = fespaceConv->GetGlobalNumbers();
  // second: get pointer to array where the start of numbering for
  // each cell is stored
  BeginIndexConv = fespaceConv->GetBeginIndex();

  // allocate memory for geometric position of nodes
  x_conv = new double[N_DOFConv]; 
  y_conv = new double[N_DOFConv]; 
  z_conv = new double[N_DOFConv]; 
   
// ########################################################################
// computing the convolution
// loop over all cells
// ########################################################################

  // get pointer to set of mesh cells which define the fe space
  Coll = fespace->GetCollection(); 
  // get number of mesh cells
  N_Cells = Coll->GetN_Cells();

  // direct a pointer to the values of the finite element vector
  // function
  // first part of the array : u1
  // second part of the array : u2
  Values = u->GetValues();           

  // store positition of each mesh cells in the collection of mesh 
  // cells: set ClipBoard to i
  for(i=0;i<N_Cells;i++)                       
  {
    cell = Coll->GetCell(i);
    cell->SetClipBoard(i);
  }

  // loop over all mesh cells to accumulate the values for the 
  // convolution 
  for(i=0;i<N_Cells;i++)
  {
    // get current mesh cell
    cell = Coll->GetCell(i);        
    // get cell diameter (longest distance between to vertices of the cell)
    hK = cell->GetDiameter();       
    delta = CharacteristicFilterWidth(hK);
 
    // ####################################################################
    // find local used finite elements on this cell
    // ####################################################################
    
    // have a look how much 1 are in Used ->   N_LocalUsedElements
    CurrentElement = fespace->GetFE3D(i, cell);
    LocalUsedElements[0] = CurrentElement;
    N_LocalUsedElements = 1;
    
    // ####################################################################
    // calculate values on original element
    // ####################################################################
    
    // From id of fe space the situation in the referenz mesh cell is 
    // known. Now, get information in the original mesh cell. 
    // Every finite element gets automatically a quadrature rule.
    // input: N_LocalUsedElements, LocalUsedElements, cell, SecondDer
    // output : N_Points - number of quadrature points
    //          xi - xi-values of quad. points in ref. cell
    //          eta - eta-values of quad. points in ref. cell
    //          weigths - weights fro quad. rule in quad. points 
    //          X - X-values of quad. points in original cell  
    //          Y - Y-values of quad. points in original cell  
    //          AbsDetjk - absolute value of Jacobian in quad. points
    //          xi, eta, weights, X, Y,  AbsDetjk are pointers to arrays

    TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, zeta,
                           weights, X, Y, Z, AbsDetjk);
 
    //cout << " Points " << N_Points ;
    // compute geometric positions of the fe nodes
    // get id of finite element in current mesh cell
    CurrentElement = fespace->GetFE3D(i, cell);
    CurrentElementConv = fespaceConv->GetFE3D(i, cell);
    // get fe from its id
    Element = TFEDatabase3D::GetFE3D(CurrentElementConv);
    // get reference transformation 
    RefTrans = Element->GetRefTransID();
    // get pointer to the nodal functionals (fe nodes) of the fe 
    // (in ref. cell)
    nf = Element->GetNodalFunctional3D();
    // get number and coordinates of local dof in ref cell
    // xi_ref, eta_ref, zeta_ref are pointers
    nf->GetPointsForAll(N_loc_dofConv, xi_ref, eta_ref, zeta_ref);
      
    // get coordinates of fe nodes in original cell
    // input: RefTrans, N_loc_dof, xi_ref, eta_ref
    // output : X_orig, Y_orig - pointers to arrays with coordinates
    //          AbsDetjk - same as above    
    TFEDatabase3D::GetOrigFromRef(RefTrans,N_loc_dofConv, 
                                  xi_ref, eta_ref, zeta_ref,
                                  X_orig, Y_orig, Z_orig, AbsDetjkConv);   
   
    // get pointer to basis functions of the current element
    BaseFunct = BaseFuncts[CurrentElement];  
    // # of basis functions, is the same as N_loc_dof
    N_ = N_BaseFunct[CurrentElement];        

    // find the part of the global index array where the information
    // for the current mesh cell are stored: DOF
    DOF = GlobalNumbers + BeginIndex[i];
    DOFConv = GlobalNumbersConv + BeginIndexConv[i];
    // copy the values of the finite element function from the global
    // array (Values) to local arrays 
    for(l=0;l<N_;l++)
    {
      FEFunctValues[l] = Values[DOF[l]];
      FEFunctValues1[l] = Values[DOF[l]+N_U];
      FEFunctValues2[l] = Values[DOF[l]+2*N_U];
    }

    // compute values of fe function in the quadrature points with
    // linear combinations of the local fe values and the basis functions
    for(k=0;k<N_Derivatives;k++)
    {
      // get fe values of the basis functions in the quad. points
      // in the current element (pointer to pointer)
      OrigFEValues = TFEDatabase3D::GetOrigElementValues(BaseFunct,
                                            NeededDerivatives[k]);
      // accumulate the values of the fe function 
      // do it for all quad. points
      for(j=0;j<N_Points;j++)
      {
        // get pointer to array where the values for quad. point j are
        // stored
        Orig = OrigFEValues[j];
        // initialize 
        value = value1 = value2 = 0;
        // compute linear combination over all local dof
        for(l=0;l<N_;l++)
        {
          // value += u1[l] * phi[l]
          value += FEFunctValues[l] * Orig[l];
          value1 += FEFunctValues1[l] * Orig[l];
          value2 += FEFunctValues2[l] * Orig[l];
        } // endfor l
        // copy result to position in array 
        FEValue_u1_QuadPoint[j] = value;
        FEValue_u2_QuadPoint[j] = value1;
        FEValue_u3_QuadPoint[j] = value2;
      } // endfor j
    } // endfor k
    
    // compute part of the convolution in all local fe nodes (dof)
    // which comes from current mesh cell 
    // loop over all local dof
    for(l=0;l<N_loc_dofConv;l++)         
    {
      // get global indes of local dof
      index = DOFConv[l];          
      // initialize
      value = value1 = value2 = 0;
      // loop over all quadrature points 
      for(j=0;j<N_Points;j++) 
      {
        // square of the distance between local dof and quad. point
        // in original cell 
        distance_sq = ((X_orig[l]-X[j])*(X_orig[l]-X[j])+
                       (Y_orig[l]-Y[j])*(Y_orig[l]-Y[j])+
                       (Z_orig[l]-Z[j])*(Z_orig[l]-Z[j])); 
        // compute value of GaussianFilter
        // input: delta, distance_sq
        g = GaussianFilter3D(delta,distance_sq);  
        // apply quadrature rule in ref. cell 
        value += g*AbsDetjk[j]*weights[j]*FEValue_u1_QuadPoint[j];
        value1 += g*AbsDetjk[j]*weights[j]*FEValue_u2_QuadPoint[j];
        value2 += g*AbsDetjk[j]*weights[j]*FEValue_u3_QuadPoint[j];
      }
      // add local result to global array
      u_conv[index] += value;
      u_conv[index+N_UConv] += value1;
      u_conv[index+2*N_UConv] += value2;
      // store global coordinates of the fe nodes (dof)
      x_conv[index] = X_orig[l];
      y_conv[index] = Y_orig[l];
      z_conv[index] = Z_orig[l];
      x_conv[index+N_UConv] = X_orig[l];
      y_conv[index+N_UConv] = Y_orig[l];
      z_conv[index+N_UConv] = Z_orig[l];
      x_conv[index+2*N_UConv] = X_orig[l];
      y_conv[index+2*N_UConv] = Y_orig[l];
      z_conv[index+2*N_UConv] = Z_orig[l];
      
      // starting to compute the integral in the neighbour cells
      // number of faces
      N_Faces=cell->GetN_Faces();
      // loop over all faces
      for(ll=0;ll<N_Faces;ll++)                           
      {
        // get pointer to edge[l]
        joint=cell->GetJoint(ll); 
        
        // if boundary edge continue 
        if (joint->GetType() == BoundaryFace)
          continue;
        
        // get pointer to neighbour cell
        neigh=cell->GetJoint(ll)->GetNeighbour(cell);
        
        // get id of neighbour
        neigh_i = neigh->GetClipBoard();
        
        // get id of finite element in neigh mesh cell wrt the
        // finite element space of the convolution
        CurrentElementNeigh = fespaceConv->GetFE3D(neigh_i, neigh);
        
        // get pointer to basis functions of the neighbour element
        BaseFunctNeigh = BaseFuncts[CurrentElementNeigh]; 
 
        // # of basis functions, is the same as N_loc_dof
        N_Neigh = N_BaseFunct[CurrentElementNeigh]; 
     
        // find the part of the global index array where the information
        // for the neigh mesh cell are stored: DOF
        DOFNeigh = GlobalNumbersConv + BeginIndexConv[neigh_i];
        
        // check if computation on neigh cell necessary
        // for all dof
        same_dof = 0;
        for (j=0;j<N_Neigh;j++)
        {
          //cout << "N " << DOFNeigh[j] << " ind " << index << endl;
          if (DOFNeigh[j]==index)
          same_dof = 1;
        }
        // if the cells have the same dof, continue
        if (same_dof==1)
          continue;
        
        ////////////////////////////////////////////////////////////////////
        
        // if not, do something! 
        // i.e. the same thing as before
        
        // get id of finite element in neigh mesh cell wrt the finite
        // element space of the velocity
        CurrentElementNeigh = fespace->GetFE3D(neigh_i, neigh);
        
        // get pointer to basis functions of the neighbour element
        BaseFunctNeigh = BaseFuncts[CurrentElementNeigh]; 
 
        // # of basis functions, is the same as N_loc_dof
        N_Neigh = N_BaseFunct[CurrentElementNeigh]; 
     
        // find the part of the global index array where the information
        // for the neigh mesh cell are stored: DOF
        DOFNeigh = GlobalNumbers + BeginIndex[neigh_i];

        // initialize auxiliary array UsedNeigh
        memset(UsedNeigh, 0, N_FEs3D*SizeOfInt);
        
        // get id of finite element space in neigh cell 
        CurrentElementNeigh = fespace->GetFE3D(neigh_i, neigh);
        // set entry in auxiliary array to 1
        UsedNeigh[CurrentElementNeigh] = 1;
       
        // have a look how much 1 are in UsedNeigh -> N_LocalUsedElementsNeigh
        N_LocalUsedElementsNeigh = 0;
        memset(LocalUsedElementsNeigh, 0, SizeOfInt*N_FEs3D);
        j = 0;
        for(k=0;k<N_FEs3D;k++)
          if(UsedNeigh[k])
          {
            LocalUsedElementsNeigh[j] = (FE3D)k;
            j++;
          }
        N_LocalUsedElementsNeigh = j;
        
        // ####################################################################
        // calculate values on original element
        // ####################################################################
        
        // From id of fe space the situation in the referenz mesh cell is 
        // known. Now, get information in the original mesh cell. 
        // Every finite element gets automatically a quadrature rule.
        // input: N_LocalUsedElementsNeigh, LocalUsedElementsNeigh, 
        //        neigh, SecondDer
        // output : N_PointsNeigh - number of quadrature points
        //          xiNeigh - xi-values of quad. points in ref. cell
        //          etaNeigh - eta-values of quad. points in ref. cell
        //          weigthsNeigh - weights fro quad. rule in quad. points 
        //          XNeigh - X-values of quad. points in original cell  
        //          YNeigh - Y-values of quad. points in original cell  
        //          AbsDetjkNeigh - absolute value of Jacobian in quad. points
        //          xiNeigh, etaNeigh, weightsNeigh, XNeigh, YNeigh,  
        //          AbsDetjkNeigh are pointers to arrays

        TFEDatabase3D::GetOrig(N_LocalUsedElementsNeigh, LocalUsedElementsNeigh, 
                               Coll, neigh, SecondDer,
                               N_PointsNeigh, xiNeigh, etaNeigh, zetaNeigh,
                               weightsNeigh,
                               XNeigh, YNeigh, ZNeigh, AbsDetjkNeigh);
    
        // copy the values of the finite element function from the global
        // array (Values) to local arrays 
        for(lll=0;lll<N_Neigh;lll++)
        {
          FEFunctValuesNeigh[lll] = Values[DOFNeigh[lll]];
          FEFunctValues1Neigh[lll] = Values[DOFNeigh[lll]+N_U];
          FEFunctValues2Neigh[lll] = Values[DOFNeigh[lll]+2*N_U];
        }
        
        // compute values of fe function in the quadrature points with
        // linear combinations of the local fe values and the basis functions
        for(k=0;k<N_Derivatives;k++)
        {
          // get fe values of the basis functions in the quad. points
          // in the current element (pointer to pointer)
          OrigFEValuesNeigh = TFEDatabase3D::GetOrigElementValues(BaseFunctNeigh,
                                                       NeededDerivatives[k]);
          // accumulate the values of the fe function 
          // do it for all quad. points
          for(j=0;j<N_PointsNeigh;j++)
          {
            // get pointer to array where the values for quad. point j are
            // stored
            Orig = OrigFEValuesNeigh[j];
            // initialize 
            value = value1 = value2 = 0;
            // compute linear combination over all local dof
            for(lll=0;lll<N_Neigh;lll++)
            {
              // value += u1[l] * phi[l]
              value += FEFunctValuesNeigh[lll] * Orig[lll];
              value1 += FEFunctValues1Neigh[lll] * Orig[lll];
              value2 += FEFunctValues2Neigh[lll] * Orig[lll];
            } // endfor lll
            // copy result to position in array 
            FEValue_u1_QuadPointNeigh[j] = value;
            FEValue_u2_QuadPointNeigh[j] = value1;
            FEValue_u3_QuadPointNeigh[j] = value2;
          } // endfor j
        } // endfor k
        
        // compute part of the convolution in all local fe nodes (dof)
        // which comes from current mesh cell 
        // loop over all local dof
        
        // initialize
        value = value1 = value2 = 0;
        // loop over all quadrature points 
        for(j=0;j<N_PointsNeigh;j++) 
        { 
          // square of the distance between local dof and quad. point
          // in original cell 
          distance_sq = ((X_orig[l]-XNeigh[j])*(X_orig[l]-XNeigh[j])+
                         (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j])+
                         (Z_orig[l]-ZNeigh[j])*(Z_orig[l]-ZNeigh[j])); 
          //cout << "X " << XNeigh[j] << " Y " << YNeigh[j] << endl;
          // compute value of GaussianFilter
          // input: delta, distance_sq
          g = GaussianFilter3D(delta,distance_sq);
          //cout << "dist " <<  sqrt(distance_sq) << endl;
          //
          // apply quadrature rule in ref. cell 
          value += g*AbsDetjkNeigh[j]*weightsNeigh[j]
                        *FEValue_u1_QuadPointNeigh[j];
          value1 += g*AbsDetjkNeigh[j]*weightsNeigh[j]
                        *FEValue_u2_QuadPointNeigh[j];
          value2 += g*AbsDetjkNeigh[j]*weightsNeigh[j]
                        *FEValue_u3_QuadPointNeigh[j];
            
        }
        // add local result to global array
        u_conv[index] += value;
        u_conv[index+N_UConv] += value1; 
        u_conv[index+2*N_UConv] += value2; 
      } // end of ll loop (edges)
    } // end of l loop (local dof)
  } // endfor i (cells)
  
 
  // for (i=0;i<N_DOFConv; i++)
  //{
  // cout << i<< " x " <<   x_conv[i] << " y  " << y_conv[i] << " z  " << z_conv[i] ;
  // cout << " u " << Values[i] << " u_conv " << u_conv[i] << endl;
  //}

  // release memory
  delete SecondDer;
  delete x_conv;
  delete y_conv;
  delete z_conv;
  TDatabase::ParamDB->INTERNAL_QUAD_RULE=0;
} // ConvoluteVelocity

void  ConvoluteVelocityFull3D(TFEVectFunct3D *u, TFEVectFunct3D *uConv)
{
  int i,j,k,l,n,m, N_UsedElements, N_LocalUsedElements;
  int N_LocalUsedElementsNeigh;
  int N_Cells, N_Points, N_Parameters, N_, N_U, N_DOF, N_loc_dofConv;
  int Used[N_FEs3D], UsedNeigh[N_FEs3D], *N_BaseFunct;
  int neigh_i, N_Neigh, ll, lll, N_PointsNeigh, N_Edges, same_dof;
  FE3D LocalUsedElements[N_FEs3D], CurrentElement, CurrentElementConv;
  FE3D CurrentElementNeigh,LocalUsedElementsNeigh[N_FEs3D] ;
  TFE3D *Element;
  TNodalFunctional3D *nf;
  BaseFunct3D BaseFunct, *BaseFuncts, BaseFunctNeigh;
  TCollection *Coll;
  TBaseCell *cell, *neigh;
  TFE3D *ele;
  double *weights, *xi, *eta, *zeta, *xi_ref, *eta_ref, *zeta_ref;
  double *weightsNeigh, *xiNeigh, *etaNeigh, *zetaNeigh;
  double X_orig[MaxN_PointsForNodal3D], Y_orig[MaxN_PointsForNodal3D];
  double Z_orig[MaxN_PointsForNodal3D];
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D], AbsDetjkConv[MaxN_QuadPoints_3D];
  double XNeigh[MaxN_QuadPoints_3D], YNeigh[MaxN_QuadPoints_3D];
  double ZNeigh[MaxN_QuadPoints_3D];
  double AbsDetjkNeigh[MaxN_QuadPoints_3D];
  RefTrans3D RefTrans;
  double *Param[MaxN_QuadPoints_3D], *aux;
  double FEValue_u1_QuadPoint[MaxN_QuadPoints_3D];
  double FEValue_u2_QuadPoint[MaxN_QuadPoints_3D];
  double FEValue_u3_QuadPoint[MaxN_QuadPoints_3D];
  double FEValue_u1_QuadPointNeigh[MaxN_QuadPoints_3D];
  double FEValue_u2_QuadPointNeigh[MaxN_QuadPoints_3D];
  double FEValue_u3_QuadPointNeigh[MaxN_QuadPoints_3D];
  double *ExactVal[MaxN_QuadPoints_3D];
  int *DOF, ActiveBound, DirichletBound, end, last, *DOFConv, *DOFNeigh;
  double **OrigFEValues, *Orig, value, value1, value2;
  double **OrigFEValuesNeigh;
  double FEFunctValues[MaxN_BaseFunctions3D];
  double FEFunctValues1[MaxN_BaseFunctions3D];
  double FEFunctValues2[MaxN_BaseFunctions3D];
  double FEFunctValuesNeigh[MaxN_BaseFunctions3D];
  double FEFunctValues1Neigh[MaxN_BaseFunctions3D];
  double FEFunctValues2Neigh[MaxN_BaseFunctions3D];
  int *GlobalNumbers, *BeginIndex, *GlobalNumbersConv, *BeginIndexConv;
  double LocError[4];
  double hK, delta;
  bool *SecondDer;
  int n_fespaces = 1, N_Derivatives = 1, index;
  int N_UConv, N_DOFConv, *conv_comp;
  double *Values, *u_conv,distance_sq,h_fine,g;
  double *u_values, *x_conv, *y_conv, *z_conv;
  TJoint *joint;
  TVertex *vertex0;
  double x_vertex, y_vertex, z_vertex;
  TFESpace3D *fespace, *fespaceConv; 

  //  TDatabase::ParamDB->INTERNAL_QUAD_RULE=1; there is something wrong !!!!
  // gives an array where the needed derivatives are described
  // defined in NavierStokes.h 
  // here we get as result {D000}
  MultiIndex3D NeededDerivatives[1] = { D000 };

  // get fe spaces of velocity u and of covoluted tensor
  fespace = u->GetFESpace3D();
  fespaceConv = uConv->GetFESpace3D();

  // gives a pointer to all available basis functions which
  // are described in the fedatabase
  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  // gives the number of available basis functions from the
  // fedatabase
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  // say that we don't need second derivatives
  // first: allocate memory for an array of length n_fespaces 
  // second: set all entries to FALSE
  SecondDer = new bool[n_fespaces];
  for(i=0;i<n_fespaces;i++)
    SecondDer[i] = FALSE;

  // get information of the numbering of the degrees of freedom
  // first: get pointer to array where a global numbering is stored
  GlobalNumbers = fespace->GetGlobalNumbers();
  // second: get pointer to array where the start of numbering for
  // each cell is stored
  BeginIndex = fespace->GetBeginIndex();
 
  // get number of dof for one velocity component N_U
  // -> total number of dof is twice N_U
  N_U = u->GetLength();
  N_DOF = 3*N_U;

  // get pointer to set of mesh cells which define the fe space
  Coll = fespace->GetCollection();
 
  // # dof in the new fe space
  N_UConv = fespaceConv->GetN_DegreesOfFreedom();
   N_DOFConv = 3*N_UConv;

  // allocate memory for values of convoluted function
  u_conv = uConv->GetValues(); 
  // initialize u_conv to 0
  memset(u_conv,0,N_DOFConv*SizeOfDouble);
  
  // allocate memory for values of convoluted function
  conv_comp = new int[N_UConv]; 
  // initialize conv_comp to 0
  memset(conv_comp,0,N_UConv*SizeOfInt);

  // get information of the numbering of the dof for the space of convolution
  // first: get pointer to array where a global numbering is stored
  GlobalNumbersConv = fespaceConv->GetGlobalNumbers();
  // second: get pointer to array where the start of numbering for
  // each cell is stored
  BeginIndexConv = fespaceConv->GetBeginIndex();

  // allocate memory for geometric position of nodes
  x_conv = new double[N_DOFConv]; 
  y_conv = new double[N_DOFConv]; 
  z_conv = new double[N_DOFConv]; 
   
// ########################################################################
// computing the convolution
// loop over all cells
// ########################################################################

  // get pointer to set of mesh cells which define the fe space
  Coll = fespace->GetCollection(); 
  // get number of mesh cells
  N_Cells = Coll->GetN_Cells();

  // direct a pointer to the values of the finite element vector
  // function
  // first part of the array : u1
  // second part of the array : u2
  Values = u->GetValues();           

  // store positition of each mesh cells in the collection of mesh 
  // cells: set ClipBoard to i
  for(i=0;i<N_Cells;i++)                       
  {
    cell = Coll->GetCell(i);
    cell->SetClipBoard(i);
  }

  // loop over all mesh cells to accumulate the values for the 
  // convolution 
  for(i=0;i<N_Cells;i++)
  {
    // get current mesh cell
    cell = Coll->GetCell(i);        
    // get cell diameter (longest distance between to vertices of the cell)
    hK = cell->GetDiameter();       
    delta =  CharacteristicFilterWidth(hK);
 
    // ####################################################################
    // find local used finite elements on this cell
    // ####################################################################
    
    // initialize auxiliary array Used
    memset(Used, 0, N_FEs3D*SizeOfInt);
 
    for(j=0;j<n_fespaces;j++)
    {
      // get id of finite element space in this cell 
      CurrentElement = fespace->GetFE3D(i, cell);
      // set entry in auxiliary array to 1
      Used[CurrentElement] = 1;
    }

    // have a look how much 1 are in Used ->   N_LocalUsedElements
    N_LocalUsedElements = 0;
    memset(LocalUsedElements, 0, SizeOfInt*N_FEs3D);
    j = 0;
    for(k=0;k<N_FEs3D;k++)
      if(Used[k])
      {
        LocalUsedElements[j] = (FE3D)k;
        j++;
      }
    N_LocalUsedElements = j;
    
    // ####################################################################
    // calculate values on original element
    // ####################################################################
    
    // From id of fe space the situation in the referenz mesh cell is 
    // known. Now, get information in the original mesh cell. 
    // Every finite element gets automatically a quadrature rule.
    // input: N_LocalUsedElements, LocalUsedElements, cell, SecondDer
    // output : N_Points - number of quadrature points
    //          xi - xi-values of quad. points in ref. cell
    //          eta - eta-values of quad. points in ref. cell
    //          weigths - weights fro quad. rule in quad. points 
    //          X - X-values of quad. points in original cell  
    //          Y - Y-values of quad. points in original cell  
    //          AbsDetjk - absolute value of Jacobian in quad. points
    //          xi, eta, weights, X, Y,  AbsDetjk are pointers to arrays

    TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                           Coll, cell, SecondDer, N_Points, xi, eta, zeta,
                           weights, X, Y, Z, AbsDetjk);
    
    // compute geometric positions of the fe nodes
    // get id of finite element in current mesh cell
    CurrentElement = fespace->GetFE3D(i, cell);
    CurrentElementConv = fespaceConv->GetFE3D(i, cell);
    // get fe from its id
    Element = TFEDatabase3D::GetFE3D(CurrentElementConv);
    // get reference transformation 
    RefTrans = Element->GetRefTransID();
    // get pointer to the nodal functionals (fe nodes) of the fe 
    // (in ref. cell)
    nf = Element->GetNodalFunctional3D();
    // get number and coordinates of local dof in ref cell
    // xi_ref, eta_ref are pointers
    nf->GetPointsForAll(N_loc_dofConv, xi_ref, eta_ref, zeta_ref);
      
    // get coordinates of fe nodes in original cell
    // input: RefTrans, N_loc_dof, xi_ref, eta_ref
    // output : X_orig, Y_orig - pointers to arrays with coordinates
    //          AbsDetjk - same as above    
    TFEDatabase3D::GetOrigFromRef(RefTrans,N_loc_dofConv, 
                                  xi_ref, eta_ref, zeta_ref,
                                  X_orig, Y_orig, Z_orig, AbsDetjkConv);   
   
    // get pointer to basis functions of the current element
    BaseFunct = BaseFuncts[CurrentElement];  
    // # of basis functions, is the same as N_loc_dof
    N_ = N_BaseFunct[CurrentElement];        

    // find the part of the global index array where the information
    // for the current mesh cell are stored: DOF
    DOF = GlobalNumbers + BeginIndex[i];
    DOFConv = GlobalNumbersConv + BeginIndexConv[i];
    // copy the values of the finite element function from the global
    // array (Values) to local arrays 
    for(l=0;l<N_;l++)
    {
      FEFunctValues[l] = Values[DOF[l]];
      FEFunctValues1[l] = Values[DOF[l]+N_U];
      FEFunctValues2[l] = Values[DOF[l]+2*N_U];
    }

    // compute values of fe function in the quadrature points with
    // linear combinations of the local fe values and the basis functions
    for(k=0;k<N_Derivatives;k++)
    {
      // get fe values of the basis functions in the quad. points
      // in the current element (pointer to pointer)
      OrigFEValues = TFEDatabase3D::GetOrigElementValues(BaseFunct,
                                            NeededDerivatives[k]);
      // accumulate the values of the fe function 
      // do it for all quad. points
      for(j=0;j<N_Points;j++)
      {
        // get pointer to array where the values for quad. point j are
        // stored
        Orig = OrigFEValues[j];
        // initialize 
        value = value1 = value2 =  0;
        // compute linear combination over all local dof
        for(l=0;l<N_;l++)
        {
          // value += u1[l] * phi[l]
          value += FEFunctValues[l] * Orig[l];
          value1 += FEFunctValues1[l] * Orig[l];
          value2 += FEFunctValues2[l] * Orig[l];
        } // endfor l
        // copy result to position in array 
        FEValue_u1_QuadPoint[j] = value;
        FEValue_u2_QuadPoint[j] = value1;
        FEValue_u3_QuadPoint[j] = value2;
      } // endfor j
    } // endfor k
    
    // compute part of the convolution in all local fe nodes (dof)
    // which comes from current mesh cell 
    // loop over all local dof
    for(l=0;l<N_loc_dofConv;l++)         
    {
      // get global index of local dof
      index = DOFConv[l];          

      // check if value is already computed
      // if yes : continue
      if (conv_comp[index])
        continue;

      // set conv_comp[index]
      conv_comp[index] = 1;
      // initialize
      value = value1 = value2 = 0;
      // loop over all quadrature points 
      for(j=0;j<N_Points;j++) 
      {
        // square of the distance between local dof and quad. point
        // in original cell 
        distance_sq = ((X_orig[l]-X[j])*(X_orig[l]-X[j])+
                       (Y_orig[l]-Y[j])*(Y_orig[l]-Y[j])+
                       (Z_orig[l]-Z[j])*(Z_orig[l]-Z[j])); 
        // compute value of GaussianFilter
        // input: delta, distance_sq
        g = GaussianFilter3D(delta,distance_sq);  
        // apply quadrature rule in ref. cell 
        value += g*AbsDetjk[j]*weights[j]*FEValue_u1_QuadPoint[j];
        value1 += g*AbsDetjk[j]*weights[j]*FEValue_u2_QuadPoint[j];
        value2 += g*AbsDetjk[j]*weights[j]*FEValue_u3_QuadPoint[j];
      }
      // add local result to global array
      u_conv[index] += value;
      u_conv[index+N_UConv] += value1;
      u_conv[index+2*N_UConv] += value2;
      // store global coordinates of the fe nodes (dof)
      x_conv[index] = X_orig[l];
      y_conv[index] = Y_orig[l];
      z_conv[index] = Z_orig[l];
      x_conv[index+N_UConv] = X_orig[l];
      y_conv[index+N_UConv] = Y_orig[l];
      z_conv[index+N_UConv] = Z_orig[l];
      x_conv[index+2*N_UConv] = X_orig[l];
      y_conv[index+2*N_UConv] = Y_orig[l];
      z_conv[index+2*N_UConv] = Z_orig[l];
      
      // starting to compute the integral in ALL the cells      
      // loop over all cells
      for(ll=0;ll<N_Cells;ll++)                           
      {
        // check if we are at the current mesh cell
        if (ll==i)
          continue;
        
        neigh = Coll->GetCell(ll); 
        // get one vertex of the neigh 
        vertex0 = neigh->GetVertex(0);
        // get coordinates of vertex
        x_vertex = vertex0->GetX();
        y_vertex = vertex0->GetY();
        z_vertex = vertex0->GetZ();
        // compute distance to current point where the convolution
        // is computed
        distance_sq = ((X_orig[l]-x_vertex)*(X_orig[l]-x_vertex)+
                       (Y_orig[l]-y_vertex)*(Y_orig[l]-y_vertex)+
                       (Z_orig[l]-z_vertex)*(Z_orig[l]-z_vertex));
        // check if neigh is close enough
        if (distance_sq > 2*delta *delta)
          continue;
        // get id of neighbour
        neigh_i = neigh->GetClipBoard();
        
        // get id of finite element in neigh mesh cell wrt the finite
        // element space of the velocity
        CurrentElementNeigh = fespace->GetFE3D(neigh_i, neigh);
        
        // get pointer to basis functions of the neighbour element
        BaseFunctNeigh = BaseFuncts[CurrentElementNeigh]; 
 
        // # of basis functions, is the same as N_loc_dof
        N_Neigh = N_BaseFunct[CurrentElementNeigh]; 
     
        // find the part of the global index array where the information
        // for the neigh mesh cell are stored: DOF
        DOFNeigh = GlobalNumbers + BeginIndex[neigh_i];

        // initialize auxiliary array UsedNeigh
        memset(UsedNeigh, 0, N_FEs3D*SizeOfInt);
        
        // set entry in auxiliary array to 1
        UsedNeigh[CurrentElementNeigh] = 1;
       
        // have a look how much 1 are in UsedNeigh -> N_LocalUsedElementsNeigh
        N_LocalUsedElementsNeigh = 0;
        memset(LocalUsedElementsNeigh, 0, SizeOfInt*N_FEs3D);
        j = 0;
        for(k=0;k<N_FEs3D;k++)
          if(UsedNeigh[k])
          {
            LocalUsedElementsNeigh[j] = (FE3D)k;
            j++;
          }
        N_LocalUsedElementsNeigh = j;
        
        // ####################################################################
        // calculate values on the original element
        // ####################################################################
        
        // From id of fe space the situation in the referenz mesh cell is 
        // known. Now, get information in the original mesh cell. 
        // Every finite element gets automatically a quadrature rule.
        // input: N_LocalUsedElementsNeigh, LocalUsedElementsNeigh, 
        //        neigh, SecondDer
        // output : N_PointsNeigh - number of quadrature points
        //          xiNeigh - xi-values of quad. points in ref. cell
        //          etaNeigh - eta-values of quad. points in ref. cell
        //          weigthsNeigh - weights fro quad. rule in quad. points 
        //          XNeigh - X-values of quad. points in original cell  
        //          YNeigh - Y-values of quad. points in original cell  
        //          AbsDetjkNeigh - absolute value of Jacobian in quad. points
        //          xiNeigh, etaNeigh, weightsNeigh, XNeigh, YNeigh,  
        //          AbsDetjkNeigh are pointers to arrays

        TFEDatabase3D::GetOrig(N_LocalUsedElementsNeigh, LocalUsedElementsNeigh, 
                               Coll, neigh, SecondDer,
                               N_PointsNeigh, xiNeigh, etaNeigh, zetaNeigh,
                               weightsNeigh,
                               XNeigh, YNeigh, ZNeigh, AbsDetjkNeigh);
    
        // copy the values of the finite element function from the global
        // array (Values) to local arrays 
        for(lll=0;lll<N_Neigh;lll++)
        {
          FEFunctValuesNeigh[lll] = Values[DOFNeigh[lll]];
          FEFunctValues1Neigh[lll] = Values[DOFNeigh[lll]+N_U];
          FEFunctValues2Neigh[lll] = Values[DOFNeigh[lll]+2*N_U];
        }
        
        // compute values of fe function in the quadrature points with
        // linear combinations of the local fe values and the basis functions
        for(k=0;k<N_Derivatives;k++)
        {
          // get fe values of the basis functions in the quad. points
          // in the current element (pointer to pointer)
          OrigFEValuesNeigh = TFEDatabase3D::GetOrigElementValues(BaseFunctNeigh,
                                                                NeededDerivatives[k]);
          // accumulate the values of the fe function 
          // do it for all quad. points
          for(j=0;j<N_PointsNeigh;j++)
          {
            // get pointer to array where the values for quad. point j are
            // stored
            Orig = OrigFEValuesNeigh[j];
            // initialize 
            value = value1 = value2 = 0;
            // compute linear combination over all local dof
            for(lll=0;lll<N_Neigh;lll++)
            {
              // value += u1[l] * phi[l]
              value += FEFunctValuesNeigh[lll] * Orig[lll];
              value1 += FEFunctValues1Neigh[lll] * Orig[lll];
              value2 += FEFunctValues2Neigh[lll] * Orig[lll];
            } // endfor lll
            // copy result to position in array 
            FEValue_u1_QuadPointNeigh[j] = value;
            FEValue_u2_QuadPointNeigh[j] = value1;
            FEValue_u3_QuadPointNeigh[j] = value2;
          } // endfor j
        } // endfor k
        
        // compute part of the convolution in all local fe nodes (dof)
        // which comes from current mesh cell 
        // loop over all local dof
        
        // initialize
        value = value1 = value2 = 0;
        // loop over all quadrature points 
        for(j=0;j<N_PointsNeigh;j++) 
        { 
          // square of the distance between local dof and quad. point
          // in original cell 
          distance_sq = ((X_orig[l]-XNeigh[j])*(X_orig[l]-XNeigh[j])+
                         (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j])+
                         (Z_orig[l]-ZNeigh[j])*(Z_orig[l]-ZNeigh[j])); 
          //cout << "X " << XNeigh[j] << " Y " << YNeigh[j] << endl;
          // compute value of GaussianFilter
          // input: delta, distance_sq
          g = GaussianFilter3D(delta,distance_sq);
          //cout << "dist " <<  sqrt(distance_sq) << endl;
          
          // apply quadrature rule in ref. cell 
          value += g*AbsDetjkNeigh[j]*weightsNeigh[j]*FEValue_u1_QuadPointNeigh[j];
          value1 += g*AbsDetjkNeigh[j]*weightsNeigh[j]*FEValue_u2_QuadPointNeigh[j];
          value2 += g*AbsDetjkNeigh[j]*weightsNeigh[j]*FEValue_u3_QuadPointNeigh[j];       
        }
        // add local `result to global array
        u_conv[index] += value;
        u_conv[index+N_UConv] += value1; 
        u_conv[index+2*N_UConv] += value2; 
      } // end of ll loop (edges)
    } // end of l loop (local dof)
  } // endfor i (cells)
  
 
/*  for (i=0;i<N_DOFConv; i++)
  {
    cout << i<< " x " <<   x_conv[i] << " y  " << y_conv[i] << " z  " << z_conv[i] ;
    cout << " u " << Values[i] << " u_conv_full " << u_conv[i] << endl;
  }
*/

  // release memory
  delete SecondDer;
  delete conv_comp;
  delete x_conv;
  delete y_conv;
  delete z_conv;
  TDatabase::ParamDB->INTERNAL_QUAD_RULE=0;
} // ConvoluteVelocityFull

// ========================================================================
// convolute D = (grad w grad w^T) by local integration
// use D_{21} = D_{12}
// input u
// output DuTensor
// ========================================================================
void  ConvoluteSymmetricTensor3D(TFEVectFunct3D *u, TFEVectFunct3D *duTensor)

{
  int i,j,k,l,n,m, N_UsedElements, N_LocalUsedElements;
  int N_LocalUsedElementsNeigh;
  int N_Cells, N_Points, N_Parameters, N_, N_U, N_DOF, N_loc_dofConv;
  int Used[N_FEs3D], UsedNeigh[N_FEs3D], *N_BaseFunct;
  int neigh_i, N_Neigh, ll, lll, N_PointsNeigh, N_Faces, same_dof;
  FE3D LocalUsedElements[N_FEs3D], CurrentElement, CurrentElementConv;
  FE3D CurrentElementNeigh,LocalUsedElementsNeigh[N_FEs3D] ;
  TFE3D *Element;
  TNodalFunctional3D *nf;
  BaseFunct3D BaseFunct, *BaseFuncts, BaseFunctNeigh;
  TCollection *Coll;
  TBaseCell *cell, *neigh;
  TFE3D *ele;
  double *weights, *xi, *eta, *zeta, *xi_ref, *eta_ref, *zeta_ref;
  double *weightsNeigh, *xiNeigh, *etaNeigh, *zetaNeigh;
  double X_orig[MaxN_PointsForNodal3D], Y_orig[MaxN_PointsForNodal3D];
  double Z_orig[MaxN_PointsForNodal3D];
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D], AbsDetjkConv[MaxN_QuadPoints_3D];
  double XNeigh[MaxN_QuadPoints_3D], YNeigh[MaxN_QuadPoints_3D];
  double ZNeigh[MaxN_QuadPoints_3D];
  double AbsDetjkNeigh[MaxN_QuadPoints_3D];
  RefTrans3D RefTrans;
  double *Param[MaxN_QuadPoints_3D], *aux;
  double FEValue_D1u1_QuadPoint[MaxN_QuadPoints_3D];
  double FEValue_D1u2_QuadPoint[MaxN_QuadPoints_3D];
  double FEValue_D1u3_QuadPoint[MaxN_QuadPoints_3D];
  double FEValue_D2u1_QuadPoint[MaxN_QuadPoints_3D];
  double FEValue_D2u2_QuadPoint[MaxN_QuadPoints_3D];
  double FEValue_D2u3_QuadPoint[MaxN_QuadPoints_3D];
  double FEValue_D3u1_QuadPoint[MaxN_QuadPoints_3D];
  double FEValue_D3u2_QuadPoint[MaxN_QuadPoints_3D];
  double FEValue_D3u3_QuadPoint[MaxN_QuadPoints_3D];
  double FEValue_D1u1_QuadPointNeigh[MaxN_QuadPoints_3D];
  double FEValue_D1u2_QuadPointNeigh[MaxN_QuadPoints_3D];
  double FEValue_D1u3_QuadPointNeigh[MaxN_QuadPoints_3D];
  double FEValue_D2u1_QuadPointNeigh[MaxN_QuadPoints_3D];
  double FEValue_D2u2_QuadPointNeigh[MaxN_QuadPoints_3D];
  double FEValue_D2u3_QuadPointNeigh[MaxN_QuadPoints_3D];
  double FEValue_D3u1_QuadPointNeigh[MaxN_QuadPoints_3D];
  double FEValue_D3u2_QuadPointNeigh[MaxN_QuadPoints_3D];
  double FEValue_D3u3_QuadPointNeigh[MaxN_QuadPoints_3D];
  double *ExactVal[MaxN_QuadPoints_3D];
  int *DOF, ActiveBound, DirichletBound, end, last, *DOFConv, *DOFNeigh;
  double **OrigFEValues, *Orig;
  double value, value1, value2;
  double value11, value12, value13, value22, value23, value33;
  double **OrigFEValuesNeigh;
  double FEFunctValues[MaxN_BaseFunctions3D];
  double FEFunctValues1[MaxN_BaseFunctions3D];
  double FEFunctValues2[MaxN_BaseFunctions3D];
  double FEFunctValuesNeigh[MaxN_BaseFunctions3D];
  double FEFunctValues1Neigh[MaxN_BaseFunctions3D];
  double FEFunctValues2Neigh[MaxN_BaseFunctions3D];
  int *GlobalNumbers, *BeginIndex, *GlobalNumbersConv, *BeginIndexConv;
  double LocError[4];
  int N_UConv, N_DOFConv, index;
  double *Values, *u_conv,distance_sq,h_fine,g, *u_values;
  double *x_conv, *y_conv, *z_conv, delta, hK;
  TJoint *joint;
  TFESpace3D *fespace, *fespaceConv; 

  bool SecondDer[1] = { FALSE };
  int N_Derivatives = 3;
  MultiIndex3D NeededDerivatives[3] = { D100, D010, D001 };
  //TDatabase::ParamDB->INTERNAL_QUAD_RULE=1; there is something wrong !!!!

  // get fe spaces of velocity u and of covoluted tensor
  fespace = u->GetFESpace3D();
  fespaceConv = duTensor->GetFESpace3D();

  // gives a pointer to all available basis functions which
  // are described in the fedatabase
  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  // gives the number of available basis functions from the
  // fedatabase
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  // get information of the numbering of the degrees of freedom
  // first: get pointer to array where a global numbering is stored
  GlobalNumbers = fespace->GetGlobalNumbers();
  // second: get pointer to array where the start of numbering for
  // each cell is stored
  BeginIndex = fespace->GetBeginIndex();
 
  // get number of dof for one velocity component N_U
  // -> total number of dof is twice N_U
  N_U = u->GetLength();
  N_DOF = 3*N_U;
  
  // get pointer to set of mesh cells which define the fe space
  Coll = fespace->GetCollection();
  // get number of mesh cells
  N_Cells = Coll->GetN_Cells();
 
  // # dof in the new fe space
  N_UConv = fespaceConv->GetN_DegreesOfFreedom();
 
  N_DOFConv = 6*N_UConv;
  // allocate memory for values of convoluted function
  u_conv = duTensor->GetValues(); 
  // initialize u_conv to 0
  memset(u_conv,0,N_DOFConv*SizeOfDouble);

  // get information of the numbering of the dof for the space of convolution
  // first: get pointer to array where a global numbering is stored
  GlobalNumbersConv = fespaceConv->GetGlobalNumbers();
  // second: get pointer to array where the start of numbering for
  // each cell is stored
  BeginIndexConv = fespaceConv->GetBeginIndex();

  // allocate memory for geometric position of nodes
  x_conv = new double[N_DOFConv]; 
  y_conv = new double[N_DOFConv]; 
  z_conv = new double[N_DOFConv]; 
   
// ########################################################################
// computing the convolution
// loop over all cells
// ########################################################################

  // direct a pointer to the values of the finite element vector
  // function
  // first part of the array : u1
  // second part of the array : u2
  Values = u->GetValues();           

  // store positition of each mesh cells in the collection of mesh 
  // cells: set ClipBoard to i
  for(i=0;i<N_Cells;i++)                       
  {
    cell = Coll->GetCell(i);
    cell->SetClipBoard(i);
  }
 
  // loop over all mesh cells to accumulate the values for the 
  // convolution 
  for(i=0;i<N_Cells;i++)
  {
    // get current mesh cell
    cell = Coll->GetCell(i);        
    hK = cell->GetDiameter();
    delta =  CharacteristicFilterWidth(hK);

    // ####################################################################
    // find local used finite elements on this cell
    // ####################################################################
    CurrentElement = fespace->GetFE3D(i, cell);
    LocalUsedElements[0] = CurrentElement;
    N_LocalUsedElements = 1;
    
    // ####################################################################
    // calculate values on original element
    // ####################################################################
    
    // From id of fe space the situation in the referenz mesh cell is 
    // known. Now, get information in the original mesh cell. 
    // Every finite element gets automatically a quadrature rule.
    // input: N_LocalUsedElements, LocalUsedElements, cell, SecondDer
    // output : N_Points - number of quadrature points
    //          xi - xi-values of quad. points in ref. cell
    //          eta - eta-values of quad. points in ref. cell
    //          weigths - weights fro quad. rule in quad. points 
    //          X - X-values of quad. points in original cell  
    //          Y - Y-values of quad. points in original cell  
    //          AbsDetjk - absolute value of Jacobian in quad. points
    //          xi, eta, weights, X, Y,  AbsDetjk are pointers to arrays

    TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                           Coll, cell, SecondDer, N_Points, 
                           xi, eta, zeta,
                           weights, X, Y, Z, AbsDetjk);
    //cout << N_Points << " ";
    // compute geometric positions of the fe nodes
    // get id of finite element in current mesh cell
    CurrentElement = fespace->GetFE3D(i, cell);
    CurrentElementConv = fespaceConv->GetFE3D(i, cell);
    // get fe from its id
    Element = TFEDatabase3D::GetFE3D(CurrentElementConv);
    // get reference transformation 
    RefTrans = Element->GetRefTransID();
    // get pointer to the nodal functionals (fe nodes) of the fe 
    // (in ref. cell)
    nf = Element->GetNodalFunctional3D();
    // get number and coordinates of local dof in ref cell
    // xi_ref, eta_ref are pointers
    nf->GetPointsForAll(N_loc_dofConv, xi_ref, eta_ref, zeta_ref);
      
    // get coordinates of fe nodes in original cell
    // input: RefTrans, N_loc_dof, xi_ref, eta_ref
    // output : X_orig, Y_orig - pointers to arrays with coordinates
    //          AbsDetjk - same as above    
    TFEDatabase3D::GetOrigFromRef(RefTrans,N_loc_dofConv, 
                                  xi_ref, eta_ref, zeta_ref,
                                  X_orig, Y_orig, Z_orig, AbsDetjkConv);   
   
    // get pointer to basis functions of the current element
    BaseFunct = BaseFuncts[CurrentElement];  
    // # of basis functions, is the same as N_loc_dof
    N_ = N_BaseFunct[CurrentElement];        

    // find the part of the global index array where the information
    // for the current mesh cell are stored: DOF
    DOF = GlobalNumbers + BeginIndex[i];
    DOFConv = GlobalNumbersConv + BeginIndexConv[i];
    // copy the values of the finite element function from the global
    // array (Values) to local arrays 
    for(l=0;l<N_;l++)
    {
      FEFunctValues[l] = Values[DOF[l]];
      FEFunctValues1[l] = Values[DOF[l]+N_U];
      FEFunctValues2[l] = Values[DOF[l]+2*N_U];
    }

    // compute values of fe function in the quadrature points with
    // linear combinations of the local fe values and the basis functions
    for(k=0;k<N_Derivatives;k++)
    {
      // get fe values of the basis functions in the quad. points
      // in the current element (pointer to pointer)
      OrigFEValues = TFEDatabase3D::GetOrigElementValues(BaseFunct,
                                            NeededDerivatives[k]);
      // accumulate the values of the fe function 
      // do it for all quad. points
      for(j=0;j<N_Points;j++)
      {
        // get pointer to array where the values for quad. point j are
        // stored
        Orig = OrigFEValues[j];
        // initialize 
        value = value1 = value2 = 0;
        // compute linear combination over all local dof
        for(l=0;l<N_;l++)
        {
          // value += u1[l] * phi[l]
          value += FEFunctValues[l] * Orig[l];
          value1 += FEFunctValues1[l] * Orig[l];
          value2 += FEFunctValues2[l] * Orig[l];
        } // endfor l
        // copy result to position in array 

        switch(k)
        {
          case 0 :
          FEValue_D1u1_QuadPoint[j] = value;
          FEValue_D1u2_QuadPoint[j] = value1;
          FEValue_D1u3_QuadPoint[j] = value2;
          break;
          case 1 :
          FEValue_D2u1_QuadPoint[j] = value;
          FEValue_D2u2_QuadPoint[j] = value1;
          FEValue_D2u3_QuadPoint[j] = value2;
          break;
          case 2 :
          FEValue_D3u1_QuadPoint[j] = value;
          FEValue_D3u2_QuadPoint[j] = value1;
          FEValue_D3u3_QuadPoint[j] = value2;
          break;
        }
      } // endfor j
    } // endfor k
    
    // compute part of the convolution in all local fe nodes (dof)
    // which comes from current mesh cell 
    // loop over all local dof
    for(l=0;l<N_loc_dofConv;l++)         
    {
      // get global indes of local dof
      index = DOFConv[l];          
      // initialize
      value11 = value12 = value13 = value22 = value23 = value33 = 0;
      // loop over all quadrature points 
      for(j=0;j<N_Points;j++) 
      {
        // square of the distance between local dof and quad. point
        // in original cell 
        distance_sq = ((X_orig[l]-X[j])*(X_orig[l]-X[j])+
                       (Y_orig[l]-Y[j])*(Y_orig[l]-Y[j])+
                       (Z_orig[l]-Z[j])*(Z_orig[l]-Z[j])); 
        // compute value of GaussianFilter
        // input: delta, distance_sq
        g = GaussianFilter3D(delta,distance_sq);  
        // apply quadrature rule in ref. cell 
        value11 += g*AbsDetjk[j]*weights[j]*
                   ( FEValue_D1u1_QuadPoint[j]*FEValue_D1u1_QuadPoint[j]
                    +FEValue_D2u1_QuadPoint[j]*FEValue_D2u1_QuadPoint[j]
                    +FEValue_D3u1_QuadPoint[j]*FEValue_D3u1_QuadPoint[j]);
        value12 += g*AbsDetjk[j]*weights[j]*
                   ( FEValue_D1u1_QuadPoint[j]*FEValue_D1u2_QuadPoint[j]
                     +FEValue_D2u1_QuadPoint[j]*FEValue_D2u2_QuadPoint[j]
                     +FEValue_D3u1_QuadPoint[j]*FEValue_D3u2_QuadPoint[j]);
        value13 += g*AbsDetjk[j]*weights[j]*
                   ( FEValue_D1u1_QuadPoint[j]*FEValue_D1u3_QuadPoint[j]
                     +FEValue_D2u1_QuadPoint[j]*FEValue_D2u3_QuadPoint[j]
                     +FEValue_D3u1_QuadPoint[j]*FEValue_D3u3_QuadPoint[j]);
        value22 += g*AbsDetjk[j]*weights[j]*
                   ( FEValue_D1u2_QuadPoint[j]*FEValue_D1u2_QuadPoint[j]
                     +FEValue_D2u2_QuadPoint[j]*FEValue_D2u2_QuadPoint[j]
                     +FEValue_D3u2_QuadPoint[j]*FEValue_D3u2_QuadPoint[j]);
        value23 += g*AbsDetjk[j]*weights[j]*
                   ( FEValue_D1u2_QuadPoint[j]*FEValue_D1u3_QuadPoint[j]
                     +FEValue_D2u2_QuadPoint[j]*FEValue_D2u3_QuadPoint[j]
                     +FEValue_D3u2_QuadPoint[j]*FEValue_D3u3_QuadPoint[j]);
        value33 += g*AbsDetjk[j]*weights[j]*
                   ( FEValue_D1u3_QuadPoint[j]*FEValue_D1u3_QuadPoint[j]
                     +FEValue_D2u3_QuadPoint[j]*FEValue_D2u3_QuadPoint[j]
                     +FEValue_D3u3_QuadPoint[j]*FEValue_D3u3_QuadPoint[j]);
      }
      // add local result to global array
      u_conv[index] += value11;
      u_conv[index+N_UConv] += value12;
      u_conv[index+2*N_UConv] += value13;
      u_conv[index+3*N_UConv] += value22;
      u_conv[index+4*N_UConv] += value23;
      u_conv[index+5*N_UConv] += value33;
      // store global coordinates of the fe nodes (dof)
      x_conv[index] = X_orig[l];
      y_conv[index] = Y_orig[l];
      z_conv[index] = Z_orig[l];
      x_conv[index+N_UConv] = X_orig[l];
      y_conv[index+N_UConv] = Y_orig[l];
      z_conv[index+N_UConv] = Z_orig[l];
      x_conv[index+2*N_UConv] = X_orig[l];
      y_conv[index+2*N_UConv] = Y_orig[l];
      z_conv[index+2*N_UConv] = Z_orig[l];
      
      // starting to compute the integral in the neighbour cells
      // number of edges
      N_Faces=cell->GetN_Faces();
      // loop over all edges
      for(ll=0;ll<N_Faces;ll++)                           
      {
        // get pointer to edge[l]
        joint=cell->GetJoint(ll); 
        
        // if boundary edge continue 
        if (joint->GetType() == BoundaryFace)
          continue;
        
        // get pointer to neighbour cell
        neigh=cell->GetJoint(ll)->GetNeighbour(cell);
        
        // get id of neighbour
        neigh_i = neigh->GetClipBoard();
        
        // get id of finite element in neigh mesh cell wrt the
        // finite element space of the convolution
        CurrentElementNeigh = fespaceConv->GetFE3D(neigh_i, neigh);
        
        // get pointer to basis functions of the neighbour element
        BaseFunctNeigh = BaseFuncts[CurrentElementNeigh]; 
 
        // # of basis functions, is the same as N_loc_dof
        N_Neigh = N_BaseFunct[CurrentElementNeigh]; 
     
        // find the part of the global index array where the information
        // for the neigh mesh cell are stored: DOF
        DOFNeigh = GlobalNumbersConv + BeginIndexConv[neigh_i];
        
        // check if computation on neigh cell necessary
        // for all dof
        same_dof = 0;
        for (j=0;j<N_Neigh;j++)
        {
          //cout << "N " << DOFNeigh[j] << " ind " << index << endl;
          if (DOFNeigh[j]==index)
          same_dof = 1;
        }
        // if the cells have the same dof, continue
        if (same_dof==1)
          continue;
        
        ////////////////////////////////////////////////////////////////////
        
        // if not, do something! 
        // i.e. the same thing as before
        
        // get id of finite element in neigh mesh cell wrt the finite
        // element space of the velocity
        CurrentElementNeigh = fespace->GetFE3D(neigh_i, neigh);
        
        // get pointer to basis functions of the neighbour element
        BaseFunctNeigh = BaseFuncts[CurrentElementNeigh]; 
 
        // # of basis functions, is the same as N_loc_dof
        N_Neigh = N_BaseFunct[CurrentElementNeigh]; 
     
        // find the part of the global index array where the information
        // for the neigh mesh cell are stored: DOF
        DOFNeigh = GlobalNumbers + BeginIndex[neigh_i];

        // initialize auxiliary array UsedNeigh
        memset(UsedNeigh, 0, N_FEs3D*SizeOfInt);
        
        // get id of finite element space in neigh cell 
        CurrentElementNeigh = fespace->GetFE3D(neigh_i, neigh);
        // set entry in auxiliary array to 1
        UsedNeigh[CurrentElementNeigh] = 1;
       
        // have a look how much 1 are in UsedNeigh -> N_LocalUsedElementsNeigh
        N_LocalUsedElementsNeigh = 1;
        LocalUsedElementsNeigh[0] = CurrentElementNeigh;
        
        // ####################################################################
        // calculate values on original element
        // ####################################################################
        
        // From id of fe space the situation in the referenz mesh cell is 
        // known. Now, get information in the original mesh cell. 
        // Every finite element gets automatically a quadrature rule.
        // input: N_LocalUsedElementsNeigh, LocalUsedElementsNeigh, 
        //        neigh, SecondDer
        // output : N_PointsNeigh - number of quadrature points
        //          xiNeigh - xi-values of quad. points in ref. cell
        //          etaNeigh - eta-values of quad. points in ref. cell
        //          weigthsNeigh - weights fro quad. rule in quad. points 
        //          XNeigh - X-values of quad. points in original cell  
        //          YNeigh - Y-values of quad. points in original cell  
        //          AbsDetjkNeigh - absolute value of Jacobian in quad. points
        //          xiNeigh, etaNeigh, weightsNeigh, XNeigh, YNeigh,  
        //          AbsDetjkNeigh are pointers to arrays

        TFEDatabase3D::GetOrig(N_LocalUsedElementsNeigh, LocalUsedElementsNeigh, 
                               Coll, neigh, SecondDer, N_PointsNeigh, 
                               xiNeigh, etaNeigh, zetaNeigh, weightsNeigh,
                               XNeigh, YNeigh, ZNeigh, AbsDetjkNeigh);
    
        // copy the values of the finite element function from the global
        // array (Values) to local arrays 
        for(lll=0;lll<N_Neigh;lll++)
        {
          FEFunctValuesNeigh[lll] = Values[DOFNeigh[lll]];
          FEFunctValues1Neigh[lll] = Values[DOFNeigh[lll]+N_U];
          FEFunctValues2Neigh[lll] = Values[DOFNeigh[lll]+2*N_U];
        }
        
        // compute values of fe function in the quadrature points with
        // linear combinations of the local fe values and the basis functions
        for(k=0;k<N_Derivatives;k++)
        {
          // get fe values of the basis functions in the quad. points
          // in the current element (pointer to pointer)
          OrigFEValuesNeigh = TFEDatabase3D::GetOrigElementValues(BaseFunctNeigh,
                                                       NeededDerivatives[k]);
          // accumulate the values of the fe function 
          // do it for all quad. points
          for(j=0;j<N_PointsNeigh;j++)
          {
            // get pointer to array where the values for quad. point j are
            // stored
            Orig = OrigFEValuesNeigh[j];
            // initialize 
            value = value1 = value2 = 0;
            // compute linear combination over all local dof
            for(lll=0;lll<N_Neigh;lll++)
            {
              // value += u1[l] * phi[l]
              value += FEFunctValuesNeigh[lll] * Orig[lll];
              value1 += FEFunctValues1Neigh[lll] * Orig[lll];
              value2 += FEFunctValues2Neigh[lll] * Orig[lll];
            } // endfor lll
            // copy result to position in array 

            switch(k)
            {
              case 0 :
                FEValue_D1u1_QuadPoint[j] = value;
                FEValue_D1u2_QuadPoint[j] = value1;
                FEValue_D1u3_QuadPoint[j] = value2;
                break;
              case 1 :
                FEValue_D2u1_QuadPoint[j] = value;
                FEValue_D2u2_QuadPoint[j] = value1;
                FEValue_D2u3_QuadPoint[j] = value2;
                break;
              case 2 :
                FEValue_D3u1_QuadPoint[j] = value;
                FEValue_D3u2_QuadPoint[j] = value1;
                FEValue_D3u3_QuadPoint[j] = value2;
                break;
            }
          } // endfor j
        } // endfor k
        
        // compute part of the convolution in all local fe nodes (dof)
        // which comes from current mesh cell 
        // loop over all local dof
        
        // initialize
        value11 = value12 = value13 = value22 = value23 = value33 = 0;
        // loop over all quadrature points 
        for(j=0;j<N_PointsNeigh;j++) 
        { 
          // square of the distance between local dof and quad. point
          // in original cell 
          switch(TDatabase::ParamDB->INTERNAL_PERIODIC_IDENTITY)
          {
            case 0:
              distance_sq = ((X_orig[l]-XNeigh[j])*(X_orig[l]-XNeigh[j])+
                             (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j])+
                             (Z_orig[l]-ZNeigh[j])*(Z_orig[l]-ZNeigh[j])); 
              break;
            case 1:
              // for mixing layer problem in 3D
              if ((fabs(X_orig[l]-XNeigh[j])> 1.5)||(fabs(Z_orig[l]-ZNeigh[j])>1.5))
              {
                 // OutPut("points "<<  X_orig[l] << " " << Y_orig[l] << " " << Z_orig[l] << 
                 //       " ; " <<  XNeigh[j] << " " << YNeigh[j] << " " << ZNeigh[j] << " ");
                if (X_orig[l]-XNeigh[j]>1.5)
                  distance_sq = ((X_orig[l]-XNeigh[j]-2)*(X_orig[l]-XNeigh[j]-2)+
                                 (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j])+
                                 (Z_orig[l]-ZNeigh[j])*(Z_orig[l]-ZNeigh[j])); 
                if (X_orig[l]-XNeigh[j]<-1.5)
                  distance_sq = ((X_orig[l]-XNeigh[j]+2)*(X_orig[l]-XNeigh[j]+2)+
                                 (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j])+
                                 (Z_orig[l]-ZNeigh[j])*(Z_orig[l]-ZNeigh[j])); 
                if (Z_orig[l]-ZNeigh[j]>1.5)
                  distance_sq = ((X_orig[l]-XNeigh[j])*(X_orig[l]-XNeigh[j])+
                                 (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j])+
                                 (Z_orig[l]-ZNeigh[j]-2)*(Z_orig[l]-ZNeigh[j]-2)); 
                if (Z_orig[l]-ZNeigh[j]<-1.5)
                  distance_sq = ((X_orig[l]-XNeigh[j])*(X_orig[l]-XNeigh[j])+
                                 (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j])+
                                 (Z_orig[l]-ZNeigh[j]+2)*(Z_orig[l]-ZNeigh[j]+2)); 
                //OutPut(distance_sq << endl);
              }
              else
                distance_sq = ((X_orig[l]-XNeigh[j])*(X_orig[l]-XNeigh[j])+
                               (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j])+
                               (Z_orig[l]-ZNeigh[j])*(Z_orig[l]-ZNeigh[j])); 
              break;
            default :
              OutPut("Not implemented for this INTERNAL_PERIODIC_IDENTITY !!!" << endl);
              exit(4711);
                
         }
          //cout << "X " << XNeigh[j] << " Y " << YNeigh[j] << endl;
          // compute value of GaussianFilter
          // input: delta, distance_sq
          g = GaussianFilter3D(delta,distance_sq);
          //cout << "dist " <<  sqrt(distance_sq) << endl;
          //
          // apply quadrature rule in ref. cell 
          value11 += g*AbsDetjk[j]*weightsNeigh[j]*
            ( FEValue_D1u1_QuadPoint[j]*FEValue_D1u1_QuadPoint[j]
              +FEValue_D2u1_QuadPoint[j]*FEValue_D2u1_QuadPoint[j]
              +FEValue_D3u1_QuadPoint[j]*FEValue_D3u1_QuadPoint[j]);
          value12 += g*AbsDetjk[j]*weightsNeigh[j]*
            ( FEValue_D1u1_QuadPoint[j]*FEValue_D1u2_QuadPoint[j]
              +FEValue_D2u1_QuadPoint[j]*FEValue_D2u2_QuadPoint[j]
              +FEValue_D3u1_QuadPoint[j]*FEValue_D3u2_QuadPoint[j]);
          value13 += g*AbsDetjk[j]*weightsNeigh[j]*
            ( FEValue_D1u1_QuadPoint[j]*FEValue_D1u3_QuadPoint[j]
              +FEValue_D2u1_QuadPoint[j]*FEValue_D2u3_QuadPoint[j]
              +FEValue_D3u1_QuadPoint[j]*FEValue_D3u3_QuadPoint[j]);
          value22 += g*AbsDetjk[j]*weightsNeigh[j]*
            ( FEValue_D1u2_QuadPoint[j]*FEValue_D1u2_QuadPoint[j]
              +FEValue_D2u2_QuadPoint[j]*FEValue_D2u2_QuadPoint[j]
              +FEValue_D3u2_QuadPoint[j]*FEValue_D3u2_QuadPoint[j]);
          value23 += g*AbsDetjk[j]*weightsNeigh[j]*
            ( FEValue_D1u2_QuadPoint[j]*FEValue_D1u3_QuadPoint[j]
              +FEValue_D2u2_QuadPoint[j]*FEValue_D2u3_QuadPoint[j]
              +FEValue_D3u2_QuadPoint[j]*FEValue_D3u3_QuadPoint[j]);
          value33 += g*AbsDetjk[j]*weightsNeigh[j]*
            ( FEValue_D1u3_QuadPoint[j]*FEValue_D1u3_QuadPoint[j]
              +FEValue_D2u3_QuadPoint[j]*FEValue_D2u3_QuadPoint[j]
              +FEValue_D3u3_QuadPoint[j]*FEValue_D3u3_QuadPoint[j]);
        }
        // add local result to global array
        u_conv[index] += value11;
        u_conv[index+N_UConv] += value12;
        u_conv[index+2*N_UConv] += value13;
        u_conv[index+3*N_UConv] += value22;
        u_conv[index+4*N_UConv] += value23;
        u_conv[index+5*N_UConv] += value33;
      } // end of ll loop (edges)
    } // end of l loop (local dof)
  } // endfor i (cells)
  
  // for (i=0;i<N_DOFConv; i++)
  // {
  //   cout << u_conv[i] << "  " << u_conv[i+N_UConv] << "  ";
  //   cout << u_conv[i+2*N_UConv] << "  " << u_conv[i+3*N_UConv] << endl;
  // }

  // release memory
  delete x_conv;
  delete y_conv;
  delete z_conv;
  TDatabase::ParamDB->INTERNAL_QUAD_RULE=0;

} // ConvoluteSymmetricTensor

// ========================================================================
// convolute D = (grad w grad w^T) by full integration
// use D_{21} = D_{12}
// input u
// output DuTensor
// ========================================================================

void  ConvoluteSymmetricTensorFull3D(TFEVectFunct3D *u, TFEVectFunct3D *duTensor)

{
  int i,j,k,l,n,m, N_UsedElements, N_LocalUsedElements;
  int N_LocalUsedElementsNeigh;
  int N_Cells, N_Points, N_Parameters, N_, N_U, N_DOF, N_loc_dofConv;
  int Used[N_FEs3D], UsedNeigh[N_FEs3D], *N_BaseFunct;
  int neigh_i, N_Neigh, ll, lll, N_PointsNeigh, N_Edges, same_dof;
  FE3D LocalUsedElements[N_FEs3D], CurrentElement, CurrentElementConv;
  FE3D CurrentElementNeigh,LocalUsedElementsNeigh[N_FEs3D] ;
  TFE3D *Element;
  TNodalFunctional3D *nf;
  BaseFunct3D BaseFunct, *BaseFuncts, BaseFunctNeigh;
  TCollection *Coll;
  TBaseCell *cell, *neigh;
  TFE3D *ele;
  double *weights, *xi, *eta, *zeta, *xi_ref, *eta_ref, *zeta_ref;
  double *weightsNeigh, *xiNeigh, *etaNeigh, *zetaNeigh;
  double X_orig[MaxN_PointsForNodal3D], Y_orig[MaxN_PointsForNodal3D];
  double Z_orig[MaxN_PointsForNodal3D];
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D], AbsDetjkConv[MaxN_QuadPoints_3D];
  double XNeigh[MaxN_QuadPoints_3D], YNeigh[MaxN_QuadPoints_3D];
  double ZNeigh[MaxN_QuadPoints_3D];
  double AbsDetjkNeigh[MaxN_QuadPoints_3D];
  RefTrans3D RefTrans;
  double *Param[MaxN_QuadPoints_3D], *aux;
  double FEValue_D1u1_QuadPoint[MaxN_QuadPoints_3D];
  double FEValue_D1u2_QuadPoint[MaxN_QuadPoints_3D];
  double FEValue_D1u3_QuadPoint[MaxN_QuadPoints_3D];
  double FEValue_D2u1_QuadPoint[MaxN_QuadPoints_3D];
  double FEValue_D2u2_QuadPoint[MaxN_QuadPoints_3D];
  double FEValue_D2u3_QuadPoint[MaxN_QuadPoints_3D];
  double FEValue_D3u1_QuadPoint[MaxN_QuadPoints_3D];
  double FEValue_D3u2_QuadPoint[MaxN_QuadPoints_3D];
  double FEValue_D3u3_QuadPoint[MaxN_QuadPoints_3D];
  double FEValue_D1u1_QuadPointNeigh[MaxN_QuadPoints_3D];
  double FEValue_D1u2_QuadPointNeigh[MaxN_QuadPoints_3D];
  double FEValue_D1u3_QuadPointNeigh[MaxN_QuadPoints_3D];
  double FEValue_D2u1_QuadPointNeigh[MaxN_QuadPoints_3D];
  double FEValue_D2u2_QuadPointNeigh[MaxN_QuadPoints_3D];
  double FEValue_D2u3_QuadPointNeigh[MaxN_QuadPoints_3D];
  double FEValue_D3u1_QuadPointNeigh[MaxN_QuadPoints_3D];
  double FEValue_D3u2_QuadPointNeigh[MaxN_QuadPoints_3D];
  double FEValue_D3u3_QuadPointNeigh[MaxN_QuadPoints_3D];
  double *ExactVal[MaxN_QuadPoints_3D];
  int *DOF, ActiveBound, DirichletBound, end, last, *DOFConv, *DOFNeigh;
  double **OrigFEValues, *Orig;
  double value, value1, value2;
  double value11, value12, value13, value22, value23, value33;
  double **OrigFEValuesNeigh;
  double FEFunctValues[MaxN_BaseFunctions3D];
  double FEFunctValues1[MaxN_BaseFunctions3D];
  double FEFunctValues2[MaxN_BaseFunctions3D];
  double FEFunctValuesNeigh[MaxN_BaseFunctions3D];
  double FEFunctValues1Neigh[MaxN_BaseFunctions3D];
  double FEFunctValues2Neigh[MaxN_BaseFunctions3D];
  int *GlobalNumbers, *BeginIndex, *GlobalNumbersConv, *BeginIndexConv;
  double LocError[4];
  int N_UConv, N_DOFConv, index, *conv_comp;
  double *Values, *u_conv,distance_sq,h_fine,g, *u_values;
  double *x_conv, *y_conv, *z_conv, delta, hK;
  TJoint *joint;
  TFESpace3D *fespace, *fespaceConv; 
  TVertex *vertex0;
  double x_vertex, y_vertex, z_vertex;

  bool SecondDer[1] = { FALSE };
  int N_Derivatives = 3;
  MultiIndex3D NeededDerivatives[3] = { D100, D010, D001 };
  TDatabase::ParamDB->INTERNAL_QUAD_RULE=1;

  // get fe spaces of velocity u and of covoluted tensor
  fespace = u->GetFESpace3D();
  fespaceConv = duTensor->GetFESpace3D();

  // gives a pointer to all available basis functions which
  // are described in the fedatabase
  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  // gives the number of available basis functions from the
  // fedatabase
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  // get information of the numbering of the degrees of freedom
  // first: get pointer to array where a global numbering is stored
  GlobalNumbers = fespace->GetGlobalNumbers();
  // second: get pointer to array where the start of numbering for
  // each cell is stored
  BeginIndex = fespace->GetBeginIndex();
 
  // get number of dof for one velocity component N_U
  // -> total number of dof is twice N_U
  N_U = u->GetLength();
  N_DOF = 3*N_U;
  
  // get pointer to set of mesh cells which define the fe space
  Coll = fespace->GetCollection();
  // get number of mesh cells
  N_Cells = Coll->GetN_Cells();
 
  // # dof in the new fe space
  N_UConv = fespaceConv->GetN_DegreesOfFreedom();
 
  N_DOFConv = 6*N_UConv;
  // allocate memory for values of convoluted function
  u_conv = duTensor->GetValues(); 
  // initialize u_conv to 0
  memset(u_conv,0,N_DOFConv*SizeOfDouble);

  // allocate memory for values of convoluted function
  conv_comp = new int[N_UConv]; 
  // initialize conv_comp to 0
  memset(conv_comp,0,N_UConv*SizeOfInt);

  // get information of the numbering of the dof for the space of convolution
  // first: get pointer to array where a global numbering is stored
  GlobalNumbersConv = fespaceConv->GetGlobalNumbers();
  // second: get pointer to array where the start of numbering for
  // each cell is stored
  BeginIndexConv = fespaceConv->GetBeginIndex();

  // allocate memory for geometric position of nodes
  x_conv = new double[N_DOFConv]; 
  y_conv = new double[N_DOFConv]; 
  z_conv = new double[N_DOFConv]; 
   
// ########################################################################
// computing the convolution
// loop over all cells
// ########################################################################

  // direct a pointer to the values of the finite element vector
  // function
  // first part of the array : u1
  // second part of the array : u2
  Values = u->GetValues();           

  // store positition of each mesh cells in the collection of mesh 
  // cells: set ClipBoard to i
  for(i=0;i<N_Cells;i++)                       
  {
    cell = Coll->GetCell(i);
    cell->SetClipBoard(i);
  }

  // loop over all mesh cells to accumulate the values for the 
  // convolution 
  for(i=0;i<N_Cells;i++)
  {
    // get current mesh cell
    cell = Coll->GetCell(i);        
    hK = cell->GetDiameter();
    delta =  CharacteristicFilterWidth(hK);

    // ####################################################################
    // find local used finite elements on this cell
    // ####################################################################
    CurrentElement = fespace->GetFE3D(i, cell);
    LocalUsedElements[0] = CurrentElement;
    N_LocalUsedElements = 1;
    
    // ####################################################################
    // calculate values on original element
    // ####################################################################
    
    // From id of fe space the situation in the referenz mesh cell is 
    // known. Now, get information in the original mesh cell. 
    // Every finite element gets automatically a quadrature rule.
    // input: N_LocalUsedElements, LocalUsedElements, cell, SecondDer
    // output : N_Points - number of quadrature points
    //          xi - xi-values of quad. points in ref. cell
    //          eta - eta-values of quad. points in ref. cell
    //          weigths - weights fro quad. rule in quad. points 
    //          X - X-values of quad. points in original cell  
    //          Y - Y-values of quad. points in original cell  
    //          AbsDetjk - absolute value of Jacobian in quad. points
    //          xi, eta, weights, X, Y,  AbsDetjk are pointers to arrays

    TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                           Coll, cell, SecondDer, N_Points, 
                           xi, eta, zeta, weights, X, Y, Z, AbsDetjk);
 
    // compute geometric positions of the fe nodes
    // get id of finite element in current mesh cell
    CurrentElement = fespace->GetFE3D(i, cell);
    CurrentElementConv = fespaceConv->GetFE3D(i, cell);
    // get fe from its id
    Element = TFEDatabase3D::GetFE3D(CurrentElementConv);
    // get reference transformation 
    RefTrans = Element->GetRefTransID();
    // get pointer to the nodal functionals (fe nodes) of the fe 
    // (in ref. cell)
    nf = Element->GetNodalFunctional3D();
    // get number and coordinates of local dof in ref cell
    // xi_ref, eta_ref are pointers
    nf->GetPointsForAll(N_loc_dofConv, xi_ref, eta_ref, zeta_ref);
      
    // get coordinates of fe nodes in original cell
    // input: RefTrans, N_loc_dof, xi_ref, eta_ref
    // output : X_orig, Y_orig - pointers to arrays with coordinates
    //          AbsDetjk - same as above    
    TFEDatabase3D::GetOrigFromRef(RefTrans,N_loc_dofConv, 
                                  xi_ref, eta_ref, zeta_ref,
                                  X_orig, Y_orig, Z_orig, AbsDetjkConv);   
   
    // get pointer to basis functions of the current element
    BaseFunct = BaseFuncts[CurrentElement];  
    // # of basis functions, is the same as N_loc_dof
    N_ = N_BaseFunct[CurrentElement];        

    // find the part of the global index array where the information
    // for the current mesh cell are stored: DOF
    DOF = GlobalNumbers + BeginIndex[i];
    DOFConv = GlobalNumbersConv + BeginIndexConv[i];
    // copy the values of the finite element function from the global
    // array (Values) to local arrays 
    for(l=0;l<N_;l++)
    {
      FEFunctValues[l] = Values[DOF[l]];
      FEFunctValues1[l] = Values[DOF[l]+N_U];
      FEFunctValues2[l] = Values[DOF[l]+2*N_U];
    }

    // compute values of fe function in the quadrature points with
    // linear combinations of the local fe values and the basis functions
    for(k=0;k<N_Derivatives;k++)
    {
      // get fe values of the basis functions in the quad. points
      // in the current element (pointer to pointer)
      OrigFEValues = TFEDatabase3D::GetOrigElementValues(BaseFunct,
                                            NeededDerivatives[k]);
      // accumulate the values of the fe function 
      // do it for all quad. points
      for(j=0;j<N_Points;j++)
      {
        // get pointer to array where the values for quad. point j are
        // stored
        Orig = OrigFEValues[j];
        // initialize 
        value = value1 = value2 =0;
        // compute linear combination over all local dof
        for(l=0;l<N_;l++)
        {
          // value += u1[l] * phi[l]
          value += FEFunctValues[l] * Orig[l];
          value1 += FEFunctValues1[l] * Orig[l];
          value2 += FEFunctValues2[l] * Orig[l];
        } // endfor l
        // copy result to position in array 

        switch(k)
        {
          case 0 :
          FEValue_D1u1_QuadPoint[j] = value;
          FEValue_D1u2_QuadPoint[j] = value1;
          FEValue_D1u3_QuadPoint[j] = value2;
          break;
          case 1 :
          FEValue_D2u1_QuadPoint[j] = value;
          FEValue_D2u2_QuadPoint[j] = value1;
          FEValue_D2u3_QuadPoint[j] = value2;
          break;
          case 2 :
          FEValue_D3u1_QuadPoint[j] = value;
          FEValue_D3u2_QuadPoint[j] = value1;
          FEValue_D3u3_QuadPoint[j] = value2;
          break;
        }
      } // endfor j
    } // endfor k
    
    // compute part of the convolution in all local fe nodes (dof)
    // which comes from current mesh cell 
    // loop over all local dof
    for(l=0;l<N_loc_dofConv;l++)         
    {
      // get global indes of local dof
      index = DOFConv[l];          
      // check if value is already computed
      // if yes : continue
      if (conv_comp[index])
        continue;

      // set conv_comp[index]
      conv_comp[index] = 1;

      // initialize
      value11 = value12 = value13 = value22 = value23 = value33 = 0;
      // loop over all quadrature points 
      for(j=0;j<N_Points;j++) 
      {
        // square of the distance between local dof and quad. point
        // in original cell 
        distance_sq = ((X_orig[l]-X[j])*(X_orig[l]-X[j])+
                       (Y_orig[l]-Y[j])*(Y_orig[l]-Y[j])+
                       (Z_orig[l]-Z[j])*(Z_orig[l]-Z[j])); 
        // compute value of GaussianFilter
        // input: delta, distance_sq
        g = GaussianFilter3D(delta,distance_sq);  
        // apply quadrature rule in ref. cell 
        value11 += g*AbsDetjk[j]*weights[j]*
                   ( FEValue_D1u1_QuadPoint[j]*FEValue_D1u1_QuadPoint[j]
                    +FEValue_D2u1_QuadPoint[j]*FEValue_D2u1_QuadPoint[j]
                    +FEValue_D3u1_QuadPoint[j]*FEValue_D3u1_QuadPoint[j]);
        value12 += g*AbsDetjk[j]*weights[j]*
                   ( FEValue_D1u1_QuadPoint[j]*FEValue_D1u2_QuadPoint[j]
                     +FEValue_D2u1_QuadPoint[j]*FEValue_D2u2_QuadPoint[j]
                     +FEValue_D3u1_QuadPoint[j]*FEValue_D3u2_QuadPoint[j]);
        value13 += g*AbsDetjk[j]*weights[j]*
                   ( FEValue_D1u1_QuadPoint[j]*FEValue_D1u3_QuadPoint[j]
                     +FEValue_D2u1_QuadPoint[j]*FEValue_D2u3_QuadPoint[j]
                     +FEValue_D3u1_QuadPoint[j]*FEValue_D3u3_QuadPoint[j]);
        value22 += g*AbsDetjk[j]*weights[j]*
                   ( FEValue_D1u2_QuadPoint[j]*FEValue_D1u2_QuadPoint[j]
                     +FEValue_D2u2_QuadPoint[j]*FEValue_D2u2_QuadPoint[j]
                     +FEValue_D3u2_QuadPoint[j]*FEValue_D3u2_QuadPoint[j]);
        value23 += g*AbsDetjk[j]*weights[j]*
                   ( FEValue_D1u2_QuadPoint[j]*FEValue_D1u3_QuadPoint[j]
                     +FEValue_D2u2_QuadPoint[j]*FEValue_D2u3_QuadPoint[j]
                     +FEValue_D3u2_QuadPoint[j]*FEValue_D3u3_QuadPoint[j]);
        value33 += g*AbsDetjk[j]*weights[j]*
                   ( FEValue_D1u3_QuadPoint[j]*FEValue_D1u3_QuadPoint[j]
                     +FEValue_D2u3_QuadPoint[j]*FEValue_D2u3_QuadPoint[j]
                     +FEValue_D3u3_QuadPoint[j]*FEValue_D3u3_QuadPoint[j]);
      }
      // add local result to global array
      u_conv[index] += value11;
      u_conv[index+N_UConv] += value12;
      u_conv[index+2*N_UConv] += value13;
      u_conv[index+3*N_UConv] += value22;
      u_conv[index+4*N_UConv] += value23;
      u_conv[index+5*N_UConv] += value33;
      // store global coordinates of the fe nodes (dof)
      x_conv[index] = X_orig[l];
      y_conv[index] = Y_orig[l];
      z_conv[index] = Z_orig[l];
      x_conv[index+N_UConv] = X_orig[l];
      y_conv[index+N_UConv] = Y_orig[l];
      z_conv[index+N_UConv] = Z_orig[l];
      x_conv[index+2*N_UConv] = X_orig[l];
      y_conv[index+2*N_UConv] = Y_orig[l];
      z_conv[index+2*N_UConv] = Z_orig[l];
      
      // starting to compute the integral in ALL the cells
      
      // loop over all cells
      for(ll=0;ll<N_Cells;ll++)                           
      {
        // check if we are at the current mesh cell
        if (ll==i)
          continue;

        neigh = Coll->GetCell(ll); 
        // get one vertex of the neigh 
        vertex0 = neigh->GetVertex(0);
        // get coordinates of vertex
        x_vertex = vertex0->GetX();
        y_vertex = vertex0->GetY();
        z_vertex = vertex0->GetZ();
        // compute distance to current point where the convolution
        // is computed
        distance_sq = ((X_orig[l]-x_vertex)*(X_orig[l]-x_vertex)+
                       (Y_orig[l]-y_vertex)*(Y_orig[l]-y_vertex)+
                       (Z_orig[l]-z_vertex)*(Z_orig[l]-z_vertex));
        // check if neigh is close enough
        if (distance_sq > 2*delta *delta)
          continue;
        // get id of neighbour
        neigh_i = neigh->GetClipBoard();

        // get id of finite element in neigh mesh cell wrt the finite
        // element space of the velocity
        CurrentElementNeigh = fespace->GetFE3D(neigh_i, neigh);
        
        // get pointer to basis functions of the neighbour element
        BaseFunctNeigh = BaseFuncts[CurrentElementNeigh]; 
 
        // # of basis functions, is the same as N_loc_dof
        N_Neigh = N_BaseFunct[CurrentElementNeigh]; 
     
        // find the part of the global index array where the information
        // for the neigh mesh cell are stored: DOF
        DOFNeigh = GlobalNumbers + BeginIndex[neigh_i];

        // initialize auxiliary array UsedNeigh
        memset(UsedNeigh, 0, N_FEs3D*SizeOfInt);
        
        // set entry in auxiliary array to 1
        UsedNeigh[CurrentElementNeigh] = 1;
       
        // have a look how much 1 are in UsedNeigh -> N_LocalUsedElementsNeigh
        N_LocalUsedElementsNeigh = 0;
        memset(LocalUsedElementsNeigh, 0, SizeOfInt*N_FEs3D);
        j = 0;
        for(k=0;k<N_FEs3D;k++)
          if(UsedNeigh[k])
          {
            LocalUsedElementsNeigh[j] = (FE3D)k;
            j++;
          }
        N_LocalUsedElementsNeigh = j;
      
        // ####################################################################
        // calculate values on original element
        // ####################################################################
        
        // From id of fe space the situation in the referenz mesh cell is 
        // known. Now, get information in the original mesh cell. 
        // Every finite element gets automatically a quadrature rule.
        // input: N_LocalUsedElementsNeigh, LocalUsedElementsNeigh, 
        //        neigh, SecondDer
        // output : N_PointsNeigh - number of quadrature points
        //          xiNeigh - xi-values of quad. points in ref. cell
        //          etaNeigh - eta-values of quad. points in ref. cell
        //          weigthsNeigh - weights fro quad. rule in quad. points 
        //          XNeigh - X-values of quad. points in original cell  
        //          YNeigh - Y-values of quad. points in original cell  
        //          AbsDetjkNeigh - absolute value of Jacobian in quad. points
        //          xiNeigh, etaNeigh, weightsNeigh, XNeigh, YNeigh,  
        //          AbsDetjkNeigh are pointers to arrays

        TFEDatabase3D::GetOrig(N_LocalUsedElementsNeigh, LocalUsedElementsNeigh, 
                               Coll, neigh, SecondDer, N_PointsNeigh, 
                               xiNeigh, etaNeigh, zetaNeigh, weightsNeigh,
                               XNeigh, YNeigh, ZNeigh, AbsDetjkNeigh);
    
        // copy the values of the finite element function from the global
        // array (Values) to local arrays 
        for(lll=0;lll<N_Neigh;lll++)
        {
          FEFunctValuesNeigh[lll] = Values[DOFNeigh[lll]];
          FEFunctValues1Neigh[lll] = Values[DOFNeigh[lll]+N_U];
          FEFunctValues2Neigh[lll] = Values[DOFNeigh[lll]+2*N_U];
        }
        
        // compute values of fe function in the quadrature points with
        // linear combinations of the local fe values and the basis functions
        for(k=0;k<N_Derivatives;k++)
        {
          // get fe values of the basis functions in the quad. points
          // in the current element (pointer to pointer)
          OrigFEValuesNeigh = TFEDatabase3D::GetOrigElementValues(BaseFunctNeigh,
                                                       NeededDerivatives[k]);
          // accumulate the values of the fe function 
          // do it for all quad. points
          for(j=0;j<N_PointsNeigh;j++)
          {
            // get pointer to array where the values for quad. point j are
            // stored
            Orig = OrigFEValuesNeigh[j];
            // initialize 
            value = value1 = value2 = 0;
            // compute linear combination over all local dof
            for(lll=0;lll<N_Neigh;lll++)
            {
              // value += u1[l] * phi[l]
              value += FEFunctValuesNeigh[lll] * Orig[lll];
              value1 += FEFunctValues1Neigh[lll] * Orig[lll];
              value2 += FEFunctValues2Neigh[lll] * Orig[lll];
            } // endfor lll
            // copy result to position in array 

            switch(k)
            {
              case 0 :
                FEValue_D1u1_QuadPoint[j] = value;
                FEValue_D1u2_QuadPoint[j] = value1;
                FEValue_D1u3_QuadPoint[j] = value2;
                break;
              case 1 :
                FEValue_D2u1_QuadPoint[j] = value;
                FEValue_D2u2_QuadPoint[j] = value1;
                FEValue_D2u3_QuadPoint[j] = value2;
                break;
              case 2 :
                FEValue_D3u1_QuadPoint[j] = value;
                FEValue_D3u2_QuadPoint[j] = value1;
                FEValue_D3u3_QuadPoint[j] = value2;
                break;
            }
          } // endfor j
        } // endfor k
        
        // compute part of the convolution in all local fe nodes (dof)
        // which comes from current mesh cell 
        // loop over all local dof
        
        // initialize
        value11 = value12 = value13 = value22 = value23 = value33 = 0;
        // loop over all quadrature points 
        for(j=0;j<N_PointsNeigh;j++) 
        { 
          // square of the distance between local dof and quad. point
          // in original cell 
          distance_sq = ((X_orig[l]-XNeigh[j])*(X_orig[l]-XNeigh[j])+
                         (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j])+
                         (Z_orig[l]-ZNeigh[j])*(Z_orig[l]-ZNeigh[j])); 
          //cout << "X " << XNeigh[j] << " Y " << YNeigh[j] << endl;
          // compute value of GaussianFilter
          // input: delta, distance_sq
          g = GaussianFilter3D(delta,distance_sq);
          //cout << "dist " <<  sqrt(distance_sq) << endl;
          //
          // apply quadrature rule in ref. cell 
          value11 += g*AbsDetjk[j]*weightsNeigh[j]*
            ( FEValue_D1u1_QuadPoint[j]*FEValue_D1u1_QuadPoint[j]
              +FEValue_D2u1_QuadPoint[j]*FEValue_D2u1_QuadPoint[j]
              +FEValue_D3u1_QuadPoint[j]*FEValue_D3u1_QuadPoint[j]);
          value12 += g*AbsDetjk[j]*weightsNeigh[j]*
            ( FEValue_D1u1_QuadPoint[j]*FEValue_D1u2_QuadPoint[j]
              +FEValue_D2u1_QuadPoint[j]*FEValue_D2u2_QuadPoint[j]
              +FEValue_D3u1_QuadPoint[j]*FEValue_D3u2_QuadPoint[j]);
          value13 += g*AbsDetjk[j]*weightsNeigh[j]*
            ( FEValue_D1u1_QuadPoint[j]*FEValue_D1u3_QuadPoint[j]
              +FEValue_D2u1_QuadPoint[j]*FEValue_D2u3_QuadPoint[j]
              +FEValue_D3u1_QuadPoint[j]*FEValue_D3u3_QuadPoint[j]);
          value22 += g*AbsDetjk[j]*weightsNeigh[j]*
            ( FEValue_D1u2_QuadPoint[j]*FEValue_D1u2_QuadPoint[j]
              +FEValue_D2u2_QuadPoint[j]*FEValue_D2u2_QuadPoint[j]
              +FEValue_D3u2_QuadPoint[j]*FEValue_D3u2_QuadPoint[j]);
          value23 += g*AbsDetjk[j]*weightsNeigh[j]*
            ( FEValue_D1u2_QuadPoint[j]*FEValue_D1u3_QuadPoint[j]
              +FEValue_D2u2_QuadPoint[j]*FEValue_D2u3_QuadPoint[j]
              +FEValue_D3u2_QuadPoint[j]*FEValue_D3u3_QuadPoint[j]);
          value33 += g*AbsDetjk[j]*weightsNeigh[j]*
            ( FEValue_D1u3_QuadPoint[j]*FEValue_D1u3_QuadPoint[j]
              +FEValue_D2u3_QuadPoint[j]*FEValue_D2u3_QuadPoint[j]
              +FEValue_D3u3_QuadPoint[j]*FEValue_D3u3_QuadPoint[j]);
        }
        // add local result to global array
        u_conv[index] += value11;
        u_conv[index+N_UConv] += value12; 
        u_conv[index+2*N_UConv] += value22; 
        u_conv[index+3*N_UConv] += value22;
        u_conv[index+4*N_UConv] += value23;
        u_conv[index+5*N_UConv] += value33;
      } // end of ll loop (edges)
    } // end of l loop (local dof)
  } // endfor i (cells)
  
  /* for (i=0;i<N_UConv; i++)
    {
      cout << x_conv[i] << " " << y_conv[i];
      cout << "    " << Values[i] << " "<< Values[i+N_U] << "     ";
      cout << u_conv[i] << "  " << u_conv[i+N_UConv] << "  ";
      cout << u_conv[i+2*N_UConv] << endl;
      }*/

  // release memory
  delete x_conv;
  delete y_conv;
  delete z_conv;
  delete conv_comp;
  TDatabase::ParamDB->INTERNAL_QUAD_RULE=0;

} // ConvoluteSymmetricTensorFull

// same as velocity convolution, only one component
void  ConvolutePressure3D(TFEFunction3D *u, TFEFunction3D *uConv)
{
  int i,j,k,l,n,m, N_UsedElements, N_LocalUsedElements;
  int N_LocalUsedElementsNeigh;
  int N_Cells, N_Points, N_Parameters, N_, N_U, N_DOF, N_loc_dofConv;
  int Used[N_FEs3D], UsedNeigh[N_FEs3D], *N_BaseFunct;
  int neigh_i, N_Neigh, ll, lll, N_PointsNeigh, N_Faces, same_dof;
  FE3D LocalUsedElements[N_FEs3D], CurrentElement, CurrentElementConv;
  FE3D CurrentElementNeigh,LocalUsedElementsNeigh[N_FEs3D] ;
  TFE3D *Element;
  TNodalFunctional3D *nf;
  BaseFunct3D BaseFunct, *BaseFuncts, BaseFunctNeigh;
  TCollection *Coll;
  TBaseCell *cell, *neigh;
  TFE3D *ele;
  double *weights, *xi, *eta, *zeta, *xi_ref, *eta_ref, *zeta_ref;
  double *weightsNeigh, *xiNeigh, *etaNeigh, *zetaNeigh;
  double X_orig[MaxN_PointsForNodal3D], Y_orig[MaxN_PointsForNodal3D];
  double Z_orig[MaxN_PointsForNodal3D];
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D], AbsDetjkConv[MaxN_QuadPoints_3D];
  double XNeigh[MaxN_QuadPoints_3D], YNeigh[MaxN_QuadPoints_3D];
  double ZNeigh[MaxN_QuadPoints_3D];
  double AbsDetjkNeigh[MaxN_QuadPoints_3D];
  RefTrans3D RefTrans;
  double *Param[MaxN_QuadPoints_3D], *aux;
  double FEValue_u1_QuadPoint[MaxN_QuadPoints_3D];
  double FEValue_u1_QuadPointNeigh[MaxN_QuadPoints_3D];
  double *ExactVal[MaxN_QuadPoints_3D];
  int *DOF, ActiveBound, DirichletBound, end, last, *DOFConv, *DOFNeigh;
  double **OrigFEValues, *Orig, value;
  double **OrigFEValuesNeigh;
  double FEFunctValues[MaxN_BaseFunctions3D];
  double FEFunctValuesNeigh[MaxN_BaseFunctions3D];
  int *GlobalNumbers, *BeginIndex, *GlobalNumbersConv, *BeginIndexConv;
  double LocError[4];
  double hK,delta;
  bool *SecondDer;
  int n_fespaces = 1, N_Derivatives = 1, index;
  int N_UConv, N_DOFConv;
  double *Values, *u_conv,distance_sq,h_fine,g, *u_values;
  double *x_conv, *y_conv, *z_conv;
  TJoint *joint;
  TFESpace3D *fespace, *fespaceConv; 
  TDatabase::ParamDB->INTERNAL_QUAD_RULE=1;

  // gives an array where the needed derivatives are described
  // defined in NavierStokes.h 
  // here we get as result {D00}
  MultiIndex3D NeededDerivatives[1] = { D000 };

  // get fe spaces of velocity u and of covoluted tensor
  fespace = u->GetFESpace3D();
  fespaceConv = uConv->GetFESpace3D();

  // gives a pointer to all available basis functions which
  // are described in the fedatabase
  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  // gives the number of available basis functions from the
  // fedatabase
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  // say that we don't need second derivatives
  // first: allocate memory for an array of length n_fespaces 
  // second: set all entries to FALSE
  SecondDer = new bool[n_fespaces];
  for(i=0;i<n_fespaces;i++)
    SecondDer[i] = FALSE;

  // get information of the numbering of the degrees of freedom
  // first: get pointer to array where a global numbering is stored
  GlobalNumbers = fespace->GetGlobalNumbers();
  // second: get pointer to array where the start of numbering for
  // each cell is stored
  BeginIndex = fespace->GetBeginIndex();
 
  // get number of dof for one velocity component N_U
  // -> total number of dof is twice N_U
  N_U = u->GetLength();
  N_DOF = N_U;
  
  // get pointer to set of mesh cells which define the fe space
  Coll = fespace->GetCollection();
 
  // # dof in the new fe space
  N_UConv = fespaceConv->GetN_DegreesOfFreedom();
 
  N_DOFConv = N_UConv;
  // allocate memory for values of convuluted function
  u_conv = uConv->GetValues(); 
  // initialize u_conv to 0
  memset(u_conv,0,N_DOFConv*SizeOfDouble);

  // get information of the numbering of the dof for the space of convolution
  // first: get pointer to array where a global numbering is stored
  GlobalNumbersConv = fespaceConv->GetGlobalNumbers();
  // second: get pointer to array where the start of numbering for
  // each cell is stored
  BeginIndexConv = fespaceConv->GetBeginIndex();

  // allocate memory for geometric position of nodes
  x_conv = new double[N_DOFConv]; 
  y_conv = new double[N_DOFConv]; 
  z_conv = new double[N_DOFConv]; 
   
// ########################################################################
// computing the convolution
// loop over all cells
// ########################################################################

  // get pointer to set of mesh cells which define the fe space
  Coll = fespace->GetCollection(); 
  // get number of mesh cells
  N_Cells = Coll->GetN_Cells();

  // direct a pointer to the values of the finite element vector
  // function
  // first part of the array : u1
  // second part of the array : u2
  Values = u->GetValues();           

  // store positition of each mesh cells in the collection of mesh 
  // cells: set ClipBoard to i
  for(i=0;i<N_Cells;i++)                       
  {
    cell = Coll->GetCell(i);
    cell->SetClipBoard(i);
  }

  // loop over all mesh cells to accumulate the values for the 
  // convolution 
  for(i=0;i<N_Cells;i++)
  {
    // get current mesh cell
    cell = Coll->GetCell(i);        
    // get cell diameter (longest distance between to vertices of the cell)
    hK = cell->GetDiameter();       
    delta = CharacteristicFilterWidth(hK);
 
    // ####################################################################
    // find local used finite elements on this cell
    // ####################################################################
    
    // have a look how much 1 are in Used ->   N_LocalUsedElements
    CurrentElement = fespace->GetFE3D(i, cell);
    LocalUsedElements[0] = CurrentElement;
    N_LocalUsedElements = 1;
    
    // ####################################################################
    // calculate values on original element
    // ####################################################################
    
    // From id of fe space the situation in the referenz mesh cell is 
    // known. Now, get information in the original mesh cell. 
    // Every finite element gets automatically a quadrature rule.
    // input: N_LocalUsedElements, LocalUsedElements, cell, SecondDer
    // output : N_Points - number of quadrature points
    //          xi - xi-values of quad. points in ref. cell
    //          eta - eta-values of quad. points in ref. cell
    //          weigths - weights fro quad. rule in quad. points 
    //          X - X-values of quad. points in original cell  
    //          Y - Y-values of quad. points in original cell  
    //          AbsDetjk - absolute value of Jacobian in quad. points
    //          xi, eta, weights, X, Y,  AbsDetjk are pointers to arrays

    TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, zeta,
                           weights, X, Y, Z, AbsDetjk);
 
    // compute geometric positions of the fe nodes
    // get id of finite element in current mesh cell
    CurrentElement = fespace->GetFE3D(i, cell);
    CurrentElementConv = fespaceConv->GetFE3D(i, cell);
    // get fe from its id
    Element = TFEDatabase3D::GetFE3D(CurrentElementConv);
    // get reference transformation 
    RefTrans = Element->GetRefTransID();
    // get pointer to the nodal functionals (fe nodes) of the fe 
    // (in ref. cell)
    nf = Element->GetNodalFunctional3D();
    // get number and coordinates of local dof in ref cell
    // xi_ref, eta_ref, zeta_ref are pointers
    nf->GetPointsForAll(N_loc_dofConv, xi_ref, eta_ref, zeta_ref);
      
    // get coordinates of fe nodes in original cell
    // input: RefTrans, N_loc_dof, xi_ref, eta_ref
    // output : X_orig, Y_orig - pointers to arrays with coordinates
    //          AbsDetjk - same as above    
    TFEDatabase3D::GetOrigFromRef(RefTrans,N_loc_dofConv, 
                                  xi_ref, eta_ref, zeta_ref,
                                  X_orig, Y_orig, Z_orig, AbsDetjkConv);   
   
    // get pointer to basis functions of the current element
    BaseFunct = BaseFuncts[CurrentElement];  
    // # of basis functions, is the same as N_loc_dof
    N_ = N_BaseFunct[CurrentElement];        

    // find the part of the global index array where the information
    // for the current mesh cell are stored: DOF
    DOF = GlobalNumbers + BeginIndex[i];
    DOFConv = GlobalNumbersConv + BeginIndexConv[i];
    // copy the values of the finite element function from the global
    // array (Values) to local arrays 
    for(l=0;l<N_;l++)
    {
      FEFunctValues[l] = Values[DOF[l]];
    }

    // compute values of fe function in the quadrature points with
    // linear combinations of the local fe values and the basis functions
    for(k=0;k<N_Derivatives;k++)
    {
      // get fe values of the basis functions in the quad. points
      // in the current element (pointer to pointer)
      OrigFEValues = TFEDatabase3D::GetOrigElementValues(BaseFunct,
                                            NeededDerivatives[k]);
      // accumulate the values of the fe function 
      // do it for all quad. points
      for(j=0;j<N_Points;j++)
      {
        // get pointer to array where the values for quad. point j are
        // stored
        Orig = OrigFEValues[j];
        // initialize 
        value = 0;
        // compute linear combination over all local dof
        for(l=0;l<N_;l++)
        {
          // value += u1[l] * phi[l]
          value += FEFunctValues[l] * Orig[l];
        } // endfor l
        // copy result to position in array 
        FEValue_u1_QuadPoint[j] = value;
      } // endfor j
    } // endfor k
    
    // compute part of the convolution in all local fe nodes (dof)
    // which comes from current mesh cell 
    // loop over all local dof
    for(l=0;l<N_loc_dofConv;l++)         
    {
      // get global indes of local dof
      index = DOFConv[l];          
      // initialize
      value = 0;
      // loop over all quadrature points 
      for(j=0;j<N_Points;j++) 
      {
        // square of the distance between local dof and quad. point
        // in original cell 
        distance_sq = ((X_orig[l]-X[j])*(X_orig[l]-X[j])+
                       (Y_orig[l]-Y[j])*(Y_orig[l]-Y[j])+
                       (Z_orig[l]-Z[j])*(Z_orig[l]-Z[j])); 
        // compute value of GaussianFilter
        // input: delta, distance_sq
        g = GaussianFilter3D(delta,distance_sq);  
        // apply quadrature rule in ref. cell 
        value += g*AbsDetjk[j]*weights[j]*FEValue_u1_QuadPoint[j];
      }
      // add local result to global array
      u_conv[index] += value;
      // store global coordinates of the fe nodes (dof)
      x_conv[index] = X_orig[l];
      y_conv[index] = Y_orig[l];
      z_conv[index] = Z_orig[l];
      
      // starting to compute the integral in the neighbour cells
      // number of faces
      N_Faces=cell->GetN_Faces();
      // loop over all faces
      for(ll=0;ll<N_Faces;ll++)                           
      {
        // get pointer to edge[l]
        joint=cell->GetJoint(ll); 
        
        // if boundary edge continue 
        if (joint->GetType() == BoundaryFace)
          continue;
        
        // get pointer to neighbour cell
        neigh=cell->GetJoint(ll)->GetNeighbour(cell);
        
        // get id of neighbour
        neigh_i = neigh->GetClipBoard();
        
        // get id of finite element in neigh mesh cell wrt the
        // finite element space of the convolution
        CurrentElementNeigh = fespaceConv->GetFE3D(neigh_i, neigh);
        
        // get pointer to basis functions of the neighbour element
        BaseFunctNeigh = BaseFuncts[CurrentElementNeigh]; 
 
        // # of basis functions, is the same as N_loc_dof
        N_Neigh = N_BaseFunct[CurrentElementNeigh]; 
     
        // find the part of the global index array where the information
        // for the neigh mesh cell are stored: DOF
        DOFNeigh = GlobalNumbersConv + BeginIndexConv[neigh_i];
        
        // check if computation on neigh cell necessary
        // for all dof
        same_dof = 0;
        for (j=0;j<N_Neigh;j++)
        {
          //cout << "N " << DOFNeigh[j] << " ind " << index << endl;
          if (DOFNeigh[j]==index)
          same_dof = 1;
        }
        // if the cells have the same dof, continue
        if (same_dof==1)
          continue;
        
        ////////////////////////////////////////////////////////////////////
        
        // if not, do something! 
        // i.e. the same thing as before
        
        // get id of finite element in neigh mesh cell wrt the finite
        // element space of the velocity
        CurrentElementNeigh = fespace->GetFE3D(neigh_i, neigh);
        
        // get pointer to basis functions of the neighbour element
        BaseFunctNeigh = BaseFuncts[CurrentElementNeigh]; 
 
        // # of basis functions, is the same as N_loc_dof
        N_Neigh = N_BaseFunct[CurrentElementNeigh]; 
     
        // find the part of the global index array where the information
        // for the neigh mesh cell are stored: DOF
        DOFNeigh = GlobalNumbers + BeginIndex[neigh_i];

        // initialize auxiliary array UsedNeigh
        memset(UsedNeigh, 0, N_FEs3D*SizeOfInt);
        
        // get id of finite element space in neigh cell 
        CurrentElementNeigh = fespace->GetFE3D(neigh_i, neigh);
        // set entry in auxiliary array to 1
        UsedNeigh[CurrentElementNeigh] = 1;
       
        // have a look how much 1 are in UsedNeigh -> N_LocalUsedElementsNeigh
        N_LocalUsedElementsNeigh = 0;
        memset(LocalUsedElementsNeigh, 0, SizeOfInt*N_FEs3D);
        j = 0;
        for(k=0;k<N_FEs3D;k++)
          if(UsedNeigh[k])
          {
            LocalUsedElementsNeigh[j] = (FE3D)k;
            j++;
          }
        N_LocalUsedElementsNeigh = j;
        
        // ####################################################################
        // calculate values on original element
        // ####################################################################
        
        // From id of fe space the situation in the referenz mesh cell is 
        // known. Now, get information in the original mesh cell. 
        // Every finite element gets automatically a quadrature rule.
        // input: N_LocalUsedElementsNeigh, LocalUsedElementsNeigh, 
        //        neigh, SecondDer
        // output : N_PointsNeigh - number of quadrature points
        //          xiNeigh - xi-values of quad. points in ref. cell
        //          etaNeigh - eta-values of quad. points in ref. cell
        //          weigthsNeigh - weights fro quad. rule in quad. points 
        //          XNeigh - X-values of quad. points in original cell  
        //          YNeigh - Y-values of quad. points in original cell  
        //          AbsDetjkNeigh - absolute value of Jacobian in quad. points
        //          xiNeigh, etaNeigh, weightsNeigh, XNeigh, YNeigh,  
        //          AbsDetjkNeigh are pointers to arrays

        TFEDatabase3D::GetOrig(N_LocalUsedElementsNeigh, LocalUsedElementsNeigh, 
                               Coll, neigh, SecondDer,
                               N_PointsNeigh, xiNeigh, etaNeigh, zetaNeigh,
                               weightsNeigh,
                               XNeigh, YNeigh, ZNeigh, AbsDetjkNeigh);
    
        // copy the values of the finite element function from the global
        // array (Values) to local arrays 
        for(lll=0;lll<N_Neigh;lll++)
        {
          FEFunctValuesNeigh[lll] = Values[DOFNeigh[lll]];
        }
        
        // compute values of fe function in the quadrature points with
        // linear combinations of the local fe values and the basis functions
        for(k=0;k<N_Derivatives;k++)
        {
          // get fe values of the basis functions in the quad. points
          // in the current element (pointer to pointer)
          OrigFEValuesNeigh = TFEDatabase3D::GetOrigElementValues(BaseFunctNeigh,
                                                       NeededDerivatives[k]);
          // accumulate the values of the fe function 
          // do it for all quad. points
          for(j=0;j<N_PointsNeigh;j++)
          {
            // get pointer to array where the values for quad. point j are
            // stored
            Orig = OrigFEValuesNeigh[j];
            // initialize 
            value = 0;
            // compute linear combination over all local dof
            for(lll=0;lll<N_Neigh;lll++)
            {
              // value += u1[l] * phi[l]
              value += FEFunctValuesNeigh[lll] * Orig[lll];
            } // endfor lll
            // copy result to position in array 
            FEValue_u1_QuadPointNeigh[j] = value;
          } // endfor j
        } // endfor k
        
        // compute part of the convolution in all local fe nodes (dof)
        // which comes from current mesh cell 
        // loop over all local dof
        
        // initialize
        value = 0;
        // loop over all quadrature points 
        for(j=0;j<N_PointsNeigh;j++) 
        { 
          // square of the distance between local dof and quad. point
          // in original cell 
          distance_sq = ((X_orig[l]-XNeigh[j])*(X_orig[l]-XNeigh[j])+
                         (Y_orig[l]-YNeigh[j])*(Y_orig[l]-YNeigh[j])+
                         (Z_orig[l]-ZNeigh[j])*(Z_orig[l]-ZNeigh[j])); 
          //cout << "X " << XNeigh[j] << " Y " << YNeigh[j] << endl;
          // compute value of GaussianFilter
          // input: delta, distance_sq
          g = GaussianFilter3D(delta,distance_sq);
          //cout << "dist " <<  sqrt(distance_sq) << endl;
          //
          // apply quadrature rule in ref. cell 
          value += g*AbsDetjkNeigh[j]*weightsNeigh[j]
                        *FEValue_u1_QuadPointNeigh[j];            
        }
        // add local result to global array
        u_conv[index] += value;
      } // end of ll loop (edges)
    } // end of l loop (local dof)
  } // endfor i (cells)
  
 
  // for (i=0;i<N_DOFConv; i++)
  //{
  // cout << i<< " x " <<   x_conv[i] << " y  " << y_conv[i] << " z  " << z_conv[i] ;
  // cout << " u " << Values[i] << " u_conv " << u_conv[i] << endl;
  //}

  // release memory
  delete SecondDer;
  delete x_conv;
  delete y_conv;
  delete z_conv;
  TDatabase::ParamDB->INTERNAL_QUAD_RULE=0;
} // ConvolutePressure

#endif
