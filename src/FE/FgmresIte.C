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
// @(#)FgmresIte.C        1.24 06/27/00
//
// Class:       TFgmresIte
// Purpose:     flexible gmres
//
// Author:      Volker John 24.10.2000
//
// History:     24.10.2000 start of implementation
//
// =======================================================================
#include <ItMethod.h>
#include <FgmresIte.h>
#include <MooNMD_Io.h>
#include <Database.h>
#include <LinAlg.h>
#include <NSE_MultiGrid.h>

#ifdef _MPI 
#include <ParFECommunicator3D.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef __3D__
    
#ifndef GEO_DIM    
#define GEO_DIM = 3;
#endif
    
#endif
    
#ifdef __2D__

#ifndef GEO_DIM    
#define GEO_DIM = 2;
#endif
    
#endif


TFgmresIte::TFgmresIte(MatVecProc *MatVec,
			DefectProc *Defect,
			TItMethod *Prec,
			int n_aux, int n_dof,
			int scalar_
#ifdef _MPI
  	#ifdef  __3D__
			       ,TParFECommunicator3D *parcomm_U, TParFECommunicator3D *parcomm_P
	#else		       
				,TParFECommunicator2D *parcomm_U, TParFECommunicator2D *parcomm_P
	#endif 
#endif
                      )
: TItMethod(MatVec, Defect, Prec, n_aux, n_dof)
{
  int i;
  double *aux;

  matvec = MatVec;
  matvecdefect = Defect;
  prec = Prec;
  N_DOF = n_dof;

  
  scalar = scalar_;
  if (scalar_)
  {
    res_norm_min = TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR;
    red_factor= TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR;
    maxit = TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR;
  }
  else
  {
    res_norm_min = TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SADDLE;
    red_factor= TDatabase::ParamDB->SC_LIN_RED_FACTOR_SADDLE;
    maxit = TDatabase::ParamDB->SC_LIN_MAXIT_SADDLE;
  }

  prec_maxit = TDatabase::ParamDB->SC_AMG_PREC_IT;
  div_factor = TDatabase::ParamDB->SC_DIV_FACTOR;
  minit = TDatabase::ParamDB->SC_MINIT;
  restart = TDatabase::ParamDB->SC_GMRES_RESTART;

  if (restart > maxit)
    restart = maxit;
  if (restart <=0 )
  {
    OutPut("WARNING: restart too small: " << restart);
    OutPut("         Set restart to 10 !!!");
    restart = 10;
  }

  if (prec_maxit<1)
  {
    OutPut("WARNING: Number of preconditioner iterations too small: "
      << prec_maxit << endl);
    OutPut("         Set number of preconditioner iterations to 1 !!!");
    prec_maxit = TDatabase::ParamDB->SC_AMG_PREC_IT = 1;
  }

  s = new double[restart+1];
  cosi = new double[restart+1];
  sn =  new double[restart+1];

  H = new double*[restart+1];
  H[0] = new double[(restart+1)*restart];
  memset(H[0],0,((restart+1)*restart)*SizeOfDouble);
  for(i=1;i<=restart;i++)
    H[i] = H[0] + i*restart;

  v = new double*[restart+1];
  v[0] = new double[(restart+1)*N_DOF];
  memset(v[0],0,((restart+1)*N_DOF)*SizeOfDouble);
  for(i=1;i<=restart;i++)
    v[i] = v[0] + i*N_DOF;

  zv = new double*[restart+1];
  zv[0] = new double[(restart+1)*N_DOF];
  memset(zv[0],0,((restart+1)*N_DOF)*SizeOfDouble);
  for(i=1;i<=restart;i++)
    zv[i] = zv[0] + i*N_DOF;

  defect =  new double[N_DOF];

  N_Aux = n_aux;
  if (n_aux>0)
  {
    AuxArray = new double* [n_aux];
    aux = new double[n_aux*N_DOF];
    for(i=0;i<n_aux;i++)
      AuxArray[i] = aux+i*N_DOF;
  }
  
#ifdef _MPI
  if(scalar)
   ParComm = parcomm_U;
  else{
    ParCommU = parcomm_U;
    ParCommP = parcomm_P;    
  }
#endif
}

TFgmresIte::~TFgmresIte()
{
  delete s;
  delete cosi;
  delete sn;
  delete H[0];
  delete H;
  delete v[0];
  delete v;
  delete zv[0];
  delete zv;
  delete defect;

  if (N_Aux>0)
  {
    delete AuxArray[0];
    delete AuxArray;
  }
}


void GeneratePlaneRotation(double dx, double dy, double *cs, double *sn)
{
  if (dy == 0.0)
  {
    cs[0] = 1.0;
    sn[0] = 0.0;
  }
  else if (fabs(dy) > fabs(dx))
  {
    double temp = dx / dy;
    sn[0] = 1.0 / sqrt( 1.0 + temp*temp );
    cs[0] = temp * sn[0];
  }
  else
  {
    double temp = dy / dx;
    cs[0] = 1.0 / sqrt( 1.0 + temp*temp );
    sn[0] = temp * cs[0];
  }
}

// rotate a vector by theta (cs = cos(theta), sn= sin(theta)
void ApplyPlaneRotation(double *dx, double *dy, double cs, double sn)
{
  double temp  =  cs * dx[0] + sn * dy[0];
  dy[0] = -sn * dx[0] + cs * dy[0];
  dx[0] = temp;
}


void UpdateGmresIterate(double *x, int Len_x, int k, double **h,
double *s, double **v)
{
  int i,j;
  /* Backsolve: */

  for (i = k; i >= 0; i--){
    s[i] /= h[i][i];
    for (j = i - 1; j >= 0; j--)
      s[j] -= h[j][i] * s[i];
  }

  for (j = 0; j <= k; j++)
    for(i=0; i<Len_x; i++)
      x[i] += v[j][i] * s[j]; // ======== COMPUTE z*Y
      
      // =========== COMMUPDATE FOR X REQUIRED?
      
}


int TFgmresIte::Iterate (TSquareMatrix **sqmat,
TMatrix **mat, double *sol,
double *rhs)
{
  int i=0, j,k,l, verbose = TDatabase::ParamDB->SC_VERBOSE;
  int maxite;
  double res, res0, reslast, t1, t2, temp,tempGlobal;
  double beta,residlast,dnorm0,end_residual;
  double eps = 1e-14;
  
  int flexible = TDatabase::ParamDB->SC_FLEXIBLE_KRYLOV_SPACE_SOLVER; 
  
// ========== 	check flexibility parameter   ========================
   // ======================== decide between parallel GMRES and sequential ===========
   // =================================================================
  
   
  if(flexible)
  {
    t1 = GetTime();

#ifdef _MPI
  TDatabase::ParamDB->time_GMRES_start = MPI_Wtime();
#else
  TDatabase::ParamDB->time_GMRES_start = GetTime();
#endif
      
#ifdef _MPI  


    
      int ii, rank, *MasterOfDof, *MasterOfDofU,*MasterOfDofP;  
      double  resglobal=0.0;
      int N_U = ParCommU->GetNDof();
      int N_P = ParCommP->GetNDof();
      
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);

#endif 

    if (verbose>1
    #ifdef _MPI
    && rank==0
    #endif
    )   
      OutPut("Entering fgmres" << endl);
    
    if ((TDatabase::ParamDB->INTERNAL_GMRES_INFO)&&(verbose==-1)){
      OutPut("fgmres maxite " << maxit << " restart " << restart << endl);
      TDatabase::ParamDB->INTERNAL_GMRES_INFO = 0;
    }

// ======== COMPUTE RESIDUAL ========

    matvecdefect(sqmat, mat, sol, rhs, v[0]); // Assume sol, rhs are updated 
    update(v[0]);
    // no need to update v for beta

// ================== COMPUTE BETA VALUE ======================    

    
    res=(Dotprod(N_DOF,v[0],v[0]));    
    res=res0=reslast=sqrt(res);     
      
                                 /* norm of residual */
    beta=res;

// =============== END OF COMPUTE BETA 

    maxite = maxit; // CONTROL PARAMETER

    // ===========	CHECK IF SOLUTION IS OBTAINED ARLEADY
    
    if (beta <= res_norm_min ){
      if ((minit==0)||(beta<1e-20)){
        if (verbose>0)
          OutPut("(no) fgmres iteration " << 0 << " " << beta << endl);
        end_residual=res;
        return(0);
      }
      else
        maxite = minit;      // IF RESIDAUL IS REACHED AS PER GMRES DO ATLEAST MINIMUM ITERATIONS
    }
    
    // ===============		GMRES ITERATION STARTS	======================================
    //========================================================================================
    //====================================================================================================
    
    if (verbose>0
#ifdef _MPI
    && rank==0
#endif
    )
      OutPut("fgmres Iteration " << 0 << " " << beta << endl);

    j=1;
    // iterate
    while (j<=maxite)
    {
      Dscal(N_DOF, 1.0/beta, v[0]); // v = r/beta
      memset(s,0,(restart+1)*SizeOfDouble);
      s[0]=beta; // store before a restart

// =================  iterate until converged or restart
      for(i=0; i<restart && j<=maxite; i++, j++){
	
        if (prec==NULL)
          memcpy(zv[i],v[i],N_DOF*SizeOfDouble); // if no preconditioner then M = I (identity matrix)
        
        
// --------THIS IS TO OBTAIN THE PRECONDITIONED RESIDUAL AND THE SOLUTION---------------------------------------        
        
        else
	{ 
	      
	      memset(zv[i],0,N_DOF*SizeOfDouble);
	  
// 	      if (prec_maxit>1)
// 		dnorm0 = sqrt(Ddot(N_DOF,v[i],v[i]))*TDatabase::ParamDB->SC_AMG_PREC_RED_FACTOR;
	      
	      for (k=0;k<prec_maxit;k++)
	      {								
	
		
	      #ifdef _MPI
	      TDatabase::ParamDB->time_GMRES_end = MPI_Wtime();
	      #else
	      TDatabase::ParamDB->time_GMRES_end = GetTime();
	      #endif

	      TDatabase::ParamDB->time_GMRES += TDatabase::ParamDB->time_GMRES_end - TDatabase::ParamDB->time_GMRES_start; 
	      
              prec->Iterate(sqmat, mat, zv[i], v[i]); // Send the updated zv into the routine	  
		
		
	      #ifdef _MPI
	      TDatabase::ParamDB->time_GMRES_start = MPI_Wtime();
	      #else
	      TDatabase::ParamDB->time_GMRES_start = GetTime();
	      #endif
					
	     if (k<prec_maxit-1)
	     {
		
		
              matvecdefect(sqmat, mat, zv[i], v[i], defect); // ------- CHECK IF THE RESIUDAL HAS DECREASED THE REQUIRED AMOUNT

              res=(Dotprod(N_DOF,defect,defect));    
              res=res0=reslast=sqrt(res);                   
  
		
//===================================DEBUG MESSAGES =================================
//=========================					===========================

cout << "==========	After "<< k+1 <<" multigrid sweep "<< res <<"	==================="<< endl; 	      
	     
//               if (res <= dnorm0)
//                 break;
	      }
	    }
        }

// --------	END OF PRECONDITIONING ---------------------------------------        
                                 
        // ------------ GET THE W VECTOR
        matvec(sqmat, mat, zv[i], defect); // send back an updated defect
	update(defect);
        // -------- NO NEED OF COMM UPDATE AS ZV IS ALREADY UPDATED ACROSS	
 
// --------- COMPUTE THE ENTRIES OF THE HESSENBERG MATRIX
	
        for(k=0; k<=i; k++)
	{ 
	  
          H[k][i]= Dotprod(N_DOF, defect,v[k]);    
	  Daxpy(N_DOF, -H[k][i], v[k], defect); // ------ COMPUTE NEW W VECTOR
	  // defect and v[k] already updated so need of further communication
	}
	
	
// --------- END OF COMPUTE THE ENTRIES OF THE HESSENBERG MATRIX


// --------------- COMPUTE THE NORM OF - W VECTOR FOR H_i+1_i

        H[i+1][i]=sqrt(Dotprod(N_DOF,defect,defect));
           
        Dcopy(N_DOF, defect, v[i+1]); 
	
        if(H[i+1][i]!=0)
          Dscal(N_DOF, 1/H[i+1][i], v[i+1]);
	// ------------------- obtain new value of V
        
	else
	{
          Error("Fehler!" << endl);
          OutPut("Fehler in here actually FgmresIte.C !!!" << endl);
        }

// --------------- COMPUTE THE NORM OF - W VECTOR FOR H_i+1_i



// ========================================================================================
// =================== OBTAINED H MATRIX, BETA. SAME ACROSS ALL PROCESSES
// ================================ COMPUTE Y VECTOR TO UPDATE SOLUTION
// ======================================


        for(k=0;k<i;k++)
          ApplyPlaneRotation(&H[k][i], &H[k+1][i], cosi[k], sn[k]);

        GeneratePlaneRotation(H[i][i], H[i+1][i], &cosi[i], &sn[i]);
        ApplyPlaneRotation(&H[i][i], &H[i+1][i], cosi[i], sn[i]);
        ApplyPlaneRotation(&s[i], &s[i+1], cosi[i], sn[i]);

        residlast=res;  // --------------- CHECK WHY IS THIS RESIDUAL BEING CHECKED UPON?

// ========= START OF STOPPING CRITERIA
    // ========= ROUTINE TERMINATES IF WE GET INTO THIS LOOP
	
        if ((j >= minit) && ((((res=fabs(s[i+1])) < res_norm_min)|| ((res=fabs(s[i+1]))<red_factor*res0)))){
          
	  memset(defect,0,N_DOF*SizeOfDouble);
          
	  UpdateGmresIterate(defect /* The w vector*/, 
			     N_DOF, 
			     i/*to mention limit of acess to H*/,
			     H /*the hessenberg matrix*/, 
			     s /*beta value*/, 
			     zv/*Required to compute final solution update*/);
                                 /* new iterate */
          
	  Daxpy(N_DOF,1.0,defect,sol); // Defect = ZmYm in previous step,
          // Add sol = sol + update
          // COMMUPDATE REQURIED? 
          
    // ======= SOLUTION IS UPDATED HERE
	  
	  t2 = GetTime();
          if (verbose>=1
#ifdef _MPI
          && rank==0
#endif
          ){
            OutPut("FGMRES ITE: " << setw(4) << j);
            OutPut(" t: " << setw(6) << t2-t1);
            OutPut(" t/cyc: " << setw(6) << (t2-t1)/j);
            OutPut(" res : "  << setw(8) << res);
            OutPut(" rate: " << pow(res/res0,1.0/j) << endl);
          }
          if (verbose==-2
#ifdef _MPI
            && rank==0
#endif
          )
          {
            OutPut("vanka fgmres ite: " << setw(4) << j);
            OutPut(" t: " << setw(6) << t2-t1);
            OutPut(" t/cyc: " << setw(6) << (t2-t1)/j);
            OutPut(" res : "  << setw(8) << res);
            OutPut(" rate: " << pow(res/res0,1.0/j) << endl);

          }
         
	  
#ifdef _MPI
  TDatabase::ParamDB->time_GMRES_end = MPI_Wtime();
#else
  TDatabase::ParamDB->time_GMRES_end = GetTime();
#endif

TDatabase::ParamDB->time_GMRES += TDatabase::ParamDB->time_GMRES_end - TDatabase::ParamDB->time_GMRES_start; 

      return(j); // RETURN THE NUMBER OF ITERATIONS THAT TOOK PLACE
	  
        }
        
// ======================== END OF STOPPING CRITERIA CASE                        
        // ===================== SOME CHECKS WRT RESIDUAL =============== 
        
        if (verbose>0 //&& //j%5==0
#ifdef _MPI
          && rank==0
#endif
        ){
          OutPut("fgmres iteration " << j << " " << res << " " << res/residlast);
          OutPut(endl);
        }
        
        if (res>div_factor*res0){
          OutPut("fgmres iteration diverges !!!" << endl);
          exit(4711);
        }

      }                            /* endfor i */

// =============================== CALCULATED UNTIL J = RESTART
// ========================= EVALUATE NEW SOL. , UPDATE V , BETA TO RESTART THE ROUTINE

      memset(v[0],0, N_DOF*SizeOfDouble);
      UpdateGmresIterate(v[0], N_DOF, restart-1, H, s, zv);
      
      Daxpy(N_DOF,1.0,v[0],sol);   

      // ============= SOL. IS UPDATED HERE
      // ------------- REUQIRE COMM UPDATE ?
      
      matvecdefect(sqmat, mat, sol, rhs, v[0]); // recieve an updated v[0]
      update(v[0]);
      
      // ========== UPDATED V FOR NEXT ITERATION
      
// --------------------- UPDATE BETA FOR NEXT ITERATION
      
      res=(Dotprod(N_DOF,v[0],v[0]));    
      res=res0=reslast=sqrt(res);     
     
      beta=res;

// ===================== BETA UPDATED FOR NEXT ITERATION 

      // -------- IS THIS CHECK REQUIRED ?
      if (verbose>1){
        if(fabs(beta-res)>0.01*beta){
          OutPut("WARNING: restart residual changed " <<  beta << " " << res);
          OutPut(endl);
        }
      }
      
    } 
// =========================== COMPLETED WHILE - ITERATIONS :: GMRES HAS TRIED ITS BEST
    
    t2 = GetTime();
 
    // ============= OUTPUTS BASED ON THE RESIDUAL DETERMINED ===================
    
    if (j>=maxite && res>res_norm_min &&res>res0*red_factor ){
      if (verbose>-1)
        OutPut("FGMRES not converged !!!" << endl);
    }

    if (j>maxite)
      j=maxite;
    
    if (verbose>-1){
      OutPut("FGMRES ITE: " << setw(4) << j);
      OutPut(" t: " << setw(6) << t2-t1);
      OutPut(" t/cyc: " << setw(6) << (t2-t1)/j);
      OutPut(" res : "  << setw(8) << res);
      OutPut(" rate: " << pow(res/res0,1.0/j) << endl);
    }
    if (verbose==-2){
      OutPut("vanka fgmres ite: " << setw(4) << j);
      OutPut(" t: " << setw(6) << t2-t1);
      OutPut(" t/cyc: " << setw(6) << (t2-t1)/j);
      OutPut(" res : "  << setw(8) << res);
      OutPut(" rate: " << pow(res/res0,1.0/j) << endl);
    }
    
    
    // ===============		GMRES ITERATION ENDS	======================================
    //========================================================================================
    //====================================================================================================

    
#ifdef _MPI
  TDatabase::ParamDB->time_GMRES_end = MPI_Wtime();
#else
  TDatabase::ParamDB->time_GMRES_end = GetTime();
#endif

TDatabase::ParamDB->time_GMRES += TDatabase::ParamDB->time_GMRES_end - TDatabase::ParamDB->time_GMRES_start; 
 
        return(j);
    // =============================== END BY RETURNING THE NUMBER OF ITERATIONS =============== 
  }
  
// ============================  END OF FLEXIBLE ROUTINE =================================


  else
  {
    double *r = new double[N_DOF];
    if (verbose>1)
      OutPut("Entering GMRES" <<endl);
    
    Dcopy(N_DOF, rhs, r); // IS THIS REQUIRED ??

    // ======= GET RESIDUAL IN r
    
    matvecdefect(sqmat, mat, sol, rhs, r);

    // ======= EVALUATE BETA
    beta=sqrt(Ddot(N_DOF, r,r));

    // =========== CONTROL CHECKS - CHECK IF RESIDUAL ALREADY REACHED AND PERFORM MINIMUM ITERATIONS
    double resid = beta;
    double start_residual=resid;
    if ((beta  <= res_norm_min )&&(minit==0)){
      if (verbose>0){
        OutPut("GMRES Iteration 0 " << beta << endl);
        end_residual=resid;
      }
      return(0);
    }
    
    if (verbose>0)
      OutPut("GMRES Iteration 0 " << beta << endl);
    if (verbose>1)
      OutPut("Entering GMRES iteration cycle" << endl);

    
    //===================================	GMRES LOOP BEGINS	====================================//
    //========================================================================================================
    //=========================================================================================================
    
   
    
    
    j=1;
    while (j<=maxit)
    {
 
      // ============= CALCULATE V[0]
      Dcopy(N_DOF, r, v[0]);
      Dscal(N_DOF, 1.0/beta, v[0]);      
      
      memset(s,0,(restart+1)*SizeOfDouble);
      s[0]=beta;
      
      for (i=0;i<restart && j<=maxit;i++, j++){
        
	// ========= IF NO PRE-COONDITIONER THEN EVALUTE A*V[i]
	if (prec==NULL)
          matvec(sqmat, mat, v[i], defect);  //
        
	// =============== GET THE UPDATED PRE-CONDITIONED RESIDUAL INTO r  
	else{
          Dcopy(N_DOF, v[i], zv[0]);
          memset(r,0.0,N_DOF*SizeOfDouble);
	  
          for (k=0;k<prec_maxit;k++){
            prec->Iterate(sqmat, mat, r, v[i]);		// ======== GET Z into r
          }                    
          matvec(sqmat, mat, r, defect);	// =========== EVALUATE W VECTOR INTO defect
        }
        
        // ==============	GET ENTRIES INTO H-MATRIX, UPDATE DEFECT
        for(k=0; k<=i; k++){
          H[k][i]=Ddot(N_DOF, defect, v[k]);
          Daxpy(N_DOF, -H[k][i], v[k], defect);
        }
                
        H[i+1][i]=sqrt(Ddot(N_DOF, defect,defect)); // =========	H_i+1_i
	
        Dcopy(N_DOF, defect, v[i+1]); 
        Dscal(N_DOF, 1.0/H[i+1][i], v[i+1]);
        // ============ 	UPDATE V FOR NEXT ITERATION
	
	for(k=0;k<i;k++)	
          ApplyPlaneRotation(&H[k][i], &H[k+1][i], cosi[k], sn[k]);		
        GeneratePlaneRotation(H[i][i], H[i+1][i], &cosi[i], &sn[i]);
        ApplyPlaneRotation(&H[i][i], &H[i+1][i], cosi[i], sn[i]);
        ApplyPlaneRotation(&s[i], &s[i+1], cosi[i], sn[i]);
        
	residlast=resid;
	
// ========= START OF STOPPING CRITERIA
    // ========= ROUTINE TERMINATES IF WE GET INTO THIS LOOP
	
        if (((resid=fabs(s[i+1])) < res_norm_min) || ((resid=fabs(s[i+1]))<red_factor*start_residual)){
	  
          memset(zv[0],0,((restart+1)*N_DOF)*SizeOfDouble);
	  
          UpdateGmresIterate(zv[0],N_DOF, i, H, s, v);
	  
	  // ===========	GET UPDATE INTO ZV[0]
          
	  if (prec==NULL) // ===============	IF NO PRE-CONDITIONER THEN PROCEED
            Daxpy(N_DOF, 1.0, zv[0], sol);
	  
          else{		// ==================		NOT SURE ??
            memset(r,0.0,N_DOF*SizeOfDouble);
            for (k=0;k<prec_maxit;k++)
              prec->Iterate(sqmat, mat, r, zv[0]);
            
            Daxpy(N_DOF, 1.0, r, sol);            /* new iterate */
          }
          
          if (verbose>0)
            OutPut("GMRES (right) : iterations " << j << " residual " <<resid << endl);
          end_residual=resid;

          return(j);
        }
        
  // ========= END OF STOPPING CRITERIA      
        
        if (verbose>0)
          OutPut("GMRES (right) Iteration " << j << " residuum " << resid << " " << resid/residlast << endl);
	
        if (resid>div_factor*start_residual){
          OutPut("GMRES (right) iteration diverges !!!" << endl);
          exit(4711);
        }
 
      }      // ============================	END OF GMRES RESTART LOOP
      
      memset(zv[0],0,N_DOF*SizeOfDouble); // =============	RESET ZV

      
      
      UpdateGmresIterate(zv[0], N_DOF, restart-1, H, s, v);
      
      if (prec==NULL)      
        Dcopy(N_DOF, zv[0],r);

      else{
        memset(r,0,N_DOF*SizeOfDouble);
        for (l=0;l<prec_maxit;l++)        
          prec->Iterate(sqmat, mat, r, zv[0]);
      }
      
      Daxpy(N_DOF, 1.0, r, sol);	
      // ============================	SOLUTION IS UPDATED HERE
      
      // ============================================	UPDATE RESIDUAL, V -, BETA FOR START OF NEXT GMRES 
      
      Dcopy(N_DOF,rhs,r);                         /* copy rhs (b) into r */
      
      
      matvecdefect(sqmat, mat, sol, rhs, r);
      beta=sqrt(Ddot(N_DOF, r, r));               /* norm of residual */

      if (verbose>1){
        if(fabs(beta-resid)>0.01*beta)
          OutPut("restart residual changed " << beta << " " << resid << endl);
      }
    }
    
      //===================================	GMRES LOOP ENDS	====================================//
    //========================================================================================================
    //=========================================================================================================

    if (verbose>0){
      OutPut("GMRES (right) : (maximal) iterations " << maxit << " residual " << resid << endl);
    }
    end_residual=resid;
    	  
    return(j);
  }
}

double TFgmresIte::Dotprod(int x, double* v1, double* v2)
{
  int *MasterOfDof,*MasterOfDofP,*MasterOfDofU;
  double sum=0.0;
  double glo_sum=0.0;
  int NDofU,NDofP;
#ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if(scalar)
    MasterOfDof = ParComm->GetMaster();
  else
  {
    MasterOfDofU = ParCommU->GetMaster();
    MasterOfDofP = ParCommP->GetMaster();
    NDofU = ParCommU->GetNDof();
    NDofP = ParCommP->GetNDof();
  }    
#endif

#ifdef _MPI
  if(scalar)
  {
    for(int ii=0;ii<N_DOF;ii++)
    {
      if(MasterOfDof[ii]==rank)
      {
	 for(int loc=0;loc<GEO_DIM;loc++)
	 {
	    sum+= v1[ii + loc*NDofU]*v2[ii + loc*NDofU];  
	 }
      }
    }
  }
  else
  {
    for(int ii=0;ii<NDofU;ii++)
    {
      if(MasterOfDofU[ii]==rank)
      {
	 for(int loc=0;loc<GEO_DIM;loc++)
	 {
	    sum+= v1[ii + loc*NDofU]*v2[ii + loc*NDofU];  
	 }
      }
    }
    
    for(int ii=0;ii<NDofP;ii++)
    {
      if(MasterOfDofP[ii]==rank)
      {
	    sum+= v1[ii + GEO_DIM*NDofU]*v2[ii + GEO_DIM*NDofU];  	 
      }
    }
  }
  
  TDatabase::ParamDB->time_communication_start =  MPI_Wtime();
  MPI_Allreduce(&sum, &glo_sum, 1, MPI_DOUBLE, MPI_SUM, TDatabase::ParamDB->Comm);  
  TDatabase::ParamDB->time_communication_end =  MPI_Wtime();
  TDatabase::ParamDB->time_communication += TDatabase::ParamDB->time_communication_end - 
					      TDatabase::ParamDB->time_communication_start;
#else

  sum = Ddot(N_DOF,v1,v2);
  glo_sum = sum;

#endif

  
  return glo_sum;
}


void TFgmresIte::update(double* &v1)
{
  int *MasterOfDof,*MasterOfDofP,*MasterOfDofU;
  double sum=0.0;
  double glo_sum=0.0;
  int NDofU,NDofP;
#ifdef _MPI
  if(scalar)
  {
    
    ParComm->CommUpdate(v1);
    
  }
  else
  {
    NDofU = ParCommU->GetNDof();
    NDofP = ParCommP->GetNDof();

    ParCommU->CommUpdate(v1);
    ParCommP->CommUpdate(v1 + GEO_DIM*NDofU);
 
  }   
#endif

}
