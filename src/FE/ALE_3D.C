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
   
#include <ALE_3D.h>

void SolveGridEquation(double **Entries, double *sol, double *rhs,
                       int *KCol, int *RowPtr, int N_DOF)
{
  int i,j,k, col, Diognal;
  double *Entries11, *Entries12, *Entries13,
          *Entries21,*Entries22, *Entries23,
	  *Entries31,*Entries32, *Entries33;
	  
  double sum1, sum2, sum3; 
  int start, end;

  double max_error, error=1.e-12;
  int iter;

  Entries11 = Entries[0];
  Entries12 = Entries[1];
  Entries13 = Entries[2];
  Entries21 = Entries[3];
  Entries22 = Entries[4];
  Entries23 = Entries[5];
  Entries31 = Entries[6];
  Entries32 = Entries[7];
  Entries33 = Entries[8];

  max_error = 1.; iter = 0;
  while(max_error>error)
  {
    max_error = 0.0; iter++;
    for(i=0;i<N_DOF;i++)
     {
      start = RowPtr[i];
      end = RowPtr[i+1];
      sum1 = rhs[i];
      sum2 = rhs[i+N_DOF];
      sum3 = rhs[i+2*N_DOF];
      
      for(k=start;k<end;k++)
       {
        col = KCol[k];
        if (col==i) Diognal = k;
        sum1 -= Entries11[k]  * sol[col]
                +Entries12[k] * sol[col+N_DOF]
                +Entries13[k] * sol[col+2*N_DOF];
    
        sum2 -= Entries21[k]  * sol[col]
                +Entries22[k] * sol[col+N_DOF]
                +Entries23[k] * sol[col+2*N_DOF];
   
        sum3 -= Entries31[k]  * sol[col]
                +Entries32[k] * sol[col+N_DOF]
                +Entries33[k] * sol[col+2*N_DOF]; 
       } // endfor k

       sol[i] += sum1/Entries11[Diognal];
       sol[i+N_DOF] += sum2/Entries22[Diognal];
       sol[i+2*N_DOF] += sum3/Entries33[Diognal];
       
      if(max_error<fabs(sum1/Entries11[Diognal])) max_error = fabs(sum1/Entries11[Diognal]);
      if(max_error<fabs(sum2/Entries22[Diognal])) max_error = fabs(sum2/Entries22[Diognal]);
      if(max_error<fabs(sum3/Entries33[Diognal])) max_error = fabs(sum3/Entries33[Diognal]);
      
//     OutPut("r1 "<<fabs(sum1/Entries11[Diognal])<< " r2 " << fabs(sum2/Entries22[Diognal]) << " r3 " << fabs(sum3/Entries33[Diognal]) << endl); 
     } // endfor i
   if(iter == 1000) break;
       cout << " max_error " << max_error <<endl;
//    exit(0);
  } // end while
OutPut("Grid Solver: Number iteration "<<iter<< " max error" << max_error << endl);
exit(0);

}


// ======================================================================
// Standard Galerkin
// ======================================================================
void GridAssemble4(double Mult, double *coeff,
                  double *param, double hK,
                  double **OrigValues, int *N_BaseFuncts,
                  double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12,**MatrixA13,
         **MatrixA21, **MatrixA22,**MatrixA23,
	 **MatrixA31, **MatrixA32,**MatrixA33;
  double val, lambda=0;
  double *MatrixRow11, *MatrixRow12,*MatrixRow13,
         *MatrixRow21, *MatrixRow22,*MatrixRow23,
	 *MatrixRow31, *MatrixRow32,*MatrixRow33;
  double ansatz100, ansatz010,ansatz001, ansatz000;
  double test100, test010,test001, test000;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,k,l, N_U;
  double c0=1., detjk;

//   double lame = TDatabase::ParamDB->LameC;

//   if(TDatabase::ParamDB->P0)
   { detjk = coeff[19]; } // see DiscreteForm3D.C    
//   else
//    { detjk = 1; }

//   varying stiffness and divergence term values
//   lame *= detjk;
  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8]; 
  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z

  for(i=0;i<N_U;i++)
  {
    MatrixRow11 = MatrixA11[i];
    MatrixRow12 = MatrixA12[i];
    MatrixRow13 = MatrixA13[i];
    MatrixRow21 = MatrixA21[i];
    MatrixRow22 = MatrixA22[i];
    MatrixRow23 = MatrixA23[i];
    MatrixRow31 = MatrixA31[i];
    MatrixRow32 = MatrixA32[i];
    MatrixRow33 = MatrixA33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val  = c0*(2*test100*ansatz100+test010*ansatz010+test001*ansatz001) +lambda*test100*ansatz100;
      MatrixRow11[j] += Mult * val/detjk;   
      val  = c0*test010*ansatz100 + lambda*test100*ansatz010;
      MatrixRow12[j] += Mult * val/detjk;  
      val  = c0*test001*ansatz100 + lambda*test100*ansatz001;
      MatrixRow13[j] += Mult * val/detjk;

      val  = c0*test100*ansatz010 + lambda*test010*ansatz100;
      MatrixRow21[j] += Mult * val/detjk;  
      val  = c0*(test100*ansatz100+2*test010*ansatz010+test001*ansatz001) + lambda*test010*ansatz010;
      MatrixRow22[j] += Mult * val/detjk; 
      val  = c0*test001*ansatz010 + lambda*test010*ansatz001;
      MatrixRow23[j] += Mult * val/detjk;
   
      val  = c0*test100*ansatz001 + lambda*test001*ansatz100;
      MatrixRow31[j] += Mult * val/detjk;     
      val  = c0*test010*ansatz001 + lambda*test001*ansatz010;
      MatrixRow32[j] += Mult * val/detjk;  
      val  = c0*(test100*ansatz100+test010*ansatz010+2*test001*ansatz001) + lambda*test001*ansatz001;
      MatrixRow33[j] += Mult * val/detjk;
      
/*      // Laplace
      val  = c0*(test100*ansatz100+test010*ansatz010+test001*ansatz001) +lambda*test100*ansatz100;
      MatrixRow11[j] += Mult * val/detjk;   
      val  =  lambda*test100*ansatz010;
      MatrixRow12[j] += Mult * val/detjk;  
      val  =   lambda*test100*ansatz001;
      MatrixRow13[j] += Mult * val/detjk;

      val  =  lambda*test010*ansatz100;
      MatrixRow21[j] += Mult * val/detjk;  
      val  = c0*(test100*ansatz100+test010*ansatz010+test001*ansatz001) + lambda*test010*ansatz010;
      MatrixRow22[j] += Mult * val/detjk; 
      val  = c0*test001*ansatz010 + lambda*test010*ansatz001;
      MatrixRow23[j] += Mult * val/detjk;
   
      val  =  lambda*test001*ansatz100;
      MatrixRow31[j] += Mult * val/detjk;     
      val  =  lambda*test001*ansatz010;
      MatrixRow32[j] += Mult * val/detjk;  
      val  = c0*(test100*ansatz100+test010*ansatz010+test001*ansatz001) + lambda*test001*ansatz001;
      MatrixRow33[j] += Mult * val/detjk;  */    
      

    } // endfor j
  } // endfor i
 }
