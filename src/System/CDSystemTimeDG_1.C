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
   
/** ************************************************************************ 
*
* @brief     dG(1) in time discretization
* @author    Sashikumaar Ganesan, 
* @date      24.12.15
* @History    
 ************************************************************************  */

#include <Database.h>
#include <string.h>
#include <CDSystemTimeDG_1.h>
#include <LinAlg.h>

TCDSystemTimeDG_1::TCDSystemTimeDG_1(TSquareMatrix2D *mat_m, TSquareMatrix2D *mat_A)
                  :TCDSystemTimeDG(2, 2, mat_m, mat_A)        
{
 int i, j, k, pos, end, len, disp1, disp2, *RowPtr, *KCol;  
 
  RowPtr = Mat_M->GetRowPtr(); 
  KCol  = Mat_M->GetKCol();  
  
   //first rwo block matrices
   pos =0;
   Sys_rowptr[0] = 0;
   for(i=0; i<N_U; i++)
    {
     len = RowPtr[i+1] - RowPtr[i];
     disp1 = Sys_rowptr[i];
     disp2 = disp1 + len;
     k = RowPtr[i];
     for(j=0; j<len; j++)
     {
      Sys_colindex[disp1+j] =  KCol[k];
      Sys_colindex[disp2+j] = N_U + KCol[k];
      pos +=2; 
      k++;
     } // for(j=RowPtr[i]; j<end; j++)  
    Sys_rowptr[i+1] = pos;
     
     /**testing */
//      for(j=Sys_rowptr[i]; j<Sys_rowptr[i+1]; j++)
//            cout << i << " Sys_colindex  "  << j << " " << Sys_colindex[j]<< endl;
//      cout<<endl;
    }// for(i=0; i<N_U; i++

    //second row block
    for(i=0; i<N_U; i++)
    {
     len = RowPtr[i+1] - RowPtr[i];
     disp1 = Sys_rowptr[N_U +i];
     disp2 = disp1+len;
     k = RowPtr[i];
     for(j=0; j<len; j++)
     {
      Sys_colindex[disp1+j] =  KCol[k];
      Sys_colindex[disp2+j] = N_U + KCol[k];
      k++;
      pos +=2;   
     } // for(j=RowPtr[i]; j<end; j++)
    Sys_rowptr[N_U + i+1] = pos;
    
     /**testing */
//      for(j=Sys_rowptr[N_U +i]; j<Sys_rowptr[N_U +i+1]; j++)
//            cout << i << " Sys_colindex  "  << j << " " << Sys_colindex[j]<< endl;
//      cout<<endl;    
    }// for(i=0; i<N_U; i++ 
  
  // dG structure is generated with increasing col entries, since M & A are in this order
  Sys_structure->SetColOrder(1);
  
}


TCDSystemTimeDG_1::~TCDSystemTimeDG_1()
{
  
}


void TCDSystemTimeDG_1::AssembleSysMat(double *Mu_old, double *Rhs)
{
  int i, j, k, len, *RowPtr, *KCol, N_Entries, N_Rhs, disp1, disp2, *Sys_colindex_aux;  
  double *Sys_Entries, *Sys_Entries_aux, *rhs_dG_aux, *M_Entries, *A_Entries, tau, *Rhs_loc;
  double T[2], m, tau_a;
  
//   N_U = Mat_M->GetN_Columns();  
  RowPtr = Mat_M->GetRowPtr(); 
  KCol  = Mat_M->GetKCol();  
  N_Entries = RowPtr[N_U];
  
  M_Entries = Mat_M->GetEntries();
  A_Entries = Mat_A->GetEntries();
  Sys_Entries = Sys_Mat->GetEntries();
  Sys_Entries_aux = Sys_Entries+(N_ColBlockMat*N_Entries); //second row block
  Sys_colindex_aux = Sys_colindex+(N_ColBlockMat*N_Entries); //second row block;  
  
  rhs_dG_aux = rhs_dG+N_U; // second rhs
  
  tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
     
  //Gauss Radau quarature
  // t_1 = t^{n-1} + dt/3, and t_2 = t^n (which is already available)
  // w_1 = 3/2,            and w_2 = 1/2
  N_Rhs = 1;
  T[0] = (TDatabase::TimeDB->CURRENTTIME - tau) + tau/3.; //t_1
  Rhs_loc = new double[N_Rhs*N_U];
  
  this->AssembleRhs(N_Rhs, T, Rhs_loc);
 
   for(i=0; i<N_U; i++)
    {
     len = RowPtr[i+1] - RowPtr[i];
     k = RowPtr[i];
     disp1 = Sys_rowptr[i];
     disp2 = disp1+len;     
      
    if(i<N_UActive)
     {         
      rhs_dG[i] = Mu_old[i] + 0.25*tau*(3.*Rhs_loc[i] + Rhs[i]);
      rhs_dG_aux[i] = 0.25*tau*(Rhs_loc[i] + Rhs[i]);

      for(j=0; j<len; j++)
       {
        tau_a = tau*A_Entries[k];
        m=M_Entries[k];
        k++;
      
        Sys_Entries[disp1+j] =  m +  tau_a;    // S11
        Sys_Entries[disp2+j] =   m + 0.5*tau_a;   // S12  
      
        Sys_Entries_aux[disp1+j] = 0.5*tau_a; // S21
        Sys_Entries_aux[disp2+j] =  0.5*m +  (1./3.)*tau_a;// S22  
       } // for(j=0; j<len; j++)
      } //  if(i<N_UActive)
     else //Dirichlet
      {      
       rhs_dG[i] =  Rhs[i];
       rhs_dG_aux[i] = 0.;

       for(j=0; j<len; j++)
        {
         if(Sys_colindex[disp1+j]==i)
          {
           Sys_Entries[disp1+j] =  1.;    // S11
           Sys_Entries[disp2+j] =  0.;    // S12        
          }
        else
         {
          Sys_Entries[disp1+j] =  0.;    // S11
          Sys_Entries[disp2+j] =  0.;    // S12   
         }
          
         if(Sys_colindex_aux[disp2+j]==(N_U+i) )
          {         
           Sys_Entries_aux[disp1+j] = 0.; // S21        
           Sys_Entries_aux[disp2+j] = 1.; // S22  
          }
         else
          {         
           Sys_Entries_aux[disp1+j] = 0.; // S21        
           Sys_Entries_aux[disp2+j] = 0.; // S22  
          }          
         }// for(j=0; j<len; j++)
       
      }// //Dirichlet 
    }// for(i=0; i<N_U; i++  
  
//       cout <<  " TCDSystemTimeDG_1 AssembleSysMat " << T[0]  << endl;
//       exit(0);
      
}

void TCDSystemTimeDG_1::AssembleALESysMat_Qp1(double *Mu_old, double *Rhs)
{
  int i, j, k, len, *RowPtr, *KCol, N_Entries, N_Rhs, disp1, disp2, *Sys_colindex_aux;  
  double *Sys_Entries, *Sys_Entries_aux, *rhs_dG_aux, *M_Entries, *A_Entries, tau;
  double *M_Qp1_Entries, *A_Qp1_Entries;
  double  m, m_Qp1, tau_a, tau_a_Qp1;
  
  RowPtr = Mat_M->GetRowPtr(); 
  KCol  = Mat_M->GetKCol();  
  N_Entries = RowPtr[N_U];
  
  M_Entries = Mat_M->GetEntries(); //t^n+1
  A_Entries = Mat_A->GetEntries(); //t^n+1
  
  if(QpMatricsAdded)
  {
   M_Qp1_Entries = Mat_M_Qp1->GetEntries(); //t^n+1/3
   A_Qp1_Entries = Mat_A_Qp1->GetEntries(); //t^n+1/3   
  }
  
  Sys_Entries = Sys_Mat->GetEntries();
  Sys_Entries_aux = Sys_Entries+(N_ColBlockMat*N_Entries); //second row block
  Sys_colindex_aux = Sys_colindex+(N_ColBlockMat*N_Entries); //second row block;  
  rhs_dG_aux = rhs_dG+N_U; // second rhs

  tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;

  //Gauss Radau quarature
  // t_1 = t^{n-1} + dt/3, and t_2 = t^n (which is already available)
  // w_1 = 3/2,            and w_2 = 1/2
  
   for(i=0; i<N_U; i++)
    {
     len = RowPtr[i+1] - RowPtr[i];
     k = RowPtr[i];
     disp1 = Sys_rowptr[i];
     disp2 = disp1+len;      
      
     if(i<N_UActive)
      {
       rhs_dG[i] = Mu_old[i] + 0.25*tau*(3.*Rhs_Qp1[i] + Rhs[i]);
       rhs_dG_aux[i] = 0.25*tau*(Rhs_Qp1[i] + Rhs[i]);
     
        for(j=0; j<len; j++)
         {

          tau_a = tau*A_Entries[k]; 
          m=M_Entries[k];

          tau_a_Qp1 = tau*A_Qp1_Entries[k];
          m_Qp1 = M_Qp1_Entries[k];

          k++;
    
          if(CONSERVATIVEALE)
           {
            Sys_Entries[disp1+j] =  m  +  0.75*tau_a_Qp1 + 0.25*tau_a;   // S11
            Sys_Entries[disp2+j] =   m +  0.25*tau_a_Qp1 + 0.25*tau_a;   // S12  
      
            Sys_Entries_aux[disp1+j] = -0.75*m_Qp1 + 0.75*m + 0.25*tau_a_Qp1     + 0.25*tau_a; // S21
            Sys_Entries_aux[disp2+j] = -0.25*m_Qp1 + 0.75*m + (1./12.)*tau_a_Qp1 + 0.25*tau_a; // S22  
           }
          else // NON-CONSERVATIVEALE
           {
            Sys_Entries[disp1+j] =                   m + 0.75*tau_a_Qp1 + 0.25*tau_a;   // S11
            Sys_Entries[disp2+j] = 0.75*m_Qp1 + 0.25*m + 0.25*(tau_a_Qp1 + tau_a);   // S12  
      
            Sys_Entries_aux[disp1+j] =                            0.25*(tau_a_Qp1 + tau_a); // S21
            Sys_Entries_aux[disp2+j] = 0.25*m_Qp1 + 0.25*m  + (1./12.)*tau_a_Qp1 + 0.25*tau_a; // S22  	     
//             cout << "  NON-CONSERVATIVEALE  AssembleALESysMat_Qp1 " <<endl; 
//             exit(0);
           }
          } // for(j=0; j<len; j++)
      } //  if(i<N_UActive)
     else //Dirichlet
      {
       rhs_dG[i] =  Rhs[i];
       rhs_dG_aux[i] = 0.;
  
        for(j=0; j<len; j++)
         {
          if(Sys_colindex[disp1+j]==i)
           {
            Sys_Entries[disp1+j] =  1.;    // S11
            Sys_Entries[disp2+j] =  0.;    // S12        
           }
          else
          {
           Sys_Entries[disp1+j] =  0.;    // S11
           Sys_Entries[disp2+j] =  0.;    // S12   
          }          
 
          if(Sys_colindex_aux[disp2+j]==(N_U+i) )
           {         
            Sys_Entries_aux[disp1+j] = 0.; // S21        
            Sys_Entries_aux[disp2+j] = 1.; // S22     
           }
          else
           {         
            Sys_Entries_aux[disp1+j] = 0.; // S21        
            Sys_Entries_aux[disp2+j] = 0.; // S22  
           } 
          
         }// for(j=0; j<len; j++)
      } // else if(i<N_UActive)
    }// for(i=0; i<N_U; i++  
  
//       cout <<  " TCDSystemTimeDG_1 AssembleALESysMat "    << endl;
//       exit(0);
      
} // AssembleALESysMat


void TCDSystemTimeDG_1::SoveTimedG(double *Sol)
{
 int i;

   //solve the dG system
   this->SovedGSystem();
  
   // u(t^n) = u(t_1) + u(t_2); in dG(1)
   for(i=0; i<N_U; i++)
    {
     Sol[i] = Sol_dG[i] + Sol_dG[N_U+i];
    }
    
//   cout <<  " TCDSystemTimeDG_1 SoveTimedG: " << Ddot(N_U, Sol, Sol)  << endl;
//   exit(0); 
}
