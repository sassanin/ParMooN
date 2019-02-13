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
   
/* 21 quadrature points + 4 vertex values */
static double NF_S_Q_Q2_2D_Xi[] = 
        { -0.77459666924148337703585307995647992, 0, 
           0.77459666924148337703585307995647992, 
           
	  1, 1, 1, 
          
	  0.77459666924148337703585307995647992, 0, 
          -0.77459666924148337703585307995647992, 
          
	  -1, -1, -1,
          
	  -0.77459666924148337703585307995647992, 0, 
           0.77459666924148337703585307995647992, 
          
	  -0.77459666924148337703585307995647992, 0, 
           0.77459666924148337703585307995647992, 
	  
	  -0.77459666924148337703585307995647992, 0, 
           0.77459666924148337703585307995647992,
	  
	  -1,1,1,-1 
	};

static double NF_S_Q_Q2_2D_Eta[] = 
        { -1, -1, -1,
 
          -0.77459666924148337703585307995647992, 0,
           0.77459666924148337703585307995647992,

           1, 1, 1,

           0.77459666924148337703585307995647992, 0,
          -0.77459666924148337703585307995647992, 

	   -0.77459666924148337703585307995647992,
	  -0.77459666924148337703585307995647992,
	  -0.77459666924148337703585307995647992,

	  0,0,0,

	  0.77459666924148337703585307995647992,
	  0.77459666924148337703585307995647992,
	  0.77459666924148337703585307995647992,

	  -1,-1,1,1
	};

/* not used */
static double NF_S_Q_Q2_2D_T[] = 
        { -1, -0.77459666924148337703585307995647992, 0,
           0.77459666924148337703585307995647992, 1 };

void NF_S_Q_Q2_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues, double *Functionals)
{
  double f0,f1,f2,f3,f4,f5,f6,f7,f8;
  static double weights[3] = { 0.5555555555555555555555555555555556,
                               0.88888888888888888888888888888888889,
                               0.5555555555555555555555555555555556 };
  /* vertex values */
  f0 = PointValues[21];
  f1 = PointValues[22];
  f2 = PointValues[23];
  f3 = PointValues[24];
   
  /* line integrals */
  f4 = ( weights[0]*PointValues[0]
                    +weights[1]*PointValues[1]
                    +weights[2]*PointValues[2]) * 0.5;
  f5 = ( weights[0]*PointValues[3]
                    +weights[1]*PointValues[4]
                    +weights[2]*PointValues[5]) * 0.5;
  f6 = ( weights[0]*PointValues[6]
                    +weights[1]*PointValues[7]
                    +weights[2]*PointValues[8]) * 0.5;
  f7 = ( weights[0]*PointValues[9]
                    +weights[1]*PointValues[10]
                    +weights[2]*PointValues[11]) * 0.5;
  
  /* cell integral */
  f8 = ( weights[0]*weights[0]*PointValues[12]
	 +weights[1]*weights[0]*PointValues[13]
	 +weights[2]*weights[0]*PointValues[14]
	 +weights[0]*weights[1]*PointValues[15]
	 +weights[1]*weights[1]*PointValues[16]
	 +weights[2]*weights[1]*PointValues[17]
	 +weights[0]*weights[2]*PointValues[18]
	 +weights[1]*weights[2]*PointValues[19]
	 +weights[2]*weights[2]*PointValues[20])*0.25;

  /* basis transformation */
  Functionals[0]=f0;
  Functionals[2]=f1;
  Functionals[8]=f2;
  Functionals[6]=f3;
  Functionals[1]=-0.25*(f0+f1)+1.5*f4;
  Functionals[5]=-0.25*(f1+f2)+1.5*f5;
  Functionals[7]=-0.25*(f2+f3)+1.5*f6;
  Functionals[3]=-0.25*(f0+f3)+1.5*f7;
  Functionals[4]=0.0625*(f0+f1+f2+f3)-0.375*(f4+f5+f6+f7)
                 +2.25*f8;
  
  /*
  Functionals[0]=f0;
  Functionals[2]=f1;
  Functionals[8]=f2;
  Functionals[6]=f3;
  Functionals[1]=0.125*(f0+f1)+0.75*f4;
  Functionals[5]=0.125*(f1+f2)+0.75*f5;
  Functionals[7]=0.125*(f2+f3)+0.75*f6;
  Functionals[3]=0.125*(f0+f3)+0.75*f7;
  Functionals[4]=0.015625*(f0+f1+f2+f3)+0.09375*(f4+f5+f6+f7)
                 +0.5625*f8;
  */
  /*
  cout<<"f"<<endl;
  cout<<f0<<" "<<f1<<" "<<f2<<endl;
  cout<<f3<<" "<<f4<<" "<<f5<<endl;
  cout<<f6<<" "<<f7<<" "<<f8<<endl<<endl;
  cout<<"Functionals"<<endl;
   
  cout<<Functionals[0]<<" "<<Functionals[1]<<" "<<Functionals[2]<<endl;
  cout<<Functionals[3]<<" "<<Functionals[4]<<" "<<Functionals[5]<<endl;
  cout<<Functionals[6]<<" "<<Functionals[7]<<" "<<Functionals[8]<<endl<<endl;
  */

}

/* not used */
void NF_S_Q_Q2_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues, double *Functionals)
{
  static double weights[3] = { 0.5555555555555555555555555555555556,
                               0.88888888888888888888888888888888889,
                               0.5555555555555555555555555555555556 };
  double s;

  s =(   weights[0]*PointValues[1]
	+weights[1]*PointValues[2]
	+weights[2]*PointValues[3])*0.5;

  Functionals[0] = PointValues[0];
  Functionals[1] = -0.25*(PointValues[0]+PointValues[4])
    +1.5*s;
  Functionals[2] = PointValues[4];
}

/*
    TNodalFunctional2D(NodalFunctional2D id,
                       int n_allfunctionals, int n_edgefunctionals,
                       int n_pointsall, int n_pointsedge,
                       double *xi, double *eta, double *t,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evaledge);
*/

TNodalFunctional2D *NF_S_Q_Q2_2D_Obj = new TNodalFunctional2D
        (NF_S_Q_Q2_2D, 9, 3, 25, 5, NF_S_Q_Q2_2D_Xi, NF_S_Q_Q2_2D_Eta,
         NF_S_Q_Q2_2D_T, NF_S_Q_Q2_2D_EvalAll, NF_S_Q_Q2_2D_EvalEdge);
