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
   

static double NF_D_L_P2_1D_Xi[] = { -0.77459666924148337703585307995647992,
                                     0,
                                     0.77459666924148337703585307995647992 };
static double NF_D_L_P2_1D_Eta[] = { 0, 0, 0 };
static double *NF_D_L_P2_1D_T = NULL;

void NF_D_L_P2_1D_EvalAll(double *PointValues, double *Functionals)
{

  Functionals[0] =0.5*(0.55555555555555560*PointValues[0]
                       + 0.88888888888888880*PointValues[1]
                       + 0.55555555555555560*PointValues[2]   );

  Functionals[1] =0.5*(-1.27459666924148337703585307995647992*PointValues[0]
                       + 0.5*PointValues[1]
                       + 0.27459666924148337703585307995647992 *PointValues[2]  );

  Functionals[2] = 0.5*(- 0.27459666924148337703585307995647992*PointValues[0]
                       + 0.5*PointValues[1]
                       + 1.27459666924148337703585307995647992*PointValues[2]  );
}

void NF_D_L_P2_1D_EvalEdge( double *PointValues, double *Functionals)
{

}

TNodalFunctional1D *NF_D_L_P2_1D_Obj = new TNodalFunctional1D
        (NF_D_L_P2_1D, 3, 0, 3, 0, NF_D_L_P2_1D_Xi, NF_D_L_P2_1D_Eta,
         NF_D_L_P2_1D_T, NF_D_L_P2_1D_EvalAll, NF_D_L_P2_1D_EvalEdge);
