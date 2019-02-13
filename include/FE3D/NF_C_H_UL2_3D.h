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
   
/*
    TNodalFunctional3D(NodalFunctional3D id,
                       int n_allfunctionals, int *n_facefunctionals,
                       int n_pointsall, int *n_pointsface,
                       double *xi, double *eta, double *zeta,
                       double **xiarray, double **etaarray,
                       double **zetaarray,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evalface);
*/

/* for all functionals */
static double NF_C_H_UL2_3D_Xi[]   = {-1, 0, 1, -1, 0, 1, -1, 0, 1,
                                      -1, 0, 1, 1, 1, 0, -1, -1,
									  -1, 0, 1, -1, 0, 1, -1, 0, 1,
									  -0.7745966692414833770358530, 0., 0.7745966692414833770358530,
		                              -0.7745966692414833770358530, 0., 0.7745966692414833770358530,
		                              -0.7745966692414833770358530, 0., 0.7745966692414833770358530,
		                              -0.7745966692414833770358530, 0., 0.7745966692414833770358530,
		                              -0.7745966692414833770358530, 0., 0.7745966692414833770358530,
		                              -0.7745966692414833770358530, 0., 0.7745966692414833770358530,
		                              -0.7745966692414833770358530, 0., 0.7745966692414833770358530,
		                              -0.7745966692414833770358530, 0., 0.7745966692414833770358530,
		                              -0.7745966692414833770358530, 0., 0.7745966692414833770358530};
		
static double NF_C_H_UL2_3D_Eta[]  = {-1, -1, -1, 0, 0, 0, 1, 1, 1,
                                      -1, -1, -1, 0, 1, 1, 1, 0,
									  -1, -1, -1, 0, 0, 0, 1, 1, 1,
									  -0.7745966692414833770358530, -0.7745966692414833770358530, -0.7745966692414833770358530,
									  0., 0., 0.,
									  0.7745966692414833770358530, 0.7745966692414833770358530, 0.7745966692414833770358530,
									  -0.7745966692414833770358530, -0.7745966692414833770358530, -0.7745966692414833770358530,
									  0., 0., 0.,
									  0.7745966692414833770358530, 0.7745966692414833770358530, 0.7745966692414833770358530,
									  -0.7745966692414833770358530, -0.7745966692414833770358530, -0.7745966692414833770358530,
									  0., 0., 0.,
									  0.7745966692414833770358530, 0.7745966692414833770358530, 0.7745966692414833770358530};
									  
static double NF_C_H_UL2_3D_Zeta[] = {-1, -1, -1, -1, -1, -1, -1, -1, -1,
                                       0, 0, 0, 0, 0, 0, 0, 0,
   									   1, 1, 1, 1, 1, 1, 1, 1, 1,
									   -0.7745966692414833770358530, -0.7745966692414833770358530, -0.7745966692414833770358530,
									   -0.7745966692414833770358530, -0.7745966692414833770358530, -0.7745966692414833770358530,
									   -0.7745966692414833770358530, -0.7745966692414833770358530, -0.7745966692414833770358530,
									    0., 0., 0.,
										0., 0., 0.,
										0., 0., 0.,
										0.7745966692414833770358530, 0.7745966692414833770358530, 0.7745966692414833770358530,
										0.7745966692414833770358530, 0.7745966692414833770358530, 0.7745966692414833770358530,
										0.7745966692414833770358530, 0.7745966692414833770358530, 0.7745966692414833770358530};
									 
static double NF_C_H_UL2_3D_W26[] = {0.1714677640603566529492459, 0.2743484224965706447187930, 0.1714677640603566529492459,
                                     0.2743484224965706447187929, 0.4389574759945130315500680, 0.2743484224965706447187929,
									 0.1714677640603566529492459, 0.2743484224965706447187930, 0.1714677640603566529492459,
									 0.2743484224965706447187929, 0.4389574759945130315500680, 0.2743484224965706447187929,
									 0.4389574759945130315500680, 0.7023319615912208504801077, 0.4389574759945130315500680,
									 0.2743484224965706447187929, 0.4389574759945130315500680, 0.2743484224965706447187929,
									 0.1714677640603566529492459, 0.2743484224965706447187930, 0.1714677640603566529492459,
									 0.2743484224965706447187929, 0.4389574759945130315500680, 0.2743484224965706447187929,
									 0.1714677640603566529492459, 0.2743484224965706447187930, 0.1714677640603566529492459};
									
static double NF_C_H_UL2_3D_W27[] = {-0.1328183589234367930445567, -0.2125093742774988688712904, -0.1328183589234367930445567,
                                     -0.2125093742774988688712903, -0.3400149988439981901940640, -0.2125093742774988688712903,
									 -0.1328183589234367930445567, -0.2125093742774988688712904, -0.1328183589234367930445567,
									 0., 0., 0.,
									 0., 0., 0.,
									 0., 0., 0.,
									 0.1328183589234367930445567, 0.2125093742774988688712904, 0.1328183589234367930445567,
									 0.2125093742774988688712903, 0.3400149988439981901940640, 0.2125093742774988688712903,
									 0.1328183589234367930445567, 0.2125093742774988688712904, 0.1328183589234367930445567};

static double NF_C_H_UL2_3D_W28[] = {-0.1328183589234367930445567, -0.2125093742774988688712904, -0.1328183589234367930445567,
                                     0., 0., 0.,
									 0.1328183589234367930445567, 0.2125093742774988688712904, 0.1328183589234367930445567,
 									 -0.2125093742774988688712903, -0.3400149988439981901940640, -0.2125093742774988688712903,
									 0., 0., 0.,
									 0.2125093742774988688712903, 0.3400149988439981901940640, 0.2125093742774988688712903,
									 -0.1328183589234367930445567, -0.2125093742774988688712904, -0.1328183589234367930445567,
									 0., 0., 0.,
									 0.1328183589234367930445567, 0.2125093742774988688712904, 0.1328183589234367930445567};

static double NF_C_H_UL2_3D_W29[] = {-0.1328183589234367930445567, 0., 0.1328183589234367930445567,
                                     -0.2125093742774988688712903, 0., 0.2125093742774988688712903,
									 -0.1328183589234367930445567, 0., 0.1328183589234367930445567,
									 -0.2125093742774988688712903, 0., 0.2125093742774988688712903,
									 -0.3400149988439981901940640, 0., 0.3400149988439981901940640,
									 -0.2125093742774988688712903, 0., 0.2125093742774988688712903,
									 -0.1328183589234367930445567, 0., 0.1328183589234367930445567,
									 -0.2125093742774988688712903, 0., 0.2125093742774988688712903,
									 -0.1328183589234367930445567, 0., 0.1328183589234367930445567};


/* face 0                                  0     1     2     3     4     5     6     7     8 */
static double NF_C_H_UL2_3D_F0_Xi[]   = { -1,    0,    1,   -1,    0,    1,   -1,    0,    1 };
static double NF_C_H_UL2_3D_F0_Eta[]  = { -1,   -1,   -1,    0,    0,    0,    1,    1,    1 };
static double NF_C_H_UL2_3D_F0_Zeta[] = { -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1 };

/* face 1                                  0     9    17     1    10    18     2    11    19 */
static double NF_C_H_UL2_3D_F1_Xi[]   = { -1,   -1,   -1,    0,    0,    0,    1,    1,    1 };
static double NF_C_H_UL2_3D_F1_Eta[]  = { -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1 };
static double NF_C_H_UL2_3D_F1_Zeta[] = { -1,    0,    1,   -1,    0,    1,   -1,    0,    1 };

/* face 2                                  2    11    19     5    12    22     8    13    25 */
static double NF_C_H_UL2_3D_F2_Xi[]   = {  1,    1,    1,    1,    1,    1,    1,    1,    1 };
static double NF_C_H_UL2_3D_F2_Eta[]  = { -1,   -1,   -1,    0,    0,    0,    1,    1,    1 };
static double NF_C_H_UL2_3D_F2_Zeta[] = { -1,    0,    1,   -1,    0,    1,   -1,    0,    1 };

/* face 3                                  8    13    25     7    14    24     6    15    23 */
static double NF_C_H_UL2_3D_F3_Xi[]   = {  1,    1,    1,    0,    0,    0,   -1,   -1,   -1 };
static double NF_C_H_UL2_3D_F3_Eta[]  = {  1,    1,    1,    1,    1,    1,    1,    1,    1 };
static double NF_C_H_UL2_3D_F3_Zeta[] = { -1,    0,    1,   -1,    0,    1,   -1,    0,    1 };

/* face 4                                  0     3     6     9    16    15    17    20    23 */
static double NF_C_H_UL2_3D_F4_Xi[]   = { -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1 };
static double NF_C_H_UL2_3D_F4_Eta[]  = { -1,    0,    1,   -1,    0,    1,   -1,    0,    1 };
static double NF_C_H_UL2_3D_F4_Zeta[] = { -1,   -1,   -1,    0,    0,    0,    1,    1,    1 };

/* face 5                                 17    20    23    18    21    24    19    22    25 */
static double NF_C_H_UL2_3D_F5_Xi[]   = { -1,   -1,   -1,    0,    0,    0,    1,    1,    1 };
static double NF_C_H_UL2_3D_F5_Eta[]  = { -1,    0,    1,   -1,    0,    1,   -1,    0,    1 };
static double NF_C_H_UL2_3D_F5_Zeta[] = {  1,    1,    1,    1,    1,    1,    1,    1,    1 };

static double *NF_C_H_UL2_3D_XiArray[6] = { 
                        NF_C_H_UL2_3D_F0_Xi,
                        NF_C_H_UL2_3D_F1_Xi,
                        NF_C_H_UL2_3D_F2_Xi,
                        NF_C_H_UL2_3D_F3_Xi,
                        NF_C_H_UL2_3D_F4_Xi,
                        NF_C_H_UL2_3D_F5_Xi };

static double *NF_C_H_UL2_3D_EtaArray[6] = { 
                        NF_C_H_UL2_3D_F0_Eta,
                        NF_C_H_UL2_3D_F1_Eta,
                        NF_C_H_UL2_3D_F2_Eta,
                        NF_C_H_UL2_3D_F3_Eta,
                        NF_C_H_UL2_3D_F4_Eta,
                        NF_C_H_UL2_3D_F5_Eta };

static double *NF_C_H_UL2_3D_ZetaArray[6] = { 
                        NF_C_H_UL2_3D_F0_Zeta,
                        NF_C_H_UL2_3D_F1_Zeta,
                        NF_C_H_UL2_3D_F2_Zeta,
                        NF_C_H_UL2_3D_F3_Zeta,
                        NF_C_H_UL2_3D_F4_Zeta,
                        NF_C_H_UL2_3D_F5_Zeta };

static double NF_C_H_UL2_3D_T[9] = { 0, 0.5, 1,   0, 0.5,   1, 0, 0.5, 1 };
static double NF_C_H_UL2_3D_S[9] = { 0,   0, 0, 0.5, 0.5, 0.5, 1,   1, 1 };

void NF_C_H_UL2_3D_EvalAll(TCollection *Coll, TBaseCell *Cell,
                          double *PointValues, double *Functionals)
{
  memcpy(Functionals, PointValues, 26*SizeOfDouble);
  
    Functionals[26] = NF_C_H_UL2_3D_W26[0]*PointValues[ 26]
					+NF_C_H_UL2_3D_W26[1]*PointValues[ 27]
					+NF_C_H_UL2_3D_W26[2]*PointValues[ 28]
					+NF_C_H_UL2_3D_W26[3]*PointValues[ 29]
					+NF_C_H_UL2_3D_W26[4]*PointValues[ 30]
					+NF_C_H_UL2_3D_W26[5]*PointValues[ 31]
					+NF_C_H_UL2_3D_W26[6]*PointValues[ 32]
					+NF_C_H_UL2_3D_W26[7]*PointValues[ 33]
					+NF_C_H_UL2_3D_W26[8]*PointValues[ 34]
					+NF_C_H_UL2_3D_W26[9]*PointValues[ 35]
					+NF_C_H_UL2_3D_W26[10]*PointValues[ 36]
					+NF_C_H_UL2_3D_W26[11]*PointValues[ 37]
					+NF_C_H_UL2_3D_W26[12]*PointValues[ 38]
					+NF_C_H_UL2_3D_W26[13]*PointValues[ 39]
					+NF_C_H_UL2_3D_W26[14]*PointValues[ 40]
					+NF_C_H_UL2_3D_W26[15]*PointValues[ 41]
					+NF_C_H_UL2_3D_W26[16]*PointValues[ 42]
					+NF_C_H_UL2_3D_W26[17]*PointValues[ 43]
					+NF_C_H_UL2_3D_W26[18]*PointValues[ 44]
					+NF_C_H_UL2_3D_W26[19]*PointValues[ 45]
					+NF_C_H_UL2_3D_W26[20]*PointValues[ 46]
					+NF_C_H_UL2_3D_W26[21]*PointValues[ 47]
					+NF_C_H_UL2_3D_W26[22]*PointValues[ 48]
					+NF_C_H_UL2_3D_W26[23]*PointValues[ 49]
					+NF_C_H_UL2_3D_W26[24]*PointValues[ 50]
					+NF_C_H_UL2_3D_W26[25]*PointValues[ 51]
					+NF_C_H_UL2_3D_W26[26]*PointValues[ 52];
    Functionals[27] = NF_C_H_UL2_3D_W27[0]*PointValues[ 26]
					+NF_C_H_UL2_3D_W27[1]*PointValues[ 27]
					+NF_C_H_UL2_3D_W27[2]*PointValues[ 28]
					+NF_C_H_UL2_3D_W27[3]*PointValues[ 29]
					+NF_C_H_UL2_3D_W27[4]*PointValues[ 30]
					+NF_C_H_UL2_3D_W27[5]*PointValues[ 31]
					+NF_C_H_UL2_3D_W27[6]*PointValues[ 32]
					+NF_C_H_UL2_3D_W27[7]*PointValues[ 33]
					+NF_C_H_UL2_3D_W27[8]*PointValues[ 34]
					+NF_C_H_UL2_3D_W27[9]*PointValues[ 35]
					+NF_C_H_UL2_3D_W27[10]*PointValues[ 36]
					+NF_C_H_UL2_3D_W27[11]*PointValues[ 37]
					+NF_C_H_UL2_3D_W27[12]*PointValues[ 38]
					+NF_C_H_UL2_3D_W27[13]*PointValues[ 39]
					+NF_C_H_UL2_3D_W27[14]*PointValues[ 40]
					+NF_C_H_UL2_3D_W27[15]*PointValues[ 41]
					+NF_C_H_UL2_3D_W27[16]*PointValues[ 42]
					+NF_C_H_UL2_3D_W27[17]*PointValues[ 43]
					+NF_C_H_UL2_3D_W27[18]*PointValues[ 44]
					+NF_C_H_UL2_3D_W27[19]*PointValues[ 45]
					+NF_C_H_UL2_3D_W27[20]*PointValues[ 46]
					+NF_C_H_UL2_3D_W27[21]*PointValues[ 47]
					+NF_C_H_UL2_3D_W27[22]*PointValues[ 48]
					+NF_C_H_UL2_3D_W27[23]*PointValues[ 49]
					+NF_C_H_UL2_3D_W27[24]*PointValues[ 50]
					+NF_C_H_UL2_3D_W27[25]*PointValues[ 51]
					+NF_C_H_UL2_3D_W27[26]*PointValues[ 52];
    Functionals[28] = NF_C_H_UL2_3D_W28[0]*PointValues[ 26]
					+NF_C_H_UL2_3D_W28[1]*PointValues[ 27]
					+NF_C_H_UL2_3D_W28[2]*PointValues[ 28]
					+NF_C_H_UL2_3D_W28[3]*PointValues[ 29]
					+NF_C_H_UL2_3D_W28[4]*PointValues[ 30]
					+NF_C_H_UL2_3D_W28[5]*PointValues[ 31]
					+NF_C_H_UL2_3D_W28[6]*PointValues[ 32]
					+NF_C_H_UL2_3D_W28[7]*PointValues[ 33]
					+NF_C_H_UL2_3D_W28[8]*PointValues[ 34]
					+NF_C_H_UL2_3D_W28[9]*PointValues[ 35]
					+NF_C_H_UL2_3D_W28[10]*PointValues[ 36]
					+NF_C_H_UL2_3D_W28[11]*PointValues[ 37]
					+NF_C_H_UL2_3D_W28[12]*PointValues[ 38]
					+NF_C_H_UL2_3D_W28[13]*PointValues[ 39]
					+NF_C_H_UL2_3D_W28[14]*PointValues[ 40]
					+NF_C_H_UL2_3D_W28[15]*PointValues[ 41]
					+NF_C_H_UL2_3D_W28[16]*PointValues[ 42]
					+NF_C_H_UL2_3D_W28[17]*PointValues[ 43]
					+NF_C_H_UL2_3D_W28[18]*PointValues[ 44]
					+NF_C_H_UL2_3D_W28[19]*PointValues[ 45]
					+NF_C_H_UL2_3D_W28[20]*PointValues[ 46]
					+NF_C_H_UL2_3D_W28[21]*PointValues[ 47]
					+NF_C_H_UL2_3D_W28[22]*PointValues[ 48]
					+NF_C_H_UL2_3D_W28[23]*PointValues[ 49]
					+NF_C_H_UL2_3D_W28[24]*PointValues[ 50]
					+NF_C_H_UL2_3D_W28[25]*PointValues[ 51]
					+NF_C_H_UL2_3D_W28[26]*PointValues[ 52];
    Functionals[29] = NF_C_H_UL2_3D_W29[0]*PointValues[ 26]
					+NF_C_H_UL2_3D_W29[1]*PointValues[ 27]
					+NF_C_H_UL2_3D_W29[2]*PointValues[ 28]
					+NF_C_H_UL2_3D_W29[3]*PointValues[ 29]
					+NF_C_H_UL2_3D_W29[4]*PointValues[ 30]
					+NF_C_H_UL2_3D_W29[5]*PointValues[ 31]
					+NF_C_H_UL2_3D_W29[6]*PointValues[ 32]
					+NF_C_H_UL2_3D_W29[7]*PointValues[ 33]
					+NF_C_H_UL2_3D_W29[8]*PointValues[ 34]
					+NF_C_H_UL2_3D_W29[9]*PointValues[ 35]
					+NF_C_H_UL2_3D_W29[10]*PointValues[ 36]
					+NF_C_H_UL2_3D_W29[11]*PointValues[ 37]
					+NF_C_H_UL2_3D_W29[12]*PointValues[ 38]
					+NF_C_H_UL2_3D_W29[13]*PointValues[ 39]
					+NF_C_H_UL2_3D_W29[14]*PointValues[ 40]
					+NF_C_H_UL2_3D_W29[15]*PointValues[ 41]
					+NF_C_H_UL2_3D_W29[16]*PointValues[ 42]
					+NF_C_H_UL2_3D_W29[17]*PointValues[ 43]
					+NF_C_H_UL2_3D_W29[18]*PointValues[ 44]
					+NF_C_H_UL2_3D_W29[19]*PointValues[ 45]
					+NF_C_H_UL2_3D_W29[20]*PointValues[ 46]
					+NF_C_H_UL2_3D_W29[21]*PointValues[ 47]
					+NF_C_H_UL2_3D_W29[22]*PointValues[ 48]
					+NF_C_H_UL2_3D_W29[23]*PointValues[ 49]
					+NF_C_H_UL2_3D_W29[24]*PointValues[ 50]
					+NF_C_H_UL2_3D_W29[25]*PointValues[ 51]
					+NF_C_H_UL2_3D_W29[26]*PointValues[ 52];
}

void NF_C_H_UL2_3D_EvalFace(TCollection *Coll, TBaseCell *Cell, int Joint, 
                           double *PointValues, double *Functionals)
{
  memcpy(Functionals, PointValues, 9*SizeOfDouble);
}

static int NF_C_H_UL2_3D_N_AllFunctionals = 30;
static int NF_C_H_UL2_3D_N_PointsAll = 53;
static int NF_C_H_UL2_3D_N_FaceFunctionals[] = { 9, 9, 9, 9, 9, 9 };
static int NF_C_H_UL2_3D_N_PointsFace[] = { 9, 9, 9, 9, 9, 9 };

TNodalFunctional3D *NF_C_H_UL2_3D_Obj = new TNodalFunctional3D
        (NF_C_H_UL2_3D, NF_C_H_UL2_3D_N_AllFunctionals,
         NF_C_H_UL2_3D_N_FaceFunctionals, NF_C_H_UL2_3D_N_PointsAll,
         NF_C_H_UL2_3D_N_PointsFace,
         NF_C_H_UL2_3D_Xi, NF_C_H_UL2_3D_Eta, NF_C_H_UL2_3D_Zeta,
         NF_C_H_UL2_3D_XiArray, NF_C_H_UL2_3D_EtaArray,
         NF_C_H_UL2_3D_ZetaArray,
         NF_C_H_UL2_3D_T, NF_C_H_UL2_3D_S,
         NF_C_H_UL2_3D_EvalAll, NF_C_H_UL2_3D_EvalFace);
