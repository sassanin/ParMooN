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
   
static double NF_D_Q_Q3_2D_Xi[16]={
   -0.86113631159405257522, -0.3399810435848562648,
    0.3399810435848562648, 0.86113631159405257522,
   -0.86113631159405257522, -0.3399810435848562648,
    0.3399810435848562648, 0.86113631159405257522,
   -0.86113631159405257522, -0.3399810435848562648,
    0.3399810435848562648, 0.86113631159405257522,
   -0.86113631159405257522, -0.3399810435848562648,
    0.3399810435848562648, 0.86113631159405257522 };

static double NF_D_Q_Q3_2D_Eta[16]={
   -0.86113631159405257522, -0.86113631159405257522,
   -0.86113631159405257522, -0.86113631159405257522,
   -0.3399810435848562648, -0.3399810435848562648,
   -0.3399810435848562648, -0.3399810435848562648,
    0.3399810435848562648, 0.3399810435848562648,
    0.3399810435848562648, 0.3399810435848562648,
    0.86113631159405257522, 0.86113631159405257522,
    0.86113631159405257522, 0.86113631159405257522 };

static double NF_D_Q_Q3_2D_Weight0[16]={
    0.03025074832140050138, 0.056712962962962962963,
    0.056712962962962962963, 0.03025074832140050138,
    0.056712962962962962963, 0.10632332575267357269,
    0.10632332575267357269, 0.056712962962962962963,
    0.056712962962962962963, 0.10632332575267357269,
    0.10632332575267357269, 0.056712962962962962963,
    0.03025074832140050138, 0.056712962962962962963,
    0.056712962962962962963, 0.03025074832140050138 };

static double NF_D_Q_Q3_2D_Weight1[16] = { 
   -0.078150053497352415165, -0.057843996998812350609,
    0.057843996998812350609, 0.078150053497352415165,
   -0.14651277523648811176, -0.10844374574041975333,
    0.10844374574041975333, 0.14651277523648811176,
   -0.14651277523648811176, -0.10844374574041975333,
    0.10844374574041975333, 0.14651277523648811176,
   -0.078150053497352415165, -0.057843996998812350609,
    0.057843996998812350609, 0.078150053497352415165 };

static double NF_D_Q_Q3_2D_Weight2[16] = { 
    0.092617751245468615539, -0.092617751245468615539,
   -0.092617751245468615539, 0.092617751245468615539,
    0.17363626976398713406, -0.17363626976398713406,
   -0.17363626976398713406, 0.17363626976398713406,
    0.17363626976398713406, -0.17363626976398713406,
   -0.17363626976398713406, 0.17363626976398713406,
    0.092617751245468615539, -0.092617751245468615539,
   -0.092617751245468615539, 0.092617751245468615539 };

static double NF_D_Q_Q3_2D_Weight3[16] = { 
   -0.064531770405098972688, 0.16345220357384786324,
   -0.16345220357384786324, 0.064531770405098972688,
   -0.12098173129587467994, 0.30643403161502865704,
   -0.30643403161502865704, 0.12098173129587467994,
   -0.12098173129587467994, 0.30643403161502865704,
   -0.30643403161502865704, 0.12098173129587467994,
   -0.064531770405098972688, 0.16345220357384786324,
   -0.16345220357384786324, 0.064531770405098972688 };

static double NF_D_Q_Q3_2D_Weight4[16] = { 
   -0.078150053497352415165, -0.14651277523648811176,
   -0.14651277523648811176, -0.078150053497352415165,
   -0.057843996998812350609, -0.10844374574041975333,
   -0.10844374574041975333, -0.057843996998812350609,
    0.057843996998812350609, 0.10844374574041975333,
    0.10844374574041975333, 0.057843996998812350609,
    0.078150053497352415165, 0.14651277523648811176,
    0.14651277523648811176, 0.078150053497352415165 };

static double NF_D_Q_Q3_2D_Weight5[16] = { 
    0.20189354645876384279, 0.14943469867024414309,
   -0.14943469867024414309, -0.20189354645876384279,
    0.14943469867024414309, 0.11060645354123615721,
   -0.11060645354123615721, -0.14943469867024414309,
   -0.14943469867024414309, -0.11060645354123615721,
    0.11060645354123615721, 0.14943469867024414309,
   -0.20189354645876384279, -0.14943469867024414309,
    0.14943469867024414309, 0.20189354645876384279 };

static double NF_D_Q_Q3_2D_Weight6[16] = { 
   -0.23926952608697493807, 0.23926952608697493807,
    0.23926952608697493807, -0.23926952608697493807,
   -0.17709912059562590953, 0.17709912059562590953,
    0.17709912059562590953, -0.17709912059562590953,
    0.17709912059562590953, -0.17709912059562590953,
   -0.17709912059562590953, 0.17709912059562590953,
    0.23926952608697493807, -0.23926952608697493807,
   -0.23926952608697493807, 0.23926952608697493807 };

static double NF_D_Q_Q3_2D_Weight7[16] = { 
    0.16671195224184350797, -0.42226388312251070245,
    0.42226388312251070245, -0.16671195224184350797,
    0.12339448578202241628, -0.3125452855751768413,
    0.3125452855751768413, -0.12339448578202241628,
   -0.12339448578202241628, 0.3125452855751768413,
   -0.3125452855751768413, 0.12339448578202241628,
   -0.16671195224184350797, 0.42226388312251070245,
   -0.42226388312251070245, 0.16671195224184350797 };

static double NF_D_Q_Q3_2D_Weight8[16] = { 
    0.092617751245468615539, 0.17363626976398713406,
    0.17363626976398713406, 0.092617751245468615539,
   -0.092617751245468615539, -0.17363626976398713406,
   -0.17363626976398713406, -0.092617751245468615539,
   -0.092617751245468615539, -0.17363626976398713406,
   -0.17363626976398713406, -0.092617751245468615539,
    0.092617751245468615539, 0.17363626976398713406,
    0.17363626976398713406, 0.092617751245468615539 };

static double NF_D_Q_Q3_2D_Weight9[16] = { 
   -0.23926952608697493807, -0.17709912059562590953,
    0.17709912059562590953, 0.23926952608697493807,
    0.23926952608697493807, 0.17709912059562590953,
   -0.17709912059562590953, -0.23926952608697493807,
    0.23926952608697493807, 0.17709912059562590953,
   -0.17709912059562590953, -0.23926952608697493807,
   -0.23926952608697493807, -0.17709912059562590953,
    0.17709912059562590953, 0.23926952608697493807 };

static double NF_D_Q_Q3_2D_Weight10[16] = { 
    0.28356481481481481481, -0.28356481481481481481,
   -0.28356481481481481481, 0.28356481481481481481,
   -0.28356481481481481481, 0.28356481481481481481,
    0.28356481481481481481, -0.28356481481481481481,
   -0.28356481481481481481, 0.28356481481481481481,
    0.28356481481481481481, -0.28356481481481481481,
    0.28356481481481481481, -0.28356481481481481481,
   -0.28356481481481481481, 0.28356481481481481481 };

static double NF_D_Q_Q3_2D_Weight11[16] = { 
   -0.19757486311771497539, 0.50043639814413277831,
   -0.50043639814413277831, 0.19757486311771497539,
    0.19757486311771497539, -0.50043639814413277831,
    0.50043639814413277831, -0.19757486311771497539,
    0.19757486311771497539, -0.50043639814413277831,
    0.50043639814413277831, -0.19757486311771497539,
   -0.19757486311771497539, 0.50043639814413277831,
   -0.50043639814413277831, 0.19757486311771497539 };

static double NF_D_Q_Q3_2D_Weight12[16] = { 
   -0.064531770405098972688, -0.12098173129587467994,
   -0.12098173129587467994, -0.064531770405098972688,
    0.16345220357384786324, 0.30643403161502865704,
    0.30643403161502865704, 0.16345220357384786324,
   -0.16345220357384786324, -0.30643403161502865704,
   -0.30643403161502865704, -0.16345220357384786324,
    0.064531770405098972688, 0.12098173129587467994,
    0.12098173129587467994, 0.064531770405098972688 };

static double NF_D_Q_Q3_2D_Weight13[16] = { 
    0.16671195224184350797, 0.12339448578202241628,
   -0.12339448578202241628, -0.16671195224184350797,
   -0.42226388312251070245, -0.3125452855751768413,
    0.3125452855751768413, 0.42226388312251070245,
    0.42226388312251070245, 0.3125452855751768413,
   -0.3125452855751768413, -0.42226388312251070245,
   -0.16671195224184350797, -0.12339448578202241628,
    0.12339448578202241628, 0.16671195224184350797 };

static double NF_D_Q_Q3_2D_Weight14[16] = { 
   -0.19757486311771497539, 0.19757486311771497539,
    0.19757486311771497539, -0.19757486311771497539,
    0.50043639814413277831, -0.50043639814413277831,
   -0.50043639814413277831, 0.50043639814413277831,
   -0.50043639814413277831, 0.50043639814413277831,
    0.50043639814413277831, -0.50043639814413277831,
    0.19757486311771497539, -0.19757486311771497539,
   -0.19757486311771497539, 0.19757486311771497539 };

static double NF_D_Q_Q3_2D_Weight15[16] = { 
    0.13766103725342861723, -0.34868096356390300054,
    0.34868096356390300054, -0.13766103725342861723,
   -0.34868096356390300054, 0.8831722960799047161,
   -0.8831722960799047161, 0.34868096356390300054,
    0.34868096356390300054, -0.8831722960799047161,
    0.8831722960799047161, -0.34868096356390300054,
   -0.13766103725342861723, 0.34868096356390300054,
   -0.34868096356390300054, 0.13766103725342861723 };

static double *NF_D_Q_Q3_2D_T = NULL;

void NF_D_Q_Q3_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  Functionals[0] =  NF_D_Q_Q3_2D_Weight0[0]*PointValues[0]
                   +NF_D_Q_Q3_2D_Weight0[1]*PointValues[1]
                   +NF_D_Q_Q3_2D_Weight0[2]*PointValues[2]
                   +NF_D_Q_Q3_2D_Weight0[3]*PointValues[3]
                   +NF_D_Q_Q3_2D_Weight0[4]*PointValues[4]
                   +NF_D_Q_Q3_2D_Weight0[5]*PointValues[5]
                   +NF_D_Q_Q3_2D_Weight0[6]*PointValues[6]
                   +NF_D_Q_Q3_2D_Weight0[7]*PointValues[7]
                   +NF_D_Q_Q3_2D_Weight0[8]*PointValues[8]
                   +NF_D_Q_Q3_2D_Weight0[9]*PointValues[9]
                   +NF_D_Q_Q3_2D_Weight0[10]*PointValues[10]
                   +NF_D_Q_Q3_2D_Weight0[11]*PointValues[11]
                   +NF_D_Q_Q3_2D_Weight0[12]*PointValues[12]
                   +NF_D_Q_Q3_2D_Weight0[13]*PointValues[13]
                   +NF_D_Q_Q3_2D_Weight0[14]*PointValues[14]
                   +NF_D_Q_Q3_2D_Weight0[15]*PointValues[15];

  Functionals[1] =  NF_D_Q_Q3_2D_Weight1[0]*PointValues[0]
                   +NF_D_Q_Q3_2D_Weight1[1]*PointValues[1]
                   +NF_D_Q_Q3_2D_Weight1[2]*PointValues[2]
                   +NF_D_Q_Q3_2D_Weight1[3]*PointValues[3]
                   +NF_D_Q_Q3_2D_Weight1[4]*PointValues[4]
                   +NF_D_Q_Q3_2D_Weight1[5]*PointValues[5]
                   +NF_D_Q_Q3_2D_Weight1[6]*PointValues[6]
                   +NF_D_Q_Q3_2D_Weight1[7]*PointValues[7]
                   +NF_D_Q_Q3_2D_Weight1[8]*PointValues[8]
                   +NF_D_Q_Q3_2D_Weight1[9]*PointValues[9]
                   +NF_D_Q_Q3_2D_Weight1[10]*PointValues[10]
                   +NF_D_Q_Q3_2D_Weight1[11]*PointValues[11]
                   +NF_D_Q_Q3_2D_Weight1[12]*PointValues[12]
                   +NF_D_Q_Q3_2D_Weight1[13]*PointValues[13]
                   +NF_D_Q_Q3_2D_Weight1[14]*PointValues[14]
                   +NF_D_Q_Q3_2D_Weight1[15]*PointValues[15];

  Functionals[2] =  NF_D_Q_Q3_2D_Weight2[0]*PointValues[0]
                   +NF_D_Q_Q3_2D_Weight2[1]*PointValues[1]
                   +NF_D_Q_Q3_2D_Weight2[2]*PointValues[2]
                   +NF_D_Q_Q3_2D_Weight2[3]*PointValues[3]
                   +NF_D_Q_Q3_2D_Weight2[4]*PointValues[4]
                   +NF_D_Q_Q3_2D_Weight2[5]*PointValues[5]
                   +NF_D_Q_Q3_2D_Weight2[6]*PointValues[6]
                   +NF_D_Q_Q3_2D_Weight2[7]*PointValues[7]
                   +NF_D_Q_Q3_2D_Weight2[8]*PointValues[8]
                   +NF_D_Q_Q3_2D_Weight2[9]*PointValues[9]
                   +NF_D_Q_Q3_2D_Weight2[10]*PointValues[10]
                   +NF_D_Q_Q3_2D_Weight2[11]*PointValues[11]
                   +NF_D_Q_Q3_2D_Weight2[12]*PointValues[12]
                   +NF_D_Q_Q3_2D_Weight2[13]*PointValues[13]
                   +NF_D_Q_Q3_2D_Weight2[14]*PointValues[14]
                   +NF_D_Q_Q3_2D_Weight2[15]*PointValues[15];

  Functionals[3] =  NF_D_Q_Q3_2D_Weight3[0]*PointValues[0]
                   +NF_D_Q_Q3_2D_Weight3[1]*PointValues[1]
                   +NF_D_Q_Q3_2D_Weight3[2]*PointValues[2]
                   +NF_D_Q_Q3_2D_Weight3[3]*PointValues[3]
                   +NF_D_Q_Q3_2D_Weight3[4]*PointValues[4]
                   +NF_D_Q_Q3_2D_Weight3[5]*PointValues[5]
                   +NF_D_Q_Q3_2D_Weight3[6]*PointValues[6]
                   +NF_D_Q_Q3_2D_Weight3[7]*PointValues[7]
                   +NF_D_Q_Q3_2D_Weight3[8]*PointValues[8]
                   +NF_D_Q_Q3_2D_Weight3[9]*PointValues[9]
                   +NF_D_Q_Q3_2D_Weight3[10]*PointValues[10]
                   +NF_D_Q_Q3_2D_Weight3[11]*PointValues[11]
                   +NF_D_Q_Q3_2D_Weight3[12]*PointValues[12]
                   +NF_D_Q_Q3_2D_Weight3[13]*PointValues[13]
                   +NF_D_Q_Q3_2D_Weight3[14]*PointValues[14]
                   +NF_D_Q_Q3_2D_Weight3[15]*PointValues[15];

  Functionals[4] =  NF_D_Q_Q3_2D_Weight4[0]*PointValues[0]
                   +NF_D_Q_Q3_2D_Weight4[1]*PointValues[1]
                   +NF_D_Q_Q3_2D_Weight4[2]*PointValues[2]
                   +NF_D_Q_Q3_2D_Weight4[3]*PointValues[3]
                   +NF_D_Q_Q3_2D_Weight4[4]*PointValues[4]
                   +NF_D_Q_Q3_2D_Weight4[5]*PointValues[5]
                   +NF_D_Q_Q3_2D_Weight4[6]*PointValues[6]
                   +NF_D_Q_Q3_2D_Weight4[7]*PointValues[7]
                   +NF_D_Q_Q3_2D_Weight4[8]*PointValues[8]
                   +NF_D_Q_Q3_2D_Weight4[9]*PointValues[9]
                   +NF_D_Q_Q3_2D_Weight4[10]*PointValues[10]
                   +NF_D_Q_Q3_2D_Weight4[11]*PointValues[11]
                   +NF_D_Q_Q3_2D_Weight4[12]*PointValues[12]
                   +NF_D_Q_Q3_2D_Weight4[13]*PointValues[13]
                   +NF_D_Q_Q3_2D_Weight4[14]*PointValues[14]
                   +NF_D_Q_Q3_2D_Weight4[15]*PointValues[15];

  Functionals[5] =  NF_D_Q_Q3_2D_Weight5[0]*PointValues[0]
                   +NF_D_Q_Q3_2D_Weight5[1]*PointValues[1]
                   +NF_D_Q_Q3_2D_Weight5[2]*PointValues[2]
                   +NF_D_Q_Q3_2D_Weight5[3]*PointValues[3]
                   +NF_D_Q_Q3_2D_Weight5[4]*PointValues[4]
                   +NF_D_Q_Q3_2D_Weight5[5]*PointValues[5]
                   +NF_D_Q_Q3_2D_Weight5[6]*PointValues[6]
                   +NF_D_Q_Q3_2D_Weight5[7]*PointValues[7]
                   +NF_D_Q_Q3_2D_Weight5[8]*PointValues[8]
                   +NF_D_Q_Q3_2D_Weight5[9]*PointValues[9]
                   +NF_D_Q_Q3_2D_Weight5[10]*PointValues[10]
                   +NF_D_Q_Q3_2D_Weight5[11]*PointValues[11]
                   +NF_D_Q_Q3_2D_Weight5[12]*PointValues[12]
                   +NF_D_Q_Q3_2D_Weight5[13]*PointValues[13]
                   +NF_D_Q_Q3_2D_Weight5[14]*PointValues[14]
                   +NF_D_Q_Q3_2D_Weight5[15]*PointValues[15];

  Functionals[6] =  NF_D_Q_Q3_2D_Weight6[0]*PointValues[0]
                   +NF_D_Q_Q3_2D_Weight6[1]*PointValues[1]
                   +NF_D_Q_Q3_2D_Weight6[2]*PointValues[2]
                   +NF_D_Q_Q3_2D_Weight6[3]*PointValues[3]
                   +NF_D_Q_Q3_2D_Weight6[4]*PointValues[4]
                   +NF_D_Q_Q3_2D_Weight6[5]*PointValues[5]
                   +NF_D_Q_Q3_2D_Weight6[6]*PointValues[6]
                   +NF_D_Q_Q3_2D_Weight6[7]*PointValues[7]
                   +NF_D_Q_Q3_2D_Weight6[8]*PointValues[8]
                   +NF_D_Q_Q3_2D_Weight6[9]*PointValues[9]
                   +NF_D_Q_Q3_2D_Weight6[10]*PointValues[10]
                   +NF_D_Q_Q3_2D_Weight6[11]*PointValues[11]
                   +NF_D_Q_Q3_2D_Weight6[12]*PointValues[12]
                   +NF_D_Q_Q3_2D_Weight6[13]*PointValues[13]
                   +NF_D_Q_Q3_2D_Weight6[14]*PointValues[14]
                   +NF_D_Q_Q3_2D_Weight6[15]*PointValues[15];

  Functionals[7] =  NF_D_Q_Q3_2D_Weight7[0]*PointValues[0]
                   +NF_D_Q_Q3_2D_Weight7[1]*PointValues[1]
                   +NF_D_Q_Q3_2D_Weight7[2]*PointValues[2]
                   +NF_D_Q_Q3_2D_Weight7[3]*PointValues[3]
                   +NF_D_Q_Q3_2D_Weight7[4]*PointValues[4]
                   +NF_D_Q_Q3_2D_Weight7[5]*PointValues[5]
                   +NF_D_Q_Q3_2D_Weight7[6]*PointValues[6]
                   +NF_D_Q_Q3_2D_Weight7[7]*PointValues[7]
                   +NF_D_Q_Q3_2D_Weight7[8]*PointValues[8]
                   +NF_D_Q_Q3_2D_Weight7[9]*PointValues[9]
                   +NF_D_Q_Q3_2D_Weight7[10]*PointValues[10]
                   +NF_D_Q_Q3_2D_Weight7[11]*PointValues[11]
                   +NF_D_Q_Q3_2D_Weight7[12]*PointValues[12]
                   +NF_D_Q_Q3_2D_Weight7[13]*PointValues[13]
                   +NF_D_Q_Q3_2D_Weight7[14]*PointValues[14]
                   +NF_D_Q_Q3_2D_Weight7[15]*PointValues[15];

  Functionals[8] =  NF_D_Q_Q3_2D_Weight8[0]*PointValues[0]
                   +NF_D_Q_Q3_2D_Weight8[1]*PointValues[1]
                   +NF_D_Q_Q3_2D_Weight8[2]*PointValues[2]
                   +NF_D_Q_Q3_2D_Weight8[3]*PointValues[3]
                   +NF_D_Q_Q3_2D_Weight8[4]*PointValues[4]
                   +NF_D_Q_Q3_2D_Weight8[5]*PointValues[5]
                   +NF_D_Q_Q3_2D_Weight8[6]*PointValues[6]
                   +NF_D_Q_Q3_2D_Weight8[7]*PointValues[7]
                   +NF_D_Q_Q3_2D_Weight8[8]*PointValues[8]
                   +NF_D_Q_Q3_2D_Weight8[9]*PointValues[9]
                   +NF_D_Q_Q3_2D_Weight8[10]*PointValues[10]
                   +NF_D_Q_Q3_2D_Weight8[11]*PointValues[11]
                   +NF_D_Q_Q3_2D_Weight8[12]*PointValues[12]
                   +NF_D_Q_Q3_2D_Weight8[13]*PointValues[13]
                   +NF_D_Q_Q3_2D_Weight8[14]*PointValues[14]
                   +NF_D_Q_Q3_2D_Weight8[15]*PointValues[15];

  Functionals[9] =  NF_D_Q_Q3_2D_Weight9[0]*PointValues[0]
                   +NF_D_Q_Q3_2D_Weight9[1]*PointValues[1]
                   +NF_D_Q_Q3_2D_Weight9[2]*PointValues[2]
                   +NF_D_Q_Q3_2D_Weight9[3]*PointValues[3]
                   +NF_D_Q_Q3_2D_Weight9[4]*PointValues[4]
                   +NF_D_Q_Q3_2D_Weight9[5]*PointValues[5]
                   +NF_D_Q_Q3_2D_Weight9[6]*PointValues[6]
                   +NF_D_Q_Q3_2D_Weight9[7]*PointValues[7]
                   +NF_D_Q_Q3_2D_Weight9[8]*PointValues[8]
                   +NF_D_Q_Q3_2D_Weight9[9]*PointValues[9]
                   +NF_D_Q_Q3_2D_Weight9[10]*PointValues[10]
                   +NF_D_Q_Q3_2D_Weight9[11]*PointValues[11]
                   +NF_D_Q_Q3_2D_Weight9[12]*PointValues[12]
                   +NF_D_Q_Q3_2D_Weight9[13]*PointValues[13]
                   +NF_D_Q_Q3_2D_Weight9[14]*PointValues[14]
                   +NF_D_Q_Q3_2D_Weight9[15]*PointValues[15];

  Functionals[10] =  NF_D_Q_Q3_2D_Weight10[0]*PointValues[0]
                   +NF_D_Q_Q3_2D_Weight10[1]*PointValues[1]
                   +NF_D_Q_Q3_2D_Weight10[2]*PointValues[2]
                   +NF_D_Q_Q3_2D_Weight10[3]*PointValues[3]
                   +NF_D_Q_Q3_2D_Weight10[4]*PointValues[4]
                   +NF_D_Q_Q3_2D_Weight10[5]*PointValues[5]
                   +NF_D_Q_Q3_2D_Weight10[6]*PointValues[6]
                   +NF_D_Q_Q3_2D_Weight10[7]*PointValues[7]
                   +NF_D_Q_Q3_2D_Weight10[8]*PointValues[8]
                   +NF_D_Q_Q3_2D_Weight10[9]*PointValues[9]
                   +NF_D_Q_Q3_2D_Weight10[10]*PointValues[10]
                   +NF_D_Q_Q3_2D_Weight10[11]*PointValues[11]
                   +NF_D_Q_Q3_2D_Weight10[12]*PointValues[12]
                   +NF_D_Q_Q3_2D_Weight10[13]*PointValues[13]
                   +NF_D_Q_Q3_2D_Weight10[14]*PointValues[14]
                   +NF_D_Q_Q3_2D_Weight10[15]*PointValues[15];

  Functionals[11] =  NF_D_Q_Q3_2D_Weight11[0]*PointValues[0]
                   +NF_D_Q_Q3_2D_Weight11[1]*PointValues[1]
                   +NF_D_Q_Q3_2D_Weight11[2]*PointValues[2]
                   +NF_D_Q_Q3_2D_Weight11[3]*PointValues[3]
                   +NF_D_Q_Q3_2D_Weight11[4]*PointValues[4]
                   +NF_D_Q_Q3_2D_Weight11[5]*PointValues[5]
                   +NF_D_Q_Q3_2D_Weight11[6]*PointValues[6]
                   +NF_D_Q_Q3_2D_Weight11[7]*PointValues[7]
                   +NF_D_Q_Q3_2D_Weight11[8]*PointValues[8]
                   +NF_D_Q_Q3_2D_Weight11[9]*PointValues[9]
                   +NF_D_Q_Q3_2D_Weight11[10]*PointValues[10]
                   +NF_D_Q_Q3_2D_Weight11[11]*PointValues[11]
                   +NF_D_Q_Q3_2D_Weight11[12]*PointValues[12]
                   +NF_D_Q_Q3_2D_Weight11[13]*PointValues[13]
                   +NF_D_Q_Q3_2D_Weight11[14]*PointValues[14]
                   +NF_D_Q_Q3_2D_Weight11[15]*PointValues[15];

  Functionals[12] =  NF_D_Q_Q3_2D_Weight12[0]*PointValues[0]
                   +NF_D_Q_Q3_2D_Weight12[1]*PointValues[1]
                   +NF_D_Q_Q3_2D_Weight12[2]*PointValues[2]
                   +NF_D_Q_Q3_2D_Weight12[3]*PointValues[3]
                   +NF_D_Q_Q3_2D_Weight12[4]*PointValues[4]
                   +NF_D_Q_Q3_2D_Weight12[5]*PointValues[5]
                   +NF_D_Q_Q3_2D_Weight12[6]*PointValues[6]
                   +NF_D_Q_Q3_2D_Weight12[7]*PointValues[7]
                   +NF_D_Q_Q3_2D_Weight12[8]*PointValues[8]
                   +NF_D_Q_Q3_2D_Weight12[9]*PointValues[9]
                   +NF_D_Q_Q3_2D_Weight12[10]*PointValues[10]
                   +NF_D_Q_Q3_2D_Weight12[11]*PointValues[11]
                   +NF_D_Q_Q3_2D_Weight12[12]*PointValues[12]
                   +NF_D_Q_Q3_2D_Weight12[13]*PointValues[13]
                   +NF_D_Q_Q3_2D_Weight12[14]*PointValues[14]
                   +NF_D_Q_Q3_2D_Weight12[15]*PointValues[15];

  Functionals[13] =  NF_D_Q_Q3_2D_Weight13[0]*PointValues[0]
                   +NF_D_Q_Q3_2D_Weight13[1]*PointValues[1]
                   +NF_D_Q_Q3_2D_Weight13[2]*PointValues[2]
                   +NF_D_Q_Q3_2D_Weight13[3]*PointValues[3]
                   +NF_D_Q_Q3_2D_Weight13[4]*PointValues[4]
                   +NF_D_Q_Q3_2D_Weight13[5]*PointValues[5]
                   +NF_D_Q_Q3_2D_Weight13[6]*PointValues[6]
                   +NF_D_Q_Q3_2D_Weight13[7]*PointValues[7]
                   +NF_D_Q_Q3_2D_Weight13[8]*PointValues[8]
                   +NF_D_Q_Q3_2D_Weight13[9]*PointValues[9]
                   +NF_D_Q_Q3_2D_Weight13[10]*PointValues[10]
                   +NF_D_Q_Q3_2D_Weight13[11]*PointValues[11]
                   +NF_D_Q_Q3_2D_Weight13[12]*PointValues[12]
                   +NF_D_Q_Q3_2D_Weight13[13]*PointValues[13]
                   +NF_D_Q_Q3_2D_Weight13[14]*PointValues[14]
                   +NF_D_Q_Q3_2D_Weight13[15]*PointValues[15];

  Functionals[14] =  NF_D_Q_Q3_2D_Weight14[0]*PointValues[0]
                   +NF_D_Q_Q3_2D_Weight14[1]*PointValues[1]
                   +NF_D_Q_Q3_2D_Weight14[2]*PointValues[2]
                   +NF_D_Q_Q3_2D_Weight14[3]*PointValues[3]
                   +NF_D_Q_Q3_2D_Weight14[4]*PointValues[4]
                   +NF_D_Q_Q3_2D_Weight14[5]*PointValues[5]
                   +NF_D_Q_Q3_2D_Weight14[6]*PointValues[6]
                   +NF_D_Q_Q3_2D_Weight14[7]*PointValues[7]
                   +NF_D_Q_Q3_2D_Weight14[8]*PointValues[8]
                   +NF_D_Q_Q3_2D_Weight14[9]*PointValues[9]
                   +NF_D_Q_Q3_2D_Weight14[10]*PointValues[10]
                   +NF_D_Q_Q3_2D_Weight14[11]*PointValues[11]
                   +NF_D_Q_Q3_2D_Weight14[12]*PointValues[12]
                   +NF_D_Q_Q3_2D_Weight14[13]*PointValues[13]
                   +NF_D_Q_Q3_2D_Weight14[14]*PointValues[14]
                   +NF_D_Q_Q3_2D_Weight14[15]*PointValues[15];

  Functionals[15] =  NF_D_Q_Q3_2D_Weight15[0]*PointValues[0]
                   +NF_D_Q_Q3_2D_Weight15[1]*PointValues[1]
                   +NF_D_Q_Q3_2D_Weight15[2]*PointValues[2]
                   +NF_D_Q_Q3_2D_Weight15[3]*PointValues[3]
                   +NF_D_Q_Q3_2D_Weight15[4]*PointValues[4]
                   +NF_D_Q_Q3_2D_Weight15[5]*PointValues[5]
                   +NF_D_Q_Q3_2D_Weight15[6]*PointValues[6]
                   +NF_D_Q_Q3_2D_Weight15[7]*PointValues[7]
                   +NF_D_Q_Q3_2D_Weight15[8]*PointValues[8]
                   +NF_D_Q_Q3_2D_Weight15[9]*PointValues[9]
                   +NF_D_Q_Q3_2D_Weight15[10]*PointValues[10]
                   +NF_D_Q_Q3_2D_Weight15[11]*PointValues[11]
                   +NF_D_Q_Q3_2D_Weight15[12]*PointValues[12]
                   +NF_D_Q_Q3_2D_Weight15[13]*PointValues[13]
                   +NF_D_Q_Q3_2D_Weight15[14]*PointValues[14]
                   +NF_D_Q_Q3_2D_Weight15[15]*PointValues[15];
}

void NF_D_Q_Q3_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
                           double *Functionals)
{
}

/*
    TNodalFunctional2D(NodalFunctional2D id,
                       int n_allfunctionals, int n_edgefunctionals,
                       int n_pointsall, int n_pointsedge,
                       double *xi, double *eta, double *t,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evaledge);
*/

TNodalFunctional2D *NF_D_Q_Q3_2D_Obj = new TNodalFunctional2D
        (NF_D_Q_Q3_2D, 16, 0, 16, 0, NF_D_Q_Q3_2D_Xi, NF_D_Q_Q3_2D_Eta,
         NF_D_Q_Q3_2D_T, NF_D_Q_Q3_2D_EvalAll, NULL);
