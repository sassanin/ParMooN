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
   
static double NF_D_Q_Q4_2D_Xi[25]={
   -0.9061798459386639928, -0.53846931010568309104, 0,
    0.53846931010568309104, 0.9061798459386639928,
   -0.9061798459386639928, -0.53846931010568309104, 0,
    0.53846931010568309104, 0.9061798459386639928,
   -0.9061798459386639928, -0.53846931010568309104, 0,
    0.53846931010568309104, 0.9061798459386639928,
   -0.9061798459386639928, -0.53846931010568309104, 0,
    0.53846931010568309104, 0.9061798459386639928,
   -0.9061798459386639928, -0.53846931010568309104, 0,
    0.53846931010568309104, 0.9061798459386639928 };
 
static double NF_D_Q_Q4_2D_Eta[25]={
   -0.9061798459386639928, -0.9061798459386639928, -0.9061798459386639928,
   -0.9061798459386639928, -0.9061798459386639928,
   -0.53846931010568309104, -0.53846931010568309104, -0.53846931010568309104,
   -0.53846931010568309104, -0.53846931010568309104,
    0, 0, 0,
    0, 0,
    0.53846931010568309104, 0.53846931010568309104, 0.53846931010568309104,
    0.53846931010568309104, 0.53846931010568309104,
    0.9061798459386639928, 0.9061798459386639928, 0.9061798459386639928,
    0.9061798459386639928, 0.9061798459386639928 };

static double NF_D_Q_Q4_2D_Weight0[25] = {
    0.014033587215607158989, 0.02835, 0.03369626809688022578,
    0.02835, 0.014033587215607158989,
    0.02835, 0.057271351055997779283, 0.068071633137687675455,
    0.057271351055997779283, 0.02835,
    0.03369626809688022578, 0.068071633137687675455, 0.080908641975308641975,
    0.068071633137687675455, 0.03369626809688022578,
    0.02835, 0.057271351055997779283, 0.068071633137687675455,
    0.057271351055997779283, 0.02835,
    0.014033587215607158989, 0.02835, 0.03369626809688022578,
    0.02835, 0.014033587215607158989 };

static double NF_D_Q_Q4_2D_Weight1[25] = {
   -0.038150861703017099767, -0.045796814824488346893, 0,
    0.045796814824488346893, 0.038150861703017099767,
   -0.077070595897083372587, -0.092516594675830526938, 0,
    0.092516594675830526938, 0.077070595897083372587,
   -0.091604637098216524605, -0.1099634560002535146, 0,
    0.1099634560002535146, 0.091604637098216524605,
   -0.077070595897083372587, -0.092516594675830526938, 0,
    0.092516594675830526938, 0.077070595897083372587,
   -0.038150861703017099767, -0.045796814824488346893, 0,
    0.045796814824488346893, 0.038150861703017099767 };

static double NF_D_Q_Q4_2D_Weight2[25] = {
    0.051344886912150381713, -0.0092245517910500994885, -0.08424067024220056445,
   -0.0092245517910500994885, 0.051344886912150381713,
    0.10372455179105009949, -0.01863501036894050517, -0.17017908284421918864,
   -0.01863501036894050517, 0.10372455179105009949,
    0.1232850195901221217, -0.022149217120986319232, -0.20227160493827160494,
   -0.022149217120986319232, 0.1232850195901221217,
    0.10372455179105009949, -0.01863501036894050517, -0.17017908284421918864,
   -0.01863501036894050517, 0.10372455179105009949,
    0.051344886912150381713, -0.0092245517910500994885, -0.08424067024220056445,
   -0.0092245517910500994885, 0.051344886912150381713 };

static double NF_D_Q_Q4_2D_Weight3[25] = {
   -0.049218852456151369207, 0.082829478484557283863, 0,
   -0.082829478484557283863, 0.049218852456151369207,
   -0.099429635893813178533, 0.16732825890914621298, 0,
   -0.16732825890914621298, 0.099429635893813178533,
   -0.11818016465090349137, 0.19888317010182638258, 0,
   -0.19888317010182638258, 0.11818016465090349137,
   -0.099429635893813178533, 0.16732825890914621298, 0,
   -0.16732825890914621298, 0.099429635893813178533,
   -0.049218852456151369207, 0.082829478484557283863, 0,
   -0.082829478484557283863, 0.049218852456151369207 };

static double NF_D_Q_Q4_2D_Weight4[25] = {
    0.031036949974581418315, -0.087899402388066799318, 0.11372490482697076201,
   -0.087899402388066799318, 0.031036949974581418315,
    0.062699402388066799318, -0.17757028330791475165, 0.22974176183969590466,
   -0.17757028330791475165, 0.062699402388066799318,
    0.074523311195148893522, -0.21105664452848222685, 0.27306666666666666667,
   -0.21105664452848222685, 0.074523311195148893522,
    0.062699402388066799318, -0.17757028330791475165, 0.22974176183969590466,
   -0.17757028330791475165, 0.062699402388066799318,
    0.031036949974581418315, -0.087899402388066799318, 0.11372490482697076201,
   -0.087899402388066799318, 0.031036949974581418315 };

static double NF_D_Q_Q4_2D_Weight5[25] = {
   -0.038150861703017099767, -0.077070595897083372587, -0.091604637098216524605,
   -0.077070595897083372587, -0.038150861703017099767,
   -0.045796814824488346893, -0.092516594675830526938, -0.1099634560002535146,
   -0.092516594675830526938, -0.045796814824488346893,
    0, 0, 0,
    0, 0,
    0.045796814824488346893, 0.092516594675830526938, 0.1099634560002535146,
    0.092516594675830526938, 0.045796814824488346893,
    0.038150861703017099767, 0.077070595897083372587, 0.091604637098216524605,
    0.077070595897083372587, 0.038150861703017099767 };

static double NF_D_Q_Q4_2D_Weight6[25] = {
    0.10371462594140193502, 0.12450045180640912036, 0,
   -0.12450045180640912036, -0.10371462594140193502,
    0.12450045180640912036, 0.14945204072526473164, 0,
   -0.14945204072526473164, -0.12450045180640912036,
    0, 0, 0,
    0, 0,
   -0.12450045180640912036, -0.14945204072526473164, 0,
    0.14945204072526473164, 0.12450045180640912036,
   -0.10371462594140193502, -0.12450045180640912036, 0,
    0.12450045180640912036, 0.10371462594140193502 };

static double NF_D_Q_Q4_2D_Weight7[25] = {
   -0.13958310513537167423, 0.025077308762601018473, 0.22901159274554131151,
    0.025077308762601018473, -0.13958310513537167423,
   -0.16755746353184382749, 0.030103143531526934242, 0.2749086400006337865,
    0.030103143531526934242, -0.16755746353184382749,
    0, 0, 0,
    0, 0,
    0.16755746353184382749, -0.030103143531526934242, -0.2749086400006337865,
   -0.030103143531526934242, 0.16755746353184382749,
    0.13958310513537167423, -0.025077308762601018473, -0.22901159274554131151,
   -0.025077308762601018473, 0.13958310513537167423 };

static double NF_D_Q_Q4_2D_Weight8[25] = {
    0.13380339640797924484, -0.22517521215694801014, 0,
    0.22517521215694801014, -0.13380339640797924484,
    0.16061942233140254033, -0.27030339640797924484, 0,
    0.27030339640797924484, -0.16061942233140254033,
    0, 0, 0,
    0, 0,
   -0.16061942233140254033, 0.27030339640797924484, 0,
   -0.27030339640797924484, 0.16061942233140254033,
   -0.13380339640797924484, 0.22517521215694801014, 0,
   -0.22517521215694801014, 0.13380339640797924484 };

static double NF_D_Q_Q4_2D_Weight9[25] = {
   -0.084375175639116632927, 0.2389580007423570182, -0.30916565020648077054,
    0.2389580007423570182, -0.084375175639116632927,
   -0.10128511184382284494, 0.28684844384425065082, -0.37112666400085561177,
    0.28684844384425065082, -0.10128511184382284494,
    0, 0, 0,
    0, 0,
    0.10128511184382284494, -0.28684844384425065082, 0.37112666400085561177,
   -0.28684844384425065082, 0.10128511184382284494,
    0.084375175639116632927, -0.2389580007423570182, 0.30916565020648077054,
   -0.2389580007423570182, 0.084375175639116632927 };

static double NF_D_Q_Q4_2D_Weight10[25] = {
    0.051344886912150381713, 0.10372455179105009949, 0.1232850195901221217,
    0.10372455179105009949, 0.051344886912150381713,
   -0.0092245517910500994885, -0.01863501036894050517, -0.022149217120986319232,
   -0.01863501036894050517, -0.0092245517910500994885,
   -0.08424067024220056445, -0.17017908284421918864, -0.20227160493827160494,
   -0.17017908284421918864, -0.08424067024220056445,
   -0.0092245517910500994885, -0.01863501036894050517, -0.022149217120986319232,
   -0.01863501036894050517, -0.0092245517910500994885,
    0.051344886912150381713, 0.10372455179105009949, 0.1232850195901221217,
    0.10372455179105009949, 0.051344886912150381713 };

static double NF_D_Q_Q4_2D_Weight11[25] = {
   -0.13958310513537167423, -0.16755746353184382749, 0,
    0.16755746353184382749, 0.13958310513537167423,
    0.025077308762601018473, 0.030103143531526934242, 0,
   -0.030103143531526934242, -0.025077308762601018473,
    0.22901159274554131151, 0.2749086400006337865, 0,
   -0.2749086400006337865, -0.22901159274554131151,
    0.025077308762601018473, 0.030103143531526934242, 0,
   -0.030103143531526934242, -0.025077308762601018473,
   -0.13958310513537167423, -0.16755746353184382749, 0,
    0.16755746353184382749, 0.13958310513537167423 };

static double NF_D_Q_Q4_2D_Weight12[25] = {
    0.18785627448765265213, -0.03375, -0.30821254897530530425,
   -0.03375, 0.18785627448765265213,
   -0.03375, 0.0060634785987671009604, 0.055373042802465798079,
    0.0060634785987671009604, -0.03375,
   -0.30821254897530530425, 0.055373042802465798079, 0.50567901234567901235,
    0.055373042802465798079, -0.30821254897530530425,
   -0.03375, 0.0060634785987671009604, 0.055373042802465798079,
    0.0060634785987671009604, -0.03375,
    0.18785627448765265213, -0.03375, -0.30821254897530530425,
   -0.03375, 0.18785627448765265213 };

static double NF_D_Q_Q4_2D_Weight13[25] = {
   -0.18007772171725312491, 0.3030494014425796084, 0,
   -0.3030494014425796084, 0.18007772171725312491,
    0.032352515903623760697, -0.054445438815296630174, 0,
    0.054445438815296630174, -0.032352515903623760697,
    0.29545041162725872843, -0.49720792525456595646, 0,
    0.49720792525456595646, -0.29545041162725872843,
    0.032352515903623760697, -0.054445438815296630174, 0,
    0.054445438815296630174, -0.032352515903623760697,
   -0.18007772171725312491, 0.3030494014425796084, 0,
   -0.3030494014425796084, 0.18007772171725312491 };

static double NF_D_Q_Q4_2D_Weight14[25] = {
    0.11355533421780251827, -0.32159880477613359864, 0.41608694111666216074,
   -0.32159880477613359864, 0.11355533421780251827,
   -0.020401195223866401364, 0.057777999115530815067, -0.074753607783328827407,
    0.057777999115530815067, -0.020401195223866401364,
   -0.1863082779878722338, 0.52764161132120556714, -0.68266666666666666667,
    0.52764161132120556714, -0.1863082779878722338,
   -0.020401195223866401364, 0.057777999115530815067, -0.074753607783328827407,
    0.057777999115530815067, -0.020401195223866401364,
    0.11355533421780251827, -0.32159880477613359864, 0.41608694111666216074,
   -0.32159880477613359864, 0.11355533421780251827 };

static double NF_D_Q_Q4_2D_Weight15[25] = {
   -0.049218852456151369207, -0.099429635893813178533, -0.11818016465090349137,
   -0.099429635893813178533, -0.049218852456151369207,
    0.082829478484557283863, 0.16732825890914621298, 0.19888317010182638258,
    0.16732825890914621298, 0.082829478484557283863,
    0, 0, 0,
    0, 0,
   -0.082829478484557283863, -0.16732825890914621298, -0.19888317010182638258,
   -0.16732825890914621298, -0.082829478484557283863,
    0.049218852456151369207, 0.099429635893813178533, 0.11818016465090349137,
    0.099429635893813178533, 0.049218852456151369207 };

static double NF_D_Q_Q4_2D_Weight16[25] = {
    0.13380339640797924484, 0.16061942233140254033, 0,
   -0.16061942233140254033, -0.13380339640797924484,
   -0.22517521215694801014, -0.27030339640797924484, 0,
    0.27030339640797924484, 0.22517521215694801014,
    0, 0, 0,
    0, 0,
    0.22517521215694801014, 0.27030339640797924484, 0,
   -0.27030339640797924484, -0.22517521215694801014,
   -0.13380339640797924484, -0.16061942233140254033, 0,
    0.16061942233140254033, 0.13380339640797924484 };

static double NF_D_Q_Q4_2D_Weight17[25] = {
   -0.18007772171725312491, 0.032352515903623760697, 0.29545041162725872843,
    0.032352515903623760697, -0.18007772171725312491,
    0.3030494014425796084, -0.054445438815296630174, -0.49720792525456595646,
   -0.054445438815296630174, 0.3030494014425796084,
    0, 0, 0,
    0, 0,
   -0.3030494014425796084, 0.054445438815296630174, 0.49720792525456595646,
    0.054445438815296630174, -0.3030494014425796084,
    0.18007772171725312491, -0.032352515903623760697, -0.29545041162725872843,
   -0.032352515903623760697, 0.18007772171725312491 };

static double NF_D_Q_Q4_2D_Weight18[25] = {
    0.17262125498505972143, -0.29050105421495461416, 0,
    0.29050105421495461416, -0.17262125498505972143,
   -0.29050105421495461416, 0.48887874501494027857, 0,
   -0.48887874501494027857, 0.29050105421495461416,
    0, 0, 0,
    0, 0,
    0.29050105421495461416, -0.48887874501494027857, 0,
    0.48887874501494027857, -0.29050105421495461416,
   -0.17262125498505972143, 0.29050105421495461416, 0,
   -0.29050105421495461416, 0.17262125498505972143 };

static double NF_D_Q_Q4_2D_Weight19[25] = {
   -0.10885335573993384866, 0.30828238358833349035, -0.39885805569679928338,
    0.30828238358833349035, -0.10885335573993384866,
    0.18318725929795338757, -0.51880260884478540818, 0.67123069909366404122,
   -0.51880260884478540818, 0.18318725929795338757,
    0, 0, 0,
    0, 0,
   -0.18318725929795338757, 0.51880260884478540818, -0.67123069909366404122,
    0.51880260884478540818, -0.18318725929795338757,
    0.10885335573993384866, -0.30828238358833349035, 0.39885805569679928338,
   -0.30828238358833349035, 0.10885335573993384866 };

static double NF_D_Q_Q4_2D_Weight20[25] = {
    0.031036949974581418315, 0.062699402388066799318, 0.074523311195148893522,
    0.062699402388066799318, 0.031036949974581418315,
   -0.087899402388066799318, -0.17757028330791475165, -0.21105664452848222685,
   -0.17757028330791475165, -0.087899402388066799318,
    0.11372490482697076201, 0.22974176183969590466, 0.27306666666666666667,
    0.22974176183969590466, 0.11372490482697076201,
   -0.087899402388066799318, -0.17757028330791475165, -0.21105664452848222685,
   -0.17757028330791475165, -0.087899402388066799318,
    0.031036949974581418315, 0.062699402388066799318, 0.074523311195148893522,
    0.062699402388066799318, 0.031036949974581418315 };

static double NF_D_Q_Q4_2D_Weight21[25] = {
   -0.084375175639116632927, -0.10128511184382284494, 0,
    0.10128511184382284494, 0.084375175639116632927,
    0.2389580007423570182, 0.28684844384425065082, 0,
   -0.28684844384425065082, -0.2389580007423570182,
   -0.30916565020648077054, -0.37112666400085561177, 0,
    0.37112666400085561177, 0.30916565020648077054,
    0.2389580007423570182, 0.28684844384425065082, 0,
   -0.28684844384425065082, -0.2389580007423570182,
   -0.084375175639116632927, -0.10128511184382284494, 0,
    0.10128511184382284494, 0.084375175639116632927 };

static double NF_D_Q_Q4_2D_Weight22[25] = {
    0.11355533421780251827, -0.020401195223866401364, -0.1863082779878722338,
   -0.020401195223866401364, 0.11355533421780251827,
   -0.32159880477613359864, 0.057777999115530815067, 0.52764161132120556714,
    0.057777999115530815067, -0.32159880477613359864,
    0.41608694111666216074, -0.074753607783328827407, -0.68266666666666666667,
   -0.074753607783328827407, 0.41608694111666216074,
   -0.32159880477613359864, 0.057777999115530815067, 0.52764161132120556714,
    0.057777999115530815067, -0.32159880477613359864,
    0.11355533421780251827, -0.020401195223866401364, -0.1863082779878722338,
   -0.020401195223866401364, 0.11355533421780251827 };

static double NF_D_Q_Q4_2D_Weight23[25] = {
   -0.10885335573993384866, 0.18318725929795338757, 0,
   -0.18318725929795338757, 0.10885335573993384866,
    0.30828238358833349035, -0.51880260884478540818, 0,
    0.51880260884478540818, -0.30828238358833349035,
   -0.39885805569679928338, 0.67123069909366404122, 0,
   -0.67123069909366404122, 0.39885805569679928338,
    0.30828238358833349035, -0.51880260884478540818, 0,
    0.51880260884478540818, -0.30828238358833349035,
   -0.10885335573993384866, 0.18318725929795338757, 0,
   -0.18318725929795338757, 0.10885335573993384866 };

static double NF_D_Q_Q4_2D_Weight24[25] = {
    0.068641912358186242182, -0.1944, 0.25151617528362751564,
   -0.1944, 0.068641912358186242182,
   -0.1944, 0.55055808764181375782, -0.71231617528362751564,
    0.55055808764181375782, -0.1944,
    0.25151617528362751564, -0.71231617528362751564, 0.9216,
   -0.71231617528362751564, 0.25151617528362751564,
   -0.1944, 0.55055808764181375782, -0.71231617528362751564,
    0.55055808764181375782, -0.1944,
    0.068641912358186242182, -0.1944, 0.25151617528362751564,
   -0.1944, 0.068641912358186242182 };

static double *NF_D_Q_Q4_2D_T = NULL;

void NF_D_Q_Q4_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  Functionals[0] =  NF_D_Q_Q4_2D_Weight0[0]*PointValues[0]
                   +NF_D_Q_Q4_2D_Weight0[1]*PointValues[1]
                   +NF_D_Q_Q4_2D_Weight0[2]*PointValues[2]
                   +NF_D_Q_Q4_2D_Weight0[3]*PointValues[3]
                   +NF_D_Q_Q4_2D_Weight0[4]*PointValues[4]
                   +NF_D_Q_Q4_2D_Weight0[5]*PointValues[5]
                   +NF_D_Q_Q4_2D_Weight0[6]*PointValues[6]
                   +NF_D_Q_Q4_2D_Weight0[7]*PointValues[7]
                   +NF_D_Q_Q4_2D_Weight0[8]*PointValues[8]
                   +NF_D_Q_Q4_2D_Weight0[9]*PointValues[9]
                   +NF_D_Q_Q4_2D_Weight0[10]*PointValues[10]
                   +NF_D_Q_Q4_2D_Weight0[11]*PointValues[11]
                   +NF_D_Q_Q4_2D_Weight0[12]*PointValues[12]
                   +NF_D_Q_Q4_2D_Weight0[13]*PointValues[13]
                   +NF_D_Q_Q4_2D_Weight0[14]*PointValues[14]
                   +NF_D_Q_Q4_2D_Weight0[15]*PointValues[15]
                   +NF_D_Q_Q4_2D_Weight0[16]*PointValues[16]
                   +NF_D_Q_Q4_2D_Weight0[17]*PointValues[17]
                   +NF_D_Q_Q4_2D_Weight0[18]*PointValues[18]
                   +NF_D_Q_Q4_2D_Weight0[19]*PointValues[19]
                   +NF_D_Q_Q4_2D_Weight0[20]*PointValues[20]
                   +NF_D_Q_Q4_2D_Weight0[21]*PointValues[21]
                   +NF_D_Q_Q4_2D_Weight0[22]*PointValues[22]
                   +NF_D_Q_Q4_2D_Weight0[23]*PointValues[23]
                   +NF_D_Q_Q4_2D_Weight0[24]*PointValues[24];

  Functionals[1] =  NF_D_Q_Q4_2D_Weight1[0]*PointValues[0]
                   +NF_D_Q_Q4_2D_Weight1[1]*PointValues[1]
                   +NF_D_Q_Q4_2D_Weight1[2]*PointValues[2]
                   +NF_D_Q_Q4_2D_Weight1[3]*PointValues[3]
                   +NF_D_Q_Q4_2D_Weight1[4]*PointValues[4]
                   +NF_D_Q_Q4_2D_Weight1[5]*PointValues[5]
                   +NF_D_Q_Q4_2D_Weight1[6]*PointValues[6]
                   +NF_D_Q_Q4_2D_Weight1[7]*PointValues[7]
                   +NF_D_Q_Q4_2D_Weight1[8]*PointValues[8]
                   +NF_D_Q_Q4_2D_Weight1[9]*PointValues[9]
                   +NF_D_Q_Q4_2D_Weight1[10]*PointValues[10]
                   +NF_D_Q_Q4_2D_Weight1[11]*PointValues[11]
                   +NF_D_Q_Q4_2D_Weight1[12]*PointValues[12]
                   +NF_D_Q_Q4_2D_Weight1[13]*PointValues[13]
                   +NF_D_Q_Q4_2D_Weight1[14]*PointValues[14]
                   +NF_D_Q_Q4_2D_Weight1[15]*PointValues[15]
                   +NF_D_Q_Q4_2D_Weight1[16]*PointValues[16]
                   +NF_D_Q_Q4_2D_Weight1[17]*PointValues[17]
                   +NF_D_Q_Q4_2D_Weight1[18]*PointValues[18]
                   +NF_D_Q_Q4_2D_Weight1[19]*PointValues[19]
                   +NF_D_Q_Q4_2D_Weight1[20]*PointValues[20]
                   +NF_D_Q_Q4_2D_Weight1[21]*PointValues[21]
                   +NF_D_Q_Q4_2D_Weight1[22]*PointValues[22]
                   +NF_D_Q_Q4_2D_Weight1[23]*PointValues[23]
                   +NF_D_Q_Q4_2D_Weight1[24]*PointValues[24];

  Functionals[2] =  NF_D_Q_Q4_2D_Weight2[0]*PointValues[0]
                   +NF_D_Q_Q4_2D_Weight2[1]*PointValues[1]
                   +NF_D_Q_Q4_2D_Weight2[2]*PointValues[2]
                   +NF_D_Q_Q4_2D_Weight2[3]*PointValues[3]
                   +NF_D_Q_Q4_2D_Weight2[4]*PointValues[4]
                   +NF_D_Q_Q4_2D_Weight2[5]*PointValues[5]
                   +NF_D_Q_Q4_2D_Weight2[6]*PointValues[6]
                   +NF_D_Q_Q4_2D_Weight2[7]*PointValues[7]
                   +NF_D_Q_Q4_2D_Weight2[8]*PointValues[8]
                   +NF_D_Q_Q4_2D_Weight2[9]*PointValues[9]
                   +NF_D_Q_Q4_2D_Weight2[10]*PointValues[10]
                   +NF_D_Q_Q4_2D_Weight2[11]*PointValues[11]
                   +NF_D_Q_Q4_2D_Weight2[12]*PointValues[12]
                   +NF_D_Q_Q4_2D_Weight2[13]*PointValues[13]
                   +NF_D_Q_Q4_2D_Weight2[14]*PointValues[14]
                   +NF_D_Q_Q4_2D_Weight2[15]*PointValues[15]
                   +NF_D_Q_Q4_2D_Weight2[16]*PointValues[16]
                   +NF_D_Q_Q4_2D_Weight2[17]*PointValues[17]
                   +NF_D_Q_Q4_2D_Weight2[18]*PointValues[18]
                   +NF_D_Q_Q4_2D_Weight2[19]*PointValues[19]
                   +NF_D_Q_Q4_2D_Weight2[20]*PointValues[20]
                   +NF_D_Q_Q4_2D_Weight2[21]*PointValues[21]
                   +NF_D_Q_Q4_2D_Weight2[22]*PointValues[22]
                   +NF_D_Q_Q4_2D_Weight2[23]*PointValues[23]
                   +NF_D_Q_Q4_2D_Weight2[24]*PointValues[24];

  Functionals[3] =  NF_D_Q_Q4_2D_Weight3[0]*PointValues[0]
                   +NF_D_Q_Q4_2D_Weight3[1]*PointValues[1]
                   +NF_D_Q_Q4_2D_Weight3[2]*PointValues[2]
                   +NF_D_Q_Q4_2D_Weight3[3]*PointValues[3]
                   +NF_D_Q_Q4_2D_Weight3[4]*PointValues[4]
                   +NF_D_Q_Q4_2D_Weight3[5]*PointValues[5]
                   +NF_D_Q_Q4_2D_Weight3[6]*PointValues[6]
                   +NF_D_Q_Q4_2D_Weight3[7]*PointValues[7]
                   +NF_D_Q_Q4_2D_Weight3[8]*PointValues[8]
                   +NF_D_Q_Q4_2D_Weight3[9]*PointValues[9]
                   +NF_D_Q_Q4_2D_Weight3[10]*PointValues[10]
                   +NF_D_Q_Q4_2D_Weight3[11]*PointValues[11]
                   +NF_D_Q_Q4_2D_Weight3[12]*PointValues[12]
                   +NF_D_Q_Q4_2D_Weight3[13]*PointValues[13]
                   +NF_D_Q_Q4_2D_Weight3[14]*PointValues[14]
                   +NF_D_Q_Q4_2D_Weight3[15]*PointValues[15]
                   +NF_D_Q_Q4_2D_Weight3[16]*PointValues[16]
                   +NF_D_Q_Q4_2D_Weight3[17]*PointValues[17]
                   +NF_D_Q_Q4_2D_Weight3[18]*PointValues[18]
                   +NF_D_Q_Q4_2D_Weight3[19]*PointValues[19]
                   +NF_D_Q_Q4_2D_Weight3[20]*PointValues[20]
                   +NF_D_Q_Q4_2D_Weight3[21]*PointValues[21]
                   +NF_D_Q_Q4_2D_Weight3[22]*PointValues[22]
                   +NF_D_Q_Q4_2D_Weight3[23]*PointValues[23]
                   +NF_D_Q_Q4_2D_Weight3[24]*PointValues[24];

  Functionals[4] =  NF_D_Q_Q4_2D_Weight4[0]*PointValues[0]
                   +NF_D_Q_Q4_2D_Weight4[1]*PointValues[1]
                   +NF_D_Q_Q4_2D_Weight4[2]*PointValues[2]
                   +NF_D_Q_Q4_2D_Weight4[3]*PointValues[3]
                   +NF_D_Q_Q4_2D_Weight4[4]*PointValues[4]
                   +NF_D_Q_Q4_2D_Weight4[5]*PointValues[5]
                   +NF_D_Q_Q4_2D_Weight4[6]*PointValues[6]
                   +NF_D_Q_Q4_2D_Weight4[7]*PointValues[7]
                   +NF_D_Q_Q4_2D_Weight4[8]*PointValues[8]
                   +NF_D_Q_Q4_2D_Weight4[9]*PointValues[9]
                   +NF_D_Q_Q4_2D_Weight4[10]*PointValues[10]
                   +NF_D_Q_Q4_2D_Weight4[11]*PointValues[11]
                   +NF_D_Q_Q4_2D_Weight4[12]*PointValues[12]
                   +NF_D_Q_Q4_2D_Weight4[13]*PointValues[13]
                   +NF_D_Q_Q4_2D_Weight4[14]*PointValues[14]
                   +NF_D_Q_Q4_2D_Weight4[15]*PointValues[15]
                   +NF_D_Q_Q4_2D_Weight4[16]*PointValues[16]
                   +NF_D_Q_Q4_2D_Weight4[17]*PointValues[17]
                   +NF_D_Q_Q4_2D_Weight4[18]*PointValues[18]
                   +NF_D_Q_Q4_2D_Weight4[19]*PointValues[19]
                   +NF_D_Q_Q4_2D_Weight4[20]*PointValues[20]
                   +NF_D_Q_Q4_2D_Weight4[21]*PointValues[21]
                   +NF_D_Q_Q4_2D_Weight4[22]*PointValues[22]
                   +NF_D_Q_Q4_2D_Weight4[23]*PointValues[23]
                   +NF_D_Q_Q4_2D_Weight4[24]*PointValues[24];

  Functionals[5] =  NF_D_Q_Q4_2D_Weight5[0]*PointValues[0]
                   +NF_D_Q_Q4_2D_Weight5[1]*PointValues[1]
                   +NF_D_Q_Q4_2D_Weight5[2]*PointValues[2]
                   +NF_D_Q_Q4_2D_Weight5[3]*PointValues[3]
                   +NF_D_Q_Q4_2D_Weight5[4]*PointValues[4]
                   +NF_D_Q_Q4_2D_Weight5[5]*PointValues[5]
                   +NF_D_Q_Q4_2D_Weight5[6]*PointValues[6]
                   +NF_D_Q_Q4_2D_Weight5[7]*PointValues[7]
                   +NF_D_Q_Q4_2D_Weight5[8]*PointValues[8]
                   +NF_D_Q_Q4_2D_Weight5[9]*PointValues[9]
                   +NF_D_Q_Q4_2D_Weight5[10]*PointValues[10]
                   +NF_D_Q_Q4_2D_Weight5[11]*PointValues[11]
                   +NF_D_Q_Q4_2D_Weight5[12]*PointValues[12]
                   +NF_D_Q_Q4_2D_Weight5[13]*PointValues[13]
                   +NF_D_Q_Q4_2D_Weight5[14]*PointValues[14]
                   +NF_D_Q_Q4_2D_Weight5[15]*PointValues[15]
                   +NF_D_Q_Q4_2D_Weight5[16]*PointValues[16]
                   +NF_D_Q_Q4_2D_Weight5[17]*PointValues[17]
                   +NF_D_Q_Q4_2D_Weight5[18]*PointValues[18]
                   +NF_D_Q_Q4_2D_Weight5[19]*PointValues[19]
                   +NF_D_Q_Q4_2D_Weight5[20]*PointValues[20]
                   +NF_D_Q_Q4_2D_Weight5[21]*PointValues[21]
                   +NF_D_Q_Q4_2D_Weight5[22]*PointValues[22]
                   +NF_D_Q_Q4_2D_Weight5[23]*PointValues[23]
                   +NF_D_Q_Q4_2D_Weight5[24]*PointValues[24];

  Functionals[6] =  NF_D_Q_Q4_2D_Weight6[0]*PointValues[0]
                   +NF_D_Q_Q4_2D_Weight6[1]*PointValues[1]
                   +NF_D_Q_Q4_2D_Weight6[2]*PointValues[2]
                   +NF_D_Q_Q4_2D_Weight6[3]*PointValues[3]
                   +NF_D_Q_Q4_2D_Weight6[4]*PointValues[4]
                   +NF_D_Q_Q4_2D_Weight6[5]*PointValues[5]
                   +NF_D_Q_Q4_2D_Weight6[6]*PointValues[6]
                   +NF_D_Q_Q4_2D_Weight6[7]*PointValues[7]
                   +NF_D_Q_Q4_2D_Weight6[8]*PointValues[8]
                   +NF_D_Q_Q4_2D_Weight6[9]*PointValues[9]
                   +NF_D_Q_Q4_2D_Weight6[10]*PointValues[10]
                   +NF_D_Q_Q4_2D_Weight6[11]*PointValues[11]
                   +NF_D_Q_Q4_2D_Weight6[12]*PointValues[12]
                   +NF_D_Q_Q4_2D_Weight6[13]*PointValues[13]
                   +NF_D_Q_Q4_2D_Weight6[14]*PointValues[14]
                   +NF_D_Q_Q4_2D_Weight6[15]*PointValues[15]
                   +NF_D_Q_Q4_2D_Weight6[16]*PointValues[16]
                   +NF_D_Q_Q4_2D_Weight6[17]*PointValues[17]
                   +NF_D_Q_Q4_2D_Weight6[18]*PointValues[18]
                   +NF_D_Q_Q4_2D_Weight6[19]*PointValues[19]
                   +NF_D_Q_Q4_2D_Weight6[20]*PointValues[20]
                   +NF_D_Q_Q4_2D_Weight6[21]*PointValues[21]
                   +NF_D_Q_Q4_2D_Weight6[22]*PointValues[22]
                   +NF_D_Q_Q4_2D_Weight6[23]*PointValues[23]
                   +NF_D_Q_Q4_2D_Weight6[24]*PointValues[24];

  Functionals[7] =  NF_D_Q_Q4_2D_Weight7[0]*PointValues[0]
                   +NF_D_Q_Q4_2D_Weight7[1]*PointValues[1]
                   +NF_D_Q_Q4_2D_Weight7[2]*PointValues[2]
                   +NF_D_Q_Q4_2D_Weight7[3]*PointValues[3]
                   +NF_D_Q_Q4_2D_Weight7[4]*PointValues[4]
                   +NF_D_Q_Q4_2D_Weight7[5]*PointValues[5]
                   +NF_D_Q_Q4_2D_Weight7[6]*PointValues[6]
                   +NF_D_Q_Q4_2D_Weight7[7]*PointValues[7]
                   +NF_D_Q_Q4_2D_Weight7[8]*PointValues[8]
                   +NF_D_Q_Q4_2D_Weight7[9]*PointValues[9]
                   +NF_D_Q_Q4_2D_Weight7[10]*PointValues[10]
                   +NF_D_Q_Q4_2D_Weight7[11]*PointValues[11]
                   +NF_D_Q_Q4_2D_Weight7[12]*PointValues[12]
                   +NF_D_Q_Q4_2D_Weight7[13]*PointValues[13]
                   +NF_D_Q_Q4_2D_Weight7[14]*PointValues[14]
                   +NF_D_Q_Q4_2D_Weight7[15]*PointValues[15]
                   +NF_D_Q_Q4_2D_Weight7[16]*PointValues[16]
                   +NF_D_Q_Q4_2D_Weight7[17]*PointValues[17]
                   +NF_D_Q_Q4_2D_Weight7[18]*PointValues[18]
                   +NF_D_Q_Q4_2D_Weight7[19]*PointValues[19]
                   +NF_D_Q_Q4_2D_Weight7[20]*PointValues[20]
                   +NF_D_Q_Q4_2D_Weight7[21]*PointValues[21]
                   +NF_D_Q_Q4_2D_Weight7[22]*PointValues[22]
                   +NF_D_Q_Q4_2D_Weight7[23]*PointValues[23]
                   +NF_D_Q_Q4_2D_Weight7[24]*PointValues[24];

  Functionals[8] =  NF_D_Q_Q4_2D_Weight8[0]*PointValues[0]
                   +NF_D_Q_Q4_2D_Weight8[1]*PointValues[1]
                   +NF_D_Q_Q4_2D_Weight8[2]*PointValues[2]
                   +NF_D_Q_Q4_2D_Weight8[3]*PointValues[3]
                   +NF_D_Q_Q4_2D_Weight8[4]*PointValues[4]
                   +NF_D_Q_Q4_2D_Weight8[5]*PointValues[5]
                   +NF_D_Q_Q4_2D_Weight8[6]*PointValues[6]
                   +NF_D_Q_Q4_2D_Weight8[7]*PointValues[7]
                   +NF_D_Q_Q4_2D_Weight8[8]*PointValues[8]
                   +NF_D_Q_Q4_2D_Weight8[9]*PointValues[9]
                   +NF_D_Q_Q4_2D_Weight8[10]*PointValues[10]
                   +NF_D_Q_Q4_2D_Weight8[11]*PointValues[11]
                   +NF_D_Q_Q4_2D_Weight8[12]*PointValues[12]
                   +NF_D_Q_Q4_2D_Weight8[13]*PointValues[13]
                   +NF_D_Q_Q4_2D_Weight8[14]*PointValues[14]
                   +NF_D_Q_Q4_2D_Weight8[15]*PointValues[15]
                   +NF_D_Q_Q4_2D_Weight8[16]*PointValues[16]
                   +NF_D_Q_Q4_2D_Weight8[17]*PointValues[17]
                   +NF_D_Q_Q4_2D_Weight8[18]*PointValues[18]
                   +NF_D_Q_Q4_2D_Weight8[19]*PointValues[19]
                   +NF_D_Q_Q4_2D_Weight8[20]*PointValues[20]
                   +NF_D_Q_Q4_2D_Weight8[21]*PointValues[21]
                   +NF_D_Q_Q4_2D_Weight8[22]*PointValues[22]
                   +NF_D_Q_Q4_2D_Weight8[23]*PointValues[23]
                   +NF_D_Q_Q4_2D_Weight8[24]*PointValues[24];

  Functionals[9] =  NF_D_Q_Q4_2D_Weight9[0]*PointValues[0]
                   +NF_D_Q_Q4_2D_Weight9[1]*PointValues[1]
                   +NF_D_Q_Q4_2D_Weight9[2]*PointValues[2]
                   +NF_D_Q_Q4_2D_Weight9[3]*PointValues[3]
                   +NF_D_Q_Q4_2D_Weight9[4]*PointValues[4]
                   +NF_D_Q_Q4_2D_Weight9[5]*PointValues[5]
                   +NF_D_Q_Q4_2D_Weight9[6]*PointValues[6]
                   +NF_D_Q_Q4_2D_Weight9[7]*PointValues[7]
                   +NF_D_Q_Q4_2D_Weight9[8]*PointValues[8]
                   +NF_D_Q_Q4_2D_Weight9[9]*PointValues[9]
                   +NF_D_Q_Q4_2D_Weight9[10]*PointValues[10]
                   +NF_D_Q_Q4_2D_Weight9[11]*PointValues[11]
                   +NF_D_Q_Q4_2D_Weight9[12]*PointValues[12]
                   +NF_D_Q_Q4_2D_Weight9[13]*PointValues[13]
                   +NF_D_Q_Q4_2D_Weight9[14]*PointValues[14]
                   +NF_D_Q_Q4_2D_Weight9[15]*PointValues[15]
                   +NF_D_Q_Q4_2D_Weight9[16]*PointValues[16]
                   +NF_D_Q_Q4_2D_Weight9[17]*PointValues[17]
                   +NF_D_Q_Q4_2D_Weight9[18]*PointValues[18]
                   +NF_D_Q_Q4_2D_Weight9[19]*PointValues[19]
                   +NF_D_Q_Q4_2D_Weight9[20]*PointValues[20]
                   +NF_D_Q_Q4_2D_Weight9[21]*PointValues[21]
                   +NF_D_Q_Q4_2D_Weight9[22]*PointValues[22]
                   +NF_D_Q_Q4_2D_Weight9[23]*PointValues[23]
                   +NF_D_Q_Q4_2D_Weight9[24]*PointValues[24];

  Functionals[10] =  NF_D_Q_Q4_2D_Weight10[0]*PointValues[0]
                   +NF_D_Q_Q4_2D_Weight10[1]*PointValues[1]
                   +NF_D_Q_Q4_2D_Weight10[2]*PointValues[2]
                   +NF_D_Q_Q4_2D_Weight10[3]*PointValues[3]
                   +NF_D_Q_Q4_2D_Weight10[4]*PointValues[4]
                   +NF_D_Q_Q4_2D_Weight10[5]*PointValues[5]
                   +NF_D_Q_Q4_2D_Weight10[6]*PointValues[6]
                   +NF_D_Q_Q4_2D_Weight10[7]*PointValues[7]
                   +NF_D_Q_Q4_2D_Weight10[8]*PointValues[8]
                   +NF_D_Q_Q4_2D_Weight10[9]*PointValues[9]
                   +NF_D_Q_Q4_2D_Weight10[10]*PointValues[10]
                   +NF_D_Q_Q4_2D_Weight10[11]*PointValues[11]
                   +NF_D_Q_Q4_2D_Weight10[12]*PointValues[12]
                   +NF_D_Q_Q4_2D_Weight10[13]*PointValues[13]
                   +NF_D_Q_Q4_2D_Weight10[14]*PointValues[14]
                   +NF_D_Q_Q4_2D_Weight10[15]*PointValues[15]
                   +NF_D_Q_Q4_2D_Weight10[16]*PointValues[16]
                   +NF_D_Q_Q4_2D_Weight10[17]*PointValues[17]
                   +NF_D_Q_Q4_2D_Weight10[18]*PointValues[18]
                   +NF_D_Q_Q4_2D_Weight10[19]*PointValues[19]
                   +NF_D_Q_Q4_2D_Weight10[20]*PointValues[20]
                   +NF_D_Q_Q4_2D_Weight10[21]*PointValues[21]
                   +NF_D_Q_Q4_2D_Weight10[22]*PointValues[22]
                   +NF_D_Q_Q4_2D_Weight10[23]*PointValues[23]
                   +NF_D_Q_Q4_2D_Weight10[24]*PointValues[24];

  Functionals[11] =  NF_D_Q_Q4_2D_Weight11[0]*PointValues[0]
                   +NF_D_Q_Q4_2D_Weight11[1]*PointValues[1]
                   +NF_D_Q_Q4_2D_Weight11[2]*PointValues[2]
                   +NF_D_Q_Q4_2D_Weight11[3]*PointValues[3]
                   +NF_D_Q_Q4_2D_Weight11[4]*PointValues[4]
                   +NF_D_Q_Q4_2D_Weight11[5]*PointValues[5]
                   +NF_D_Q_Q4_2D_Weight11[6]*PointValues[6]
                   +NF_D_Q_Q4_2D_Weight11[7]*PointValues[7]
                   +NF_D_Q_Q4_2D_Weight11[8]*PointValues[8]
                   +NF_D_Q_Q4_2D_Weight11[9]*PointValues[9]
                   +NF_D_Q_Q4_2D_Weight11[10]*PointValues[10]
                   +NF_D_Q_Q4_2D_Weight11[11]*PointValues[11]
                   +NF_D_Q_Q4_2D_Weight11[12]*PointValues[12]
                   +NF_D_Q_Q4_2D_Weight11[13]*PointValues[13]
                   +NF_D_Q_Q4_2D_Weight11[14]*PointValues[14]
                   +NF_D_Q_Q4_2D_Weight11[15]*PointValues[15]
                   +NF_D_Q_Q4_2D_Weight11[16]*PointValues[16]
                   +NF_D_Q_Q4_2D_Weight11[17]*PointValues[17]
                   +NF_D_Q_Q4_2D_Weight11[18]*PointValues[18]
                   +NF_D_Q_Q4_2D_Weight11[19]*PointValues[19]
                   +NF_D_Q_Q4_2D_Weight11[20]*PointValues[20]
                   +NF_D_Q_Q4_2D_Weight11[21]*PointValues[21]
                   +NF_D_Q_Q4_2D_Weight11[22]*PointValues[22]
                   +NF_D_Q_Q4_2D_Weight11[23]*PointValues[23]
                   +NF_D_Q_Q4_2D_Weight11[24]*PointValues[24];

  Functionals[12] =  NF_D_Q_Q4_2D_Weight12[0]*PointValues[0]
                   +NF_D_Q_Q4_2D_Weight12[1]*PointValues[1]
                   +NF_D_Q_Q4_2D_Weight12[2]*PointValues[2]
                   +NF_D_Q_Q4_2D_Weight12[3]*PointValues[3]
                   +NF_D_Q_Q4_2D_Weight12[4]*PointValues[4]
                   +NF_D_Q_Q4_2D_Weight12[5]*PointValues[5]
                   +NF_D_Q_Q4_2D_Weight12[6]*PointValues[6]
                   +NF_D_Q_Q4_2D_Weight12[7]*PointValues[7]
                   +NF_D_Q_Q4_2D_Weight12[8]*PointValues[8]
                   +NF_D_Q_Q4_2D_Weight12[9]*PointValues[9]
                   +NF_D_Q_Q4_2D_Weight12[10]*PointValues[10]
                   +NF_D_Q_Q4_2D_Weight12[11]*PointValues[11]
                   +NF_D_Q_Q4_2D_Weight12[12]*PointValues[12]
                   +NF_D_Q_Q4_2D_Weight12[13]*PointValues[13]
                   +NF_D_Q_Q4_2D_Weight12[14]*PointValues[14]
                   +NF_D_Q_Q4_2D_Weight12[15]*PointValues[15]
                   +NF_D_Q_Q4_2D_Weight12[16]*PointValues[16]
                   +NF_D_Q_Q4_2D_Weight12[17]*PointValues[17]
                   +NF_D_Q_Q4_2D_Weight12[18]*PointValues[18]
                   +NF_D_Q_Q4_2D_Weight12[19]*PointValues[19]
                   +NF_D_Q_Q4_2D_Weight12[20]*PointValues[20]
                   +NF_D_Q_Q4_2D_Weight12[21]*PointValues[21]
                   +NF_D_Q_Q4_2D_Weight12[22]*PointValues[22]
                   +NF_D_Q_Q4_2D_Weight12[23]*PointValues[23]
                   +NF_D_Q_Q4_2D_Weight12[24]*PointValues[24];

  Functionals[13] =  NF_D_Q_Q4_2D_Weight13[0]*PointValues[0]
                   +NF_D_Q_Q4_2D_Weight13[1]*PointValues[1]
                   +NF_D_Q_Q4_2D_Weight13[2]*PointValues[2]
                   +NF_D_Q_Q4_2D_Weight13[3]*PointValues[3]
                   +NF_D_Q_Q4_2D_Weight13[4]*PointValues[4]
                   +NF_D_Q_Q4_2D_Weight13[5]*PointValues[5]
                   +NF_D_Q_Q4_2D_Weight13[6]*PointValues[6]
                   +NF_D_Q_Q4_2D_Weight13[7]*PointValues[7]
                   +NF_D_Q_Q4_2D_Weight13[8]*PointValues[8]
                   +NF_D_Q_Q4_2D_Weight13[9]*PointValues[9]
                   +NF_D_Q_Q4_2D_Weight13[10]*PointValues[10]
                   +NF_D_Q_Q4_2D_Weight13[11]*PointValues[11]
                   +NF_D_Q_Q4_2D_Weight13[12]*PointValues[12]
                   +NF_D_Q_Q4_2D_Weight13[13]*PointValues[13]
                   +NF_D_Q_Q4_2D_Weight13[14]*PointValues[14]
                   +NF_D_Q_Q4_2D_Weight13[15]*PointValues[15]
                   +NF_D_Q_Q4_2D_Weight13[16]*PointValues[16]
                   +NF_D_Q_Q4_2D_Weight13[17]*PointValues[17]
                   +NF_D_Q_Q4_2D_Weight13[18]*PointValues[18]
                   +NF_D_Q_Q4_2D_Weight13[19]*PointValues[19]
                   +NF_D_Q_Q4_2D_Weight13[20]*PointValues[20]
                   +NF_D_Q_Q4_2D_Weight13[21]*PointValues[21]
                   +NF_D_Q_Q4_2D_Weight13[22]*PointValues[22]
                   +NF_D_Q_Q4_2D_Weight13[23]*PointValues[23]
                   +NF_D_Q_Q4_2D_Weight13[24]*PointValues[24];

  Functionals[14] =  NF_D_Q_Q4_2D_Weight14[0]*PointValues[0]
                   +NF_D_Q_Q4_2D_Weight14[1]*PointValues[1]
                   +NF_D_Q_Q4_2D_Weight14[2]*PointValues[2]
                   +NF_D_Q_Q4_2D_Weight14[3]*PointValues[3]
                   +NF_D_Q_Q4_2D_Weight14[4]*PointValues[4]
                   +NF_D_Q_Q4_2D_Weight14[5]*PointValues[5]
                   +NF_D_Q_Q4_2D_Weight14[6]*PointValues[6]
                   +NF_D_Q_Q4_2D_Weight14[7]*PointValues[7]
                   +NF_D_Q_Q4_2D_Weight14[8]*PointValues[8]
                   +NF_D_Q_Q4_2D_Weight14[9]*PointValues[9]
                   +NF_D_Q_Q4_2D_Weight14[10]*PointValues[10]
                   +NF_D_Q_Q4_2D_Weight14[11]*PointValues[11]
                   +NF_D_Q_Q4_2D_Weight14[12]*PointValues[12]
                   +NF_D_Q_Q4_2D_Weight14[13]*PointValues[13]
                   +NF_D_Q_Q4_2D_Weight14[14]*PointValues[14]
                   +NF_D_Q_Q4_2D_Weight14[15]*PointValues[15]
                   +NF_D_Q_Q4_2D_Weight14[16]*PointValues[16]
                   +NF_D_Q_Q4_2D_Weight14[17]*PointValues[17]
                   +NF_D_Q_Q4_2D_Weight14[18]*PointValues[18]
                   +NF_D_Q_Q4_2D_Weight14[19]*PointValues[19]
                   +NF_D_Q_Q4_2D_Weight14[20]*PointValues[20]
                   +NF_D_Q_Q4_2D_Weight14[21]*PointValues[21]
                   +NF_D_Q_Q4_2D_Weight14[22]*PointValues[22]
                   +NF_D_Q_Q4_2D_Weight14[23]*PointValues[23]
                   +NF_D_Q_Q4_2D_Weight14[24]*PointValues[24];

  Functionals[15] =  NF_D_Q_Q4_2D_Weight15[0]*PointValues[0]
                   +NF_D_Q_Q4_2D_Weight15[1]*PointValues[1]
                   +NF_D_Q_Q4_2D_Weight15[2]*PointValues[2]
                   +NF_D_Q_Q4_2D_Weight15[3]*PointValues[3]
                   +NF_D_Q_Q4_2D_Weight15[4]*PointValues[4]
                   +NF_D_Q_Q4_2D_Weight15[5]*PointValues[5]
                   +NF_D_Q_Q4_2D_Weight15[6]*PointValues[6]
                   +NF_D_Q_Q4_2D_Weight15[7]*PointValues[7]
                   +NF_D_Q_Q4_2D_Weight15[8]*PointValues[8]
                   +NF_D_Q_Q4_2D_Weight15[9]*PointValues[9]
                   +NF_D_Q_Q4_2D_Weight15[10]*PointValues[10]
                   +NF_D_Q_Q4_2D_Weight15[11]*PointValues[11]
                   +NF_D_Q_Q4_2D_Weight15[12]*PointValues[12]
                   +NF_D_Q_Q4_2D_Weight15[13]*PointValues[13]
                   +NF_D_Q_Q4_2D_Weight15[14]*PointValues[14]
                   +NF_D_Q_Q4_2D_Weight15[15]*PointValues[15]
                   +NF_D_Q_Q4_2D_Weight15[16]*PointValues[16]
                   +NF_D_Q_Q4_2D_Weight15[17]*PointValues[17]
                   +NF_D_Q_Q4_2D_Weight15[18]*PointValues[18]
                   +NF_D_Q_Q4_2D_Weight15[19]*PointValues[19]
                   +NF_D_Q_Q4_2D_Weight15[20]*PointValues[20]
                   +NF_D_Q_Q4_2D_Weight15[21]*PointValues[21]
                   +NF_D_Q_Q4_2D_Weight15[22]*PointValues[22]
                   +NF_D_Q_Q4_2D_Weight15[23]*PointValues[23]
                   +NF_D_Q_Q4_2D_Weight15[24]*PointValues[24];

  Functionals[16] =  NF_D_Q_Q4_2D_Weight16[0]*PointValues[0]
                   +NF_D_Q_Q4_2D_Weight16[1]*PointValues[1]
                   +NF_D_Q_Q4_2D_Weight16[2]*PointValues[2]
                   +NF_D_Q_Q4_2D_Weight16[3]*PointValues[3]
                   +NF_D_Q_Q4_2D_Weight16[4]*PointValues[4]
                   +NF_D_Q_Q4_2D_Weight16[5]*PointValues[5]
                   +NF_D_Q_Q4_2D_Weight16[6]*PointValues[6]
                   +NF_D_Q_Q4_2D_Weight16[7]*PointValues[7]
                   +NF_D_Q_Q4_2D_Weight16[8]*PointValues[8]
                   +NF_D_Q_Q4_2D_Weight16[9]*PointValues[9]
                   +NF_D_Q_Q4_2D_Weight16[10]*PointValues[10]
                   +NF_D_Q_Q4_2D_Weight16[11]*PointValues[11]
                   +NF_D_Q_Q4_2D_Weight16[12]*PointValues[12]
                   +NF_D_Q_Q4_2D_Weight16[13]*PointValues[13]
                   +NF_D_Q_Q4_2D_Weight16[14]*PointValues[14]
                   +NF_D_Q_Q4_2D_Weight16[15]*PointValues[15]
                   +NF_D_Q_Q4_2D_Weight16[16]*PointValues[16]
                   +NF_D_Q_Q4_2D_Weight16[17]*PointValues[17]
                   +NF_D_Q_Q4_2D_Weight16[18]*PointValues[18]
                   +NF_D_Q_Q4_2D_Weight16[19]*PointValues[19]
                   +NF_D_Q_Q4_2D_Weight16[20]*PointValues[20]
                   +NF_D_Q_Q4_2D_Weight16[21]*PointValues[21]
                   +NF_D_Q_Q4_2D_Weight16[22]*PointValues[22]
                   +NF_D_Q_Q4_2D_Weight16[23]*PointValues[23]
                   +NF_D_Q_Q4_2D_Weight16[24]*PointValues[24];

  Functionals[17] =  NF_D_Q_Q4_2D_Weight17[0]*PointValues[0]
                   +NF_D_Q_Q4_2D_Weight17[1]*PointValues[1]
                   +NF_D_Q_Q4_2D_Weight17[2]*PointValues[2]
                   +NF_D_Q_Q4_2D_Weight17[3]*PointValues[3]
                   +NF_D_Q_Q4_2D_Weight17[4]*PointValues[4]
                   +NF_D_Q_Q4_2D_Weight17[5]*PointValues[5]
                   +NF_D_Q_Q4_2D_Weight17[6]*PointValues[6]
                   +NF_D_Q_Q4_2D_Weight17[7]*PointValues[7]
                   +NF_D_Q_Q4_2D_Weight17[8]*PointValues[8]
                   +NF_D_Q_Q4_2D_Weight17[9]*PointValues[9]
                   +NF_D_Q_Q4_2D_Weight17[10]*PointValues[10]
                   +NF_D_Q_Q4_2D_Weight17[11]*PointValues[11]
                   +NF_D_Q_Q4_2D_Weight17[12]*PointValues[12]
                   +NF_D_Q_Q4_2D_Weight17[13]*PointValues[13]
                   +NF_D_Q_Q4_2D_Weight17[14]*PointValues[14]
                   +NF_D_Q_Q4_2D_Weight17[15]*PointValues[15]
                   +NF_D_Q_Q4_2D_Weight17[16]*PointValues[16]
                   +NF_D_Q_Q4_2D_Weight17[17]*PointValues[17]
                   +NF_D_Q_Q4_2D_Weight17[18]*PointValues[18]
                   +NF_D_Q_Q4_2D_Weight17[19]*PointValues[19]
                   +NF_D_Q_Q4_2D_Weight17[20]*PointValues[20]
                   +NF_D_Q_Q4_2D_Weight17[21]*PointValues[21]
                   +NF_D_Q_Q4_2D_Weight17[22]*PointValues[22]
                   +NF_D_Q_Q4_2D_Weight17[23]*PointValues[23]
                   +NF_D_Q_Q4_2D_Weight17[24]*PointValues[24];

  Functionals[18] =  NF_D_Q_Q4_2D_Weight18[0]*PointValues[0]
                   +NF_D_Q_Q4_2D_Weight18[1]*PointValues[1]
                   +NF_D_Q_Q4_2D_Weight18[2]*PointValues[2]
                   +NF_D_Q_Q4_2D_Weight18[3]*PointValues[3]
                   +NF_D_Q_Q4_2D_Weight18[4]*PointValues[4]
                   +NF_D_Q_Q4_2D_Weight18[5]*PointValues[5]
                   +NF_D_Q_Q4_2D_Weight18[6]*PointValues[6]
                   +NF_D_Q_Q4_2D_Weight18[7]*PointValues[7]
                   +NF_D_Q_Q4_2D_Weight18[8]*PointValues[8]
                   +NF_D_Q_Q4_2D_Weight18[9]*PointValues[9]
                   +NF_D_Q_Q4_2D_Weight18[10]*PointValues[10]
                   +NF_D_Q_Q4_2D_Weight18[11]*PointValues[11]
                   +NF_D_Q_Q4_2D_Weight18[12]*PointValues[12]
                   +NF_D_Q_Q4_2D_Weight18[13]*PointValues[13]
                   +NF_D_Q_Q4_2D_Weight18[14]*PointValues[14]
                   +NF_D_Q_Q4_2D_Weight18[15]*PointValues[15]
                   +NF_D_Q_Q4_2D_Weight18[16]*PointValues[16]
                   +NF_D_Q_Q4_2D_Weight18[17]*PointValues[17]
                   +NF_D_Q_Q4_2D_Weight18[18]*PointValues[18]
                   +NF_D_Q_Q4_2D_Weight18[19]*PointValues[19]
                   +NF_D_Q_Q4_2D_Weight18[20]*PointValues[20]
                   +NF_D_Q_Q4_2D_Weight18[21]*PointValues[21]
                   +NF_D_Q_Q4_2D_Weight18[22]*PointValues[22]
                   +NF_D_Q_Q4_2D_Weight18[23]*PointValues[23]
                   +NF_D_Q_Q4_2D_Weight18[24]*PointValues[24];

  Functionals[19] =  NF_D_Q_Q4_2D_Weight19[0]*PointValues[0]
                   +NF_D_Q_Q4_2D_Weight19[1]*PointValues[1]
                   +NF_D_Q_Q4_2D_Weight19[2]*PointValues[2]
                   +NF_D_Q_Q4_2D_Weight19[3]*PointValues[3]
                   +NF_D_Q_Q4_2D_Weight19[4]*PointValues[4]
                   +NF_D_Q_Q4_2D_Weight19[5]*PointValues[5]
                   +NF_D_Q_Q4_2D_Weight19[6]*PointValues[6]
                   +NF_D_Q_Q4_2D_Weight19[7]*PointValues[7]
                   +NF_D_Q_Q4_2D_Weight19[8]*PointValues[8]
                   +NF_D_Q_Q4_2D_Weight19[9]*PointValues[9]
                   +NF_D_Q_Q4_2D_Weight19[10]*PointValues[10]
                   +NF_D_Q_Q4_2D_Weight19[11]*PointValues[11]
                   +NF_D_Q_Q4_2D_Weight19[12]*PointValues[12]
                   +NF_D_Q_Q4_2D_Weight19[13]*PointValues[13]
                   +NF_D_Q_Q4_2D_Weight19[14]*PointValues[14]
                   +NF_D_Q_Q4_2D_Weight19[15]*PointValues[15]
                   +NF_D_Q_Q4_2D_Weight19[16]*PointValues[16]
                   +NF_D_Q_Q4_2D_Weight19[17]*PointValues[17]
                   +NF_D_Q_Q4_2D_Weight19[18]*PointValues[18]
                   +NF_D_Q_Q4_2D_Weight19[19]*PointValues[19]
                   +NF_D_Q_Q4_2D_Weight19[20]*PointValues[20]
                   +NF_D_Q_Q4_2D_Weight19[21]*PointValues[21]
                   +NF_D_Q_Q4_2D_Weight19[22]*PointValues[22]
                   +NF_D_Q_Q4_2D_Weight19[23]*PointValues[23]
                   +NF_D_Q_Q4_2D_Weight19[24]*PointValues[24];

  Functionals[20] =  NF_D_Q_Q4_2D_Weight20[0]*PointValues[0]
                   +NF_D_Q_Q4_2D_Weight20[1]*PointValues[1]
                   +NF_D_Q_Q4_2D_Weight20[2]*PointValues[2]
                   +NF_D_Q_Q4_2D_Weight20[3]*PointValues[3]
                   +NF_D_Q_Q4_2D_Weight20[4]*PointValues[4]
                   +NF_D_Q_Q4_2D_Weight20[5]*PointValues[5]
                   +NF_D_Q_Q4_2D_Weight20[6]*PointValues[6]
                   +NF_D_Q_Q4_2D_Weight20[7]*PointValues[7]
                   +NF_D_Q_Q4_2D_Weight20[8]*PointValues[8]
                   +NF_D_Q_Q4_2D_Weight20[9]*PointValues[9]
                   +NF_D_Q_Q4_2D_Weight20[10]*PointValues[10]
                   +NF_D_Q_Q4_2D_Weight20[11]*PointValues[11]
                   +NF_D_Q_Q4_2D_Weight20[12]*PointValues[12]
                   +NF_D_Q_Q4_2D_Weight20[13]*PointValues[13]
                   +NF_D_Q_Q4_2D_Weight20[14]*PointValues[14]
                   +NF_D_Q_Q4_2D_Weight20[15]*PointValues[15]
                   +NF_D_Q_Q4_2D_Weight20[16]*PointValues[16]
                   +NF_D_Q_Q4_2D_Weight20[17]*PointValues[17]
                   +NF_D_Q_Q4_2D_Weight20[18]*PointValues[18]
                   +NF_D_Q_Q4_2D_Weight20[19]*PointValues[19]
                   +NF_D_Q_Q4_2D_Weight20[20]*PointValues[20]
                   +NF_D_Q_Q4_2D_Weight20[21]*PointValues[21]
                   +NF_D_Q_Q4_2D_Weight20[22]*PointValues[22]
                   +NF_D_Q_Q4_2D_Weight20[23]*PointValues[23]
                   +NF_D_Q_Q4_2D_Weight20[24]*PointValues[24];

  Functionals[21] =  NF_D_Q_Q4_2D_Weight21[0]*PointValues[0]
                   +NF_D_Q_Q4_2D_Weight21[1]*PointValues[1]
                   +NF_D_Q_Q4_2D_Weight21[2]*PointValues[2]
                   +NF_D_Q_Q4_2D_Weight21[3]*PointValues[3]
                   +NF_D_Q_Q4_2D_Weight21[4]*PointValues[4]
                   +NF_D_Q_Q4_2D_Weight21[5]*PointValues[5]
                   +NF_D_Q_Q4_2D_Weight21[6]*PointValues[6]
                   +NF_D_Q_Q4_2D_Weight21[7]*PointValues[7]
                   +NF_D_Q_Q4_2D_Weight21[8]*PointValues[8]
                   +NF_D_Q_Q4_2D_Weight21[9]*PointValues[9]
                   +NF_D_Q_Q4_2D_Weight21[10]*PointValues[10]
                   +NF_D_Q_Q4_2D_Weight21[11]*PointValues[11]
                   +NF_D_Q_Q4_2D_Weight21[12]*PointValues[12]
                   +NF_D_Q_Q4_2D_Weight21[13]*PointValues[13]
                   +NF_D_Q_Q4_2D_Weight21[14]*PointValues[14]
                   +NF_D_Q_Q4_2D_Weight21[15]*PointValues[15]
                   +NF_D_Q_Q4_2D_Weight21[16]*PointValues[16]
                   +NF_D_Q_Q4_2D_Weight21[17]*PointValues[17]
                   +NF_D_Q_Q4_2D_Weight21[18]*PointValues[18]
                   +NF_D_Q_Q4_2D_Weight21[19]*PointValues[19]
                   +NF_D_Q_Q4_2D_Weight21[20]*PointValues[20]
                   +NF_D_Q_Q4_2D_Weight21[21]*PointValues[21]
                   +NF_D_Q_Q4_2D_Weight21[22]*PointValues[22]
                   +NF_D_Q_Q4_2D_Weight21[23]*PointValues[23]
                   +NF_D_Q_Q4_2D_Weight21[24]*PointValues[24];

  Functionals[22] =  NF_D_Q_Q4_2D_Weight22[0]*PointValues[0]
                   +NF_D_Q_Q4_2D_Weight22[1]*PointValues[1]
                   +NF_D_Q_Q4_2D_Weight22[2]*PointValues[2]
                   +NF_D_Q_Q4_2D_Weight22[3]*PointValues[3]
                   +NF_D_Q_Q4_2D_Weight22[4]*PointValues[4]
                   +NF_D_Q_Q4_2D_Weight22[5]*PointValues[5]
                   +NF_D_Q_Q4_2D_Weight22[6]*PointValues[6]
                   +NF_D_Q_Q4_2D_Weight22[7]*PointValues[7]
                   +NF_D_Q_Q4_2D_Weight22[8]*PointValues[8]
                   +NF_D_Q_Q4_2D_Weight22[9]*PointValues[9]
                   +NF_D_Q_Q4_2D_Weight22[10]*PointValues[10]
                   +NF_D_Q_Q4_2D_Weight22[11]*PointValues[11]
                   +NF_D_Q_Q4_2D_Weight22[12]*PointValues[12]
                   +NF_D_Q_Q4_2D_Weight22[13]*PointValues[13]
                   +NF_D_Q_Q4_2D_Weight22[14]*PointValues[14]
                   +NF_D_Q_Q4_2D_Weight22[15]*PointValues[15]
                   +NF_D_Q_Q4_2D_Weight22[16]*PointValues[16]
                   +NF_D_Q_Q4_2D_Weight22[17]*PointValues[17]
                   +NF_D_Q_Q4_2D_Weight22[18]*PointValues[18]
                   +NF_D_Q_Q4_2D_Weight22[19]*PointValues[19]
                   +NF_D_Q_Q4_2D_Weight22[20]*PointValues[20]
                   +NF_D_Q_Q4_2D_Weight22[21]*PointValues[21]
                   +NF_D_Q_Q4_2D_Weight22[22]*PointValues[22]
                   +NF_D_Q_Q4_2D_Weight22[23]*PointValues[23]
                   +NF_D_Q_Q4_2D_Weight22[24]*PointValues[24];

  Functionals[23] =  NF_D_Q_Q4_2D_Weight23[0]*PointValues[0]
                   +NF_D_Q_Q4_2D_Weight23[1]*PointValues[1]
                   +NF_D_Q_Q4_2D_Weight23[2]*PointValues[2]
                   +NF_D_Q_Q4_2D_Weight23[3]*PointValues[3]
                   +NF_D_Q_Q4_2D_Weight23[4]*PointValues[4]
                   +NF_D_Q_Q4_2D_Weight23[5]*PointValues[5]
                   +NF_D_Q_Q4_2D_Weight23[6]*PointValues[6]
                   +NF_D_Q_Q4_2D_Weight23[7]*PointValues[7]
                   +NF_D_Q_Q4_2D_Weight23[8]*PointValues[8]
                   +NF_D_Q_Q4_2D_Weight23[9]*PointValues[9]
                   +NF_D_Q_Q4_2D_Weight23[10]*PointValues[10]
                   +NF_D_Q_Q4_2D_Weight23[11]*PointValues[11]
                   +NF_D_Q_Q4_2D_Weight23[12]*PointValues[12]
                   +NF_D_Q_Q4_2D_Weight23[13]*PointValues[13]
                   +NF_D_Q_Q4_2D_Weight23[14]*PointValues[14]
                   +NF_D_Q_Q4_2D_Weight23[15]*PointValues[15]
                   +NF_D_Q_Q4_2D_Weight23[16]*PointValues[16]
                   +NF_D_Q_Q4_2D_Weight23[17]*PointValues[17]
                   +NF_D_Q_Q4_2D_Weight23[18]*PointValues[18]
                   +NF_D_Q_Q4_2D_Weight23[19]*PointValues[19]
                   +NF_D_Q_Q4_2D_Weight23[20]*PointValues[20]
                   +NF_D_Q_Q4_2D_Weight23[21]*PointValues[21]
                   +NF_D_Q_Q4_2D_Weight23[22]*PointValues[22]
                   +NF_D_Q_Q4_2D_Weight23[23]*PointValues[23]
                   +NF_D_Q_Q4_2D_Weight23[24]*PointValues[24];

  Functionals[24] =  NF_D_Q_Q4_2D_Weight24[0]*PointValues[0]
                   +NF_D_Q_Q4_2D_Weight24[1]*PointValues[1]
                   +NF_D_Q_Q4_2D_Weight24[2]*PointValues[2]
                   +NF_D_Q_Q4_2D_Weight24[3]*PointValues[3]
                   +NF_D_Q_Q4_2D_Weight24[4]*PointValues[4]
                   +NF_D_Q_Q4_2D_Weight24[5]*PointValues[5]
                   +NF_D_Q_Q4_2D_Weight24[6]*PointValues[6]
                   +NF_D_Q_Q4_2D_Weight24[7]*PointValues[7]
                   +NF_D_Q_Q4_2D_Weight24[8]*PointValues[8]
                   +NF_D_Q_Q4_2D_Weight24[9]*PointValues[9]
                   +NF_D_Q_Q4_2D_Weight24[10]*PointValues[10]
                   +NF_D_Q_Q4_2D_Weight24[11]*PointValues[11]
                   +NF_D_Q_Q4_2D_Weight24[12]*PointValues[12]
                   +NF_D_Q_Q4_2D_Weight24[13]*PointValues[13]
                   +NF_D_Q_Q4_2D_Weight24[14]*PointValues[14]
                   +NF_D_Q_Q4_2D_Weight24[15]*PointValues[15]
                   +NF_D_Q_Q4_2D_Weight24[16]*PointValues[16]
                   +NF_D_Q_Q4_2D_Weight24[17]*PointValues[17]
                   +NF_D_Q_Q4_2D_Weight24[18]*PointValues[18]
                   +NF_D_Q_Q4_2D_Weight24[19]*PointValues[19]
                   +NF_D_Q_Q4_2D_Weight24[20]*PointValues[20]
                   +NF_D_Q_Q4_2D_Weight24[21]*PointValues[21]
                   +NF_D_Q_Q4_2D_Weight24[22]*PointValues[22]
                   +NF_D_Q_Q4_2D_Weight24[23]*PointValues[23]
                   +NF_D_Q_Q4_2D_Weight24[24]*PointValues[24];

}

void NF_D_Q_Q4_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
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

TNodalFunctional2D *NF_D_Q_Q4_2D_Obj = new TNodalFunctional2D
        (NF_D_Q_Q4_2D, 25, 0, 25, 0, NF_D_Q_Q4_2D_Xi, NF_D_Q_Q4_2D_Eta,
         NF_D_Q_Q4_2D_T, NF_D_Q_Q4_2D_EvalAll, NULL);
