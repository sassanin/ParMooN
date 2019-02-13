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
   
const int nGrid = 94;
const int nGridL = 100;
double l_grid[nGridL] = {2.500000e-06,1.356050e-05,3.068150e-05,4.780300e-05,6.492450e-05,8.204550e-05,9.916650e-05,1.162875e-04,1.334090e-04,1.505305e-04,1.676515e-04,1.847725e-04,2.018940e-04,2.190155e-04,2.361365e-04,2.532575e-04,2.703785e-04,2.875000e-04,3.046215e-04,3.217425e-04,3.388635e-04,3.559845e-04,3.731060e-04,3.902275e-04,4.073485e-04,4.244695e-04,4.415910e-04,4.587125e-04,4.758335e-04,4.929545e-04,5.100755e-04,5.271970e-04,5.443185e-04,5.614395e-04,5.785605e-04,5.956815e-04,6.128030e-04,6.299245e-04,6.470455e-04,6.641665e-04,6.812875e-04,6.984090e-04,7.155305e-04,7.326515e-04,7.497725e-04,7.668940e-04,7.840155e-04,8.011365e-04,8.182575e-04,8.353785e-04,8.525000e-04,8.696215e-04,8.867425e-04,9.038635e-04,9.209845e-04,9.381060e-04,9.552275e-04,9.723485e-04,9.894695e-04,1.006591e-03,1.023712e-03,1.040833e-03,1.057955e-03,1.075075e-03,1.092197e-03,1.109318e-03,1.126439e-03,1.143561e-03,1.160682e-03,1.177803e-03,1.194925e-03,1.212045e-03,1.229167e-03,1.246287e-03,1.263409e-03,1.280530e-03,1.297651e-03,1.314773e-03,1.331894e-03,1.349016e-03,1.366137e-03,1.383257e-03,1.400379e-03,1.417500e-03,1.434622e-03,1.451742e-03,1.468863e-03,1.485985e-03,1.503106e-03,1.520228e-03,1.537348e-03,1.554470e-03,1.571591e-03,1.588713e-03,1.605834e-03,1.622954e-03,1.640075e-03,1.657197e-03,1.674318e-03,1.691440e-03 };
double l_input[nGridL] = {0.000000e+00,1.539144e+09,2.082543e+09,1.056886e+09,1.032820e+09,9.174217e+08,7.339395e+08,6.383108e+08,5.630089e+08,4.776040e+08,3.952118e+08,3.244646e+08,2.668657e+08,2.186269e+08,1.774119e+08,1.430176e+08,1.142823e+08,9.072941e+07,7.094795e+07,5.432124e+07,4.189605e+07,3.280686e+07,2.562972e+07,2.013453e+07,1.620424e+07,1.315460e+07,1.040263e+07,7.846039e+06,5.728553e+06,4.495501e+06,3.800058e+06,3.005911e+06,2.229857e+06,1.675695e+06,1.307405e+06,1.003372e+06,7.486895e+05,6.805768e+05,7.435860e+05,7.812756e+05,7.396855e+05,6.082531e+05,3.731120e+05,1.582579e+05,8.976509e+04,1.072953e+05,9.867716e+04,5.302092e+04,3.384872e+04,7.221184e+04,1.380859e+05,1.713877e+05,1.450725e+05,9.890296e+04,6.954686e+04,3.619885e+04,1.297526e+04,1.931888e+04,5.546763e+04,1.066047e+05,1.070156e+05,5.457797e+04,1.540227e+04,3.078500e+03,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00 };

double m_grid[nGrid] = {0, 1.16415e-010, 2.32829e-010, 4.65659e-010, 9.31319e-010, 1.86264e-009, 3.72527e-009, 7.45055e-009, 1.49011e-008, 2.98022e-008, 5.96047e-008, 8.94068e-008, 1.19209e-007, 1.49012e-007, 1.78814e-007, 2.08616e-007, 2.38419e-007, 2.98023e-007, 3.57628e-007, 4.17232e-007, 4.76837e-007, 5.96046e-007, 7.15255e-007, 8.34465e-007, 9.53674e-007, 1.19209e-006, 1.43051e-006, 1.66893e-006, 1.90735e-006, 2.38419e-006, 2.86102e-006, 3.33786e-006, 3.8147e-006, 4.76837e-006, 5.72205e-006, 6.67572e-006, 7.6294e-006, 8.58307e-006, 9.53675e-006, 1.04904e-005, 1.14441e-005, 1.23978e-005, 1.33514e-005, 1.43051e-005, 1.52588e-005, 1.90735e-005, 2.28882e-005, 2.67029e-005, 3.05176e-005, 3.8147e-005, 4.57764e-005, 6.10352e-005, 9.15528e-005, 0.00012207, 0.000152588, 0.000183105, 0.000213623, 0.000244141, 0.000305176, 0.000366211, 0.000427246, 0.000488281, 0.000549317, 0.000610351, 0.000671387, 0.000732422, 0.000854492, 0.000976562, 0.0012207, 0.00146484, 0.00170898, 0.00195313, 0.00244141, 0.00292969, 0.00341797, 0.00390625, 0.00488281, 0.00585937, 0.00683594, 0.0078125, 0.00976563, 0.0117187, 0.0136719, 0.015625, 0.0195312, 0.0234375, 0.0273437, 0.03125, 0.046875, 0.0625, 0.125, 0.25, 0.5, 1};
double m_input[nGrid];
//const int nGrid = 9;
//double m_grid[nGrid] = {0 , 0,0625 , 0,125 , 0,1875 , 0,25 , 0,375 , 0,5 , 0,75 ,1,000000e+00 };
