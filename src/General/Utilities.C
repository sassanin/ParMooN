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
   
/** rearrange bytes in int array */
void SwapIntArray(int *intarray, int length)
{
  char *array, c;
  int i;

  for(i=0;i<length;i++)
  {
    array = (char *)(intarray+i);
    c = array[0];
    array[0] = array[3];
    array[3] = c;

    c = array[1];
    array[1] = array[2];
    array[2] = c;
  }
} // SwapIntArray

/** rearrange bytes in float array */
void SwapFloatArray(float *floatarray, int length)
{
  char *array, c;
  int i;

  for(i=0;i<length;i++)
  {
    array = (char *)(floatarray+i);
    c = array[0];
    array[0] = array[3];
    array[3] = c;

    c = array[1];
    array[1] = array[2];
    array[2] = c;
  }
} // SwapFloatArray

/** rearrange bytes in double array */
void SwapDoubleArray(double *doublearray, int length)
{
  char *array, c;
  int i;

  for(i=0;i<length;i++)
  {
    array = (char *)(doublearray+i);
    c = array[0];
    array[0] = array[7];
    array[7] = c;

    c = array[1];
    array[1] = array[6];
    array[6] = c;

    c = array[2];
    array[2] = array[5];
    array[5] = c;

    c = array[3];
    array[3] = array[4];
    array[4] = c;
  }
} // SwapDoubleArray
