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
   
/****************************************************************************************
*                                                                                       *
*                          compute_median_volume_fraction.cpp                           *
*                         ------------------------------------                          *
*                                                                                       *
****************************************************************************************/


/*
   This routine generates a data file and a MATLAB file for the illustration of the median volume fraction
   d_p_50

*/

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
using namespace std;

#include "compute_dp_50.h"
// #include "pbe_analysis.h"

/*
   the arguments are the file name of the data file and the name of the new data file for MATLAB
   the best choice is the first part of the origin file name - the routine creates a new second part
   and displaces the "." by "_" in the name string  - to avoid problems with MATLAB and Windows

*/


int main( int argc, char *argv[] )
{
  // cout precision for the tests of the routine
  cout.precision(12);
  // test of the correct number of parameters
  if ( argc != 3 )
  {
    cerr << " wrong call of the function ! " << endl;
    cerr << " correct call :  file name _ name of the new data file " << endl;
    cerr << " the routine expect 3 arguments " << endl;
    return 1;
  }

  cout << " start of compute_median_volume_fraction " << endl;

  // open the data file
  ifstream data_file(argv[1]);

  if ( !data_file )
  {
    cerr << " Error : file defect " << endl;
    return 2;
  }

  // file name for the MATLAB data file with the d_p_50 values over time
  string name = argv[2];
  string end  = "_mvf.txt";

  // change all "." in the file name into a "_" => to avoid problems with windows
  int pos = 1;
  while ( pos>0 )
  {
    pos = name.find(".");
    if ( pos != -1 )
    {
      name.replace(pos,1,"_");
    }
  }

  // create the new data file
  ofstream data_file_mvf((name+end).c_str());

  char read_line[200];
  // line_1.length() == 68 for all entries in the test file -> check for other files
  string line_1, line_2, time;
  string PSD = "PSD";
  int counter = 0;
  int counter_1, N = 65;
  double A_4[N], A_5[N];
  double dp_50;

  while( data_file.getline(read_line,200))
  {
    line_1 = read_line;

    if ( (line_1[0]==PSD[0]) && (line_1[1]==PSD[1]) && (line_1[2]==PSD[2]) )
    {
      if ( counter < N )
      {
        time = line_1;
        // einlesen und kürzen der 4. und 5. Spalte
        line_2 = line_1;
        line_1.erase(0,43);
        line_1.erase(12,24);
        line_2.erase(0,56);
        // umwandeln der strings in double Zahlen
        stringstream str_1, str_2;
        str_1 << line_1;
        str_2 << line_2;
        str_1 >> A_4[counter];
        str_2 >> A_5[counter];
      }
      counter++;
      if ( counter == N )
      {
        time.erase(0,4);
        counter_1 = time.length();
        time.erase(12,counter_1);
        dp_50 = compute_dp_50(N, A_4, A_5);
        data_file_mvf << time << "  " << dp_50 << endl;
        counter = 0;
      }

    }

    line_1 = "";
    line_2 = "";
    time = "";
  }

  data_file_mvf.close();

  cout << " end of compute_median_volume_fraction " << endl;

}
