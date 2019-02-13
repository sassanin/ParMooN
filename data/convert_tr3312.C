#include<stdlib.h>
#include<malloc.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using std::cout;
using std::endl;

#define DOUBLE double

/*++++++++++++++++++++++++++++++++*/
int main(int argc, char *argv[])
{ 
  std::ifstream fp_in; 
  std::ofstream fp_out;
  char *name = "tr2360", *outfile = "sold_grid4.dat";
  char line[200], *blank = " ", *help; 
  int i, n, *flag, ind[3];
  double *x, *y;
 
  fp_in.open(name,std::ios::in);      /* try to open file */
  if (fp_in.fail())              /* file does not exists */
  {
    exit(1);
  }
  fp_out.open(outfile,std::ios::out);      /* open file */
 // read number of nodes 
  fp_in.getline(line,200,'\n');
  n = (int) strtod(line,NULL);
  cout << "n " << n << endl;
  // allocate arrays for coordinates
  x = new double[n];
  y = new double[n];
  flag = new int[n];

  fp_out <<"Grobgitter" << endl;
  fp_out <<"Parametrisierung PARXC, PARYC, TMAXC" << endl;
  fp_out <<"    3312    1721   0    3    1         NEL NVT NMT NVE NBCT" << endl;
  fp_out <<"DCORVG"<<endl;
  for (i=0;i<n;i++)
  {    
    fp_in.getline(line,200,'\n');
    x[i] = strtod(line,NULL);
    help = strstr(line,blank);
    help++;
    y[i] = strtod(help,NULL);
    help = strstr(help,blank);
    help++;
    flag[i] = (int)strtod(help,NULL);
    if (flag[i]==1)
    {
      //  cout << x[i] << " " << y[i] << "       ";
      if (fabs(y[i]) < 1e-7)
       {
           y[i] = 0;
       }
       if ((fabs(1-x[i]) < 1e-7)&&(fabs(y[i]) > 1e-7))
       {
           x[i] = y[i] + 1;
           y[i] = 0;
       }
       if ((fabs(1-y[i]) < 1e-7)&&(fabs(1-x[i]) > 1e-7))
       {
           x[i] = 3 - x[i];
           y[i] = 0;
       }
       if ((fabs(x[i]) < 1e-7)&&(fabs(y[i]) > 1e-7))
       {
           x[i] = 4 - y[i];
           y[i] = 0;
       }
       //cout << x[i] << " " << y[i] << endl;
    }  
  fp_out <<  x[i] << " " << y[i] << endl;
  }
  fp_out << "KVERT" << endl;
  while (!fp_in.eof())
  {
    fp_in.getline(line,200,'\n');
    ind[0] = (int) strtod(line,NULL);
    help = strstr(line,blank);
    help++;
    ind[1] = (int) strtod(help,NULL);
    help = strstr(help,blank);
    help++;
    ind[2] = (int)strtod(help,NULL);
    if (ind[0] > -1) 
    fp_out << ind[0]+1 << " " << ind[1]+1 << " " << ind[2]+1 << " "  << endl;
  }
  fp_out << "KNPR" << endl;  
  for (i=0;i<n;i++)
  {
    fp_out << flag[i] << " ";
    if( (i+1)%10==0)
      fp_out <<  endl;
  }
  fp_out << "KMM" << endl;  
  fp_out << "     1    6" << endl;
  fp_in.close();                 /* close files */
  fp_out.close();                 /* close files */
  return 0;
} 
