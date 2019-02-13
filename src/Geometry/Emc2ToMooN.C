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
   
// #include<stdlib.h>
// #include<malloc.h>
// #include<stdio.h>
// #include<math.h>
// #include<string.h>
// #include <fstream.h>
// #include <strstream.h>
// #include <iomanip.h>
// 
// #define DOUBLE double

// //#include"Benchmark.h"
// #include"dc.h"
// 
// /*===========================*/
// /*        open file          */
// /*===========================*/
// 
// int file_open(char *argv[], char *name,int *ftq,fstream &fp_in,fstream &fp_out)
// { 
//   int n;
//   char *ct=".am_fmt",*ct1=".ftq",*geo=".GEO",help[200];
// 
//   strcpy(name,argv[1]);       
//   *ftq=0;                         /* initialize ftq flag */
//   fp_in.open(name,ios::in);       /* try to open file */
//   if (!fp_in.fail())               /* file exists */
//     {
//       if (strstr(name,ct1))*ftq=1;/* '.ftq' is in file name */
//     }
//   else                            /* file does not exist */
//     {
//       strcat(name,ct);            /* try to open name.am_fmt */
//       fp_in.open(name,ios::in);
//       if(fp_in.fail())            /* name.am_fmt does not exist */
//         {
//           *ftq=1;
//           strcpy(name,argv[1]);   /* try to open name.ftq */
//           fp_in.open(name,ios::in);
//           if(fp_in.fail())        /* name.ftq does not exist */
//             {                     
//               cout << "file dosen't exist = " << name << endl;
//               exit(1);             /* exit program */
//             }
//         }
//     }
// 
//   cout << "input file "<< name << endl; 
//   n=strlen(name);
//   if (*ftq==1)                    /* cut suffix from input file name */
//     name[n-4]='\0'; 
//   else 
//     name[n-7]='\0';
// 
//   strcat(name,geo);               /* build output file name */ 
//   fp_out.open(name,ios::out);
//   if(fp_out.fail())
//     { 
//       cout << "cannot open output file " << name << endl;
//       exit(1);
//     }
//   cout << "output file " << name << endl; 
//  
//   return(0);
// }
// 
// /*=======================================*/
// /*   get number of nodes and elements    */
// /*=======================================*/
// 
// void get_dimension_am_fmt(fstream &fp_in,int *nr_node,int *nr_ele)
// { 
//   char help[20];        
//   
//   /* get first line and write to buffer  help(string) */  
//   fp_in >> *nr_node;
//   fp_in >> *nr_ele;
//   fp_in >> help;
//   fp_in >> help;
//   
//   cout<< "nr. of elements : " <<  *nr_ele << " nr. of nodes : "
//       << *nr_node << endl;  
//   return;
// }
// 
// void get_dimension_ftq(fstream &fp_in,int *nr_node,int *nr_ele)
// { 
//   int nr_tri,nr_qua;
// 
//   fp_in >> *nr_node; 
//   fp_in >> *nr_ele;
//   fp_in >> nr_tri;
//   fp_in >> nr_qua;
// 
//   cout << "nr. of elements : " << *nr_ele <<
//     " nr. of nodes : " << *nr_node <<
//     " nr. of triang : " << nr_tri <<
//     " nr. of quad : " << nr_qua << endl;
//   return;
// }
// 
// /*++++++++++++++++++++++++++++++++*/
// int main(int argc, char *argv[])
// { 
//   fstream fp_in,fp_out;  
//   char name[250],help[200],in_name[50];
//   int *knpr,k,nr_nodes,nr_elem,v,maxnodes,nod,boundary;
//   long int i;
//   double *knoten[2],param;
//   int *elemente[4],bnd_seg;   
//   int ftq,nodes;
// 
//   cout << "This is the converter from EMC2 to MooNMD" << endl;
//   cout << "!!!  It is assumed that all boundaries are marked in" << endl;
//   cout << "EMC2 with the numbers corrsponding to the included" << endl;
//   cout << "header file. Otherwise, the description of boundary" << endl;
//   cout << "points will be wrong in the output file  !!!" << endl;
// 
//   nod=4; 
//   if (file_open(argv,name,&ftq,fp_in,fp_out)!=0)     /* open file */ 
//     return(1);
// 
//   if (ftq==0)                                          /* am_fmt file */
//     {
//       get_dimension_am_fmt(fp_in,&nr_nodes,&nr_elem); /* get dimension of coarse grid */
//       knpr=new int[nr_nodes*sizeof(int)];               /* allocate arrays to hold data */
//       knoten[0]=new double[nr_nodes*sizeof(double)];       /* of the coarse grid */
//       knoten[1]=new double[nr_nodes*sizeof(double)];
//       elemente[0]=new int[nr_elem*sizeof(int)];
//       elemente[1]=new int[nr_elem*sizeof(int)];
//       elemente[2]=new int[nr_elem*sizeof(int)];
//   
//       maxnodes=3;
// 
//       k=nr_elem/2;      
//       for (i=1;i<=k;i++)       /* read elements */
//         {
//           fp_in >> elemente[0][2*i-2];
//           fp_in >> elemente[1][2*i-2];
//           fp_in >> elemente[2][2*i-2];
//           fp_in >> elemente[0][2*i-1];
//           fp_in >> elemente[1][2*i-1];
//           fp_in >> elemente[2][2*i-1];
//         }
//       if(nr_elem%2==1)         /* odd number of elements */ 
//         {
//           fp_in >> elemente[0][nr_elem-1];
//           fp_in >> elemente[1][nr_elem-1];
//           fp_in >> elemente[2][nr_elem-1];
//         }          
//       
//       k=nr_nodes/2;
//       for (i=1;i<=k;i++)
//         {
//           fp_in >> knoten[0][2*i-2];
//           fp_in >> knoten[1][2*i-2];
//           fp_in >> knoten[0][2*i-1];
//           fp_in >> knoten[1][2*i-1];
//         }
//       if(nr_nodes%2==1)
//         {
//           fp_in >> knoten[0][nr_nodes-1]; /* odd number of nodes */
//           fp_in >> knoten[1][nr_nodes-1];
//         }
//       for (i=0;i<nr_elem;i++)            /* read boundary information */
//         fp_in >> v;                      /* skip some lines */
//       for (i=0;i<nr_nodes;i++)            /* read boundary information */
//         fp_in >> knpr[i];
//     }                  /* end am_fmt */
//   else                 /* ftp file */
//     {
//       get_dimension_ftq(fp_in,&nr_nodes,&nr_elem);  /* read data of grid */                            
//       knpr=new int[nr_nodes*sizeof(int)];             /* allocate vectors to store data */
//       knoten[0]=new double[nr_nodes*sizeof(double)];
//       knoten[1]=new double[nr_nodes*sizeof(double)];
//       elemente[0]=new int[nr_elem*sizeof(int)];
//       elemente[1]=new int[nr_elem*sizeof(int)];
//       elemente[2]=new int[nr_elem*sizeof(int)];
//       elemente[3]=new int[nr_elem*sizeof(int)];
//       maxnodes=3;                                    /* maximal nodes of elements */
//       for (i=0;i<nr_elem;i++)                       /* read elements */
//         {
//           fp_in >> nodes;
//           if (nodes>maxnodes) maxnodes=nodes;        /* update maxnodes */
//           fp_in >> elemente[0][i];
//           fp_in >> elemente[1][i];
//           fp_in >> elemente[2][i];
//           if (nodes==3)
//             elemente[3][i]= 0;
//           else
//             fp_in >> elemente[3][i];
//           fp_in >> v;                                      /* next entry not needed */
//         }
//       for (i=0;i<nr_nodes;i++)                          /* read nodes */                 
//               {
//           fp_in >> knoten[0][i];
//           fp_in >> knoten[1][i];
//           fp_in >> knpr[i];
//               }
//     } /* end ftp */
//  
//    boundary=0;                      /* change description of boundary points */
// 
//    for (i=0;i<nr_nodes;i++)         /* due to boundary description in input file */
//      {
//        if(knpr[i]!=0)               /* boundary point */
//          {
//            ComputeBoundary(knpr[i],knoten[0][i],knoten[1][i],&bnd_seg,&param);
//            knpr[i]=bnd_seg;
//            if(knpr[i]>boundary) 
//              boundary=knpr[i];
//            knoten[0][i] = param;    /* change description */
//            knoten[1][i] = 0;        
//          }
//      }
//    
//    fp_out << "coarse grid" << endl; /* write data to output file */
//    fp_out << "Parametrisierung PARXC PARYC TMAXC"<< endl;
//    fp_out << "      " << nr_elem << "  " << nr_nodes << "   0   "
//           << nod << "     " << boundary  
//           << "       NEL NVT NMT NVE NBCT" << endl;
//    fp_out << "DCORVG" << endl;  
//    
//    for(i=0;i<nr_nodes;i++)
//      fp_out << "   "  << knoten[0][i] << "   " << knoten[1][i] << endl;
//    
//    fp_out << "KVERT" << endl;; 
//  
//    for(i=0;i<nr_elem;i++)
//      { 
//        if (ftq==0) 
//          fp_out << "   " << elemente[0][i] << "   " << elemente[1][i] << "   " 
//                 << elemente[2][i] << "   0" << endl;
//        else
//          fp_out << "   " << elemente[0][i] << "   " << elemente[1][i] << "   " 
//                 << elemente[2][i] << "   " << elemente[3 ][i] << endl;
//      }
//    fp_out <<  "KNPR" << endl;
//    for(i=0;i<nr_nodes;i++)
//     {
//       fp_out << knpr[i] << "  "; 
//       if((i+1)%15==0)
//         fp_out << endl; 
//     } 
//   fp_out << endl << "KMM" << endl; 
//   fp_out << "1 2 3 4" << endl; 
//   fp_in.close();                 /* close files */
//   fp_out.close();
//   return 0;
// } 
// 
// 
// 
// 
