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
   
#include <Database.h>
#include <SquareStructure2D.h>
#include <SquareMatrix2D.h>
#include <FEFunction2D.h>
#include <LinAlg.h>
//#include "../../Examples/TNSE_2D/Bulk_Fallung_Stat_copy.h"

//calculate the correspondent coordinates of a node in 3d grid  
void ComputeCoordinates(int i, int N, int Nz, double *x, double *y, double *z)
{
    double zmin = 1e-9/TDatabase::ParamDB->REACTOR_P2, zmax = 1;
    int j;
    
    
    j = i%N;
    x[0] = 1.0*j/(N-1);

    j = (int) (i/(N*N));
    z[0] = zmin + (zmax - zmin)*j/(Nz - 1);

    j = i - j*N*N;
    j = (int)(j/N);
    y[0] = 1.0*j/(N-1);
}

//assemble the matrix corresponding to the 3d finite-difference-method (7-point)
void Build_3D_FDM_Matrix(double h, TFEFunction2D *velocity1, TFEFunction2D *velocity2, 
                         TFEFunction2D *concent_C, double *f_old, 
			 double *velo1, double *velo2)
{
    int i, n, N, Nz, N_entries, N2, *row_ptr, *col_ptr, j, k, ii, BoundaryComponent;
    double x, y, z, hx, hz, values[3], u1_boundary_value, BoundaryParam, B_c_C;
    double *entries_x, *entries_y, *entries_z;
    double *velocity1_array, *velocity2_array, *concent_C_array, *G_c_C;
    double K1_x, K2_x, K3_x, K1_y, K2_y, K3_y, K1_z, K2_z, K3_z;
    double B1_x, B2_0_x, B2_1_x, B3_x, B1_y, B2_y, B3_y, B1_z, B2_0_z, B2_1_z, B3_z;
    double *result_vector_x, *result_vector_y, *result_vector_z;
    double velocity1_array_val, velocity2_array_val, concent_C_array_val, G_c_C_val;
    int very_first = 0;
    double *f_new, factor_G, *derx_val;
    BoundCond bc;
    double deltat = TDatabase::TimeDB->TIMESTEPLENGTH;
    double curr_t = TDatabase::TimeDB->CURRENTTIME;

    if (fabs (velo1[0] + 4711) < 1e-6)
    {
	very_first++;
	OutPut("very first computation of f" << endl);
    }
    //model constants (later from data file)
    double k_g = 1e-7;
    double C_g = 45.98;
    double c_C_infty = 1.37e-4;
    double l_infty = 1;
    double u_infty = TDatabase::ParamDB->P6;
    double d_p_max = TDatabase::ParamDB->REACTOR_P2;
    double k_nuc = 1e24; //k_nuc_infty = 1e+24
    double C_2 = 7.2e-9;
    double d_p_0 = 1e-6;
    double zmin = d_p_0/d_p_max;
    double zmax = 1;
    double f_infty;

    factor_G = k_g*c_C_infty*l_infty/(u_infty*d_p_max);
    f_infty = u_infty/(C_g*k_g*d_p_max*d_p_max*d_p_max*l_infty);

    //number of points in x direction
    //N = (int)sqrt(coll->GetN_Cells()+1e-4)+1;
    N = (int)(sqrt(2.0)/h+1e-4) + 1;
    //number of the layers in z direction
    Nz = TDatabase::ParamDB->N_CELL_LAYERS + 1;
    //number of unknowns in 2D (on one layer)
    N2 = N*N; //the grid is equidistant
    //OutPut("N " << N << " Nz " << Nz << endl);

    // number of unknowns in 3D
    n = N*N*Nz;
    row_ptr = new int[n+1];

    // count number of matrix N_entries
    N_entries = 0;
    for (i=0;i<n;i++)
    {
	// the node itself
	N_entries++;
	// left neighbour (x direction)
	if (i%N!=0)
	    N_entries++;
	// right neighbour (x direction)
	if ((i+1)%N!=0)
	    N_entries++;
	// lower neighbour (y direction)
	if (((int)(i/N2)==(int)((i-N)/N2))&&(i-N>=0))
	    N_entries++;
	// upper neighbour (y direction)
	if ((int)(i/N2)==(int)((i+N)/N2))
	    N_entries++;
	// lower neighbour (z direction)
	if (i>=N2)
	    N_entries++;
	// upper neighbour (z direction)
	if (i<=(Nz-1)*(N2-1))
	    N_entries++;
    }
    //OutPut("N_entries " << N_entries << endl);
 
    /*col_ptr = new int[N_entries];
   
    j = 0;
    for (i=0;i<n;i++)
    {
	row_ptr[i] = j;
	// the node itself
	col_ptr[j] = i;
	j++;
	// left neighbour (x direction)
	if (i%N!=0)
	{
	   col_ptr[j] = i-1;
	   j++;
	}
	// right neighbour (x direction)
	if ((i+1)%N!=0)
	{
	    col_ptr[j] = i+1;
	    j++;
	}
	// lower neighbour (y direction)
	if (((int)(i/N2)==(int)((i-N)/N2))&&(i-N>=0))
	{
	    col_ptr[j] = i-N;
	    j++;
	}
	// upper neighbour (y direction)
	if ((int)(i/N2)==(int)((i+N)/N2))
	{
	    col_ptr[j] = i+N;
	    j++;
	}
	// lower neighbour (z direction)
	if (i>=N2)
	{
	    col_ptr[j] = i - N2;
	    j++;
	}
	// upper neighbour (z direction)
	if (i<=(Nz-1)*(N2-1))
	{
	    col_ptr[j] = i+N2;
	    j++;
	}
    }
    */
    
    concent_C_array = new double[N2];
    memset(concent_C_array, 0, N2*SizeOfDouble);
    G_c_C = new double[n];
    memset(G_c_C, 0, n*SizeOfDouble);
    derx_val = new double[n];
    memset(derx_val, 0, n*SizeOfDouble);
    f_new = new double[n];
    memset(f_new, 0, n*SizeOfDouble);
           
    // discretization of PBE with FDM
    hx = 1.0/(N-1);
    hz = (zmax-zmin)/(Nz-1);
       
    for (i=0;i<n;i++)
    {
        ComputeCoordinates(i,N,Nz,&x,&y,&z);
	//if (curr_t > 6)
	//    OutPut(i << " " << x << " " << y << " " << z << endl);
        if (i < N2) //node is on the first layer (z=0)
	{
	    if (very_first)
	    {
		velocity1->FindGradient(x,y,values);//u1
		velo1[i] = values[0];
	    }
	    velocity1_array_val = velo1[i];
	    if (very_first)
	    {
		velocity2->FindGradient(x,y,values);//u2
		velo2[i] = values[0];
	    }
	    velocity2_array_val = velo2[i];

          concent_C->FindGradient(x,y,values);//c_C
          concent_C_array[i] = values[0];
          concent_C_array_val = values[0];

          //G_c_C[i] = (concent_C_array[i] -1)*(k_g*c_C_infty*l_infty)/(u_infty*d_p_max); //G(c_C)
	  G_c_C[i] = factor_G*(concent_C_array_val- 1);
	  G_c_C_val = G_c_C[i];  
	}
        else //node is not on the first layer
	{
          ii = i - N2*((int)(i/N2)); //the corresponding node on the first layer
          
          velocity1_array_val = velo1[ii];
          velocity2_array_val = velo2[ii];
          concent_C_array_val = concent_C_array[ii];
	  G_c_C[i] = factor_G*(concent_C_array_val-1);
	  G_c_C_val = G_c_C[i];  
	}
	
        // compute the coefficients corresponding to the 3d finite-difference-method (7-point)
	// simple upwind scheme

        if (velocity1_array_val >= 0)
	{
	    if (i%N!=0)
		derx_val[i] += velocity1_array_val*(f_old[i]-f_old[i-1])/hx;
	}
	else
	{
	    if ((i+1)%N!=0)
		derx_val[i] += velocity1_array_val*(f_old[i+1]-f_old[i])/hx;
	}

        if (velocity2_array_val >= 0)
	{
	    if (((int)(i/N2)==(int)((i-N)/N2))&&(i-N>=0))
		derx_val[i] += velocity2_array_val*(f_old[i]-f_old[i-N])/hx;
	}
	else
	{
	    if (((int)(i/N2)==(int)((i-N)/N2))&&(i-N>=0))
		derx_val[i] += velocity2_array_val*(f_old[i+N]-f_old[i])/hx;
	}

        if (factor_G >= 0)
	{
	    if (i>=N2)
	    derx_val[i] += factor_G *(f_old[i]-f_old[i-N2])/hz;
	}
	else
	{
	    if ((int)(i/N2)==(int)((i+N)/N2))
	    derx_val[i] += factor_G *(f_old[i+N2]-f_old[i])/hz;
	}	    
	f_new[i] = f_old[i] - deltat * derx_val[i];
	// set Dirichlet boundary conditions 
	// bottom face 
	if (i<N2)	
	{
	    G_c_C_val = k_g*(concent_C_array_val-c_C_infty*exp(C_2/z));
	    // compute G*n, n=(0,0,-1);
	    if (G_c_C_val > 1e-10)
	    {
		// compute rate of nucleation
		B_c_C = k_nuc*(concent_C_array_val - exp(C_2/d_p_0))*(concent_C_array_val - exp(C_2/d_p_0))
		    *(concent_C_array_val - exp(C_2/d_p_0))*(concent_C_array_val - exp(C_2/d_p_0))
		    *(concent_C_array_val - exp(C_2/d_p_0));  
		if (B_c_C < 0)
		    B_c_C = 0;
		//else
		//    OutPut(" bo " <<  B_c_C << " " << G_c_C_val);
		f_new[i] = B_c_C/ (G_c_C_val*f_infty);
            } 
	}
    }//endfor

    Dcopy(n, f_new, f_old);

    //OutPut("norm of derx_val " << sqrt(Ddot(n,derx_val,derx_val)) << endl);
    //OutPut("norm of f_new " << sqrt(Ddot(n,f_new,f_new)) << endl);
    /*OutPut("norm of result_x " << sqrt(Ddot(n,result_vector_x,result_vector_x)) << endl);
    OutPut("norm of result_y " << sqrt(Ddot(n,result_vector_y,result_vector_y)) << endl);
    OutPut("norm of result_z " << sqrt(Ddot(n,result_vector_z,result_vector_z)) << endl);
    OutPut("norm of G_c_C " << sqrt(Ddot(n,G_c_C,G_c_C)) << endl);*/

    /*double max_c_C=-1;
    for (i=0;i<N2;i++)
    {
	if (concent_C_array[i] > max_c_C)
	    max_c_C= concent_C_array[i];
    }
    OutPut("maximal_c_C(2) " <<  max_c_C << endl);
    
    double max_c_C=-1e9;
    for (i=0;i<N2;i++)
    {
	if (G_c_C[i] > max_c_C)
	    max_c_C= G_c_C[i];
    }
    OutPut("maximal_G_c_C " <<  max_c_C << endl);
    */

    delete concent_C_array;
    delete G_c_C;
    delete f_new;
    delete derx_val;
}

void Evalute_f_at_outflow(int n_dof_FDM, int Nx, int Nz, double *f)
{
    int k, i0, i;
    double xout, x, y, z;
    double zmin = 1e-9/TDatabase::ParamDB->REACTOR_P2, zmax = 1;

    // compute centre of outflow
    k = (int)TDatabase::ParamDB->P9;
    xout = (k+1)/32.0;
    //OutPut("center of outflow " << xout << endl);

    // compute index for value of f on z=0
    i0 = (int) (xout*(Nx-1));
    // check
    ComputeCoordinates(i0, Nx, Nz, &x, &y, &z);
    //OutPut(i0 << " " << x << " " << y << " " << z << endl);
    if (fabs(xout -x)>1e-6)
    {
	OutPut("Wrong index in Evalute_f_at_outflow " << endl);
	exit(1);
    }
    for (i=0;i<Nz;i++)
	OutPut( "property " << TDatabase::TimeDB->CURRENTTIME << " " <<
		zmin+i*(zmax-zmin)/(Nz-1) << " " << f[i0+i*Nx*Nx] << endl);
}

void Integral_For_Particle_Increase_Term(TFESpace2D *fespace, TFEFunction2D *fefct,
					 int n_dof_FDM, int Nx, int Nz, double *f)
{
    int i,j, N2, ii, N_Cells, N_Edges, index, k, l;
    int *GlobalNumbers, *BeginIndex, *DOF, local;
    double zmin = 1e-9/TDatabase::ParamDB->REACTOR_P2, zmax = 1, hz, z, val, hx, x, y;
    double *integral;
    TCollection *Coll;
    TBaseCell *cell;
    int *indextest;

    hz = (zmax-zmin)/(Nz-1);
    hx = 1.0/(Nx-1);
    N2 = Nx*Nx; 

    indextest = new int[N2];
    memset(indextest,0,N2*SizeOfInt);
    
    GlobalNumbers = fespace->GetGlobalNumbers();
    BeginIndex = fespace->GetBeginIndex();
    integral = fefct->GetValues();
    Coll = fespace->GetCollection(); // all spaces use same Coll
    N_Cells = Coll->GetN_Cells();
    for(i=0;i<N_Cells;i++)
    {
	cell = Coll->GetCell(i);
	N_Edges=cell->GetN_Edges();
	DOF = GlobalNumbers + BeginIndex[i];
	// loop over all vertices
	for (j=0;j<N_Edges;j++)
	{
	    // corners of mesh cell
	    x = cell->GetVertex(j)->GetX();
	    y = cell->GetVertex(j)->GetY();
	    // corresponding index of array f 
	    index = (int)(y/hx+1e-8);
	    index *=Nx;
	    index += (int)(x/hx+1e-8);

	    indextest[index]++;
	    // compute integral
	    val = zmin*zmin*f[index]/2;
	    for (k=1;k<Nz;k++)
	    {
		z = zmin + k*hz;
		ii = index+k*N2;
		val += z*z*f[ii];
	    }
	    // f should be zero at zmax, can be deleted
	    ii = index+Nz*N2;
	    val += zmax*zmax*f[ii]/2;
	    val *= hz;
	    // insert the value into the fe function
	    // compute first local index
	    switch(j)
	    {
		case 0:
		case 1:
		    local = j;
		    break;
		case 2:
		    local = 3;
		    break;
		case 3:
		    local = 2;
		    break;
	    }
	    // l - global index
	    l = DOF[local];
	    integral[l] = val;
	}
    }
    for (i=0;i<N2;i++)
    {
	if (!indextest[i])
	{
	    OutPut("Index " << i << " not found " << endl);
	    exit(1);
	}
    }
    OutPut("integral " << sqrt(Ddot(N2,integral,integral)) << endl);

    delete indextest;
}
