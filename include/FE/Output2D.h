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
   
// =======================================================================
// @(#)Output2D.h        1.3 11/15/99
// 
// Class:       TOutput2D
// Purpose:     store given data and realize output
//
// Author:      Gunar Matthies (21.08.1998)
//
// History:     start of implementation 21.08.1998 (Gunar Matthies)
//        :     parallel vtk output 26.09.09 (Sashikumaar Ganesan)
// =======================================================================

#ifndef __OUTPUT2D__
#define __OUTPUT2D__

#include <FEVectFunct2D.h>
#include <Domain.h>

#ifdef _MPI
#  include "mpi.h"
#endif

/** store given data and realize output */
class TOutput2D
{
  protected:
    /** collection for all spaces and functions */
    TCollection *Coll;

    /** number of stored FESpace */
    int N_FESpaces;

    /** maximal storage for FESpaces */
    int MaxN_FESpaces;

    /** array of stored FESpaces */
    TFESpace2D **FESpaceArray;

    /** number of stored scalar variables = TFEFunction */
    int N_ScalarVar;

    /** maximal storage for FEFunction */
    int MaxN_ScalarVar;

    /** array of stored scalar variables */
    TFEFunction2D **FEFunctionArray;

    /** number of stored vector-valued variables = TFEFunction */
    int N_VectorVar;

    /** maximal storage for FEVectFunct */
    int MaxN_VectorVar;

    /** array of stored vector-valued variables */
    TFEVectFunct2D **FEVectFunctArray;

    /** number of stored paramters */
    int N_Parameters;

    /** maximal storage for parameters */
    int MaxN_Parameters;

    /** values of parameters */
    double *ParameterValues;

    /** description for parameters */
    const char **ParameterDescription;

    /** corresponding domain */
    TDomain *Domain;

    /** add a FESpace into this output object (internal use) */
    int AddFESpace(TFESpace2D *fespace);

    /** internal data storage **/
    struct TOutputData
    {
      int N_Nodes;
      int N_Data;
      TVertex **Nodes;
      int *ConList;

      double **FEFuncValues; // first two for coords

      enum CELLTYPE {TRIA=2, QUAD} Type;

      TOutputData() : N_Nodes(0), N_Data(0), Nodes(0),
		      ConList(0), FEFuncValues(0) {};
      ~TOutputData();
    };

    TOutputData *Data;

  protected:
    void ComputeOutputData();
    void ComputeFEValues();

  public:
    /** constructor: maximum number of these things */
    TOutput2D(int maxn_fespaces, int maxn_scalar, int maxn_vect, 
              int maxn_parameters, TDomain *domain);

    /** destructor: freeing all allocated space */
    ~TOutput2D();

    /** add a FEFunction into this output object */
    int AddFEFunction(TFEFunction2D *fefunction);

    /** add a FEVectFunct into this output object */
    int AddFEVectFunct(TFEVectFunct2D *fevectfunct);

    /** add parameter into this output object */
    int AddParameter(double value, const char *descr);

    /** write stored data. This calls the other Write* functions. */
    int Write(std::string basename, int i=1, double t=0.);
    
    /** write stored data into a grape file */
    int WriteGrape(const char *name);

    /** write stored data into a gunplot file */
    int WriteGnuplot(const char *name);

    /** write stored data into a VTK file */
    int WriteVtk(const char *name);
    
    /** write a discontinuous function into a VTK file */
    void WriteVtkDiscontinuous(const char *fileName, 
                               int N_LocVertices, TVertex **Vertices);

    /** write stored data into a MATLAB file */
    int WriteMatlab(const char *name);

    /** write stored data into a MATLAB file */
    int WriteMatlabOld(const char *name);

    /** write scalar function into a gunplot file */
    int WriteGNU_iso(const char *name, int scalar);

    /** output for GMV file */
    int WriteGMV(const char *name);    

    /** write stored PARALLEL data into a pvtu and vtu files (XML or paraview) */
    int Write_ParVTK(
#ifdef _MPI
                          MPI_Comm comm,
#endif
                          int img, char *subID);
    int WriteAsciiPlt(const char *filename);
    int WriteBinaryPlt(const char *filename);
};

#endif
