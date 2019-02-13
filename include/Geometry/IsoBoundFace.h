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
// @(#)IsoBoundFace.h        1.1 08/12/99
// 
// Class:       TIsoBoundFace
// Purpose:     face on a boundary component with additional vertices
//              for isoparametric reference transformation
//
// Author:      Gunar Matthies  06.08.99
//
// =======================================================================

#ifndef __ISOBOUNDFACE__
#define __ISOBOUNDFACE__

#include <BoundFace.h>

/** Face on a boundary component */
class TIsoBoundFace : public TBoundFace
{
  protected:
    /** number of additional vertices */
    int N_Vertices;

    /** array of all additional vertices */
    TVertex **Vertices;

    /** array of all additional refvertices */
    TVertex **RefVertices;

    /** array of all additional vertices */
    double *IsoParam1;
    double *IsoParam2;


  public:
    // Constructors
    /** initialize the Face with the boundary component bdcomp and the
        paramter of starting and end point t\_0, t\_1 */
    TIsoBoundFace(TBoundComp3D *bdcomp, double *param1, double *param2);

    TIsoBoundFace(TBoundComp3D *bdcomp);

    // Methods
    /** check whether the refinement pattern on both side patch,
        dummy here: there is no neighbour */
    virtual int CheckMatchingRef(TBaseCell *Me, int J_i,
                  struct StoreGeom &Tmp);

    /** create a new instance of this class */
    virtual TJoint *NewInst(double T_0, double T_1, TBaseCell *Me);
    virtual TJoint *NewInst();

    /** return number of additional vertices */
    int GetN_Vertices()
    { return N_Vertices; }

    TVertex **GetVertices()
    { return Vertices; }

    TVertex **GetRefVertices()
    { return RefVertices; }

    void SetVertices(int n_vertices, TVertex **vertices);

    void GenVert(int N_NewVert,const int N_V, double **LinComb,
                            double *X, double *Y, double *Z );
    void GenHexaVert(int N_NewVert,const int N_V, double **LinComb,
                            double *X, double *Y, double *Z );
    void GenerateVertices(int n_vertices);

    /** return both isoparameters arrays */
    void GetParameters(double *OutIsoParam1, double *OutIsoParam2);


};

#endif
