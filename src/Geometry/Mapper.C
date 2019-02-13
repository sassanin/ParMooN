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
// @(#)Mapper.C        1.3 11/15/99
// 
// Class:       TMapper
// Purpose:     mapper for geometric objects
//
// Author:      Volker Behns  30.07.97
//
// =======================================================================

#include <Mapper.h>

#ifndef NULL
#define NULL 0
#endif

// MapTriReg0: map tri, reg, (0,0)
static const int DatMapTriReg0RefVerts[] = {0,2,1,5,4,3};
static const int DatMapTriReg0RefEdges[] = {5,4,3,2,1,0,7,6,8};
static const int DatMapTriReg0RefFaces[] = {0,2,1,3};
static const int DatMapTriReg0OrigVerts[] = {0,2,1};
static const int DatMapTriReg0OrigEdges[] = {2,1,0};

// MapTriReg1: map tri, reg, (0,1)
static const int DatMapTriReg1RefVerts[] = {1,0,2,3,5,4};
static const int DatMapTriReg1RefEdges[] = {1,0,5,4,3,2,8,7,6};
static const int DatMapTriReg1RefFaces[] = {1,0,2,3};
static const int DatMapTriReg1OrigVerts[] = {1,0,2};
static const int DatMapTriReg1OrigEdges[] = {0,2,1};

// MapTriReg2: map tri, reg, (0,2)
static const int DatMapTriReg2RefVerts[] = {2,1,0,4,3,5};
static const int DatMapTriReg2RefEdges[] = {3,2,1,0,5,4,6,8,7};
static const int DatMapTriReg2RefFaces[] = {2,1,0,3};
static const int DatMapTriReg2OrigVerts[] = {2,1,0};
static const int DatMapTriReg2OrigEdges[] = {1,0,2};

// MapQuadReg0: map quad, reg, (0,0)
static const int DatMapQuadReg0RefVerts[] = {0,3,2,1,7,6,5,4,8};
static const int DatMapQuadReg0RefEdges[] = {7,6,5,4,3,2,1,0,11,10,9,8};
static const int DatMapQuadReg0RefFaces[] = {0,3,2,1};
static const int DatMapQuadReg0OrigVerts[] = {0,3,2,1};
static const int DatMapQuadReg0OrigEdges[] = {3,2,1,0};

// MapQuadReg1: map quad, reg, (0,1)
static const int DatMapQuadReg1RefVerts[] = {1,0,3,2,4,7,6,5,8};
static const int DatMapQuadReg1RefEdges[] = {1,0,7,6,5,4,3,2,8,11,10,9};
static const int DatMapQuadReg1RefFaces[] = {1,0,3,2};
static const int DatMapQuadReg1OrigVerts[] = {1,0,3,2};
static const int DatMapQuadReg1OrigEdges[] = {0,3,2,1};

// MapQuadReg2: map quad, reg, (0,2)
static const int DatMapQuadReg2RefVerts[] = {2,1,0,3,5,4,7,6,8};
static const int DatMapQuadReg2RefEdges[] = {3,2,1,0,7,6,5,4,9,8,11,10};
static const int DatMapQuadReg2RefFaces[] = {2,1,0,3};
static const int DatMapQuadReg2OrigVerts[] = {2,1,0,3};
static const int DatMapQuadReg2OrigEdges[] = {1,0,3,2};

// MapQuadReg3: map quad, reg, (0,3)
static const int DatMapQuadReg3RefVerts[] = {3,2,1,0,6,5,4,7,8};
static const int DatMapQuadReg3RefEdges[] = {5,4,3,2,1,0,7,6,10,9,8,11};
static const int DatMapQuadReg3RefFaces[] = {3,2,1,0};
static const int DatMapQuadReg3OrigVerts[] = {3,2,1,0};
static const int DatMapQuadReg3OrigEdges[] = {2,1,0,3};

// MapTriBis0 with MapType 0 (=> Neighb is TriBis2)
static const int DatMapTriBis00RefVerts[] = {0,2,1,3};
static const int DatMapTriBis00RefEdges[] = {2,1,0,3,4};
static const int DatMapTriBis00RefFaces[] = {1,0};
static const int DatMapTriBis00OrigVerts[] = {0,2,1};
static const int DatMapTriBis00OrigEdges[] = {2,1,0};

// MapTriBis0 with MapType 1 (=> Neighb is TriBis0)
static const int DatMapTriBis01RefVerts[] = {1,0,2,3};
static const int DatMapTriBis01RefEdges[] = {1,0,3,2,4};
static const int DatMapTriBis01RefFaces[] = {1,0};
static const int DatMapTriBis01OrigVerts[] = {1,0,2};
static const int DatMapTriBis01OrigEdges[] = {0,2,1};

// MapTriBis0 with MapType 2 (=> Neighb is TriBis1)
static const int DatMapTriBis02RefVerts[] = {2,1,0,3};
static const int DatMapTriBis02RefEdges[] = {2,1,0,3,4};
static const int DatMapTriBis02RefFaces[] = {1,0};
static const int DatMapTriBis02OrigVerts[] = {2,1,0};
static const int DatMapTriBis02OrigEdges[] = {1,0,2};

// MapTriBis1 with MapType 0 (=> Neighb is TriBis1)
static const int DatMapTriBis10RefVerts[] = {0,2,1,3};
static const int DatMapTriBis10RefEdges[] = {3,2,1,0,4};
static const int DatMapTriBis10RefFaces[] = {1,0};
static const int DatMapTriBis10OrigVerts[] = {0,2,1};
static const int DatMapTriBis10OrigEdges[] = {2,1,0};

// MapTriBis1 with MapType 1 (=> Neighb is TriBis2)
static const int DatMapTriBis11RefVerts[] = {1,0,2,3};
static const int DatMapTriBis11RefEdges[] = {3,2,1,0,4};
static const int DatMapTriBis11RefFaces[] = {1,0};
static const int DatMapTriBis11OrigVerts[] = {1,0,2};
static const int DatMapTriBis11OrigEdges[] = {0,2,1};

// MapTriBis1 with MapType 2 (=> Neighb is TriBis0)
static const int DatMapTriBis12RefVerts[] = {2,1,0,3};
static const int DatMapTriBis12RefEdges[] = {2,1,0,3,4};
static const int DatMapTriBis12RefFaces[] = {1,0};
static const int DatMapTriBis12OrigVerts[] = {2,1,0};
static const int DatMapTriBis12OrigEdges[] = {1,0,2};

// MapTriBis2 with MapType 0 (=> Neighb is TriBis0)
static const int DatMapTriBis20RefVerts[] = {0,2,1,3};
static const int DatMapTriBis20RefEdges[] = {2,1,0,3,4};
static const int DatMapTriBis20RefFaces[] = {1,0};
static const int DatMapTriBis20OrigVerts[] = {0,2,1};
static const int DatMapTriBis20OrigEdges[] = {2,1,0};

// MapTriBis2 with MapType 1 (=> Neighb is TriBis1)
static const int DatMapTriBis21RefVerts[] = {1,0,2,3};
static const int DatMapTriBis21RefEdges[] = {3,2,1,0,4};
static const int DatMapTriBis21RefFaces[] = {1,0};
static const int DatMapTriBis21OrigVerts[] = {1,0,2};
static const int DatMapTriBis21OrigEdges[] = {0,2,1};

// MapTriBis2 with MapType 2 (=> Neighb is TriBis2)
static const int DatMapTriBis22RefVerts[] = {2,1,0,3};
static const int DatMapTriBis22RefEdges[] = {3,2,1,0,4};
static const int DatMapTriBis22RefFaces[] = {1,0};
static const int DatMapTriBis22OrigVerts[] = {2,1,0};
static const int DatMapTriBis22OrigEdges[] = {1,0,2};

/*
 * Double Bisections
 */
// MapTriBis01 with MapType 0 (=> Neighb is TriBis21)
static const int DatMapTriBis010RefVerts[] = {0,2,1,3,4};
static const int DatMapTriBis010RefEdges[] = {4,3,2,1,0,5,6};
static const int DatMapTriBis010RefFaces[] = {0,1,2};
static const int DatMapTriBis010OrigVerts[] = {0,2,1};
static const int DatMapTriBis010OrigEdges[] = {2,1,0};

// MapTriBis01 with MapType 1 (=> Neighb is TriBis02)
static const int DatMapTriBis011RefVerts[] = {1,0,2,3,4};
static const int DatMapTriBis011RefEdges[] = {1,0,4,3,2,5,6};
static const int DatMapTriBis011RefFaces[] = {2,0,1};
static const int DatMapTriBis011OrigVerts[] = {1,0,2};
static const int DatMapTriBis011OrigEdges[] = {0,2,1};

// MapTriBis01 with MapType 2 (=> Neighb is TriBis10)
static const int DatMapTriBis012RefVerts[] = {2,1,0,3,4};
static const int DatMapTriBis012RefEdges[] = {3,2,1,0,4,5,6};
static const int DatMapTriBis012RefFaces[] = {2,1,0};
static const int DatMapTriBis012OrigVerts[] = {2,1,0};
static const int DatMapTriBis012OrigEdges[] = {1,0,2};


// MapTriBis02 with MapType 0 (=> Neighb is TriBis20)
static const int DatMapTriBis020RefVerts[] = {0,2,1,3,4};
static const int DatMapTriBis020RefEdges[] = {4,3,2,1,0,5,6};
static const int DatMapTriBis020RefFaces[] = {0,1,2};
static const int DatMapTriBis020OrigVerts[] = {0,2,1};
static const int DatMapTriBis020OrigEdges[] = {2,1,0};

// MapTriBis02 with MapType 1 (=> Neighb is TriBis01)
static const int DatMapTriBis021RefVerts[] = {1,0,2,3,4};
static const int DatMapTriBis021RefEdges[] = {1,0,4,3,2,5,6};
static const int DatMapTriBis021RefFaces[] = {1,2,0};
static const int DatMapTriBis021OrigVerts[] = {1,0,2};
static const int DatMapTriBis021OrigEdges[] = {0,2,1};

// MapTriBis02 with MapType 2 (=> Neighb is TriBis12)
static const int DatMapTriBis022RefVerts[] = {2,1,0,3,4};
static const int DatMapTriBis022RefEdges[] = {2,1,0,4,3,5,6};
static const int DatMapTriBis022RefFaces[] = {1,2,0};
static const int DatMapTriBis022OrigVerts[] = {2,1,0};
static const int DatMapTriBis022OrigEdges[] = {1,0,2};


// MapTriBis10 with MapType 0 (=> Neighb is TriBis12)
static const int DatMapTriBis100RefVerts[] = {0,2,1,3,4};
static const int DatMapTriBis100RefEdges[] = {4,3,2,1,0,5,6};
static const int DatMapTriBis100RefFaces[] = {2,1,0};
static const int DatMapTriBis100OrigVerts[] = {0,2,1};
static const int DatMapTriBis100OrigEdges[] = {2,1,0};

// MapTriBis10 with MapType 1 (=> Neighb is TriBis20)
static const int DatMapTriBis101RefVerts[] = {1,0,2,3,4};
static const int DatMapTriBis101RefEdges[] = {1,0,4,3,2,5,6};
static const int DatMapTriBis101RefFaces[] = {1,0,2};
static const int DatMapTriBis101OrigVerts[] = {1,0,2};
static const int DatMapTriBis101OrigEdges[] = {0,2,1};

// MapTriBis10 with MapType 2 (=> Neighb is TriBis01)
static const int DatMapTriBis102RefVerts[] = {2,1,0,3,4};
static const int DatMapTriBis102RefEdges[] = {3,2,1,0,4,5,6};
static const int DatMapTriBis102RefFaces[] = {2,1,0};
static const int DatMapTriBis102OrigVerts[] = {2,1,0};
static const int DatMapTriBis102OrigEdges[] = {1,0,2};


// MapTriBis12 with MapType 0 (=> Neighb is TriBis10)
static const int DatMapTriBis120RefVerts[] = {0,2,1,3,4};
static const int DatMapTriBis120RefEdges[] = {4,3,2,1,0,5,6};
static const int DatMapTriBis120RefFaces[] = {2,1,0};
static const int DatMapTriBis120OrigVerts[] = {0,2,1};
static const int DatMapTriBis120OrigEdges[] = {2,1,0};

// MapTriBis12 with MapType 1 (=> Neighb is TriBis21)
static const int DatMapTriBis121RefVerts[] = {1,0,2,3,4};
static const int DatMapTriBis121RefEdges[] = {0,4,3,2,1,5,6};
static const int DatMapTriBis121RefFaces[] = {0,1,2};
static const int DatMapTriBis121OrigVerts[] = {1,0,2};
static const int DatMapTriBis121OrigEdges[] = {0,2,1};

// MapTriBis12 with MapType 2 (=> Neighb is TriBis02)
static const int DatMapTriBis122RefVerts[] = {2,1,0,3,4};
static const int DatMapTriBis122RefEdges[] = {2,1,0,4,3,5,6};
static const int DatMapTriBis122RefFaces[] = {2,0,1};
static const int DatMapTriBis122OrigVerts[] = {2,1,0};
static const int DatMapTriBis122OrigEdges[] = {1,0,2};


// MapTriBis20 with MapType 0 (=> Neighb is TriBis02)
static const int DatMapTriBis200RefVerts[] = {0,2,1,3,4};
static const int DatMapTriBis200RefEdges[] = {4,3,2,1,0,5,6};
static const int DatMapTriBis200RefFaces[] = {0,1,2};
static const int DatMapTriBis200OrigVerts[] = {0,2,1};
static const int DatMapTriBis200OrigEdges[] = {2,1,0};

// MapTriBis20 with MapType 1 (=> Neighb is TriBis10)
static const int DatMapTriBis201RefVerts[] = {1,0,2,3,4};
static const int DatMapTriBis201RefEdges[] = {1,0,4,3,2,5,6};
static const int DatMapTriBis201RefFaces[] = {1,0,2};
static const int DatMapTriBis201OrigVerts[] = {1,0,2};
static const int DatMapTriBis201OrigEdges[] = {0,2,1};

// MapTriBis20 with MapType 2 (=> Neighb is TriBis21)
static const int DatMapTriBis202RefVerts[] = {2,1,0,3,4};
static const int DatMapTriBis202RefEdges[] = {2,1,0,4,3,5,6};
static const int DatMapTriBis202RefFaces[] = {1,2,0};
static const int DatMapTriBis202OrigVerts[] = {2,1,0};
static const int DatMapTriBis202OrigEdges[] = {1,0,2};


// MapTriBis21 with MapType 0 (=> Neighb is TriBis01)
static const int DatMapTriBis210RefVerts[] = {0,2,1,3,4};
static const int DatMapTriBis210RefEdges[] = {4,3,2,1,0,5,6};
static const int DatMapTriBis210RefFaces[] = {0,1,2};
static const int DatMapTriBis210OrigVerts[] = {0,2,1};
static const int DatMapTriBis210OrigEdges[] = {2,1,0};

// MapTriBis21 with MapType 1 (=> Neighb is TriBis12)
static const int DatMapTriBis211RefVerts[] = {1,0,2,3,4};
static const int DatMapTriBis211RefEdges[] = {0,4,3,2,1,5,6};
static const int DatMapTriBis211RefFaces[] = {0,1,2};
static const int DatMapTriBis211OrigVerts[] = {1,0,2};
static const int DatMapTriBis211OrigEdges[] = {0,2,1};

// MapTriBis21 with MapType 2 (=> Neighb is TriBis20)
static const int DatMapTriBis212RefVerts[] = {2,1,0,3,4};
static const int DatMapTriBis212RefEdges[] = {2,1,0,4,3,5,6};
static const int DatMapTriBis212RefFaces[] = {2,0,1};
static const int DatMapTriBis212OrigVerts[] = {2,1,0};
static const int DatMapTriBis212OrigEdges[] = {1,0,2};


//Constructor
TMapper::TMapper(Mapper which)
{
  switch (which)
  {
  //
  // maps for triangles
  //
    case MapTriReg0: // map tri, reg, (0,0)
         MapRefVerts = (const int *) DatMapTriReg0RefVerts;
         MapRefEdges = (const int *) DatMapTriReg0RefEdges;
         MapRefFaces = (const int *) DatMapTriReg0RefFaces;

         MapOrigVerts = (const int *) DatMapTriReg0OrigVerts;
         MapOrigEdges = (const int *) DatMapTriReg0OrigEdges;

         break;

    case MapTriReg1: // map tri, reg, (0,1)
         MapRefVerts = (const int *) DatMapTriReg1RefVerts;
         MapRefEdges = (const int *) DatMapTriReg1RefEdges;
         MapRefFaces = (const int *) DatMapTriReg1RefFaces;

         MapOrigVerts = (const int *) DatMapTriReg1OrigVerts;
         MapOrigEdges = (const int *) DatMapTriReg1OrigEdges;

         break;

    case MapTriReg2: // map tri, reg, (0,2)
         MapRefVerts = (const int *) DatMapTriReg2RefVerts;
         MapRefEdges = (const int *) DatMapTriReg2RefEdges;
         MapRefFaces = (const int *) DatMapTriReg2RefFaces;

         MapOrigVerts = (const int *) DatMapTriReg2OrigVerts;
         MapOrigEdges = (const int *) DatMapTriReg2OrigEdges;

         break;

    case MapQuadReg0: // map quad, reg, (0,0)
         MapRefVerts = (const int *) DatMapQuadReg0RefVerts;
         MapRefEdges = (const int *) DatMapQuadReg0RefEdges;
         MapRefFaces = (const int *) DatMapQuadReg0RefFaces;

         MapOrigVerts = (const int *) DatMapQuadReg0OrigVerts;
         MapOrigEdges = (const int *) DatMapQuadReg0OrigEdges;

         break;

    case MapQuadReg1: // map quad, reg, (0,1)
         MapRefVerts = (const int *) DatMapQuadReg1RefVerts;
         MapRefEdges = (const int *) DatMapQuadReg1RefEdges;
         MapRefFaces = (const int *) DatMapQuadReg1RefFaces;

         MapOrigVerts = (const int *) DatMapQuadReg1OrigVerts;
         MapOrigEdges = (const int *) DatMapQuadReg1OrigEdges;

         break;

    case MapQuadReg2: // map quad, reg, (0,2)
         MapRefVerts = (const int *) DatMapQuadReg2RefVerts;
         MapRefEdges = (const int *) DatMapQuadReg2RefEdges;
         MapRefFaces = (const int *) DatMapQuadReg2RefFaces;

         MapOrigVerts = (const int *) DatMapQuadReg2OrigVerts;
         MapOrigEdges = (const int *) DatMapQuadReg2OrigEdges;

         break;

    case MapQuadReg3: // map quad, reg, (0,3)
         MapRefVerts = (const int *) DatMapQuadReg3RefVerts;
         MapRefEdges = (const int *) DatMapQuadReg3RefEdges;
         MapRefFaces = (const int *) DatMapQuadReg3RefFaces;

         MapOrigVerts = (const int *) DatMapQuadReg3OrigVerts;
         MapOrigEdges = (const int *) DatMapQuadReg3OrigEdges;

         break;

    case MapTriBis00:
         MapRefVerts = (const int *) DatMapTriBis00RefVerts;
         MapRefEdges = (const int *) DatMapTriBis00RefEdges;
         MapRefFaces = (const int *) DatMapTriBis00RefFaces;

         MapOrigVerts = (const int *) DatMapTriBis00OrigVerts;
         MapOrigEdges = (const int *) DatMapTriBis00OrigEdges;

         break;

    case MapTriBis01:
         MapRefVerts = (const int *) DatMapTriBis01RefVerts;
         MapRefEdges = (const int *) DatMapTriBis01RefEdges;
         MapRefFaces = (const int *) DatMapTriBis01RefFaces;

         MapOrigVerts = (const int *) DatMapTriBis01OrigVerts;
         MapOrigEdges = (const int *) DatMapTriBis01OrigEdges;

         break;

    case MapTriBis02:
         MapRefVerts = (const int *) DatMapTriBis02RefVerts;
         MapRefEdges = (const int *) DatMapTriBis02RefEdges;
         MapRefFaces = (const int *) DatMapTriBis02RefFaces;

         MapOrigVerts = (const int *) DatMapTriBis02OrigVerts;
         MapOrigEdges = (const int *) DatMapTriBis02OrigEdges;

         break;

    case MapTriBis10:
         MapRefVerts = (const int *) DatMapTriBis10RefVerts;
         MapRefEdges = (const int *) DatMapTriBis10RefEdges;
         MapRefFaces = (const int *) DatMapTriBis10RefFaces;

         MapOrigVerts = (const int *) DatMapTriBis10OrigVerts;
         MapOrigEdges = (const int *) DatMapTriBis10OrigEdges;

         break;

    case MapTriBis11:
         MapRefVerts = (const int *) DatMapTriBis11RefVerts;
         MapRefEdges = (const int *) DatMapTriBis11RefEdges;
         MapRefFaces = (const int *) DatMapTriBis11RefFaces;

         MapOrigVerts = (const int *) DatMapTriBis11OrigVerts;
         MapOrigEdges = (const int *) DatMapTriBis11OrigEdges;

         break;

    case MapTriBis12:
         MapRefVerts = (const int *) DatMapTriBis12RefVerts;
         MapRefEdges = (const int *) DatMapTriBis12RefEdges;
         MapRefFaces = (const int *) DatMapTriBis12RefFaces;

         MapOrigVerts = (const int *) DatMapTriBis12OrigVerts;
         MapOrigEdges = (const int *) DatMapTriBis12OrigEdges;

         break;

    case MapTriBis20:
         MapRefVerts = (const int *) DatMapTriBis20RefVerts;
         MapRefEdges = (const int *) DatMapTriBis20RefEdges;
         MapRefFaces = (const int *) DatMapTriBis20RefFaces;

         MapOrigVerts = (const int *) DatMapTriBis20OrigVerts;
         MapOrigEdges = (const int *) DatMapTriBis20OrigEdges;

         break;

    case MapTriBis21:
         MapRefVerts = (const int *) DatMapTriBis21RefVerts;
         MapRefEdges = (const int *) DatMapTriBis21RefEdges;
         MapRefFaces = (const int *) DatMapTriBis21RefFaces;

         MapOrigVerts = (const int *) DatMapTriBis21OrigVerts;
         MapOrigEdges = (const int *) DatMapTriBis21OrigEdges;

         break;

    case MapTriBis22:
         MapRefVerts = (const int *) DatMapTriBis22RefVerts;
         MapRefEdges = (const int *) DatMapTriBis22RefEdges;
         MapRefFaces = (const int *) DatMapTriBis22RefFaces;

         MapOrigVerts = (const int *) DatMapTriBis22OrigVerts;
         MapOrigEdges = (const int *) DatMapTriBis22OrigEdges;

         break;

    case MapTriBis010:
         MapRefVerts = (const int *) DatMapTriBis010RefVerts;
         MapRefEdges = (const int *) DatMapTriBis010RefEdges;
         MapRefFaces = (const int *) DatMapTriBis010RefFaces;

         MapOrigVerts = (const int *) DatMapTriBis010OrigVerts;
         MapOrigEdges = (const int *) DatMapTriBis010OrigEdges;

         break;
     
    case MapTriBis011:
         MapRefVerts = (const int *) DatMapTriBis011RefVerts;
         MapRefEdges = (const int *) DatMapTriBis011RefEdges;
         MapRefFaces = (const int *) DatMapTriBis011RefFaces;

         MapOrigVerts = (const int *) DatMapTriBis011OrigVerts;
         MapOrigEdges = (const int *) DatMapTriBis011OrigEdges;

         break;
     
    case MapTriBis012:
         MapRefVerts = (const int *) DatMapTriBis012RefVerts;
         MapRefEdges = (const int *) DatMapTriBis012RefEdges;
         MapRefFaces = (const int *) DatMapTriBis012RefFaces;

         MapOrigVerts = (const int *) DatMapTriBis012OrigVerts;
         MapOrigEdges = (const int *) DatMapTriBis012OrigEdges;

         break;
     
    case MapTriBis020:
         MapRefVerts = (const int *) DatMapTriBis020RefVerts;
         MapRefEdges = (const int *) DatMapTriBis020RefEdges;
         MapRefFaces = (const int *) DatMapTriBis020RefFaces;

         MapOrigVerts = (const int *) DatMapTriBis020OrigVerts;
         MapOrigEdges = (const int *) DatMapTriBis020OrigEdges;

         break;
     
    case MapTriBis021:
         MapRefVerts = (const int *) DatMapTriBis021RefVerts;
         MapRefEdges = (const int *) DatMapTriBis021RefEdges;
         MapRefFaces = (const int *) DatMapTriBis021RefFaces;

         MapOrigVerts = (const int *) DatMapTriBis021OrigVerts;
         MapOrigEdges = (const int *) DatMapTriBis021OrigEdges;

         break;
     
    case MapTriBis022:
         MapRefVerts = (const int *) DatMapTriBis022RefVerts;
         MapRefEdges = (const int *) DatMapTriBis022RefEdges;
         MapRefFaces = (const int *) DatMapTriBis022RefFaces;

         MapOrigVerts = (const int *) DatMapTriBis022OrigVerts;
         MapOrigEdges = (const int *) DatMapTriBis022OrigEdges;

         break;
     
    case MapTriBis100:
         MapRefVerts = (const int *) DatMapTriBis100RefVerts;
         MapRefEdges = (const int *) DatMapTriBis100RefEdges;
         MapRefFaces = (const int *) DatMapTriBis100RefFaces;

         MapOrigVerts = (const int *) DatMapTriBis100OrigVerts;
         MapOrigEdges = (const int *) DatMapTriBis100OrigEdges;

         break;
     
    case MapTriBis101:
         MapRefVerts = (const int *) DatMapTriBis101RefVerts;
         MapRefEdges = (const int *) DatMapTriBis101RefEdges;
         MapRefFaces = (const int *) DatMapTriBis101RefFaces;

         MapOrigVerts = (const int *) DatMapTriBis101OrigVerts;
         MapOrigEdges = (const int *) DatMapTriBis101OrigEdges;

         break;
     
    case MapTriBis102:
         MapRefVerts = (const int *) DatMapTriBis102RefVerts;
         MapRefEdges = (const int *) DatMapTriBis102RefEdges;
         MapRefFaces = (const int *) DatMapTriBis102RefFaces;

         MapOrigVerts = (const int *) DatMapTriBis102OrigVerts;
         MapOrigEdges = (const int *) DatMapTriBis102OrigEdges;

         break;
     
    case MapTriBis120:
         MapRefVerts = (const int *) DatMapTriBis120RefVerts;
         MapRefEdges = (const int *) DatMapTriBis120RefEdges;
         MapRefFaces = (const int *) DatMapTriBis120RefFaces;

         MapOrigVerts = (const int *) DatMapTriBis120OrigVerts;
         MapOrigEdges = (const int *) DatMapTriBis120OrigEdges;

         break;
     
    case MapTriBis121:
         MapRefVerts = (const int *) DatMapTriBis121RefVerts;
         MapRefEdges = (const int *) DatMapTriBis121RefEdges;
         MapRefFaces = (const int *) DatMapTriBis121RefFaces;

         MapOrigVerts = (const int *) DatMapTriBis121OrigVerts;
         MapOrigEdges = (const int *) DatMapTriBis121OrigEdges;

         break;
     
    case MapTriBis122:
         MapRefVerts = (const int *) DatMapTriBis122RefVerts;
         MapRefEdges = (const int *) DatMapTriBis122RefEdges;
         MapRefFaces = (const int *) DatMapTriBis122RefFaces;

         MapOrigVerts = (const int *) DatMapTriBis122OrigVerts;
         MapOrigEdges = (const int *) DatMapTriBis122OrigEdges;

         break;
     
    case MapTriBis200:
         MapRefVerts = (const int *) DatMapTriBis200RefVerts;
         MapRefEdges = (const int *) DatMapTriBis200RefEdges;
         MapRefFaces = (const int *) DatMapTriBis200RefFaces;

         MapOrigVerts = (const int *) DatMapTriBis200OrigVerts;
         MapOrigEdges = (const int *) DatMapTriBis200OrigEdges;

         break;
     
    case MapTriBis201:
         MapRefVerts = (const int *) DatMapTriBis201RefVerts;
         MapRefEdges = (const int *) DatMapTriBis201RefEdges;
         MapRefFaces = (const int *) DatMapTriBis201RefFaces;

         MapOrigVerts = (const int *) DatMapTriBis201OrigVerts;
         MapOrigEdges = (const int *) DatMapTriBis201OrigEdges;

         break;
     
    case MapTriBis202:
         MapRefVerts = (const int *) DatMapTriBis202RefVerts;
         MapRefEdges = (const int *) DatMapTriBis202RefEdges;
         MapRefFaces = (const int *) DatMapTriBis202RefFaces;

         MapOrigVerts = (const int *) DatMapTriBis202OrigVerts;
         MapOrigEdges = (const int *) DatMapTriBis202OrigEdges;

         break;
     
    case MapTriBis210:
         MapRefVerts = (const int *) DatMapTriBis210RefVerts;
         MapRefEdges = (const int *) DatMapTriBis210RefEdges;
         MapRefFaces = (const int *) DatMapTriBis210RefFaces;

         MapOrigVerts = (const int *) DatMapTriBis210OrigVerts;
         MapOrigEdges = (const int *) DatMapTriBis210OrigEdges;

         break;
     
    case MapTriBis211:
         MapRefVerts = (const int *) DatMapTriBis211RefVerts;
         MapRefEdges = (const int *) DatMapTriBis211RefEdges;
         MapRefFaces = (const int *) DatMapTriBis211RefFaces;

         MapOrigVerts = (const int *) DatMapTriBis211OrigVerts;
         MapOrigEdges = (const int *) DatMapTriBis211OrigEdges;

         break;
     
    case MapTriBis212:
         MapRefVerts = (const int *) DatMapTriBis212RefVerts;
         MapRefEdges = (const int *) DatMapTriBis212RefEdges;
         MapRefFaces = (const int *) DatMapTriBis212RefFaces;

         MapOrigVerts = (const int *) DatMapTriBis212OrigVerts;
         MapOrigEdges = (const int *) DatMapTriBis212OrigEdges;

         break;

    default:
         MapRefVerts = NULL;
         MapRefEdges = NULL;
         MapRefFaces = NULL;

         MapOrigVerts = NULL;
         MapOrigEdges = NULL;
  }
}
