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
// @(#)RefHexaRegDesc.C        1.4 11/15/99
//
// Class:       TRefHexaRegDesc
// Purpose:     refinement descriptor for regular refinement of a
//              hexahedron
//
// Author:      Volker Behns  31.07.97
//
// =======================================================================

#ifndef __3D__
#define __3D__
#endif

#include <RefHexaRegDesc.h>

static const Shapes DatChildType[] = {Hexahedron, Hexahedron, Hexahedron,
                                      Hexahedron, Hexahedron, Hexahedron,
                                      Hexahedron, Hexahedron };

static const Refinements DatFaceType[] = { QuadReg, QuadReg, QuadReg,
                                           QuadReg, QuadReg, QuadReg };

static const Refinements DatEdgeType[] = { LineReg, LineReg, LineReg, LineReg, 
                                           LineReg, LineReg, LineReg, LineReg, 
                                           LineReg, LineReg, LineReg, LineReg };

static const int DatChildVertex[][REFHEXAREGMAXN_VpC] =
         { { 0,  8, 12, 11, 15, 16, 26, 24}, { 1,  9, 12,  8, 13, 19, 26, 16},
           { 2, 10, 12,  9, 17, 22, 26, 19}, { 3, 11, 12, 10, 20, 24, 26, 22},
           { 4, 23, 25, 14, 15, 24, 26, 16}, { 5, 14, 25, 18, 13, 16, 26, 19},
           { 6, 18, 25, 21, 17, 19, 26, 22}, { 7, 21, 25, 23, 20, 22, 26, 24}};

static const int DatVertexChild[][REFHEXAREGMAXN_CpV] =
       { {0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {0,1},
         {1,2}, {2,3}, {0,3}, {0,1,2,3}, {1,5}, {4,5}, {0,4}, {0,1,4,5}, {2,6},
         {5,6}, {1,2,5,6}, {3,7}, {6,7}, {2,3,6,7}, {4,7}, {0,3,4,7}, {4,5,6,7},
         {0,1,2,3,4,5,6,7}};

static const int DatVertexChildIndex[][REFHEXAREGMAXN_CpV] =
       { {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {1,3},
         {1,3}, {1,3}, {3,1}, {2,2,2,2}, {4,4}, {3,1}, {4,4}, {5,7,7,5}, {4,4},
         {3,1}, {5,7,7,5}, {4,4}, {3,1}, {5,7,7,5}, {1,3}, {7,5,5,7}, {2,2,2,2},
         {6,6,6,6,6,6,6,6}};

static const int DatVertexChildLen[] =
                 { 1, 1, 1, 1, 1, 1, 1, 1, 2,
                   2, 2, 2, 4, 2, 2, 2, 4, 2,
                   2, 4, 2, 2, 4, 2, 4, 4, 8};

static const int DatChildEdge[][REFHEXAREGMAXN_EpC] =
                 { { 0,  8, 11,  7, 17, 18, 52, 40, 21, 48, 51, 41}, 
                   { 2,  9,  8,  1, 12, 26, 52, 18, 29, 49, 48, 19},
                   { 4, 10,  9,  3, 22, 34, 52, 26, 37, 50, 49, 27}, 
                   { 6, 11, 10,  5, 30, 40, 52, 34, 43, 51, 50, 35},
                   {38, 47, 44, 15, 16, 42, 53, 20, 41, 51, 48, 21}, 
                   {14, 44, 45, 25, 13, 20, 53, 28, 19, 48, 49, 29},
                   {24, 45, 46, 33, 23, 28, 53, 36, 27, 49, 50, 37}, 
                   {32, 46, 47, 39, 31, 36, 53, 42, 35, 50, 51, 43}};

static const int DatEdgeChild[][REFHEXAREGMAXN_CpE] =
               { {0}, {1}, {1}, {2}, {2}, {3}, {3}, {0}, {0,1}, {1,2}, 
                 {2,3}, {0,3}, {1}, {5}, {5}, {4}, {4}, {0}, {0,1}, {1,5}, 
                 {4,5}, {0,4}, {2}, {6}, {6}, {5}, {1,2}, {2,6}, {5,6}, {1,5}, 
                 {3}, {7}, {7}, {6}, {2,3}, {3,7}, {6,7}, {2,6}, {4}, {7}, 
                 {0,3}, {0,4}, {4,7}, {3,7}, {4,5}, {5,6}, {6,7}, {4,7}, 
                 {0,1,4,5}, {1,2,5,6}, {2,3,6,7}, {0,3,4,7}, 
                 {0,1,2,3}, {4,5,6,7}};

static const int DatEdgeChildIndex[][REFHEXAREGMAXN_CpE] = 
                 { {0}, {3}, {0}, {3}, {0}, {3}, {0}, {3}, {1,2}, {1,2},
                   {1,2}, {2,1}, {4}, {4}, {0}, {3}, {4}, {4}, {5,7}, {11,8},
                   {7,5}, {8,11}, {4}, {4}, {0}, {3}, {5,7}, {11,8}, {7,5}, {8,11},
                   {4}, {4}, {0}, {3}, {5,7}, {11,8}, {7,5}, {8,11}, {0}, {3},
                   {7,5}, {11,8}, {5,7}, {8,11}, {2,1}, {2,1}, {2,1}, {1,2},
                   {9,10,10,9}, {9,10,10,9}, {9,10,10,9}, {10,9,9,10},
                   {6,6,6,6}, {6,6,6,6}};

static const int DatEdgeChildLen[] =
                 { 1, 1, 1, 1, 1, 1, 1, 1, 2, 2,
                   2, 2, 1, 1, 1, 1, 1, 1, 2, 2,
                   2, 2, 1, 1, 1, 1, 2, 2, 2, 2,
                   1, 1, 1, 1, 2, 2, 2, 2, 1, 1,
                   2, 2, 2, 2, 2, 2, 2, 2,
                   4, 4, 4, 4, 4, 4};

static const int DatChildFace[][REFHEXAREGMAXN_FpC] = 
                 { {0,4,24,28,17,32},  {1,8,29,24,5,33},   {2,12,25,29,9,34}, 
                   {3,16,28,25,13,35}, {20,18,31,27,7,32}, {21,6,27,30,11,33},
                   {22,10,30,26,15,34}, {23,14,26,31,19,35}};

static const int DatFaceChild[][REFHEXAREGMAXN_CpF] = 
                 { {0}, {1}, {2}, {3}, {0}, {1}, {5}, {4}, {1}, {2}, 
                   {6}, {5}, {2}, {3}, {7}, {6}, {3}, {0}, {4}, {7}, 
                   {4}, {5}, {6}, {7}, {0,1}, {2,3}, {6,7}, {4,5}, {0,3}, {1,2}, 
                   {5,6}, {4,7}, {0,4}, {1,5}, {2,6}, {3,7}};

static const int DatFaceChildIndex[][REFHEXAREGMAXN_CpF] =
                 { {0}, {0}, {0}, {0}, {1}, {4}, {1}, {4}, {1}, 
                   {4}, {1}, {4}, {1}, {4}, {1}, {4}, {1}, {4}, 
                   {1}, {4}, {0}, {0}, {0}, {0}, {2,3}, {2,3}, {3,2}, 
                   {3,2}, {3,2}, {2,3}, {3,2}, {2,3}, {5,5}, {5,5}, {5,5}, {5,5}};

static const int DatFaceChildLen[] =
                 {  1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 2, 2, 2,
                    2, 2, 2, 2, 2, 2, 2, 2, 2};

static const int DatEdgeVertex[][2] = 
        {  {0,8},  {8,1},  {1,9},  {9,2}, {2,10}, {10,3}, {3,11},
          {11,0}, {8,12}, {9,12},{10,12},{11,12}, {1,13}, {13,5},
          {5,14}, {14,4}, {4,15}, {15,0}, {8,16},{13,16},{14,16},
         {15,16}, {2,17}, {17,6}, {6,18}, {18,5}, {9,19},{17,19},
         {18,19},{13,19}, {3,20}, {20,7}, {7,21}, {21,6},{10,22},
         {20,22},{21,22},{17,22}, {4,23}, {23,7},{11,24},{15,24},
         {23,24},{20,24},{14,25},{18,25},{21,25},{23,25},{16,26},
         {19,26},{22,26},{24,26},{12,26},{25,26}};

static const int DatVertexEdge[][REFHEXAREGMAXN_EpV] =
    { {0,7,17}, {1,2,12}, {3,4,22}, {5,6,30}, {15,16,38}, {13,14,25}, 
      {23,24,33}, {31,32,39}, {0,1,8,18}, {2,3,9,26}, {4,5,10,34}, {6,7,11,40},
      {8,9,10,11,52}, {12,13,19,29}, {14,15,20,44}, {16,17,21,41},
      {18,19,20,21,48}, {22,23,27,37}, {24,25,28,45}, {26,27,28,29,49},
      {30,31,35,43}, {32,33,36,46}, {34,35,36,37,50}, {38,39,42,47},
      {40,41,42,43,51}, {44,45,46,47,53}, {48,49,50,51,52,53}};

static const int DatVertexEdgeIndex[][REFHEXAREGMAXN_EpV] =
        { {0,1,1}, {1,0,0}, {1,0,0}, {1,0,0}, {1,0,0}, {1,0,1}, 
          {1,0,1}, {1,0,1}, {1,0,0,0}, {1,0,0,0}, {1,0,0,0}, {1,0,0,0},
          {1,1,1,1,0}, {1,0,0,0}, {1,0,0,0}, {1,0,0,0}, {1,1,1,1,0},
          {1,0,0,0}, {1,0,0,0}, {1,1,1,1,0}, {1,0,0,0}, {1,0,0,0}, {1,1,1,1,0},
          {1,0,0,0}, {1,1,1,1,0}, {1,1,1,1,0}, {1,1,1,1,1,1}};

static const int DatVertexEdgeLen[] = 
        {  3, 3, 3, 3, 3, 3,
           3, 3, 4, 4, 4, 4,
           5, 4, 4, 4, 5,
           4, 4, 5, 4, 4, 5,
           4, 5, 5, 6};

static const int DatFaceVertex[][REFHEXAREGMAXN_VpF] = 
        { {0,8,12,11},   {1,9,12,8},    {2,10,12,9},   {3,11,12,10},
          {0,8,16,15},   {1,13,16,8},   {5,14,16,13},  {4,15,16,14},
          {1,9,19,13},   {2,17,19,9},   {6,18,19,17},  {5,13,19,18},
          {2,10,22,17},  {3,20,22,10},  {7,21,22,20},  {6,17,22,21},
          {3,11,24,20},  {0,15,24,11},  {4,23,24,15},  {7,20,24,23},
          {4,14,25,23},  {5,18,25,14},  {6,21,25,18},  {7,23,25,21},
          {8,12,26,16},  {10,22,26,12}, {21,25,26,22}, {14,16,26,25},
          {11,12,26,24}, {9,19,26,12},  {18,25,26,19}, {23,24,26,25},
          {15,16,26,24}, {13,19,26,16}, {17,22,26,19}, {20,24,26,22}};

static const int DatVertexFace[][REFHEXAREGMAXN_FpV] =
        { {0,4,17}, {1,5,8}, {2,9,12}, {3,13,16}, {7,18,20}, {6,11,21}, 
          {10,15,22}, {14,19,23}, {0,1,4,5,24}, {1,2,8,9,29}, 
          {2,3,12,13,25}, {0,3,16,17,28}, {0,1,2,3,24,25,28,29},
          {5,6,8,11,33}, {6,7,20,21,27}, {4,7,17,18,32}, {4,5,6,7,24,27,32,33},
          {9,10,12,15,34}, {10,11,21,22,30}, {8,9,10,11,29,30,33,34},
          {13,14,16,19,35}, {14,15,22,23,26}, {12,13,14,15,25,26,34,35},
          {18,19,20,23,31}, {16,17,18,19,28,31,32,35}, {20,21,22,23,26,27,30,31},
          {24,25,26,27,28,29,30,31,32,33,34,35}};

static const int DatVertexFaceIndex[][REFHEXAREGMAXN_FpV] =
        { {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, 
          {0,0,0}, {0,0,0}, {1,3,1,3,0}, {1,3,1,3,0}, 
          {1,3,1,3,0}, {3,1,1,3,0}, {2,2,2,2,1,3,1,3}, 
          {1,3,3,1,0}, {1,3,1,3,0}, {3,1,1,3,0}, {2,2,2,2,3,1,1,3}, 
          {1,3,3,1,0}, {1,3,1,3,0}, {2,2,2,2,1,3,1,3},
          {1,3,3,1,0}, {1,3,1,3,0}, {2,2,2,2,1,3,1,3}, 
          {1,3,3,1,0}, {2,2,2,2,3,1,3,1}, {2,2,2,2,1,3,1,3}, 
          {2,2,2,2,2,2,2,2,2,2,2,2}};

static const int DatVertexFaceLen[] =
        {  3, 3, 3, 3, 3, 3,
           3, 3, 5, 5, 5, 5, 8,
           5, 5, 5, 8, 5, 5, 8,
           5, 5, 8, 5, 8, 8, 12};

static const int DatFaceEdge[][REFHEXAREGMAXN_EpF] =
      { {0,8,11,7}, {2,9,8,1}, {4,10,9,3}, {6,11,10,5}, /* 0-3 */
        {0,18,21,17}, {12,19,18,1}, {14,20,19,13}, {16,21,20,15}, /* 4-7 */
        {2,26,29,12}, {22,27,26,3}, {24,28,27,23}, {13,29,28,25}, /* 8-11 */
        {4,34,37,22}, {30,35,34,5}, {32,36,35,31}, {23,37,36,33}, /* 12-15 */
        {6,40,43,30}, {17,41,40,7}, {38,42,41,16}, {31,43,42,39}, /* 16-19 */
        {15,44,47,38}, {25,45,44,14}, {33,46,45,24}, {39,47,46,32}, /* 20-23 */
        {8,52,48,18}, {34,50,52,10}, {46,53,50,36}, {20,48,53,44}, /* 24-27 */
        {11,52,51,40}, {26,49,52,9}, {45,53,49,28}, {42,51,53,47}, /* 28-31 */
        {21,48,51,41}, {29,49,48,19}, {37,50,49,27}, {43,51,50,35}}; /* 32-35 */

static const int DatEdgeFace[][REFHEXAREGMAXN_FpE] =
      { {0,4}, {1,5}, {1,8}, {2,9}, {2,12}, {3,13}, {3,16}, {0,17},
        {0,1,24}, {1,2,29}, {2,3,25}, {0,3,28}, {5,8}, {6,11}, {6,21},
        {7,20}, {7,18}, {4,17}, {4,5,24}, {5,6,33}, {6,7,27}, {4,7,32},
        {9,12}, {10,15}, {10,22}, {11,21}, {8,9,29}, {9,10,34}, {10,11,30},
        {8,11,33}, {13,16}, {14,19}, {14,23}, {15,22}, {12,13,25}, {13,14,35},
        {14,15,26}, {12,15,34}, {18,20}, {19,23}, {16,17,28}, {17,18,32},
        {18,19,31}, {16,19,35}, {20,21,27}, {21,22,30}, {22,23,26}, {20,23,31},
        {24,27,32,33}, {29,30,33,34}, {25,26,34,35},
        {28,31,32,35}, {24,25,28,29}, {26,27,30,31}};

static const int DatEdgeFaceIndex[][REFHEXAREGMAXN_FpE] =
        { {0,0}, {3,3}, {0,0}, {3,3}, {0,0}, {3,3}, {0,0}, {3,3},  /* 0-7 */
          {1,2,0}, {1,2,3}, {1,2,3}, {2,1,0}, {0,3}, {3,0}, {0,3}, /* 8 -14 */
          {3,0}, {0,3}, {3,0}, {1,2,3}, {1,2,3}, {1,2,0}, {2,1,0}, /* 15-21 */
          {0,3}, {3,0}, {0,3}, {3,0}, {1,2,0}, {1,2,3}, {1,2,3},   /* 22-28 */
          {2,1,0}, {0,3}, {3,0}, {0,3}, {3,0}, {1,2,0}, {1,2,3},   /* 29-35 */
          {1,2,3}, {2,1,0}, {0,3}, {3,0}, {1,2,3}, {1,2,3},        /* 36-41 */
          {1,2,0}, {2,1,0}, {1,2,3}, {1,2,0}, {1,2,0}, {2,1,3},    /* 42-47 */
          {2,1,1,2}, {1,2,1,2}, {1,2,1,2},   /* 48-50 */         
          {2,1,2,1}, {1,2,1,2}, {1,2,1,2}};  /* 51-53 */

static const int DatEdgeFaceLen[] =
        {  2, 2, 2, 2, 2, 2, 2, 2,
           3, 3, 3, 3, 2, 2, 2,
           2, 2, 2, 3, 3, 3, 3,
           2, 2, 2, 2, 3, 3, 3,
           3, 2, 2, 2, 2, 3, 3,
           3, 3, 2, 2, 3, 3,
           3, 3, 3, 3, 3, 3,
           4, 4, 4, 4, 4, 4};

static const int DatNewVertexEqOldVertex[] = { 0, 1, 2, 3, 4, 5, 6, 7};
static const int DatNewVertexEqOldVertexIndex[] = { 0, 1, 2, 3, 4, 5, 6, 7};

static const int DatInteriorVertexOfCell[] = {26};

static const double DatOldFaceNewVertexPos[][REFHEXAREGMAXN_nVpoF]
                                            [REFHEXAREGMAXN_oVpoF] =
   { { {1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1}, {0.5,0.5,0,0},
       {0,0.5,0.5,0}, {0,0,0.5,0.5}, {0.5,0,0,0.5}, {0.25, 0.25, 0.25, 0.25} },
     { {1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1}, {0.5,0.5,0,0},
       {0,0.5,0.5,0}, {0,0,0.5,0.5}, {0.5,0,0,0.5}, {0.25, 0.25, 0.25, 0.25} },
     { {1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1}, {0.5,0.5,0,0},
       {0,0.5,0.5,0}, {0,0,0.5,0.5}, {0.5,0,0,0.5}, {0.25, 0.25, 0.25, 0.25} },
     { {1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1}, {0.5,0.5,0,0},
       {0,0.5,0.5,0}, {0,0,0.5,0.5}, {0.5,0,0,0.5}, {0.25, 0.25, 0.25, 0.25} },
     { {1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1}, {0.5,0.5,0,0},
       {0,0.5,0.5,0}, {0,0,0.5,0.5}, {0.5,0,0,0.5}, {0.25, 0.25, 0.25, 0.25} },
     { {1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1}, {0.5,0.5,0,0},
       {0,0.5,0.5,0}, {0,0,0.5,0.5}, {0.5,0,0,0.5}, {0.25, 0.25, 0.25, 0.25} }};


static const double DatPositionOfIntVert[][HEXAN_V] =
        { { 0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125 } };

static const int DatInteriorEdgeOfCell[] = 
      {  48, 49, 50, 51, 52, 53 };

static const int DatInteriorFaceOfCell[] = 
      {  24, 25, 26, 27, 28, 29,
         30, 31, 32 ,33, 34, 35};

static const int DatInteriorVertexOfEdge[][REFHEXAREGMAXN_iVpE] =
      { {8}, {9}, {10}, {11}, {15}, {13}, 
        {17}, {20}, {14}, {18}, {21}, {23}};
static const int DatInteriorVertexOfEdgeLen[]=
      {  1, 1, 1, 1, 1, 1, 
         1, 1, 1, 1, 1, 1};

static const int DatInteriorVertexOfFace[][REFHEXAREGMAXN_iVpF] =
      { {12}, {16}, {19}, {22}, {24}, {25}};
static const int DatInteriorVertexOfFaceLen[] =
      { 1, 1, 1, 1, 1, 1};

static const int DatInteriorEdgeOfFace[][REFHEXAREGMAXN_iEpF] =
      { {8,9,10,11}, {18,19,20,21}, {26,27,28,29},
        {34,35,36,37}, {40,41,42,43}, {44,45,46,47}};
static const int DatInteriorEdgeOfFaceLen[] =
      { 4, 4, 4, 4, 4, 4};

static const int DatOldEdgeNewVertex[][REFHEXAREGMAXN_nVpoE] =
      { {0,8,1}, {1,9,2}, {2,10,3}, {3,11,0}, {0,15,4}, {1,13,5}, 
        {2,17,6}, {3,20,7}, {4,14,5}, {5,18,6}, {6,21,7}, {7,23,4}};
static const int DatOldEdgeNewVertexLen[] =
      { 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};

static const int DatOldEdgeNewEdge[][REFHEXAREGMAXN_nEpoE] =
      {{0,1}, {2,3}, {4,5}, {6,7}, {17,16}, {12,13}, {22,23}, 
       {30,31}, {15,14}, {25,24}, {33,32}, {39,38}};
static const int DatOldEdgeNewEdgeLen[] =
      { 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};

static const int DatOldFaceNewVertex[][REFHEXAREGMAXN_nVpoF] =
      { {0,1,2,3,8,9,10,11,12}, {0,4,5,1,15,14,13,8,16},
        {1,5,6,2,13,18,17,9,19}, {2,6,7,3,17,21,20,10,22},
        {0,3,7,4,11,20,23,15,24}, {4,7,6,5,23,21,18,14,25}};

static const int DatOldFaceNewVertexLen[] =
      { 9, 9, 9, 9, 9, 9};

static const int DatOldFaceNewEdge[][REFHEXAREGMAXN_nEpoF] =
  { {0,1,2,3,4,5,6,7,8,9,10,11}, {17,16,15,14,13,12,1,0,21,20,19,18},
    {12,13,25,24,23,22,3,2,29,28,27,26}, {22,23,33,32,31,30,5,4,37,36,35,34},
    {7,6,30,31,39,38,16,17,40,43,42,41}, {38,39,32,33,24,25,14,15,47,46,45,44}};

static const int DatOldFaceNewEdgeLen[] =
      { 12, 12, 12, 12, 12, 12};

static const int DatOldFaceNewFace[][REFHEXAREGMAXN_nFpoF] =
      { {0,1,2,3}, {4,7,6,5}, {8,11,10,9}, 
        {12,15,14,13}, {17,16,19,18}, {20,23,22,21}};

static const int DatOldFaceNewFaceLen[] =
      { 4, 4, 4, 4, 4, 4};

static const int DatNewEdgeOldEdge[] = 
      {  0,  0,  1,  1,  2,  2,  3,  3, -1, -1,
        -1, -1,  5,  5,  8,  8,  4,  4, -1, -1,
        -1, -1,  6,  6,  9,  9, -1, -1, -1, -1,
         7,  7, 10, 10, -1, -1, -1, -1, 11, 11,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1,-1, -1, -1};
                
static const int DatNewFaceOldFace[] =
      { 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };

static const int DatOldFaceNewLocFace[][HEXAN_F] =
       { { 0, 1, -1, -1, 4, -1}, { 0, 4, 1, -1, -1, -1}, 
         { 0, -1, 4, 1, -1, -1}, { 0, -1, -1, 4, 1, -1},
         { -1, 4, -1, -1, 1, 0}, { -1, 1, 4, -1, -1, 0},
         { -1, -1, 1, 4, -1, 0}, { -1, -1, -1, 1, 4, 0}};

static const int DatChildTwistIndex[] = 
       { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

// Constructor
TRefHexaRegDesc::TRefHexaRegDesc(TShapeDesc *shape) : TRefDesc(shape)
{

  Type = HexaReg;

  //set all numbers
  N_Vertices = 27;
  N_Edges = 54;
  N_Faces = 36;
  N_Children = 8;
  N_InnerVertices = 1;
  N_NewVertEqOldVert = 8;
  N_InnerEdges = 6;
  N_InnerFaces = 12;

  // initialize all dimension values
  MaxN_VpC = REFHEXAREGMAXN_VpC;
  MaxN_CpV = REFHEXAREGMAXN_CpV;
  MaxN_EpC = REFHEXAREGMAXN_EpC;
  MaxN_CpE = REFHEXAREGMAXN_CpE;
  MaxN_EpV = REFHEXAREGMAXN_EpV;
  MaxN_EpF = REFHEXAREGMAXN_EpF;
  MaxN_FpE = REFHEXAREGMAXN_FpE;
  MaxN_VpF = REFHEXAREGMAXN_VpF;
  MaxN_FpV = REFHEXAREGMAXN_FpV;
  MaxN_FpC = REFHEXAREGMAXN_FpC;
  MaxN_CpF = REFHEXAREGMAXN_CpF;
  MaxN_iVpE = REFHEXAREGMAXN_iVpE;
  MaxN_iVpF = REFHEXAREGMAXN_iVpF;
  MaxN_iEpF = REFHEXAREGMAXN_iEpF;
  MaxN_nVpoE = REFHEXAREGMAXN_nVpoE;
  MaxN_nEpoE = REFHEXAREGMAXN_nEpoE;
  MaxN_nVpoF = REFHEXAREGMAXN_nVpoF;
  MaxN_oVpoF = REFHEXAREGMAXN_oVpoF;
  MaxN_nEpoF = REFHEXAREGMAXN_nEpoF;
  MaxN_nFpoF = REFHEXAREGMAXN_nFpoF;

  // initialize all pointers
  ChildType = (const Shapes *) DatChildType;
  FaceType = (const Refinements *) DatFaceType;
  EdgeType = (const Refinements *) DatEdgeType;

  ChildVertex = (const int *) DatChildVertex;
  VertexChild = (const int *) DatVertexChild;
  VertexChildIndex = (const int *) DatVertexChildIndex;
  VertexChildLen = (const int *) DatVertexChildLen;

  ChildFace = (const int *) DatChildFace;
  FaceChild = (const int *) DatFaceChild;
  FaceChildIndex = (const int *) DatFaceChildIndex;
  FaceChildLen = (const int *) DatFaceChildLen;

  ChildEdge = (const int *) DatChildEdge;
  EdgeChild = (const int *) DatEdgeChild;
  EdgeChildIndex = (const int *) DatEdgeChildIndex;
  EdgeChildLen = (const int *) DatEdgeChildLen;

  EdgeVertex = (const int *) DatEdgeVertex;
  VertexEdge = (const int *) DatVertexEdge;
  VertexEdgeIndex = (const int *) DatVertexEdgeIndex;
  VertexEdgeLen = (const int *) DatVertexEdgeLen;

  FaceVertex = (const int *) DatFaceVertex;
  VertexFace = (const int *) DatVertexFace;
  VertexFaceIndex = (const int *) DatVertexFaceIndex;
  VertexFaceLen = (const int *) DatVertexFaceLen;

  FaceEdge = (const int *) DatFaceEdge;
  EdgeFace = (const int *) DatEdgeFace;
  EdgeFaceIndex = (const int *) DatEdgeFaceIndex;
  EdgeFaceLen = (const int *) DatEdgeFaceLen;

  NewVertexEqOldVertex = (const int *) DatNewVertexEqOldVertex;
  NewVertexEqOldVertexIndex = (const int *) DatNewVertexEqOldVertexIndex;

  InteriorVertexOfCell = (const int *) DatInteriorVertexOfCell;
  PositionOfIntVert = (const double *) DatPositionOfIntVert;

  InteriorEdgeOfCell = (const int *) DatInteriorEdgeOfCell;
  InteriorFaceOfCell = (const int *) DatInteriorFaceOfCell;
  InteriorVertexOfEdge = (const int *) DatInteriorVertexOfEdge;
  InteriorVertexOfEdgeLen = (const int *) DatInteriorVertexOfEdgeLen;
  InteriorVertexOfFace = (const int *) DatInteriorVertexOfFace;
  InteriorVertexOfFaceLen = (const int *) DatInteriorVertexOfFaceLen;
  InteriorEdgeOfFace = (const int *) DatInteriorEdgeOfFace;
  InteriorEdgeOfFaceLen = (const int *) DatInteriorEdgeOfFaceLen;
                
  OldEdgeNewVertex = (const int *) DatOldEdgeNewVertex;
  OldEdgeNewVertexLen = (const int *) DatOldEdgeNewVertexLen;
 
  OldEdgeNewEdge = (const int *) DatOldEdgeNewEdge;
  OldEdgeNewEdgeLen = (const int *) DatOldEdgeNewEdgeLen;
  NewEdgeOldEdge = (const int *) DatNewEdgeOldEdge;
  NewFaceOldFace = (const int *) DatNewFaceOldFace;

  OldFaceNewVertex = (const int *) DatOldFaceNewVertex;
  OldFaceNewVertexPos = (const double *) DatOldFaceNewVertexPos;
  OldFaceNewVertexLen = (const int *) DatOldFaceNewVertexLen;
  OldFaceNewEdge = (const int *) DatOldFaceNewEdge;
  OldFaceNewEdgeLen = (const int *) DatOldFaceNewEdgeLen;
  OldFaceNewFace = (const int *) DatOldFaceNewFace;
  OldFaceNewFaceLen = (const int *) DatOldFaceNewFaceLen;

  OldFaceNewLocFace = (const int *) DatOldFaceNewLocFace;
  ChildTwistIndex = (const int *) DatChildTwistIndex;
}

// Methods
