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
   

#ifndef _FEPC_EASY_HELPER
#define _FEPC_EASY_HELPER

#include "funktion.h"

func_p
func_multi(func_p f, func_p g);

func_p
func_div(func_p f, func_p g);

folgen_vektor_p
folgen_vektor_multi(folgen_vektor_p f, folgen_vektor_p g, int step);

folgen_vektor_p
folgen_vektor_div(folgen_vektor_p f, folgen_vektor_p g, int step);

fepc_real_t folge_norm2(folge_p folge, int step);

folge_p
folge_multi(folge_p f, folge_p g, int step);

folge_p
folge_div(folge_p f, folge_p g, int step);

vec_real_p get_mean_points(vec_p v, vec_p grad, int step);

#endif
