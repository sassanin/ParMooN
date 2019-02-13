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
   
#ifndef __FALTUNG_H
#define __FALTUNG_H

#include "faltung_hilfe.h"

/* Berechnet die Projektion der Faltung der Funktionen f und g auf die Gitterstruktur von w.
 Von f und g sind die Gitterstrukturen sowie die dazugehoerigen Eintraege gegeben. Von w ist
 nur die Gitterstruktur gegeben.

 Eingabewerte: Funktionen f,g mit Gitterstruktur und Eintraegen (dh. zugehoerigen Folgen)
							Funktion w mit Gitterstruktur ohne Eintraege
 Ausgabewert:	Faltung von f und g, projeziert auf die Gitterstruktur von w */


/* Berechnung der FEPC-Faltung */
void
faltung_fepc(func_p0_p f, func_p0_p g, func_p0_p w, fepc_real_t h, func_p0_p back);

void
faltung_fepc(func_p1_p f, func_p1_p g, func_p1_p w, fepc_real_t h, func_p1_p back);

#endif
