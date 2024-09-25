/*

 $Id: parm.h,v 1.3 2009/05/08 23:17:35 rhuey Exp $

 AutoGrid 

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.

 AutoGrid is a Trade Mark of The Scripps Research Institute.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 */

#define C3      0
#define C2      1
#define C1      2
#define Cac     3
#define Cpl     4
    
#define N3pl    5
#define Nox     6
#define N3      7
#define Ntr     8
#define Npl     9
#define N1     10
#define Nam    11

#define O3     12
#define O2     13
#define Om     14

#define S3pl   15
#define S3     16
#define S2     17
#define Sac    18
#define Sox    19
#define S      20

#define HC     21 
#define H      22 

#define P3     23
#define Pac    24 
#define Pox    25

#define B      26
#define Bac    27
#define Box    28

#define Al     29
#define As     30
#define Be     31
#define Br     32
#define Ca     33
#define Cl     34
#define Cu     35  
#define Fl     36  
#define Fe     37 
#define Ge     38  
#define I      39 
#define K      40 
#define Li     41 
#define Mg     42  
#define Mn     43
#define Na     44
#define Ni     45
#define Pb     46
#define Si     47
#define Zn     48

#define num_atom_types     49
#define num_hbnd_types     16


extern float get_Rij(int type1, int type2);
extern float get_epsij(int type1, int type2);

