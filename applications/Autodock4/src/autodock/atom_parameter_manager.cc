/*

 $Id: atom_parameter_manager.cc,v 1.5 2010/10/01 22:51:39 mp Exp $

 AutoDock 

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.

 AutoDock is a Trade Mark of The Scripps Research Institute.

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

#define MAXKEY (256*256)

#include <stdlib.h>
#include <string.h>
#include "structs.h" // needed for ParameterEntry structure
#include "atom_parameter_manager.h"

typedef ParameterEntry PE;

static PE *dictionary[MAXKEY];

static unsigned int 
hash(const char key[]) {
    switch (strlen(key)) {
        case 0: return 0;
        case 1: return (unsigned int)key[0];
        default: return (unsigned int)key[0] + 256*(unsigned int)key[1];
    }
}

void 
apm_enter(const char key[], const PE& value) {
    if (dictionary[hash(key)] == NULL) {
        dictionary[hash(key)] = (PE *) calloc(1, sizeof(PE));
    }
    *(dictionary[hash(key)]) = value;  // this replaces, as well as inserts
    return;
}

PE * 
apm_find(const char key[]) {
    return dictionary[hash(key)];
}
