#!/bin/csh -f
#
# $Id: paramdat2h.csh,v 1.10 2011/03/08 04:18:37 mp Exp $
# 
# AutoDock 
# 
# Copyright (C) 1989-2007,  Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson, 
# All Rights Reserved.
# 
# AutoDock is a Trade Mark of The Scripps Research Institute.
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

echo 'const char *param_string_4_0[MAX_LINES] = {'
egrep -v '^#|^$' $1 | sed 's/\(.*\)$/"\1\\n", /'
echo ' };'
echo 'const char *param_string_4_1[MAX_LINES] = {'
egrep -v '^#|^$' $2 | sed 's/\(.*\)$/"\1\\n", /'
echo ' };'
echo '// EOF'
