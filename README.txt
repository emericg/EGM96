This archive contains translation to C of an EGM96 implementation.
The original files are located at 
http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/egm96.html

Files EGM96 and CORCOEF are just filtered with sed -e"s/D/E/g"
EGM96 and CORRCOEF from that site.

f477.c is rewritten in C f477.f. once I was interested in C source for EGM96 
and couldn't find it with http://www.google.com. here it is.

Terms of use for original software are not clear enough, as you can
see in f477.c (it keeps original notes); this translation 
(i.e. the minimal part of work that can be honestly 
licensed by D.Ineiev) is distributed under the next conditions:

(C) 2006 D.Ineiev

  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.

