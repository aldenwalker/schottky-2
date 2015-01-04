/*	schottky.cc

	Version 0.01 March 15 2014

	Copyright Danny Calegari

	released under the terms of the GNU GPL, version 3.0
	see attached licence for more details
*/

// standard libraries to include

#include <iostream>
#include <cstdlib>
#include <cmath>

#include "ifs.h"
#include "ifs_gui.h"

#define TWOPI 6.28318530718


int main(int argc, char *argv[]) { 
  char c;
  int w = 1024;
  cpx I(0,1);
  int mode = 0;
  
  if (argc > 1) {
    if (std::string(argv[1]) == "-a" && argc == 5) {
      double a1 = atof(argv[2]);
      double a2 = atof(argv[3]);
      int trap_depth = atoi(argv[4]);
      IFSGui G;
      G.mand_trap_depth = trap_depth;
      G.mand_ll = (1.0/sqrt(2.0))*cpx( cos(a2), sin(a1) );
      G.mand_ur = (1.0/sqrt(2.0))*cpx( cos(a1), sin(a2) );
      G.find_traps_along_circle_in_window(2, false);
      return 0;
    }
  }
  
  //std::cout << sizeof(long long int) <<"\n";
  
//  std::cout << "enter 'i' for IFS or 'm' for mandelbrot (n for new GUI):";
//  std::cin >> c;
//   if (c == 'n') {
//     IFSGui G;
//     G.launch(BOTH, cpx(-0.58, 0.33));
//   } else {
//     mode = (c=='i' ? 0 : 1);
//     ifs G(cos(TWOPI/3.0)+I*sin(TWOPI/3.0), 0.5, w, mode);      // default value : Sierpinski triangle
//     G.user_interface();
//   }
  
  IFSGui G;
  G.launch(BOTH, cpx(-0.58, 0.33));
  
  return 0;
}
