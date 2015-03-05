/*	schottky.cc

	Copyright Danny Calegari and Alden Walker

	released under the terms of the GNU GPL, version 3.0
	see attached licence for more details
*/


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
    
    if (std::string(argv[1]) == "-4db" && argc == 7) {
      cpx center( atof(argv[2]), atof(argv[3]) );
      double rad = atof(argv[4]);
      double step = atof(argv[5]);
      int depth = atof(argv[6]);
      int num_boxes = (2.0*rad)/step;
      std::vector<std::vector<std::vector<std::vector<bool> > > > res(0);
      for (double a1=center.real()-rad; a1 < center.real()+rad; a1 += step) {
        res.push_back(std::vector<std::vector<std::vector<bool> > >(0));
        for (double b1=center.imag()-rad; b1 < center.imag()+rad; b1 += step) {
        res.back().push_back(std::vector<std::vector<bool> >(0));
          for (double a2=center.real()-rad; a2 < center.real()+rad; a2 += step) {
          res.back().back().push_back(std::vector<bool>(0));
            for (double b2=center.imag()-rad; b2 < center.imag()+rad; b2 += step) {
              ifs G( cpx(a1,b1), cpx(a2,b2), 1, 0);
              int diff;
              //std::cout << "Testing " << cpx(a1,b1) << " " << cpx(a2,b2) << "\n";
              res.back().back().back().push_back(G.is_connected(depth, diff));
            }
          }
        }
      }
      for (int i=0; i<(int)res.size(); ++i) {
        for (int j=0; j<(int)res[i].size(); ++j) {
          for (int k=0; k<(int)res[i][j].size(); ++k) {
            for (int ell=0; ell<(int)res[i][j][k].size(); ++ell) {
              std::cout << "{" << i << "," << j << "," << k << "," << ell << "," << res[i][j][k][ell] << "},";
            }
          }
        }
      }
      return 0;
    }
    
              
  }
  
  //std::cout << sizeof(long long int) <<"\n";
  
 std::cout << "enter 'i' for IFS or 'm' for mandelbrot (n for new GUI):";
 std::cin >> c;
  if (c == 'n') {
    IFSGui G;
    G.launch(BOTH, cpx(-0.58, 0.33));
  } else {
    mode = (c=='i' ? 0 : 1);
    ifs G(cos(TWOPI/3.0)+I*sin(TWOPI/3.0), 0.5, w, mode);      // default value : Sierpinski triangle
    G.user_interface();
  }
  
  //IFSGui G;
  //G.launch(BOTH, cpx(-0.58, 0.33));
  
  return 0;
}
