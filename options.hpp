/*
   worms: a simple worm code

   Copyright (C) 2013-2021 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

// default & command line options

#include <iostream>
#include <string>

struct options {
  unsigned int L;
  double T;
  double H;
  unsigned int sweeps;
  unsigned int therm;
  bool valid;

  options(unsigned int argc, char *argv[], unsigned int L_def, double T_def) :
    L(L_def), T(T_def), H(0), sweeps(1 << 16), therm(sweeps >> 3), valid(true) {
    for (unsigned int i = 1; i < argc; ++i) {
      switch (argv[i][0]) {
      case '-' :
        switch (argv[i][1]) {
        case 'L' :
          if (++i == argc) { usage(); return; }
          L = std::atoi(argv[i]); break;
        case 'T' :
          if (++i == argc) { usage(); return; }
          T = std::atof(argv[i]); break;
        case 'H' :
          if (++i == argc) { usage(); return; }
          H = std::atof(argv[i]); break;
        case 'm' :
          if (++i == argc) { usage(); return; }
          therm = std::atoi(argv[i]); break;
        case 'n' :
          if (++i == argc) { usage(); return; }
          sweeps = std::atoi(argv[i]); break;
        case 'h' :
          usage(std::cout); return;
        default :
          usage(); return;
        }
        break;
      default :
        usage(); return;
      }
    }
    if (T <= 0 || sweeps == 0) {
      std::cerr << "invalid parameter\n"; usage(); return;
    }
    std::cout << "System Linear Size     = " << L << '\n'
              << "Temperature            = " << T << '\n'
              << "Magnetic Field         = " << H << '\n'
              << "MCS for Thermalization = " << therm << '\n'
              << "MCS for Measurement    = " << sweeps << '\n';
  }

  void usage(std::ostream& os = std::cerr) {
    os << "[command line options]\n"
       << "  -L int    System Linear Size\n"
       << "  -T double Temperature\n"
       << "  -H double Magnetic Field\n"
       << "  -m int    MCS for Thermalization\n"
       << "  -n int    MCS for Measurement\n"
       << "  -h        this help\n";
    valid = false;
  }
};
