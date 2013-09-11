// default & command line options

#include <iostream>
#include <string>
#include <boost/lexical_cast.hpp>

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
          L = boost::lexical_cast<unsigned int>(argv[i]); break;
        case 'T' :
          if (++i == argc) { usage(); return; }
          T = boost::lexical_cast<double>(argv[i]); break;
        case 'H' :
          if (++i == argc) { usage(); return; }
          H = boost::lexical_cast<double>(argv[i]); break;
        case 'm' :
          if (++i == argc) { usage(); return; }
          therm = boost::lexical_cast<unsigned int>(argv[i]); break;
        case 'n' :
          if (++i == argc) { usage(); return; }
          sweeps = boost::lexical_cast<unsigned int>(argv[i]); break;
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
