// default & command line options

#include <iostream>
#include <string>
#include <boost/lexical_cast.hpp>

struct options {
  unsigned int length;
  double temperature;
  unsigned int sweeps;
  unsigned int therm;
  bool valid;

  options(unsigned int argc, char *argv[], unsigned int len_def, double temp_def) :
    length(len_def), temperature(temp_def), sweeps(1 << 16), therm(sweeps >> 3), valid(true) {
    for (unsigned int i = 1; i < argc; ++i) {
      switch (argv[i][0]) {
      case '-' :
        switch (argv[i][1]) {
        case 'l' :
          if (++i == argc) { usage(); return; }
          length = boost::lexical_cast<unsigned int>(argv[i]); break;
        case 't' :
          if (++i == argc) { usage(); return; }
          temperature = boost::lexical_cast<double>(argv[i]); break;
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
    if (temperature <= 0 || sweeps == 0) {
      std::cerr << "invalid parameter\n"; usage(); return;
    }
    std::cout << "System Linear Size     = " << length << '\n'
              << "Temperature            = " << temperature << '\n'
              << "MCS for Thermalization = " << therm << '\n'
              << "MCS for Measurement    = " << sweeps << '\n';
  }

  void usage(std::ostream& os = std::cerr) {
    os << "[command line options]\n"
       << "  -l int    System Linear Size\n"
       << "  -t double Temperature\n"
       << "  -m int    MCS for Thermalization\n"
       << "  -n int    MCS for Measurement\n"
       << "  -h        this help\n";
    valid = false;
  }
};
