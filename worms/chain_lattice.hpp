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

#pragma once

class chain_lattice {
public:
  chain_lattice(unsigned int L) : length_(L) {}
  unsigned int coordination_num() const { return 2; }
  unsigned int num_sites() const { return length_; }
  unsigned int num_bonds() const { return num_sites(); }
  unsigned int source(unsigned int b) const { return b; }
  unsigned int target(unsigned int b) const { return (b == length_-1) ? 0 : b+1; }
  double phase(unsigned int s) const { return (s & 1) ? 1.0 : -1.0; }
private:
  unsigned int length_;
};
