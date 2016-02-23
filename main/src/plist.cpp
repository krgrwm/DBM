#include "plist.h"

PList::PList(const int size) : poslist(size), plist(size) {}

void PList::append(const int i, const Pos &pos, const double p) {
  this->poslist[i] = pos;
  this->plist[i] = p;
}

int PList::size() const {
  return this->plist.size();
}

double PList::p(const int i) const {
  return this->plist[i];
}

Pos PList::pos(const int i) const {
  return this->poslist[i];
}

PosVal PList::at(const int i) const {
  return PosVal(this->poslist[i], this->plist[i]);
}

const std::vector<double> PList::get_plist() {
  return this->plist;
}
