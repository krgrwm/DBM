#include "plist.h"

PList::PList(const int size) : plist(size) {}

void PList::append(const int i, const PosVal pv) {
  this->plist[i] = pv;
}

int PList::size() const {
  return this->plist.size();
}

double PList::p(const int i) const {
  return this->plist[i].second;
}

Pos PList::pos(const int i) const {
  return this->plist[i].first;
}

PosVal PList::at(const int i) const {
  return this->plist[i];
}
