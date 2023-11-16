//
//  Copyright (C) 2013-2014 Greg Landrum and NextMove Software
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <cstring>
#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/FileParseException.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/FileParserUtils.h>
#include <RDGeneral/LocaleSwitcher.h>
#include "ProximityBonds.h"

#include <GraphMol/MonomerInfo.h>

// PDBWriter support multiple "flavors" of PDB output
// flavor & 1 : Ignore atoms in alternate conformations and dummy atoms
// flavor & 2 : Read each MODEL into a separate molecule.


/* Nawrocki: Reading HIN files. */

namespace RDKit {

namespace {

Atom *HINAtomFromSymbol(const char *symb) {
  PRECONDITION(symb, "bad char ptr");
  if (symb[0] == 'D' && !symb[1]) {
    auto *result = new Atom(1);
    result->setIsotope(2);
    return result;
  } else if (symb[0] == 'T' && !symb[1]) {
    auto *result = new Atom(1);
    result->setIsotope(3);
    return result;
  }
  int elemno = PeriodicTable::getTable()->getAtomicNumber(symb);
  printf("Element number: %i\n",elemno);
  return elemno > 0 ? new Atom(elemno) : (Atom *)nullptr;
}

void HINAtomLine(RWMol *mol, const char *ptr, unsigned int len,
                 unsigned int flavor, std::map<int, Atom *> &amap) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(ptr, "bad char ptr");
  std::string tmp;

  Atom *atom = (Atom *)nullptr;
  char symb[3];

  symb[0] = ptr[12];
  printf("Symbol: %s\n",symb);
  atom = HINAtomFromSymbol(symb);
  mol->addAtom(atom, true, true);
                 }

void parseHINBlock(RWMol *&mol, const char *str, bool sanitize, bool removeHs,
                   unsigned int flavor, bool proximityBonding) {
  PRECONDITION(str, "bad char ptr");
  std::map<int, Atom *> amap;
  std::map<Bond *, int> bmap;
  Utils::LocaleSwitcher ls;
  bool multi_conformer = false;
  int conformer_atmidx = 0;
  Conformer *conf = nullptr;

  while (*str) {
    unsigned int len;
    const char *next = str + 1;
    for (;;) {
      if (*next == '\r') {
        len = (unsigned int)(next - str);
        if (next[1] == '\n') {
          next += 2;
        } else {
          next++;
        }
        break;
      } else if (*next == '\n') {
        len = (unsigned int)(next - str);
        next++;
        break;
      } else if (*next == '\0') {
        len = (unsigned int)(next - str);
        break;
      }
      next++;
    }

    // ATOM records

    if (str[0] == 'a' && str[1] == 't' && str[2] == 'o' && str[3] == 'm' &&
        str[4] == ' ' && str[5] == ' ') {
      if (!multi_conformer) {
        if (!mol) {
          mol = new RWMol();
        }
        printf(str);
        HINAtomLine(mol, str, len, flavor, amap);
      }
      // HETATM records
    }
    str = next;
  }
                   }


} // namespace

RWMol *HINBlockToMol(const char *str, bool sanitize, bool removeHs,
                     unsigned int flavor, bool proximityBonding) {
  RWMol *mol = nullptr;
  try {
    parseHINBlock(mol, str, sanitize, removeHs, flavor, proximityBonding);
  } catch (...) {
    delete mol;
    throw;
  }

  return mol;
}

RWMol *HINDataStreamToMol(std::istream *inStream, bool sanitize, bool removeHs,
                          unsigned int flavor, bool proximityBonding) {
  PRECONDITION(inStream, "bad stream");
  std::string buffer;
  while (!inStream->eof() && !inStream->fail()) {
    std::string line;
    std::getline(*inStream, line);
    buffer += line;
    buffer += '\n';
    auto ptr = line.c_str();
    // Check for END

    printf("%s\n",ptr);

    if (ptr[0] == 'e' && ptr[1] == 'n' && ptr[2] == 'd' && ptr[3] == 'm' && ptr[4] == 'o' && ptr[5] == 'l')
         {
      printf("endmol found - breake\n");
      break;
    }
  }
  return HINBlockToMol(buffer.c_str(), sanitize, removeHs, flavor,
                       proximityBonding);
}

RWMol *HINFileToMol(const std::string &fileName, bool sanitize, bool removeHs,
                    unsigned int flavor, bool proximityBonding) {
  std::ifstream ifs(fileName.c_str(), std::ios_base::binary);
  if (!ifs || ifs.bad()) {
    std::ostringstream errout;
    errout << "Bad input file " << fileName;
    throw BadFileException(errout.str());
  }

  return HINDataStreamToMol(static_cast<std::istream *>(&ifs), sanitize,
                            removeHs, flavor, proximityBonding);
}
}  // namespace RDKit
