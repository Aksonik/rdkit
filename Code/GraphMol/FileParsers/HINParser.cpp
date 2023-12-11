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

/* Parsing HIN files by Grzegorz Nawrocki. */

namespace RDKit {

namespace {

    /* pairBondType - an atom pairwise array (atom numbers from zero) for storing bond types. */

    const int rows=100;
    const int cols=100;

    std::vector<std::vector<std::string>> pairBondType(rows, std::vector<std::string>(cols));

constexpr int BCNAM(char A, char B, char C) { return (A << 16) | (B << 8) | C; }

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

  int serialno;
  tmp = std::string(ptr + 6, 5);
  try {
    serialno = FileParserUtils::toInt(tmp);
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Non-integer PDB serial number " << tmp;
    throw FileParseException(errout.str());
  }

  Atom *atom = (Atom *)nullptr;
  char symb[3];

  /* The HIN format is flexible, i.e. order of columns matters, not their exact postion.
     Let's check which colum starts where. */

  std::string currstr=" ";
  std::string prevstr=" ";
  int colnum=0;
  int colarr[35][2]={0};

  std::cout << "Parse HIN line: " << std::string(ptr,len) << std::endl;

  for (int charnum=0; charnum<len; charnum+=1){
    currstr=std::string(ptr+charnum,1);
    if(charnum>0){
      prevstr=std::string(ptr+charnum-1,1);
    }

    if((prevstr==" ") & (currstr!=" ")){
      colnum+=1;
      colarr[colnum][0]=charnum;
    }
    else if((prevstr!=" ") & (currstr==" ")){
      colarr[colnum][1]=charnum;
      std::cout << "column " << colnum << " starts " << colarr[colnum][0] << " ends " << colarr[colnum][1] << std::endl;
    }
  }

  /*
     column:
     1  - keyword, e.g. atom
     2  - index (starts with 1)
     3  - label
     4  - element, e.g. C
     5  - type (according to a FF)
     6  - flat, e.g. h for heteroatom
     7  - partial charge
     8  - x
     9  - y
     10 - z
     11 - number of bonds
     12 - index of a bonded atom
     13 - bond type (s, d, t, a)
     etc.
  */

  if(ptr[colarr[4][0]]>='A' && ptr[colarr[4][0]]<='Z'){
    symb[0]=ptr[colarr[4][0]];
    if(ptr[colarr[4][0]+1]>='a' && ptr[colarr[4][0]+1]<='z'){
      symb[1]=ptr[colarr[4][0]+1];
      symb[2]='\0';
    }
    else{
      symb[1]='\0';
    }
  }
  printf("Symbol: %s\n",symb);
  atom = HINAtomFromSymbol(symb);

  if (!atom) {
    std::ostringstream errout;
    errout << "Cannot determine element for PDB atom #" << serialno;
    throw FileParseException(errout.str());
  }

  mol->addAtom(atom, true, true);
  amap[serialno] = atom;
  printf("Atom serial number: %i\n",serialno);

  if (len >= colarr[8][1]) {
    RDGeom::Point3D pos;
    try {
      pos.x = FileParserUtils::toDouble(std::string(ptr + colarr[8][0], colarr[8][1]-colarr[8][0]));
      if (len >= colarr[9][1]) {
        pos.y = FileParserUtils::toDouble(std::string(ptr + colarr[9][0], colarr[9][1]-colarr[9][0]));
      }
      if (len >= colarr[10][1]) {
        pos.z = FileParserUtils::toDouble(std::string(ptr + colarr[10][0], colarr[10][1]-colarr[10][0]));
      }
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Problem with coordinates for HIN atom #" << serialno;
      throw FileParseException(errout.str());
    }

    Conformer *conf;
    if (!mol->getNumConformers()) {
      conf = new RDKit::Conformer(mol->getNumAtoms());
      conf->set3D(pos.z != 0.0);
      conf->setId(0);
      mol->addConformer(conf, false);
    } else {
      conf = &mol->getConformer();
      if (pos.z != 0.0) {
        conf->set3D(true);
      }
    }
    printf("Coor: %f %f %f\n",pos.x,pos.y,pos.z);
    printf("Index: %i\n",atom->getIdx());
    conf->setAtomPos(atom->getIdx(), pos);
  }

  printf("Len: %i\n",len);
  if (len >= 0) {
    double charge = 0.000;
    charge=FileParserUtils::toDouble(std::string(ptr + colarr[7][0], colarr[7][1]-colarr[7][0]));

    if (charge != 0) {
      atom->setFormalCharge(charge);
      printf("Formal charge: %f\n",charge);
    }
  }

  int bondNum=std::stoi(std::string(ptr+colarr[11][0],1));
  std::cout << "Number of bonds: " << bondNum << std::endl;

  int atomA=std::stoi(std::string(ptr+colarr[2][0],colarr[2][1]-colarr[2][0]));
  int atomB=0;

  for (int atomBndx=12;atomBndx<12+2*bondNum;atomBndx+=2){
    int atomB=std::stoi(std::string(ptr+colarr[atomBndx][0],colarr[atomBndx][1]-colarr[atomBndx][0]));

    std::string bondABType=std::string(ptr+colarr[atomBndx+1][0],colarr[atomBndx+1][1]-colarr[atomBndx+1][0]);
    char bondABTypeChar=bondABType[0];

    std::cout << "Bond type: " << atomA << "-" << atomB << " " << bondABType << std::endl;

    pairBondType[atomA-1][atomB-1]=bondABTypeChar;
  }

  tmp=std::string(ptr+colarr[5][0],colarr[5][1]-colarr[5][0]);
  std::cout << "Atom type: " << tmp << std::endl;
  AtomPDBResidueInfo *info = new AtomPDBResidueInfo(tmp, serialno);
  atom->setMonomerInfo(info);

                 }

void PDBBondLine(RWMol *mol, const char *ptr, unsigned int len,
                 std::map<int, Atom *> &amap, std::map<Bond *, int> &bmap) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(ptr, "bad char ptr");

  for (int i=0; i<rows; i+=1){
    for (int j=i+1; j<cols; j+=1){
        if(pairBondType[i][j]=="s"){
           Bond *bond=mol->getBondBetweenAtoms(amap[i+1]->getIdx(),amap[j+1]->getIdx());
           bond=new Bond(Bond::SINGLE);
           bond->setOwningMol(mol);           
           bond->setBeginAtom(amap[i+1]);
           bond->setEndAtom(amap[j+1]);
           mol->addBond(bond,true);
          }
        else if(pairBondType[i][j]=="d"){
           Bond *bond=mol->getBondBetweenAtoms(amap[i+1]->getIdx(),amap[j+1]->getIdx());
           bond=new Bond(Bond::DOUBLE);
           bond->setOwningMol(mol);           
           bond->setBeginAtom(amap[i+1]);
           bond->setEndAtom(amap[j+1]);
           mol->addBond(bond,true);
          }
        else if(pairBondType[i][j]=="t"){
           Bond *bond=mol->getBondBetweenAtoms(amap[i+1]->getIdx(),amap[j+1]->getIdx());
           bond=new Bond(Bond::TRIPLE);
           bond->setOwningMol(mol);           
           bond->setBeginAtom(amap[i+1]);
           bond->setEndAtom(amap[j+1]);
           mol->addBond(bond,true);
          }
        else if(pairBondType[i][j]=="a"){
           Bond *bond=mol->getBondBetweenAtoms(amap[i+1]->getIdx(),amap[j+1]->getIdx());
           bond=new Bond(Bond::AROMATIC);
           bond->setOwningMol(mol);           
           bond->setBeginAtom(amap[i+1]);
           bond->setEndAtom(amap[j+1]);
           mol->addBond(bond,true);
          }
      }
  }
}

void PDBTitleLine(RWMol *mol, const char *ptr, unsigned int len) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(ptr, "bad char ptr");
  std::string title;
  while (ptr[len - 1] == ' ') {
    len--;
  }
  if (ptr[len - 1] == ';') {
    len--;
  }
  if (len > 21 && !strncmp(ptr + 10, " MOLECULE: ", 11)) {
    title = std::string(ptr + 21, len - 21);
  } else if (len > 10) {
    title = std::string(ptr + 10, len - 10);
  }
  if (!title.empty()) {
    mol->setProp(common_properties::_Name, title);
  }
}

void PDBConformerLine(RWMol *mol, const char *ptr, unsigned int len,
                      Conformer *&conf, int &conformer_atmidx) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(ptr, "bad char ptr");

  if (len >= 38) {
    RDGeom::Point3D pos;
    try {
      pos.x = FileParserUtils::toDouble(std::string(ptr + 30, 8));
      if (len >= 46) {
        pos.y = FileParserUtils::toDouble(std::string(ptr + 38, 8));
      }
      if (len >= 54) {
        pos.z = FileParserUtils::toDouble(std::string(ptr + 46, 8));
      }
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Problem with multi-conformer coordinates";
      throw FileParseException(errout.str());
    }

    if (conformer_atmidx == 0) {
      conf = new RDKit::Conformer(mol->getNumAtoms());
      conf->setId(mol->getNumConformers());
      conf->set3D(pos.z != 0.0);
      mol->addConformer(conf, false);
    } else if (pos.z != 0.0) {
      conf->set3D(true);
    }

    if (conformer_atmidx < rdcast<int>(mol->getNumAtoms())) {
      conf->setAtomPos(conformer_atmidx, pos);
      conformer_atmidx++;
    }
  }
}

// This function determines whether a standard atom name in
// in a recognized PDB amino acid should be chiral or not.
// This is used to avoid chirality on VAL.CG and LEU.CG.
bool StandardPDBChiralAtom(const char *resnam, const char *atmnam) {
  switch (BCNAM(resnam[0], resnam[1], resnam[2])) {
    case BCNAM('G', 'L', 'Y'):
      return false;
    case BCNAM('I', 'L', 'E'):
    case BCNAM('T', 'H', 'R'):
      // Alpha and beta carbons (" CA " and " CB ").
      return atmnam[0] == ' ' && atmnam[1] == 'C' &&
             (atmnam[2] == 'A' || atmnam[2] == 'B') && atmnam[3] == ' ';
    case BCNAM('A', 'L', 'A'):
    case BCNAM('A', 'R', 'G'):
    case BCNAM('A', 'S', 'N'):
    case BCNAM('A', 'S', 'P'):
    case BCNAM('C', 'Y', 'S'):
    case BCNAM('G', 'L', 'N'):
    case BCNAM('G', 'L', 'U'):
    case BCNAM('H', 'I', 'S'):
    case BCNAM('L', 'E', 'U'):
    case BCNAM('L', 'Y', 'S'):
    case BCNAM('M', 'E', 'T'):
    case BCNAM('P', 'H', 'E'):
    case BCNAM('P', 'R', 'O'):
    case BCNAM('S', 'E', 'R'):
    case BCNAM('T', 'R', 'P'):
    case BCNAM('T', 'Y', 'R'):
    case BCNAM('V', 'A', 'L'):
      return atmnam[0] == ' ' && atmnam[1] == 'C' && atmnam[2] == 'A' &&
             atmnam[3] == ' ';
  }
  return false;
}

void StandardPDBResidueChirality(RWMol *mol) {
  for (ROMol::AtomIterator atomIt = mol->beginAtoms();
       atomIt != mol->endAtoms(); ++atomIt) {
    Atom *atom = *atomIt;
    if (atom->getChiralTag() != Atom::CHI_UNSPECIFIED) {
      auto *info = (AtomPDBResidueInfo *)atom->getMonomerInfo();
      if (info && info->getMonomerType() == AtomMonomerInfo::PDBRESIDUE &&
          !info->getIsHeteroAtom() &&
          !StandardPDBChiralAtom(info->getResidueName().c_str(),
                                 info->getName().c_str())) {
        atom->setChiralTag(Atom::CHI_UNSPECIFIED);
        if (atom->hasProp(common_properties::_CIPCode)) {
          atom->clearProp(common_properties::_CIPCode);
        }
      }
    }
  }
}

void BasicPDBCleanup(RWMol &mol) {
  ROMol::VERTEX_ITER atBegin, atEnd;
  boost::tie(atBegin, atEnd) = mol.getVertices();
  while (atBegin != atEnd) {
    Atom *atom = mol[*atBegin];
    atom->calcExplicitValence(false);

    // correct four-valent neutral N -> N+
    // This was github #1029
    if (atom->getAtomicNum() == 7 && atom->getFormalCharge() == 0 &&
        atom->getExplicitValence() == 4) {
      atom->setFormalCharge(1);
    }
    ++atBegin;
  }
}

void initPairBondType(){

    /* Initialize the bond type array with "x", i.e. no bond. */
    
    for (int i=0; i<rows; ++i) {
        for (int j=0; j<cols; ++j) {
            pairBondType[i][j]="x";
        }
    }
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

  initPairBondType();

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

    // atom records
    if (str[0] == 'a' && str[1] == 't' && str[2] == 'o' && str[3] == 'm' &&
        str[4] == ' ' && str[5] == ' ') {
      if (!multi_conformer) {
        if (!mol) {
          mol = new RWMol();
        }
        HINAtomLine(mol, str, len, flavor, amap);
      } else {
        PDBConformerLine(mol, str, len, conf, conformer_atmidx);
      }
      // HETATM records
    } else if (str[0] == 'H' && str[1] == 'E' && str[2] == 'T' &&
               str[3] == 'A' && str[4] == 'T' && str[5] == 'M') {
      if (!multi_conformer) {
        if (!mol) {
          mol = new RWMol();
        }
        HINAtomLine(mol, str, len, flavor, amap);
      } else {
        PDBConformerLine(mol, str, len, conf, conformer_atmidx);
      }
      // CONECT records
    } else if (str[0] == 'C' && str[1] == 'O' && str[2] == 'N' &&
               str[3] == 'E' && str[4] == 'C' && str[5] == 'T') {
      if (mol && !multi_conformer) {
        PDBBondLine(mol, str, len, amap, bmap);
      }
    } else if (str[0] == 'e' && str[1] == 'n' && str[2] == 'd' &&
               str[3] == 'm' && str[4] == 'o' && str[5] == 'l') {
      if (mol && !multi_conformer) {
        PDBBondLine(mol, str, len, amap, bmap);
      }
      // COMPND records
    } else if (str[0] == 'C' && str[1] == 'O' && str[2] == 'M' &&
               str[3] == 'P' && str[4] == 'N' && str[5] == 'D') {
      if (!mol) {
        mol = new RWMol();
      }
      if (len > 10 &&
          (str[9] == ' ' || !strncmp(str + 9, "2 MOLECULE: ", 12))) {
        PDBTitleLine(mol, str, len);
      }
      // HEADER records
    } else if (str[0] == 'H' && str[1] == 'E' && str[2] == 'A' &&
               str[3] == 'D' && str[4] == 'E' && str[5] == 'R') {
      if (!mol) {
        mol = new RWMol();
      }
      PDBTitleLine(mol, str, len < 50 ? len : 50);
      // ENDMDL records
    } else if (str[0] == 'E' && str[1] == 'N' && str[2] == 'D' &&
               str[3] == 'M' && str[4] == 'D' && str[5] == 'L') {
      if (!mol) {
        break;
      }
      multi_conformer = true;
      conformer_atmidx = 0;
      conf = nullptr;
    }
    str = next;
  }

  if (!mol) {
    return;
  }

  if (proximityBonding) {
    ConnectTheDots(mol, ctdIGNORE_H_H_CONTACTS);
  }
  // flavor & 8 doesn't encode double bonds
  if (proximityBonding || flavor & 8) {
    StandardPDBResidueBondOrders(mol);
  }

  BasicPDBCleanup(*mol);

  if (sanitize) {
    if (removeHs) {
      MolOps::removeHs(*mol, false, false);
    } else {
      MolOps::sanitizeMol(*mol);
    }
  } else {
    // we need some properties for the chiral setup
    mol->updatePropertyCache(false);
  }

  /* Set tetrahedral chirality from 3D co-ordinates */
  MolOps::assignChiralTypesFrom3D(*mol);
  StandardPDBResidueChirality(mol);
}
}  // namespace

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

RWMol *HINBlockToMol(const std::string &str, bool sanitize, bool removeHs,
                     unsigned int flavor, bool proximityBonding) {
  return HINBlockToMol(str.c_str(), sanitize, removeHs, flavor,
                       proximityBonding);
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

    printf("HIN line: %s\n",ptr);

    // Check for END
    if (ptr[0] == 'e' && ptr[1] == 'n' && ptr[2] == 'd' && ptr[3] == 'm' && ptr[4] == 'o' && ptr[5] == 'l')
         {
      break;
    }
  }
  return HINBlockToMol(buffer.c_str(), sanitize, removeHs, flavor,
                       proximityBonding);
}
RWMol *HINDataStreamToMol(std::istream &inStream, bool sanitize, bool removeHs,
                          unsigned int flavor, bool proximityBonding) {
  return HINDataStreamToMol(&inStream, sanitize, removeHs, flavor,
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
