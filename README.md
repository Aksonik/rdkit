# RDKit

#### Branch: reading-HIN-format

It has an additional function to read a molecule written in a HIN format.  
It is coded in a similar way as Chem.rdmolfiles.MolFromPDBFile to read PDB file.

import rdkit  
print(rdkit.__version__)

from rdkit import Chem  
from rdkit.Chem import Descriptors

mol=Chem.rdmolfiles.MolFromHINFile("Ethane.hin")  
mass=Descriptors.MolWt(mol)  
print(mass)

print(mol.GetNumAtoms())

for conf in mol.GetConformers():  
     print(conf.GetPositions())

distMat = Chem.Get3DDistanceMatrix(supplHIN)  
print(distMat[0,2])

for atom in mol.GetAtoms():  
  atom_charge = atom.GetFormalCharge()  
  print("Atom (symbol, atomic numer, valence(imp), valence(exp), charge:",atom.GetSymbol(),atom.GetAtomicNum(),atom.GetImplicitValence(),atom.GetExplicitValence(),atom_charge)  
  print("Neighbors:",[x.GetAtomicNum() for x in atom.GetNeighbors()])  

print("Number of bonds:",mol.GetNumBonds())
  
for bond in mol.GetBonds():  
  aid1=bond.GetBeginAtomIdx()  
  aid2=bond.GetEndAtomIdx()  
  
  atom=mol.GetAtomWithIdx(aid1)  
  aid1_symbol=atom.GetSymbol()

  atom=mol.GetAtomWithIdx(aid2)  
  aid2_symbol=atom.GetSymbol()
