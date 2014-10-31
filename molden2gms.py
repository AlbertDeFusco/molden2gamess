#!/usr/bin/env python
import sys
import math

# Store calculation parameters
# derived from the contents of the
# Molden file
class Parameters(object):
  def setUnits(self,units):
    self.units=units
  def setTitle(self,title):
    self.title=title

# Atoms have coordinates and a basis set
class Atom(object):
  def __init__(self,name,number):
    self.name=name
    self.number=number
    #a basis set is a list of
    #contracted gaussian primitives
    self.basis=list()

  def setCoord(self,coord):
    self.coord=coord

  def addGTOs(self,shell,prims):
    self.basis.append( (shell,prims) )

  def printAtom(self):
    print "%-8s %5.1f %15.8f %15.8f %15.8f" % (
        self.name,
        self.number,
        self.coord[0],self.coord[1],self.coord[2] )

  def printBasis(self):
    for contraction in self.basis:
      print "%2s %5d" % (contraction[0].upper(),len(contraction[1]))
      for num,GTO in enumerate(contraction[1]):
        if contraction[0] == "l":
          print "%4d %20.10f %20.10f %20.10f" % (num+1,GTO[0],GTO[1],GTO[2])
        else:
          print "%4d %20.10f %20.10f" % (num+1,GTO[0],GTO[1])
    print

class Orbital(object):
  def __init__(self):
    self.alpha=list()
    self.beta=list()

# Extract a Molden format group
# to a list for further processing.
# Each item of the list is a list
# of text columns.
def grabGroup(inputFile,groupName):
  group=list()
  for num,line in enumerate(inputFile):
    if "["+groupName+"]" in line:
      group.append(line.split())
      for line2 in inputFile[num+1:]:
        if "[" in line2:
          break
        else:
          group.append(line2.split())

  return group


# return a list of Atoms
def readCoord(params,inputFile):
  moldenAtoms = grabGroup(inputFile,"Atoms")

  if moldenAtoms[0][1] == "AU":
    params.setUnits("bohr")
  elif moldenAtoms[0][1] == "Angs":
    params.setUnits("angs")

  atoms=list()
  for atom in moldenAtoms[1:]:
    thisAtom = Atom(atom[0].upper(), float(atom[2]))
    thisAtom.setCoord([float(atom[3]),float(atom[4]),float(atom[5])])
    atoms.append(thisAtom)

  return atoms

# add the basis set to the list of Atoms
def readBasis(atoms,inputFile):
  moldenBasis = grabGroup(inputFile,"GTO")
  #search for a basis set for each atom
  for iAtom,atom in enumerate(atoms):
    for num,line in enumerate(moldenBasis[1:]):
      if line != [] and ( line[0] == str(iAtom+1) and line[1]=='0' ):
        #basis found, now read shells until a blank line is found
        for num2,line2 in enumerate(moldenBasis[num+1:]):
          if line2==[]:
            break
          elif line2[0] == 'sp':
            shell="l"
            length=int(moldenBasis[num+num2+1][1])
            firstGTO=num+num2+2
            lastGTO=num+num2+2+length
            GTOs=[(float(exp),float(coeffS),float(coeffP)) for exp,coeffS,coeffP in moldenBasis[firstGTO:lastGTO]]
            atom.addGTOs(shell,GTOs)
          elif line2[0] in ['s','p','d']:
            shell=moldenBasis[num+num2+1][0]
            length=int(moldenBasis[num+num2+1][1])
            firstGTO=num+num2+2
            lastGTO=num+num2+2+length
            GTOs=[(float(exp),float(coeff)) for exp,coeff in moldenBasis[firstGTO:lastGTO]]
            atom.addGTOs(shell,GTOs)
          else:
            continue

def printData(params,atoms):
  print " $DATA"
  print params.title
  print "C1"
  for atom in atoms:
    atom.printAtom()
    atom.printBasis()
  print " $END"

def printContrl(params):
  print " $CONTRL coord=unique units=%4s $END" % params.units





if __name__ == "__main__":
  try:
    moldenFile = sys.argv[1]
  except:
    print "Usage %s <molden-file>" % sys.argv[0]
    sys.exit(1)

  #read the whole file into a list
  molden = open(moldenFile).readlines()

  params=Parameters()
  try:
    params.setTitle(grabGroup(inputFile,"Title")[1][0])
  except:
    params.setTitle("molden2gamess")

  atoms=readCoord(params,molden)
  readBasis(atoms,molden)

  printContrl(params)
  printData(params,atoms)


