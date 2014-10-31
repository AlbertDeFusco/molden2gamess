#!/usr/bin/env python
import sys
import math

#%%%%%%%%%%%
#
# three clases are used to generate a GAMESS input
# Parameters
#   just units and title
#
# Atoms
#   stores coordinates and basis set
#
# Orbitals
#   will store a $VEC keyword



# Store calculation parameters
# derived from the contents of the
# Molden file
class Parameters(object):
  def setUnits(self,units):
    self.units=units
  def setTitle(self,title):
    self.title=title
  def setNumMOs(self,numMOs):
    self.numMOs=numMOs

# Atoms have coordinates and a basis set
class Atom(object):
  def __init__(self,name,number):
    self.name=name
    self.number=number
    #a basis set is a list of
    #contracted gaussian primitives
    self.basis=list()

  #adding two atoms together sums
  #the number of primitive GTOs
  def __add__(self,other):
    try:
      return len(self.basis)+len(other.basis)
    except:
      return 0.0

  def __radd__(self,other):
    try:
      return other + len(self.basis)
    except:
      return self

  def setCoord(self,coord):
    self.coord=coord

  def addGTOs(self,shell,prims):
    self.basis.append( (shell,prims) )

  def printAtom(self):
    print "%-8s %5.1f %20.15f %20.15f %20.15f" % (
        self.name,
        self.number,
        self.coord[0],self.coord[1],self.coord[2] )

  def printBasis(self):
    for shell in self.basis:
      print "%2s %5d" % (shell[0].upper(),len(shell[1]))
      for num,GTO in enumerate(shell[1]):
        if shell[0] == "l":
          print "%4d %20.10f %20.10f %20.10f" % (num+1,GTO[0],GTO[1],GTO[2])
        else:
          print "%4d %20.10f %20.10f" % (num+1,GTO[0],GTO[1])
    print

  def numGTOs(self):
    self.gtoMap = {}
    numGTOs=0
    for shell in self.basis:
      if shell[0] == 'l':
        numPrim=4
      elif shell[0] == 's':
        numPrim=1
      elif shell[0] == 'p':
        numPrim=3
      elif shell[0] == 'd':
        numPrim=6

      if shell[0] in self.gtoMap.keys():
        self.gtoMap[shell[0]].append( (numGTOs,numGTOs+numPrim) )
      else:
        self.gtoMap[shell[0]] = [ (numGTOs,numGTOs+numPrim) ]

      numGTOs=numGTOs+numPrim

    return numGTOs


class Orbitals(object):
  def __init__(self):
    #Separate lists of alpha and beta orbitals are
    #kept. If orbitals are doubly occupied they are
    #stored in alpha and beta is empty
    #Each entry in alpha or beta begins with the MO
    #energy followed by the coefficients
    self.alpha=list()
    self.beta=list()

  def addAlpha(self,moEnergy,moCoeffs):
    self.alpha.append( (moEnergy,moCoeffs) )
  def addBeta(self,moEnergy,moCoeffs):
    self.beta.append( (moEnergy,moCoeffs) )

  def printMOEnergies(self):
    if(self.beta == []):
      print "! MO Energies"
      alphaE = [ ("%10.5f " % orb[0]) for orb in self.alpha]
      print(" ".join(alphaE))
    else:
      print "! Alpha MO Energies"
      alphaE = [ ("%10.5f " % orb[0]) for orb in self.alpha]
      print(" ".join(alphaE))
      print "! Beta MO Energies"
      betaE = [ ("%10.5f " % orb[0]) for orb in self.beta]
      print(" ".join(betaE))

  def printMOs(self):
    self.printMOCoefficients(self.alpha)
    if self.beta != []:
      self.printMOCoefficients(self.beta)

  def printMOCoefficients(self,orbs,orbWidth=5):
    for num,orb in enumerate(orbs):
      numLines=len(orb[1])/orbWidth
      remainder=len(orb[1])%orbWidth
      for line in range(numLines):
        head= "%2d%3d" % ((num+1)%100,line+1)
        s = line*orbWidth
        e = s + orbWidth
        thisLine = [ ("%15.8E" % myOrb) for myOrb in orb[1][s:e] ]
        print(head+"".join(thisLine))
      if remainder>0:
        head= "%2d%3d" % ((num+1)%100,line+2)
        thisLine = [ ("%15.8E" % myOrb) for myOrb in orb[1][e:e+remainder] ]
        print(head+"".join(thisLine))

  def makeMap(self,atoms):
    #generate a look-up table of shell type and primitive GTOs
    #in the MO coefficients
    self.gtoMap = {}
    totalGTOs=0
    for num,atom in enumerate(atoms):
      nGTOs = atom.numGTOs()
      for shell in atom.gtoMap.keys():
        for entry in atom.gtoMap[shell]:
          start=entry[0]+totalGTOs
          end=entry[1]+totalGTOs
          if shell in self.gtoMap.keys():
            self.gtoMap[shell].append( (start,end) )
          else:
            self.gtoMap[shell] = [ (start,end) ]
      totalGTOs=totalGTOs+nGTOs

  def fixD(self):
    for num,orb in enumerate(self.alpha):
      for s,e in self.gtoMap['d']:
        orb[1][s:e] = [c*math.sqrt(3) for c in orb[1][s:e]]


#%%%%%%%%%%%%
#
# these routines read the Molden file
# and populate the objects
#
#


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

def readOrbitals(atoms,params,inputFile):
  #each orbital will have exactly numGTOs entries
  numGTOs=sum( [atom.numGTOs() for atom in atoms] )

  gamessOrbitals = Orbitals()
  gamessOrbitals.makeMap(atoms)


  count=0
  moldenOrbitals = grabGroup(inputFile,"MO")
  string=moldenOrbitals[1][0]
  for num,line in enumerate(moldenOrbitals):
    if line[0] == string:
      # I'm not certain how many different options I can have
      if string == "Sym=":
        energy = float(moldenOrbitals[num+1][1])
        spin   = moldenOrbitals[num+2][1]
        orb    = [float(coeff[1]) for coeff in moldenOrbitals[num+4:num+4+numGTOs]]
      elif string == "Ene=":
        energy = float(moldenOrbitals[num][1])
        spin   = moldenOrbitals[num+1][1]
        orb    = [float(coeff[1]) for coeff in moldenOrbitals[num+3:num+3+numGTOs]]

      if spin == "Alpha":
        gamessOrbitals.addAlpha(energy,orb)
      if spin == "Beta":
        gamessOrbitals.addBeta(energy,orb)

      count=count+1
  params.setNumMOs(len(gamessOrbitals.alpha))

  gamessOrbitals.fixD()
  return gamessOrbitals

#%%%%%%%%%%%%
#
# Each group in the GAMESS input
# gets its own function
#
#

def printData(params,atoms):
  print " $DATA"
  print params.title
  print "C1"
  for atom in atoms:
    atom.printAtom()
    atom.printBasis()
  print " $END"

def printContrl(params):
  print " $CONTRL ispher=XX scftyp=XX coord=unique units=%4s $END" % params.units

def printGuess(params):
  print " $GUESS guess=moread norb=%-d prtmo=.f. purify=.t. $END" % params.numMOs


def printVec(orbitals):
  print " $VEC"
  orbitals.printMOs()
  print " $END"



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
  orbitals=readOrbitals(atoms,params,molden)

  printContrl(params)
  printGuess(params)
  printData(params,atoms)
  printVec(orbitals)

  print >> sys.stderr, "REMINDER: fix SCFTYP and ISPHER"

