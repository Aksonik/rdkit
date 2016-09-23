# copyright 2000 greg landrum
""" Generic file manipulation stuff

"""
from __future__ import print_function
import numpy
import string, re


class ReFile:
  """convenience class for dealing with files with comments

  blank (all whitespace) lines, and lines beginning with comment
    characters are skipped.

  anything following a comment character on a line is stripped off
  """

  def readline(self):
    """ read the next line and return it.

    return '' on EOF

    """
    result = ''
    while result == '':
      inLine = self.inFile.readline()
      if inLine == '':
        return ''
      result = string.strip(self.regExp.split(inLine)[0])
    return result

  def readlines(self):
    """ return a list of all the lines left in the file

    return [] if there are none

    """
    res = []
    inLines = self.inFile.readlines()
    for line in inLines:
      result = string.strip(self.regExp.split(line)[0])
      if result != '':
        res.append(result)

    return res

  def rewind(self):
    """ rewinds the file (seeks to the beginning)

    """
    self.inFile.seek(0)

  def __init__(self, fileName, mode='r', comment=r'#', trailer=r'\n'):
    if trailer is not None and trailer != '':
      comment = comment + r'|' + trailer
    self.regExp = re.compile(comment)
    self.inFile = open(fileName, mode)


def ReadDataFile(fileName, comment=r'#', depVarCol=0, dataType=numpy.float):
  """ read in the data file and return a tuple of two Numeric arrays:
  (independant variables, dependant variables).

  **ARGUMENTS:**

  - fileName: the fileName

  - comment: the comment character for the file

  - depVarcol: the column number containing the dependant variable

  - dataType: the Numeric short-hand for the data type

  RETURNS:

   a tuple of two Numeric arrays:

    (independant variables, dependant variables).
  
  """
  inFile = ReFile(fileName)
  dataLines = inFile.readlines()
  nPts = len(dataLines)

  if dataType in [numpy.float, numpy.float32, numpy.float64]:
    _convfunc = float
  else:
    _convfunc = int

  nIndVars = len(string.split(dataLines[0])) - 1
  indVarMat = numpy.zeros((nPts, nIndVars), dataType)
  depVarVect = numpy.zeros(nPts, dataType)
  for i in range(nPts):
    splitLine = string.split(dataLines[i])
    depVarVect[i] = _convfunc(splitLine[depVarCol])
    del splitLine[depVarCol]
    indVarMat[i, :] = map(_convfunc, splitLine)

  return indVarMat, depVarVect


if __name__ == '__main__':
  import sys

  fileN = sys.argv[1]
  iV, dV = ReadDataFile(fileN)
  print('iV:', iV)
  print('dV:', dV)
