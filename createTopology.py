from time import sleep
import decimal
from itertools import combinations

def extract_numbers (line):
	lineString = line.replace ("\t", ",").replace (" ", ",").replace ("\n", "")
	for item in lineString.split (","):
		try:
			yield decimal.Decimal (item)
		except:
			pass

def readAtomInfo (inputFileName, atomInfo):
	with open ("atomEntries.testing", "r") as inputFile:
		for line in inputFile:
			lineArray = list (extract_numbers (line))
			atomInfo.append ({'sino': int (lineArray[0]), 'molType': int (lineArray[1]), 'atomType': int (lineArray[2]), 'charge': int (lineArray[3]), 'x': float (lineArray[4]), 'y': float (lineArray[5]), 'z': float (lineArray[6]), 'bondAtom1': int (lineArray[7]), 'bondAtom2': int (lineArray[8]), 'bondAtom3': int (lineArray[9]), 'bondAtom4': int (lineArray[10])})

	return atomInfo

def createBonds (atomInfo, bondInfo):
	lineArray = []
	sino_bond = 1
	bondTypeArr = []
	bondTypeArr.append({'atom1': 0, 'atom2': 0, 'bondType': 0})

	def bondCheck (atom1, atom2, bondInfo):
		for bond in bondInfo:
			if ((atom1 == bond['bondAtom1'] or atom1 == bond['bondAtom2']) and (atom2 == bond['bondAtom1'] or atom2 == bond['bondAtom2'])):
				return 0

		return 1

	def findBondType (atom1, atom2, bondTypeArr):
		returnBondType = 0
		if (atom1 < atom2):
			ascAtom1 = atom1
			ascAtom2 = atom2
		else:
			ascAtom1 = atom2
			ascAtom2 = atom1

		for items in bondTypeArr:
			if (ascAtom1 == items ['atom1'] and ascAtom2 == items ['atom2']):
				returnBondType = items ['bondType']

		if (returnBondType == 0):
			existingLength = len (bondTypeArr)
			bondTypeArr.append ({'atom1': ascAtom1, 'atom2': ascAtom2, 'bondType': existingLength})
			returnBondType = existingLength

		return returnBondType, bondTypeArr

	def addBond (atomLine, atomInfo, dictString, bondInfo, sino_bond, bondTypeArr):
		if (atomLine [dictString] and bondCheck (atomLine ['sino'], atomLine [dictString], bondInfo)):
			try:
				secondAtomType = atomInfo [int (atomLine [dictString]) - 1]['atomType']
				bondType, bondTypeArr = findBondType (atomLine ['atomType'], secondAtomType, bondTypeArr)
				if (atomLine ['sino'] < atomLine [dictString]):
					bondInfo.append ({'sino': sino_bond, 'bondType': bondType, 'bondAtom1': atomLine ['sino'], 'bondAtom2': atomLine [dictString], 'bondAtom1Type': atomLine ['atomType'], 'bondAtom2Type': secondAtomType})
					sino_bond += 1
				else:
					bondInfo.append ({'sino': sino_bond, 'bondType': bondType, 'bondAtom1': atomLine [dictString], 'bondAtom2': atomLine ['sino'], 'bondAtom1Type': secondAtomType, 'bondAtom2Type': atomLine ['atomType']})
					sino_bond += 1
			except:
				pass

		return bondInfo, bondTypeArr, sino_bond

	for atomLine in atomInfo:
		bondInfo, bondTypeArr, sino_bond = addBond (atomLine, atomInfo, 'bondAtom1', bondInfo, sino_bond, bondTypeArr)
		bondInfo, bondTypeArr, sino_bond = addBond (atomLine, atomInfo, 'bondAtom2', bondInfo, sino_bond, bondTypeArr)
		bondInfo, bondTypeArr, sino_bond = addBond (atomLine, atomInfo, 'bondAtom3', bondInfo, sino_bond, bondTypeArr)
		bondInfo, bondTypeArr, sino_bond = addBond (atomLine, atomInfo, 'bondAtom4', bondInfo, sino_bond, bondTypeArr)

	return bondInfo

def createAngles (atomInfo, angleInfo, bondInfo):
	# angleInfo contains sino, angleType, angleAtom1, angleAtom2, angleAtom3, angleAtom1Type, angleAtom2Type, angleAtom3Type
	angleTypeArr = []
	sino_angle = 1

	# Initializing connectedAtoms dict of lists
	connectedAtoms = {}
	for atomLine in atomInfo:
		connectedAtoms [atomLine ['sino']] = []

	angleTypeArr.append ({'atom1': 0, 'atom2': 0, 'atom3': 0, 'angleType': 0})

	def findConnectedAtoms (concernedBondAtom, primaryConnect, bondInfo, connectedAtoms):
		if (primaryConnect not in connectedAtoms [concernedBondAtom]):
			connectedAtoms [concernedBondAtom].append (primaryConnect)
		if (concernedBondAtom not in connectedAtoms [primaryConnect]):
			connectedAtoms [primaryConnect].append (concernedBondAtom)

		for bondLine in bondInfo:
			if ((bondLine ['bondAtom1'] == concernedBondAtom) and (bondLine ['bondAtom2'] not in connectedAtoms [concernedBondAtom])):
				connectedAtoms [concernedBondAtom].append (bondLine ['bondAtom2'])
				connectedAtoms [concernedBondAtom].sort ()
			if ((bondLine ['bondAtom2'] == concernedBondAtom) and (bondLine ['bondAtom1'] not in connectedAtoms [concernedBondAtom])):
				connectedAtoms [concernedBondAtom].append (bondLine ['bondAtom1'])
				connectedAtoms [concernedBondAtom].sort ()

		return connectedAtoms

	for bondLine in bondInfo:
		connectedAtoms = findConnectedAtoms (bondLine ['bondAtom1'], bondLine ['bondAtom2'], bondInfo, connectedAtoms)
	
	# print (connectedAtoms)

	# Creating angleInfo
	for atom in connectedAtoms:
		if (len (connectedAtoms [atom]) > 1):
			comb = combinations (connectedAtoms [atom], 2)
			for i in list (comb):
				firstAtomType = atomInfo [int (i [0]) - 1]['atomType']
				secondAtomType = atomInfo [int (atom) - 1]['atomType']
				thirdAtomType = atomInfo [int (i [1]) - 1]['atomType']
				angleInfo.append ({'sino': sino_angle, 'angleType': 0, 'angleAtom1': i[0], 'angleAtom2': atom, 'angleAtom3': i[1], 'angleAtom1Type': firstAtomType, 'angleAtom2Type': secondAtomType, 'angleAtom3Type': thirdAtomType})
			# print (atom, connectedAtoms [atom])

	for angle in angleInfo:
		print (angle)
		sleep (1)

	return angleInfo

def main():
	atomInfo = []
	bondInfo = []
	angleInfo = []
	dihedralInfo = []

	atomInfo = readAtomInfo ("atomEntries.testing", atomInfo)
	bondInfo = createBonds (atomInfo, bondInfo)
	angleInfo = createAngles (atomInfo, angleInfo, bondInfo)
	# dihedralInfo = createDihedrals (dihedralInfo, angleInfo, bondInfo)

	# for bond in bondInfo:
	# 	print (bond)

if __name__ == '__main__':
	main()