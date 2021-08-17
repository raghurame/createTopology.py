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
	atomTypeArr = []
	with open ("atomEntries.testing", "r") as inputFile:
		for line in inputFile:
			lineArray = list (extract_numbers (line))
			atomInfo.append ({'sino': int (lineArray[0]), 'molType': int (lineArray[1]), 'atomType': int (lineArray[2]), 'charge': int (lineArray[3]), 'x': float (lineArray[4]), 'y': float (lineArray[5]), 'z': float (lineArray[6]), 'bondAtom1': int (lineArray[7]), 'bondAtom2': int (lineArray[8]), 'bondAtom3': int (lineArray[9]), 'bondAtom4': int (lineArray[10])})

	for atomLine in atomInfo:
		if (atomLine ['atomType'] not in atomTypeArr):
			atomTypeArr.append (atomLine ['atomType'])

	return atomInfo, atomTypeArr

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

	return bondInfo, bondTypeArr

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
	
	def findAngleType (firstAtomType, secondAtomType, thirdAtomType, angleTypeArr):
		if (firstAtomType <= thirdAtomType):
			ascFirstAtomType = firstAtomType
			ascThirdAtomType = thirdAtomType
		else:
			ascFirstAtomType = thirdAtomType
			ascThirdAtomType = firstAtomType

		assignedType = 0

		for types in angleTypeArr:
			if (ascFirstAtomType == types ['atom1'] and secondAtomType == types ['atom2'] and ascThirdAtomType == types ['atom3']):
				assignedType = types ['angleType']

		if (assignedType == 0):
			assignedType = len (angleTypeArr)
			angleTypeArr.append ({'atom1': ascFirstAtomType, 'atom2': secondAtomType, 'atom3': ascThirdAtomType, 'angleType': len (angleTypeArr)})

		return assignedType, angleTypeArr

	# Creating angleInfo
	for atom in connectedAtoms:
		if (len (connectedAtoms [atom]) > 1):
			comb = combinations (connectedAtoms [atom], 2)
			for i in list (comb):
				firstAtomType = atomInfo [int (i [0]) - 1]['atomType']
				secondAtomType = atomInfo [int (atom) - 1]['atomType']
				thirdAtomType = atomInfo [int (i [1]) - 1]['atomType']
				angleType, angleTypeArr = findAngleType (firstAtomType, secondAtomType, thirdAtomType, angleTypeArr)
				angleInfo.append ({'sino': sino_angle, 'angleType': angleType, 'angleAtom1': i [0], 'angleAtom2': atom, 'angleAtom3': i [1], 'angleAtom1Type': firstAtomType, 'angleAtom2Type': secondAtomType, 'angleAtom3Type': thirdAtomType})
				sino_angle += 1

	return angleInfo, angleTypeArr

def createDihedrals (atomInfo, dihedralInfo, angleInfo, bondInfo):
	# dihInfo contains sino, dihType, dihAtom1, dihAtom2, dihAtom3, dihAtom4, dihAtom1Type, dihAtom2Type, dihAtom3Type, dihAtom4Type
	dihTypeArr = []
	sino_dih = 1

	dihTypeArr.append ({'atom1': 0, 'atom2': 0, 'atom3': 0, 'atom4': 0, 'dihType': 0})

	def findDihType (firstAtomType, secondAtomType, thirdAtomType, fourthAtomType, dihTypeArr):
		if (firstAtomType < fourthAtomType):
			ascFirstAtomType = firstAtomType
			ascSecondAtomType = secondAtomType
			ascThirdAtomType = thirdAtomType
			ascFourthAtomType = fourthAtomType
		else:
			ascFirstAtomType = fourthAtomType
			ascSecondAtomType = thirdAtomType
			ascThirdAtomType = secondAtomType
			ascFourthAtomType = firstAtomType

		assignedType = 0

		for types in dihTypeArr:
			if (ascFirstAtomType == types ['atom1'] and ascSecondAtomType == types ['atom2'] and ascThirdAtomType == types ['atom3'] and ascFourthAtomType == types ['atom4']):
				assignedType = types ['dihType']

		if (assignedType == 0):
			assignedType = len (dihTypeArr)
			dihTypeArr.append ({'atom1': ascFirstAtomType, 'atom2': ascSecondAtomType, 'atom3': ascThirdAtomType, 'atom4': ascFourthAtomType, 'dihType': len (dihTypeArr)})

		return assignedType, dihTypeArr

	for x in range (0, 360, 1):
		for y in range (x + 1, 360, 1):
			if (angleInfo [x]['angleAtom2'] == angleInfo [y]['angleAtom1'] and angleInfo [x]['angleAtom3'] == angleInfo [y]['angleAtom2']):
				firstAtomType = atomInfo [int (angleInfo [x]['angleAtom1']) - 1]['atomType']
				secondAtomType = atomInfo [int (angleInfo [x]['angleAtom2']) - 1]['atomType']
				thirdAtomType = atomInfo [int (angleInfo [x]['angleAtom3']) - 1]['atomType']
				fourthAtomType = atomInfo [int (angleInfo [y]['angleAtom3']) - 1]['atomType']
				dihType, dihTypeArr = findDihType (firstAtomType, secondAtomType, thirdAtomType, fourthAtomType, dihTypeArr)
				dihedralInfo.append ({'sino': sino_dih, 'dihType': dihType, 'dihAtom1': angleInfo [x]['angleAtom1'], 'dihAtom2': angleInfo [x]['angleAtom2'], 'dihAtom3': angleInfo [x]['angleAtom3'], 'dihAtom4': angleInfo [y]['angleAtom3'], 'dihAtom1Type': firstAtomType, 'dihAtom2Type': secondAtomType, 'dihAtom3Type': thirdAtomType, 'dihAtom4Type': fourthAtomType})
				sino_dih += 1

	return dihedralInfo, dihTypeArr

def createImpropers (atomInfo, bondInfo, angleInfo, dihedralInfo, improperInfo):
	improperTypeArr = []
	improperTypeArr.append ({'atom1': 0, 'atom2': 0, 'atom3': 0, 'atom4': 0, 'improperType': 0})

	return improperInfo, improperTypeArr

def printDataFile (atomInfo, atomTypeArr, bondInfo, bondTypeArr, angleInfo, angleTypeArr, dihedralInfo, dihTypeArr, improperInfo, improperTypeArr):

	nAtoms = len (atomInfo)
	nBonds = len (bondInfo)
	nAngles = len (angleInfo)
	nDihedrals = len (dihedralInfo)
	nImpropers = len (improperInfo)

	nAtomTypes = len (atomTypeArr)
	nBondTypes = len (bondTypeArr) - 1
	nAngleTypes = len (angleTypeArr) - 1
	nDihedralTypes = len (dihTypeArr) - 1
	nImproperTypes = 0

	coords_x = []
	coords_y = []
	coords_z = []

	with open ("atomEntries.output", "w") as file:
		for atomLine in atomInfo:
			file.write ("\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format (atomLine ['sino'], atomLine ['molType'], atomLine ['atomType'], atomLine ['charge'], atomLine ['x'], atomLine ['y'], atomLine ['z']))
			coords_x.append (atomLine ['x'])
			coords_y.append (atomLine ['y'])
			coords_z.append (atomLine ['z'])

	coords_xlo = min (coords_x) - 0.5 * min (coords_x)
	coords_xhi = max (coords_x) + 0.5 * max (coords_x)
	coords_ylo = min (coords_y) - 0.5 * min (coords_y)
	coords_yhi = max (coords_y) + 0.5 * max (coords_y)
	coords_zlo = min (coords_z) - 0.5 * min (coords_z)
	coords_zhi = max (coords_z) + 0.5 * max (coords_z)

	with open ("bondEntries.output", "w") as file:
		for bondLine in bondInfo:
			file.write ("\t{}\t{}\t{}\t{}\n".format (bondLine ['sino'], bondLine ['bondType'], bondLine ['bondAtom1'], bondLine ['bondAtom2']))

	with open ("angleEntries.output", "w") as file:
		for angleLine in angleInfo:
			file.write ("\t{}\t{}\t{}\t{}\t{}\n".format (angleLine ['sino'], angleLine ['angleType'], angleLine ['angleAtom1'], angleLine ['angleAtom2'], angleLine ['angleAtom3']))

	with open ("dihedralEntries.output", "w") as file:
		for dihedralLine in dihedralInfo:
			file.write ("\t{}\t{}\t{}\t{}\t{}\t{}\n".format (dihedralLine ['sino'], dihedralLine ['dihType'], dihedralLine ['dihAtom1'], dihedralLine ['dihAtom2'], dihedralLine ['dihAtom3'], dihedralLine ['dihAtom4']))

	with open ("output.data", "w") as file:
		file.write ("Created by you v1.8.1 on today, this month, this year, current time.\n\n\t{}\tatoms\n\t{}\tbonds\n\t{}\tangles\n\t{}\tdihedrals\n\t{}\timpropers\n\n\t{} atom types\n\t{} bond types\n\t{} angle types\n\t{} dihedral types\n\t{} improper types\n\n\t{}\t{}\txlo xhi\n\t{}\t{}\tylo yhi\n\t{}\t{}\tzlo zhi\n\nMasses\n\n\t1\t13.0907\t#CG311 CH\n\t2\t14.1707\t#CG321 CH2\n\t3\t15.2507\t#CG331 CH3\n\nAtoms\n\n".format (nAtoms, nBonds, nAngles, nDihedrals, nImpropers, nAtomTypes, nBondTypes, nAngleTypes, nDihedralTypes, nImproperTypes, coords_xlo, coords_xhi, coords_ylo, coords_yhi, coords_zlo, coords_zhi))

		for atomLine in atomInfo:
			file.write ("\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format (atomLine ['sino'], atomLine ['molType'], atomLine ['atomType'], atomLine ['charge'], atomLine ['x'], atomLine ['y'], atomLine ['z']))

		file.write ("\nBonds\n\n")
		for bondLine in bondInfo:
			file.write ("\t{}\t{}\t{}\t{}\n".format (bondLine ['sino'], bondLine ['bondType'], bondLine ['bondAtom1'], bondLine ['bondAtom2']))

		file.write ("\nAngles\n\n")
		for angleLine in angleInfo:
			file.write ("\t{}\t{}\t{}\t{}\t{}\n".format (angleLine ['sino'], angleLine ['angleType'], angleLine ['angleAtom1'], angleLine ['angleAtom2'], angleLine ['angleAtom3']))

		file.write ("\nDihedrals\n\n")
		for dihedralLine in dihedralInfo:
			file.write ("\t{}\t{}\t{}\t{}\t{}\t{}\n".format (dihedralLine ['sino'], dihedralLine ['dihType'], dihedralLine ['dihAtom1'], dihedralLine ['dihAtom2'], dihedralLine ['dihAtom3'], dihedralLine ['dihAtom4']))

		# file.write ("\nImpropers\n\n")



def main():
	atomInfo = []
	bondInfo = []
	angleInfo = []
	dihedralInfo = []
	improperInfo = []

	atomInfo, atomTypeArr = readAtomInfo ("atomEntries.testing", atomInfo)
	bondInfo, bondTypeArr = createBonds (atomInfo, bondInfo)
	angleInfo, angleTypeArr = createAngles (atomInfo, angleInfo, bondInfo)
	dihedralInfo, dihTypeArr = createDihedrals (atomInfo, dihedralInfo, angleInfo, bondInfo)
	improperInfo, improperTypeArr = createImpropers (atomInfo, bondInfo, angleInfo, dihedralInfo, improperInfo)

	printDataFile (atomInfo, atomTypeArr, bondInfo, bondTypeArr, angleInfo, angleTypeArr, dihedralInfo, dihTypeArr, improperInfo, improperTypeArr)

	# for angle in angleInfo:
	# 	print (angle)
	# 	sleep (1)

	# for dihedralLine in dihedralInfo:
	# 	print (dihedralLine)
	# 	sleep (1)

if __name__ == '__main__':
	main()