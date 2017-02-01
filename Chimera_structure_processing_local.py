# Call this script from cygwin as: ./Chimera\ 1.10.1/bin/chimera.exe --start "Attribute Calculator" --start "Reply Log" --script "./chimeraStructureProcessForVOiD_6-9-2016.py"

import os
import chimera
import re
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
from chimera import replyobj # for emitting status messages
from chimera import selection
from chimera import dialogs

# Gather the names of .pdb files in the specified folder
directoryName = "./Protein_Structures_For_Maisonneuve_2009"

# To get all pdb files in the specified directory:
#file_names = [fn for fn in os.listdir(directoryName) if fn.endswith(".pdb")]

# To specify individual pdb files:
file_names = ["./Ecoli_I-TASSER/E14278.pdb"]


### EXAMPLE OF DEFINING A DICTIONARY
# Dictionaries:
# Amino acid 3-to-1 letter code dictionary.
aaDict = dict({'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'})
# Residue number to 1 letter code dictionary.
resDict = dict({})

# Loop through the files, opening, processing, and closing each in turn.
for fn in file_names:

        # Open file in Chimera.
        replyobj.status("Processing " + fn) # show what file we're working on
#        rc("open " + directoryName + "/" + fn)
        rc("open " + fn)

        # Parse input file name:
#        fn = re.sub(r directoryName, '', fn)
        fn = re.sub(r'./', '', fn)
        fn = re.sub(r'.pdb', '', fn)
        fn = re.sub(r'allStruct', '', fn)

        # Output filenames:
        areaSASFile = './chimeraVOiD_areaSAS_' + fn + '.txt'
        f_areaSAS = open(str(areaSASFile), 'w')
        depthFile = './chimeraVOiD_depth_' + fn + '.txt'
        f_depth = open(str(depthFile), 'w')
        distMatrixFile = './chimeraVOiD_distMatrix_' + fn + '.txt'
        f_distMatrix = open(distMatrixFile, 'w')
        rlbcoorFile = './' + fn + '.rlbcoor'
        f_rlbcoor = open(rlbcoorFile, 'w')
        ssFile = './chimeraVOiD_SS_' + fn + '.txt'
        f_SS = open(str(ssFile), 'w')
        hydroFile = './chimeraVOiD_hydro_' + fn + '.txt'
        f_hydro = open(str(hydroFile), 'w')

        # Clean out all non-peptide atoms.
        rc("select protein")
        rc("select invert")
        rc("delete selected")

        # Delete all alternative locations for residues.
        rc("select #0:@.B")
        rc("delete selected")

        # Initialize areaSAS for each atom of protein, otherwise throws error for any atom without areaSAS attribute defined.
        rc("select protein")
        atoms = selection.currentAtoms()
        for a in atoms:
                a.areaSAS = 0

        # Generate surface to protein.
        rc("select protein")
        rc("split")
        rc("surface allComponents false")
#        rc("surface probeRadius 2.8")
        rc("~select")

        # Select appropriate RKPT(W) atoms that are carbonylatable. R,K,P come from Requena 2003. T infered based on conversion of threonine to 2-amino-3-ketobutyric acid as noted in Levine 2001 and elsewhere. W comes from Taylor 2003.
        # 'R_CD', 'K_CE', 'P_CD', 'T_CB', 'W_CD1'
#        rc("select :arg@cd|:lys@ce|:pro@cd|:thr@cb|:trp@cd1")
        rc("select :arg@cd|:lys@ce|:pro@cd|:thr@cb")
        atoms = selection.currentAtoms()



        # COMPUTE SOLVENT ACCESSIBLE SURFACE AREA FOR EACH ATOM AND SS ASSIGNMENT AND HYDROPHOBICITY OF RESIDUE INCLUDING ATOM.
        atomsList = []
        areaSASDict = dict({})
        helixDict = dict({})
        sheetDict = dict({})
        hydroDict = dict({})

        for a in atoms:
                areaSAS = a.areaSAS
                residue = a.residue
                type = residue.type
                resHelix = residue.isHelix
                resSheet = residue.isSheet
                resHydro = residue.kdHydrophobicity
                resNum = str(residue)
                resNum = re.sub(r"#.+ .+ ", "", resNum)

### EXAMPLE OF LOOKING UP A KEY AND RETURNING A VALUE FROM A DICTIONARY
                resDict[resNum] = str(aaDict[type])
                name = a.name
                altLoc = a.altLoc
                if altLoc == "":
                        string = str(residue) + " " + str(aaDict[type]) + " " + str(name) + altLoc
                else:
                        string = str(residue) + " " + str(aaDict[type]) + " " + str(name) + '.' + altLoc
                string = re.sub(r" ARG ", ":", string)
                string = re.sub(r" LYS ", ":", string)
                string = re.sub(r" PRO ", ":", string)
                string = re.sub(r" THR ", ":", string)
                string = re.sub(r" TRP ", ":", string)
                areaSASDict[string] = areaSAS
                helixDict[string] = resHelix
                sheetDict[string] = resSheet
                hydroDict[string] = resHydro
                atomsList.append(string)
        atomsList = sorted(atomsList)



        # COMPUTE DEPTH OF EACH ATOM (i.e. minimum distance to surface)
        depthDict = dict({})
        minDepth = 1000000
        rc("select #0")
        rc("~select #0:@")
#        rc("measure distance :arg@cd|:lys@ce|:pro@cd|:thr@cb|:trp@cd1 selection multiple true show true")
        rc("measure distance :arg@cd|:lys@ce|:pro@cd|:thr@cb selection multiple true show true")
        r = dialogs.find('reply')
        text = r.text.get('1.0', 'end')
        f_temp = open('./chimeraVOiD_temp.txt', 'w')
        f_temp.write(text)
        f_temp.close()

### EXAMPLE OF OPENING A TEXT FILE, ITERATING THROUGH THE LINES, AND SOME REGULAR EXPRESSION OPERATIONS.
        with open('./chimeraVOiD_temp.txt','r') as f_temp:
                for line in f_temp:
                        line = line.rstrip()
                        if line[:21] == 'minimum distance from':
                                line = line.replace('minimum distance from ', '')
                                line = re.sub(r" to #0:\? = ", "\t", line)
                                resNum = line
                                resNum = re.sub(r"#.+:", "", resNum)
                                resNum = re.sub(r"@.+", "", resNum)
                                line = re.sub(r'@', (" " + resDict[resNum] + " "), line)
                                string = re.sub(r"\t.+", "", line)
                                depth = re.sub(r".+\t", "", line)
                                depthDict[string] = depth
                                if float(depth) < minDepth:
                                        minDepth = float(depth)
        os.remove('chimeraVOiD_temp.txt')



        # COMPUTE INTERATOMIC DISTANCES. (Notes that the order of atoms selected is not deterministic and so varies across runs.)
        distMatrix = [[0 for x in range(len(atoms))] for x in range(len(atoms))]
        atomsIndex = dict({})
        resTypeChainIndex = dict({})
        xformCoordIndex = dict({})
        counter = 0
        for a1 in atoms:
                residue = a1.residue
                type = residue.type
                name = a1.name
                altLoc = a1.altLoc
                if altLoc == "":
                        string = str(residue) + " " + str(aaDict[type]) + " " + str(name) + altLoc
                else:
                        string = str(residue) + " " + str(aaDict[type]) + " " + str(name) + '.' + altLoc
                string = re.sub(r" ARG ", ":", string)
                string = re.sub(r" LYS ", ":", string)
                string = re.sub(r" PRO ", ":", string)
                string = re.sub(r" THR ", ":", string)
                string = re.sub(r" TRP ", ":", string)

                atomsIndex[string] = counter
                counter = counter + 1

        counter1 = 0
        for a1 in atoms:
                residue = a1.residue
                type = residue.type
                name = a1.name
                altLoc = a1.altLoc
                if altLoc == "":
                        string = str(residue) + " " + str(aaDict[type]) + " " + str(name) + altLoc
                else:
                        string = str(residue) + " " + str(aaDict[type]) + " " + str(name) + '.' + altLoc
                string = re.sub(r" ARG ", ":", string)
                string = re.sub(r" LYS ", ":", string)
                string = re.sub(r" PRO ", ":", string)
                string = re.sub(r" THR ", ":", string)
                string = re.sub(r" TRP ", ":", string)


                # For rlbcoor data.
                resTypeChain = str(residue) + " A"
                resTypeChain = re.sub(r"#0 ", "",resTypeChain)
                resTypeChainIndex[string] =  str(resTypeChain)
                xformCoordIndex[string] = str(a1.xformCoord())
                counter2 = 0
                for a2 in atoms:
                        distance = chimera.distance(a1.xformCoord(), a2.xformCoord())
                        distMatrix[counter2][counter1] = str(distance)
                        counter2 = counter2 + 1
                counter1 = counter1 + 1


### EXAMPLES OF WRITING VARIABLES TO OUTPUT FILES.
        # PRINT OUT ALL RESULTS.
        f_distMatrix.write("\t")
        for a1 in atomsList:
                atom1 = str(a1)
                f_distMatrix.write(atom1)
                f_distMatrix.write("\t")
        f_distMatrix.write("\n")

        for a1 in atomsList:
                atom1 = str(a1)
                areaSAS = areaSASDict[atom1]
                depth = depthDict[atom1]
                helix = helixDict[atom1]
                sheet = sheetDict[atom1]
                hydro = hydroDict[atom1]

                f_areaSAS.write(atom1)
                f_areaSAS.write("\t")
                f_areaSAS.write(str(areaSAS))
                f_areaSAS.write("\n")

                f_SS.write(atom1)
                f_SS.write("\t")
                f_SS.write(str(helix))
                f_SS.write("\t")
                f_SS.write(str(sheet))
                f_SS.write("\n")

                f_hydro.write(atom1)
                f_hydro.write("\t")
                f_hydro.write(str(hydro))
                f_hydro.write("\n")

                f_depth.write(atom1)
                f_depth.write("\t")
                f_depth.write(str(depth))
                f_depth.write("\n")

                f_distMatrix.write(atom1)
                f_distMatrix.write("\t")

                for a2 in atomsList:
                        atom2 = str(a2)
                        a1Index = atomsIndex[atom1]
                        a2Index = atomsIndex[atom2]
                        distance = distMatrix[a2Index][a1Index]
                        f_distMatrix.write(str(distance))
                        f_distMatrix.write("\t")
                f_distMatrix.write("\n")

        f_rlbcoor.write("\n")
        f_rlbcoor.write("\n")
        f_rlbcoor.write(" loop_\n")
        f_rlbcoor.write("   _relibase_surface_pseudocentre_id\n")
        f_rlbcoor.write("   _relibase_surface_pseudocentre_name\n")
        f_rlbcoor.write("   _relibase_surface_pseudocentre_number\n")
        f_rlbcoor.write("   _relibase_surface_pseudocentre_residue_type\n")
        f_rlbcoor.write("   _relibase_surface_pseudocentre_residue_number\n")
        f_rlbcoor.write("   _relibase_surface_pseudocentre_chainid\n")
        f_rlbcoor.write("   _relibase_surface_pseudocentre_backbone\n")
        f_rlbcoor.write("   _relibase_surface_pseudocentre_x\n")
        f_rlbcoor.write("   _relibase_surface_pseudocentre_y\n")
        f_rlbcoor.write("   _relibase_surface_pseudocentre_z\n")
        counter2 = 0
        for a1 in atomsList:
                atom1 = str(a1)
                depth = depthDict[atom1]
                label = ''
                num = 0
                backbone = 0
                coords = xformCoordIndex[atom1]

                xyzCoords = coords.split(' ',2)

                if xyzCoords[0].find('.') > -1:
                        xDecimal = xyzCoords[0].split('.',1)[1]
                        xAddZeros = 4 - len(xDecimal)
                        for i in range(xAddZeros):
                                xyzCoords[0] = xyzCoords[0] + "0"
                else:
                        xyzCoords[0] = xyzCoords[0] + ".0000"
                if xyzCoords[1].find('.') > -1:
                        yDecimal = xyzCoords[1].split('.',1)[1]
                        yAddZeros = 4 - len(yDecimal)
                        for i in range(yAddZeros):
                                xyzCoords[1] = xyzCoords[1] + "0"
                else:
                        xyzCoords[1] = xyzCoords[1] + ".0000"
                if xyzCoords[2].find('.') > -1:
                        zDecimal = xyzCoords[2].split('.',1)[1]
                        zAddZeros = 4 - len(zDecimal)
                        for i in range(zAddZeros):
                                xyzCoords[2] = xyzCoords[2] + "0"
                else:
                        xyzCoords[2] = xyzCoords[2] + ".0000"

                if float(depth) == float(minDepth):
                        label = 'Donor'
                        num = 0
                else:
                        label = 'Acceptor'
                        num = 1
                string = str(num) + " " + str(label) + " " + str(counter2) + " " + str(resTypeChainIndex[str(atom1)]) + " " + str(backbone) + " " + str(xyzCoords[0]) + " " + str(xyzCoords[1]) + " " + str(xyzCoords[2]) + "\n"
                f_rlbcoor.write(string)

                counter2 = counter2 + 1


        f_distMatrix.close()
        f_areaSAS.close()
        f_temp.close()
        f_depth.close()
        f_SS.close()
        f_hydro.close()

        rc("close all")
rc("stop now")
