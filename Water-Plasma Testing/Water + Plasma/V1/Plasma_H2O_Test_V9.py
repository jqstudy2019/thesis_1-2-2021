########################################################################################################################
#   Water_Place_Test_v6.py ----- A code to build the .data file for a LAMMPS Simulation
#
#       8/4/2020
#
#       General Description -   Work On adding elements to the simulation ( Water Plasma CytoSkeleton Etc)
#
#                               ALL Non Essential Code Removed and/or commented out
#
########################################################################################################################

##      Import Extra Modules / Tools        ##

import numpy as np  # Numpy is used for Array Creation and Manipulation Mainly
import random as random
from copy import deepcopy  # DeepCopy Is used to Copy Arrays


########################################################################################################################
##      Self Defined Functions      ##

def find_atom_id(x_find, y_find, z_find):
    for iV3, posV3 in enumerate(molecule_1):
        if (molecule_1[posV3].get('X') == x_find) \
                and (molecule_1[posV3].get('Y') == y_find) \
                and (molecule_1[posV3].get('Z') == z_find):
            return molecule_1[posV3].get('Atom_ID')


def find_atom_type(atom_id):
    for iV5, posV5 in enumerate(ATOMS):
        if ATOMS[posV5].get('Atom_ID') == atom_id:
            return ATOMS[posV5].get('Atom_Type')


def packing_roll(coeff):
    temp_roll = random.randint(1, 100)
    if temp_roll in range(0, coeff):
        return 'Place Atom'
    else:
        return 'Leave Empty'


########################################################################################################################

###      Set Variables Used For Later       ###


##     Output File Options      ##

filename = input("Enter a name for the Output File (remember .data at the end!) :")  # String Name of File to be Created

##      Simulation Box Settings     ##

sim_box_side_length = 100  # Each Side of the Simulation Box Will have sides of these lengths

SIM_SPACE = np.zeros((sim_box_side_length, sim_box_side_length, sim_box_side_length))
# Creates an Array Called SIM_SPACE (is a 3-D array (x,y,z of each simulation space position)

x0 = int(np.floor(SIM_SPACE.shape[2] / 2))  # Finds the Center of the Sim Space ( x0, y0, z0 )
y0 = int(np.floor(SIM_SPACE.shape[1] / 2))
z0 = int(np.floor(SIM_SPACE.shape[0] / 2))

##      Molecule Geometry Creation Settings         ##

membrane_atom_type = 1  # LAMMPS Track molecules by their ID number ; this is the number of molecules created

membrane_molecule_type = 1

molecule_1_thickness = 15  # The thickness is an extra range past the radius to include in the molecules

molecule_1_radius = 10  # Radius of each molecule

num_atoms_molecule_1 = 0  # Placeholder for total number of atoms created as part of the molecule

molecule_1_bond_length = 2  # Functions Similar to molecule radius but on the bond scale

molecule_1 = {}  # initializes the dictionary for molecule 1

molecule_1_bonds = {}  # dictionary for bonds

num_molecule_1_bonds = 0

molecule_1_bond_type = 1


##      Water Variables     ##

num_atoms_water = 0

water_atom_type = 2

water_molecule_type = 2

water_offset = 50

water_packing_coeff = 5


##      Plasma Variables

num_atoms_plasma = 0

plasma_atom_type = 3

plasma_molecule_type = 3

plasma_offset = 25

plasma_packing_coeff = 50


##       LAMMPS Variables       ##

num_atoms_total = 0

total_atom_types = 3

total_bond_types = 1

##         Program Performance Tracking         ##

positions_checked = 0  # Keeps Count of How many points are iterated through ; easy way to check it iterating correctly

##        EXPERIMENTAL       ##

ATOMS = {}

########################################################################################################################

###     Creates 3-D Set of Points      ###

# Iterates Each Position of the sim (x,y,z

for x in range(0, sim_box_side_length):
    for y in range(0, sim_box_side_length):
        for z in range(0, sim_box_side_length):
            deb = (molecule_1_radius ** 2) - ((x0 - x) ** 2) - ((y0 - y) ** 2) - ((z0 - z) ** 2)
            positions_checked = positions_checked + 1
            if (deb >= 0) and (deb <= molecule_1_thickness):
                num_atoms_molecule_1 = num_atoms_molecule_1 + 1
                SIM_SPACE[x, y, z] = membrane_atom_type
                num_atoms_total = num_atoms_total + 1
                ATOMS[num_atoms_total] = {'Atom_ID': num_atoms_total,
                                          'Molecule_ID': membrane_molecule_type,
                                          'Atom_Type': membrane_atom_type,
                                          'X': x,
                                          'Y': y,
                                          'Z': z,
                                          'Positions String': '{} {} {}'.format(x, y, z)}

                for x1 in range(x - molecule_1_bond_length, x + molecule_1_bond_length):
                    for y1 in range(y - molecule_1_bond_length, y + molecule_1_bond_length):
                        for z1 in range(z - molecule_1_bond_length, z + molecule_1_bond_length):
                            Atom_ID_2 = find_atom_id(x1, y1, z1)
                            if Atom_ID_2 is not None:
                                if Atom_ID_2 != num_atoms_molecule_1:
                                    Atom_Type_2 = find_atom_type(Atom_ID_2)
                                    if Atom_Type_2 == membrane_atom_type:
                                        num_molecule_1_bonds = num_molecule_1_bonds + 1
                                        molecule_1_bonds[num_molecule_1_bonds] = {'Bond_ID': num_molecule_1_bonds,
                                                                                  'Bond_Type': molecule_1_bond_type,
                                                                                  'atom1': num_atoms_molecule_1,
                                                                                  'atom2': Atom_ID_2}

########################################################################################################################

###         Plasma Placement Section         ###

for x2 in range(0, sim_box_side_length):
    for y2 in range(0, sim_box_side_length):
        for z2 in range(0, sim_box_side_length):
            if SIM_SPACE[x2, y2, z2] == 0:
                debV2 = (molecule_1_radius ** 2) - ((x0 - x2) ** 2) - ((y0 - y2) ** 2) - ((z0 - z2) ** 2)
                if debV2 >= plasma_offset:
                    roll = packing_roll(plasma_packing_coeff)
                    if roll == 'Place Atom':
                        num_atoms_plasma = num_atoms_plasma + 1
                        SIM_SPACE[x2, y2, z2] = plasma_atom_type
                        num_atoms_total = num_atoms_total + 1
                        ATOMS[num_atoms_total] = {'Atom_ID': num_atoms_total,
                                                  'Molecule_ID': plasma_molecule_type,
                                                  'Atom_Type': plasma_atom_type,
                                                  'X': x2,
                                                  'Y': y2,
                                                  'Z': z2,
                                                  'Positions String': '{} {} {}'.format(x2, y2, z2)}


########################################################################################################################

###         Water Placement Section         ###

for x3 in range(0, sim_box_side_length):
    for y3 in range(0, sim_box_side_length):
        for z3 in range(0, sim_box_side_length):
            if SIM_SPACE[x3, y3, z3] == 0:
                debV2 = (molecule_1_radius ** 2) - ((x0 - x3) ** 2) - ((y0 - y3) ** 2) - ((z0 - z3) ** 2)
                if debV2 <= water_offset:
                    roll = packing_roll(water_packing_coeff)
                    if roll == 'Place Atom':
                        num_atoms_water = num_atoms_water + 1
                        SIM_SPACE[x3, y3, z3] = water_atom_type
                        num_atoms_total = num_atoms_total + 1
                        ATOMS[num_atoms_total] = {'Atom_ID': num_atoms_total,
                                                  'Molecule_ID': water_molecule_type,
                                                  'Atom_Type': water_atom_type,
                                                  'X': x3,
                                                  'Y': y3,
                                                  'Z': z3,
                                                  'Positions String': '{} {} {}'.format(x3, y3, z3)}

########################################################################################################################

###     Prints Sim Info To Make Sure Its Running Correct        ###

print("Number of Points in Simulation Space Iterated thorough: {} \n".format(positions_checked))

print("Number of Water Atoms: {} ".format(num_atoms_water))

print("Number of Plasma Atoms: {} ".format(num_atoms_plasma))

print("Number of Plasma Atoms: {} ".format(num_atoms_molecule_1))

########################################################################################################################

###     Write LAMMPS Data File      ###

with open(filename, 'w+') as fdata:  # opens a text file named a for the 'filename' variable
    fdata.write('{}\n\n'.format(filename))  # First line is a comment line

    ##     Header of Data File     ##

    #   Atoms Header #

    fdata.write('{} atoms\n'.format(num_atoms_total))  # Specify number of atoms
    fdata.write('{} atom types\n'.format(total_atom_types))

    # Bonds Header  #
    fdata.write('{} bonds\n'.format(num_molecule_1_bonds))  # Specify number of atoms
    fdata.write('{} bond types\n'.format(total_bond_types))

    #   specify box dimensions      #

    fdata.write('{} {} xlo xhi\n'.format(0.0, sim_box_side_length))  # Writes X Position
    fdata.write('{} {} ylo yhi\n'.format(0.0, sim_box_side_length))  # Writes Y Position
    fdata.write('{} {} zlo zhi\n'.format(0.0, sim_box_side_length))  # Writes Z Position

    fdata.write('\n')  # Skips the Next line for Formatting

    ##      Data File Body      ##

    #       Atoms section       #
    fdata.write('Atoms\n\n')
    for i, pos in enumerate(ATOMS):
        fdata.write('{} {} {} {} \n'.format(ATOMS[pos].get('Atom_ID'),
                                            ATOMS[pos].get('Molecule_ID'),
                                            ATOMS[pos].get('Atom_Type'),
                                            ATOMS[pos].get('Positions String')))
    fdata.write('\n')

    #       Bonds Section       #
    fdata.write('Bonds\n\n')
    for iV2, posV2 in enumerate(molecule_1_bonds):
        fdata.write('{} {} {} {} \n'.format(molecule_1_bonds[posV2].get('Bond_ID'),
                                            molecule_1_bonds[posV2].get('Bond_Type'),
                                            molecule_1_bonds[posV2].get('atom1'),
                                            molecule_1_bonds[posV2].get('atom2')))

########################################################################################################################
###         Finished!!!     ###

print(' Data File Created Successfully!!! ; File name => {} '.format(filename))

########################################################################################################################
