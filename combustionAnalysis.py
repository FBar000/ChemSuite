from array import array
import copy
from math import floor
from getMoles import *
from BalChemEq.methods import *
import numpy as np
from decimal import Decimal, ROUND_DOWN



def moles(atom, term, t_mass): 
    """
    Get the moles of an atom in a compound of a given mass.

    Arguments:
        atom (string) The symbol of the atom to be counted.
        term (string) The formula for the compound to count from.
        t_mass (float) The mass of the compound to count from.
    Return:
        moles_of_atom (float) The moles of the atom present in the compound.
    """
    return t_mass / Decimal(getMolarMass(term)) * Decimal(count(atom, (term, 1)))

def check(raw_number): 
    """
    Returns true if float is within 10^-3 of an integer. False otherwise.

    Arguments: 
        raw_number (Decimal): A Decimal type
    """
    return abs(raw_number.quantize(Decimal('0.')) - raw_number) < Decimal('0.05')

def checkD(raw_array):
    """
    Returns true if items in array are within 10^-3 of an integer. False otherwise.

    Arguments: 
        raw_array (ndarray type=Decimal): A Decimal type.
    """
    for dec in raw_array:
        if not check(dec):
            return False
    return True

def checkZero(raw_number):
    """
    Returns true if float is within 10^-3 of 0. False otherwise.

    Arguments: 
        raw_number (Decimal): A Decimal type
    """
    return check(raw_number-Decimal('0'))


def combustionAnalysisFindUnknown(unknown_MolarMass, unknown_compound_mass, product_information):
    """
    Get the composition of an unknown reactant in a combustion reaction with excess O2

    Arguments:
        unknown_MolarMass (Decimal): The molar mass of the unknown compound.
        unknown_compound_mass (Decimal): The mass of the sample of the unknown compound.
        product_information (dict={string:Decimal}): A dict with the products and their masses.
    Return:
        solution (dict={string:int}): The composition of the unknown compound.
    """
    # Molar mass of unknown compound
    molar_mass = unknown_MolarMass
    # Sig figs assumed to unkown_compound_mass
    sf = Decimal('0.0001')
    uk_mantissa = Decimal('0.01')
    # A dictionary of the present atoms
    atoms = findAtoms("".join(product_information.keys()))
    atoms.remove("O")
    # A dict with {atom: moles}
    info = {}
    for atm in atoms:
        tmp = 0
        for term, mass in product_information.items():
            tmp += moles(atm, term, mass).quantize(sf)
        info[atm] = tmp
    # Find Moles Oxygen
    tmp = 0
    for atom, mole in info.items():
        tmp += (mole * Decimal(getMolarMass(atom))).quantize(sf) 
    tmp2 = (unknown_compound_mass - tmp).quantize(uk_mantissa) / Decimal(getMolarMass("O"))
    # Add "O" to list only if nonzero
    if not checkZero(tmp2):
        info["O"] = tmp2
        atoms.append("O")
    # Add information into a vector
    info_vector = np.fromiter(list(info.values()), dtype=Decimal)
    info_vector = info_vector / np.amin(info_vector)
    # Make the vector all integers
    adder = copy.copy(info_vector)
    while not checkD(info_vector):
        info_vector += adder
    # Round to integer
    aux = Decimal('0')
    for i in range(info_vector.size):
        info_vector[i] = info_vector[i].quantize(aux) 
    # Extract info from vector
    solution = {}
    simple_molarmass = 0
    for i in range(len(atoms)):
        atom = atoms[i]
        mole = int(info_vector[i])
        solution[atom] = mole
        simple_molarmass += Decimal(getMolarMass(atom)) * mole
    # Scale to have accurate molar mass
    factor = round(molar_mass / simple_molarmass)
    for i in solution:
        solution[i] *= factor
    # Done
    return solution

if __name__ == '__main__':

    X_mass = Decimal('5.50')
    X_mm = Decimal('78')
    product_information = {
        "CO2": Decimal('18.62'),
        "H2O": Decimal('3.81')
    }
    print(combustionAnalysisFindUnknown(X_mm, X_mass,product_information))