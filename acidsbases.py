"""
This module contains scripts for calculating pH/OH values, acid/base constants, etc.
"""

import getMoles
from decimal import Decimal

def calcTitrantMass(pH, Ka, V, M, X):
    """
    Returns the mass of the base titrant needed to achieve a given pH from volume of an acid solution with a given molarity.

    This assumes that adding a substance does not change the volume of the solution.

    Args: 
        pH  (Decimal):  The target pH
        Ka  (Decimal):  The acid's equilibrium constant 
        V   (Decimal):  The acid's volume.   (Liters)
        M   (Decimal):  The acid's molarity. (Molal)
        X   (string):   The base's chemical formula.

    Return:
        m (Decimal):    The mass of base required.
    """
    cH = 10**(pH)
    molar_mass = Decimal(getMoles.getMolarMass(X))
    return cH * M * V * molar_mass * Ka

if __name__ == '__main__':

    givens = {
        "pH": 5.31,
        "Ka": 1.8*10**(-5),
        "V": 0.125,
        "M": 1.3,
        "X": "NaCH3CO2"
    }

    for i in givens:
        if type(givens[i]) == float:
            givens[i] = Decimal(givens[i])

    print(calcTitrantMass(
                        givens["pH"],
                        givens["Ka"],
                        givens["V"],
                        givens["M"],
                        givens["X"]))