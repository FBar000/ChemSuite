"""
This module contains scripts for calculating pH/OH values, acid/base constants, etc.
"""

import getMoles
from decimal import Decimal

def calcBaseTitrantMass(pH, Ka, V, M, X):
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

def calcAcidTitrantMass(pH, Kb, V, M, X):
    """
    Returns the mass of the acid titrant needed to achieve a given pH from volume of an acid solution with a given molarity.

    This assumes that adding a substance does not change the volume of the solution.

    Args: 
        pH  (Decimal):  The target pH
        Kb  (Decimal):  The base's equilibrium constant 
        M   (Decimal):  The base's molarity. (Molal)
        V   (Decimal):  The base's volume.   (Liters)
        X   (string):   The acid's chemical formula.

    Return:
        m (Decimal):    The mass of acid required.
    """
    Ka = Decimal(10**(-14)) / Kb
    molar_mass = Decimal(getMoles.getMolarMass(X))
    cH = 10**(pH)
    return molar_mass * M * V / (cH * Ka) 




if __name__ == '__main__':

    givens = {
        "pH": 3.66,
        "K": 4.5*10**(-4),
        "V": 0.3,
        "M": 0.7,
        "X": "KNO2"
    }

    for i in givens:
        if type(givens[i]) == float:
            givens[i] = Decimal(givens[i])

    print(calcAcidTitrantMass(
                        givens["pH"],
                        givens["K"],
                        givens["V"],
                        givens["M"],
                        givens["X"]))
    
    print(calcBaseTitrantMass(
                        givens["pH"],
                        givens["K"],
                        givens["V"],
                        givens["M"],
                        givens["X"]))