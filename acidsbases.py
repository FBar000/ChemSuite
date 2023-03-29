"""
This module contains scripts for calculating pH/OH values, acid/base constants, etc.
"""

import math
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

def titrateWASB(Va, Ma, Vb, Mb, pKa):
    """
    Calculate the pH of a weak acid after titration with a strong base.

    Args (Decimals):
        Va, Vb: Volumes of acid and base    (Liters)
        Ma, Mb: Masses of acid and base     (Molal)
        pKa: Power of acid constant  

    Return: 
        pH (Decimal)
    """
    return pKa + Decimal(math.log((Ma * Va) / (Mb * Vb), 10))

def EQPpH(V, M, M2, pKa):
    """
    Calculate the pH of a solution at it's equivalence point.

    Args (Decimals):
        V: Volume of analyte    (Liters)
        M: Molarity of analyte  (Molal)
        M2: Molarity of titrant  (Molal)
        pKa: Analyte constant

    Return:
        pH (Decimal)
    """
    return pKa + Decimal(math.log(M / M2, 10))

if __name__ == '__main__':

    givens = {
        "V": 0.22,
        "M": 0.5901,
        "M2": 135,
        "pKa": 4.82
    }

    for i in givens:
        givens[i]=Decimal(givens[i])
        
    givens['pKa'] = Decimal(10^(-14)) / givens['pKa']

    print(EQPpH(**givens))