"""
This module contains scripts for calculating pH/OH values, acid/base constants, etc.
"""

import math
import getMoles
from decimal import Decimal
from BalChemEq import BCE
import numpy as np

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

    HA + B -> HB + A

    Args (Decimals):
        Va, Vb: Volumes of acid and base    (Liters)
        Ma, Mb: Masses of acid and base     (Molal)
        pKa: Power of acid constant  

    Return: 
        pH (Decimal)
    """
    acid_moles = Va * Ma
    base_moles = Vb * Mb
    if acid_moles >= base_moles:
        x = acid_moles / base_moles
        log = Decimal(math.log(1 / (x - 1), 10)) # Algebra simplifies to this
        pH = pKa + log     # HH Equation
    else: 
        givens = {"A": acid_moles,
                  "H": Decimal(0),
                  "HA": Decimal(0)}
        equation = "A + H : HA"
        sols = ICEDiagram(givens, equation, pKa)
        pH = -Decimal(math.log(sols["H"], 10))
    return pH


def EQPpH(M, M2, pKa):
    """
    Calculate the pH of a solution at it's equivalence point. 
    Assumes acid/base are both monoprotic.

    Args (Decimals):
        V: Volume of analyte    (Liters)
        M: Molarity of analyte  (Molal)
        M2: Molarity of titrant  (Molal)
        pKa: Analyte constant

    Return:
        pH (Decimal)
    """
    print(M)
    print(M2)
    print(pKa)
    x = Decimal(math.log(M**(-1) + M2**(-1), 10))
    return Decimal(0.5) * (pKa + x)

def ICEDiagram(molarities, equation, K):
    """
    Returns the final row of an ICE diagram.
    Assumes two reactants, =<two products

    Args:
        molarities (dict[Decimal]): a dictionary of the initial molarities of the species. Absentees assumed zero.
        equation (string): an unbalanced chemical equation for the reaction
        K (Decimal): The equilibrium constant
    Return:
        final_molarities (dict[Decimal]): the equilibrium concentrations
    """
    coefficients = BCE.solutionCoefficients(equation)
    # decimalize everyting
    decimalize(molarities)
    decimalize(coefficients)
    coefficients["null"] = 0
    # assign indices to equation terms
    term_list = equation.split(" ")
    while len(term_list) < 4:
        term_list.append("null")
    # easy reference to coefficients
    a, b, c, d = coefficients[term_list[0]], coefficients[term_list[1]], coefficients[term_list[2]], coefficients[term_list[3]]
    # easy reference to molarities
    A, B, C, D = molarities[term_list[0]], molarities[term_list[1]], molarities[term_list[2]], molarities[term_list[3]]
    # first coef. in quadratic
    X1 = determinantStep(a, c, b*K, d)
    # second coef.
    X2 = determinantStep(a, b, -A, -B) * K - determinantStep(c, d, -C, -D)
    X2*= Decimal(-1)
    # third coef.
    X3 = determinantStep(A, C, B*K, D)
    # Change
    change = quadraticPosRoot(X1, X2, X3)
    # Adjust signs of coefficients
    coefficients[term_list[2]], coefficients[term_list[3]] = -coefficients[term_list[2]], -coefficients[term_list[3]]
    # Pack solutions
    final_molarities={}
    for i in molarities:
        final_molarities[i] = molarities[i] - change * coefficients[i]
    return final_molarities

def determinantStep(a, b, c, d):
    "For use in quadratic solver only"
    return a*d - b*c

def decimalize(dict):
    """
    Convert dict values into Decimal types.

    Return None.
    """
    for i in dict:
        dict[i] = Decimal(dict[i])

def quadraticPosRoot(a, b, c):
    return (-b + math.sqrt(b**2 - 4*a*c)) / 2*a



if __name__ == '__main__':

    givens = {
        "M1": 0.65,
        "V1": 119.5,
        "M2": 0.88,
        "V2": 100.6,
        "pK": 4.89 
    }
    decimalize(givens)
    print(titrateWASB(**givens))