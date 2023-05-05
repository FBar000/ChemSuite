"""
Scripts to solve ALEKS Ch16 problems
"""


import math
import re
import getMoles
from decimal import Decimal
from BalChemEq import BCE
import numpy as np
import pyromat as pm


def equationEntropy(equation):
    """
    Find the entropy of a chemical equation
    
    Arg:
        equation (str)

    Return: 
        entropy (float)
    """
    reactants, products = BCE.processEquation(equation)
    p_entropy = sum([molecule[1] * getEntropy(molecule[0]) for molecule in products[0]])
    r_entropy = sum([molecule[1] * getEntropy(molecule[0]) for molecule in reactants[0]])
    entropy = p_entropy - r_entropy
    return entropy


def getEntropy(molecule):
    return
