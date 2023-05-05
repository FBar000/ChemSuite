"""
Scripts to solve ALEKS Ch16 problems
"""



from BalChemEq import BCE
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
    rxn_entropy = p_entropy - r_entropy
    return rxn_entropy


def getEntropy(molecule, temperature=273.15, pressure=1):
    """
    Get entropy (kJ/mol) for a molecule at a given temperature (K) and pressure (atm)
    Arg: 
        molecule (str)

    Return:
        entropy (float)
    """
    name = '.'.join('ig', molecule)
    gas = pm.get(name)
    entropy = gas.h(temperature, pressure)
    return entropy
