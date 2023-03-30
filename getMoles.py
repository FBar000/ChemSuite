import mendeleev
# from methods import *
from BalChemEq.methods import *
from BalChemEq import BCE

def getMolarMass(base_substance, target_substance=''):
    """
    Get grams of target_substance per mole of substance_formula. 
    
    target_substance defaults to substance_formula,
    
    Argument:
        base_substance (string): the formual of the chemical substance
        target_substance (string): a subset of the substance formula
    Return:
        const (float): grams of target_substance per mole of substance_foruma
    """
    # Default to molar mass
    if target_substance == '':
        target_substance = base_substance
    const = 0
    List =  findAtoms(target_substance)
    for atom in List:
        amount_of_atom = count(atom, (base_substance, 1))
        atomic_weight = mendeleev.element(atom).atomic_weight
        const += amount_of_atom * atomic_weight
    return const


def massToMoles(substance_formula, mass):
    """
    Returns the moles of the substance.

    Arguments:
        substance_formula (string): A chemical formula (molecular formula)
        mass (float): The mass (in grams) of the substance
    Returns:
        moles (float): The moles of the element 
    """
    const =  getMolarMass(substance_formula)
    moles = mass / const
    return moles


def molesToMass(substance_formula, moles):
    """
    Returns the moles of the substance.

    Arguments:
        substance_formula (string): A chemical formula (molecular formula)
        mass (float): The moles of the substance
    Returns:
        moles (float): The moles of the element 
    """
    const = getMolarMass(substance_formula)
    moles = moles * const
    return moles


def massToMolarity(substance_formula, mass, volume):
    """
    Get molarity of a solution from the dissolved substance's mass and volume

    Arguments:
        substance_formula (string): Solute formula
        mass (string): Mass of solute
        volume (float): Volume of solvent
    """
    return massToMoles(substance_formula, mass) / volume

def molarityToMass(substance_formula, moles, volume):
    """
    Inverse of massToMolarity
    """
    return volume * molesToMass(substance_formula, moles)


def percentMass(base_substance, target_substance):
    """
    Percent mass target_substance of base_substance

    Argument:
        base_substance (string): the formual of the chemical substance
        target_substance (string): a subset of the substance formula
    Return:
        (float): percent mass
    """
    mass_part = getMolarMass(base_substance, target_substance)
    mass_whole = getMolarMass(base_substance)
    return mass_part / mass_whole


def massCompToMoleComp(mass_composition):
    """
    Return a dictionary with the moles of each element in a compound of 100g.

    Arg:
        mass_composition (dict): A mapping between elements and their percentage in the compound
    Returns
        mole_comp (dict): a dict with the percent mole composition of each element
    """
    mole_comp = {}
    mole_total = 0
    for key in mass_composition:
        tmp = 100 * mass_composition[key] / getMolarMass(key)
        mole_comp[key] = tmp
        mole_total += tmp
    for key in mass_composition:
        mole_comp[key] /= mole_total
    return mole_comp

def moleRatio(term1, term2, chemicalEquationSolution):
    """
    Compute the ratio of moles of term1 to term2 in the equation whose solution is given.

    (E.g. The ratio of H2 to O2 in "2 H2 + O2 : H2O" is 2 / 1)

    Arg:
        term1 (string): A term in the chemical equation.
        term2 (string): A term in the chemical equation.
        chemicalEquationSolution (dict) A dictionary with the solution to a chemical equation.
    Returns:
        ratio (float): The ratio of term1 to term2 in the reaction described by the solution
    """
    return chemicalEquationSolution[term1] / chemicalEquationSolution[term2]

def productsFromMoles(equation, reactant_amounts):
    """
    Compute the moles of the products produced in a chemical reaction.

    Arguments:
        equation (string): An equation that represents the chemical reaction.
        reactant_amounts (dict): A dictionary with the amount (moles) of each reactant.
    Return 
        product_amounts (dict): A dictionary with the amount (moles) of each product.
    """
    # Solve chemical equation
    coef = BCE.solutionCoefficients(equation)
    # Find limiting reactant
    react_amt = len(reactant_amounts)
    react_terms = list(reactant_amounts.keys())[:react_amt]
    lim_reagent, moles = react_terms[0], reactant_amounts[react_terms[0]]
    for key in react_terms[1:]:
        mole_equivalence = reactant_amounts[key] * moleRatio(lim_reagent, key, coef)
        if mole_equivalence < moles:
            lim_reagent = key
            moles = reactant_amounts[key]
    # Get product keys
    coef = BCE.solutionCoefficients(equation)
    prod_key = list(coef.keys())[react_amt:]
    # Get moles of product
    product_amounts = {}
    for key in prod_key:
        product_amounts[key] = moles * moleRatio(key, lim_reagent, coef)
    return product_amounts
        
def productsFromMass(equation, reactant_amounts):
    """
    Compute the mass of the products produced in a chemical reaction between two quantities of reactants.

    Arguments:
        equation (string): An equation that represents the chemical reaction.
        reactant_amounts (dict): A dictionary with the amount (grams) of each reactant.
    Return 
        product_amounts (dict): A dictionary with the amount (grams) of each product.
    """
    reactant_moles = {}
    for key, value in reactant_amounts.items():
        reactant_moles[key] = massToMoles(key, value)
    solution_moles = productsFromMoles(equation, reactant_moles)
    for key, value in solution_moles.items():
        solution_moles[key] = molesToMass(key, value)
    return solution_moles
    
def reactantsFromMoles(equation, product_amounts):
    """
    Compute the moles of the reactants required for a specific amount of a product.

    Arguments:
        equation (string): An equation that represents the chemical reaction.
        product_amounts (dict): A dictionary with the amount (moles) of a reactant.
    Return 
        reactant_amounts (dict): A dictionary with the amount (moles) of each product.
    """
    # Solve chemical equation
    coef = BCE.solutionCoefficients(equation)
    # Find limiting product
    prod_amt = len(product_amounts)
    prod_terms = list(reversed(list(product_amounts.keys())))[:prod_amt]
    lim_reagent, moles = prod_terms[0], product_amounts[prod_terms[0]]
    for key in prod_terms[1:]:
        mole_equivalence = product_amounts[key] * moleRatio(lim_reagent, key, coef)
        if mole_equivalence < moles:
            lim_reagent = key
            moles = product_amounts[key]
    # Get reactant keys
    coef = BCE.solutionCoefficients(equation)
    prod_key = list(reversed(list(coef.keys())))[prod_amt:]
    # Get moles of product
    reactant_amount = {}
    for key in prod_key:
        reactant_amount[key] = moles * moleRatio(key, lim_reagent, coef)
    return reactant_amount


def calcTheoreticalProduct(equation, target_product, product_amounts):
    """
    Calculate the mass of the target_product in the reaction described by equation, given the product_amounts.

    Arguments:
        equation (str): An equation representing the chemical reaction (unbalanced).
        target_product (str): Formula of the compound to be found.
        product_amounts (str): A dictionary with formula, mass (g) pairs for the products
    Return:
        theoretical_mass (float): The theoretical mass of the target_product
    """
    sols = BCE.solutionCoefficients(equation)
    theoretical_mass = np.inf
    for product in product_amounts:
        tmp = massToMoles(product, product_amounts[product]) * sols[target_product] / sols[product] * getMolarMass(target_product)
        if tmp < theoretical_mass:
            theoretical_mass= tmp
    return theoretical_mass

def getLimitingReagent(equation, reactant_amts):
    """
    Find the limiting reactant.
    
    Arguments:
        equation (str): The unbalanced chemical equation.
        reactant_amts (dict): A dictionary with formula, mole pair values for the reactants.
    Return:
        limit_reagent: The formula of the limiting reagent.
    """
    sols = BCE.solutionCoefficients(equation)
    lim, lim_store = "", np.inf
    base_key = list(sols.keys())[0]
    for ingredient in reactant_amts:
        tmp = reactant_amts[ingredient] * sols[base_key] / sols[ingredient]
        if tmp < lim_store:
            lim = ingredient
            lim_store = reactant_amts[ingredient]
    return lim
