from BalChemEq import BCE
import getMoles 


def findConcentrationFromPrecipitate(equation, sample, mass_obtained):
    """
    Find the molarity (g / L) of a sample given the mass of precipitate formed in a reaction.

    Arguments:
        equation (str): The unbalanced chemical equation
        sample (tuple): A tuple containing the formula of the solute and the volume of solution (L)
        mass_obtained (tuple): A tuple containing the symbol for the precipitate and its mass (g)
    Return:
        molarity (float): The molarity (g / L) of the sample
    """
    target_substance = sample[0]
    sol = BCE.solutionCoefficients(equation)
    target_mass = (
                getMoles.massToMoles(mass_obtained[0], mass_obtained[1]) * 
                sol[target_substance] / sol[mass_obtained[0]] * 
                getMoles.getMolarMass(target_substance)
                )
    molarity = target_mass / sample[1]
    return molarity

def getNeutralizingMass(equation, target, present):
    """
    Find the mass of a base/acid needed to neutralize an acid/base.

    Arguments:
        equation (str): The unbalanced chemical equation for the neutralization reaction.
        target (str): The formula for the acid/base to find.
        present (tuple): A tuple containing the formula and moles of the base/acid present.
    Return:
        y_mass (float): The mass g required to neutralize x.
    
    
    """
    x_eq = present[0]
    x_mol = present[1]
    sol = BCE.solutionCoefficients(equation)
    y_mass = x_mol * sol[target] / sol[x_eq] * getMoles.getMolarMass(target)
    return y_mass

def getTitrantVolume(equation, analyte, titrant):
    """
    Get the volume of a titrant required to neutralize the analyte.

    Arguments:
        equation (str): Unbalanced chemical equation.
        analyte (tuple): Formula and mass (g) of the analyte.
        titrant (tuple): Formula and molarity of the titrant.

    Return:
        titrant_volume (float): The volume of the titrant (L).
    """
    sols = BCE.solutionCoefficients(equation)
    titrant_volume = getMoles.massToMoles(analyte[0],analyte[1]) * sols[titrant[0]] / sols[analyte[0]] * (1 / titrant[1])
    return titrant_volume

def findIonMolarity(equation, solute, solvent):
    """
    Find the molarity of a cation of the solute in a solution after a reaction.

    Arguments:
        equation (str): Unbalanced chemical equation.
        solute (tuple): Formula and mass (g) of the solute.
        titrant (tuple): Formula, volume, and molarity of the solvent.

    Return:
        ratio (float): The concentration of unreacted solute cations.
    """
    BCE.writeBCESteps(equation)
    sols = BCE.solutionCoefficients(equation)

    reactants = {
        solute[0]: getMoles.massToMoles(solute[0],solute[1]),
        solvent[0]: solvent[1]*solvent[2]
    }

    limiting = getMoles.getLimitingReagent(equation, reactants)
    if limiting == solute[0]:
        print(0)
    else:
        unreacted = getMoles.massToMoles(solute[0], solute[1]) - (sols[solute[0]] / sols[limiting]) * reactants[limiting]
        print(unreacted / solvent[1])

    print(getMoles.massToMoles(solute[0], solute[1]) / 0.1)


def massPercentAnalyte(eq_analyte, grams_solution, titrant, mole_ratio):
    """
    Find the mass  percentage of an analyte in a sample.

    Arguments:
        eq_analyte (str): The chemical equation for the analyte.
        grams_solution (float): The mass of the solution in grams.
        titrant (tuple): A tuple with the volume (L) and molarity (M) of the titrant solution.
        mole_ratio (float): The ratio of moles analyte to the titrant in the balanced equation.
    Return:
        mass_percentage (float): The ratio of the grams of the analyte to the grams of the solution.
    """
    grams_solution = 21
    liters_titrant, molarity_titrant = titrant
    moles_titrant = liters_titrant * molarity_titrant
    moles_analyte = moles_titrant * mole_ratio
    grams_analyte = getMoles.molesToMass(eq_analyte, moles_analyte)
    return grams_analyte / grams_solution

def massPercentCation(eq_target, eq_ppt, grams_ppt, target_to_sample, grams_sample, mole_ratio):
    """
    Find the mass  percentage of an analyte in a sample.

    Arguments:
        eq_target (str): The chemical equatio of the target ion.
        eq_ppt (str): The chemical equation of the precipitate.
        grams_ppt (float): The mass of the precipitate in grams.
        titrant_to_sample (float): The ratio of moles of the target to one mole of the analyte.
        mole_ratio (float): The mole ratio of analyte to the titrant in the balanced equation.
    Return:
        mass_percentage (float): The ratio of the grams of the analyte to the grams of the solution.
    """
    moles_ppt = getMoles.massToMoles(eq_ppt, grams_ppt)
    moles_analyte = moles_ppt * mole_ratio
    moles_ppt = moles_analyte * target_to_sample
    grams_ppt = getMoles.molesToMass(eq_target, moles_ppt)
    return grams_ppt / grams_sample

