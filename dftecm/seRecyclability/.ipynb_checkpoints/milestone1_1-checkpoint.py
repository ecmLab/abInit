from pymatgen.ext.matproj import MPRester
from pymatgen.analysis.reaction_calculator import ComputedReaction

# Initialize Materials Project API
api_key = "H9GD9SySZmkHxPyxI3T5ylATO9Z8XFAd"  # Replace with your Materials Project API Key
mpr = MPRester(api_key)

# Define the materials and their hydrolysis reactions
materials = {
    "Li3PS4": {
        "reactants": ["mp-696444", "mp-697111"],  # Li3PS4 + H2O
        "products": ["mp-697113", "mp-697112"],  # Li3PO4 + H2S
    },
    "Li3YCl6": {
        "reactants": ["mp-XXXXX", "mp-697111"],  # Li3YCl6 + H2O (Placeholder ID)
        "products": ["mp-YYYYY", "mp-697113"],  # YCl3 + LiOH (Placeholder ID)
    },
    "LLZO": {
        "reactants": ["mp-XXXXX", "mp-697111"],  # LLZO + H2O (Placeholder ID)
        "products": ["mp-YYYYY", "mp-YYYYY"],  # Products Placeholder ID
    },
    "Na3SbS4": {
        "reactants": ["mp-XXXXX", "mp-697111"],  # Na3SbS4 + H2O (Placeholder ID)
        "products": ["mp-YYYYY", "mp-YYYYY"],  # Products Placeholder ID
    },
    "Li3InCl6": {
        "reactants": ["mp-XXXXX", "mp-697111"],  # Li3InCl6 + H2O (Placeholder ID)
        "products": ["mp-YYYYY", "mp-YYYYY"],  # Products Placeholder ID
    },
}

# Function to calculate hydrolysis reaction energy
def calculate_hydrolysis_energy(material_key, material_data):
    reactants = [mpr.get_entry_by_material_id(mid) for mid in material_data["reactants"]]
    products = [mpr.get_entry_by_material_id(mid) for mid in material_data["products"]]
    reaction = ComputedReaction(reactants, products)
    return reaction.calculated_reaction_energy

# Iterate over materials to compute hydrolysis reaction energies
for material, data in materials.items():
    try:
        energy = calculate_hydrolysis_energy(material, data)
        print(f"Hydrolysis Reaction Energy for {material}: {energy:.3f} eV")
    except Exception as e:
        print(f"Error calculating energy for {material}: {e}")
