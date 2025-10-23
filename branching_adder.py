#!/usr/bin/env python3
"""
Extract branching ratios from JSON file and apply to OpenMC depletion chain.

Simple standalone script - just set the file paths below and run.
"""

import json
import openmc.deplete
from collections import defaultdict

# ============================================================================
# CONFIGURE FILE PATHS HERE
# ============================================================================
INPUT_JSON_FILE = "branching_ratios_14MeV.json"
INPUT_CHAIN_FILE = "JENDL_chain.xml"
OUTPUT_CHAIN_FILE = "JENDL_chain_with_branching_isomers.xml"

# ============================================================================
# Element symbols mapping (Z=1 to Z=100)
# ============================================================================
ELEMENT_SYMBOLS = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O',
    9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P',
    16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti',
    23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu',
    30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr',
    37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc',
    44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
    51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La',
    58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd',
    65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu',
    72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt',
    79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At',
    86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa', 92: 'U',
    93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es',
    100: 'Fm'
}

# ============================================================================
# Complete MT number to OpenMC reaction type mapping
# ============================================================================
MT_TO_REACTION = {
    # Basic reactions
    'MT_102': '(n,gamma)',     # Radiative capture
    'MT_103': '(n,p)',         # (n,proton)
    'MT_104': '(n,d)',         # (n,deuteron)
    'MT_105': '(n,t)',         # (n,triton)
    'MT_106': '(n,3He)',       # (n,helium-3)
    'MT_107': '(n,a)',         # (n,alpha)

    # Multiple neutron emission
    'MT_16': '(n,2n)',         # (n,2n)
    'MT_17': '(n,3n)',         # (n,3n)
    'MT_37': '(n,4n)',         # (n,4n)

    # Neutron + charged particle
    'MT_22': '(n,na)',         # (n,n+alpha)
    'MT_28': '(n,np)',         # (n,n+proton)
    'MT_32': '(n,nd)',         # (n,n+deuteron)
    'MT_33': '(n,nt)',         # (n,n+triton)
    'MT_44': '(n,n2p)',        # (n,n+2protons)
    'MT_45': '(n,npa)',        # (n,n+proton+alpha)

    # Multiple neutrons + charged particle
    'MT_11': '(n,2nd)',        # (n,2n+deuteron)
    'MT_24': '(n,2na)',        # (n,2n+alpha)
    'MT_30': '(n,2n2a)',       # (n,2n+2alpha)
    'MT_34': '(n,n3He)',       # (n,n+helium-3)
    'MT_41': '(n,2np)',        # (n,2n+proton)

    # Multiple charged particles
    'MT_108': '(n,2a)',        # (n,2alpha)
    'MT_109': '(n,3a)',        # (n,3alpha)
    'MT_111': '(n,2p)',        # (n,2proton)
    'MT_112': '(n,pa)',        # (n,proton+alpha)
    'MT_115': '(n,pd)',        # (n,proton+deuteron)
    'MT_116': '(n,pt)',        # (n,proton+triton)
    'MT_117': '(n,da)',        # (n,deuteron+alpha)

    # Complex reactions
    'MT_155': '(n,ta)',        # (n,triton+alpha)

    # Special reactions
    'MT_197': '(n,3p)',        # (n,3protons)
}


def parse_parent_notation(parent_str):
    """Convert '(Z, A, m)' string to OpenMC nuclide name."""
    parent_str = parent_str.strip('()')
    parts = [int(x.strip()) for x in parent_str.split(',')]

    # Handle both (Z, A) and (Z, A, m) formats
    if len(parts) == 2:
        z, a = parts
        m = 0  # Default to ground state
    elif len(parts) == 3:
        z, a, m = parts
    else:
        raise ValueError(f"Invalid parent notation: {parent_str}")

    element = ELEMENT_SYMBOLS.get(z)
    if not element:
        raise ValueError(f"Unknown element with Z={z}")

    nuclide = f"{element}{a}"

    # Add metastable state suffix if isomer
    if m == 1:
        nuclide += "_m1"
    elif m == 2:
        nuclide += "_m2"
    elif m == 3:
        nuclide += "_m3"
    elif m > 3:
        nuclide += f"_m{m}"

    return nuclide


def get_ground_state_q_value(nuclide, reaction_type, ground_state_product):
    """
    Extract Q-value from existing ground state reaction.

    Parameters:
    - nuclide: openmc.deplete.Nuclide object
    - reaction_type: str like "(n,2n)"
    - ground_state_product: str like "Kr79"

    Returns:
    - Q-value in eV, or None if not found
    """
    for reaction in nuclide.reactions:
        if reaction.type == reaction_type and reaction.target == ground_state_product:
            return reaction.Q
    return None


def zap_to_nuclide(zap, isomer_state):
    """Convert ZAP number and isomer state to OpenMC nuclide name."""
    z = int(zap // 1000)
    a = int(zap % 1000)

    element = ELEMENT_SYMBOLS.get(z)
    if not element:
        raise ValueError(f"Unknown element with Z={z}")

    nuclide = f"{element}{a}"

    # Add metastable state suffix if needed
    if isomer_state == "m1":
        nuclide += "_m1"
    elif isomer_state == "m2":
        nuclide += "_m2"
    elif isomer_state == "m3":
        nuclide += "_m3"
    elif isomer_state and isomer_state.startswith('m'):
        nuclide += f"_{isomer_state}"

    return nuclide


def add_missing_isomeric_products(chain, branching_data):
    """
    Add isomeric products to chain reactions that don't already have them.
    Uses same Q-value as ground state reaction.

    Parameters:
    - chain: openmc.deplete.Chain object
    - branching_data: dict from extract_branching_ratios()

    Returns:
    - dict with statistics
    """
    print(f"\n{'='*70}")
    print(f"Adding Missing Isomeric Products")
    print(f"{'='*70}")

    stats = {
        'total_isomers_to_add': 0,
        'successfully_added': 0,
        'already_exists': 0,
        'no_q_value': 0,
        'parent_missing': 0,
        'failed': 0
    }

    for reaction_type, parent_data in sorted(branching_data.items()):
        for parent, products in parent_data.items():
            # Check if parent exists
            if parent not in chain.nuclide_dict:
                stats['parent_missing'] += len([p for p in products.values() if p['is_isomer']])
                continue

            parent_nuclide = chain.nuclides[chain.nuclide_dict[parent]]

            # Find ground state product and its Q-value
            ground_state_product = None
            ground_state_q = None

            for product_name, product_data in products.items():
                if not product_data['is_isomer']:
                    ground_state_product = product_name
                    ground_state_q = get_ground_state_q_value(
                        parent_nuclide, reaction_type, ground_state_product
                    )
                    break

            # If no ground state in JSON, try to find any existing reaction of this type
            if ground_state_q is None:
                for reaction in parent_nuclide.reactions:
                    if reaction.type == reaction_type and reaction.Q is not None:
                        ground_state_q = reaction.Q
                        break

            # Add isomeric products
            for product_name, product_data in products.items():
                if not product_data['is_isomer']:
                    continue  # Skip ground state

                stats['total_isomers_to_add'] += 1

                # Check if this isomeric product already exists in reactions
                already_exists = False
                for reaction in parent_nuclide.reactions:
                    if reaction.type == reaction_type and reaction.target == product_name:
                        already_exists = True
                        break

                if already_exists:
                    stats['already_exists'] += 1
                    continue

                # Need Q-value to add reaction
                if ground_state_q is None:
                    print(f"  Warning: No Q-value for {parent} {reaction_type} → {product_name}")
                    stats['no_q_value'] += 1
                    continue

                # Add the isomeric product using same Q as ground state
                try:
                    parent_nuclide.add_reaction(
                        type=reaction_type,
                        target=product_name,
                        Q=ground_state_q,  # Use same Q as ground state for now
                        branching_ratio=0.0  # Will be set later by set_branch_ratios
                    )
                    print(f"  Added: {parent} {reaction_type} → {product_name} (Q={ground_state_q/1e6:.2f} MeV)")
                    stats['successfully_added'] += 1
                except Exception as e:
                    print(f"  Error adding {parent} {reaction_type} → {product_name}: {e}")
                    stats['failed'] += 1

    # Print summary
    print(f"\n{'='*70}")
    print(f"Isomeric Product Addition Summary:")
    print(f"{'='*70}")
    print(f"Total isomers to add: {stats['total_isomers_to_add']}")
    print(f"Successfully added: {stats['successfully_added']}")
    print(f"Already existed: {stats['already_exists']}")
    print(f"Skipped (no Q-value): {stats['no_q_value']}")
    print(f"Skipped (parent missing): {stats['parent_missing']}")
    print(f"Failed: {stats['failed']}")

    return stats


def extract_branching_ratios(json_file_path):
    """Extract branching ratios from JSON and format for OpenMC."""

    print(f"\nLoading branching ratio data from: {json_file_path}")
    with open(json_file_path, 'r') as f:
        data = json.load(f)

    # Organize results by reaction type
    branching_by_reaction = defaultdict(lambda: defaultdict(dict))

    # Statistics
    stats = {
        'total_parents': 0,
        'total_reactions': 0,
        'reactions_with_branching': 0,
        'zero_sum_skipped': 0,
        'sum_errors': [],
        'unknown_mt': set(),
        'parse_errors': []
    }

    for parent_key, reactions in data.items():
        if not reactions:
            continue

        stats['total_parents'] += 1

        try:
            parent_nuclide = parse_parent_notation(parent_key)
        except ValueError as e:
            stats['parse_errors'].append(f"Parent {parent_key}: {e}")
            continue

        for mt_key, products in reactions.items():
            if mt_key not in MT_TO_REACTION:
                stats['unknown_mt'].add(mt_key)
                continue

            reaction_type = MT_TO_REACTION[mt_key]
            stats['total_reactions'] += 1

            if len(products) > 1:
                stats['reactions_with_branching'] += 1

            # Collect products in temporary dict first
            temp_products = {}
            total_br = 0.0

            for product in products:
                zap = product['zap']
                isomer = product.get('matched_isomer', 'm=0')
                br = product['branching_ratio_14MeV']
                elfs_ev = product.get('elfs_ev', 0.0)

                # Convert isomer notation
                if isomer == 'm=0':
                    isomer_state = None
                elif isomer.startswith('m'):
                    isomer_state = isomer
                else:
                    isomer_state = None

                try:
                    product_nuclide = zap_to_nuclide(zap, isomer_state)
                except ValueError as e:
                    stats['parse_errors'].append(f"Product ZAP {zap}: {e}")
                    continue

                # Store in temporary dict, summing branching ratios for duplicate products
                if product_nuclide in temp_products:
                    # Product already exists - add to its branching ratio
                    temp_products[product_nuclide]['branching_ratio'] += br
                else:
                    # New product - create entry
                    temp_products[product_nuclide] = {
                        'branching_ratio': br,
                        'elfs_ev': elfs_ev,
                        'is_isomer': (isomer_state is not None)
                    }
                total_br += br

            # Skip reactions with zero branching ratio sum
            if total_br < 0.001:  # Effectively zero
                stats['zero_sum_skipped'] += 1
                continue

            # Add products to main dictionary only if sum is non-zero
            for product_nuclide, product_data in temp_products.items():
                branching_by_reaction[reaction_type][parent_nuclide][product_nuclide] = product_data

            # Validate sum for non-zero reactions
            if abs(total_br - 1.0) > 0.01:
                stats['sum_errors'].append({
                    'parent': parent_nuclide,
                    'reaction': reaction_type,
                    'sum': total_br
                })

    # Print statistics
    print(f"\n{'='*70}")
    print(f"Extraction Statistics")
    print(f"{'='*70}")
    print(f"Total parent nuclides: {stats['total_parents']}")
    print(f"Total reactions: {stats['total_reactions']}")
    print(f"Reactions with isomer branching: {stats['reactions_with_branching']}")
    print(f"Reactions skipped (zero sum): {stats['zero_sum_skipped']}")
    print(f"Unique reaction types: {len(branching_by_reaction)}")

    if stats['unknown_mt']:
        print(f"\nWarning: {len(stats['unknown_mt'])} unknown MT numbers found:")
        print(f"  {sorted(stats['unknown_mt'])}")

    if stats['parse_errors']:
        print(f"\nWarning: {len(stats['parse_errors'])} parsing errors")
        for error in stats['parse_errors'][:5]:
            print(f"  {error}")

    if stats['sum_errors']:
        print(f"\nWarning: {len(stats['sum_errors'])} reactions don't sum to 1.0:")
        # Sort by largest difference from 1.0 first
        sorted_errors = sorted(stats['sum_errors'],
                              key=lambda x: abs(x['sum'] - 1.0),
                              reverse=True)
        for error in sorted_errors[:10]:
            diff = abs(error['sum'] - 1.0)
            print(f"  {error['parent']} {error['reaction']}: sum = {error['sum']:.4f} (diff = {diff:.4f})")

    print(f"\n{'='*70}")
    print("Reactions by Type:")
    print(f"{'='*70}")
    for reaction_type in sorted(branching_by_reaction.keys()):
        parent_count = len(branching_by_reaction[reaction_type])
        print(f"  {reaction_type:20s}: {parent_count:4d} parent nuclides")

    return dict(branching_by_reaction)


def apply_branching_ratios_to_chain(chain_file, branching_data, output_file):
    """Apply branching ratios to chain and save."""

    print(f"\n{'='*70}")
    print(f"Loading chain: {chain_file}")
    print(f"{'='*70}")

    chain = openmc.deplete.Chain.from_xml(chain_file)
    print(f"Chain loaded with {len(chain.nuclides)} nuclides")

    # Add missing isomeric products before setting branching ratios
    add_stats = add_missing_isomeric_products(chain, branching_data)

    print(f"\n{'='*70}")
    print(f"Applying Branching Ratios")
    print(f"{'='*70}")

    total_applied = 0
    total_skipped_parent = 0
    total_skipped_products = 0
    total_failed = 0

    for reaction_type, parent_data in sorted(branching_data.items()):
        print(f"\n{reaction_type}:")

        applied = 0
        skipped_parent = 0
        skipped_products = 0
        failed = 0

        for parent, products in parent_data.items():
            # Check if parent exists in chain
            if parent not in chain.nuclide_dict:
                skipped_parent += 1
                continue

            # Extract just the branching ratios (products is now dict of dicts)
            products_br = {
                product_name: product_data['branching_ratio']
                for product_name, product_data in products.items()
            }

            # Check if all products exist
            missing_products = [p for p in products_br if p not in chain.nuclide_dict]
            if missing_products:
                skipped_products += 1
                continue

            try:
                chain.set_branch_ratios(
                    {parent: products_br},
                    reaction=reaction_type,
                    strict=False
                )
                applied += 1
            except Exception as e:
                print(f"  Error for {parent}: {e}")
                failed += 1

        print(f"  Applied: {applied}")
        if skipped_parent > 0:
            print(f"  Skipped (parent not in chain): {skipped_parent}")
        if skipped_products > 0:
            print(f"  Skipped (products not in chain): {skipped_products}")
        if failed > 0:
            print(f"  Failed: {failed}")

        total_applied += applied
        total_skipped_parent += skipped_parent
        total_skipped_products += skipped_products
        total_failed += failed

    print(f"\n{'='*70}")
    print(f"Overall Summary:")
    print(f"{'='*70}")
    print(f"Total applied: {total_applied}")
    print(f"Total skipped (parent not in chain): {total_skipped_parent}")
    print(f"Total skipped (products not in chain): {total_skipped_products}")
    print(f"Total failed: {total_failed}")

    # Export
    print(f"\n{'='*70}")
    print(f"Saving modified chain to: {output_file}")
    print(f"{'='*70}")

    chain.export_to_xml(output_file)
    print(f"\nSuccess! Modified chain saved.\n")

    return chain


# ============================================================================
# MAIN SCRIPT EXECUTION
# ============================================================================
if __name__ == "__main__":
    print("\n" + "="*70)
    print("OpenMC Branching Ratio Extraction and Application")
    print("="*70)

    # Extract branching ratios from JSON
    branching_data = extract_branching_ratios(INPUT_JSON_FILE)

    # Apply to chain file
    chain = apply_branching_ratios_to_chain(
        INPUT_CHAIN_FILE,
        branching_data,
        OUTPUT_CHAIN_FILE
    )
