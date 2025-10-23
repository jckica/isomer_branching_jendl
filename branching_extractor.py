from endf_parserpy import EndfParserPy
import pprint
import os
from typing import Dict, List, Set, Tuple, Optional, Union
import re
import numpy as np
parser = EndfParserPy()
path_to_files = "neutronics_workflow\Code\workspace\JENDL\jendl5-n\jendl5-n"
path_to_file = r"neutronics_workflow\Code\workspace\JENDL\jendl5-n\jendl5-n\n_011-Na-023.dat"
path_to_decay_files = "neutronics_workflow\Code\workspace\JENDL\jendl5-dec_upd5\jendl5-dec_upd5"
# All MT numbers from OpenMC REACTIONS dictionary
list_of_mt_to_extract = [
    # Individual MT numbers
    11, 16, 17, 22, 23, 24, 25, 28, 29, 30, 32, 33, 34, 35, 36, 37, 41, 42, 44, 45,
    102, 103, 104, 105, 106, 107, 108, 109, 111, 112, 113, 114, 115, 116, 117, 152, 153, 154, 155, 156, 157, 158, 159,
    160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179,
    180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200,
]
def parse_filename(filename: str) -> Optional[Tuple[int, int, int]]:
        """
        Parse JENDL filename to extract Z, A, and metastable state
        
        Args:
            filename: JENDL data filename (e.g., 'n_092-U-238.dat' or 'dec-092-U-238m1.dat')
            
        Returns:
            Tuple of (Z, A, metastable_state) or None if parsing fails 
        """
        # Pattern for neutron files: n_ZZZ-Symbol-AAA.dat or n_ZZZ-Symbol-AAAm?.dat
        neutron_pattern = r'n_(\d{3})-([A-Z][a-z]?)-(\d{3})(m\d+)?\.dat'
        
        # Pattern for decay files: dec-ZZZ-Symbol-AAA.dat or dec-ZZZ-Symbol-AAAm?.dat  
        decay_pattern = r'dec-(\d{3})-([A-Z][a-z]?)-(\d{3})(m\d+)?\.dat'
        
        for pattern in [neutron_pattern, decay_pattern]:
            match = re.match(pattern, filename)
            if match:
                z = int(match.group(1))
                symbol = match.group(2)
                a = int(match.group(3))
                metastable = match.group(4)
                
                # Parse metastable state
                m = 0
                if metastable:
                    m = int(metastable[1:])  # Extract number after 'm'
                
                return z, a, m
        
        return None
def get_za_combinations_as_tuples(folder_path: str) -> List[Tuple[Tuple[int, int], Tuple[int, ...]]]:
    """
    Parse a folder of files and return a list of tuples where each tuple contains
    (Z, A) and all isomeric states for that nuclide.
    
    Args:
        folder_path: Path to the folder containing JENDL data files
        
    Returns:
        List of tuples where each tuple is ((Z, A), (m1, m2, m3, ...))
        representing all isomeric states for each nuclide
    """
    za_combinations = {}
    
    try:
        files = os.listdir(folder_path)
    except FileNotFoundError:
        print(f"Folder not found: {folder_path}")
        return []
    
    for filename in files:
        if not filename.endswith('.dat'):
            continue
            
        parsed = parse_filename(filename)
        if parsed is not None:
            z, a, m = parsed
            za_key = (z, a)
            
            if za_key not in za_combinations:
                za_combinations[za_key] = []
            
            if m not in za_combinations[za_key]:
                za_combinations[za_key].append(m)
    
    # Sort the metastable states for each ZA combination
    for za_key in za_combinations:
        za_combinations[za_key].sort()
    
    # Convert to list of tuples
    result = []
    for (z, a), metastable_states in za_combinations.items():
        result.append(((z, a), tuple(metastable_states)))
    
    # Sort by Z, then A for consistent ordering
    result.sort(key=lambda x: (x[0][0], x[0][1]))
    
    return result

def get_decay_filename_for_isomer(z: int, a: int, m: int) -> str:
    """
    Build decay filename from Z, A, metastable state.
    Uses same pattern as parse_filename() in reverse.

    Args:
        z: Atomic number
        a: Mass number
        m: Metastable state (0, 1, 2, ...)

    Returns:
        Decay filename (e.g., "dec-030-Zn-061m1.dat")
    """
    ELEMENTS = {
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
        100: 'Fm', 101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db',
        106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt', 110: 'Ds', 111: 'Rg',
        112: 'Cn', 113: 'Nh', 114: 'Fl', 115: 'Mc', 116: 'Lv', 117: 'Ts', 118: 'Og'
    }

    symbol = ELEMENTS.get(z, 'Xx')
    if m == 0:
        return f"dec-{z:03d}-{symbol}-{a:03d}.dat"
    else:
        return f"dec-{z:03d}-{symbol}-{a:03d}m{m}.dat"

def extract_elis_from_decay_file(decay_file_path: str) -> Optional[float]:
    """
    Extract ELIS (excitation energy of isomeric state) from decay file MF=1, MT=451.

    Args:
        decay_file_path: Full path to decay file

    Returns:
        ELIS in eV, or None if not found or file doesn't exist
    """
    if not os.path.exists(decay_file_path):
        return None

    try:
        parser_local = EndfParserPy()
        decay_data = parser_local.parsefile(decay_file_path)

        # Check for MF=1, MT=451 (General Information)
        if 1 not in decay_data or 451 not in decay_data[1]:
            return None

        mt451 = decay_data[1][451]
        elis = mt451.get('ELIS', None)

        if elis is not None:
            return float(elis)

        return None
    except Exception as e:
        print(f"  Warning: Error reading decay file {os.path.basename(decay_file_path)}: {e}")
        return None

def match_elfs_to_isomer(elfs_ev: float, zap: float, decay_za_combinations: List,
                         path_to_decay_files: str, tolerance_ev: float = 1000.0) -> dict:
    """
    Match ELFS energy to isomeric state by comparing against ELIS energies.

    Args:
        elfs_ev: Excitation energy from parent MF=8 (eV)
        zap: Product nuclide identifier
        decay_za_combinations: List of available decay file (Z,A) and isomers
        path_to_decay_files: Path to decay file directory
        tolerance_ev: Energy matching tolerance (default ±1 keV)

    Returns:
        Dictionary with matching results:
        {
            'matched_isomer': 'm=0' or 'm1' or 'm2' or None,
            'elis_ev': float or None,
            'energy_diff_ev': float or None,
            'match_quality': 'exact' | 'within_tolerance' | 'no_match' | 'no_decay_file'
        }
    """
    za = int(zap)
    z = za // 1000
    a = za % 1000

    # Ground state always matches ELFS=0.0
    if abs(elfs_ev) < 0.1:  # Effectively zero
        return {
            'matched_isomer': 'm=0',
            'elis_ev': 0.0,
            'energy_diff_ev': 0.0,
            'match_quality': 'exact'
        }

    # Find available isomers for this (Z,A) in decay library
    available_isomers = []
    for (decay_z, decay_a), isomers in decay_za_combinations:
        if decay_z == z and decay_a == a:
            available_isomers = list(isomers)
            break

    if not available_isomers:
        # No decay files available - assign to ground state
        return {
            'matched_isomer': 'm=0',
            'elis_ev': None,
            'energy_diff_ev': None,
            'match_quality': 'no_decay_file'
        }

    # Try matching against each isomer's ELIS energy
    best_match = None
    best_diff = float('inf')

    for m in available_isomers:
        if m == 0:
            continue  # Skip ground state (already handled)

        # Build decay filename
        decay_filename = get_decay_filename_for_isomer(z, a, m)
        decay_filepath = os.path.join(path_to_decay_files, decay_filename)

        # Extract ELIS energy
        elis_energy = extract_elis_from_decay_file(decay_filepath)

        if elis_energy is None:
            continue

        # Calculate energy difference
        diff = abs(elfs_ev - elis_energy)

        if diff < best_diff:
            best_diff = diff
            best_match = {
                'matched_isomer': f'm{m}',
                'elis_ev': elis_energy,
                'energy_diff_ev': diff,
                'match_quality': 'exact' if diff < 1.0 else ('within_tolerance' if diff <= tolerance_ev else 'poor_match')
            }

    # If no match within tolerance, assign to ground state
    if best_match is None or best_diff > tolerance_ev:
        return {
            'matched_isomer': 'm=0',
            'elis_ev': best_match['elis_ev'] if best_match else None,
            'energy_diff_ev': best_diff if best_match else None,
            'match_quality': 'no_match'
        }

    return best_match

def endf_interpolate(x_new, x, y, nbt, interp_codes):
    """
    Interpolate using ENDF-6 interpolation schemes (vectorized implementation).

    Args:
        x_new: array of x values to interpolate at
        x: array of tabulated x values
        y: array of tabulated y values
        nbt: array of interpolation region boundaries (indices)
        interp_codes: array of interpolation scheme codes for each region
            1 = histogram (constant)
            2 = linear-linear
            3 = linear-log (y linear in ln(x))
            4 = log-linear (ln(y) linear in x)
            5 = log-log (ln(y) linear in ln(x))
            6 = Gamow (not implemented, falls back to linear)

    Returns:
        y_new: interpolated y values
    """
    x = np.asarray(x)
    y = np.asarray(y)
    x_new = np.asarray(x_new)
    nbt = np.asarray(nbt, dtype=int)
    interp_codes = np.asarray(interp_codes, dtype=int)

    y_new = np.zeros_like(x_new, dtype=float)

    # Handle out of range values
    out_of_range = (x_new < x[0]) | (x_new > x[-1])
    y_new[out_of_range] = 0.0

    # Work only with in-range values
    in_range = ~out_of_range
    if not np.any(in_range):
        return y_new

    x_valid = x_new[in_range]

    # Find which interpolation region each x_valid belongs to
    boundary_x = x[nbt - 1]  # NBT is 1-indexed, convert to 0-indexed
    region_indices = np.searchsorted(boundary_x, x_valid)
    region_indices = np.clip(region_indices, 0, len(nbt) - 1)

    # Find bracketing points for all x_valid at once
    idx = np.searchsorted(x, x_valid) - 1
    idx = np.clip(idx, 0, len(x) - 2)

    # Extract bracketing values (now arrays)
    x1 = x[idx]
    x2 = x[idx + 1]
    y1 = y[idx]
    y2 = y[idx + 1]

    # Get interpolation code for each point
    point_interp_codes = interp_codes[region_indices]

    # Initialize result array for valid points
    y_valid = np.zeros_like(x_valid)

    # Handle exact matches
    exact_lower = (x_valid == x1)
    exact_upper = (x_valid == x2)
    y_valid[exact_lower] = y1[exact_lower]
    y_valid[exact_upper] = y2[exact_upper]

    # Points that need interpolation (not exact matches)
    needs_interp = ~(exact_lower | exact_upper)

    # Process each interpolation scheme
    for code in np.unique(point_interp_codes):
        # Mask for points using this interpolation code that need interpolation
        mask = needs_interp & (point_interp_codes == code)
        if not np.any(mask):
            continue

        x_m = x_valid[mask]
        x1_m = x1[mask]
        x2_m = x2[mask]
        y1_m = y1[mask]
        y2_m = y2[mask]

        if code == 1:  # Histogram
            y_valid[mask] = y1_m

        elif code == 2:  # Linear-linear
            frac = (x_m - x1_m) / (x2_m - x1_m)
            y_valid[mask] = y1_m + frac * (y2_m - y1_m)

        elif code == 3:  # Linear-log (y linear in ln(x))
            # Check for valid log arguments
            valid_log = (x1_m > 0) & (x2_m > 0)

            # Linear-log for valid points
            if np.any(valid_log):
                submask = mask.copy()
                submask[mask] = valid_log
                x_s = x_valid[submask]
                x1_s = x1[submask]
                x2_s = x2[submask]
                y1_s = y1[submask]
                y2_s = y2[submask]

                frac = (np.log(x_s) - np.log(x1_s)) / (np.log(x2_s) - np.log(x1_s))
                y_valid[submask] = y1_s + frac * (y2_s - y1_s)

            # Fallback to linear for invalid points
            if np.any(~valid_log):
                submask = mask.copy()
                submask[mask] = ~valid_log
                x_s = x_valid[submask]
                x1_s = x1[submask]
                x2_s = x2[submask]
                y1_s = y1[submask]
                y2_s = y2[submask]

                frac = (x_s - x1_s) / (x2_s - x1_s)
                y_valid[submask] = y1_s + frac * (y2_s - y1_s)

        elif code == 4:  # Log-linear (ln(y) linear in x)
            # Check for valid log arguments
            valid_log = (y1_m > 0) & (y2_m > 0)

            # Log-linear for valid points
            if np.any(valid_log):
                submask = mask.copy()
                submask[mask] = valid_log
                x_s = x_valid[submask]
                x1_s = x1[submask]
                x2_s = x2[submask]
                y1_s = y1[submask]
                y2_s = y2[submask]

                frac = (x_s - x1_s) / (x2_s - x1_s)
                y_valid[submask] = y1_s * np.exp(frac * np.log(y2_s / y1_s))

            # Fallback to linear for invalid points
            if np.any(~valid_log):
                submask = mask.copy()
                submask[mask] = ~valid_log
                x_s = x_valid[submask]
                x1_s = x1[submask]
                x2_s = x2[submask]
                y1_s = y1[submask]
                y2_s = y2[submask]

                frac = (x_s - x1_s) / (x2_s - x1_s)
                y_valid[submask] = y1_s + frac * (y2_s - y1_s)

        elif code == 5:  # Log-log
            # Check for valid log arguments
            valid_log = (x1_m > 0) & (x2_m > 0) & (y1_m > 0) & (y2_m > 0)

            # Log-log for valid points
            if np.any(valid_log):
                submask = mask.copy()
                submask[mask] = valid_log
                x_s = x_valid[submask]
                x1_s = x1[submask]
                x2_s = x2[submask]
                y1_s = y1[submask]
                y2_s = y2[submask]

                frac = (np.log(x_s) - np.log(x1_s)) / (np.log(x2_s) - np.log(x1_s))
                y_valid[submask] = y1_s * np.exp(frac * np.log(y2_s / y1_s))

            # Fallback to linear for invalid points
            if np.any(~valid_log):
                submask = mask.copy()
                submask[mask] = ~valid_log
                x_s = x_valid[submask]
                x1_s = x1[submask]
                x2_s = x2[submask]
                y1_s = y1[submask]
                y2_s = y2[submask]

                frac = (x_s - x1_s) / (x2_s - x1_s)
                y_valid[submask] = y1_s + frac * (y2_s - y1_s)

        else:  # Code 6 (Gamow) and unknown codes - fallback to linear
            frac = (x_m - x1_m) / (x2_m - x1_m)
            y_valid[mask] = y1_m + frac * (y2_m - y1_m)

    # Place valid results back into output array
    y_new[in_range] = y_valid

    return y_new

def branching_ratios(lfs_data):
    """
    lfs_data: dict {LFS: (energies, xs, nbt, interp_codes)}
    energies in eV, xs in barns or multiplicities
    nbt: interpolation region boundaries
    interp_codes: interpolation scheme codes
    """
    # build union grid
    union_E = np.unique(np.concatenate([E for E, _, _, _ in lfs_data.values()]))

    # Ensure 14 MeV (14,000,000 eV) is included if within the energy range
    min_E = np.min(union_E)
    max_E = np.max(union_E)
    fusion_energy = 14_000_000  # 14 MeV in eV

    if min_E <= fusion_energy <= max_E and fusion_energy not in union_E:
        union_E = np.append(union_E, fusion_energy)
        union_E = np.sort(union_E)

    # interpolate each LFS using ENDF-aware interpolation
    interp_xs = {}
    for LFS, (E, xs, nbt, interp_codes) in lfs_data.items():
        interp_xs[LFS] = endf_interpolate(union_E, E, xs, nbt, interp_codes)

    # compute totals
    total = sum(interp_xs.values())

    # branching ratios - avoid division by zero
    br = {}
    for LFS in interp_xs:
        # Use np.divide with where to avoid warnings, set to 0 where total is 0
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = np.divide(interp_xs[LFS], total)
            # Where total is zero or very small, set branching ratio to 0
            ratio = np.where(np.abs(total) < 1e-30, 0.0, ratio)
        br[LFS] = ratio

    return union_E, br

def process_mf9_data(data, mt, debug=False):
    """
    Process MF=9 data and create lfs_data dictionary for branching ratio calculation
    MF=9 contains multiplicities (Y) instead of cross sections (sigma)
    """
    try:
        subsections = data[9][mt]['subsection']
        num_states = len(subsections)

        lfs_data = {}
        for state in range(1, num_states+1):
            # Debug: Print available keys in first subsection
            if debug and state == 1:
                print(f"\n=== DEBUG: MF9 MT{mt} Subsection Structure ===")
                print(f"Available keys in subsection[{state}]: {list(subsections[state].keys())}")
                for key in subsections[state].keys():
                    value = subsections[state][key]
                    if isinstance(value, (list, np.ndarray)):
                        print(f"  {key}: {type(value).__name__} with {len(value)} elements")
                        if key in ['NBT', 'INT']:
                            print(f"    Values: {value}")
                    else:
                        print(f"  {key}: {value} ({type(value).__name__})")
                print("=" * 50)

            lfs = subsections[state]['LFS']
            energies = np.array(subsections[state]['E'])
            multiplicities = np.array(subsections[state]['Y'])  # MF=9 uses Y (multiplicities)
            nbt = subsections[state].get('NBT', [len(energies)])  # Default: single region
            interp_codes = subsections[state].get('INT', [2])  # Default: linear interpolation

            lfs_data[lfs] = (energies, multiplicities, nbt, interp_codes)

        return lfs_data
    except Exception as e:
        print(f"Error processing MF=9 data for MT {mt}: {e}")
        return None

def process_mf10_data(data, mt, debug=False):
    """
    Process MF=10 data and create lfs_data dictionary for branching ratio calculation
    MF=10 contains cross sections (sigma)
    """
    try:
        subsections = data[10][mt]['subsection']
        num_states = len(subsections)

        lfs_data = {}
        for state in range(1, num_states+1):
            # Debug: Print available keys in first subsection
            if debug and state == 1:
                print(f"\n=== DEBUG: MF10 MT{mt} Subsection Structure ===")
                print(f"Available keys in subsection[{state}]: {list(subsections[state].keys())}")
                for key in subsections[state].keys():
                    value = subsections[state][key]
                    if isinstance(value, (list, np.ndarray)):
                        print(f"  {key}: {type(value).__name__} with {len(value)} elements")
                        if key in ['NBT', 'INT']:
                            print(f"    Values: {value}")
                    else:
                        print(f"  {key}: {value} ({type(value).__name__})")
                print("=" * 50)

            lfs = subsections[state]['LFS']
            energies = np.array(subsections[state]['E'])
            cross_sections = np.array(subsections[state]['sigma'])  # MF=10 uses sigma (cross sections)
            nbt = subsections[state].get('NBT', [len(energies)])  # Default: single region
            interp_codes = subsections[state].get('INT', [2])  # Default: linear interpolation

            lfs_data[lfs] = (energies, cross_sections, nbt, interp_codes)

        return lfs_data
    except Exception as e:
        print(f"Error processing MF=10 data for MT {mt}: {e}")
        return None

decay_za_combinations = get_za_combinations_as_tuples(path_to_decay_files)
print(decay_za_combinations)

# Process all files in the folder
branching_ratios_dict = {}

# Track all failed energy matches
failed_matches = []

try:
    files = os.listdir(path_to_files)
except FileNotFoundError:
    print(f"Folder not found: {path_to_files}")
    exit()

# Filter for .dat files only
dat_files = [f for f in files if f.endswith('.dat')]
print(f"Found {len(dat_files)} .dat files to process")

for file_idx, filename in enumerate(dat_files):
    file_path = os.path.join(path_to_files, filename)
    print(f"\nProcessing file {file_idx + 1}/{len(dat_files)}: {filename}")
    
    try:
        data = parser.parsefile(file_path)
    except Exception as e:
        print(f"Error parsing file {filename}: {e}")
        continue
    
    # Parse the filename to get ZA
    parsed_filename = parse_filename(filename)
    if not parsed_filename:
        print(f"Could not parse filename: {filename}")
        continue
    
    z, a, m = parsed_filename
    za_key = (z, a, m)  # Include metastable state to distinguish isomeric parents

    # Initialize dictionary entry for this nuclide if it doesn't exist
    if za_key not in branching_ratios_dict:
        branching_ratios_dict[za_key] = {}
    else:
        print(f"  Warning: Duplicate parent state ({z}, {a}, {m}) - data may be overwritten")
    
    for mt in list_of_mt_to_extract:
        # Check if this MT number exists in the data
        lfs_data = None

        # Enable debug mode for first file only (set to True to see interpolation schemes)
        debug_mode = False  # Change to True to see NBT/INT values

        # Try to extract MF=9 data first
        if 9 in data and mt in data[9]:
            print(f"  File: {filename}, MT: {mt} - Found MF=9 data")
            lfs_data = process_mf9_data(data, mt, debug=debug_mode)

        # If no MF=9, try MF=10
        elif 10 in data and mt in data[10]:
            print(f"  File: {filename}, MT: {mt} - Found MF=10 data")
            lfs_data = process_mf10_data(data, mt, debug=debug_mode)
            
        # If neither exists, continue to next MT
        else:
            continue
        
        # Process the data if we found any
        if lfs_data is not None:
            print(f"  Processing MT {mt} with {len(lfs_data)} states")

            # Print interpolation schemes being used
            interp_scheme_names = {
                1: "Histogram",
                2: "Linear-Linear",
                3: "Linear-Log",
                4: "Log-Linear",
                5: "Log-Log",
                6: "Gamow"
            }
            for lfs, (E, xs, nbt, interp_codes) in lfs_data.items():
                scheme_names = [interp_scheme_names.get(code, f"Unknown({code})") for code in interp_codes]
                print(f"    LFS {lfs}:")
                print(f"      NBT: {nbt}")
                print(f"      INT: {interp_codes}")
                if len(scheme_names) == 1:
                    print(f"      Scheme: {scheme_names[0]}")
                else:
                    print(f"      Schemes: {', '.join(scheme_names)}")

            br = branching_ratios(lfs_data)
            fusion_energy = 14_000_000  # 14 MeV in eV
            union_E, br_dict = br
            
            # Find the index closest to 14 MeV
            energy_idx = np.argmin(np.abs(union_E - fusion_energy))
            actual_energy = union_E[energy_idx]
            
            # Create list of tuples (LFS, branching_ratio)
            lfs_branching_tuples = []
            for lfs, br_values in br_dict.items():
                lfs_branching_tuples.append((lfs, br_values[energy_idx]))
            
            # Store in dictionary with MT as additional identifier
            mt_key = f"MT_{mt}"

            # Extract ELFS and match to isomers if MF=8 data is available
            if 8 in data and mt in data[8]:
                mt8_subsections = data[8][mt].get('subsection', {})
                enhanced_lfs_data = []

                for (lfs, br_value) in lfs_branching_tuples:
                    # Find matching LFS in MF=8 subsections
                    elfs_ev = None
                    zap = None
                    for sub_key, sub_data in mt8_subsections.items():
                        if sub_data.get('LFS') == lfs:
                            elfs_ev = sub_data.get('ELFS')
                            zap = sub_data.get('ZAP')
                            break

                    # Build enhanced entry
                    lfs_entry = {
                        'lfs': lfs,
                        'branching_ratio_14MeV': float(br_value)
                    }

                    # Perform energy matching if we have the data
                    if elfs_ev is not None and zap is not None:
                        lfs_entry['elfs_ev'] = float(elfs_ev)
                        lfs_entry['zap'] = float(zap)

                        matched_isomer_info = match_elfs_to_isomer(
                            elfs_ev, zap, decay_za_combinations,
                            path_to_decay_files, tolerance_ev=1000.0
                        )
                        lfs_entry.update(matched_isomer_info)

                        # Track failed matches
                        if lfs_entry.get('match_quality') in ['no_match', 'no_decay_file']:
                            failed_matches.append({
                                'parent_z': z,
                                'parent_a': a,
                                'mt': mt,
                                'lfs': lfs_entry['lfs'],
                                'elfs_ev': lfs_entry.get('elfs_ev'),
                                'zap': lfs_entry.get('zap'),
                                'matched_isomer': lfs_entry.get('matched_isomer'),
                                'elis_ev': lfs_entry.get('elis_ev'),
                                'energy_diff_ev': lfs_entry.get('energy_diff_ev'),
                                'match_quality': lfs_entry['match_quality']
                            })

                    enhanced_lfs_data.append(lfs_entry)

                branching_ratios_dict[za_key][mt_key] = enhanced_lfs_data

                # Print matching results
                print(f"  MT {mt}: Energy {actual_energy/1e6:.3f} MeV")
                for entry in enhanced_lfs_data:
                    if 'matched_isomer' in entry:
                        elis_str = f"{entry.get('elis_ev'):.1f}" if entry.get('elis_ev') is not None else "None"
                        diff_str = f"{entry.get('energy_diff_ev'):.1f}" if entry.get('energy_diff_ev') is not None else "None"
                        print(f"    LFS {entry['lfs']}: BR={entry['branching_ratio_14MeV']:.4f}, "
                              f"ELFS={entry.get('elfs_ev', 0):.1f} eV → {entry['matched_isomer']} "
                              f"(ELIS={elis_str} eV, diff={diff_str} eV, {entry['match_quality']})")
                    else:
                        print(f"    LFS {entry['lfs']}: BR={entry['branching_ratio_14MeV']:.4f}")
            else:
                # No MF=8 data - store simple format
                branching_ratios_dict[za_key][mt_key] = [
                    {'lfs': lfs, 'branching_ratio_14MeV': float(br_value)}
                    for lfs, br_value in lfs_branching_tuples
                ]
                print(f"  MT {mt}: Energy {actual_energy/1e6:.3f} MeV")
                print(f"  Branching ratios: {lfs_branching_tuples}")
        else:
            print(f"  Failed to process data for MT {mt}")

# Save JSON files immediately after main loop completes
import json

# Convert tuples to lists for JSON serialization
def convert_tuples_to_lists(obj):
    """Recursively convert tuples to lists for JSON serialization"""
    if isinstance(obj, tuple):
        return list(obj)
    elif isinstance(obj, dict):
        # Convert tuple keys to strings and tuple values to lists
        converted_dict = {}
        for key, value in obj.items():
            if isinstance(key, tuple):
                # Convert tuple key to string representation
                # Handle 3-element tuples (Z, A, m) for parent nuclides
                if len(key) == 3:
                    converted_key = f"({key[0]}, {key[1]}, {key[2]})"
                else:
                    converted_key = str(key)
            else:
                converted_key = key
            converted_dict[converted_key] = convert_tuples_to_lists(value)
        return converted_dict
    elif isinstance(obj, list):
        return [convert_tuples_to_lists(item) for item in obj]
    else:
        return obj

# Convert the dictionary for JSON serialization
json_serializable_dict = convert_tuples_to_lists(branching_ratios_dict)

# Save to JSON file
output_filename = "branching_ratios_14MeV.json"
with open(output_filename, 'w') as f:
    json.dump(json_serializable_dict, f, indent=2)

# Save failed matches to a separate JSON file
failed_matches_filename = "failed_isomer_matches_14MeV.json"
with open(failed_matches_filename, 'w') as f:
    json.dump(failed_matches, f, indent=2)

# Print summary of failed matches
print("\n" + "=" * 80)
print("ENERGY MATCHING SUMMARY")
print("=" * 80)
print(f"\nTotal failed matches: {len(failed_matches)}")

if failed_matches:
    # Categorize by failure type
    no_match = [f for f in failed_matches if f['match_quality'] == 'no_match']
    no_decay = [f for f in failed_matches if f['match_quality'] == 'no_decay_file']

    print(f"  - Energy mismatch (>1000 eV tolerance): {len(no_match)}")
    print(f"  - Missing decay file: {len(no_decay)}")

    print("\nDetailed Failed Matches:")
    print("-" * 80)
    for idx, f in enumerate(failed_matches, 1):
        za_str = f"({f['parent_z']}, {f['parent_a']})"
        zap = int(f['zap']) if f['zap'] else 0
        prod_z = zap // 1000
        prod_a = zap % 1000

        print(f"\n{idx}. Parent {za_str} MT={f['mt']} LFS={f['lfs']} → Product ({prod_z}, {prod_a})")
        print(f"   ELFS: {f['elfs_ev']:.1f} eV")
        print(f"   Assigned to: {f['matched_isomer']}")

        if f['match_quality'] == 'no_match':
            elis_str = f"{f['elis_ev']:.1f}" if f['elis_ev'] is not None else "None"
            diff_str = f"{f['energy_diff_ev']:.1f}" if f['energy_diff_ev'] is not None else "None"
            print(f"   Closest ELIS: {elis_str} eV")
            print(f"   Energy difference: {diff_str} eV (exceeds 1000 eV tolerance)")
        elif f['match_quality'] == 'no_decay_file':
            print(f"   Reason: No decay file available for product nuclide")
else:
    print("\nAll LFS states successfully matched to isomers!")

print("\n" + "=" * 80)

print(f"\nFinal branching ratios dictionary:")
pprint.pprint(branching_ratios_dict)

# Print parent nuclide statistics
print("\n" + "=" * 80)
print("PARENT NUCLIDE SUMMARY")
print("=" * 80)


ground_state_parents = 0
isomeric_parents = 0
metastable_breakdown = {}

for (z, a, m) in branching_ratios_dict.keys():
    if m == 0:
        ground_state_parents += 1
    else:
        isomeric_parents += 1
        metastable_label = f"m{m}"
        metastable_breakdown[metastable_label] = metastable_breakdown.get(metastable_label, 0) + 1

print(f"\nTotal parent nuclides: {len(branching_ratios_dict)}")
print(f"  Ground state (m=0): {ground_state_parents}")
print(f"  Isomeric states (m>0): {isomeric_parents}")

if metastable_breakdown:
    print(f"\nIsomeric state breakdown:")
    for state in sorted(metastable_breakdown.keys(), key=lambda x: int(x[1:])):
        print(f"  {state}: {metastable_breakdown[state]} parent nuclides")

print("\n" + "=" * 80)

print(f"\nBranching ratios dictionary saved to: {output_filename}")
print(f"Failed isomer matches saved to: {failed_matches_filename}")