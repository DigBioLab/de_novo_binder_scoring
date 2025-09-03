import os
import sys
import json
from pymol import cmd, stored
import re

# Define non-polar residues
NON_POLAR_RESIDUES = ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO', 'GLY']

def calculate_sasa(pdb_file, chain):
    """Calculate SASA for a given chain in a PDB file."""
    cmd.delete("all")
    cmd.load(pdb_file, "complex")
    cmd.set("dot_solvent", 1)
    sasa = cmd.get_area(f"complex and chain {chain}")
    cmd.delete("all")
    return sasa

def find_interface(pdb_file, cA='A', cB='B', cutoff=1.0):
    """
    Find interface residues between chains A and B, count non-polar residues,
    and calculate the percentage of non-polar residues in the interface.
    """
    cmd.delete("all")
    cmd.load(pdb_file, "complex")
    cmd.set("dot_solvent", 1)
    # Calculate total area for the complex and then for each chain
    cmd.get_area("complex", load_b=1)
    cmd.alter("complex", "q=b")
    cmd.extract("chA", f"complex and chain {cA}")
    cmd.extract("chB", f"complex and chain {cB}")
    cmd.get_area("chA", load_b=1)
    cmd.get_area("chB", load_b=1)
    cmd.alter("chA or chB", "b=b-q")
    stored.residues = []
    cmd.iterate("chA or chB", "stored.residues.append((model, chain,resi, resn, b))")
    
    unique_interface = {}
    nonpolar_count = 0
    binder_res=[]
    target_res=[]
    for model,chain,resi, resn, diff in stored.residues:
        if abs(diff) >= cutoff:
            unique_interface[(resi, resn)]=resn
            if resn.upper() in NON_POLAR_RESIDUES:
                nonpolar_count += 1
            if chain==cA:
                binder_res.append(f"{resn}_{resi}")
            if chain==cB:
                target_res.append(f"{resn}_{resi}")
                
    #print(unique_interface)
    interface = list(unique_interface.items())
    nonpolar_count = sum(1 for (_, resn) in interface if resn.upper() in NON_POLAR_RESIDUES)
    num_interface = len(interface)
    percent_nonpolar = (nonpolar_count / num_interface * 100) if num_interface else 0
    
    #drop duplicates
    residues_target=list(set(target_res))
    residues_binder=list(set(binder_res))

    residues_target=":".join(target_res)
    residues_binder=":".join(binder_res)

    return interface, num_interface, nonpolar_count, percent_nonpolar,residues_target,residues_binder


def get_secondary_structure(pdb_file):
    """
    Load a PDB file into PyMOL, assign secondary structure using DSSP, and return a dictionary 
    mapping each residue (by chain and residue number) to its secondary structure code.
    
    :param pdb_file: Path to the PDB file.
    :return: Dictionary with keys as (chain, resi) tuples and values as secondary structure codes.
    """
    
    # Load the PDB file into an object called "structure"
    cmd.delete("all")
    cmd.load(pdb_file, "structure")
    
    # Run DSSP to assign secondary structure
    cmd.dss("structure")
    
    # Create a dictionary of secondary structure assignments using only the CA atoms.
    # (CA atoms are typically used because they have DSSP-assigned secondary structure.)
    
    
    res_ss_dict = {"A":{},"B":{}}

    modelA = cmd.get_model("structure and (name CA and chain A)")
    modelB = cmd.get_model("structure and (name CA and chain B)")
    
    for atom in modelA.atom:
        # Use a key that combines chain and residue number.
        key = f"{atom.resn}_{atom.resi}"
        # Retrieve the secondary structure code (e.g., 'H' for helix, 'E' for sheet).
        # DSSP may leave coil regions as an empty string, so we convert those to "C".
        ss = (atom.ss or "").strip()
        if not ss:
            ss = "C"
        res_ss_dict["A"][key] = ss

    for atom in modelB.atom:
        # Use a key that combines chain and residue number.
        key = f"{atom.resn}_{atom.resi}"
        # Retrieve the secondary structure code (e.g., 'H' for helix, 'E' for sheet).
        # DSSP may leave coil regions as an empty string, so we convert those to "C".
        ss = (atom.ss or "").strip()
        if not ss:
            ss = "C"
        res_ss_dict["B"][key] = ss
    
    return res_ss_dict

def compute_ss_percentages(residues_str, chain_ss_dict):
    """
    Given a colon-separated string of residue IDs (that match keys in chain_ss_dict)
    and the secondary structure dictionary for that chain,
    compute the percentages of helix, sheet, and loop.
    
    We assume that 'H' (or h) indicates helix, 'E' (or e) indicates sheet,
    and any other secondary structure code is considered loop.
    """
    # Split the string on colon; filter out any empty strings.
    res_ids = [res.strip() for res in residues_str.split(":") if res.strip()]
    total = len(res_ids)
    helix_count = 0
    sheet_count = 0
    loop_count = 0

    
    for res in res_ids:
        # Look up the secondary structure for this residue if available.
        ss = chain_ss_dict.get(res, None)
        #print(set(val for val in chain_ss_dict.values()))
        if ss:
            ss = ss.upper()
            if ss in ["H"]:
                helix_count += 1 # CHECK THIS
            elif ss in ["S"]:
                sheet_count += 1
            elif ss in ["L"]:
                loop_count += 1
            else:
                print(f"{ss} not found in H,L,S possible codes")    
        else:
            print(f"{res} not found in dict")
            # Optionally, if the residue isn't found, you could decide to count it as loop
            # or simply ignore it. Here we ignore missing residues.
            continue
    

    #print(f"length of interface residues {len(res_ids)}, helix {helix_count}, {sheet_count}, {loop_count} sum {helix_count+sheet_count+loop_count}")
    if not helix_count+sheet_count+loop_count ==total:
        print(f"{helix_count}+{sheet_count}+{loop_count} TOTAL {total}")
        raise ValueError

    if total > 0:
        helix_percent = 100.0 * helix_count / total
        sheet_percent = 100.0 * sheet_count / total
        loop_percent  = 100.0 * loop_count / total
    else:
        helix_percent = sheet_percent = loop_percent = 0.0
    #print(f"{helix_percent} {sheet_percent} {loop_percent} total {helix_percent+loop_percent+sheet_percent}")
    if not round(helix_percent+loop_percent+sheet_percent,1)==100.0:
        if total>0:
            print(f"The secondary structure percentages do not sum to 100 percent even though there is a res string {residues_str} and chainSS dict {chain_ss_dict}")
            raise ValueError
        else: 
            print("There are no interface residues by the pymol SASA change estimation")

    return helix_percent, sheet_percent, loop_percent

def count_hydrogen_bonds(pdb_file):
    """Count hydrogen bonds between chain A and chain B."""
    cmd.delete("all")
    cmd.load(pdb_file, "complex")
    cmd.h_add("complex")
    cmd.select("donors_A", "(elem n or elem o) and (neighbor hydro) and chain A")
    cmd.select("acceptors_A", "(elem o or (elem n and not neighbor hydro)) and chain A")
    cmd.select("donors_B", "(elem n or elem o) and (neighbor hydro) and chain B")
    cmd.select("acceptors_B", "(elem o or (elem n and not neighbor hydro)) and chain B")
    hbonds_AB = cmd.find_pairs("donors_A", "acceptors_B", cutoff=3.2, angle=45)
    hbonds_BA = cmd.find_pairs("donors_B", "acceptors_A", cutoff=3.2, angle=45)
    total_hbonds = len(hbonds_AB) + len(hbonds_BA)
    return total_hbonds

def process_pdb_directory(pdb_dir, prefix, file_handle):
    """
    Process all PDB files in a directory:
      - Calculate SASA for binder (chain A) and target (chain B)
      - Compute delta SASA (unbound target SASA - bound target SASA)
      - Identify interface residues, count non-polar ones and compute their percentage
      - Count hydrogen bonds
      - Calculate the secondary structure percentages for the interface residues:
          For binder (chain A): reported as paratope_helix, paratope_sheet, paratope_loop.
          For target (chain B): reported as epitope_helix, epitope_sheet, epitope_loop.
    Write the metrics for each model to the CSV file.
    """
    for pdb in os.listdir(pdb_dir):
        if pdb.endswith(".pdb"):
            pdb_path = os.path.join(pdb_dir, pdb)
            # Use only the pdb file's basename as the model id.
            binder_id = os.path.splitext(pdb)[0]
            
            
            # Bound state: calculate SASA for binder (chain A) and target (chain B)
            sasa_binder = calculate_sasa(pdb_path, 'A')
            sasa_target = calculate_sasa(pdb_path, 'B')
            #print(sasa_binder)
            #print(sasa_target)
            
            # Compute delta SASA: difference between unbound target SASA and bound target SASA
            cmd.delete("all")
            cmd.load(pdb_path, "target_only")
            # Remove binder (chain A) to leave unbound target (chain B)
            cmd.remove("target_only and chain A")
            unbound_sasa_target = cmd.get_area("target_only and chain B")
            delta_sasa = unbound_sasa_target - sasa_target
            
            # Get interface metrics:
            # interface, num_interface, nonpolar_count, percent_nonpolar,
            # residues_target, residues_binder are returned.
            interface, num_interface, nonpolar_count, percent_nonpolar, residues_target, residues_binder = find_interface(pdb_path, cA='A', cB='B')
            #print(interface)
       
            residues_binder = str(residues_binder)
            residues_target = str(residues_target)

            # Count hydrogen bonds between target and binder
            num_hbonds = count_hydrogen_bonds(pdb_path)

            # Get the secondary structure dictionary.
            ss_dict = get_secondary_structure(pdb_path)
            # ss_dict is assumed to be of the format:
            # {"A": {THRELETTER_RESNUM: SS_code, ...}, "B": {THRELETTER_RESNUM: SS_code, ...}}
            

            # For our calculations, chain A is binder (paratope) and chain B is target (epitope).
            binder_ss_dict = ss_dict.get("A", {})
            target_ss_dict = ss_dict.get("B", {})

            #print(binder_ss_dict)
            #print(target_ss_dict)

            binder_ss_string=":".join([f"{key}={val}" for key,val in binder_ss_dict.items()])
            target_ss_string=":".join([f"{key}={val}" for key,val in target_ss_dict.items()])



            # Calculate secondary structure percentages for the interface residues.
            paratope_helix, paratope_sheet, paratope_loop = compute_ss_percentages(residues_binder, binder_ss_dict)
            epitope_helix, epitope_sheet, epitope_loop = compute_ss_percentages(residues_target, target_ss_dict)
            
            # Write out the metrics using the prefix in the column names.
            # The order of values is: model, {prefix}_SASA_binder, {prefix}_SASA_target, {prefix}_delta_SASA,
            # {prefix}_interface_residues, {prefix}_nonpolar_count, {prefix}_percent_nonpolar, {prefix}_hbonds,
            # {prefix}_res_target, {prefix}_res_binder, {prefix}_ss_dict,
            # paratope_helix, paratope_sheet, paratope_loop, epitope_helix, epitope_sheet, epitope_loop.
            file_handle.write(f"{binder_id},{sasa_binder},{sasa_target},{delta_sasa},{num_interface},{nonpolar_count},"
                              f"{percent_nonpolar},{num_hbonds},{residues_target},{residues_binder},{binder_ss_string},{target_ss_string},"
                              f"{paratope_helix},{paratope_sheet},{paratope_loop},{epitope_helix},{epitope_sheet},{epitope_loop}\n")

# Read the JSON file containing the dictionary mapping prefixes to PDB directories.
with open("pymol_files/pdb_dirs.json") as f:
    pdb_dict = json.load(f)
print(f"pdb_dict found: {pdb_dict}")

# For each prefix, create a separate CSV file with header columns incorporating the prefix.
for prefix, pdb_dir in pdb_dict.items():
    output_file = f"pymol_files/pymol_metrics_{prefix}.csv"
    header = (
        f"binder_id,{prefix}_pymol_SASA_binder,{prefix}_pymol_SASA_target,{prefix}_pymol_delta_SASA,"
        f"{prefix}_pymol_num_interface_residues,{prefix}_pymol_nonpolar_count,{prefix}_pymol_percent_nonpolar,{prefix}_pymol_hbonds,"
        f"{prefix}_pymol_intf_res_target,{prefix}_pymol_intf_res_binder,{prefix}_pymol_binder_ss,{prefix}_pymol_target_ss,"
        f"{prefix}_pymol_paratope_helix,{prefix}_pymol_paratope_sheet,{prefix}_pymol_paratope_loop,{prefix}_pymol_epitope_helix,{prefix}_pymol_epitope_sheet,{prefix}_pymol_epitope_loop\n"
    )
    with open(output_file, "w") as f:
        f.write(header)
        process_pdb_directory(pdb_dir, prefix, f)
    print(f"Metrics saved to {output_file}")