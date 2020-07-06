"""
Handles functions that require quantities on the grid
"""

import numpy as np 
from opt_einsum import contract

def compute_orbital_properties(wfn, orbitals="all"):
    """
    Gets orbital information
    
    Parameters
    ----------
    wfn : psi4.core.RHF, psi4.core.UHF
        wavefunction object
    orbitals : List
        List of orbitals requested. 
        Postive integers   ->Alpha orbitals
        Negative integers  ->Beta orbitals
        
    Returns
    ---------
    indsa0, indsb0 : List
        List with indices of orbitals requested
    
    labelsa, labelsb : List
        List with orbital labels
    """
    
    nirrep = wfn.nirrep()
    nmopi = wfn.nmopi()
    Ca = wfn.Ca_subset("AO", "ALL").np
    Cb = wfn.Cb_subset("AO", "ALL").np

    orb_a_info = []
    orb_b_info = []
    nalpha = 0
    nbeta  = 0

    for h in range(nirrep):
        nalpha += wfn.nalphapi()[h]
        nbeta  += wfn.nbetapi()[h]
        for i in range(nmopi[h]):
            orb_a_info.append((wfn.epsilon_a().nph[h][i], i, h))
            orb_b_info.append((wfn.epsilon_b().nph[h][i], i, h))

    orb_a_info.sort()
    orb_b_info.sort()
    
    if orbitals is "all":
        indsa0 = [ i for i in range(len(Ca))]
        indsb0 = [ i for i in range(len(Cb))]

    if orbitals is not "all":
        indsa0 = [ i-1      for i in orbitals if i > 0]
        indsb0 = [ abs(i)-1 for i in orbitals if i < 0]
    
    labelsa = []
    labelsb = []

    ct = wfn.basisset().molecule().point_group().char_table()

    for j in range(len(indsa0)):
        i = orb_a_info[indsa0[0]][1]
        h = orb_a_info[indsa0[0]][2]
        labelsa.append(str(j+1)+"-"+ct.gamma(h).symbol())

    for j in range(len(indsb0)):
        i = orb_b_info[indsb0[0]][1]
        h = orb_b_info[indsb0[0]][2]
        labelsb.append(str(j+1)+"-"+ct.gamma(h).symbol())
    
    return indsa0, indsb0, labelsa, labelsb

def orbitals_on_grid(C, blocks, points_func):
    #C_on_grid = []

    points_func.set_pointers(C)
    C_np = C.clone().np
    nbf = len(C_np)

    orbitals_r = {str(i_orb) : [] for i_orb in range(nbf)}

    for block in blocks:
        points_func.compute_points(block)
        npoints = block.npoints()
        lpos = np.array(block.functions_local_to_global())
        w = np.array(block.w())
        phi = np.array(points_func.basis_values()["PHI"])[:npoints, :lpos.shape[0]]

        for i_orb in range(nbf):
            if len(lpos) != 0:
                #C_local = C_np[i_orb, lpos]
                C_local = C_np[lpos, i_orb]
                orb = contract('m, pm -> p', C_local, phi)
                orbitals_r[str(i_orb)].append(orb)        

        # if lpos.shape[0] != 0:
        #     C_local = C_np[(lpos[:, None], lpos)]
        #     orb = contract('nm, pm -> np', C_local.T, phi)
        #     C_on_grid.append(orb)

    return orbitals_r
