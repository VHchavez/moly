"""
Handles functions that require quantities on the grid
"""

import numpy as np 

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

    for i in range(len(indsa0)):
        i = orb_a_info[indsa0[0]][1]
        h = orb_a_info[indsa0[0]][2]
        labelsa.append(str(i+1)+"-"+ct.gamma(h).symbol())

    for i in range(len(indsb0)):
        i = orb_b_info[indsb0[0]][1]
        h = orb_b_info[indsb0[0]][2]
        labelsb.append(str(i+1)+"-"+ct.gamma(h).symbol())
    
    return indsa0, indsb0, labelsa, labelsb
    