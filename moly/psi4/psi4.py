"""

Functions that involve the use of psi4 objects

"""

import numpy as np
import plotly.graph_objects as go


def orbital_explorer(wfn):
    """
    
    Simple plot with orbital's energies
    
    Parameters
    ----------
    wfn : psi4.core.UHF, psi4.core.RHF
        wavefunction object from psi4 calculation
    
    """
    
    vira = wfn.epsilon_a().np
    virb = wfn.epsilon_b().np
    occa = wfn.epsilon_a().np[:wfn.nalpha()]
    occb = wfn.epsilon_b().np[:wfn.nbeta()]

    xalpha = np.zeros_like(vira) - 1.0 
    xbeta  = np.zeros_like(virb) + 1.0 
    lbound = np.zeros_like(vira) - 2.0
    rbound = np.zeros_like(virb) + 2.0
    
    fig = go.Figure()
    ball_radius = 0.01

    #VIRTUAL ALPHA ORBITAL
    fig.add_trace(go.Scattergl(x=xalpha, y=vira,
                        mode='markers',
                        name='Virtual', 
                        marker_symbol='line-ew', 
                        marker_line_width=1))

    #VIRTUAL BETA ORBITAL
    fig.add_trace(go.Scattergl(x=xbeta, y=virb,
                        mode='markers',
                        name='Virtual', 
                        marker_symbol='line-ew', 
                        marker_line_width=1))

    #OCCUPIED ALPHA ORBITAL
    fig.add_trace(go.Scattergl(x=xalpha, y=occa + ball_radius, 
                             mode='markers',
                             name='Occupied'))

    #OCCUPIED BETA ORBITAL
    fig.add_trace(go.Scattergl(x=xbeta, y=occb + ball_radius,
                            mode='markers',
                            name='Occupied '))

    fig.add_trace(go.Scattergl(x=lbound, y=vira,
                        mode='markers',
                        name='markers',
                        marker = {"size" : 0,
                                  "opacity" : 0}))

    fig.add_trace(go.Scattergl(x=rbound, y=virb,
                        mode='markers',
                        name='markers',
                        marker = {"size" : 0,
                                  "opacity" : 0}))


    fig.update_layout(template = "plotly_white",
                      showlegend = False,
                      yaxis_title = "Energy (Hartrees)",
                      xaxis = {"tickmode"   : 'array',
                               "tickvals"   : [-5.0, -1.0, 1.0, 5.0],
                               "ticktext"   : ["", "Alpha", "Beta", ""],
                               "fixedrange" : True})


    fig.show()
    
    return