<p align="center">
<br>
<img src="media/title.png" alt="moly" height=300> <br><br>
Molecular visualization in Jupyter.<br><br>
<a href="https://travis-ci.com/VHChavez/moly"><img src="https://travis-ci.com/VHChavez/moly.svg?branch=master" /></a>  
<a href="https://lgtm.com/projects/g/VHchavez/moly/context:python"><img src="https://img.shields.io/lgtm/grade/python/g/VHchavez/moly.svg?logo=lgtm&logoWidth=18" /></a>  
<a href="https://opensource.org/licenses/BSD-3-Clause"><img src="https://img.shields.io/badge/License-BSD%203--Clause-blue.svg" /></a>
<br>
</p>

---

<br>
(Package is under development!)


### Features:  
Geometry  
Volumes (Density, MOs, ESP)    

### Supports:
xyz files  
Psi4 geometries  
QCElemental molecules  
Cube files  

### Installation 
* pip:
    ```
    pip install moly
     ```
* conda:
    ```
    coming soon!
    ```
      
<br>

 
* ### Basic Geometry
*Define a figure and add elements to it. These can be molecues:*
 
 ```
 import moly
 caffeine = moly.molecule_factory("xyz", file='/caffeine.xyz')
 fig = moly.Figure()
 fig.add_molecule(caffeine)
 fig.show()
 ```
 
#### Produces
<img src="/media/caffeine.png" alt="caffeine" height=300> <br>

<br>

* ### Basic Cube file 
*Geometry and volumentric information can extracted from cube files.*
 ```
formal = moly.molecule_factory("Cube", file='orbitals.cube')
fig = moly.Figure()
fig.add_molecule(formal)
fig.add_blob(iso=0.05)
fig.show()
 ```
 #### Produces
  <img src="/media/formaldehyde.png" alt="formal" height=300> <br>
  
  
* ### Basic Layering
*Geometries can be brought from different sources and be thrown in the same figure*
 ```
#Implementation with QCArchive

import qcportal as ptl
client = ptl.FractalClient()

#Get molecule from QCArchive
ds = client.get_collection("ReactionDataset", "S22")
dimers = ds.get_molecules()
ammonia_dimer = dimers.loc['Ammonia Dimer', 'molecule'][0]

#Molecules
ammonia_dimer = moly.molecule_factory("QC", molecule=ammonia_dimer)
bucky = moly.molecule_factory("xyz", file='bucky.xyz')

#Figure, resolution and surface material can be changed
fig = moly.Figure(resolution=(800,800), surface="shiny")
fig.add_molecule(ammonia_dimer)
fig.add_molecule(bucky)
fig.show()
 ```
 
  #### Produces
  <img src="/media/bucky.png" alt="bucky" height=300> <br>
 

 

#### Copyright
Copyright (c) 2020, VH Chavez


##### Acknowledgements
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
