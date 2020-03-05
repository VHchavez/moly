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
Volumes from Cube Files    

### Supports:
xyz files  
Psi4 geometries  
QCElemental molecules  

### Installation 
* git:
    ```
    git clone https://github.com/VHchavez/moly.git
    cd moly
    pip install .
     ```
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
*Define a figure and add molecules to it.*
 
 ```
 import moly
 fig = moly.Figure()
 molecule = moly.Molecule.from_file("caffeine.xyz")
 fig.add_molecule("caffeine", molecule)
 fig.show()
 ```
 
<br>
  
* ### Basic Layering
*Geometries can be brought from different sources and be thrown in the same figure*
 ```
#Molecules from QCArchive

import qcportal as ptl
client = ptl.FractalClient()

#Get molecule from QCArchive
ds = client.get_collection("ReactionDataset", "S22")
dimers = ds.get_molecules()
ammonia_dimer = dimers.loc['Ammonia Dimer', 'molecule'][0]

#Different surfaces are available. 
#Resolution can be increased if saving figure is desired. 
fig = moly.Figure(figsize=(800,800), surface="shiny")
fig.add_molecule("dimer", ammonia_dimer)
fig.add_molecule("bucky ball", moly.Molecue.from_file("bucky.xyz"))
fig.show()
 ```

 

#### Copyright
Copyright (c) 2020, VH Chavez


##### Acknowledgements
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
