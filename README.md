# VSP2AVL
Create AVL geometry files using OpenVSP's DegenGeom export files

## Introduction
VSP2AVL is intended to provide a fast way to convert OpenVSP (Open Vehicle Sketch Pad) aircraft geometry into an AVL (Athena Vortex Lattice) file compatible with AVL aerodynamic analysis. This allows for rapid iteration on aircraft design using OpenVSP and AVL together.

In its current state, VSP2AVL can convert basic lifting and control surfaces into AVL geometry format. VSP2AVL *cannot* accurately convert wings which have been rotated about the z-axis past around 45 degrees.

VSP2AVL has limited capability to convert bodies into AVL, and—at the moment—it can only convert them into AVL *bodies*, not cruciform lifting surfaces, which is more recommended.

## Usage
1. With a .vsp3 file open in OpenVSP, export the model geometry using `Analysis > DegenGeom > .csv > Execute`. 

2. Export airfoil geometry data into the same folder as the `.csv` using `File > Export > Airfoil Points (.dat) > OK`. This airfoil coordinate data will be modified and pointed to in the `.avl` file.

3. In `VSP2AVL.py`, define `filepath` as the file path to DegenGeom `.csv` file.

4. Add relevant AVL reference information: `Sref`, `Cref`, `Bref`, `Xref`, `Yref`, `Zref` and `mach_number`. Specify model geometry tolerance as `tolerance`. Units of tolerance are the same as units in `DegenGeom.csv`. If non-lifting-surface bodies are desired to be converted to AVL, set `write_bodies` to `True`. Set lifting-surface vortex resolution in with `vortices_per_unit_length`. Running `VSP2AVL.py` will create a new file of name `filename.avl` which can be used with AVL. 

## Current Work
* Bug testing body modeling 

* Implementing cruciform lifting surface body modeling

* Reading airfoil data from DegenGeom instead of requiring airfoil geometry export from VSP

## Future Steps
Once the program is working as intended and has been thoroughly tested, I would like to get this implemented in with OpenVSP to allow for even easier geometry conversion.
