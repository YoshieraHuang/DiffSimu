# DiffSimu

`DiffSimu` is a diffraction simulation for polycrystal and single crystal

# Dependencies

## Public Python modules
* `numpy`[w](http://www.numpy.org/)
* `scipy`[w](http://www.scipy.org)
* `matplotlib`[w](https://matplotlib.org)
* `multiprocessing` for multi processing
* `logging`
* `itertools`
* `re` support for regex
* `functools`
* `pyquaternion` for rotation calculation

## Private Python modules
(may be included in this project)
* `Vec` in which class `Vector` supports vecter operations

## Modules

There are 8 main modules and 2 auxiliary modules.

8 main modules:
1. `lattice` : *abbr.* `LTTC`. contains several attributes involved in **Crystallography** like
    * `LatticeParameter()`
    * `Element()`
    * `Atom()`
    * `index()`
    * `Familyindex()`

2. `xray`: *abbr.* `XR`. contains class `Xray` which can model monochromatic and polychromatic **x-ray**
3. `peak`: *abbr.* `PK`. contains class `Peak` which can model **diffraction peak**. Both 1-D diffraction peak in Bragg diffraction and 2-D diffrraction peak in Laue diffraction can be modeled by `Peak`
4. `polyXtal`: *abbr.* `PX`. contains class `PolyXtal` which can model **polycrystal**
5. `singleXtal`: *abbr.* `SX`. contains class `SingleXtal` which can model **single crystal**
6. `simu1d`: *abbr.* `S1D`. contains class `Porfile1D` which simulates profile of **Bragg diffraction**
7. `simu2d`: *abbr.* `S2D`. contains class `Pattern2D` which simulates pattern of **Laue diffraction**
8. `detector`: *abbr.* `DET`. contains class `Detector` which models **detector**. According to the geometry of diffraction and detector, diffraction pattern on detector can be simulated

2 auxiliary modules:
1. `LUT`: "Look-Up Table". Datebase in diffraction simulation
2. `myfunctools`: *abbr.* `FT`. Useful function tools

## Class `LatticeParameter`
Lattice paramter (a, b, c, alpha , beta , gamma). Matrices that contains basis vectors in direct space and reciprocal space can calculated
### Initialization
* from lattice parameters (6-element list)
```python
import lattice as LTTC
LP = LTTC.LatticeParameters(a,b,c,alpha,beta,gamma)
```
* from direct matrix (3*3 list)
```python
LP = LTTC.LatticeParameters(DirectMatrix)
```

**NOTE**:
* angles must be in *degrees*
* if the lattice parameters are not legal, `ValueError` will be raised

### Attributes
* `a`, `b`, `c`, `alpha`, `beta`, `gamma`: `float`
    - Lattice parameters

* `direct_matrix`: `np.ndarray`, read-only
    - Matrix contains basis vectors in real space. Note that the vectors are vertical

* `reciprocal_matrix`: `np.ndarray`, read-only
    - Matrix contains basis vectors in reciprocal space

### Methods
* `show(fout = None)`: No return
    - print information in this `LatticeParameter` to `fout`
    - `fout`: an instance of `File`


## Class `FracCoor`

`FracCoor` is 3-dimensional fractional coordinates which is used to describe postions of atom in unit cell. `FracCoor` inherits from `Vector`, which allows to use `FracCoor` as 3D vector

### Initializaiton

* from 3-dimensional list, Note that if dimension of list is not 3, or even one coordinate is negative or exceeds 1, error will be raise
```python
import lattice as LTTC
FC = LTTC.FracCoor((0.5,0.5,0.5))
```

### Attributes

* `x`, `y`, `z`: `float`
    - x-, y- and z- coordinates.
    - is same as `FC[0]`, `FC[1]` and `FC[2]`


## Class `Element`

`Element` is a class for atom element.

### Initialization
* from name:
```python
import lattice as LTTC
Cu = LTTC.Element('Cu')
O = LTTC.Element('O')
```

* Input name will be regularized. Expressions like `Element('cu')` and `Element('o')` are allowed
* Then coefficients of scattering factor of this element will be searched in dictionary `LUT.dict_element` and stored into attribute `sc_factor_coeff`

### Attributes
* `name`: `str`
    - regularized name of element
* `sc_factor_coeff`

### Methods
* `sc_factor(tth, wl)`: `complex`
    - **scattering factor** when x-ray with wavelength `wl` be scattered to angle `tth`
    - `wl`: wavelength in *Angstrom*
    - `tth`: $2\theta$, scattering angle in *radian*

### Class `Atom`

`Atom` can describe atom in unit cell.

### Initialization
* from `Element` and `FracCoor`
```python
import lattice as LTTC
atom1 = LTTC.Atom('Cu', (0,0,0))
Ta = LTTC.Element('Ta')
FC = LTTC.FracCoor((0.5,0.5,0.5))
atom2 = LTTC.Atom(Cu, (0,0.5,0.5))
atom3 = LTTC.Atom('Cu', FC)
atom4 = LTTC.Atom(Cu, FC)
```

### Attributes
* `element`: `Element`
    - element of atom
* `coor`: `FracCoor`
    - coordinates of atom in unit cell

## Class `index`

`index` is Miller index of lattice plane. inherit from `Vector`

### Initialization
* from 3D integer list:
```python
import lattice as LTTC
hkl = LTTC.index((1,1,1))
```

### Attributes
* `h`, `k`, `l`: `float`
    - three integer elements of Miller index
    - same as `index()[0]`, `index()[1]` and `index()[2]`

* `str`: `string`
    - traditional representation of `index`, *hkl*

## Class `Familyindex`

`Familyindex` is a subclass of `index` and represents a family of Miller index. Negative index is not allowed in `Familyindex` 

### Initialization
* like `index`

### New Attributes
* `multiplicity`: `integer`
    - Multiplicity of this family. Number of `sons()`

### New Method
* `sons()`: list of `index`
    - sons of this family. No duplicates.

### Class `Lattice`

`Lattice` is a class of crystal lattice that contains lattice parameters and atoms

### Initialization
* from material symbol:
```python
import lattice as LTTC
Cu_lattice = LTTC.Lattice(material = 'Cu')
```
    * lattice parameters and atoms will be searched at `LUT.dict_material` and added to `Lattice`
* from structure symbol
```python
import lattice as LTTC
Cu_lattice = Lattice(structure = 'fcc', args = (3.615, 'Cu'))
```
    * structure symbol is used to search strucure informations in `LUT.dict_structure` and `args` is used to complete structure information 

### Attributes
* `LP`:
    - `LatticeParameter` of this `Lattice`
* `atoms`:
    - list of `Atom` in `Lattice`

* `reciprocal_matrix`:
    - same as `Lattice.LP.reciprocal_matrix`

### Methods
* `Add_latticeparameters(latticeparam)`:
    - assign `latticeparam` to `Lattice.LP`
    - `latticeparam`: `LatticeParameters`

* `Add_atom(atom)`:
    - add `atom` to `Lattice.atoms`
    - `atom`: `Atom`

* `strct_factor(hkls)`: Generator
    - return **structure factor**(complex) in `Lattice` for each `index` in `hkls`. If structure factor is zero, return `None`
    - For simplification, not all coefficients of scattering factor is considered. It's only recommamded to confirm whether `index` is structurally extinct or not.
    - `hkls`: list or generator of `index`

* `sc_factor(hkl, tth, wl)`:
    - return **scattering factor** in `Lattice` for `hkl`. if scattering factor is zero , return `None`
    - Whole scattering factor is considered.
    - `hkl`: `index`
    - `tth`: $2\theta$, scattering angle in *radian*
    - `wl`: wavelength in *Angstrom*

* `isextinct(hkls)`: Generator
    - for each `index` in `hkls`, if `index` is structurally extinct, return `True`, else return `False`

* `vec_in_lattice(hkls, rcp_matrix = None)`: Generator
    - for each `index` in `hkls`, return **reciprocal vectors** in `Lattice`. if `index` is `None`, return `None`
    - Reciprocal vector is calculated by `Lattice.reciprocal_matrix`. if `rcp_matrix` is assigned, `rcp_matrix` will be used to calculate reciprocal vector.
* `Vec_in_lattice(hkls, rcp_matrix = None)`:
    - return *list* of **reciprocal vectors** from generator `vec_in_lattice`
    - If `hkls` has only one `index`, corresponding reciprocal vector, not a list, will be returned

* `d_spacing(hkls)`: Generator
    - return *d_spacing* in `Lattice` for each `index` in `hkls`

* `D_spacing(hkls)`:
    - list version of `d_spacing()`

* `isperpendicular(hkl1, hkl2)`:
    - If `hkl1` and `hkl2` is perpendicular in `Lattice`, return `True`, else return `False`
    - `hkl1`, `hkl2`: `index`

* `isparallel(hkl1,hkl2)`:
    - If `hkl1` and `hkl2` is parallel in `Lattice`, return `True`, else return `False`
    - `hkl1`, `hkl2`: `index`

* `show()`:
    - print abstract view of `Lattice`

## Functions in `Lattice`
* `Gen_hklfamilies(hklrange = (10,10,10), lattice = None)`: Generator
    - Generate a series of `Familyindex`
    - `hklrange`: maximum of `index`
    - `lattice`: instance of `Lattice`
        + if `lattice` is not `None`, `Familyindex` that is structurally extinct in `lattice` will be filtered off
* `Gen_hkls(hklfamilies)`: Generator
    - Generate `index` belonging to `hklfamilies`

## Class `Xray`

`Xray` describes x-ray. Monochromatic and Polychromatic x-ray both can be modeled

### Initialization
* from wavelength or energy or intensity:
    - from wavelength in *Angstrom*.  If intensity is not given, intensity will be set to 1

```python
import xray as XR
xr1 = XR.Xray(wavelength = 0.5) # Mono

import numpy as np
xr2 = XR.Xray(wavelength = np.linspace(0.4, 0.6, 100)) #Poly

import scipy.signal
xr3 = XR.Xray(wavelength = np.linspace(0.4, 0.6, 100), intensity = scipy.signal.gaussian(100, 0.05)) #Poly and Gaussian intensity
```

* from energy in *keV*

```python
import xray as XR
xr1 = XR.Xray(energy = 20) # Mono

import numpy as np
xr2 = XR.Xray(energy = np.linspace(20, 25, 100)) # Poly
```
* from file
    - Data in file must have 2 columns, first is wavelength or energy, other intensity. If first column is in wavelength, `islambda` should be `True`.If energy, `islambda` should be `False`. `islambda` 

```python
import xray as XR
xr1 = XR.Xray(filename = 'spectrum_wl.txt', islambda = True) # wavelength
xr2 = XR.Xray(filename = 'spectrum_energy.txt', islambda = False) # energy
```
