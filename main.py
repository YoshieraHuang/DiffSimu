import lattice as LTTC
import singleXtal as SX
import polyXtal as PX
import simu1d as S1D
import simu2d as S2D
import xray as XR
import detector as DTCT
from Vec import Vector

# l = LTTC.Lattice(material = 'Si')# lattice
# # Single Crystal
# sx = SX.SingleXtal(l, z = (0,0,1), x = (1,0,0))
# sx.Calc_rcp_space((5,5,5))
# sx.Save_rcp_space('rcp.dat')
# # Projection
# sx.Project_rcp_space(proj_vec = (1,1,1), vecx = (1,-1,0), origin = (0,0,0))
# sx.Save_proj_rcp_space('rcp_proj.dat')

# Poly Crystal
## Create Lattice
#l = LTTC.Lattice(structure = 'bcc', args = (3.304, 'Ta'))
### Create lattice parameters
x = 3.304
alpha = 90
l = LTTC.Lattice()
ltts_p = LTTC.LatticeParameter((x,x,x,alpha,alpha,alpha))
l.Add_latticeparameters(ltts_p)
### Create atoms and add them to lattice
l.Add_atom(l.Atoms_from_structure('bcc','Ta'))
# or
#l.Add_atom(LTTC.Atom('Ta', (0,0,0)))
#l.Add_atom(LTTC.Atom('Ta', (0.5,0.5,0.5)))
## Create PolyCrystal
px = PX.PolyXtal(l)
# ## Calculate RCP Space
hkls = (LTTC.Familyindex((1,0,0)), LTTC.Familyindex((0,0,2)), LTTC.Familyindex((1,0,1)))
px.Calc_rcp_space(hkls = hkls)
px.Save_proj_rcp_space('rcp.dat', res_gamma = 1)
px.Calc_rcp_space((5,5,5))
px.Save_proj_rcp_space('rcp3.dat', res_gamma = 1)
# px.Calc_rcp_space((8,8,8))
# px.Save_simple_rcp_space('rcp.dat')
# px.Save_rcp_space('rcp2.dat', res_tth =2, res_gamma = 2)
## Do Diffraction
# xr = XR.Xray(wavelength = 1.54)
# p = S1D.Profile1D(px = px, xray = xr)
# p.Calc(range_2th = (30,90), precision = 1)
# p.show()
# ## Project to a detector
# inc = Vector(0,0,-1)
# p.Add_geometry(inc = inc, x = Vector(1,0,0))
# det1 = DTCT.Detector(normal = inc , size = (400, 400), poni = ('c','c'), dist = 100)
# # det1.rotate_by(axis = (1,1,0), degree = 60)
# # det1.set(poni = (-50,'c'))
# det1.Calc_peaks(p)
# det1.show()