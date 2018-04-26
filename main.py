import lattice as LTTC
import singleXtal as SX
import polyXtal as PX
import simu1d as S1D
import simu2d as S2D
import xray as XR

l = LTTC.Lattice(material = 'Si')# lattice
# Single Crystal
sx = SX.SingleXtal(l, z = (0,0,1), x = (1,0,0))
sx.Calc_rcp_space((15,15,15))
sx.Save_rcp_space('rcp.dat') 

# Poly Crystal
px = PX.PolyXtal(l)
px.Calc_rcp_space((8,8,8))
px.Save_simple_rcp_space('rcp.dat')
px.Save_rcp_space('rcp2.dat', res_tth =2, res_gamma = 2)