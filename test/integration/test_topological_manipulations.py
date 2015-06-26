import OpenPNM
import pytest
import scipy as sp
from OpenPNM.Utilities import topology
ctrl = OpenPNM.Base.Controller()
topo = topology()


def test_subdivide():
    pn = OpenPNM.Network.Cubic(shape=[5, 5, 5],
                               spacing=0.001,
                               name='micro_net')
    pn['pore.micro'] = True
    nano_pores = [2, 13, 14, 15]
    pn.subdivide(pores=nano_pores, shape=[4, 4, 4], labels='nano')
    assert pn.Np == (125+4*64-4)
    assert pn.Nt == (300+(4*144)-16+15*16+16)
    ctrl.export(network=pn, filename='nano')


def test_clone_and_trim():
    ctrl.clear()
    pn = OpenPNM.Network.Cubic(shape=[5, 5, 5], name='net')
    geom = OpenPNM.Geometry.GenericGeometry(network=pn, name='geo1')
    geom.set_locations(pores=pn.Ps, throats=pn.Ts)
    assert sorted(list(ctrl.keys())) == ['geo1', 'net']
    pn2 = ctrl.clone_simulation(pn, name='clone')
    assert sorted(list(ctrl.keys())) == ['geo1', 'geo1_clone', 'net',
                                         'net_clone']
    topo.trim(network=pn2, pores=pn2.pores('top'))


def test_trim_extend():
    pn = OpenPNM.Network.Cubic(shape=[5, 5, 5])
    assert sp.all(sp.in1d(pn.find_neighbor_pores(pores=0), [1, 5, 25]))
    assert [pn.Np, pn.Nt] == [125, 300]
    pn.trim(pores=[0])
    assert sp.all(sp.in1d(pn.find_neighbor_pores(pores=0), [1, 5, 25]))
    assert [pn.Np, pn.Nt] == [124, 297]
    pn.extend(pore_coords=[0, 0, 0], throat_conns=[[124, 0]])
    assert [pn.Np, pn.Nt] == [125, 298]
    assert sp.all(sp.in1d(pn.find_neighbor_pores(pores=0), [1, 5, 25, 124]))


def test_stitch():
    ctrl = OpenPNM.Base.Controller()
    [Nx, Ny, Nz] = [10, 10, 10]
    pn = OpenPNM.Network.Cubic(shape=[Nx, Ny, Nz])
    pn2 = OpenPNM.Network.Cubic(shape=[Nx, Ny, Nz])
    pn2['pore.coords'][:, 2] += Nz
    pn.stitch(donor=pn2,
              P_network=pn.pores('top'),
              P_donor=pn2.pores('bottom'),
              len_max=1,
              method='nearest')
    assert pn.Np == 2*pn2.Np  # Ensure number of pores doubled
    assert pn.Nt == (2*pn2.Nt + Nx*Ny)  # Ensure correct number of new throats
    assert pn2 not in ctrl.values()  # Donor Network is removed from Controller
    # Reuse the donor Network in another stitch
    pn2['pore.coords'][:, 2] -= 2*Nz
    pn.stitch(donor=pn2,
              P_network=pn.pores('bottom'),
              P_donor=pn2.pores('top'),
              len_max=1,
              method='nearest')
    assert pn.Np == 3*pn2.Np  # Ensure number of pores increased again
    assert pn.Nt == (3*pn2.Nt + 2*Nx*Ny)  # Ensure correct num of new throats
    ctrl.clear()