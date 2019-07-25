import openpnm as op
import scipy as sp
import matplotlib.pyplot as plt

# network
sp.random.seed(17)
pn = op.network.Cubic(shape=[40, 40], spacing=1e-4)

# geometry
geo = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)

# phase
gas = op.phases.Air(network=pn)

# physics
phys = op.physics.Standard(network=pn, phase=gas, geometry=geo)
gas['pore.concentration'] = 0
# physics: sink
phys['pore.sinkA'] = -1e10
phys['pore.sinkb'] = 1
mod_sink = op.models.physics.generic_source_term.power_law
phys.add_model(propname='pore.sink', A1='pore.sinkA', A2='pore.sinkb',
               X='pore.concentration', model=mod_sink)
# physics: src
phys['pore.srcA'] = 1e-10
phys['pore.srcb'] = 1
mod_src = op.models.physics.generic_source_term.power_law
phys.add_model(propname='pore.source', A1='pore.srcA', A2='pore.srcb',
               X='pore.concentration', model=mod_src)

# algorithm
rx = op.algorithms.FickianDiffusion(network=pn)
rx.setup(phase=gas)

# comment the following line on and off to see problem
rx.set_source(propname='pore.source', pores=[410, 610, 810])
rx.set_source(propname='pore.sink', pores=[430, 630, 830])

rx.set_value_BC(values=1, pores=pn.pores('front'))

rx.run()

# plot
im = sp.reshape(rx['pore.concentration'], [40, 40]).T
fig = plt.figure(figsize=[10, 10])
plt.imshow(im)
plt.colorbar()
