import em1d
import numpy as np

nx = 120
box = 4 * np.pi
dt = 0.1
tmax = 50.0

ndump = 10

d = em1d.Density()
d.n = 1.0

ppc = 500
ufl = np.array([0.2,    0.0,  0.0], dtype=np.float32)
uth = np.array([0.001,0.001,0.001], dtype=np.float32)

right = em1d.Species( "right", -1.0, ppc, ufl, uth, d )

ufl[0] = -ufl[0]
left  = em1d.Species( "left", -1.0, ppc, ufl, uth, d )

sim = em1d.Simulation( nx, box, dt, tmax, [right,left])

def report( sim, n ):
    if ( n % ndump == 0 ):
        right.report(em1d.Species.particles)
        left.report(em1d.Species.particles)
        sim.emf.report(em1d.EMF.efld,0)

sim.run(report)

