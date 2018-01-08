import em1d
import numpy as np

nx = 120
box = 4 * np.pi
dt = 0.1
tmax = 50.0

ndump = 10

ppc = 500
ufl = [0.2,    0.0,  0.0]
uth = [0.001,0.001,0.001]

right = em1d.Species( "right", -1.0, ppc, ufl = ufl, uth = uth )

ufl[0] = -ufl[0]
left  = em1d.Species( "left", -1.0, ppc, ufl = ufl, uth = uth )

def rep( sim ):
    if ( sim.n % ndump == 0 ):
        right.report(em1d.SpeciesDiag.particles)
        left.report(em1d.SpeciesDiag.particles)
        sim.emf.report(em1d.EMF.efld,0)

sim = em1d.Simulation( nx, box, dt, [right,left], report = rep )

sim.run(tmax)

