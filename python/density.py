import em1d
import numpy as np

# Custom density profile
def custom_n0(x):
    return 1.0 + 0.5*np.sin(x/np.pi)*np.sin(x/np.pi)

# Time step
dt = 0.019
tmax = 10.0

# Simulation box
nx  = 64
box = 20.0

# Diagnostic frequency
ndump = 100

# Background plasma
#density = em1d.Density( type = em1d.Density.uniform )
#density = em1d.Density( type = em1d.Density.step, start = 17.5 )
#density = em1d.Density( type = em1d.Density.slab, start = 17.5, end = 22.5 )
#density = em1d.Density( type = em1d.Density.linear_ramp, start = 17.5, end = 22.5, ramp = [1.0,2.0] )
density = em1d.Density( type = em1d.Density.custom, custom = custom_n0 )

print( "initializing electrons ")

electrons = em1d.Species( "electrons", -1.0, 128, density = density )

# Diagnostics
def rep( sim ):
    if ( sim.n % ndump == 0 ):

        # Charge density
        electrons.report( em1d.Species.charge )

        # raw particle data
        electrons.report( em1d.Species.particles )

# Initialize simulation data
print( "initializing sim ")
sim = em1d.Simulation( nx, box, dt, electrons, report = rep )

# Set moving window
sim.set_moving_window()

# Run simulation
sim.run( tmax )
