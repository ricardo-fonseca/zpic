import em1d
import numpy as np

# Time step
dt = 0.019
tmax = 22.8

# Simulation box
nx  = 1000
box = 20.0

# Diagnostic frequency
ndump = 50

# Background plasma

ppc = 128 # Particles per cell

electrons = em1d.Species( "electrons", -1.0, ppc,
	                      density = em1d.Density( type = em1d.Density.step, start = 20.0))

# Initialize simulation data
sim = em1d.Simulation( nx, box, dt, tmax, electrons )

# Add laser pulse
sim.add_laser( em1d.Laser( start = 17.0, fwhm = 2.0, a0 = 2.0, omega0 = 10.0, polarization = np.pi/2 ))

# Set moving window
sim.set_moving_window()

# Set current smoothing
sim.set_smooth( em1d.Smooth(xtype = em1d.Smooth.compensated, xlevel = 4) )

# Simulation reports
def report( sim ):
    if ( sim.n % ndump == 0 ):
        # All electric field components
        sim.emf.report( em1d.EMF.efld, 0 )
        sim.emf.report( em1d.EMF.efld, 1 )
        sim.emf.report( em1d.EMF.efld, 2 )

        # Charge density
        electrons.report( em1d.Species.charge )

        # x1u1 phasespace
        electrons.report( em1d.phasespace( em1d.Species.x1, em1d.Species.u1),
                          pha_nx = [1024,512], pha_range = [[0.0,20.0],[-2.0,2.0]])

# Run the simulation
sim.run(report)
