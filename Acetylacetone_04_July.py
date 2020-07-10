#!/usr/bin/env python
# coding: utf-8

# In[265]:


from simtk.openmm.app import *
from openforcefield.utils import *


# In[266]:


help (get_data_file_path)


# In[267]:


pdb_file_path = get_data_file_path ('Acetylacetone.pdb')
pdbfile = PDBFile (pdb_file_path)

# I got the pdb file from
# https://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/P2D
# renamed it, then I manually moved it to the data directory of openforcefield.


# In[268]:


from openforcefield.topology import Molecule
acetylacetone = Molecule.from_smiles('O=C(C)CC(=O)C')
print(acetylacetone)

# I got the SMILES-code from:
# https://en.wikipedia.org/wiki/Acetylacetone


# In[269]:


print (pdb_file_path)


# In[270]:


from openforcefield.topology import *
from openforcefield.typing.engines.smirnoff import ForceField

omm_topology =  pdbfile.topology
off_topology = Topology.from_openmm(omm_topology, unique_molecules=[acetylacetone])

forcefield = ForceField('openff-1.0.0.offxml')

system = forcefield.create_openmm_system(off_topology)


# In[ ]:





# In[271]:


from simtk import *

# Langevin Dynamics:
time_step = 2*unit.femtoseconds
temperature = 300*unit.kelvin
friction = 1/unit.picosecond # collision rate
integrator = openmm.LangevinIntegrator(temperature, friction, time_step)

# Lenghth of simulation:
num_steps = 1

trj_freq = 1 # number of steps per written trajectory frame
data_freq = 1 #number of steps per fritten simulation statistics

# Set up OpenMM simulation:
simulation = openmm.app.Simulation (omm_topology, system, integrator)

# initial positions:
positions = pdbfile.getPositions()
simulation.context.setPositions(positions)

# Randomize the velocities from a Boltzmann distribution at a given temperature.
simulation.context.setVelocitiesToTemperature(temperature)

# MINIMIZE ENERGY: I deleted this line to get the non-minimized values:
simulation.minimizeEnergy()
minimized_coords = simulation.context.getState(getPositions=True).getPositions()
minimized_forces = simulation.context.getState(getForces=True).getForces()

# Configure information in output files:
pdb_reporter = openmm.app.PDBReporter('trajectory.pdb', trj_freq)
state_data_reporter = openmm.app.StateDataReporter ('data.csv', data_freq, step=True, potentialEnergy=True, temperature=True, density=True)
simulation.reporters.append (pdb_reporter)
simulation.reporters.append (state_data_reporter)


# In[ ]:





# In[272]:


import time

print("Starting simulation")
start = time.process_time()

# Run the simulation
simulation.step(num_steps)

end = time.process_time()
print("Elapsed time %.2f seconds" % (end-start))
print("Hurra!")


# In[273]:


ff_applied_parameters = forcefield.label_molecules(off_topology)[0]
ff_values=[]
ff_valuefile = open ('ff_valuefile.txt', 'w+')

for atoms, bonds in ff_applied_parameters['Bonds'].items():
    ff_valuefile.write (F'{atoms},{bonds}')
    ff_valuefile.write ('\n')
    
ff_valuefile.close()


# In[274]:


import numpy

ff_valuefile = open ('ff_valuefile.txt', 'r')

ff_valuetab = numpy.genfromtxt(fname=ff_valuefile, dtype='unicode', delimiter='Lenght')

ff_valuefile.close

man_atom_a = [0,1,1,2,2,2,3,3,3,4,4,6,6,6]
man_atom_b = [1,2,3,7,8,9,4,10,11,5,6,12,13,14]


# In[ ]:





# In[275]:


import numpy
import os

mincoorfile = os.path.join('trajectory.pdb')

mincoor = numpy.genfromtxt(fname=mincoorfile, skip_header=2, skip_footer=2, dtype='unicode')

print (mincoor)


# In[276]:


xyz_data = mincoor[:,6:9]
xyz_data = xyz_data.astype(numpy.float)
print (xyz_data)

xa = xyz_data[man_atom_a,0]
xb = xyz_data[man_atom_b,0]
ya = xyz_data[man_atom_a,1]
yb = xyz_data[man_atom_b,1]
za = xyz_data[man_atom_a,2]
zb = xyz_data[man_atom_b,2]


# In[277]:


bondlength_min=numpy.sqrt(((xa-xb)**2)+((ya-yb)**2)+((za-zb)**2))

print (bondlength_min)    


# In[ ]:




