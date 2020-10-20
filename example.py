import numpy as np
import matplotlib as plt

import openmc
import openmc.mgxs as mgxs

temperatures = [922, 972, 1022, 1072, 1122]

#define the fuel material
msre_fuel = openmc.Material(1, "molten-salt-fuel", temperature[0])
msre_fuel.add_nuclide('Li7', 0.1090, 'wo')
msre_fuel.add_nuclide('Li6', 5e-6, 'wo')
msre_fuel.add_nuclide('F19', 0.6680, 'wo')
msre_fuel.add_nuclide('Be9', 0.0627, 'wo')
msre_fuel.add_nuclide('U235', 0.0167, 'wo')
msre_fuel.add_nuclide('U238', 0.0344, 'wo')
msre_fuel.set_density('g/cm3', 2.146) #see Moltres paper

#define the moderator material
msre_moderator = openmc.Material(2, "moderator", temperature[0])
msre_moderator.add_nuclide('C12', 1, 'wo')
msre_moderator.set_density('g/cm3', 1.86) #see Moltres paper
#msre_moderator.add_s_alpha_beta('c_Graphite')

#export the materials in an xml file
msre_materials = openmc.Materials([msre_fuel, msre_moderator])
msre_materials.export_to_xml()


#define the geometry
R = 72.5 #[cm], radius of the reactor core
H = 151.75 #[cm], height of the reactor core
num_segments = 14
pitch = R / num_segments
x = 0.2379522
fuel_rad = x * pitch
mod_rad = pitch - fuel_rad
dw = 0.1 #[cm], OpenMC can't construct 2D geometry so we just make a very thin region

#bounds
bounds = [0, 0, 0, R, 0, H]
top = openmc.ZPlane(z0=H, boundary_type='vacuum')
bottom = openmc.ZPlane(z0=0, boundary_type='vacuum')
inside = openmc.XPlane(x0=0, boundary_type='vacuum')
outside = openmc.XPlane(x0=R, boundary_type='vacuum')
slicea = openmc.YPlane(y0=-dw/2, boundary_type='reflective')
sliceb = openmc.YPlane(y0=+dw/2, boundary_type='reflective')
core = +slicea & -sliceb & -top & +bottom & +inside & -outside

#set up regions for fuel and moderator
fuel_coord_inner = 0
fuel_coord_outer = fuel_rad
fuel_channel_inner = openmc.XPlane(x0 = fuel_coord_inner)
fuel_channel_outer = openmc.XPlane(x0 = fuel_coord_outer)
fuel_channel = +fuel_channel_inner & -fuel_channel_outer & core
fuel_region = fuel_channel 

mod_coord_inner = fuel_rad
mod_coord_outer = fuel_rad + 2 * mod_rad

mod_channel_inner = openmc.XPlane(mod_coord_inner)
mod_channel_outer = openmc.XPlane(mod_coord_outer)
mod_channel = +mod_channel_inner & -mod_channel_outer & core
moderator_region = mod_channel



for i in range(1,13):
    fuel_coord_inner = mod_coord_outer
    fuel_coord_outer = fuel_coord_inner + 2*fuel_rad
    fuel_channel_inner = openmc.XPlane(x0 = fuel_coord_inner)
    fuel_channel_outer = openmc.XPlane(x0 = fuel_coord_outer)
    fuel_channel = +fuel_channel_inner & -fuel_channel_outer & core
    fuel_region = fuel_region | fuel_channel 

    mod_coord_inner = fuel_coord_outer
    mod_coord_outer = fuel_coord_outer +  2 * mod_rad

    mod_channel_inner = openmc.XPlane(mod_coord_inner)
    mod_channel_outer = openmc.XPlane(mod_coord_outer)
    mod_channel = +mod_channel_inner & -mod_channel_outer & core
    moderator_region = moderator_region | mod_channel


fuel_cell = openmc.Cell(1)
fuel_cell.region = fuel_region
fuel_cell.fill = msre_fuel

moderator_cell = openmc.Cell(2)
moderator_cell.region = moderator_region
moderator_cell.fill = msre_moderator

root = openmc.Universe(cells=(fuel_cell, moderator_cell))
geom = openmc.Geometry(root)
geom.export_to_xml()

#following the MGXS example on the OpenMC website from here on out
#define the simulation settings
batches = 50
inactie = 10
partices = 2500

settings_file = openmc.Settings()
settings_file.batches = batches
settings_file.inactitve = inactive
settings_file.particles = particles
settings_file.output = {'tallies' : True}

uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
settings_file.source = openmc.Source(space=uniform_dist)

settings.export_to_xml()


#Create our energy groups
neutron_groups = mgxs.EnergyGroups()
neutron_groups.group_edges = np.array([0., 0.625, 20.0e6]) #2 neutron groups 

delayed_groups = list(range(1,7))



# Instantiate a few different sections
chi_prompt = mgxs.Chi(domain=cell, groups=energy_groups, by_nuclide=True, prompt=True)
prompt_nu_fission = mgxs.FissionXS(domain=cell, groups=energy_groups, by_nuclide=True, nu=True, prompt=True)
chi_delayed = mgxs.ChiDelayed(domain=cell, energy_groups=energy_groups, by_nuclide=True)
delayed_nu_fission = mgxs.DelayedNuFissionXS(domain=cell, energy_groups=energy_groups, delayed_groups=delayed_groups, by_nuclide=True)
beta = mgxs.Beta(domain=cell, energy_groups=energy_groups, delayed_groups=delayed_groups, by_nuclide=True)
decay_rate = mgxs.DecayRate(domain=cell, energy_groups=one_group, delayed_groups=delayed_groups, by_nuclide=True)

chi_prompt.nuclides = ['U235', 'Pu239']
prompt_nu_fission.nuclides = ['U235', 'Pu239']
chi_delayed.nuclides = ['U235', 'Pu239']
delayed_nu_fission.nuclides = ['U235', 'Pu239']
beta.nuclides = ['U235', 'Pu239']
decay_rate.nuclides = ['U235', 'Pu239']


