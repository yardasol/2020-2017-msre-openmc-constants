import numpy as np
import matplotlib as plt

import openmc
import openmc.mgxs as mgxs

ttemperatures = [922, 972, 1022, 1072, 1122]

###################
#### MATERIALS ####
###################

# Create the fuel Material
fuel_material = openmc.Material(1, "molten_salt_fuel")
fuel_material.add_nuclide('Li7', 0.1090, 'wo')
fuel_material.add_nuclide('Li6', 5e-6, 'wo')
fuel_material.add_nuclide('F19', 0.6680, 'wo')
fuel_material.add_nuclide('Be9', 0.0627, 'wo')
fuel_material.add_nuclide('U235', 0.0167, 'wo')
fuel_material.add_nuclide('U238', 0.0344, 'wo')
fuel_material.set_density('g/cm3', 2.146) # see Moltres paper

# Create the moderator Material
moderator_material = openmc.Material(2, "moderator")
moderator_material.add_element('C', 1, 'wo')
moderator_material.set_density('g/cm3', 1.86) # see Moltres paper
moderator_material.add_s_alpha_beta('c_Graphite')

# Create the Materials and export to XML
msre_materials = openmc.Materials([fuel_material, moderator_material])
msre_materials.export_to_xml()


##################
#### GEOMETRY ####
##################

# Reactor_core parameters
R = 72.5 # [cm], radius of the reactor reactor_core
H = 151.75 # [cm], height of the reactor reactor_core
num_segments = 14 # fuel and moderator segments
pitch = R / num_segments
x = 0.2379522 # what is this???
fuel_radius = x * pitch
moderator_radius = pitch - fuel_radius

# Create the reactor core Region
reactor_core_top_surface = openmc.ZPlane(z0=H, boundary_type='vacuum')
reactor_core_bottom_surface = openmc.ZPlane(z0=0, boundary_type='vacuum')
reactor_core_left_surface = openmc.XPlane(x0=0, boundary_type='reflective')
reactor_core_right_surface = openmc.XPlane(x0=R, boundary_type='vacuum')
reactor_core_region = -reactor_core_top_surface & +reactor_core_bottom_surface & +reactor_core_left_surface & -reactor_core_right_surface


# Create the reactor core Universe
reactor_core_universe = openmc.Universe(name='reactor_core')

# Create a fuel channel Region
fuel_channel_left_surface_coord = 0
fuel_channel_right_surface_coord = fuel_radius #might need to change this to be the full diameter
fuel_channel_left_surface = openmc.XPlane(x0 = fuel_channel_left_surface_coord)
fuel_channel_right_surface = openmc.XPlane(x0 = fuel_channel_right_surface_coord)
fuel_channel_region = +fuel_channel_left_surface & -fuel_channel_right_surface & reactor_core_region

# Create a fuel channel Cell
fuel_channel_cell = openmc.Cell(name='fuel_channel_0') 
fuel_channel_cell.fill = fuel_material
fuel_channel_cell.region = fuel_channel_region

# Create a moderator Region
moderator_channel_left_surface_coord = fuel_radius
moderator_channel_right_surface_coord = fuel_radius + moderator_radius#2 * moderator_radius
moderator_channel_left_surface = openmc.XPlane(moderator_channel_left_surface_coord)
moderator_channel_right_surface = openmc.XPlane(moderator_channel_right_surface_coord)
moderator_channel_region = +moderator_channel_left_surface & -moderator_channel_right_surface & reactor_core_region

# Create a moderator Cell
moderator_channel_cell = openmc.Cell(name='moderator_channel_0') 
moderator_channel_cell.fill = moderator_material
moderator_channel_cell.region = moderator_channel_region

# Add the fuel and moderator channel Cells to the reactor core Universe
reactor_core_universe.add_cell(fuel_channel_cell)
reactor_core_universe.add_cell(moderator_channel_cell)
cells = []
cells.append(fuel_channel_cell)
cells.append(moderator_channel_cell)

# Repeat above for the remaining 13 channels
for i in range(1,14):
    fuel_channel_region = fuel_channel_region.translate([pitch, 0, 0])
    
    fuel_channel_cell = openmc.Cell(name=f'fuel_channel_{i}') 
    fuel_channel_cell.fill = fuel_material
    fuel_channel_cell.region = fuel_channel_region
    
    moderator_channel_region = moderator_channel_region.translate([pitch, 0, 0])
    
    moderator_channel_cell = openmc.Cell(name=f'moderator_channel_{i}') 
    moderator_channel_cell.fill = moderator_material
    moderator_channel_cell.region = moderator_channel_region
    
    reactor_core_universe.add_cell(fuel_channel_cell)
    reactor_core_universe.add_cell(moderator_channel_cell)
    cells.append(fuel_channel_cell)
    cells.append(moderator_channel_cell)

#Create the root Cell
root_cell = openmc.Cell(name = 'root_cell', fill=reactor_core_universe)
root_cell.region = reactor_core_region

# Create the root Universe
root_universe = openmc.Universe(universe_id = 0, name='root_universe')
root_universe.add_cell(root_cell)

# Create the Geometry and export to XML
geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()


##################
#### SETTINGS ####
##################

# OpenMC simulation parameters
batches = 50
inactive = 10
particles = 2500

# Instatiate a Settings Object
settings_file = openmc.Settings()
settings_file.batches = batches
settings_file.inactive = inactive
settings_file.particles = particles
settings_file.output = {'tallies' : True}

# Create an initial uniform spatial source distribution over fissionable zones
dist_bounds = [0, 0, 0, R, H, R]
uniform_dist = openmc.stats.Box(dist_bounds[:3], dist_bounds[3:], only_fissionable=True)
settings_file.source = openmc.Source(space=uniform_dist)

# Export Settings to XML
settings_file.export_to_xml()


#######################
#### ENERGY GROUPS ####
#######################

# Instantiate a 2-group EnergyGroups object
neutron_groups = mgxs.EnergyGroups()
neutron_groups.group_edges = np.array([0., 0.625, 20.0e6]) #are these right? 

# Initialize a 2-energy-group and 6-delayed-group MGXS Library
mgxs_lib = mgxs.Library(geometry)
mgxs_lib.energy_groups = neutron_groups
mgxs_lib.num_delayed_groups = 6


# Specify multi-group cross section types to compute
mgxs_lib.mgxs_types = ['total', 'transport', 'nu-scatter matrix', 'kappa-fission', 'inverse-velocity', 'chi-prompt',
                      'prompt-nu-fission', 'chi-delayed', 'delayed-nu-fission', 'beta']

# Specify a "cell" domain type for the cross section tally filters
mgxs_lib.domain_type = 'cell'

# Specify the cell domain over which to compute multi-group cross sections
mgxs_lib.domain = cells

# Construct all tallies needed forthe multi-group cross section library
mgxs_lib.build_library()

# Create a "tallies.xml" file for the MGXS Library
tallies_file = openmc.Tallies()
mgxs_lib.add_to_tallies_file(tallies_file, merge=True)

# Instantiate a current tally
cell_filter = openmc.CellFilter(cells)
current_tally = openmc.Tally(name='current_tally')
current_tally.scores = ['current']
current_tally.filters = [cell_filter]

# Add current tally to the tallies file
tallies_file.append(current_tally)

# Export to "tallies.xml"
tallies_file.export_to_xml()


###############
#### PLOTS ####
###############

# Instantiate a Plot object
plot = openmc.Plot()

# Plot parameters
plot.basis = 'xz'
plot.origin = (R/2, 0, H/2)
plot.width = (1.5*R, 1.5*H)
plot.pixels = (int(10*R), int(10*H))
plot.color_by = 'material'

# Create a Plots object and export to XML
plots = openmc.Plots([plot])
plots.export_to_xml()

# Plot the geometry
openmc.plot_geometry()



# Run OpenMC
openmc.run()
