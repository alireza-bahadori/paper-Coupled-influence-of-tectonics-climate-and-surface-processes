#!/usr/bin/env python
# coding: utf-8

# # 4-D simulation of Bahadori et al. (2022): Coupled influence of tectonics, climate, and surface processes on landscape evolution in southwestern North America

# In[1]:


import UWGeodynamics as GEO
import numpy as np
import math
import underworld.function as fn
from UWGeodynamics import non_dimensionalise as nd
import os 


# ## Define variables for model setup

# In[3]:


# this scale coeff mulifplies the length and with of the domain by 10 and we will have 100km for each one degree.
scale_coef = 10

#  define spacing for making your mesh in x and y.
n = 5
#  define spacing for making your mesh in z direction in km.
nz = 4

# spacing for making 3D structures like topo, crust, etc.
n_2 = 10

# spacing for making badlands model. 
n_3 = 1

dfs = int(n)
dfs_2 = int(n_2)
dfs_3 = int(n_3)
nz = int(nz)

# maximum_elevation_of_3D_model_from_sealevel km
elsl = 6
# upper crust basal depth km
ucd = -10.0
# middle crust basal depth km
mcd = -25.0
# mantle_asthenosphere_thickness_for_3D_model km
mat = 70

# air and sticky air density coefficient.
factor = 1000

# define laterally vaying temperature at the base.
heat_add = (mat)*5.0


# ## Mesh resolution

# In[12]:


resx = 42
resy = 34
resz = 49
Res_mesh = resx, resy, resz
model_top = int(elsl)


# In[6]:


print("grid resolution is: resx: {}, resy: {}, resz: {}".format(resx,resy,resz))


# # Working with units
# 
# Note that this is not an obligation and you can use values without units 
# 
# 
# The geodynamics module enables usage of a broad range of units using a *UnitRegistry*. You can attach a unit to any value that you need to define. A value with a units attached becomes a *Quantity* python object. The geodynamics module take care of the conversion internally so you may use any units that you think are appropriate. You can also mix them.
# 
# The module will also try to work out a scaling of the values to help the computation process. The user can chose to alter the way his or her values are scaled or can rely on the default options.
# 
# To use the units system, you can link the unit registry as follow:

# In[7]:


u = GEO.UnitRegistry


# ## Scaling

# In[8]:


model_length  = 1. * u.kilometer
surfaceTemp   = 273.15 * u.degK # 0 * u.degC
baseModelTemp = 1633.15 * u.degK # 1360 * u.degC
bodyforce = (3300 * u.kilogram / u.metre**3 * 9.81 * u.meter / u.second**2) 

KL = model_length
Kt = 1. * u.year
KM = bodyforce * KL**2 * Kt**2
KT = (baseModelTemp - surfaceTemp)

GEO.scaling_coefficients["[length]"] = KL
GEO.scaling_coefficients["[time]"] = Kt
GEO.scaling_coefficients["[mass]"]= KM
GEO.scaling_coefficients["[temperature]"] = KT


# # Model setup
# 
# The starting model is 2100 x 1700 km long and ~200 km thick.

# In[9]:


Model = GEO.Model(elementRes=Res_mesh,
                  minCoord=(0.0 * scale_coef * u.kilometer, 0.0 * scale_coef * u.kilometer, (-1*((resz*nz)-model_top)) * u.kilometer),  
                  maxCoord=(resx*n * scale_coef * u.kilometer, resy*n * scale_coef * u.kilometer, model_top * u.kilometer), 
                  gravity=(0.0, 0.0, -9.81 * u.meter / u.second**2))

print(Model.mesh.data)
# print(len(Model.mesh.data))


# # Output directory for geodynamic model

# In[10]:


Model.outputDir="Output_WUS_3D"


# # Creating a synthetic 3D lithosphere structure
# 
# Note that this simulation starts at step #5 which is 35.5Ma. The initial 3-D lithosphere structure has been already produced and provided in directory "Output_WUS_3D". To restart the simulations from 35.5Ma we just create a synthetic lithosphere structure here.

# In[ ]:


air = Model.add_material(name="Air", shape=GEO.shapes.Layer3D(top=Model.top, bottom=5.0 * u.kilometer))
water = Model.add_material(name="Water", shape=GEO.shapes.Layer3D(top=4.0 * u.kilometer, bottom=3.0 * u.kilometer))
ml = Model.add_material(name="MantleLithosphere", shape=GEO.shapes.Layer3D(top=-60.0 * u.kilometer, bottom=-120.0 * u.kilometer))
lc = Model.add_material(name="LowerCrust", shape=GEO.shapes.Layer3D(top=mcd * u.kilometer, bottom=-60.0 * u.kilometer))
mc = Model.add_material(name="MiddleCrust", shape=GEO.shapes.Layer3D(top=ucd * u.kilometer, bottom=mcd * u.kilometer))
uc = Model.add_material(name="UpperCrust", shape=GEO.shapes.Layer3D(top=-3.0 * u.kilometer, bottom=ucd * u.kilometer))
oc = Model.add_material(name="OceanicCruct", shape=GEO.shapes.Layer3D(top=2.0 * u.kilometer, bottom=1.0 * u.kilometer))
topo = Model.add_material(name="Topography", shape=GEO.shapes.Layer3D(top=-2.0 * u.kilometer, bottom=-3.0 * u.kilometer))
sair = Model.add_material(name="StickyAir", shape=GEO.shapes.Layer3D(top=air.bottom, bottom=4.0 * u.kilometer))
ma = Model.add_material(name="MantleAsthenosphere", shape=GEO.shapes.Layer3D(top=ml.bottom, bottom=Model.bottom))
sediment = Model.add_material(name="Sediment")
wz = Model.add_material(name="Trench", shape=GEO.shapes.Layer3D(top=water.bottom, bottom=2.0 * u.kilometer))
lc_CP = Model.add_material(name="LowerCrustCP", shape=GEO.shapes.Layer3D(top=oc.bottom, bottom=0.0 * u.kilometer))
mc_CP = Model.add_material(name="MiddleCrustCP", shape=GEO.shapes.Layer3D(top=lc_CP.bottom, bottom=-1.0 * u.kilometer))
ml_o = Model.add_material(name="MantleLithosphereOceanic", shape=GEO.shapes.Layer3D(top=mc_CP.bottom, bottom=-2.0 * u.kilometer))


# # Define model parameters

# In[65]:


air.density = 1*factor * u.kilogram / u.metre**3
sair.density = 1*factor * u.kilogram / u.metre**3

water.density = 1*factor * u.kilogram / u.metre**3

# thermalExpansivity: r_tE
r_tE = 3e-5 / u.kelvin
# reference_temperature: r_tem
r_tem = 298 * u.kelvin

ml.density  = GEO.LinearDensity(3330. * u.kilogram / u.metre**3,thermalExpansivity = r_tE,
                                reference_temperature = r_tem)   

uc.density = GEO.LinearDensity(2700. * u.kilogram / u.metre**3,thermalExpansivity = r_tE,
                                reference_temperature = r_tem) 

mc.density = GEO.LinearDensity(2900. * u.kilogram / u.metre**3,thermalExpansivity = r_tE,
                                reference_temperature = r_tem) 

lc.density = GEO.LinearDensity(2900. * u.kilogram / u.metre**3,thermalExpansivity = r_tE,
                                reference_temperature = r_tem) 

topo.density = GEO.LinearDensity(2500. * u.kilogram / u.metre**3,thermalExpansivity = r_tE,
                                reference_temperature = r_tem) 

ml_o.density  = GEO.LinearDensity(3380. * u.kilogram / u.metre**3,thermalExpansivity = r_tE,
                                reference_temperature = r_tem)

oc.density = GEO.LinearDensity(2750. * u.kilogram / u.metre**3,thermalExpansivity = r_tE,
                                reference_temperature = r_tem)

wz.density = GEO.LinearDensity(2700. * u.kilogram / u.metre**3,thermalExpansivity = r_tE,
                                reference_temperature = r_tem)

lc_CP.density = GEO.LinearDensity(2900. * u.kilogram / u.metre**3,thermalExpansivity = r_tE,
                                reference_temperature = r_tem)

mc_CP.density = GEO.LinearDensity(2900. * u.kilogram / u.metre**3,thermalExpansivity = r_tE,
                                reference_temperature = r_tem)

sediment.density = 2300. * u.kilogram / u.metre**3

upper_mantle_density = 3400.0 * u.kilogram / u.metre**3

ma.density  = GEO.LinearDensity(upper_mantle_density,thermalExpansivity = r_tE,
                                reference_temperature = r_tem)     

air.radiogenicHeatProd = 0. * u.microwatt / u.meter**3
sair.radiogenicHeatProd = 0. * u.microwatt / u.meter**3
ma.radiogenicHeatProd  = 0.022 * u.microwatt / u.meter**3  
ml.radiogenicHeatProd  = 0.022 * u.microwatt / u.meter**3  
uc.radiogenicHeatProd = 0.9 * u.microwatt / u.meter**3  
mc.radiogenicHeatProd = 0.9 * u.microwatt / u.meter**3  
lc.radiogenicHeatProd = 0.9 * u.microwatt / u.meter**3  
topo.radiogenicHeatProd = 0.5 * u.microwatt / u.meter**3  
sediment.radiogenicHeatProd = 0.6 * u.microwatt / u.meter**3
water.radiogenicHeatProd = 0.5 * u.microwatt / u.meter**3
ml_o.radiogenicHeatProd  = 0.022 * u.microwatt / u.meter**3
oc.radiogenicHeatProd = 0.5 * u.microwatt / u.meter**3
wz.radiogenicHeatProd = 0.5 * u.microwatt / u.meter**3
lc_CP.radiogenicHeatProd = 0.9 * u.microwatt / u.meter**3
mc_CP.radiogenicHeatProd = 0.9 * u.microwatt / u.meter**3

Model.diffusivity = 1.0e-6 * u.metre**2 / u.second 
Model.capacity    = 1000. * u.joule / (u.kelvin * u.kilogram)

air.diffusivity = 1.0e-6 * u.metre**2 / u.second
sair.diffusivity = 1.0e-6 * u.metre**2 / u.second
thermal_diffusivity_upper_crust = 8.3e-7 * u.metre**2 / u.second 
thermal_diffusivity_lower_crust = 6.7e-7 * u.metre**2 / u.second 
thermal_diffusivity_mantle = 1.0e-6 * u.metre**2 / u.second 

ml.diffusivity  = thermal_diffusivity_mantle
uc.diffusivity = thermal_diffusivity_upper_crust
mc.diffusivity = thermal_diffusivity_lower_crust
lc.diffusivity = thermal_diffusivity_lower_crust
topo.diffusivity = thermal_diffusivity_upper_crust
sediment.diffusivity = 1.0e-6 * u.metre**2 / u.second  
ma.diffusivity  = thermal_diffusivity_mantle
ml_o.diffusivity  = thermal_diffusivity_mantle
water.diffusivity = thermal_diffusivity_mantle
oc.diffusivity = thermal_diffusivity_mantle
wz.diffusivity = thermal_diffusivity_mantle
lc_CP.diffusivity = thermal_diffusivity_lower_crust
mc_CP.diffusivity = thermal_diffusivity_lower_crust

air.capacity = 100. * u.joule / (u.kelvin * u.kilogram)
sair.capacity = 100. * u.joule / (u.kelvin * u.kilogram)
ml.capacity  = 1000. * u.joule / (u.kelvin * u.kilogram)
uc.capacity = 1000. * u.joule / (u.kelvin * u.kilogram) 
mc.capacity = 1000. * u.joule / (u.kelvin * u.kilogram) 
lc.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)
topo.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)
sediment.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)
ma.capacity  = 1000. * u.joule / (u.kelvin * u.kilogram)
water.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)
ml_o.capacity  = 1000. * u.joule / (u.kelvin * u.kilogram)
oc.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)
wz.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)
lc_CP.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)
mc_CP.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)


# # Thermal Boundary Conditions
# 
# Absolute temperatures at the top and bottom of the model are defined here.

# In[66]:


import pandas as pd
import numpy as np

with open("UWG_temperature/UWG_temp_35.5_Ma.txt", "r") as input_crust:  
    crust = pd.read_csv(input_crust, dtype=float, delimiter=",",skiprows=0)
    df_crust = pd.DataFrame(crust)
    df_crust = np.array(df_crust)

    x = df_crust[:,0]
    y = df_crust[:,1]
    z = df_crust[:,2]
    f = len(x)
    q1 = x.reshape(f,1)
    q2 = y.reshape(f,1)
    q3 = z.reshape(f,1)
    index1 = q1
    temp = z + 274.15

    temp_data0 = []
    coor_temp = []

    for data in range(len(index1)):
        temp_data = temp[data]

        x_coor = x[data]*scale_coef
        y_coor = y[data]*scale_coef
        coor_temp = (float(x_coor),float(y_coor),float(temp_data))
        temp_data0.append(coor_temp)

    max_y = int((resy*dfs*scale_coef)+(1*dfs*scale_coef))
    
    temp_data1 = []

    for data in range(0,max_y,dfs*scale_coef):
        for node in range(len(temp_data0)):
            val1 = temp_data0[node]
            if round(val1[1]) == data:
                if round(val1[0]) % (dfs*scale_coef) == 0:
                    temp_data1.append(val1)                                   

    values1 = []
    for data in range(len(temp_data1)):
        value = temp_data1[data]
        values = value[2]
        values1.append(values)

values10 = []
for data in range(len(values1)):

    value = values1[data] + heat_add 
    values10.append(value)          

base = int((resx+1) * (resy+1))
nodes_all = []
for node in range(0,base,1):
 
    nodes_all.append(node)

temp_at_nodes = []
for data in range(len(nodes_all)):
    temp = ([nodes_all[data]], values10[data] * u.degK)
    temp_at_nodes.append(temp)


Model.set_temperatureBCs(top= 273.15 * u.degK,
                     materials=[(air, 273.15 * u.degK),
                                (sair, 273.15 * u.degK)],
                     nodeSets=(temp_at_nodes))

Model.set_heatFlowBCs(bottom=(-200 * u.microwatt / u.metre**2,ma))


# # Kinematic Boundary Conditions
# 
# We use free slip on the top, left, right, back, and front of the model.

# In[27]:


Model.set_velocityBCs(top=[None, None, 0. * u.centimeter / u.year],
                        left=[0.0 * u.centimeter / u.year, None, None],
                        right=[0.0 * u.centimeter / u.year, None, None],
                        back=[None, 0.0 * u.centimeter / u.year, None],
                        front=[None, 0.0 * u.centimeter / u.year, None])


# # Stress Boundary Conditions
# 
# We use a traction field from the mantle convection model at the base of the model.

# In[28]:


import pandas as pd
import numpy as np

with open("UWG_traction/traction_vector_uwg_35.5_Ma.txt", "r") as input_crust:
    crust = pd.read_csv(input_crust, dtype=float, delimiter=",",skiprows=0)
    df_crust = pd.DataFrame(crust)
    df_crust = np.array(df_crust)

    x = df_crust[:,0]
    y = df_crust[:,1]
    k1 = df_crust[:,2]
    k2 = df_crust[:,3]
    k3 = df_crust[:,4]
    f = len(x)
    q1 = x.reshape(f,1)
    q2 = y.reshape(f,1)
    q3 = k1.reshape(f,1)
    q4 = k2.reshape(f,1)
    q5 = k3.reshape(f,1)
    index1 = q1

    tracton_data0 = []
    coor_traction = []

    for data in range(len(index1)):
        x_coor = x[data]*scale_coef
        y_coor = y[data]*scale_coef
        coor_traction = (float(x_coor),float(y_coor),float(k1[data]),float(k2[data]),float(k3[data]*-1))
        tracton_data0.append(coor_traction)

    max_y = int((resy*dfs*scale_coef)+(1*dfs*scale_coef))

    tracton_data1 = []

    for data in range(0,max_y,dfs*scale_coef):
        for node in range(len(tracton_data0)):
            val1 = tracton_data0[node]
            if round(val1[1]) == data:
                if round(val1[0]) % (dfs*scale_coef) == 0:
                    tracton_data1.append(val1)

    values1 = []
    for data in range(len(tracton_data1)):
        value = tracton_data1[data]
        values = value[2],value[3],value[4]
        values1.append(values)
            
base = int((resx+1) * (resy+1))
nodes_all = []
for node in range(0,base,1):
    nodes_all.append(node)

traction_at_nodes = []
for data in range(len(nodes_all)):
    tract = ([nodes_all[data]], [values1[data][0] * u.pascal, values1[data][1] * u.pascal, values1[data][2] * u.pascal])
    traction_at_nodes.append(tract)

Model.set_stressBCs(top=[None, None, None], bottom=[0. * u.pascal, 0. * u.pascal, 0. * u.pascal], nodeSets=(traction_at_nodes))


# # Define Viscosities
# 
# The crust and the mantle have a visco-plastic rheology with a
# temperature and stress dependent viscosity for stresses below the yield stress,
# and a depth dependent plastic branch above it.
# We use power law relationships between strain rate and stress
# to describe dislocation creep. 
# 
# The viscosity varies with the temperature and stress. Please see Mthods section of the manuscript for more information.

# In[40]:


rh=GEO.ViscousCreepRegistry()


# In[41]:


Model.minViscosity = 1e18 * u.pascal * u.second
Model.maxViscosity = 1e24 * u.pascal * u.second

air.viscosity      = 1.0e19 * u.pascal * u.second
sair.viscosity     = 2.5e19 * u.pascal * u.second

uc.viscosity = 0.0015 * rh.Wet_Quartz_Dislocation_Gleason_and_Tullis_1995

topo.viscosity = 0.0015 * rh.Wet_Quartz_Dislocation_Gleason_and_Tullis_1995

lc.viscosity  = 10.0 * GEO.ViscousCreep(name='Wet Diorite (Carter and Tsenn, 1987)',
                                 preExponentialFactor=3.2e-2/u.megapascal ** 2.4 /u.second,
                                 stressExponent=2.4,
                                 activationVolume=0*u.meter ** 3 / u.mole,activationEnergy=212 * u.kilojoules/u.mole,
                                 f=1.0)  # Wet Diorite (Carter and Tsenn, 1987) eta _{0}=3.98e16

mc.viscosity  = 0.5 * GEO.ViscousCreep(name='Wet Diorite (Carter and Tsenn, 1987)',
                                 preExponentialFactor=3.2e-2/u.megapascal ** 2.4 /u.second,
                                 stressExponent=2.4,
                                 activationVolume=0*u.meter ** 3 / u.mole,activationEnergy=212 * u.kilojoules/u.mole,
                                 f=1.0)  # Wet Diorite (Carter and Tsenn, 1987) eta _{0}=3.98e16

ml.viscosity = 0.4 * GEO.ViscousCreep(name='Wet Olivine (Burov, 2011)',
                                 preExponentialFactor=4.8/u.megapascal ** 3.0 /u.second,
                                 stressExponent=3.0,
                                 activationVolume=0.,activationEnergy=502 * u.kilojoules/u.mole,
                                 f=1.0)  # Dry Olivine (Burov, 2011) eta _{0}=1.97e17

ma.viscosity = 3.0 * rh.Dry_Olivine_Dislocation_Karato_and_Wu_1993
     
sediment.viscosity = rh.Wet_Quartz_Dislocation_Gleason_and_Tullis_1995

water.viscosity      = 2.5e19 * u.pascal * u.second
oc.viscosity = 0.0015 * rh.Wet_Quartz_Dislocation_Gleason_and_Tullis_1995
wz.viscosity = 1.0 * GEO.ViscousCreep(name='Dry Maryland Diabase(strong) (Burov, 2011)',
                                 preExponentialFactor=8/u.megapascal ** 4.7 /u.second,
                                 stressExponent=4.7,
                                 activationVolume=0*u.meter ** 3 / u.mole,activationEnergy=485 * u.kilojoules/u.mole,
                                 f=1.0)  # Dry Maryland Diabase(strong) (Burov, 2011) eta _{0}=3.98e16

ml_o.viscosity = 0.001 * GEO.ViscousCreep(name='Dry Maryland Diabase(strong) (Burov, 2011)',
                                 preExponentialFactor=8/u.megapascal ** 4.7 /u.second,
                                 stressExponent=4.7,
                                 activationVolume=0*u.meter ** 3 / u.mole,activationEnergy=485 * u.kilojoules/u.mole,
                                 f=1.0)  # Dry Maryland Diabase(strong) (Burov, 2011) eta _{0}=3.98e16

lc_CP.viscosity  = 1.0 * GEO.ViscousCreep(name='Dry Maryland Diabase(strong) (Burov, 2011)',
                                 preExponentialFactor=8/u.megapascal ** 4.7 /u.second,
                                 stressExponent=4.7,
                                 activationVolume=0*u.meter ** 3 / u.mole,activationEnergy=485 * u.kilojoules/u.mole,
                                 f=1.0)  # Dry Maryland Diabase(strong) (Burov, 2011) eta _{0}=3.98e16 #ref: Geophysical and tectonic framework of the eastern
                                         # Basin and Range-Colorado Plateau-Rocky Mountain transition

mc_CP.viscosity = 1.5 * rh.Wet_Quartz_Dislocation_Gleason_and_Tullis_1995


# # Define Plasticity
# 
# Plastic behavior is assigned using the same approach as for viscosities.

# In[ ]:


pl = GEO.PlasticityRegistry()


# In[42]:


topo.plasticity = GEO.DruckerPrager(name="Topo",
                                                cohesion=1. * u.megapascal,
                                                cohesionAfterSoftening=0.2 * u.megapascal,
                                                frictionCoefficient=0.1,
                                                frictionAfterSoftening=0.01,
                                                epsilon1=0.0, epsilon2=0.2)

uc.plasticity = GEO.DruckerPrager(name="Upper Crust",
                                                cohesion=1. * u.megapascal,
                                                cohesionAfterSoftening=0.2 * u.megapascal,
                                                frictionCoefficient=0.1,
                                                frictionAfterSoftening=0.01,
                                                epsilon1=0.0, epsilon2=0.2)

mc.plasticity = GEO.DruckerPrager(name="Middle Crust",
                                                cohesion=10. * u.megapascal,
                                                cohesionAfterSoftening=2. * u.megapascal,
                                                frictionCoefficient=0.1,
                                                frictionAfterSoftening=0.01,
                                                epsilon1=0.0, epsilon2=0.5)

lc.plasticity = GEO.DruckerPrager(name="Lower Crust",
                                                cohesion=10. * u.megapascal,
                                                cohesionAfterSoftening=2. * u.megapascal,
                                                frictionCoefficient=0.1,
                                                frictionAfterSoftening=0.01,
                                                epsilon1=0.0, epsilon2=0.5)

ml.plasticity = GEO.DruckerPrager(name="Upper Mantle",
                                           cohesion=10. * u.megapascal,
                                           cohesionAfterSoftening=2. * u.megapascal,
                                           frictionCoefficient=0.1,
                                           frictionAfterSoftening=0.01,
                                           epsilon1=0.0, epsilon2=0.5)

ma.plasticity = GEO.DruckerPrager(name="Mantle Asthenosphere",
                                                cohesion=10. * u.megapascal,
                                                cohesionAfterSoftening=2. * u.megapascal,
                                                frictionCoefficient=0.1,
                                                frictionAfterSoftening=0.01,
                                                epsilon1=0.0, epsilon2=0.5)

sediment.plasticity = pl.Huismans_et_al_2011_Crust
sediment.plasticity.epsilon1 = 0.01
sediment.plasticity.epsilon2 = 1.0

ml_o.plasticity = GEO.DruckerPrager(name="Upper Mantle",
                                           cohesion=10. * u.megapascal,
                                           cohesionAfterSoftening=2. * u.megapascal,
                                           frictionCoefficient=0.1,
                                           frictionAfterSoftening=0.01,
                                           epsilon1=0.0, epsilon2=0.5)

oc.plasticity = GEO.DruckerPrager(name="Oceanic Crust",
                                                cohesion=5. * u.megapascal,
                                                cohesionAfterSoftening=1. * u.megapascal,
                                                frictionCoefficient=0.05,
                                                frictionAfterSoftening=0.01,
                                                epsilon1=0.0, epsilon2=0.5)

water.plasticity = GEO.DruckerPrager(name="Oceanic Crust",
                                                cohesion=5. * u.megapascal,
                                                cohesionAfterSoftening=1. * u.megapascal,
                                                frictionCoefficient=0.05,
                                                frictionAfterSoftening=0.01,
                                                epsilon1=0.0, epsilon2=0.5)

wz.plasticity = GEO.DruckerPrager(name="Weak Zone",
                                                cohesion=5. * u.megapascal,
                                                cohesionAfterSoftening=1. * u.megapascal,
                                                frictionCoefficient=0.05,
                                                frictionAfterSoftening=0.01,
                                                epsilon1=0.0, epsilon2=0.5)

lc_CP.plasticity = GEO.DruckerPrager(name="Lower Crust CP",
                                                cohesion=15. * u.megapascal,
                                                cohesionAfterSoftening=3. * u.megapascal,
                                                frictionCoefficient=0.44,
                                                frictionAfterSoftening=0.088,
                                                epsilon1=0.0, epsilon2=0.5)

mc_CP.plasticity = GEO.DruckerPrager(name="Middle Crust CP",
                                           cohesion=15. * u.megapascal,
                                           cohesionAfterSoftening=3. * u.megapascal,
                                           frictionCoefficient=0.44,
                                           frictionAfterSoftening=0.088,
                                           epsilon1=0.0, epsilon2=0.5)


# ## Melt

# In our experiments, the viscosity decreases linearly by 2 orders of magnitude when the melt fraction increases from 15 to 30%. When the melt fraction is 15%, the viscosity of the melted crust is that of the non-melted surrounding; when the melt fraction is 30%, its viscosity is a thousand times lower than in surrounding material. Please see Methods section of the manuscript for more information.

# In[44]:


solidii = GEO.SolidusRegistry()
crust_solidus = solidii.Crustal_Solidus

liquidii = GEO.LiquidusRegistry()
crust_liquidus = liquidii.Crustal_Liquidus


lc.add_melt_modifier(crust_solidus, crust_liquidus, 
                         latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
                         meltFraction=0.,
                         meltFractionLimit=0.3,
                         meltExpansion=0.13, 
                         viscosityChangeX1 = 0.15,
                         viscosityChangeX2 = 0.30,
                         viscosityChange = 1e-2
                        )  

mc.add_melt_modifier(crust_solidus, crust_liquidus, 
                         latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
                         meltFraction=0.,
                         meltFractionLimit=0.3,
                         meltExpansion=0.13, 
                         viscosityChangeX1 = 0.15,
                         viscosityChangeX2 = 0.30,
                         viscosityChange = 1e-2
                        )  


uc.add_melt_modifier(crust_solidus, crust_liquidus, 
                         latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
                         meltFraction=0.,
                         meltFractionLimit=0.3,
                         meltExpansion=0.13, 
                         viscosityChangeX1 = 0.15,
                         viscosityChangeX2 = 0.30,
                         viscosityChange = 1e-2
                        )  

topo.add_melt_modifier(crust_solidus, crust_liquidus, 
                         latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
                         meltFraction=0.,
                         meltFractionLimit=0.3,
                         meltExpansion=0.13, 
                         viscosityChangeX1 = 0.15,
                         viscosityChangeX2 = 0.30,
                         viscosityChange = 1e-2
                        ) 

wz.add_melt_modifier(crust_solidus, crust_liquidus,
                         latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
                         meltFraction=0.,
                         meltFractionLimit=0.3,
                         meltExpansion=0.13,
                         viscosityChangeX1 = 0.15,
                         viscosityChangeX2 = 0.30,
                         viscosityChange = 1e-2
                        )


# # Compute initial condition

# In[45]:


Model.init_model(defaultStrainRate=1e-14 * 1/u.second)


# In[46]:


Model.update_melt_fraction()


# ## Surface Processes Model

# In[60]:


Model.surfaceProcesses = GEO.surfaceProcesses.Badlands(airIndex=[air.index], sedimentIndex=sediment.index,
                                          XML="maped_rain/wus_test1.xml", resolution= dfs_3 * scale_coef * u.kilometer, surfElevation=0.,
                                          checkpoint_interval=100000 * u.years, restartFolder="outbdls")


# ## Solver options

# In[ ]:


Model.solver.set_inner_method("lu")

Model.solver.set_penalty(1e1)

GEO.rcParams["initial.nonlinear.max.iterations"] = 30

GEO.rcParams["nonlinear.max.iterations"] = 30

GEO.rcParams["initial.nonlinear.tolerance"]= 2e-2 

GEO.rcParams["nonlinear.tolerance"]= 2e-2

GEO.rcParams["temperature.SIunits"]= u.degC

GEO.rcParams["velocityField.SIunits"]= u.millimeter / u.year

GEO.rcParams["default.outputs"]= ["temperature","pressureField","strainRateField","velocityField","projStressField","projTimeField","projMaterialField","projViscosityField","projMeltField","projPlasticStrain","projDensityField","projStressTensor"]


# ## Run Geodynamics Model

# In[ ]:


Model.run_for(500000.0 * u.years, restartStep=-1, checkpoint_interval=100000.0*u.years,dt=100000.0*u.years, restartDir="Output_WUS_3D")

