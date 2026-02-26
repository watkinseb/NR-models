# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 09:01:26 2025

@author: Erik Watkins

"""

import numpy as np
from refl1d.names import Parameter, SLD, Slab, FitProblem, load4, Experiment, Mixture

## =============================================================
path='C:/Users/esw/Documents/BIOPHYSICS/SALB/Decane/reduced/Oct25/ReScale/'


Lipid = 'DPhPC'
Run = 4

SysErr = 0.00
NBox =2 

Symmetric =False
SameThickness = False
SamePhi = False
SameRough= True
SameSLD= False
FixScale = True

# =============================================================================
# DEFINE PARAMETER STARTING VALUES AND RANGES
# =============================================================================

TotalThickness = 50

AverageSLD = 2.0
MinSLD = -0.4
MaxSLD = 6.56

AveragePhi = 1.0
MinPhi = 0.0
MaxPhi = 1.0

AverageRough = 5.0
MinRough = 3.0
MaxRough = 10.0

ThicknessValue = np.full(NBox, TotalThickness/NBox)
ThicknessRange = [[6,40] for _ in range(NBox)]

RoughValue = np.full(NBox, AverageRough)
RoughRange = [[MinRough, MaxRough] for _ in range(NBox)]

PhiValue = np.full(NBox, AveragePhi)
PhiRange = [[MinPhi, MaxPhi] for _ in range(NBox)]

SLDValue = np.full(NBox, AverageSLD)
SLDRange = [[MinSLD, MaxSLD] for _ in range(NBox)]


## =========================================================================

if Lipid == 'DPhPC':
    if Run == 0:
        #DPhPC 
        probe_d2o = load4(path+'215447_rescale.txt', back_reflectivity=True, name='D2O')
        probe_h2o = load4(path+'215450_rescale.txt', back_reflectivity=True, name='H2O')
        probe_smw = load4(path+'215453_rescale.txt', back_reflectivity=True, name='SMW')
    
    if Run == 1:
        #DPhPC dDecane 5:1
        probe_d2o = load4(path+'215460_rescale.txt', back_reflectivity=True, name='D2O')
        probe_h2o = load4(path+'215463_rescale.txt', back_reflectivity=True, name='H2O')
        probe_smw = load4(path+'215466_rescale.txt', back_reflectivity=True, name='SMW')
        
    if Run ==2:
        # #DPhPC dDecane 2:1
        probe_d2o = load4(path+'215549_rescale.txt', back_reflectivity=True, name='D2O')
        probe_h2o = load4(path+'215552_rescale.txt', back_reflectivity=True, name='H2O')
        probe_smw = load4(path+'215555_rescale.txt', back_reflectivity=True, name='SMW')
     
    if Run == 3:
        #DPhPC dDecane 1:1
        probe_d2o = load4(path+'215469_rescale.txt', back_reflectivity=True, name='D2O')
        probe_h2o = load4(path+'215472_rescale.txt', back_reflectivity=True, name='H2O')
        probe_smw = load4(path+'215475_rescale.txt', back_reflectivity=True, name='SMW')
        
    if Run == 4:
        #DPhPC dDecane 1:5
        probe_d2o = load4(path+'215478_rescale.txt', back_reflectivity=True, name='D2O')
        probe_h2o = load4(path+'215481_rescale.txt', back_reflectivity=True, name='H2O')
        probe_smw = load4(path+'215484_rescale.txt', back_reflectivity=True, name='SMW')


if Lipid == 'eggPC':
    if Run == 0:
        #eggPC 
        probe_d2o = load4(path+'213975_rescale.txt', back_reflectivity=True, name='D2O')
        probe_h2o = load4(path+'213978_rescale.txt', back_reflectivity=True, name='H2O')
        probe_smw = load4(path+'213981_rescale.txt', back_reflectivity=True, name='SMW')

    if Run == 1:
        #eggPC dDecane 5:1
        probe_d2o = load4(path+'215672_rescale.txt', back_reflectivity=True, name='D2O')
        probe_h2o = load4(path+'215675_rescale.txt', back_reflectivity=True, name='H2O')
        probe_smw = load4(path+'215678_rescale.txt', back_reflectivity=True, name='SMW')

    if Run == 2:
        #eggPC dDecane 1:1
        probe_d2o = load4(path+'215083_rescale.txt', back_reflectivity=True, name='D2O')
        probe_h2o = load4(path+'215086_rescale.txt', back_reflectivity=True, name='H2O')
        probe_smw = load4(path+'215089_rescale.txt', back_reflectivity=True, name='SMW')
        
    if Run == 3:
        # eggPC dDecane 1:5. 
        # fits with asymetric 2 box, not sym 5box
        probe_d2o = load4(path+'215095_rescale.txt', back_reflectivity=True, name='D2O')
        probe_h2o = load4(path+'215107_rescale.txt', back_reflectivity=True, name='H2O')
        probe_smw = load4(path+'215110_rescale.txt', back_reflectivity=True, name='SMW')


if Lipid == 'DOPC':
    if Run == 0:
        #DOPC 
        probe_d2o = load4(path+'215540_rescale.txt', back_reflectivity=True, name='D2O')
        probe_h2o = load4(path+'215543_rescale.txt', back_reflectivity=True, name='H2O')
        probe_smw = load4(path+'215546_rescale.txt', back_reflectivity=True, name='SMW')

    if Run == 1:
        #DOPC dDecane 5:1
        probe_d2o = load4(path+'215681_rescale.txt', back_reflectivity=True, name='D2O')
        probe_h2o = load4(path+'215684_rescale.txt', back_reflectivity=True, name='H2O')
        probe_smw = load4(path+'215687_rescale.txt', back_reflectivity=True, name='SMW')

    if Run == 2:
        #DOPC dDecane 1:1
        probe_d2o = load4(path+'215690_rescale.txt', back_reflectivity=True, name='D2O')
        probe_h2o = load4(path+'215693_rescale.txt', back_reflectivity=True, name='H2O')
        probe_smw = load4(path+'215696_rescale.txt', back_reflectivity=True, name='SMW')

        
    if Run == 3:
        # DOPC dDecane 1:5. 
        probe_d2o = load4(path+'215095_rescale.txt', back_reflectivity=True, name='D2O')
        probe_h2o = load4(path+'215107_rescale.txt', back_reflectivity=True, name='H2O')
        probe_smw = load4(path+'215110_rescale.txt', back_reflectivity=True, name='SMW')





# =============================================================================
# LOAD DATA
# =============================================================================

#clean up nans
for p in [probe_d2o, probe_h2o, probe_smw]:
    Q, R, dR, dQ = p.Q, p.R, p.dR, p.dQ
    mask = np.isfinite(Q) & np.isfinite(R) & np.isfinite(dR) & np.isfinite(dQ)
    p.Q = Q[mask]
    p.R = R[mask]
    p.dR = dR[mask]
    p.dQ = dQ[mask]
    p.dR = p.dR + p.R*SysErr
    
# Backgrounds (one for each model)
probe_d2o.background.range(0, 1e-5)
probe_h2o.background.range(0, 1e-5)
probe_smw.background.range(0, 1e-5)

if not FixScale:
    #Intensity 
    probe_d2o.intensity.range(0.9,1.1)
    probe_h2o.intensity.range(0.8,1.2)
    probe_smw.intensity.range(0.8,1.2)
        
# =============================================================================
# MATERIALS
# =============================================================================
D2O_SLD = Parameter(name='D2O SLD', value=6.3).range(5.3,6.36)
D2O = SLD(name='D2O', rho=D2O_SLD, irho=0.0)
H2O = SLD(name='H2O', rho=-0.56, irho=0.0)
SMW = SLD(name='SMW', rho=(0.38*D2O_SLD) + (0.62*-0.56), irho=0.0)
                 
SiOx = SLD(name='siox', rho=3.5, irho=0.0000)
Si = SLD(name='silicon', rho=2.07, irho=0.0000)     


# =============================================================================
# SUBSTRATE STRUCTURE
# =============================================================================
SiOx_thickness = Parameter(name='SiOx thickness', value=15).range(10, 25)
SiOx_rough = Parameter(name ='SiOx roughness', value=3).range(2, 5)
Si_rough = Parameter(name ='Si roughness', value=3)#.range(2, 5)
SiOx_rough = Si_rough

SiOx_layer = Slab(material=SiOx, thickness=SiOx_thickness, interface=SiOx_rough)
Si_substrate = Slab(material=Si, thickness=0.0000, interface=Si_rough)


# =============================================================================
# FUNCTION TO BUILD N PARAMETERS
# =============================================================================

def make_N_pars(N, Name, Value, Range):
    """
    Inputs are number of layers, N, Par name and initial values and ranges
    Creates roughness parameters for each of the N boxes
    """
    Par = []

    for n in range(N):
        pName = Name+'_'+str(n+1)
        v = Parameter(name=pName, value=Value[n]).range(Range[n][0], Range[n][1])
        Par.append(v)

    return Par

def make_N_SLDs(N, Name, Value, Range):
    sld = []
  
    for n in range(N):
        pName = Name+'_'+str(n+1)
        p = Parameter(name=pName, value=Value[n]).range(Range[n][0], Range[n][1])
        v = SLD(name=pName, rho=p)
        sld.append(v)

    return sld


def sum_params(seq):
    s = 0*seq[0]  # start with a Parameter-like zero
    for x in seq:
        s = s + x
    return s

# =============================================================================
# BUILD ALL EXPERIMENTS
# =============================================================================
## === Critical edge sampling ===
probe_d2o.critical_edge(substrate=Si, surface=D2O)


Rough = make_N_pars(NBox, 'Roughness', RoughValue, RoughRange)
Phi = make_N_pars(NBox, 'Volume Fraction', PhiValue, PhiRange)
Thickness = make_N_pars(NBox, 'Thickness', ThicknessValue, ThicknessRange)
SLD = make_N_SLDs(NBox, 'SLD', SLDValue, SLDRange)



if SameThickness:
    for i in range(1,NBox):
        Thickness[i]=Thickness[0]
        
if SamePhi:
    for i in range(1,NBox):
        Phi[i]=Phi[0]

if SameRough:
    for i in range(0,NBox):
        Rough[i]=Si_rough
        
if SameSLD:
    for i in range(1,NBox):
        SLD[i]=SLD[0]            
        
if Symmetric:
    if NBox % 2 == 1: 
        # odd Nbox
        for i in range((NBox//2)):
            j = -(i + 1)
            SLD[j].rho=SLD[i].rho
            Thickness[j]=Thickness[i]
            Phi[j]=Phi[i]
            Rough[j-1]=Rough[i]

    else:
        # even Nbox, so two center boxes are the same 
        for i in range((NBox//2)):
            j = -(i + 1)
            SLD[j].rho=SLD[i].rho
            Thickness[j]=Thickness[i]
            Phi[j]=Phi[i]
            Rough[j-1]=Rough[i]
            
            
SLD_layers_D2O = []
SLD_layers_H2O = []
SLD_layers_SMW = []

Layers_D2O = []
Layers_H2O = []
Layers_SMW = []

NBox_D = Si_substrate | SiOx_layer
NBox_H = Si_substrate | SiOx_layer
NBox_S = Si_substrate | SiOx_layer
for n in range(NBox):

    #calculate the SLD in each layer
    d = Mixture.byvolume(D2O, SLD[n], Phi[n]*100)
    h = Mixture.byvolume(H2O, SLD[n], Phi[n]*100)
    s = Mixture.byvolume(SMW, SLD[n], Phi[n]*100)
    
    
    SLD_layers_D2O.append(d)
    SLD_layers_H2O.append(h)   
    SLD_layers_SMW.append(s)   
    
    # layers
    lay_d = Slab(d, Thickness[n], Rough[n])
    lay_h = Slab(h, Thickness[n], Rough[n])
    lay_s = Slab(s, Thickness[n], Rough[n])
    
    NBox_D |= lay_d
    NBox_H |= lay_h
    NBox_S |= lay_s
    

NBox_D |= D2O
NBox_H |= H2O
NBox_S |= SMW

expt_d2o = Experiment(probe=probe_d2o, sample=NBox_D, name='D2O')
expt_h2o = Experiment(probe=probe_h2o, sample=NBox_H, name='H2O')
expt_smw = Experiment(probe=probe_smw, sample=NBox_S, name='SMW')



step = True
STEPSIZE=0.25

problem = FitProblem([expt_d2o, expt_h2o, expt_smw])

