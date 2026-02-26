"""BUILD AN OILY BILAYER OUT OF VOLUME FRACTION BOXES"""

import numpy as np
from refl1d.names import Parameter, SLD, Slab, FitProblem, load4, Mixture
from molgroups.refl1d_interface import (Substrate, VolumeFractionBox,
                                        MolgroupsLayer,
                                        MolgroupsExperiment,
                                        make_samples)

## =========================================================================
path='C:/Users/esw/Documents/BIOPHYSICS/SLB-Decane/data/SALB/'

Lipid = 'eggPC'
Run = 1

FixScale = False

Symmetric =True
FixTailExtension = True
OneBilayerRoughness = True        # all bilayer components have the same roughness (and its bigger than Si)

AsymRoughness = False

SameRoughness = False             # everything matches Si roughness
FixTails = True                   # use Lit values for tail volume
FixHeads = True                   # use Lit values for head volume

SysErr = 0.00

# assymetric fits have an issue with outer leaflet roughness smearing into the inner leaflet
# for small oil layers, need to fix all the internal roughnesses and only allow the outer tail/HG/water interfaces to change

## =========================================================================



if Run == 0:
    #eggPC 
    probe_d2o = load4(path+'213975_rescale.txt', back_reflectivity=True, name='D2O')
    probe_h2o = load4(path+'213978_rescale.txt', back_reflectivity=True, name='H2O')
    probe_smw = load4(path+'213981_rescale.txt', back_reflectivity=True, name='SMW')
    NoOil = True

if Run == 1:
    #eggPC dDecane 5:1
    probe_d2o = load4(path+'215672_rescale.txt', back_reflectivity=True, name='D2O')
    probe_h2o = load4(path+'215675_rescale.txt', back_reflectivity=True, name='H2O')
    probe_smw = load4(path+'215678_rescale.txt', back_reflectivity=True, name='SMW')
    NoOil = False
    dOil = True

if Run == 2:
    #eggPC dDecane 1:1
    probe_d2o = load4(path+'215083_rescale.txt', back_reflectivity=True, name='D2O')
    probe_h2o = load4(path+'215086_rescale.txt', back_reflectivity=True, name='H2O')
    probe_smw = load4(path+'215089_rescale.txt', back_reflectivity=True, name='SMW')
    NoOil = False
    dOil = True
    
if Run == 3:
    # eggPC dDecane 1:5. 
    probe_d2o = load4(path+'215095_rescale.txt', back_reflectivity=True, name='D2O')
    probe_h2o = load4(path+'215107_rescale.txt', back_reflectivity=True, name='H2O')
    probe_smw = load4(path+'215110_rescale.txt', back_reflectivity=True, name='SMW')
    NoOil = False
    dOil = True

if Run == 4:
    # eggPC hDecane 1:1. 
    probe_d2o = load4(path+'214039_rescale.txt', back_reflectivity=True, name='D2O')
    probe_h2o = load4(path+'214045_rescale.txt', back_reflectivity=True, name='H2O')
    probe_smw = load4(path+'214082_rescale.txt', back_reflectivity=True, name='SMW')
    NoOil = False
    dOil = False



# =============================================================================
# LOAD DATA
# =============================================================================
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
    # Intensity 
    probe_d2o.intensity.range(0.9,1.1)
    probe_h2o.intensity.range(0.8,1.2)
    probe_smw.intensity.range(0.8,1.2)
        
# =============================================================================



### === Component volumes ====================================================

V_CH = Parameter(name='CH volume', value=25.0, fixed=True)
V_CH2 = Parameter(name='CH2 volume', value=27.4, fixed=True)
V_CH3 = Parameter(name='CH3 volume', value=54.6, fixed=True)



if Lipid == 'DPhPC':
    N_CH, N_CH2, N_CH3 = 8, 20, 10
    #V_Tails = N_CH * V_CH + N_CH2 * V_CH2 + N_CH3 * V_CH3
    #from Tristam-Nagle 2010 CPL: V_tail = 1095, DHH = 36.4, A=80.5, VH=331
    V_Tails = Parameter(name='Tail volume ', value=1095.0, fixed=True)
    Max_Tail_Extension = 20.0      
    
if Lipid == 'eggPC':
    # eggPC, approximate as 16:0 and 18:1 chains
    N_CH, N_CH2, N_CH3 = 2, 28, 2
    #V_Tails = N_CH * V_CH + N_CH2 * V_CH2 + N_CH3 * V_CH3
    #from Tristam-Nagle 2010 CPL: V_tail = 1095, DHH = 36.4, A=80.5, VH=331
    V_Tails = Parameter(name='Tail volume ', value=1095.0, fixed=True)
    Max_Tail_Extension = 21.0      
    
if Lipid == 'DOPC':
    # DOPC, x2 and 18:1 chains
    N_CH, N_CH2, N_CH3 = 4, 28, 2
    #V_Tails = N_CH * V_CH + N_CH2 * V_CH2 + N_CH3 * V_CH3
    V_Tails = Parameter(name='Tail volume ', value=926.0, fixed=True)
    Max_Tail_Extension = 21.5

if FixHeads: V_PC = Parameter(name='PC volume', value=330.0, fixed=True)
if not FixHeads: V_PC = Parameter(name='PC volume', value=330.0).range(300,360)

if FixTails: V_Tails_adjust = Parameter(name='Tail volume adjustment', value=1.0)     
if not FixTails: V_Tails_adjust = Parameter(name='Tail volume adjustment', value=1.0).range(0.9,1.1)

# scattering lengths of relevant atoms
bH = -3.739; bD = 6.671; bC = 6.646; bO = 5.803; bN = 9.36; bP = 5.13

#HG is everything up to and including the carbonyls: C10H18NO8P
bHG = bC*10 + bH*18 +bN + bO*8 + bP
bTails = (N_CH + N_CH2 + N_CH3)*bC + (N_CH * bH + N_CH2 * bH * 2 + N_CH3 * bH * 3)



# =============================================================================
# MATERIALS
# =============================================================================
D2O_SLD = Parameter(name='D2O SLD', value=6.3).range(5.3,6.36)
D2O = SLD(name='D2O', rho=D2O_SLD, irho=0.0)
H2O = SLD(name='H2O', rho=-0.56, irho=0.0)
SMW = SLD(name='SMW', rho=(0.38*D2O_SLD) + (0.62*-0.56), irho=0.0)
                 
SiOx = SLD(name='siox', rho=3.5, irho=0.0000)
Si = SLD(name='silicon', rho=2.07, irho=0.0000)  
   
SiOx_SLD = Parameter(name='SiOx SLD', value=3.5)#.range(3,4)
dOil_SLD = Parameter(name='dOil SLD', value=6.56)#.range(3,6.56)
hOil_SLD = Parameter(name='hOil SLD', value=-0.48)#.range(0,6.5)

HG_SLD = Parameter(name='HG', value=bHG/V_PC*10, fixed=True)
Tail_SLD = Parameter(name='Tails', value=bTails/(V_Tails*V_Tails_adjust)*10, fixed=True)



## === Substrate parameters ==============================================
Si_rough = Parameter(name ='Si roughness', value=3).range(2, 4)
SiOx_VF = Parameter(name='SiOx VF', value=1, fixed = True)
SiOx_thickness = Parameter(name='SiOx thickness', value=15).range(10, 30)
SUB = Substrate(rho=Si.rho, sigma=Si_rough, name='substrate')


## === Bilayer parameters ==============================================
if Lipid == 'DPhPC':
    inAPM = Parameter(name='APMin', value=80.0).range(70,90)
    outAPM = Parameter(name='APMout', value=80.0).range(70,90)
    
if Lipid == 'DOPC':
    inAPM = Parameter(name='APMin', value=68.0).range(60,75)
    outAPM = Parameter(name='APMout', value=68.0).range(60,75)
    
if Lipid == 'eggPC':
    inAPM = Parameter(name='APMin', value=68.0).range(60,75)
    outAPM = Parameter(name='APMout', value=68.0).range(60,75)    
    
    
VF_Bilayer = Parameter(name='VF bilayer', value=1).range(0.8,1)
VF_inHG = Parameter(name='VF inner HG', value=0.5).range(0.4,0.6)
VF_outHG = Parameter(name='VF outer HG', value=0.5).range(0.4,0.6)

# these parameters allow the tail thickness to go from min to max extension
inTail_frac_extension = Parameter(name="inTail frac extension", value=0.2).range(0,1)
outTail_frac_extension = Parameter(name="outTail frac extension", value=0.2).range(0,1)

if FixTailExtension:
    inTail_frac_extension = Parameter(name="inTail frac extension", value=0, fixed=True)
    outTail_frac_extension = Parameter(name="outTail frac extension", value=0, fixed=True)



## === Roughness parameters ==============================================   
if not SameRoughness:
    inBilayer_rough_scale = Parameter(name='inBilayer roughness scale', value=1).range(1,2.5)
    outBilayer_rough_scale = Parameter(name='outBilayer roughness scale', value=1).range(1,2.5)
    inOil_rough_scale = Parameter(name='inOil roughness scale', value=1).range(1,2.5)
    outOil_rough_scale = Parameter(name='outOil roughness scale', value=1).range(1,2.5)

if OneBilayerRoughness: 
    outBilayer_rough_scale = Parameter(name='outBilayer roughness scale', value=1, fixed=True)
    inOil_rough_scale = Parameter(name='inOil roughness scale', value=1, fixed=True)
    outOil_rough_scale = Parameter(name='outOil roughness scale', value=1, fixed=True)

if SameRoughness:
    inBilayer_rough_scale = Parameter(name='inBilayer roughness scale', value=1, fixed=True)
    outBilayer_rough_scale = Parameter(name='outBilayer roughness scale', value=1, fixed=True)
    inOil_rough_scale = Parameter(name='inOil roughness scale', value=1, fixed=True)
    outOil_rough_scale = Parameter(name='outOil roughness scale', value=1, fixed=True)
    
    
    
    
if Symmetric:
    outAPM = inAPM
    VF_outHG = VF_inHG
    outTail_frac_extension = inTail_frac_extension

    
    if not AsymRoughness:
        outBilayer_rough_scale = inBilayer_rough_scale
        outOil_rough_scale = inOil_rough_scale
        inBilayer_rough = Parameter(name='inBilayer roughness', value=inBilayer_rough_scale*Si_rough, fixed=True)
        inOil_rough = Parameter(name='inOil roughness', value=inOil_rough_scale*inBilayer_rough, fixed=True)
        outOil_rough = Parameter(name='outOil roughness', value=inOil_rough, fixed=True)
        outBilayer_rough = Parameter(name='inBilayer roughness', value=inBilayer_rough, fixed=True)
    if AsymRoughness:
        inBilayer_rough = Parameter(name='inBilayer roughness', value=inBilayer_rough_scale*Si_rough, fixed=True)
        inOil_rough = Parameter(name='inOil roughness', value=inOil_rough_scale*inBilayer_rough, fixed=True)
        outOil_rough = Parameter(name='outOil roughness', value=outOil_rough_scale*inOil_rough, fixed=True)
        outBilayer_rough = Parameter(name='inBilayer roughness', value=outBilayer_rough_scale*outOil_rough, fixed=True)   
        

if not Symmetric:
    inBilayer_rough = Parameter(name='inBilayer roughness', value=inBilayer_rough_scale*Si_rough, fixed=True)
    inOil_rough = Parameter(name='inOil roughness', value=inOil_rough_scale*inBilayer_rough, fixed=True)
    outOil_rough = Parameter(name='outOil roughness', value=outOil_rough_scale*inOil_rough, fixed=True)
    outBilayer_rough = Parameter(name='inBilayer roughness', value=outBilayer_rough_scale*outOil_rough, fixed=True)

## === Thickness parameters ==============================================

Min_Bilayer_Separation = 2.0 * np.sqrt(Si_rough**2 + inBilayer_rough**2)   # non conformal version
Min_Bilayer_Separation = 2.0 * inBilayer_rough-Si_rough                    # fully conformal version

Offset = Parameter(name='Bilayer separation',value=0).range(0,10)
Water_thickness = Parameter(name='Water thickness', value=Min_Bilayer_Separation+Offset, fixed =True)

inHG_thickness = Parameter(name="inHG Thickness", value=V_PC / (VF_inHG * inAPM), fixed=True)
outHG_thickness = Parameter(name="outHG Thickness", value=V_PC / (VF_outHG * outAPM), fixed=True)

Oil_thickness = Parameter(name="Oil Thickness", value=1).range(0,20)

if Run == 0: Oil_thickness = Parameter(name="Oil Thickness", value=0, fixed=True)

## === Mixed Oil/HC Tail region parameters ==============================================
# this is the minimum tail extension, if there is no oil in them
Min_inTail_Extension = (V_Tails * V_Tails_adjust) / inAPM
Min_outTail_Extension = (V_Tails * V_Tails_adjust) / outAPM

inZ = Min_inTail_Extension + inTail_frac_extension * (Max_Tail_Extension - Min_inTail_Extension)
outZ = Min_outTail_Extension + outTail_frac_extension * (Max_Tail_Extension - Min_outTail_Extension)

inTail_thickness = Parameter(name="inTail Thickness", value=inZ, fixed=True)
outTail_thickness = Parameter(name="outTail Thickness", value=outZ, fixed=True)

inTail_OilVF = Parameter(name="inTail OilVF", value=(inZ - Min_inTail_Extension)/inZ , fixed=True)
outTail_OilVF = Parameter(name="outTail OilVF", value=(outZ - Min_outTail_Extension)/outZ, fixed=True)

inTail_SLD =  Parameter(name="inTail SLD", value=(inTail_OilVF * dOil_SLD) + (1-inTail_OilVF) * Tail_SLD , fixed=True)
outTail_SLD =  Parameter(name="outTail SLD", value=(outTail_OilVF * dOil_SLD) + (1-outTail_OilVF) * Tail_SLD, fixed=True)




## === Define the boxes  ==============================================
StartZ=SUB.overlap.value
oxide = VolumeFractionBox(name='oxide',       
                     z=StartZ+SiOx_thickness/2,
                     rhoH=SiOx_SLD,
                     rhoD=SiOx_SLD,
                     volume_fraction=1,
                     length=SiOx_thickness,
                     sigma_bottom=Si_rough,
                     sigma_top=Si_rough)

inHG = VolumeFractionBox(name='inHG',       
                     z=StartZ+SiOx_thickness+Water_thickness+inHG_thickness/2,
                     rhoH=HG_SLD,
                     rhoD=HG_SLD,
                     volume_fraction=VF_inHG*VF_Bilayer,
                     length=inHG_thickness,
                     sigma_bottom=inBilayer_rough,
                     sigma_top=inBilayer_rough)

inTail = VolumeFractionBox(name='inTail',       
                     z=StartZ+SiOx_thickness+Water_thickness+inHG_thickness+inTail_thickness/2,
                     rhoH=inTail_SLD,
                     rhoD=inTail_SLD,
                     volume_fraction=VF_Bilayer,
                     length=inTail_thickness,
                     sigma_bottom=inBilayer_rough,
                     sigma_top=inOil_rough)

if not NoOil:
    if dOil:
        Oil = VolumeFractionBox(name='dOil',       
                             z=StartZ+SiOx_thickness+Water_thickness+inHG_thickness+inTail_thickness+Oil_thickness/2,
                             rhoH=dOil_SLD,
                             rhoD=dOil_SLD,
                             volume_fraction=VF_Bilayer,
                             length=Oil_thickness,
                             sigma_bottom=inOil_rough,
                             sigma_top=outOil_rough)
    else:
        Oil = VolumeFractionBox(name='hOil',       
                             z=StartZ+SiOx_thickness+Water_thickness+inHG_thickness+inTail_thickness+Oil_thickness/2,
                             rhoH=hOil_SLD,
                             rhoD=hOil_SLD,
                             volume_fraction=VF_Bilayer,
                             length=Oil_thickness,
                             sigma_bottom=inOil_rough,
                             sigma_top=outOil_rough)

outTail = VolumeFractionBox(name='outTail',       
                     z=StartZ+SiOx_thickness+Water_thickness+inHG_thickness+inTail_thickness+Oil_thickness+outTail_thickness/2,
                     rhoH=outTail_SLD,
                     rhoD=outTail_SLD,
                     volume_fraction=VF_Bilayer,
                     length=outTail_thickness,
                     sigma_bottom=outOil_rough,
                     sigma_top=outBilayer_rough)

outHG = VolumeFractionBox(name='outHG',       
                     z=StartZ+SiOx_thickness+Water_thickness+inHG_thickness+inTail_thickness+Oil_thickness+outTail_thickness+outHG_thickness/2,
                     rhoH=HG_SLD,
                     rhoD=HG_SLD,
                     volume_fraction=VF_outHG*VF_Bilayer,
                     length=outHG_thickness,
                     sigma_bottom=outBilayer_rough,
                     sigma_top=outBilayer_rough)








## == Sample layer stack ==
layer_Si = Slab(material=Si, thickness=0.0000, interface=Si_rough)
layer_SiOx = Slab(material=SiOx, thickness=SiOx_thickness, interface=Si_rough)


# Use the bilayer definition function to generate the bilayer SLD profile, passing in the relevant parameters.
if not NoOil:
    TotalThickness=StartZ+SiOx_thickness+Water_thickness+inHG_thickness+inTail_thickness+Oil_thickness+outTail_thickness+outHG_thickness+outBilayer_rough*3
    sample_d2o, sample_h2o, sample_smw = make_samples(layer_template=MolgroupsLayer(base_group=SUB,
                                                                        overlay_groups=[oxide],
                                                                        add_groups=[inHG,inTail,Oil,outTail,outHG],   
                                                                        thickness=TotalThickness,
                                                                        name='Bilayer'),
                                          contrasts=[D2O, H2O, SMW], 
                                          substrate = layer_Si )

if NoOil:
    TotalThickness=StartZ+SiOx_thickness+Water_thickness+inHG_thickness+inTail_thickness+outTail_thickness+outHG_thickness+outBilayer_rough*3
    sample_d2o, sample_h2o, sample_smw = make_samples(layer_template=MolgroupsLayer(base_group=SUB,
                                                                        overlay_groups=[oxide],
                                                                        add_groups=[inHG,inTail,outTail,outHG],   
                                                                        thickness=TotalThickness,
                                                                        name='Bilayer'),
                                          contrasts=[D2O, H2O, SMW], 
                                          substrate = layer_Si )

## === Critical edge sampling ===
probe_d2o.critical_edge(substrate=Si, surface=D2O)


step = False
STEPSIZE=0.25

model_d2o = MolgroupsExperiment(sample=sample_d2o, probe=probe_d2o, dz=STEPSIZE, step_interfaces = step)
model_h2o = MolgroupsExperiment(sample=sample_h2o, probe=probe_h2o, dz=STEPSIZE, step_interfaces = step)
model_smw = MolgroupsExperiment(sample=sample_smw, probe=probe_smw, dz=STEPSIZE, step_interfaces = step)


problem = FitProblem([model_d2o, model_h2o, model_smw])

