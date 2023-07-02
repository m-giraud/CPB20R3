

import sys; 
CPBdir = "../../../CPlantBox"
sys.path.append(CPBdir+"/src");
sys.path.append(CPBdir);
sys.path.append("../../..");sys.path.append(".."); 
sys.path.append(CPBdir+"/src/python_modules");
sys.path.append("../build-cmake/cpp/python_binding/") # dumux python binding
sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../modules/") # python wrappers 

#import matplotlib
#matplotlib.use('AGG') #otherwise "import matplotlib.pyplot" hangs

from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
from functional.phloem_flux import PhloemFluxPython  # Python hybrid solver
#from Leuning_speedup import Leuning #about 0.7 for both
#from photosynthesis_cython import Leuning
import plantbox as pb
#from plantbox import Photosynthesis as ph
#import vtk_plot as vp
import math
import os
import numpy as np
import visualisation.vtk_plot as vp
#import matplotlib.pyplot as plt
from datetime import datetime, timedelta



def theta2H(vg,theta):#(-) to cm
    thetar =vg[0]# 0.059
    thetas = vg[1]#0.445
    alpha = vg[2]#0.00644
    n = vg[3]#1.503
    nrev = 1/(1-1/n)
    H =-(((( (thetas - thetar)/(theta - thetar))**nrev) - 1)**(1/n))/alpha
    return(H)#cm

def sinusoidal(t):
    return (np.sin(np.pi*t*2)+1)/2 

def qair2rh(qair, press,TK):
    T0 = 273.16
    RH = 26.3 * press * qair  /(exp((17.67*(TK - T0))/(TK- 29.65)))
    return RH



def div0(a, b, c):        
    return np.divide(a, b, out=np.full(len(a), c), where=b!=0)
    
def div0f(a, b, c):    
    if b != 0:
        return a/b
    else:
        return a/c
        



def setKrKx_xylem(TairC, RH,r): #inC
    #mg/cm3
    hPa2cm = 1.0197
    dEauPure = (999.83952 + TairC * (16.952577 + TairC * 
        (- 0.0079905127 + TairC * (- 0.000046241757 + TairC * 
        (0.00000010584601 + TairC * (- 0.00000000028103006)))))) /  (1 + 0.016887236 * TairC)
    siPhi = (30 - TairC) / (91 + TairC)
    siEnne=0
    mu =  pow(10, (- 0.114 + (siPhi * (1.1 + 43.1 * pow(siEnne, 1.25) )))) 
    mu = mu /(24*60*60)/100/1000; #//mPa s to hPa d, 1.11837e-10 hPa d for pure water at 293.15K
    mu = mu * hPa2cm #hPa d to cmh2o d 

    #number of vascular bundles
    VascBundle_leaf = 32
    VascBundle_stem = 52
    VascBundle_root = 1 #valid for all root type
            
    #radius of xylem type^4 * number per bundle
    rad_x_l_1   = (0.0015 **4) * 2; rad_x_l_2   = (0.0005 **4) * 2   
    rad_x_s_1   = (0.0017 **4) * 3; rad_x_s_2   = (0.0008 **4) * 1     
    rad_x_r0_1  = (0.0015 **4) * 4    
    rad_x_r12_1 = (0.00041**4) * 4; rad_x_r12_2 = (0.00087**4) * 1
    rad_x_r3_1  = (0.00068**4) * 1      

    # axial conductivity [cm^3/day]  
    betaXylX = 1#0.1      
    kz_l  = VascBundle_leaf *(rad_x_l_1 + rad_x_l_2)    *np.pi /(mu * 8)  * betaXylX
    kz_s  = VascBundle_stem *(rad_x_s_1 + rad_x_s_2)    *np.pi /(mu * 8) * betaXylX 
    kz_r0 = VascBundle_root * rad_x_r0_1                *np.pi /(mu * 8) * betaXylX  
    kz_r1 = VascBundle_root *(rad_x_r12_1 + rad_x_r12_2)*np.pi /(mu * 8)  * betaXylX
    kz_r2 = VascBundle_root *(rad_x_r12_1 + rad_x_r12_2)*np.pi /(mu * 8)  * betaXylX 
    kz_r3 = VascBundle_root * rad_x_r3_1                *np.pi /(mu * 8) * betaXylX

    #radial conductivity [1/day],0.00014 #
    betaXyl = 1#0.1#0.1
    kr_l  = 3.83e-5 * hPa2cm * betaXyl# init: 3.83e-4 cm/d/hPa
    kr_s  = 0.#1.e-20  * hPa2cm # set to almost 0
    kr_r0 =6.37e-5 * hPa2cm * betaXyl
    kr_r1 =7.9e-5  * hPa2cm * betaXyl
    kr_r2 =7.9e-5  * hPa2cm * betaXyl
    kr_r3 =6.8e-5  * hPa2cm * betaXyl
    l_kr = 0.8 #cm
    r.setKr([[kr_r0,kr_r1,kr_r2,kr_r0],[kr_s,kr_s ],[kr_l]], kr_length_=l_kr) 
    r.setKx([[kz_r0,kz_r1,kz_r2,kz_r0],[kz_s,kz_s ],[kz_l]])
    
    
    Rgaz=8.314 #J K-1 mol-1 = cm^3*MPa/K/mol
    rho_h2o = dEauPure/1000#g/cm3
    Mh2o = 18.05 #g/mol
    MPa2hPa = 10000
    hPa2cm = 1/0.9806806
    #log(-) * (cm^3*MPa/K/mol) * (K) *(g/cm3)/ (g/mol) * (hPa/MPa) * (cm/hPa) =  cm                      
    #p_a = np.log(RH) * Rgaz * rho_h2o * (TairC + 273.15)/Mh2o * MPa2hPa * hPa2cm
    #done withint solve photosynthesis
    #r.psi_air = p_a #*MPa2hPa #used only with xylem
    return r

    
def setKrKx_phloem(r): #inC

    #number of vascular bundles
    VascBundle_leaf = 32
    VascBundle_stem = 52
    VascBundle_root = 1 #valid for all root type
            
    #numPerBundle
    numL = 18
    numS = 21
    numr0 = 33
    numr1 = 25
    numr2 = 25
    numr3 = 1
        
    #radius of phloem type^4 * number per bundle
    rad_s_l   = numL* (0.00025 **4)# * 2; rad_x_l_2   = (0.0005 **4) * 2   
    rad_s_s   = numS *(0.00019 **4) #* 3; rad_x_s_2   = (0.0008 **4) * 1     
    rad_s_r0  = numr0 *(0.00039 **4) #* 4    
    rad_s_r12 = numr1*(0.00035**4) #* 4; rad_x_r12_2 = (0.00087**4) * 1
    rad_s_r3  = numr3 *(0.00068**4) #* 1      

    # axial conductivity [cm^3/day] , mu is added later as it evolves with CST  
    beta = 0.9 #Thompson 2003a
    kz_l   = VascBundle_leaf * rad_s_l   * np.pi /8 * beta  
    kz_s   = VascBundle_stem * rad_s_s   * np.pi /8 * beta
    kz_r0  = VascBundle_root * rad_s_r0  * np.pi /8 * beta
    kz_r12 = VascBundle_root * rad_s_r12 * np.pi /8 * beta
    kz_r3  = VascBundle_root * rad_s_r3  * np.pi /8 * beta
    
    #print([[kz_r0,kz_r12,kz_r12,kz_r3],[kz_s,kz_s ],[kz_l]])
    #raise Exception
    #radial conductivity [1/day],
    kr_l  = 0.#3.83e-4 * hPa2cm# init: 3.83e-4 cm/d/hPa
    kr_s  = 0.#1.e-20  * hPa2cm # set to almost 0
    # kr_r0 = 1e-1
    # kr_r1 = 1e-1
    # kr_r2 = 1e-1
    # kr_r3 = 1e-1
    kr_r0 = 5e-2
    kr_r1 = 5e-2
    kr_r2 = 5e-2
    kr_r3 = 5e-2
    l_kr = 0.8 #0.8 #cm
    
    r.setKr_st([[kr_r0,kr_r1 ,kr_r2 ,kr_r0],[kr_s,kr_s ],[kr_l]] , kr_length_= l_kr)
    r.setKx_st([[kz_r0,kz_r12,kz_r12,kz_r0],[kz_s,kz_s ],[kz_l]])
    
    a_ST = [[0.00039,0.00035,0.00035,0.00039 ],[0.00019,0.00019],[0.00025]]
    Across_s_l   = numL*VascBundle_leaf *(a_ST[2][0]**2)*np.pi# (0.00025 **2)# * 2; rad_x_l_2   = (0.0005 **4) * 2   
    Across_s_s   = numS *VascBundle_stem * (a_ST[1][0]**2)*np.pi#(0.00019 **2) #* 3; rad_x_s_2   = (0.0008 **4) * 1     
    Across_s_r0  = numr0 *VascBundle_root * (a_ST[0][0]**2)*np.pi#(0.00039 **2) #* 4    
    Across_s_r12 = numr1*VascBundle_root * (a_ST[0][1]**2)*np.pi#(0.00035**2) #* 4; rad_x_r12_2 = (0.00087**4) * 1
    Across_s_r3  =  numr3 *VascBundle_root *(a_ST[0][2]**2)*np.pi# (0.00068**2) #* 1    
    
    Perimeter_s_l   = numL*VascBundle_leaf *(a_ST[2][0])* 2 * np.pi# (0.00025 **2)# * 2; rad_x_l_2   = (0.0005 **4) * 2   
    Perimeter_s_s   = numS *VascBundle_stem * (a_ST[1][0])* 2 * np.pi#(0.00019 **2) #* 3; rad_x_s_2   = (0.0008 **4) * 1     
    Perimeter_s_r0  = numr0 *VascBundle_root * (a_ST[0][0])* 2 * np.pi#(0.00039 **2) #* 4    
    Perimeter_s_r12 = numr1*VascBundle_root * (a_ST[0][1])* 2 * np.pi#(0.00035**2) #* 4; rad_x_r12_2 = (0.00087**4) * 1
    Perimeter_s_r3  =  numr3 *VascBundle_root *(a_ST[0][2])* 2 * np.pi# (0.00068**2) #* 1  
    r.setAcross_st([[Across_s_r0,Across_s_r12,Across_s_r12,Across_s_r0],[Across_s_s,Across_s_s],[Across_s_l]])
    return r

def resistance2conductance(resistance,r, weatherX):
    resistance = resistance* (1/100) #[s/m] * [m/cm] = [s/cm]
    resistance = resistance * r.R_ph * weatherX["TairK"] / r.Patm # [s/cm] * [K] * [hPa cm3 K−1 mmol−1] * [hPa] = [s] * [cm2 mmol−1]
    resistance = resistance * (1000) * (1/10000)# [s cm2 mmol−1] * [mmol/mol] * [m2/cm2] = [s m2 mol−1]
    return 1/resistance

def write_file_array(name, data, directoryN):
    name2 = 'results'+ directoryN+ name+ '.txt'
    with open(name2, 'a') as log:
        log.write(','.join([num for num in map(str, data)])  +'\n')

def write_file_float(name, data, directoryN):
    name2 = 'results' + directoryN+  name+ '.txt'
    with open(name2, 'a') as log:
        log.write(repr( data)  +'\n')

def weather(simDuration, hp, condition, simStartSim):
    vgSoil = [0.059, 0.45, 0.00644, 1.503, 1]
    loam = [0.08, 0.43, 0.04, 1.6, 50]
    Qnigh = 0; Qday = 960e-6 #458*2.1
    if ((condition == "wet") or (simDuration <= simStartSim) or (simDuration > simStartSim +7)):
        Tnigh = 15.8; Tday = 22
        RHday = 0.6; RHnigh = 0.88
        Pair = 1010.00 #hPa
        thetaInit = 0.4#
        cs = 350e-6
    elif condition == "dry":
        Tnigh = 20.7; Tday = 30.27
        RHday = 0.44; RHnigh = 0.78
        Pair = 1070.00 #hPa
        thetaInit = 28/100   
        cs = 350e-6
    else:
        print("condition",condition)
        raise Exception("condition not recognised")

    coefhours = sinusoidal(simDuration)
    RH_ = RHnigh + (RHday - RHnigh) * coefhours
    TairC_ = Tnigh + (Tday - Tnigh) * coefhours
    Q_ = Qnigh + (Qday - Qnigh) * coefhours
    es =  6.112 * np.exp((17.67 * TairC_)/(TairC_ + 243.5))
    ea = es*RH_
    assert ea < es
    
    assert ((RH_ > 0) and(RH_ < 1))
    bl_thickness = 1/1000 #1mm * m_per_mm
    diffusivity= 2.5e-5#m2/sfor 25*C
    rbl =bl_thickness/diffusivity #s/m 13
    
    Kcanopymean = 1e-1 # m2/s
    meanCanopyL = (2/3) * hp /2
    rcanopy = meanCanopyL/Kcanopymean
    windSpeed = 2 #m/s
    zmzh = 2 #m
    karman = 0.41 #[-]

    rair = 1
    if hp > 0:
        rair = np.log((zmzh - (2/3)*hp)/(0.123*hp)) * np.log((zmzh - (2/3)*hp)/(0.1*hp)) / (karman*karman*windSpeed)

    pmean = theta2H(vgSoil, thetaInit)

    weatherVar = {'TairC' : TairC_,'TairK' : TairC_ + 273.15,'Pair':Pair,"es":es,
                    'Qlight': Q_,'rbl':rbl,'rcanopy':rcanopy,'rair':rair,"ea":ea,
                    'cs':cs, 'RH':RH_, 'p_mean':pmean, 'vg':loam}
    print("Env variables at", round(simDuration//1),"d",round((simDuration%1)*24),"hrs :\n", weatherVar)
    return weatherVar

def setFunctionalParameters(r, weatherInit):
    r.g0 = 8e-6
    r.VcmaxrefChl1 =1.28
    r.VcmaxrefChl2 = 8.33
    r.a1 = 0.6/0.4
    r.a3 = 1.5
    r.alpha = 0.4
    r.theta = 0.6
    r.k_meso = 1e-3
    r.setKrm2([[2e-5]])
    r.setKrm1([[10e-2]])
    r.setRhoSucrose([[0.51],[0.65],[0.56]])
    rootFact = 2
    r.setRmax_st([[2.4*rootFact,1.5*rootFact,0.6*rootFact,2.4*rootFact],[2.,2.],[8.]])
    
    r.KMrm = 0.1
    r.sameVolume_meso_st = False
    r.sameVolume_meso_seg = True
    r.withInitVal =True
    r.initValST = 0.
    r.initValMeso = 0.
    r.beta_loading = 0.6
    r.Vmaxloading = 0.05 #mmol/d
    r.Mloading = 0.2
    r.Gr_Y = 0.8
    r.CSTimin = 0.4
    r.surfMeso=0.0025
    r.leafGrowthZone = 2 # cm
    r.StemGrowthPerPhytomer = True # 
    r.psi_osmo_proto = -10000*1.0197 #schopfer2006
    r.fwr = 0#0.1
    r.sh = 4e-4
    r.gm=0.01
    
    r.limMaxErr = 1/100
    r.maxLoop = 10000
    r.minLoop=900
    
    r.C_targ = r.CSTimin
    r.C_targMesophyll = r.CSTimin
    r.k_S_ST = 5/25 #daudet2002
    r.k_S_Mesophyll = 5/25 #daudet2002
   


    r.cs = weatherInit["cs"]

    r.expression = 6
    r.update_viscosity = True
    r.solver = 1
    r.atol = 1e-10
    r.rtol = 1e-6
    
    SPAD= 41.0
    chl_ = (0.114 *(SPAD**2)+ 7.39 *SPAD+ 10.6)/10
    r.Chl = np.array( [chl_]) 
    r.Csoil = 1e-4
    return r

def getVTPOutput(r, C_ST, psiXyl, vtpindx, directoryN):
    ana = pb.SegmentAnalyser(r.plant.mappedSegments())
    cutoff = 1e-15 #is get value too small, makes paraview crash
    C_ST_p = C_ST
    C_ST_p[abs(C_ST_p) < cutoff] = 0
    psiXyl_p = psiXyl 
    psiXyl_p[abs(psiXyl_p) < cutoff] = 0


    ana.addData("CST", C_ST_p)
    
    ana.addData("psi_Xyl",psiXyl_p)
    ana.write("results"+directoryN+"plot_"+str(int(simStartSim))+str(condition)+"at"+ str(vtpindx) +".vtp", 
              ["organType", "subType",
               "CST", "psi_Xyl"]) 

    vp.plot_plant(r.plant,p_name = ["xylem pressure (cm)","sucrose concentration (mmol/cm3)"],
                        vals =[  psiXyl_p, C_ST_p], 
                        filename = "results"+ directoryN +"plotpsi_"+str(int(simStartSim))+str(condition)+"at"+ str(vtpindx), 
                  range_ = [0,5000])
