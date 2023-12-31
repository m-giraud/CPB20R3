""" water movement within the root (static soil) """
#directoryN = "/7to14dry/"
import sys; 
#CPBdir = "../../../CPlantBox"
CPBdir = "../../../../cpb3101/CPlantBox"
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
from phloem_flux import PhloemFluxPython  # Python hybrid solver
#from Leuning_speedup import Leuning #about 0.7 for both
#from photosynthesis_cython import Leuning
import plantbox as pb
#from plantbox import Photosynthesis as ph
#import vtk_plot as vp
import math
import os
import numpy as np
import vtk_plot as vp
#import matplotlib.pyplot as plt
from datetime import datetime, timedelta



isCluster = (os.environ['HOME'] == '/home/m.giraud')

#       qr, qs, alpha, n, ks (in [cm/d])
#vgSoil = [0.059, 0.45, 0.00644, 1.503, 1]
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

def qair2rh(qair, es_, press):
    e =qair * press / (0.378 * qair + 0.622)
    rh = e / es_
    rh=max(min(rh, 1.),0.)
    return rh



def div0(a, b, c):        
    return np.divide(a, b, out=np.full(len(a), c), where=b!=0)
    
def div0f(a, b, c):    
    if b != 0:
        return a/b
    else:
        return a/c
        



def setKrKx_xylem(r, TairC, RH): #inC
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
    kz_l  = VascBundle_leaf *(rad_x_l_1 + rad_x_l_2)    *np.pi /(mu * 8)  
    kz_s  = VascBundle_stem *(rad_x_s_1 + rad_x_s_2)    *np.pi /(mu * 8) 
    kz_r0 = VascBundle_root * rad_x_r0_1                *np.pi /(mu * 8)  
    kz_r1 = VascBundle_root *(rad_x_r12_1 + rad_x_r12_2)*np.pi /(mu * 8) 
    kz_r2 = VascBundle_root *(rad_x_r12_1 + rad_x_r12_2)*np.pi /(mu * 8)  
    kz_r3 = VascBundle_root * rad_x_r3_1                *np.pi /(mu * 8) # 4.32e-1

    #radial conductivity [1/day],
    kr_l  = 3.83e-4 * hPa2cm# init: 3.83e-4 cm/d/hPa
    kr_s  = 0.#1.e-20  * hPa2cm # set to almost 0
    kr_r0 = 6.37e-5 * hPa2cm 
    kr_r1 = 7.9e-5  * hPa2cm 
    kr_r2 = 7.9e-5  * hPa2cm  
    kr_r3 = 6.8e-5  * hPa2cm 
    l_kr = 0.8 #cm
    r.setKr([[kr_r0,kr_r1,kr_r2,kr_r0],[kr_s,kr_s ],[kr_l]], kr_length_=l_kr) 
    r.setKx([[kz_r0,kz_r1,kz_r2,kz_r0],[kz_s,kz_s ],[kz_l]])
    
    
    Rgaz=8.314 #J K-1 mol-1 = cm^3*MPa/K/mol
    rho_h2o = dEauPure/1000#g/cm3
    Mh2o = 18.05 #g/mol
    MPa2hPa = 10000
    hPa2cm = 1/0.9806806
    #log(-) * (cm^3*MPa/K/mol) * (K) *(g/cm3)/ (g/mol) * (hPa/MPa) * (cm/hPa) =  cm                      
    p_a = np.log(RH) * Rgaz * rho_h2o * (TairC + 273.15)/Mh2o * MPa2hPa * hPa2cm

    r.psi_air = p_a #*MPa2hPa #used only with xylem
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
    l_kr = 0.8 #cm
    
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
    #print(a_ST[2][0],a_ST[1][0],a_ST[0][0],a_ST[0][1],a_ST[0][2])
    #r.a_ST = a_ST #to check for water equilibrium assumption
    #tot surface/np.pi of sieve tube  (np.pi added after)
    #r.a_ST_eqs = [[rad_s_r0,rad_s_r12,rad_s_r12,rad_s_r0],[rad_s_s,rad_s_s],[rad_s_l]]
    r.setAcross_st([[Across_s_r0,Across_s_r12,Across_s_r12,Across_s_r0],[Across_s_s,Across_s_s],[Across_s_l]])
    return r
    
def addParams(r,weatherInit):
    
    r.g0 = 8e-3
    r.VcmaxrefChl1 =1.28#/2
    r.VcmaxrefChl2 = 8.33#/2
    r.a1 = 0.6/0.4#0.7/0.3#0.6/0.4 #ci/(cs - ci) for ci = 0.6*cs
    r.a3 = 1.5
    r.alpha = 0.4#0.2#/2
    r.theta = 0.6#0.9#/2
    r.k_meso = 1e-3#1e-4
    r.setKrm2([[2e-5]])
    r.setKrm1([[10e-2]])#([[2.5e-2]])
    r.setRhoSucrose([[0.51],[0.65],[0.56]])#0.51
    r.setRmax_st([[14.4,9.0,6.0,14.4],[5.,5.],[15.]])#*6 for roots, *1 for stem, *24/14*1.5 for leaves
    #r.setRmax_st([[12,9.0,6.0,12],[5.,5.],[15.]])
    r.KMrm = 0.1#VERY IMPORTANT TO KEEP IT HIGH
    #r.exud_k = np.array([2.4e-4])#*10#*(1e-1)
    #r.k_gr = 1#0
    r.sameVolume_meso_st = False
    r.sameVolume_meso_seg = True
    r.withInitVal =True
    r.initValST = 0.#0.6#0.0
    r.initValMeso = 0.#0.9#0.0
    r.beta_loading = 0.6
    r.Vmaxloading = 0.05 #mmol/d, needed mean loading rate:  0.3788921068507634
    r.Mloading = 0.2
    r.Gr_Y = 0.8
    r.CSTimin = 0.4
    r.surfMeso=0.0025

    r.cs = weatherInit["cs"]

    #r.r_forPhloem(24/14*1.5, 4)
    #r.r_forPhloem(24/14, 3)
    #r.r_forPhloem(6, 2) #because roots will have high C_ST less often
    r.expression = 6
    r.update_viscosity = True
    r.solver = 1
    r.atol = 1e-12
    r.rtol = 1e-8
    #r.doNewtonRaphson = False;r.doOldEq = False
    SPAD= 41.0
    chl_ = (0.114 *(SPAD**2)+ 7.39 *SPAD+ 10.6)/10
    r.Chl = np.array( [chl_]) 
    r.Csoil = 1e-4
    return r
    
""" Parameters """
def launchUQR(directoryN,simInit, condition,forPlants):
    def write_file_array(name, data):
        name2 = 'results'+ directoryN+ name+ '_15pm.txt'
        with open(name2, 'a') as log:
            log.write(','.join([num for num in map(str, data)])  +'\n')

    def write_file_float(name, data):
        name2 = 'results' + directoryN+  name+ '_15pm.txt'
        with open(name2, 'a') as log:
            log.write(repr( data)  +'\n')
            
    def weather(simDuration):
        vgSoil = [0.059, 0.45, 0.00644, 1.503, 1]
        Qmin = 0; Qmax = 960e-6 #458*2.1
        if condition == "dry":
            Tmin = 20.7; Tmax = 30.27
            specificHumidity = 0.0111
            Pair = 1070.00 #hPa
            thetaInit = 20/100
        elif condition == "wet":
            Tmin = 15.8; Tmax = 22
            specificHumidity = 0.0097
            Pair = 1010.00 #hPa
            thetaInit = 30/100
            
        else:
            print("condition",condition)
            raise Exception("condition not recognised")

        coefhours = sinusoidal(simDuration)
        TairC_ = Tmin + (Tmax - Tmin) * coefhours
        Q_ = Qmin + (Qmax - Qmin) * coefhours
        cs = 350e-6 #co2 paartial pressure at leaf surface (mol mol-1)
        #RH = 0.5 # relative humidity
        es =  6.112 * np.exp((17.67 * TairC_)/(TairC_ + 243.5))
        RH = qair2rh(specificHumidity, es, Pair)
        ea = es*RH

        pmean = theta2H(vgSoil, thetaInit)

        weatherVar = {'TairC' : TairC_,
                        'Qlight': Q_,"es":es,"ea":ea,
                        'cs':cs, 'RH':RH, 'p_mean':pmean, 'vg':vgSoil}
        print("Env variables at", round(simDuration//1),"d",round((simDuration%1)*24),"hrs :\n", weatherVar)
        return weatherVar

    weatherInit = weather(0)
    #simInit = 7
    simDuration = simInit # [day] init simtime
    simMax =simInit+20
    depth = 60
    dt = 1/24 #1h
    verbose = True

    # plant system 
    pl = pb.MappedPlant(seednum = 2) #pb.MappedRootSystem() #pb.MappedPlant()
    pl2 = pb.MappedPlant(seednum = 2) #pb.MappedRootSystem() #pb.MappedPlant()
    path = CPBdir+"/modelparameter/plant/"
    name = "Triticum_aestivum_adapted_2023"#"Triticum_aestivum_test_2021"#"Triticum_aestivum_adapted_2021"#"wheat_uqr15" #"manyleaves"## "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010

    pl.readParameters(path + name + ".xml")
    pl2.readParameters(path + name + ".xml")
    
    class My_Tropism(pb.Tropism):
        """ User tropism 2: depending on root age use plagio- or gravitropism """

        def __init__(self, rs, n, sigma, age):
            super(My_Tropism, self).__init__(rs)
            self.exo = pb.Exotropism(rs, 0., 0.)
            self.gravi = pb.Gravitropism(rs, 0., 0.)
            self.setTropismParameter(n, sigma)

        def tropismObjective(self, pos, old, a, b, dx, root):
            return 1 * self.exo.tropismObjective(pos, old, a, b, dx, root)+0.1 *self.gravi.tropismObjective(pos, old, a, b, dx, root)
        
    if forPlants == "mix":
        for p in pl.getOrganRandomParameter(pb.root):
            p.theta =  0#p.theta * 4
            p.tropismT = 1
            p.tropismS = p.tropismS*10
            p.tropismN = p.tropismN*10
        for p in pl2.getOrganRandomParameter(pb.root):
            if p.subType ==1:
                p.theta =  1.5#p.theta * 4
                p.tropismT = 2
                p.tropismS = p.tropismS*10
                p.tropismN = p.tropismN*10
            else:
                p.theta =  0
            
        
        #mytropism2 = My_Tropism(pl2, 1.5, 0.5, 5)  # after 5 days switch from plagio- to gravitropism
        #pl2.setTropism(mytropism2, 2, 1)  # 4 for base roots, -1 for all root types
        #pl2.setTropism(mytropism2, 2, 2)  # 4 for base roots, -1 for all root types
        #pl2.setTropism(mytropism2, 2, 3)  # 4 for base roots, -1 for all root types
        #pl2.setTropism(mytropism2, 2, 4)  # 4 for base roots, -1 for all root types
        #for p in pl2.getOrganRandomParameter(pb.root):
         #   p.tropismT = 5
          #  p.theta = 1.5708 
           # p.tropismS = p.tropismS/10
            #p.tropismN = p.tropismN*10
    if forPlants == "deep":
        pass
    if forPlants == "shallow":
        for p in pl.getOrganRandomParameter(pb.root):
            p.theta =  p.theta * 2
            p.tropismS = p.tropismS*100
            p.tropismN = p.tropismN*100
        for p in pl2.getOrganRandomParameter(pb.root):
            p.theta =  p.theta * 2
            p.tropismS = p.tropismS*100
            p.tropismN = p.tropismN*100


    #raise Exception
    sdf = pb.SDF_PlantBox(np.Inf, np.Inf, depth )
    layers = depth; soilvolume = (depth / layers) * 3 * 12

    pl.setGeometry(sdf) # creates soil space to stop roots from growing out of the soil
    pl2.setGeometry(sdf) # creates soil space to stop roots from growing out of the soil


    pl.getOrganRandomParameter(pb.seed)[0].seedPos = pb.Vector3d(1.5,0, -5) 
    pl2.getOrganRandomParameter(pb.seed)[0].seedPos = pb.Vector3d(-1.5,0, -5) 
    
    pl.initialize(verbose = True)#, stochastic = False)
    pl2.initialize(verbose = True)#, stochastic = False)
    
    if forPlants == "no":
        for p in pl.getOrganRandomParameter(pb.root):
            p.theta =  0#p.theta * 4
            p.tropismT = 1
            p.tropismS = p.tropismS*10
            p.tropismN = p.tropismN*10
        
        mytropism2 = My_Tropism(pl2, 1.5, 0.5, 5)  # after 5 days switch from plagio- to gravitropism
        pl2.setTropism(mytropism2, 2, 1)  # 4 for base roots, -1 for all root types
        pl2.setTropism(mytropism2, 2, 2)  # 4 for base roots, -1 for all root types
        pl2.setTropism(mytropism2, 2, 3)  # 4 for base roots, -1 for all root types
        pl2.setTropism(mytropism2, 2, 4)  # 4 for base roots, -1 for all root types
        
    print(pl.getOrganRandomParameter(pb.seed)[0].seedPos)
    print(pl2.getOrganRandomParameter(pb.seed)[0].seedPos)
    
    pl.simulate(simDuration, False)#, "outputpm15.txt")
    pl2.simulate(simDuration, False)#, "outputpm15.txt")
    
    ö=0

    beginning = datetime.now()





    while simDuration <=  100: 

        print('simDuration:',simDuration )

        ana = pb.SegmentAnalyser(pl.mappedSegments())
        ana2 = pb.SegmentAnalyser(pl2.mappedSegments())

        ana.write("results"+directoryN+"p1p_"+ str(ö) +".vtp", 
                  ["organType", "subType"]) 
        ana2.write("results"+directoryN+"p2p_"+ str(ö) +".vtp", 
                   ["organType", "subType"]) 

        rl_ = ana.distribution("length", 0., -depth, layers, True)                   
        rl_ = np.array(rl_)/ soilvolume  # convert to density
        rl_2 = ana2.distribution("length", 0., -depth, layers, True)                   
        rl_2 = np.array(rl_2)/ soilvolume  # convert to density 
        ö +=1

        write_file_array("rl_1", rl_)
        write_file_array("rl_2", rl_2)
        

        verbose_simulate = False
        pl.simulate(dt, verbose_simulate)
        pl2.simulate(dt,  verbose_simulate)
        simDuration += dt


    end = datetime.now()
    print(end - beginning)
    

if __name__ == '__main__':    
    directoryN = "/"+os.path.basename(__file__)[:-3]+"/"

    main_dir=os.environ['PWD']#dir of the file
    results_dir = main_dir +"/results"+directoryN
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    else:
        import shutil
        shutil.rmtree(results_dir)
        os.makedirs(results_dir)

    launchUQR(directoryN,7, "dry")