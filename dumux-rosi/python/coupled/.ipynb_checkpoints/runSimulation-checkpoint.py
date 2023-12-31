

import sys; 
CPBdir = "../../../CPlantBox"
sys.path.append(CPBdir+"/src");
sys.path.append(CPBdir);
sys.path.append("../../..");sys.path.append(".."); 
sys.path.append(CPBdir+"/src/python_modules");
sys.path.append("../build-cmake/cpp/python_binding/") # dumux python binding
sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../modules/") # python wrappers 

from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
from functional.phloem_flux import PhloemFluxPython  # Python hybrid solver
from runSimulation_modules import *

import plantbox as pb
import math
import os
import numpy as np
import visualisation.vtk_plot as vp

from datetime import datetime, timedelta



    
""" Parameters """
def launchUQR(directoryN_,simInit,simStartSim, condition,spellDuration, simEnd):
    global directoryN
    directoryN = directoryN_ 

    weatherInit = weather(0,0,condition, simStartSim)
    simDuration = simInit # [day] init simtime
    
    simMax = simEnd
    depth = 60
    dt = 1/24 
    verbose = True

    # plant system 
    pl = pb.MappedPlant(seednum = 2) 
    path = CPBdir+"/modelparameter/structural/plant/"
    name = "Triticum_aestivum_test_2021"

    pl.readParameters(path + name + ".xml")
    sdf = pb.SDF_PlantBox(np.Inf, np.Inf, depth )

    pl.setGeometry(sdf) # creates soil space to stop roots from growing out of the soil


    pl.initialize(verbose = True)
    pl.simulate(simDuration, False)
    
    
    """ Coupling to soil """
    min_b = [-3./2, -12./2, -41.]#distance between wheat plants
    max_b = [3./2, 12./2, 0.]
    rez = 1
    cell_number = [int(6*rez), int(24*rez), int(40*rez)]#1cm3? 
    layers = depth; soilvolume = (depth / layers) * 3 * 12
    k_soil = []
    initial = weatherInit["p_mean"]#mean matric potential [cm] pressure head
    
    p_mean = initial
    p_bot = p_mean + depth/2
    p_top = p_mean - depth/2
    sx = np.linspace(p_top, p_bot, depth)
    picker = lambda x,y,z : max(int(np.floor(-z)),-1) 
    sx_static_bu = sx    
    pl.setSoilGrid(picker)  # maps segment


    """ Parameters phloem and photosynthesis """
    r = PhloemFluxPython(pl,psiXylInit = min(sx),ciInit = weatherInit["cs"]*0.5) 

    r = setKrKx_phloem(r)
    r = setFunctionalParameters(r, weatherInit)

    """ for post processing """
    vtpIndx=0 
    beginning = datetime.now()
    
    
    orgs_all = r.plant.getOrgans(-1, True)
    volOrg = np.array([org.orgVolume(-1,False) for org in orgs_all]) 
    ot_orgs = np.array([org.organType() for org in orgs_all])
    st_orgs = np.array([org.getParameter("subType") for org in orgs_all])
    lenOrg = np.array([org.getLength(False) for org in orgs_all]) 
    
    Q_in = 0
    Q_ST_init = np.array([])
    Q_meso_init  = np.array([])
    Q_Rmbu      = np.array([0.])

    Q_Grbu      = np.array([0.])
    
    Q_Exudbu    = np.array([0.])
    Q_Rmmaxbu   = np.array([0.])
    Q_Grmaxbu   = np.array([0.])
    Q_STbu      = np.array([0.])
    Q_S_Mesophyllbu      = np.array([0.])
    Q_S_STbu      = np.array([0.])
    Q_mesobu    = np.array([0.])
    dynamic_soil_first = True
    dtVTP = 10
    dtPrint = 10
    
    while simDuration < simMax: 

        print('simDuration:',simDuration )

        
        hp = max([tempnode[2] for tempnode in r.get_nodes()]) /100 #canopy height
        

        weatherX = weather(simDuration, hp, condition, simStartSim)
        r.Patm = weatherX["Pair"]
        
        ##resistances
        r.g_bl = resistance2conductance(weatherX["rbl"],r, weatherX) / r.a2_bl
        r.g_canopy = resistance2conductance(weatherX["rcanopy"],r, weatherX) / r.a2_canopy
        r.g_air = resistance2conductance(weatherX["rair"],r, weatherX) / r.a2_air

        r.Qlight = weatherX["Qlight"]


        r = setKrKx_xylem(weatherX["TairC"], weatherX["RH"],r)
        
        dynamic_soil = ((simDuration > simStartSim) and (simDuration <= simStartSim +spellDuration) and (condition == "dry"))
        
        if dynamic_soil and dynamic_soil_first:
            dynamic_soil_first = False
            s = RichardsWrapper(RichardsSP())
            s.initialize()
            periodic = True
            s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
            s.setHomogeneousIC(weatherX["p_mean"], True)  # cm pressure head, equilibrium
            s.setTopBC("noFlux")
            s.setBotBC("fleeFlow")
            s.setVGParameters([weatherX['vg']])
            s.initializeProblem()
            s.setCriticalPressure(-15000)
            sx = s.getSolutionHead()  # inital condition, solverbase.py
            picker = lambda x, y, z: s.pick([x, y, z])    
            pl.setSoilGrid(picker)  # maps segment
            dtWater = 1/24/60
            
        if not dynamic_soil:
            sx = sx_static_bu
            picker = lambda x,y,z : max(int(np.floor(-z)),-1) 
            pl.setSoilGrid(picker)  # maps segment
            dtWater = dt
        
        
            
        dtWatertot = 0
        Evdt =0
        Andt = 0
        while dtWatertot < dt:
            dtWatertot += dtWater
            
            """ photosynthesis """
            r.solve_photosynthesis(sim_time_ = simDuration, sxx_=sx, cells_ = True,ea_ = weatherX["ea"],es_=weatherX["es"],
                verbose_ = False, doLog_ = False,TairC_= weatherX["TairC"],outputDir_= "./results"+directoryN)

            """ dumux """   
            fluxesSoil = r.soilFluxes(simDuration, r.psiXyl, sx, approx=False)
            
            if dynamic_soil:
                s.setSource(fluxesSoil.copy())  # richards.py 
                s.solve(dtWater)
                sx = s.getSolutionHead()  # richards.py  

                min_sx, min_rx, max_sx, max_rx = np.min(sx), np.min(r.psiXyl), np.max(sx), np.max(r.psiXyl)
                n = round((simDuration- simInit)/(simMax-simInit) * 100.)

                print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], [{:g}, {:g}] cm soil [{:g}, {:g}] cm root at {:g} days {:g}"
                        .format(min_sx, max_sx, min_rx, max_rx, s.simTime, r.psiXyl[0]))

            Andt += np.array(r.Ag4Phloem)*dtWater
            Evdt += np.sum(r.Ev)*dtWater
        
        
        """ sucrose balance """
        
        r.Ag4Phloem = Andt/dt
        
        verbose_phloem = True
        
        filename = "results"+ directoryN +"inPM.txt"  
        try:
            r.startPM(simDuration, simDuration + dt, 1, ( weatherX["TairC"]  +273.15) , verbose_phloem, filename)
        except:
            r.doTroubleshooting = True
            r.startPM(simDuration, simDuration + dt, 1, ( weatherX["TairC"]  +273.15) , verbose_phloem, filename)
            
        try:
            os.remove(filename)
        except OSError:
            pass    
        
        """ get results """
        
        errLeuning = sum(r.outputFlux)
        fluxes = np.array(r.outputFlux)
        
        Nt = len(r.plant.nodes) 
        if r.withInitVal and (len(Q_ST_init) ==0) :
            Q_ST_init = np.array(r.Q_init[0:Nt])
            Q_meso_init = np.array(r.Q_init[Nt:(Nt*2)])
            
        Q_ST    = np.array(r.Q_out[0:Nt])
        Q_meso  = np.array(r.Q_out[Nt:(Nt*2)])
        Q_Rm    = np.array(r.Q_out[(Nt*2):(Nt*3)])
        Q_Exud  = np.array(r.Q_out[(Nt*3):(Nt*4)])
        Q_Gr    = np.array(r.Q_out[(Nt*4):(Nt*5)])


        C_ST    = np.array(r.C_ST)
        Q_S_Mesophyll   = np.array(r.Q_out[(Nt*7):(Nt*8)])
        Q_S_ST   = np.array(r.Q_out[(Nt*8):(Nt*9)])
        Fl      = np.array(r.Fl)
        volST   = np.array(r.vol_ST)
        volMeso   = np.array(r.vol_Meso)
        C_S_Mesophyll   = Q_S_Mesophyll/volMeso
        C_S_ST   = Q_S_ST/volST
        C_meso  = Q_meso/volMeso
        Q_in   += sum(np.array(r.AgPhl)*dt)
        Q_out   = Q_Rm + Q_Exud + Q_Gr
        error   = sum(Q_ST +Q_S_ST+ Q_meso + Q_out + Q_S_Mesophyll)- Q_in - sum(Q_ST_init)  - sum(Q_meso_init)
        

        Q_Rmmax       = np.array(r.Q_out[(Nt*5):(Nt*6)])
        Q_Grmax       = np.array(r.Q_out[(Nt*6):(Nt*7)])

        Q_ST_i        = Q_ST      - Q_STbu
        Q_S_Mesophyll_i       = Q_S_Mesophyll     - Q_S_Mesophyllbu
        Q_S_ST_i       = Q_S_ST     - Q_S_STbu
        Q_Rm_i        = Q_Rm      - Q_Rmbu
        Q_Gr_i        = Q_Gr      - Q_Grbu

        Q_Exud_i      = Q_Exud    - Q_Exudbu
        Q_meso_i      = Q_meso    - Q_mesobu

        Q_Rmmax_i     = Q_Rmmax   - Q_Rmmaxbu
        Q_Grmax_i     = Q_Grmax   - Q_Grmaxbu

        Q_out_i       = Q_Rm_i    + Q_Exud_i      + Q_Gr_i
        Q_outmax_i    = Q_Rmmax_i + Q_Exud_i   + Q_Grmax_i


        if verbose :
            print("\n\n\n\t\tat ", int(np.floor(simDuration)),"d", int((simDuration%1)*24),"h",  round(r.Qlight *1e6),"mumol m-2 s-1")
            print("Error in Suc_balance:\n\tabs (mmol) {:5.2e}\trel (-) {:5.2e}".format(error, div0f(error,Q_in, 1.)))
            #print("Error in growth:\n\tabs (mmol) {:5.2e}\trel (-) {:5.2e}".format(errorGri, relErrorGri))
            print("Error in photos:\n\tabs (cm3/day) {:5.2e}".format(errLeuning))
            print("water fluxes (cm3/day):\n\ttrans {:5.2e}\tminExud {:5.2e}\tmaxExud {:5.2e}".format(sum(fluxesSoil.values()), min(fluxesSoil.values()), max(fluxesSoil.values())))
            print("C_ST (mmol ml-1):\n\tmean {:.2e}\tmin  {:5.2e} at {:d} segs \tmax  {:5.2e}".format(np.mean(C_ST), min(C_ST), len(np.where(C_ST == min(C_ST) )[0]), max(C_ST)))        
            print("C_me (mmol ml-1):\n\tmean {:.2e}\tmin  {:5.2e}\tmax  {:5.2e}".format(np.mean(C_meso), min(C_meso), max(C_meso)))        
            print('Q_X (mmol Suc): \n\tST   {:.2e}\tmeso {:5.2e}\tin   {:5.2e}'.format(sum(Q_ST), sum(Q_meso), Q_in))
            print('\tRm   {:.2e}\tGr   {:.2e}\tExud {:5.2e}'.format(sum(Q_Rm), sum(Q_Gr), sum(Q_Exud)))
            print("aggregated sink satisfaction at last time step (%) :\n\ttot  {:5.1f}\n\tRm   {:5.1f}\tGr   {:5.1f}".format(
                sum(Q_out_i)/sum(Q_outmax_i)*100,sum(Q_Rm_i)/sum(Q_Rmmax_i)*100, 
                 div0f(sum(Q_Gr_i), sum(Q_Grmax_i), 1.)*100))
            print("aggregated sink repartition at last time step (%) :\n\tRm   {:5.1f}\tGr   {:5.1f}\tExud {:5.1f}".format(sum(Q_Rm_i)/sum(Q_out_i)*100, sum(Q_Gr_i)/sum(Q_out_i)*100,sum(Q_Exud_i)/sum(Q_out_i)*100))
            print("aggregated sink repartition (%) :\n\tRm   {:5.1f}\tGr   {:5.1f}\tExud {:5.1f}".format(sum(Q_Rm)/sum(Q_out)*100, 
                 sum(Q_Gr)/sum(Q_out)*100,sum(Q_Exud)/sum(Q_out)*100))
            print("aggregated sink repartition for max (%) :\n\tRm   {:5.1f}\tGr   {:5.1f}\tExud {:5.1f}".format(sum(Q_Rmmax_i)/sum(Q_outmax_i)*100, 
                 sum(Q_Grmax_i)/sum(Q_outmax_i)*100,sum(Q_Exud_i)/sum(Q_outmax_i)*100))
            print("abs val for max :\n\tRm   {:5.5f}\tGr   {:5.5f}".format(sum(Q_Rmmax_i), 
                 sum(Q_Grmax_i)))
            print("Starch:\n\tQ_S_ST {:5.2e}, C_S_ST {:5.2e}".format(sum(Q_S_ST), np.mean(C_S_ST)))
            print("\tQ_S_Mesophyll {:5.2e}, C_S_Mesophyll {:5.2e}".format(sum(Q_S_Mesophyll), np.mean(C_S_Mesophyll)))
            
        if min(C_ST) < 0.0:
            print("min(C_ST) < 0", min(C_ST),np.mean(C_ST),max(C_ST))
            raise Exception
        
        ots = np.concatenate((np.array([0]), r.get_organ_types()))          
        leavesSegs = np.where(ots[1:] ==4)
        fluxes_leaves = fluxes[leavesSegs]
        if (min(r.Ev) < 0) or (min(r.Jw) < 0) or (min(fluxes_leaves)<0):
            print("leaf looses water", min(r.Ev),min(r.Jw), min(fluxes_leaves))
            raise Exception
                  
        
        if (dtVTP >=0.29/24): #get results as VTP
            getVTPOutput(r, C_ST, np.array(r.psiXyl), vtpIndx, directoryN,  simStartSim, condition)
            dtVTP = 0
        dtVTP += dt
        vtpIndx +=1
        
        dtPrint += dt
        
        if dtPrint >= 0.9/24:
            dtPrint = 0
            write_file_float("time", simDuration, directoryN)
            write_file_array("ots", ots, directoryN)
            write_file_array("C_ST", C_ST, directoryN)
            write_file_array("psiXyl", r.psiXyl, directoryN)
            write_file_array("Q_Rm_dot", Q_Rm_i/dt, directoryN)
            write_file_array("Q_Exud_dot", Q_Exud_i/dt, directoryN)
            write_file_array("Q_Gr_dot", Q_Gr_i/dt, directoryN)
            write_file_array("Q_Rm", Q_Rm, directoryN)
            write_file_array("Q_Exud", Q_Exud, directoryN)
            write_file_array("Q_Gr", Q_Gr, directoryN)
            write_file_array("volOrg",  volOrg, directoryN) #with epsilon
            write_file_array("st_orgs", st_orgs, directoryN)
            write_file_array("ot_orgs", ot_orgs, directoryN)
            write_file_array("len_orgs", lenOrg, directoryN)
            write_file_float("transrate",Evdt, directoryN)
            write_file_array("AgPhl", np.array(r.AgPhl), directoryN)
            write_file_array("fw", r.fw, directoryN)
            
                
        
        
        """ plant growth """
        r.plant.simulate(dt, False)
        
        
        """ for next time step """
        orgs = r.plant.getOrgans(-1, True)
        orgs_all = r.plant.getOrgans(-1, True)
        
        Ntbu = Nt
        Nt = len(r.plant.nodes)

        volOrg = np.array([org.orgVolume(-1,False) for org in orgs_all])
        ot_orgs = np.array([org.organType() for org in orgs_all])
        st_orgs = np.array([org.getParameter("subType") for org in orgs_all])
        lenOrg = np.array([org.getLength(False) for org in orgs_all]) 

        
        Q_Rmbu       =   np.concatenate((Q_Rm, np.full(Nt - Ntbu, 0.)))
        Q_Grbu       =   np.concatenate((Q_Gr, np.full(Nt - Ntbu, 0.))) 
        Q_Exudbu     =   np.concatenate((Q_Exud, np.full(Nt - Ntbu, 0.))) 

        simDuration += dt
        
        
    print("simDuration", simDuration, "d")
    end = datetime.now()
    print(end - beginning)

if __name__ == '__main__':
    scenario = sys.argv[1]
    if scenario == "baseline":
        spellStart = 11 #not used
        condition = "wet"
    elif scenario == "earlyDry":
        spellStart = 11 
        condition = "dry"
    elif scenario == "lateDry":
        spellStart = 18
        condition = "dry"
    else :
        print("scenario", scenario,"not recognised")
        raise Exception
                  
    simInit = 10
    simEnd = 25
    spellDuration =   7
    spellEnd = spellStart + spellDuration
                  
    directoryN = "/"+scenario+"/"
    print(directoryN)

    main_dir=os.environ['PWD']#dir of the file
    results_dir = main_dir +"/results"+directoryN
    
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    else:
        import shutil
        shutil.rmtree(results_dir)
        os.makedirs(results_dir)


    launchUQR(directoryN,simInit,spellStart, condition,spellDuration,simEnd)