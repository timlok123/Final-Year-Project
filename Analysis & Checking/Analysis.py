"""
Analysis.py 

The following groups all the code that are useful for analysis

Lastest Update on 23 Jan 2024 by Justin Chau 
"""

def LoadNumpyFile(FolderName,FileName):
    """
    This function will load the numpy file in Data directory with the name 
    given.

    @param  FolderName(s): the name of the file being stored 
    @param  Name         : the name of the numpy array 
    @return array(np)    : the numpy array loaded
    """
    import os 
    import numpy as np

    current_dir = os.path.dirname(os.path.abspath(__file__))

    # Go up one directory and then navigate to the "Checking h transition" folder
    target_dir = os.path.join(current_dir, '..', 'Data', FolderName)

    return np.load(os.path.join(target_dir, f'{FileName}.npy'))

def LoadNumpyFileMultipleFolder(FileName,*FolderNames):
    """
    This function will load the numpy file in Data directory with the name 
    given.

    @param  FolderName(s): the name of the file being stored 
    @param  Name:       the name of the numpy array 
    @return array(np):  the numpy array loaded
    """
    import os 
    import numpy as np

    current_dir = os.path.dirname(os.path.abspath(__file__))

    # Go up one directory and then navigate to the "Checking h transition" folder
    target_dir = os.path.join(current_dir, '..', 'Data', *FolderNames)

    return np.load(os.path.join(target_dir, f'{FileName}.npy'))

def checking_h_transition_magnetization():
    """
    This function is programmed to check the h_c with |m| value. This function will also 
    plot the value of ED. 

    Steps:
    1. Load the data of Metropolis & Wolff algorithm
    2. Calculate the mean and standard deviation of the graph 
    3. Fitting the data to the curve 
    4. Plot graph and save the plot in "Data\Checking h transition\Absolute m\"
        (You can comment out plt.savefig if you just want to see the plot"
    5. Plot ED value stored in "Data\Checking h transition\Absolute m\ED result"

    """
    import os 
    import numpy as np
    import matplotlib.pyplot as plt

    ED_h_array = np.arange(0,2.05,0.05) 
    Simualtion_h_array=np.arange(0.5,1.6,0.1)

    #1. Load the data 
    Simualtion_data_Metro = LoadNumpyFileMultipleFolder(
                                        "L12_ImL48_tau0.1_10thTrial", 
                                          "Checking h transition",
                                          "Absolute m"
                                          )
    Simualtion_data_Wolff = LoadNumpyFileMultipleFolder(
                                        "L12_ImL48_tau0.1_WolffUpdateBFSListTest_testTrial", 
                                          "Checking h transition",
                                          "Absolute m"
                                          )
    Simualtion_data_Wolff60 = LoadNumpyFileMultipleFolder(
                                        "L12_ImL60_tau0.1_WolffUpdateBFSListTest_10_compare", 
                                          "Checking h transition",
                                          "Absolute m"
                                          )
    Simualtion_data_Wolff72 = LoadNumpyFileMultipleFolder(
                                        "L12_ImL72_tau0.1_WolffUpdateBFSListTest_10_compare", 
                                          "Checking h transition",
                                          "Absolute m"
                                          )
    
    ED_data = LoadNumpyFileMultipleFolder(
                                        "ED_result_average_mz_list", 
                                          "Checking h transition",
                                          "Absolute m",
                                          "ED result"
                                          )
    
    plt.figure(figsize=(10,10))
    """
    plt.errorbar(Simualtion_h_array, 
                 np.mean(Simualtion_data_Metro, axis=1),
                 yerr=np.std(Simualtion_data_Metro, axis=1),
                 fmt='bo-', 
                 ecolor='black',
                 label="L=12,ImL=48 simulation data (Metro)")"""
    
    plt.errorbar(Simualtion_h_array, 
                 np.mean(Simualtion_data_Wolff, axis=1),
                 yerr=np.std(Simualtion_data_Wolff, axis=1),
                 fmt='g^-', 
                 ecolor='green',
                 label="L=12,ImL=48 simulation (Wolff)")
    
    plt.errorbar(Simualtion_h_array, 
                 np.mean(Simualtion_data_Wolff60, axis=1),
                 yerr=np.std(Simualtion_data_Wolff60, axis=1),
                 fmt='b^-', 
                 ecolor='blue',
                 label="L=12,ImL=60 simulation (Wolff)")
    
    plt.errorbar(Simualtion_h_array, 
                 np.mean(Simualtion_data_Wolff72, axis=1),
                 yerr=np.std(Simualtion_data_Wolff72, axis=1),
                 fmt='y^-', 
                 ecolor='yellow',
                 label="L=12,ImL=72 simulation (Wolff)")
    
    #plt.plot(ED_h_array,ED_data,label="ED result")
    
    plt.xlabel("h")
    plt.ylabel("$<|m_z|>$")
    plt.title("<|m_z|> versus h of 1D-transverse field Ising model")
    plt.legend()
    plt.show()

def plot_RE_versus_l_over_L_3_size():

    import numpy as np
    import matplotlib.pyplot as plt 

    #1. Load the data 
    EntropyList = np.zeros((2,4))
    for i in range(1,5):
        #Size8_SubSize4_ImL24_Iteration2000_GI_3rdSpin1Trial
        IndependentToGlued = LoadNumpyFile("VaryingSubSystemSize",
                                           f"Size8_SubSize{i}_ImL24_Iteration2000_IG_4thSpin1Trial")
        GluedToIndependent = LoadNumpyFile("VaryingSubSystemSize",
                                           f"Size8_SubSize{i}_ImL24_Iteration2000_GI_4thSpin1Trial")

        print()
        #Calculate the mean and error 
        #print(f"Connected Site {i},")
        #print(f"IG means is :{np.mean(IndependentToGlued)}")
        #print(f"GI means is :{np.mean(GluedToIndependent)}")

        EntropyList[0][i-1] = -np.log(np.mean(IndependentToGlued)/np.mean(GluedToIndependent))
        EntropyList[1][i-1] = np.sqrt(np.abs(np.std(GluedToIndependent)/np.mean(GluedToIndependent) 
                             - np.std(IndependentToGlued)/np.mean(IndependentToGlued)))

    
    EntropyListL8Im36 = np.zeros((2,4))
    for i in range(1,5):

        IndependentToGlued = LoadNumpyFile("VaryingSubSystemSize",
                                           f"Size8_SubSize{i}_ImL36_Iteration2000_IG_4thSpin1Trial")
        GluedToIndependent = LoadNumpyFile("VaryingSubSystemSize",
                                           f"Size8_SubSize{i}_ImL36_Iteration2000_GI_4thSpin1Trial")

        #Calculate the mean and error 
        EntropyListL8Im36[0][i-1] = -np.log(np.mean(IndependentToGlued)/np.mean(GluedToIndependent))
        EntropyListL8Im36[1][i-1] = np.sqrt(np.abs(np.std(GluedToIndependent)/np.mean(GluedToIndependent) 
                                                   - np.std(IndependentToGlued)/np.mean(IndependentToGlued)))
        
    
    EntropyListL8Im48 = np.zeros((2,4))
    for i in range(1,5):

        IndependentToGlued = LoadNumpyFile("VaryingSubSystemSize",
                                           f"Size8_SubSize{i}_ImL48_Iteration2000_IG_4thSpin1Trial")
        GluedToIndependent = LoadNumpyFile("VaryingSubSystemSize",
                                           f"Size8_SubSize{i}_ImL48_Iteration2000_GI_4thSpin1Trial")

        #Calculate the mean and error 
        EntropyListL8Im48[0][i-1] = -np.log(np.mean(IndependentToGlued)/np.mean(GluedToIndependent))
        EntropyListL8Im48[1][i-1] = np.sqrt(np.abs(np.std(GluedToIndependent)/np.mean(GluedToIndependent) 
                                                   - np.std(IndependentToGlued)/np.mean(IndependentToGlued)))
        
    
    plt.errorbar([i/8 for i in range(1,5)],
                 EntropyList[0],
                 yerr=EntropyList[1],
                 label="L=8,ImL=24")
    print(EntropyList[0])

    """
    plt.errorbar([i/8 for i in range(1,5)],
                 EntropyListL8Im36[0],
                 yerr=EntropyListL8Im36[1],
                 label="L=8,ImL=36")

    plt.errorbar([i/8 for i in range(1,5)],
                 EntropyListL8Im48[0],
                 yerr=EntropyListL8Im48[1],
                 label="L=8,ImL=48")
    """
    
    plt.xlabel("(Subsystem Size)/(System Size)",size=15)
    plt.ylabel("Renyi Entropy",size=15)
    plt.legend()
    plt.xlim(0,1)
    plt.show()

def checking_h_transition_EnergyPerSite():
    """
    This function is programmed to check the h_c with E/N value 

    Steps:
    1. Load the data 
    2. Calculate the mean and standard deviation of the graph 
    3. Fitting the data to the curve 
    4. Plot graph and save the plot in Data/'Checking h transition'
        (You can comment out plt.savefig if you just want to see the plot)

    """

    import os 
    import numpy as np
    import matplotlib.pyplot as plt 

    discrete_h = np.arange(0.3,1.9,0.1)

    #1. Load the data 
    data_array  = LoadNumpyFileMultipleFolder(
        "L16_ImL32_tau0.1_Energy_3rdTrial",
        "Checking h transition",
        "Energy Per Site")
    
    #2. Calculate the mean 
    mean_data  = np.mean(data_array,axis=1)

    #3. Fitting
    from scipy.optimize import curve_fit
    popt       = curve_fit(QuadraticFittingFunction,discrete_h,mean_data)[0]
    popt_cubic = curve_fit(CubicFittingFunction,discrete_h,mean_data)[0]

    # Calculate the slope 
    # FDD 
    slope_FDD = np.zeros(len(mean_data)-1)
    for i in range(0,len(mean_data)-1):
        slope_FDD[i] = (mean_data[i+1] - mean_data[i])/(0.1)
    
    print(mean_data)
    print(slope_FDD)

    #4. Graph plotting 
    plot_h = np.linspace(0.1,2.3,1000)

    slopeQuad = lambda h: 2*popt[0]*h + popt[1]
    slopeCubic = lambda h: 3*popt[0]*h**2 + 2*popt[1]*h + popt[2]

    plt.plot(plot_h,QuadraticFittingFunction(plot_h,*popt),label="(Quadratic fitting)")
    plt.plot(plot_h,CubicFittingFunction(plot_h,*popt_cubic),label="(Cubic fitting)")
    #plt.plot(plot_h,slopeQuad(plot_h),label="Slope of quadratic fitting")
    #plt.plot(plot_h,slopeCubic(plot_h),label="Slope of cubic fitting")

    plt.errorbar(discrete_h, 
                 mean_data,
                 yerr=np.std(data_array , axis=1),
                 fmt='o', 
                 ecolor='black',
                 label="L=4 (simulation)")

    plt.legend(prop={'size': 10})
    plt.xlabel("h",fontsize=15)
    plt.ylabel("Energy per site",fontsize=15)

    current_dir = os.path.dirname(os.path.abspath(__file__))
    target_dir = os.path.join(current_dir, '..', 'Data', 'Checking h transition','Energy Per Site')
    #plt.savefig(os.path.join(target_dir,"Energy_per_site_versus_h.png"))

    plt.show()
 
def checking_h_transition_BinderCumlant():
    """
    This function is programmed to check the h_c with Binder Cumlant

    Steps:
    1. Load the data 
    2. Calculate the mean and standard deviation of the graph 
    #3. Fitting the data to the curve 
    4. Plot graph and save the plot in Data/'Checking h transition'
        (You can comment out plt.savefig if you just want to see the plot)

    """
    import os 
    import numpy as np
    import matplotlib.pyplot as plt 
    
    discrete_h = np.arange(0.96,1.04,0.02)

    #1. Load the data 
    data_L4  = LoadNumpyFileMultipleFolder("L4_ImL40_tau0.1_7thTrial","Checking h transition","Absolute m")
    data_L8  = LoadNumpyFileMultipleFolder("L8_ImL40_tau0.1_7thTrial","Checking h transition","Absolute m")
    data_L16 = LoadNumpyFileMultipleFolder("L16_ImL40_tau0.1_7thTrial","Checking h transition","Absolute m")

    #2. Calculate the <|m|> and <m^2>
    mean_square_m_L4 = np.mean(np.square(data_L4),axis=1)
    mean_abs_m_L4 =  np.mean(data_L4,axis=1)

    mean_square_m_L8 = np.mean(np.square(data_L8),axis=1)
    mean_abs_m_L8 =  np.mean(data_L8,axis=1)

    BinderCumulant_L4 = 3/2*(1- 1/3*(mean_square_m_L4/(mean_abs_m_L4)**2))
    BinderCumulant_L8 = 3/2*(1- 1/3*(mean_square_m_L8/(mean_abs_m_L8)**2))

    #4. Plot the graph 

    """
    plt.errorbar(discrete_h, 
                 BinderCumulant_L4,
                 yerr=np.std(data_L4, axis=1),
                 fmt='o', 
                 ecolor='black',
                 label="L=4 (simulation)")
    """
    
    plt.scatter(discrete_h, 
                 BinderCumulant_L4,
                 color='black',
                 label="L=4 (simulation)")

    plt.scatter(discrete_h, 
                 BinderCumulant_L8,
                 color='green',
                 label="L=8 (simulation)")
    
    plt.ylabel("Binder Cumulant")
    plt.xlabel("h")
    plt.legend()
    plt.show()

def plot_RE_versus_l_over_L(
        L=8,
        ImL= 80,
        NoOfBin=10,
        TrialName="new"
        ):

    import numpy as np
    import matplotlib.pyplot as plt 

    #1. Load the data 
    EntropyList = np.zeros((2,int(L/2)))
    for i in range(0,int(L/2)):
        IndependentToGlued = LoadNumpyFile("VaryingSubSystemSize",
                                           #Size8_SubSize4_ImL80_IG_testWithJackTrial
                                           f"Size{L}_SubSize{i+1}_ImL{ImL}_IG_{TrialName}")
        GluedToIndependent = LoadNumpyFile("VaryingSubSystemSize",
                                           f"Size{L}_SubSize{i+1}_ImL{ImL}_GI_{TrialName}")
        
        # Calculate the mean and error 
        EntropyList[0][i] = -np.log(np.mean(IndependentToGlued)/np.mean(GluedToIndependent))
        EntropyList[1][i] = np.sqrt(np.abs(np.std(GluedToIndependent)/np.mean(GluedToIndependent) 
                             - np.std(IndependentToGlued)/np.mean(IndependentToGlued)))


    # 3. Plot the graph
    plt.errorbar(np.arange(1,int(L/2)+1)/(L),
                 EntropyList[0],
                 yerr=EntropyList[1]/np.sqrt(NoOfBin),
                 label=f"L={L},ImL={ImL} Simulation data")
    plt.plot(np.arange(1,int(L/2)+1)/(L),[0.3704,0.4710,0.5232,0.5525],label="ED result")
    print(EntropyList[0])
    print(np.arange(1,int(L/2)+1)/(L))

    plt.xlabel("(Subsystem Size)/(System Size)",size=15)
    plt.ylabel("Renyi Entropy",size=15)
    plt.legend()
    plt.xlim(0,1)
    plt.show()

def FittingFunctionFiniteTemperature(SystemSize,CentralCharge,ConstantAtBack):
    import numpy as np
    beta = SystemSize
    return (CentralCharge/6)*(1+1/2)*np.log(SystemSize/np.pi*np.sinh(np.pi/2)) + ConstantAtBack

def FittingFunctionZeroTemperature(SubSizeOverSystemSize,ConstantAtBack):
    import numpy as np
    SystemSize = 12 
    return (0.5/3)*(1+1/2)*np.log(SystemSize/np.pi*np.sin(np.pi*SubSizeOverSystemSize)) + ConstantAtBack

def temp():
    import numpy as np 
    from scipy.optimize import curve_fit
    import matplotlib.pyplot as plt 
    L = 32

    """
    IndependentToGlued = LoadNumpyFile("RenyiEntropyCalculation",
                                       f"Size{L}_SubSize{int(L/2)}_ImL{L*10}_IG_WolffUpdateBFSList_test1")
    GluedToIndependent = LoadNumpyFile("RenyiEntropyCalculation",
                                        f"Size{L}_SubSize{int(L/2)}_ImL{L*10}_GI_WolffUpdateBFSList_test1")
    
    #Calculate the mean and error 
    mean = -np.log(np.mean(IndependentToGlued)/np.mean(GluedToIndependent))
    error = np.sqrt(np.abs(np.std(GluedToIndependent)/np.mean(GluedToIndependent) 
                    - np.std(IndependentToGlued)/np.mean(IndependentToGlued)))/np.sqrt(20)
    
    print(f"For system size {L},")
    print(f"Mean is : {mean}")
    print(f"Error is : {error}")"""

    # Fitting 
    Size = [4,8,12,16,32]
    RenyiEntropyList = [0.5370265337472249,
                        0.6082755475232792,
                        0.6474785609459838,
                        0.6674295389875499,
                        0.707323166753652]
    popt, pcov = curve_fit(FittingFunction, Size, RenyiEntropyList)
    popt2, pcov2 = curve_fit(FittingFunctionFiniteTemperature,Size,RenyiEntropyList)

    # Generate x values for the fit function
    x_fit = np.linspace(min(Size), max(Size), 1000)
    
    # Generate y values for the fit function
    y_fit = FittingFunction(x_fit, *popt)
    y_fit2 = FittingFunctionFiniteTemperature(x_fit, *popt2)
    
    # Plot the data points
    plt.scatter(Size, RenyiEntropyList, label='Data')
    
    # Plot the fit function
    #plt.plot(x_fit, y_fit, 'r')
    plt.plot(x_fit, y_fit, 'r', label='Fit: a=%5.3f, b=%5.3f' % tuple(popt))
    plt.plot(x_fit, y_fit2, 'r', label='Fit: a=%5.3f, b=%5.3f' % tuple(popt2))

    print(f"a = {popt[0]}, b = {popt[1]}")
    print(f"a = {popt2[0]}, b = {popt2[1]}")
    
    plt.xlabel('Size')
    plt.ylabel('Renyi Entropy')
    plt.legend()
    plt.show()

def temp2():
    import numpy as np 
    from scipy.optimize import curve_fit
    import matplotlib.pyplot as plt 

    L = 8
    A = 4

    IndependentToGlued = LoadNumpyFile("RenyiEntropyCalculation",
                                       f"Size{L}_SubSize{A}_ImL{L*10}_IG_WolffUpdateBFSList_test2")
    GluedToIndependent = LoadNumpyFile("RenyiEntropyCalculation",
                                        f"Size{L}_SubSize{A}_ImL{L*10}_GI_WolffUpdateBFSList_test2")
    
    #Calculate the mean and error 
    mean = -np.log(np.mean(IndependentToGlued)/np.mean(GluedToIndependent))
    error = np.sqrt(np.abs(np.std(GluedToIndependent)/np.mean(GluedToIndependent) 
                    - np.std(IndependentToGlued)/np.mean(IndependentToGlued)))/np.sqrt(20)
    
    print(f"For system size {L},subsystem size is {A},")
    print(f"Mean is : {mean}")
    print(f"Error is : {error}")

    Size = [1,2,3,4]
    RenyiEntropyList = [0.382,0.49,0.552,0.596] 

def temp3():
    import matplotlib.pyplot as plt 
    import numpy as np
    from scipy.optimize import curve_fit

    # 1. Import the data 
    L = 32
    EE_simulation = np.array([0.35989698,0.46442355,0.50605405,0.57493231,0.59673231,0.62534444,
  0.6288317,0.64556529,0.66469825,0.70768307,0.68851669,0.73452671,
  0.70935035,0.73668448,0.75009631,0.71653473])
    EE_error = np.array([0.01004245,0.01403966,0.01300118,0.01955747,0.01542957,0.02351392,
  0.02794488,0.0152281 ,0.04193116,0.02521341,0.02545573,0.02768648,
  0.0330812,0.05010195,0.01755952,0.03231457])
    SubSystemSize = np.arange(1,17)
    sizeplot = np.linspace(0,0.98)

    #2. Fitting to the formula 
    popt_simulation, pcov = curve_fit(FittingFunctionZeroTemperature, SubSystemSize/L, EE_simulation)

    print(popt_simulation)
    plt.plot(sizeplot,FittingFunctionZeroTemperature(sizeplot,*popt_simulation),label="Fitting")
    plt.plot(SubSystemSize/L,EE_simulation,label="Simulation result (Wolff)")
    plt.legend()
    plt.xlim(0,1)
    #plt.xscale("log")
    plt.show()

    SubSystemSizePlot = np.linspace(SubSystemSize[0],SubSystemSize[-1],1000)
    plt.title("L=32")
    plt.errorbar(np.log((L/np.pi)*np.sin(np.pi*SubSystemSize/L)),
             EE_simulation,
             yerr=EE_error,
             label="Simulation result (Wolff)")
    plt.plot(np.log(L/np.pi*np.sin(np.pi*SubSystemSizePlot/L)),
                    FittingFunctionZeroTemperature(SubSystemSizePlot/L,*popt_simulation),
                    label="Fitting - fix c = 1/2")
    plt.legend()
    plt.show()

def RenyiEntaglementEntropyZeroTemperatureOneHalf(SystemSize,CentralCharge,ConstantAtBack):
    """
    The formula is 
    S = c/(3*eta)*(1+1/alpha)*log(eta*L/pi*a * sin(pi*x/L))
    
    c     - central charge 
    eta   - eta is set to 1 for periodic boundary
    alpha - RenyiEntaglement Entropy order, set it to 2 here. 
    a     - the spacing between the site. Set it to 1 
    L     - System size
    x     - Subsystem size

    Since I set x/L = 1/2, the sin term at the back would be 1. 
    """

    import numpy as np
    return (CentralCharge/4)*np.log((SystemSize)/(np.pi*1)) + ConstantAtBack

def RenyiEntaglementEntropyZeroTemperatureOneHalf2nd(SystemSize,ConstantAtBack):
    """
    The formula is 
    S = c/(3*eta)*(1+1/alpha)*log(eta*L/pi*a * sin(pi*x/L))
    
    c     - central charge 
    eta   - eta is set to 1 for periodic boundary
    alpha - RenyiEntaglement Entropy order, set it to 2 here. 
    a     - the spacing between the site. Set it to 1 
    L     - System size
    x     - Subsystem size

    Since I set x/L = 1/2, the sin term at the back would be 1. 
    """

    import numpy as np
    return (0.5/4)*np.log((SystemSize)/(np.pi*1)) + ConstantAtBack

def RenyiEntaglementEntropyZeroTemperatureOneHalfTemp(SubsystemSizeOverSystemSize,ConstantAtBack):
    """
    
    c     - central charge 
    eta   - eta is set to 1 for periodic boundary
    alpha - RenyiEntaglement Entropy order, set it to 2 here. 
    a     - the spacing between the site. Set it to 1 
    L     - System size
    x     - Subsystem size

    Since I set x/L = 1/2, the sin term at the back would be 1. 
    """

    import numpy as np
    return (0.5/4)*np.log((SystemSize)/(np.pi*1)) + ConstantAtBack

def Fitting_Renyi_Entropy_On_Different_System_Size(
        SavePlot=False,
        PlotInLogScale = True):
    """
    This function is to load the data on Renyi Entropy of different system size 
    and obtain central charge through fitting. A graph of RRE versus SystemSize will be 
    plot as well 

    I will use Renyi Entanglement Entropy formula for zero temperature for fitting
    in this function. Refer to function RenyiEntaglementEntropyZeroTemperatureOneHalf(). 

    """

    import matplotlib.pyplot as plt 
    import numpy as np
    from scipy.optimize import curve_fit

    # 1. Loading data
    SystemSize = [4,8,12,16,32]
    Data_array = np.zeros((2,len(SystemSize)))
    
    for count,L in enumerate(SystemSize):
        IndependentToGlued = LoadNumpyFileMultipleFolder(
                             f"Size{L}_SubSize{int(L/2)}_ImL{L*10}_IG_WolffUpdateBFSList_test",
                            "RenyiEntropyCalculation",
                            "Different system size")
        GluedToIndependent = LoadNumpyFileMultipleFolder(
                             f"Size{L}_SubSize{int(L/2)}_ImL{L*10}_GI_WolffUpdateBFSList_test",
                            "RenyiEntropyCalculation",
                            "Different system size")
        
        Data_array[0][count] = -np.log(np.mean(IndependentToGlued)/np.mean(GluedToIndependent))
        Data_array[1][count] = np.sqrt(np.abs(np.std(GluedToIndependent)/np.mean(GluedToIndependent) 
                               - np.std(IndependentToGlued)/np.mean(IndependentToGlued)))/np.sqrt(20)

    # 2. Fitting
    popt_simulation, pcov = curve_fit(RenyiEntaglementEntropyZeroTemperatureOneHalf, 
                                      SystemSize,Data_array[0])
    popt_simulation2, pcov = curve_fit(RenyiEntaglementEntropyZeroTemperatureOneHalf2nd, 
                                      SystemSize,Data_array[0])
    
    # 3. Plot the graph 
    SystemSizePlot = np.linspace(min(SystemSize), max(SystemSize), 1000)
    print(popt_simulation)

    plt.plot(SystemSizePlot,
             RenyiEntaglementEntropyZeroTemperatureOneHalf(SystemSizePlot,*popt_simulation),
             label = f"Fitting curve. Central charge = {popt_simulation[0]}")
    plt.plot(SystemSizePlot,
             RenyiEntaglementEntropyZeroTemperatureOneHalf2nd(SystemSizePlot,*popt_simulation2),
             label = f"Fitting curve (fixed central charge as 1/2).")
    
    plt.errorbar(SystemSize, 
                 Data_array[0],
                 yerr=Data_array[1],
                 fmt='g^-', 
                 ecolor='green',
                 label="Simulation data (Wolff)")
    
    plt.title("Second Order Renyi Entanglement Entropy on different system size")
    plt.xlabel("System Size(L)")
    
    plt.ylabel(r"Second Order Renyi Entanglement Entropy($S_2$)")
    plt.legend()

    #if PlotInLogScale:
    if True:
        plt.xscale("log")

    if SavePlot:
        plt.savefig("Second Order Renyi Entanglement Entropy on different system size")
    
    plt.show()

    #Implement the save function later

def Comparing_the_ED_and_Simulation_result(
        SavePlot=False
        ):
    """
    This function is to load the data on Renyi Entropy of L=8, A = 1,2,3,4 of Exact 
    diagonalization and simulation result. It is used to validate my Wolff algorithm.

    The no of bin is 50, each bin = 10000 update. 
    """
    import matplotlib.pyplot as plt 
    import numpy as np
    from scipy.optimize import curve_fit

    # 1. Loading data (simulation)
    # The constant is set here 
    L=8
    NoOfBins=50 

    SubSystemSize = np.arange(1,5)
    SimulationDataArray = np.zeros((2,len(SubSystemSize )))
    
    for count,A in enumerate(SubSystemSize):
        IndependentToGlued = LoadNumpyFileMultipleFolder(
                             f"Size{L}_SubSize{A}_ImL{L*10}_IG_WolffUpdateBFSList_test",
                            "RenyiEntropyCalculation",
                            "SameSystemSizeL8")
        GluedToIndependent = LoadNumpyFileMultipleFolder(
                             f"Size{L}_SubSize{A}_ImL{L*10}_GI_WolffUpdateBFSList_test",
                            "RenyiEntropyCalculation",
                            "SameSystemSizeL8")
    
        SimulationDataArray[0][count] = -np.log(np.mean(IndependentToGlued)/np.mean(GluedToIndependent))
        SimulationDataArray[1][count] = np.sqrt(np.abs(np.std(GluedToIndependent)/np.mean(GluedToIndependent) 
                               - np.std(IndependentToGlued)/np.mean(IndependentToGlued)))/np.sqrt(NoOfBins)
    
    EDData = [0.3704,0.4710,0.5232,0.5525] # From Leo Wang's code 

    #2. Plot the data  
    plt.errorbar(SubSystemSize/L, 
                 SimulationDataArray[0],
                 yerr=SimulationDataArray[1],
                 fmt='g^-', 
                 ecolor='green',
                 label="Simulation data (Wolff)")
    
    plt.scatter(SubSystemSize/L,
                EDData,
                label="ED result")
    
    plt.title("Second Order Renyi Entanglement Entropy - L = 8")
    plt.xlabel("SubsystemSize/System Size(l/L)")
    plt.ylabel(r"Second Order Renyi Entanglement Entropy($S_2$)")
    plt.legend()

    if SavePlot:
        plt.savefig("Second Order Renyi Entanglement Entropy - L = 8")
    
    plt.show()


    return 0 

def plot_RE_versus_l_over_L_Size_8_12(
        SavePlot=False
        ):

    import numpy as np
    import matplotlib.pyplot as plt 

    NoOfBin = 20

    #1. Load the data of L = 8 
    EntropyListL8 = np.zeros((2,4))
    for count,A in enumerate(range(1,int(8/2)+1)):
        IndependentToGlued = LoadNumpyFileMultipleFolder(
                             f"Size8_SubSize{A}_ImL80_IG_WolffUpdateBFSList_test",
                            "RenyiEntropyCalculation",
                            "SameSystemSizeL8")

        GluedToIndependent = LoadNumpyFileMultipleFolder(
                             f"Size8_SubSize{A}_ImL80_GI_WolffUpdateBFSList_test",
                            "RenyiEntropyCalculation",
                            "SameSystemSizeL8")

        
        # Calculate the mean and error 
        EntropyListL8[0][count] = -np.log(np.mean(IndependentToGlued)/np.mean(GluedToIndependent))
        EntropyListL8[1][count] = np.sqrt(np.abs(np.std(GluedToIndependent)/np.mean(GluedToIndependent) 
                             - np.std(IndependentToGlued)/np.mean(IndependentToGlued)))/np.sqrt(NoOfBin)
    
    # 2. Load the data of L = 12 
    EntropyListL12 = np.zeros((2,6))       
    for count,A in enumerate(range(1,int(12/2)+1)):
        IndependentToGlued = LoadNumpyFileMultipleFolder(
                             f"Size12_SubSize{A}_ImL120_IG_WolffUpdateBFSList_30_test",
                            "RenyiEntropyCalculation",
                            "SameSystemSizeL12")

        GluedToIndependent = LoadNumpyFileMultipleFolder(
                             f"Size12_SubSize{A}_ImL120_GI_WolffUpdateBFSList_30_test",
                            "RenyiEntropyCalculation",
                            "SameSystemSizeL12")
        
        # Calculate the mean and error 
        EntropyListL12[0][count] = -np.log(np.mean(IndependentToGlued)/np.mean(GluedToIndependent))
        EntropyListL12[1][count] = np.sqrt(np.abs(np.std(GluedToIndependent)/np.mean(GluedToIndependent) 
                             - np.std(IndependentToGlued)/np.mean(IndependentToGlued)))/np.sqrt(NoOfBin)
        

    print(EntropyListL12)
    # 3. Plot the graph
    plt.errorbar(np.arange(1,int(8/2)+1)/(8),
                 EntropyListL8[0],
                 yerr=EntropyListL8[1],
                 fmt="r*--",
                 ecolor="red",
                 label=f"L=8,ImL=80 Simulation data(Wolff)")
    plt.plot(np.arange(1,int(8/2)+1)/(8),[0.3704,0.4710,0.5232,0.5525],label="L=8, ED result")

    plt.errorbar(np.arange(1,int(12/2)+1)/(12),
             EntropyListL12[0],
             yerr=EntropyListL12[1],
             fmt="g^-",
             ecolor="green",
             label=f"L=12,ImL=120 Simulation data(Wolff)")

    plt.xlabel("(Subsystem Size)/(System Size)",size=15)
    plt.ylabel("Renyi Entropy",size=15)
    plt.legend()
    plt.xlim(0,1)
    if SavePlot:
        plt.savefig()
    plt.show()

def plot_RE_versus_l_over_L_Size_8_12_modify(
        SavePlot=False
        ):

    import numpy as np
    import matplotlib.pyplot as plt 

    NoOfBin = 20
    """
    # 2. Load the data of L = 12 
    EntropyListL12 = np.zeros((2,6))       
    for count,A in enumerate(range(1,int(12/2)+1)):
        IndependentToGlued = LoadNumpyFileMultipleFolder(
                             f"Size12_SubSize{A}_ImL120_IG_WolffUpdateBFSList_30_test",
                            "RenyiEntropyCalculation",
                            "SameSystemSizeL12")

        GluedToIndependent = LoadNumpyFileMultipleFolder(
                             f"Size12_SubSize{A}_ImL120_GI_WolffUpdateBFSList_30_test",
                            "RenyiEntropyCalculation",
                            "SameSystemSizeL12")
        
        # Calculate the mean and error 
        EntropyListL12[0][count] = -np.log(np.mean(IndependentToGlued)/np.mean(GluedToIndependent))
        EntropyListL12[1][count] = np.sqrt(np.abs(np.std(GluedToIndependent)/np.mean(GluedToIndependent) 
                             - np.std(IndependentToGlued)/np.mean(IndependentToGlued)))/np.sqrt(NoOfBin)"""
    
    EntropyListL12 = LoadRenyiEntanglementEntropyDataWolffHalfSystemSize(
                    FolderName="SameSystemSizeL12",
                    SystemSize=12,
                    NoOfBin = 30)

   
    # 3. Plot the graph
    plt.errorbar(np.arange(1,int(12/2)+1)/(12),
             EntropyListL12[0],
             yerr=EntropyListL12[1],
             fmt="g^-",
             ecolor="green",
             label=f"L=12,ImL=120 Simulation data(Wolff)")

    plt.xlabel("(Subsystem Size)/(System Size)",size=15)
    plt.ylabel("Renyi Entropy",size=15)
    plt.legend()
    plt.xlim(0,1)
    if SavePlot:
        plt.savefig()
    plt.show()

def LoadRenyiEntanglementEntropyDataWolffHalfSystemSize(
        FolderName="SameSystemSizeL12",
        SystemSize=12,
        NoOfBin = 30 
):
    """
    (fill out the documentation later)
    
    """

    import numpy as np
    import matplotlib.pyplot as plt 

    #NoOfBin = 30
    
    EntropyArray = np.zeros((2,int(SystemSize/2)))       
    for count,A in enumerate(range(1,int(SystemSize/2)+1)):

        try: 
            IndependentToGlued = LoadNumpyFileMultipleFolder(
                                 f"Size{SystemSize}_SubSize{A}_ImL{SystemSize*10}_IG_WolffUpdateBFSList_{NoOfBin}_test",
                                "RenyiEntropyCalculation",
                                FolderName)
    
            GluedToIndependent = LoadNumpyFileMultipleFolder(
                                 f"Size{SystemSize}_SubSize{A}_ImL{SystemSize*10}_GI_WolffUpdateBFSList_{NoOfBin}_test",
                                "RenyiEntropyCalculation",
                                FolderName)
            
            # Calculate the mean and error 
            EntropyArray[0][count] = -np.log(np.mean(IndependentToGlued)/np.mean(GluedToIndependent))
            EntropyArray[1][count] = np.sqrt(np.abs(np.std(GluedToIndependent)/np.mean(GluedToIndependent) 
                                 - np.std(IndependentToGlued)/np.mean(IndependentToGlued)))/np.sqrt(NoOfBin)
            
        except:
            print("Size{SystemSize}_SubSize{A}_ImL{SystemSize*10}_GI_WolffUpdateBFSList_{NoOfBin}_test is missing")
        


    return EntropyArray


def main():

    temp3()
    
    return 0 


if __name__ == '__main__':
    import time
    start_time = time.time()
    main() 
    end_time = time.time()
    print(f"The code took {(end_time - start_time)/60} min(s) to run.")