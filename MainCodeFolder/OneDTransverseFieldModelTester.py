from OneDTransverseFieldModel import *

def test_MetropolisUpdateNtimes():
    """
    This function is used to test the method MetropolisUpdateNtimes in OneDTransverseFieldModel
    I set a 4x4 system, update it 10 times, and see whether 
    1. Correct index obtain? 
    2. Correct delta E? 
    3. Correct flipping?
    """
    iteration = 100 
    L = 4 
    Im_L = 4 
    a = OneDTFIM(4,4)
    #a.MetropolisUpdateNtimes(10) 

    print(a.TFIM)

def testInitializeNeighbour():
    ImL = 3 
    SysL = 4 
    a = OneDTFIM(SysL,ImL)
    a.InitializeNeigbour()
    #print(a.ImaginaryTimeAxisLength)
    print(a.neighbour[1][1])
    #print(a.neighbour[3][2])

    for count,s in enumerate(a.neighbour[1][1]):
        x,y = s 
        print("count is "+str(count))
        print("x is {x}".format(x=x))
        #print("y is " +y)
    
    b = 2 
    if (b == 1) and  \
       (b == 2):
        print("Yes")
    
    else:
        print("No")

def arrayOperationtest():
    a = [1,2]
    d = [0,0]
    b = np.array([[1,2],[3,4]])
    c = np.zeros((2*2,2))
    c[0] = a 
    c[1] = d
    
    #print(np.max(np.nonzero(c))) # this could return the last non-zero elements 
    #print(a in b)
    #print(type(b[0]))

    PointToTop = 2 
    ImL = 4 
    ReL = 3
    Pocket = np.zeros((ImL*ReL,2))
    PocketPointToTop = 3
    Pocket[0] = [1,1]
    Pocket[1] = [2,3]
    Pocket[2] = [2,4]
    randomIndex = np.random.randint(0,PocketPointToTop)
    #print(randomIndex)
    s = Pocket[randomIndex]
    #print(Pocket)
    #print(s)

    a = np.array(
        [[1,2,3,4],
         [5,6,7,8],
         [9,10,11,12]]
    )

    row_averages = np.mean(a, axis=1)
    print(row_averages)

def testWolffUpdate():

    Real, Im = 4,4
    a = OneDTFIM(Real,Im)
    a.InitializeNeigbour()

    print("Configuration before update: ")
    print(a)
    
    a.WolffUpdateBFSList()

    print("Configuration after update: ")
    print(a)

def test_GetAbsolutemagentization():

    Real, Im = 4,4
    a = OneDTFIM(Real,Im)
    a.get_variable("TFIM")
    print(a.GetAbsolutemagentization())

    return 0 

def test_get_set_h():

    Real, Im = 4,4
    a = OneDTFIM(Real,Im)
    print(f"The magnetic strength is: {a.Geth_MagneticFieldStrength()}")

    a.Seth_MagneticFieldStrength(0.8)
    print(f"The magnetic strength is: {a.Geth_MagneticFieldStrength()}")
    
    a.Seth_MagneticFieldStrength(0.9)
    print(f"The magnetic strength is: {a.Geth_MagneticFieldStrength()}")

    a.Seth_MagneticFieldStrength(1.1)
    print(f"The magnetic strength is: {a.Geth_MagneticFieldStrength()}")

    return 0 

def CheckPhaseTransitionPoint_abs_m(
        Real=12,
        Im=48,
        thermalization = 5000,
        iteration = 10000,
        bin = 1000,
        h_array=np.arange(0.5,1.6,0.1), 
        UpdateMethod = "MetropolisUpdate",
        NoTrial = "new"
        ):
    """
    This function will plot |m| versus h graph, with an aim to find 
    the phase transition point. The data and graph will be saved to 
    Data folder 
    """
    import os 
    
    print(f"Real :{Real}, Im :{Im} simulation")
    print(f"The updating methods used is {UpdateMethod}")

    TFIM = OneDTFIM(Real,Im)
    TFIM.SetDetlaTau(0.1)
    data_array = np.zeros((len(h_array),int(iteration/bin)))
    NoOfBin = int(iteration/bin)
    
    # Simulation 
    for index, h in enumerate(h_array):

        TFIM.Seth_MagneticFieldStrength(h)

        #1. Thermalization 
        for i in range(thermalization):
            getattr(TFIM, UpdateMethod)()
        
        print(f"Completed Thermalization at h = {TFIM.Geth_MagneticFieldStrength()}")

        #2. Collect data
        for i in range(int(iteration/bin)):
            absolute_m_list = []

            for j in range(bin):
                getattr(TFIM, UpdateMethod)()
                absolute_m_list.append(TFIM.GetAbsolutemagnetization())
            
            #3. Bin and take sd & mean 
            data_array[index][i] = np.mean(absolute_m_list)
        
        print(f" h = {TFIM.Geth_MagneticFieldStrength()} is completed.")
        print(f" Mean is {np.mean(data_array[index])}")
        print(f" Sd is {np.std(data_array[index])}")
    
    # Plot graph
    plt.figure(figsize=(10,10))
    plt.errorbar(h_array, 
                 np.mean(data_array, axis=1),
                 yerr=np.std(data_array, axis=1)/np.sqrt(NoOfBin),
                 fmt='o', 
                 ecolor='black')
    
    plt.xlabel("h")
    plt.ylabel("|m|")

    # Saving plot and data 
    current_dir = os.path.dirname(os.path.abspath(__file__))
    target_dir = os.path.join(current_dir, '..', 'Data', 'Checking h transition','Absolute m')
    plt.savefig(os.path.join(target_dir,
                             f"L{Real}_ImL{Im}_tau0.1_{UpdateMethod}_{NoOfBin}_{NoTrial}.png"))
    np.save(os.path.join(target_dir, 
                         f'L{Real}_ImL{Im}_tau0.1_{UpdateMethod}_{NoOfBin}_{NoTrial}.npy'),data_array)

def test_savingNumpyArray():
    import os

    a = np.array([1,2,3]) 
    name1 = "TestArray1"

    # Get the current file's directory
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # Go up one directory and then navigate to the "Checking h transition" folder
    target_dir = os.path.join(current_dir, '..', 'Data', 'Checking h transition')

    # Save the numpy array in the target directory
    np.save(os.path.join(target_dir, f'{name1}.npy'), a)

    print("The data is saved")

def test_GetEnergyPerSite():

    Real, Im = 3,3
    a = OneDTFIM(Real,Im)
    #a.__InitializeNeigbour__()
    a.get_variable("TFIM")
    a.get_variable("J_x")
    a.get_variable("h")
    a.get_variable("J_y")
    print(a.GetEnergyPerSite())
    
    return 0      

def CheckPhaseTransitionPoint_EnergyPerSite(
        Real=4,
        Im=40,
        h_array=np.arange(0.1,2.2,0.1),
        NoTrial="new"):
    """
    This function will plot Energy per sites (E/N) versus h graph, 
    with an aim to find the phase transition point. The data and graph 
    will be saved to the Data folder 
    """
    import os 
    
    print("This is Energy Versus h simulation of 1D transverse field Ising model")
    print(f"Real :{Real}, Im :{Im} ")

    TFIM = OneDTFIM(Real,Im)
    TFIM.SetDetlaTau(0.1)
    TFIM.get_all_variables()

    iteration = 10000
    bin = 100

    data_array = np.zeros((len(h_array),int(iteration/bin)))
    
    # Simulation 
    for index, h in enumerate(h_array):
        TFIM.Seth_MagneticFieldStrength(h)
        #1. Thermalization 
        for i in range(10000):
            TFIM.MetropolisUpdate()
        print(f"Completed Thermalization at h = {TFIM.Geth_MagneticFieldStrength()}")

        #2. Collect data
        for i in range(int(iteration/bin)):
            EnergyPerSite_list = []
            for j in range(bin):
                TFIM.MetropolisUpdate()
                EnergyPerSite_list.append(TFIM.GetEnergyPerSite())
            
            #3. Bin and take sd & mean 
            data_array[index][i] = np.mean(EnergyPerSite_list)
        
        print(f" h = {TFIM.Geth_MagneticFieldStrength()} is completed.")
        print(f" Mean is {np.mean(data_array[index])}")
        print(f" Sd is {np.std(data_array[index])}")
    
    # Plot graph
    plt.figure(figsize=(10,10))
    plt.errorbar(h_array, 
                 np.mean(data_array, axis=1),
                 yerr=np.std(data_array, axis=1),
                 fmt='o', 
                 ecolor='black')
    
    plt.xlabel("h",fontsize=15)
    plt.ylabel("Energy Per Site",fontsize=15)
   
    # Saving the plot and data 
    current_dir = os.path.dirname(os.path.abspath(__file__))
    target_dir = os.path.join(current_dir, '..', 'Data', 'Checking h transition','Energy Per Site')
    plt.savefig(os.path.join(target_dir,f"L{Real}_ImL{Im}_tau0.1_Energy_{NoTrial}Trial.png"))
    np.save(os.path.join(target_dir, f'L{Real}_ImL{Im}_tau0.1_Energy_{NoTrial}Trial.npy'),data_array)

    #plt.show()

def test_msquare_and_mpowerfour():

    Real, Im = 4,5
    a = OneDTFIM(Real,Im)
    a.get_variable("TFIM")
    print("|m| is ",a.GetAbsolutemagnetization())
    print("|m|^2 is ",a.GetAbsolutemagnetizationSquare())
    
def CheckPhaseTransitionPoint_BinderCumulant(
        Real=4,
        Im=16,
        h_array=np.arange(0.7,1.4,0.1),
        NoTrial="new"):
    """
    This function will plot Binder Cumulant versus h graph, 
    with an aim to find the phase transition point. The data and graph 
    will be saved to the Data folder 
    """
    import os 
    
    print("This is Energy Versus h simulation of 1D transverse field Ising model")
    print(f"Real :{Real}, Im :{Im} ")

    TFIM = OneDTFIM(Real,Im)
    TFIM.SetDetlaTau(0.1)
    TFIM.get_all_variables()

    iteration = 10000
    bin = 100

    data_array_m2 = np.zeros((len(h_array),int(iteration/bin)))
    data_array_m4 = np.zeros((len(h_array),int(iteration/bin)))
    
    # Simulation 
    for index, h in enumerate(h_array):

        TFIM.Seth_MagneticFieldStrength(h)

        #1. Thermalization 
        for i in range(10000):
            TFIM.MetropolisUpdate()
        
        print(f"Completed Thermalization at h = {TFIM.Geth_MagneticFieldStrength()}")

        #2. Collect data
        for i in range(int(iteration/bin)):
            m_square_list = []
            m_power_4_list = []

            for j in range(bin):
                TFIM.MetropolisUpdate()
                m_square_list.append(TFIM.Get_magentizationSquare())
                m_power_4_list.append(TFIM.Get_magentizationPowerFour())
            
            #3. Bin and take sd & mean 
            data_array_m2[index][i] = np.mean(m_square_list)
            data_array_m4[index][i] = np.mean(m_power_4_list)
        
        print(f" h = {TFIM.Geth_MagneticFieldStrength()} is completed.")
        print(f" Mean of m^2 is {np.mean(data_array_m2[index])}")
        print(f" Mean of m^4 is {np.mean(data_array_m4[index])}")
    
    BinderCumulant = 1 - 1/3*(np.mean(data_array_m4, axis=1)/np.square(np.mean(data_array_m2, axis=1)))
    
    # Plot graph
    plt.figure(figsize=(10,10))
    plt.errorbar(h_array, 
                 BinderCumulant,
                 yerr=np.std(data_array_m4, axis=1),   # change the standard deviation later 
                 fmt='o', 
                 ecolor='black')
    
    plt.xlabel("h",fontsize=15)
    plt.ylabel("BinderRatio",fontsize=15)
   
    # Saving the plot and data 
    current_dir = os.path.dirname(os.path.abspath(__file__))
    target_dir = os.path.join(current_dir, '..', 'Data', 'Checking h transition','BinderCumulant and BinderRatio')
    plt.savefig(os.path.join(target_dir,f"L{Real}_ImL{Im}_tau0.1_BCversush_{NoTrial}Trial.png"))
    np.save(os.path.join(target_dir, f'L{Real}_ImL{Im}_tau0.1_m2_{NoTrial}Trial.npy'),data_array_m2)
    np.save(os.path.join(target_dir, f'L{Real}_ImL{Im}_tau0.1_m4_{NoTrial}Trial.npy'),data_array_m4)

def CheckingThermalizationValue(
        Real = 12,
        Im = 24,
        iteration = 10000,
        h = 1.2
):
    """
    This function is used to group the code that checks thermalization value.

    @param Real: the size/length of the 1D tranverse field Ising Model 
    @param Im: the imaginary time axis length of 1D tranverse field Ising Model  
    @param iteration: the iteration you want to try out 
    @param h: the magentic field strength
    
    """

    import matplotlib.pyplot as plt 
    import numpy as np

    TFIM = OneDTFIM(Real,Im)
    TFIM.SetDetlaTau(0.1)
    TFIM.Seth_MagneticFieldStrength(h)
    
    m_array = np.zeros(iteration)
    for i in range(iteration):
        TFIM.WolffUpdateBFSList()
        m_array[i] = TFIM.GetAbsolutemagentization()
        
    plt.plot(np.arange(0,iteration),m_array)
    plt.show()

def CheckPhaseTransitionPoint_susceptibility(
        Real=12,
        Im=24,
        thermalization = 10000,
        iteration = 10000,
        bin = 100,
        h_array=np.arange(0.5,1.6,0.1)
    ):
    """
    This function will plot sigma_mz versus h graph, with an aim to
    find the phase transition point. The data and graph will be saved to 
    Data folder. 

    @param Real: the size/length of the 1D tranverse field Ising Model 
    @param Im: the imaginary time axis length of 1D tranverse field Ising Model   
    @param thermalization: no of times of doing thermalization 
    @param iteration: the total number of times for collecting data
    @param bin: bin after the input of time 
    @param h_array: the array that specifies the range to check 

    """

    import os 
    
    print(f"Real :{Real}, Im :{Im} simulation")

    TFIM = OneDTFIM(Real,Im)
    TFIM.SetDetlaTau(0.1)
    bin = 100

    data_array = np.zeros((len(h_array),int(iteration/bin)))
    
    # Simulation 
    for index, h in enumerate(h_array):

        TFIM.Seth_MagneticFieldStrength(h)

        #1. Thermalization 
        for i in range(thermalization):
            TFIM.MetropolisUpdate()
        
        print(f"Completed Thermalization at h = {TFIM.Geth_MagneticFieldStrength()}")

        #2. Collect data
        for i in range(int(iteration/bin)):
            absolute_m_list = []
            absolute_msquare_list = []

            for j in range(bin):
                TFIM.MetropolisUpdate()
                absolute_m_list.append(TFIM.GetAbsolutemagnetization())
                absolute_msquare_list.append(TFIM.GetAbsolutemagnetizationSquare())
            
            #3. Bin and take sd & mean 
            data_array[index][i] = np.mean(absolute_msquare_list) - np.square(np.mean(absolute_m_list))
        
        print(f" h = {TFIM.Geth_MagneticFieldStrength()} is completed.")
        print(f" Mean is {np.mean(data_array[index])}")
        print(f" Sd is {np.std(data_array[index])}")
    
    # Plot graph
    plt.figure(figsize=(10,10))
    plt.errorbar(h_array, 
                 np.mean(data_array, axis=1),
                 yerr=np.std(data_array, axis=1),
                 fmt='o', 
                 ecolor='black')
    
    plt.xlabel("h")
    plt.ylabel("$\sigma^{2}_{m_z}$")

    # Save the graph 
    current_dir = os.path.dirname(os.path.abspath(__file__))
    target_dir = os.path.join(current_dir, '..', 'Data', 'Checking h transition','Absolute m')
    plt.savefig(os.path.join(target_dir,f"L{Real}_ImL{Im}_tau0.1_Susceptibility_10thTrial.png"))
    np.save(os.path.join(target_dir, f'L{Real}_ImL{Im}_tau0.1_Susceptibility_10thTrial.npy'),data_array)

def main():

    CheckPhaseTransitionPoint_abs_m(Real=12,
                                    Im=72,
                                    UpdateMethod="WolffUpdateBFSListTest",
                                    NoTrial="compare")
    
    return 0 


if __name__ == '__main__':
    import time
    start_time = time.time()
    main() 
    end_time = time.time()
    print(f"The code took {(end_time - start_time)/60} mins to run.")