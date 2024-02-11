from OneDEnsemble import *
import numpy as np 

def test_inheritance():
    """
    To see whether all variables can be successfully inherited 

    """
    Real, Im, Sub= 4,4,2
    IndependentSystem = OneDEnsemble(Real,Im,Sub)
    IndependentSystem.__InitializeNeigbour__()
    IndependentSystem.SetCheckingRange(int(Real/2))
    IndependentSystem.get_all_variables()

def test_set_neighbour():
    Real, Im, Sub= 4,8,2
    a = OneDEnsemble(Real,Im,Sub)
    
    # Old way to initialize the neighbour 
    a.old_InitializeNeigbour()
    a.print_boundary_neigbour(4)
    a.print_boundary_neigbour(7)
    a.print_boundary_neigbour(3)
    a.print_boundary_neigbour(0)

    print("New way to initialize the neighbour, ")
    a.__InitializeNeigbour__()
    a.print_boundary_neigbour(4)
    a.print_boundary_neigbour(7)
    a.print_boundary_neigbour(3)
    a.print_boundary_neigbour(0)

def test_CheckingBoundaryCondition():
    Real, Im, Sub= 4,4,2
    a = OneDEnsemble(Real,Im,Sub)
    flag = False

    while(not flag):
        a = OneDEnsemble(Real,Im,Sub)
        print(a)
        flag = a.CheckingBoundaryCondition()
    
    #print(a)

    return 0 

def test_CountingBCTrueWolff():
    Real, Im, Connected= 10,100,5
    Iteration = 100000
    binIteration = 1000
    GluedSystem = OneDEnsemble(Real,Im,Connected)
    GluedSystem.__InitializeNeigbour__()
    GluedSystem.SetDetlaTau(0.1)

    GluedToIndependentList = GluedSystem.CountingBCTrueWolff(Iteration,binIteration)

    print("GluedToIndependentList is : ",GluedToIndependentList)
    print("The mean is :",np.mean(GluedToIndependentList))
    print("The errorbar is :",np.std(GluedToIndependentList)/np.sqrt(len(GluedToIndependentList)))

    return GluedToIndependentList 

def test_CountingBCTrueWolffIndependent():
    Real, Im, Connected= 10,100,0 
    
    Iteration = 100000
    binIteration = 1000
    
    IndependentSystem = OneDEnsemble(Real,Im,Connected)
    IndependentSystem.__InitializeNeigbour__() 
    IndependentSystem.SetDetlaTau(0.1)
  
    IndependentSystem.SetCheckingRange(5)

    IndependentToGluedList = IndependentSystem.CountingBCTrueWolff(Iteration,binIteration)

    print("IndependentToGluedList is : ",IndependentToGluedList )
    print("The mean is :",np.mean(IndependentToGluedList))
    print("The errorbar is :",np.std(IndependentToGluedList)/np.sqrt(len(IndependentToGluedList)))
    
    return IndependentToGluedList

def WolffTest(): 
    
    FinalResultList = [] 

    GGArray = test_CountingBCTrueWolff()
    IGArray  = test_CountingBCTrueWolffIndependent()
    RatioArray = IGArray/GGArray
    print("The ratio array is: ",RatioArray)
    print("The mean of the ratio array is: ",np.mean(RatioArray))
    FinalResultList.append((-1)*np.log(np.mean(RatioArray)))

    print("The FinalResultList is :",FinalResultList)
    print("The means of FinalResultList is :",np.mean(FinalResultList))

def test_MetropolisUpdate():
    Real, Im, Connected= 4,6,2
    GluedSystem = OneDEnsemble(Real,Im,Connected)
    GluedSystem.__InitializeNeigbour__()
    GluedSystem.SetDetlaTau(0.1)
    GluedSystem.MetropolisUpdate()

    Real, Im, Connected= 4,6,0
    IndependentSystem = OneDEnsemble(Real,Im,Connected)
    IndependentSystem.__InitializeNeigbour__()
    IndependentSystem.SetDetlaTau(0.1)
    IndependentSystem.SetCheckingRange(2)
    IndependentSystem.MetropolisUpdate()

def MetroTest():

    Real, Im, Connected= 6,60,0
    Iteration = 100
    IndependentSystem = OneDEnsemble(Real,Im,Connected)
    IndependentSystem.__InitializeNeigbour__()
    IndependentSystem.SetDetlaTau(0.1)
    IndependentSystem.SetCheckingRange(int(Real/2))

    IndependentToGluedList = IndependentSystem.CountingBCTrueMetro(Iteration)

    Real, Im, Connected= 6,60,3
    Iteration = 100
    GluedSystem = OneDEnsemble(Real,Im,Connected)
    GluedSystem.__InitializeNeigbour__()
    GluedSystem.SetDetlaTau(0.1)

    GluedToIndependentList = GluedSystem.CountingBCTrueMetro(Iteration)

    print(GluedToIndependentList)
    print(IndependentToGluedList)
    print(IndependentToGluedList/GluedToIndependentList)

    print(np.mean(IndependentToGluedList))
    print(np.mean(GluedToIndependentList))
    print(np.mean(IndependentToGluedList/GluedToIndependentList)) 

    #Add np.save(name of file, numpy array)
    np.save(f"Size{Real}SubSize{Connected}Iteration{Iteration}GITrial",GluedToIndependentList)
    np.save(f"Size{Real}SubSize{Connected}Iteration{Iteration}IGTrial",IndependentToGluedList)

def MetroVaryConnected(
        Real=4,
        Connected=2,
        Im=40,
        NoTrial="new",
        iteration=100,
        binIteration=1000):
    
    import os
    current_dir = os.path.dirname(os.path.abspath(__file__))
    target_dir = os.path.join(current_dir, '..', 'Data', 'VaryingSubSystemSize')
    
    GluedSystem = OneDEnsemble(Real,Im,Connected)
    GluedSystem.__InitializeNeigbour__()
    GluedSystem.SetDetlaTau(0.1)
    GluedSystem.Seth_MagneticFieldStrength(1)

    print("Running the simulation of glued ensebmle, ")
    GluedToIndependentList = GluedSystem.CountingBCTrueMetro(binIteration=binIteration,TotalNoOfUpdate=iteration)
    np.save(os.path.join(target_dir, 
                         f"Size{Real}_SubSize{Connected}_ImL{Im}_GI_{NoTrial}Trial.npy"),
                       GluedToIndependentList)
    

    print("Running the simulation of Independent ensemble, ")
    IndependentSystem = OneDEnsemble(Real,Im,0)
    IndependentSystem.__InitializeNeigbour__()
    IndependentSystem.SetDetlaTau(0.1)
    IndependentSystem.SetCheckingRange(Connected)
    IndependentSystem.Seth_MagneticFieldStrength(1)

    IndependentToGluedList = IndependentSystem.CountingBCTrueMetro(binIteration=binIteration,TotalNoOfUpdate=iteration)
    np.save(os.path.join(target_dir,
                         f"Size{Real}_SubSize{Connected}_ImL{Im}_IG_{NoTrial}Trial.npy"),
                         IndependentToGluedList)

    print("Mean of Independent To Glued is: ", np.mean(IndependentToGluedList))
    print("Mean of Glued To Independent is: ", np.mean(GluedToIndependentList))
    print("The 2nd order Renyi entropy is :", -np.log(np.mean(IndependentToGluedList)/
                                                       np.mean(GluedToIndependentList)))
    print("The errorbar value is :", np.sqrt(np.abs(np.std(GluedToIndependentList)/np.mean(GluedToIndependentList) 
                                                   - np.std(IndependentToGluedList)/np.mean(IndependentToGluedList)))/np.sqrt(iteration))
    

def MetroVaryImL(Real,Im):
    import os

    # Set Connected as int(Real/2) 
    Connected = int(Real/2) 

    Iteration = 1000
    GluedSystem = OneDEnsemble(Real,Im,Connected)
    GluedSystem.__InitializeNeigbour__()

    #Delta_tau is fixed to be 0.1
    GluedSystem.SetDetlaTau(0.1)

    GluedToIndependentList = GluedSystem.CountingBCTrueMetro(Iteration)

    Iteration = 1000
    IndependentSystem = OneDEnsemble(Real,Im,0)
    IndependentSystem.__InitializeNeigbour__()
    IndependentSystem.SetDetlaTau(0.1)
    IndependentSystem.SetCheckingRange(Connected)

    IndependentToGluedList = IndependentSystem.CountingBCTrueMetro(Iteration)

    print(IndependentToGluedList)
    print(GluedToIndependentList)

    print("Mean of Independent To Glued is: ",np.mean(IndependentToGluedList))
    print("Mean of Glued To Independent is: ",np.mean(GluedToIndependentList))
    print("The 2nd order Renyi entropy is :", -np.log(np.mean(IndependentToGluedList)/
                                                      np.mean(GluedToIndependentList)))

    current_dir = os.path.dirname(os.path.abspath(__file__))
    target_dir = os.path.join(current_dir, '..', 'Data', 'VaryingImL')
    np.save(os.path.join(target_dir,
                         f"Size{Real}_ImL{Im}_Iteration{Iteration}_GI_Jacktrial.npy"),
                         GluedToIndependentList)
    np.save(os.path.join(target_dir,
                          f"Size{Real}_ImL{Im}_Iteration{Iteration}_IG_Jacktrial.npy"),
                          IndependentToGluedList)

def loadingNumpyArray():
    # Example
    GIArraySub1 = np.load("Data\Size4\Size4SubSize2Iteration100GITrial.npy")
    IGArraySub1 = np.load("Data\Size4\Size4SubSize2Iteration100IGTrial.npy")

    print(GIArraySub1)
    print("Mean of GI is: ",np.mean(GIArraySub1))
    print(IGArraySub1)
    print("Mean of IG is: ",np.mean(IGArraySub1))

    print("The 2nd order Renyi entropy is :", -np.log(np.mean(IGArraySub1)/np.mean(GIArraySub1)))
    import os

    a = np.array([1,2,3]) 
    name1 = "TestArray1"
    np.save(os.path.join('..', 'Data\VaryingImL', f'{name1}.npy'),a)

def savingNumpyArray():

    import os 
    test = np.array([1,2,3,4])

    current_dir = os.path.dirname(os.path.abspath(__file__))

    target_dir = os.path.join(current_dir, '..', 'Data', 'VaryingSubSystemSize')
    np.save(os.path.join(target_dir, 
                         "test.npy"),
                         test) 

def RenyiEntropyCalculation(
        Real=4,
        Connected=2,
        Im=40,
        NoTrial="test",
        iteration=100,
        binIteration=10000,
        UpdatingMethod="MetropolisUpdate"):
    
    """
    This function is basically modified from MetroVaryConnected() function 
    and it can choose what CountingMethod do I use to calculate the Renyi Entropy. 

    @param Real           : The real axis length     
    @param Connected      : Connected site of the system. Smaller than Real. 
    @param Im             : The imaginary time axis length 
    @param NoTrial        : The index put on the final numpy file
    @param iteration      : The number of (binIteration) you would do 
    @param binIteration   : The number of iterations of update you would do 
    @param UpdatingMethod : The updating method you used to update the system 

    """
    
    import os
    current_dir = os.path.dirname(os.path.abspath(__file__))
    target_dir = os.path.join(current_dir, '..', 'Data', 'RenyiEntropyCalculation')
    
    GluedSystem = OneDEnsemble(Real,Im,Connected)
    GluedSystem.__InitializeNeigbour__()
    GluedSystem.SetDetlaTau(0.1)
    GluedSystem.Seth_MagneticFieldStrength(1)

    print("Running the simulation of glued ensebmle, ")

    GluedToIndependentList = GluedSystem.CountingBCTrue(binIteration=binIteration,
                                                        TotalNoOfUpdate=iteration,
                                                        UpdatingMethod=UpdatingMethod)

    np.save(os.path.join(target_dir, 
                         f"Size{Real}_SubSize{Connected}_ImL{Im}_GI_{UpdatingMethod}_{NoTrial}.npy"),
                       GluedToIndependentList)
    

    print("Running the simulation of Independent ensemble, ")
    IndependentSystem = OneDEnsemble(Real,Im,0)
    IndependentSystem.__InitializeNeigbour__()
    IndependentSystem.SetDetlaTau(0.1)
    IndependentSystem.SetCheckingRange(Connected)
    IndependentSystem.Seth_MagneticFieldStrength(1)

    IndependentToGluedList = IndependentSystem.CountingBCTrue(binIteration=binIteration,
                                                              TotalNoOfUpdate=iteration,
                                                              UpdatingMethod=UpdatingMethod)
    np.save(os.path.join(target_dir,
                         f"Size{Real}_SubSize{Connected}_ImL{Im}_IG_{UpdatingMethod}_{NoTrial}.npy"),
                         IndependentToGluedList)

    print("Mean of Independent To Glued is: ", np.mean(IndependentToGluedList))
    print("Mean of Glued To Independent is: ", np.mean(GluedToIndependentList))
    print("The 2nd order Renyi entropy is :", -np.log(np.mean(IndependentToGluedList)/
                                                       np.mean(GluedToIndependentList)))
    print("The errorbar value is :", np.sqrt(np.abs(np.std(GluedToIndependentList)/np.mean(GluedToIndependentList) 
                                                   - np.std(IndependentToGluedList)/np.mean(IndependentToGluedList)))/np.sqrt(iteration))

def RenyiEntropyCalculation_GI(
        Real=4,
        Connected=2,
        Im=40,
        NoTrial="test",
        iteration=100,
        binIteration=10000,
        UpdatingMethod="MetropolisUpdate"):
    
    """
    (Fill out the description later)

    @param Real           : The real axis length     
    @param Connected      : Connected site of the system. Smaller than Real. 
    @param Im             : The imaginary time axis length 
    @param NoTrial        : The index put on the final numpy file
    @param iteration      : The number of (binIteration) you would do 
    @param binIteration   : The number of iterations of update you would do 
    @param UpdatingMethod : The updating method you used to update the system 

    """
    
    import os
    current_dir = os.path.dirname(os.path.abspath(__file__))
    target_dir = os.path.join(current_dir, '..', 'Data', 'RenyiEntropyCalculation')

    print("Running the simulation of glued ensebmle, ")
    
    GluedSystem = OneDEnsemble(Real,Im,Connected)
    GluedSystem.__InitializeNeigbour__()
    GluedSystem.SetDetlaTau(0.1)
    GluedSystem.Seth_MagneticFieldStrength(1)

    GluedToIndependentList = GluedSystem.CountingBCTrue(binIteration=binIteration,
                                                        TotalNoOfUpdate=iteration,
                                                        UpdatingMethod=UpdatingMethod)

    np.save(os.path.join(target_dir, 
                         f"Size{Real}_SubSize{Connected}_ImL{Im}_GI_{UpdatingMethod}_{iteration}_{NoTrial}.npy"),
                       GluedToIndependentList)

def RenyiEntropyCalculation_IG(
        Real=4,
        Connected=2,
        Im=40,
        NoTrial="test",
        iteration=100,
        binIteration=10000,
        UpdatingMethod="MetropolisUpdate"):
    
    """
    (Fill out the description later)

    @param Real           : The real axis length     
    @param Connected      : Connected site of the system. Smaller than Real. 
    @param Im             : The imaginary time axis length 
    @param NoTrial        : The index put on the final numpy file
    @param iteration      : The number of (binIteration) you would do 
    @param binIteration   : The number of iterations of update you would do 
    @param UpdatingMethod : The updating method you used to update the system 

    """
    
    import os
    current_dir = os.path.dirname(os.path.abspath(__file__))
    target_dir = os.path.join(current_dir, '..', 'Data', 'RenyiEntropyCalculation')

    print("Running the simulation of Independent ensemble, ")
    
    IndependentSystem = OneDEnsemble(Real,Im,0)
    IndependentSystem.__InitializeNeigbour__()
    IndependentSystem.SetDetlaTau(0.1)
    IndependentSystem.SetCheckingRange(Connected)
    IndependentSystem.Seth_MagneticFieldStrength(1)

    IndependentToGluedList = IndependentSystem.CountingBCTrue(binIteration=binIteration,
                                                              TotalNoOfUpdate=iteration,
                                                              UpdatingMethod=UpdatingMethod)
    np.save(os.path.join(target_dir,
                         f"Size{Real}_SubSize{Connected}_ImL{Im}_IG_{UpdatingMethod}_{iteration}_{NoTrial}.npy"),
                         IndependentToGluedList)
    
def main():

    #RenyiEntropyCalculation(Real=4,Connected=2,Im=40)
    size = 16
    Im,Connected = size*10, 8
    flagGI,flagIG = True,False
    #flagGI,flagIG = False,True

    if flagGI: 
        RenyiEntropyCalculation_GI(Real=size,
                            Connected=Connected,
                            Im=Im, 
                            UpdatingMethod="WolffUpdateBFSList", 
                            iteration=50,
                            binIteration=10000)
                        
    if flagIG: 
        RenyiEntropyCalculation_IG(Real=size,
                            Connected=Connected,
                            Im=Im,
                            UpdatingMethod="WolffUpdateBFSList",
                            iteration=50,
                            binIteration=10000)
    
    return 0  

if __name__ == '__main__': 
    import time
    start_time = time.time()

    main() 

    end_time = time.time()
    print(f"The code took {(end_time - start_time)/60} min(s) to run.")

