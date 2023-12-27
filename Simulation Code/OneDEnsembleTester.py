from OneDEnsemble import *
import numpy as np 

def test_inheritance():
    """
    To see whether all variables can be successfully inherited 

    """
    Real, Im, Sub= 32,40,16
    IndependentSystem = OneDEnsemble(Real,Im,Sub)
    IndependentSystem.__InitializeNeigbour__()
    IndependentSystem.setCheckingRange(int(Real/2))
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

def MetroVaryConnected(Real,Connected):
    import os

    Im = int(Real/0.1)*2
    Iteration = 1000
    GluedSystem = OneDEnsemble(Real,Im,Connected)
    GluedSystem.__InitializeNeigbour__()
    GluedSystem.SetDetlaTau(0.1)

    print("Running the simulation of glued ensebmle, ")
    GluedToIndependentList = GluedSystem.CountingBCTrueMetro(Iteration)

    print("Running the simulation of Independent ensemble, ")
    Iteration = 1000
    IndependentSystem = OneDEnsemble(Real,Im,0)
    IndependentSystem.__InitializeNeigbour__()
    IndependentSystem.SetDetlaTau(0.1)
    IndependentSystem.SetCheckingRange(Connected)

    IndependentToGluedList = IndependentSystem.CountingBCTrueMetro(Iteration)

    print(IndependentToGluedList)
    print(GluedToIndependentList)

    print("Mean of Independent To Glued is: ", np.mean(IndependentToGluedList))
    print("Mean of Glued To Independent is: ", np.mean(GluedToIndependentList))
    print("The 2nd order Renyi entropy is :", -np.log(np.mean(IndependentToGluedList)/
                                                       np.mean(GluedToIndependentList)))
    
    np.save(os.path.join('Data','VaryingSubSystemSize', 
                         f"Size{Real}_SubSize{Connected}_Iteration{Iteration}_GI.npy"),
                         GluedToIndependentList)
    np.save(os.path.join('Data','VaryingSubSystemSize', 
                         f"Size{Real}_SubSize{Connected}_Iteration{Iteration}_IG.npy"),
                         IndependentToGluedList)

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

    #Add np.save(name of file, numpy array)
    np.save(os.path.join('Data','VaryingImL', 
                         f"Size{Real}_ImL{Im}_Iteration{Iteration}_GI.npy"),
                         GluedToIndependentList)
    np.save(os.path.join('Data','VaryingImL',
                          f"Size{Real}_ImL{Im}_Iteration{Iteration}_IG.npy"),
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
def savingNumpyArray():
    import os

    a = np.array([1,2,3]) 
    name1 = "TestArray1"
    np.save(os.path.join('..', 'Data\VaryingImL', f'{name1}.npy'),a)


def main():

    #MetroVaryImL(4,80)
    #MetroVaryConnected(4,1)
    #MetroVaryConnected(4,2)
    #MetroVaryConnected(4,3)
    #MetroVaryConnected(6,2)
    #MetroVaryConnected(6,3)
    MetroVaryConnected(6,5)

    return 0 


if __name__ == '__main__':
    main()