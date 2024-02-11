import numpy as np
import os

def ImportVaryingImLData():

    os.chdir("Data\VaryingImL")

    data_path = "Data\VaryingImL"

    IndependentToGluedImL40 = np.load("Size4_ImL40_Iteration1000_IG.npy")
    GluedToIndependentImL40 = np.load("Size4_ImL40_Iteration1000_GI.npy")

    IndependentToGluedImL80 = np.load("Size4_ImL80_Iteration1000_IG.npy")
    GluedToIndependentImL80 = np.load("Size4_ImL80_Iteration1000_GI.npy")

    IndependentToGluedImL100 = np.load("Size4_ImL100_Iteration1000_IG.npy")
    GluedToIndependentImL100 = np.load("Size4_ImL100_Iteration1000_GI.npy")

    """
    IndependentToGluedImL120 = np.load(os.path.join(data_path, "Size4_ImL120_Iteration1000_IG.npy"))
    GluedToIndependentImL120 = np.load(os.path.join(data_path, "Size4_ImL120_Iteration1000_GI.npy"))
    
    """

    print("Renyi Entropy for SystemSize = 4, ")

    print("ImL = 40 :")
    print(IndependentToGluedImL40)
    print(GluedToIndependentImL40)
    print("Mean of IndependentToGlued is: ",np.mean(IndependentToGluedImL40))
    print("Mean of GluedToIndependent is: ",np.mean(GluedToIndependentImL40))
    print("The 2nd order Renyi entropy of ImL 40 is :", -np.log(np.mean(IndependentToGluedImL40)/
                                                      np.mean(GluedToIndependentImL40)))
    print()
    
    print("ImL = 80 :")
    print(IndependentToGluedImL80)
    print(GluedToIndependentImL80)
    print("Mean of IndependentToGlued is: ",np.mean(IndependentToGluedImL80))
    print("Mean of GluedToIndependent is: ",np.mean(GluedToIndependentImL80))
    print("The 2nd order Renyi entropy of ImL 80 is :", -np.log(np.mean(IndependentToGluedImL80)/
                                                      np.mean(GluedToIndependentImL80)))
    print()
    
    print("ImL = 100 :")
    print(IndependentToGluedImL100)
    print(GluedToIndependentImL100)
    print("Mean of IndependentToGlued is: ",np.mean(IndependentToGluedImL100))
    print(np.std(IndependentToGluedImL100))
    print("Mean of GluedToIndependent is: ",np.mean(GluedToIndependentImL100))
    print(np.std(GluedToIndependentImL100))

    print("The 2nd order Renyi entropy of ImL 100 is :", -np.log(np.mean(IndependentToGluedImL100)/
                                                      np.mean(GluedToIndependentImL100)))
    print()

    """
    print("ImL = 120 :")
    print(IndependentToGluedImL120)
    print(GluedToIndependentImL120)
    print("Mean of IndependentToGlued is: ",np.mean(IndependentToGluedImL120))
    print("Mean of GluedToIndependent is: ",np.mean(GluedToIndependentImL120))
    print("The 2nd order Renyi entropy is :", -np.log(np.mean(IndependentToGluedImL120)/
                                                      np.mean(GluedToIndependentImL120)))
    """
    

ImportVaryingImLData()
