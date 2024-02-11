import numpy as np
import matplotlib.pyplot as plt 
import os

#1. Import the data & calculate the value 
os.chdir("Data\VaryingSubSystemSize")

EntropyListL4 = []
for i in range(1,3):
    IndependentToGlued = np.load(f"Size4_SubSize{i}_Iteration1000_IG.npy")
    GluedToIndependent = np.load(f"Size4_SubSize{i}_Iteration1000_GI.npy") 
    EntropyListL4.append(-np.log(np.mean(IndependentToGlued)/np.mean(GluedToIndependent)))

EntropyListL6 = []
for i in range(2,5):
    IndependentToGlued = np.load(f"Size6_SubSize{i}_Iteration1000_IG.npy")
    GluedToIndependent = np.load(f"Size6_SubSize{i}_Iteration1000_GI.npy") 
    EntropyListL6.append(-np.log(np.mean(IndependentToGlued)/np.mean(GluedToIndependent)))

#2. Plot the graph 
plt.figure(figsize=(10,10))
plt.plot([i for i in range(1,4)],EntropyListL4,label="Size = 4")
plt.plot([i for i in range(2,5)],EntropyListL6,label="Size = 6")


plt.xlabel("Subsystem Size",size=15)
plt.ylabel("Renyi Entropy",size=15)
plt.legend()
plt.show()