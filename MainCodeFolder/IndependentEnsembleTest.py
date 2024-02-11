from IndependentEnsemble import *
import numpy as np

def TestBasic():
    a = 1 
    a+=1
    print(a)
    print(True and True)

def TestArrayOperation():
    test = np.random.choice([1,-1],size=(5,5))
    test1 = test[4][0:3]
    print(test)
    print(test1)
    print(type(test1))
    print(test1.shape)

    a = np.array([0,0,1,1,1])
    b = np.array([0,0,1,1,1])
    print((a==b).all())
    print(type((a==b).all()))

def TestMethod():

    TestIndependentEnsemble = IndependentEnsemble(4,0.025,0.5)
    #TestIndependentEnsemble.thermalization(10)
    #TestIndependentEnsemble.CheckingBoundaryCondition()
    iteration = 10000
    count_array = np.zeros(10)
    for i in range(10):
        count = TestIndependentEnsemble.CountIndependentToGlued()
        count_array[i] = count
    
    print(count_array)
    print(np.mean(count_array))

    return 0 
def main():
    
    TestMethod()
    return 0 

if __name__ == '__main__':
    main()