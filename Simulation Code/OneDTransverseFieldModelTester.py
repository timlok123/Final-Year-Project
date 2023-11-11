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
    print(randomIndex)
    s = Pocket[randomIndex]
    print(Pocket)
    print(s)

def testWolffUpdate():
    Real, Im = 4,4
    a = OneDTFIM(Real,Im)
    a.InitializeNeigbour()

    print("Configuration before update: ")
    print(a)
    
    a.WolffUpdateBFSList()

    print("Configuration after update: ")
    print(a)

def main():
    
    testWolffUpdate()


    return 0 


if __name__ == '__main__':
    main()


