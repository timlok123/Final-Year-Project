from OneDEnsemble import *

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
    Real, Im, Sub= 4,10,2
    a = OneDEnsemble(Real,Im,Sub)
    flag = False

    
    while(not flag):
        a = OneDEnsemble(Real,Im,Sub)
        print(a)
        flag = a.CheckingBoundaryCondition()
    
    #print(a)

    return 0 

def test_CountingBCTrueWolff():
    Real, Im, Sub= 32,40,16
    Iteration = 200000
    GluedSystem = OneDEnsemble(Real,Im,Sub)
    GluedSystem.__InitializeNeigbour__()
    GluedCount = GluedSystem.CountingBCTrueWolff(Iteration)

    print(GluedCount)

def test_CountingBCTrueWolffIndependent():
    Real, Im, Sub= 32,40,0
    Iteration = 200000
    IndependentSystem = OneDEnsemble(Real,Im,Sub)
    IndependentSystem.__InitializeNeigbour__()
    IndependentSystem.setCheckingRange(int(Real/2))
    IndependentCount = IndependentSystem.CountingBCTrueWolff(Iteration)

    print(IndependentCount)


def main():

    test_CountingBCTrueWolff()
    test_CountingBCTrueWolffIndependent()
    #test_inheritance()

    return 0 


if __name__ == '__main__':
    main()