"""
OneDEnsemble.py 
This code file is used to model Glued ensemble and Independent ensemble 
Lastest Update on 10 Nov 2023 by Justin Chau 

"""

from OneDTransverseFieldModel import OneDTFIM
import numpy as np

class OneDEnsemble(OneDTFIM):

    """
    This class is used to initialize and manage the ensemble model. 
    """

    def __init__(self,RealAxisLength,ImaginaryTimeAxisLength,SubSystemSize):
        """
        This function is to initialize the 1D tranverse field Ising Model. 
        The 1D TFIM is stored as instance 

        @param RealAxisLength: the size/length of the 1D tranverse field Ising Model 
        @param ImaginaryTimeAxisLength: the imaginary time axis length of 1D tranverse 
                                         field Ising Model 
        @param SubSystemSize: the size of the subsystem 

        """

        self.SubSystemSize = SubSystemSize 

        super().__init__(RealAxisLength,ImaginaryTimeAxisLength)

        if (self.SubSystemSize!=0):
            self.CheckingRange = SubSystemSize
        elif (self.SubSystemSize==0):
            print("Remember to set the CheckingRange by using setCheckingRange()")
            print("It is set to be 0 first")
            self.CheckingRange = 0
        else:
            print("Invalid CheckingRange")
            self.CheckingRange = 0
            
    
    def __InitializeNeigbour__(self):
        """
        This function is set to find the coordinate of neighbour of every sites and stored 
        it in neighbour. 

        The methods in this class overrides the same methods in OneDTFIM. 

        Structure:
        [ReL][ImL][Orientation][xy]

        - Orientation
        {
            "Left":  0 
            "Right": 1 
            "Up":    2 
            "Down"   3 
        }

        - xy
        {
           "x":0 
           "y":1 

        }
        
        """

        super().__InitializeNeigbour__()
        
        #1. Change the boundary condition at ImL_index = ImaginaryTimeAxisLength/2 & ImaginaryTimeAxisLength - 1 
        for ReL in range(self.SubSystemSize, self.SystemSize):
            self.neighbour[ReL][int(self.ImaginaryTimeAxisLength/2)][3][1] = self.ImaginaryTimeAxisLength - 1 
            self.neighbour[ReL][self.ImaginaryTimeAxisLength-1][2][1] = int(self.ImaginaryTimeAxisLength/2)

        #2. Change the boundary condition at ImL_index = ImaginaryTimeAxisLength/2-1 & 0  
        for ReL in range(self.SubSystemSize, self.SystemSize):
            self.neighbour[ReL][int(self.ImaginaryTimeAxisLength/2)-1][2][1] = 0 
            self.neighbour[ReL][0][3][1] = int(self.ImaginaryTimeAxisLength/2) - 1 

    def SetCheckingRange(self,CR):
        """
        This is the setter of variables checking range 
        """
        self.CheckingRange = CR
        print("The Checking Range is now set to :",self.CheckingRange)
    
    def MetropolisUpdate(self):
        """
        This is an overriding methods. This function is to update the 1D transverse field Ising 
        model ensemble with Metropolis algorithm. It will change the self.TFIM (np.array) 
        afterwards.

        Instead of calculating the index directly, it will take the the neighbour from the 
        neighbour array. 

        """

        for i in range(self.TotalMCsteps):
            #1. Randomly choose Im-Length layer and SystemSize number 
            SystemSizeIndexRandom = np.random.randint(0,self.SystemSize)
            ImTimeAxisIndexRandom = np.random.randint(0,self.ImaginaryTimeAxisLength)
            SpinAtThisSite = self.TFIM[SystemSizeIndexRandom][ImTimeAxisIndexRandom]

            #2. Calculate the energy of its neighbour 
            SumOfSpinInSameLayer = 0 
            SumOfSpinInDiffLayer = 0 

            LeftReLCoor  = self.neighbour[SystemSizeIndexRandom][ImTimeAxisIndexRandom][0][0]
            RightReLCoor = self.neighbour[SystemSizeIndexRandom][ImTimeAxisIndexRandom][1][0]
            UpImLCoor    = self.neighbour[SystemSizeIndexRandom][ImTimeAxisIndexRandom][2][1]
            DownImLCoor  = self.neighbour[SystemSizeIndexRandom][ImTimeAxisIndexRandom][3][1]
            
            SumOfSpinInSameLayer = self.TFIM[LeftReLCoor][ImTimeAxisIndexRandom] + self.TFIM[RightReLCoor][ImTimeAxisIndexRandom]
            SumOfSpinInDiffLayer = self.TFIM[SystemSizeIndexRandom][UpImLCoor] + self.TFIM[SystemSizeIndexRandom][DownImLCoor]

            #3. Determine whether to flip or not  
            DeltaETotal = 2*SpinAtThisSite*(self.J_x*SumOfSpinInSameLayer + self.J_y*SumOfSpinInDiffLayer)
            weight = np.exp(-DeltaETotal*self.beta_cl)

            if (DeltaETotal <= 0) or (np.random.rand() < weight):
                self.TFIM[SystemSizeIndexRandom][ImTimeAxisIndexRandom] *= -1

    def CheckingBoundaryCondition(self):
        """
        Checking whether the index in checking range in layer 
        i. ImaginaryTimeAxisLength-1 and int(ImaginaryTimeAxisLength/2)-1
        ii. 0 and int(self.ImaginaryTimeAxisLength/2)

        are the same. They have to be both same to return true 
        
        """

        """
        for ReL_index in range(0,self.CheckingRange):
            if (self.TFIM[ReL_index][0] != \
                 self.TFIM[ReL_index][int(self.ImaginaryTimeAxisLength/2)]):
                return False
    
        return True """

        half_length = int(self.ImaginaryTimeAxisLength / 2)
        return np.all(self.TFIM[:self.CheckingRange, 0] == self.TFIM[:self.CheckingRange, half_length])
        
    
    def CountingBCTrueMetro(
            self,
            ThermalizationIteration=10000,
            TotalNoOfUpdate=10,
            binIteration=20000
            ):
        """
        This function will repeatedly call Metropolis algorithm to update the configuration for 
        TotalNoOfUpdate times, and check how many times CheckingBoundaryCondition() 
        would return true. 

        @param TotalNoOfUpdate: Total no of updating the configurations by using Wolff 
        
        """

        print("The RealAxisLength is: ",self.SystemSize)
        print("The CheckingRange is: " ,self.CheckingRange)
        print("The ImaginaryTimeAxisLength is: ", self.ImaginaryTimeAxisLength)

        count = 0 
        RatioList = [] 
   
        #1. Thermalization 
        for i in range(ThermalizationIteration):
            self.MetropolisUpdate()

        for iteration in range(TotalNoOfUpdate):

            #2. Taking data 
            for i in range(binIteration):
                self.MetropolisUpdate()
                if (self.CheckingBoundaryCondition()):
                    count += 1

            #3. Record result 
            RatioList.append(count/binIteration)
            count = 0 
            print(f"{iteration+1} iteration(s) completed.In 1 iteration, it will do {binIteration} MC steps update")
            
        return np.array(RatioList)

    def CountingBCTrue(
            self,
            ThermalizationIteration=10000,
            TotalNoOfUpdate=10,
            binIteration=20000,
            UpdatingMethod="MetropolisUpdate"):
        """
        This function will repeatedly call the selected algo to update the configuration for 
        TotalNoOfUpdate times, and check how many times CheckingBoundaryCondition() 
        would return true. 

        (Fill out the documentation later)
        
        """

        print("The RealAxisLength is: ",self.SystemSize)
        print("The CheckingRange is: " ,self.CheckingRange)
        print("The ImaginaryTimeAxisLength is: ", self.ImaginaryTimeAxisLength)
        print("The algorithm selected for updating the TFIM is: ",UpdatingMethod)

        count = 0 
        RatioArray = np.zeros(TotalNoOfUpdate)
   
        #1. Thermalization 
        for i in range(ThermalizationIteration):
            getattr(self, UpdatingMethod)()

        for iteration in range(TotalNoOfUpdate):

            #2. Taking data 
            for i in range(binIteration):
                getattr(self, UpdatingMethod)()
                if (self.CheckingBoundaryCondition()):
                    count += 1

            #3. Record result 
            RatioArray[iteration] = count/binIteration
            count = 0 
            print((f"{iteration+1} iteration(s) completed."
                f"In 1 iteration, it will do {binIteration} steps of updating"))
            
        return RatioArray 

    """
    Below are the function used for debugging 
    """
    
    def old_InitializeNeigbour(self):
        super().__InitializeNeigbour__()