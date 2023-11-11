"""
OneDEnsemble.py 
This code file is used to model Glued ensemble and Independent ensemble 
Lastest Update on 10 Nov 2023 by Justin Chau 

"""

from OneDTransverseFieldModel import OneDTFIM

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
        super().__init__(RealAxisLength,ImaginaryTimeAxisLength)

        self.SubSystemSize = SubSystemSize 
        
        if (self.SubSystemSize!=0):
            self.CheckingRange = SubSystemSize
        else:
            print("Remember to set the CheckingRange by using setCheckingRange()")
            print("It is set to be half of RealAxisLength first")
            self.CheckingRange = int(self.SystemSize/2)
    
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

    def setCheckingRange(self,CR):
        """
        This is the setter of variables checking range 
        """
        self.CheckingRange = CR
         
    
    def CheckingBoundaryCondition(self):
        """
        Checking whether the index in checking range in layer 
        i. ImaginaryTimeAxisLength-1 and int(ImaginaryTimeAxisLength/2)-1
        ii. 0 and int(self.ImaginaryTimeAxisLength/2)

        are the same. They have to be both same to return true 
        
        """

        for ReL_index in range(0,self.CheckingRange):
            if ((self.TFIM[ReL_index][self.ImaginaryTimeAxisLength-1] != \
                self.TFIM[ReL_index][int(self.ImaginaryTimeAxisLength/2)-1]) or \
                (self.TFIM[ReL_index][0] != \
                 self.TFIM[ReL_index][int(self.ImaginaryTimeAxisLength/2)])):
                 return False
        
        return True 

    def CountingBCTrueWolff(self,TotalNoOfUpdate,bins=1):
        """
        This function will repeatedly call Wolff algo to update the configuration for 
        TotalNoOfUpdate, and check how many times how many time 
        CheckingBoundaryCondition() would return true. 

        @param TotalNoOfUpdate: Total no of updating the configurations by using Wolff 
        @param bins: Divide TotalNoOfUpdate to be no of bins. Take average every bins time. 
        
        """

        count = 0 
        for i in range(0,TotalNoOfUpdate):
            self.WolffUpdateBFSList()
            if (self.CheckingBoundaryCondition()):
                count += 1
            if (i%10000 == 0) and (i!=0):
                print(f"Iteration {i} completed.")
            
        return count 
        
    """
    Below are the function used for debugging 
    """
    
    def old_InitializeNeigbour(self):
        super().__InitializeNeigbour__()