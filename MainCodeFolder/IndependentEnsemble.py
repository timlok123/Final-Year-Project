from OneDTransverseFieldModel import OneDTFIM

class IndependentEnsemble: 
    """
    This class is used to initialize and manage IndependentEnsemble
    (2 disconnected 1D transverse field Ising model)
    """

    def __init__(self, SystemSize, Temperature, SubSystemSizeProp):
        """
        This function is to initialize the IndependentEnsemble. 

        @param SystemSize: the size/length of the 1D tranverse field Ising Model 
        @param Temperature: the temperature of the system 
        @param SubSystemSizeProp(A): the proportion of sub-system sizes in the system 

        The following variables are initialized, 
        - SystemSize
        - Temperature
        - ImaginaryTimeAxisLength
        - SubSystemSize
        - ThermalizationIteration 
        - SystemI
        - SystemII 
        
        """

        self.SystemSize = SystemSize
        self.Temperature = Temperature
        self.SubSystemSizeProp = SubSystemSizeProp
        self.ThermalizationIteration = 10000

        """
        1. Calculate the ImaginaryTimeAxisLength by using 
           - using the relation 1/(Temperature) = SystemSize*ImaginaryTimeAxisLength

           Note: We want set the ImaginaryTimeAxisLength to be 10 
           => set temperature to be 0.025 
        """
        #self.ImaginaryTimeAxisLength = int(1/(SystemSize*Temperature))
        self.ImaginaryTimeAxisLength = 40 # Fix the ImaginaryTimeAxisLength to be 40  first 
       
        """
        2. Calculate the SubSystemSize
        """
        self.SubSystemSize = int(SystemSize*SubSystemSizeProp)

        """
        3. Initialize the 2 independent 1D tranverse field Ising model 
        """
        self.SystemI = OneDTFIM(self.SystemSize,self.ImaginaryTimeAxisLength)
        self.SystemII = OneDTFIM(self.SystemSize,self.ImaginaryTimeAxisLength)
    
    def thermalization(self):
        """
        This function is to thermalize both subsystem 

        @param iteration: No of thermalization iteration  
        """
        self.SystemI.MetropolisUpdateNtimes(self.ThermalizationIteration)
        self.SystemII.MetropolisUpdateNtimes(self.ThermalizationIteration)
    
    def CheckingBoundaryCondition(self):
        """
        This function would check the boundary condition of the independent ensemble 
        If it acts like the glued ensemble, the function will return True 

        @return flag(boolean): true if it acts like glued ensemble
        
        """
        #(a==b).all()

        LayerIm_LMinus1SubI = self.SystemI.TFIM[self.ImaginaryTimeAxisLength-1][0:self.SubSystemSize] 
        LayerIm_LMinus1SubII = self.SystemII.TFIM[self.ImaginaryTimeAxisLength-1][0:self.SubSystemSize] 
        LayerZeroSubI = self.SystemI.TFIM[0][0:self.SubSystemSize]
        LayerZeroSubII= self.SystemII.TFIM[0][0:self.SubSystemSize]

        LayerIm_LMinus1Same = (LayerIm_LMinus1SubI == LayerIm_LMinus1SubII).all()
        LayerZeroSame = (LayerZeroSubI == LayerZeroSubII).all()

        return (LayerIm_LMinus1Same and LayerZeroSame)

    def CountIndependentToGlued(self,total_iteration=10000,iteration=160):
        """
        This function is to return the count of indepenent ensemble acts like 
        glued ensemble when (i%iteration steps) == 0 
        """

        """
        1. Thermalization 
           The number of iteration is specified by the instance variable ThermalizationIteration
        """
        self.thermalization() 

        """
        2. Keep on updating the subsystemI and subsystemII with methods thermalization().    
           Checking whether the independent ensemble acts like the glued ensemble 
           every iterations times. 
        """
        count = 0 

        for i in range(0,total_iteration):
            self.SystemI.MetropolisUpdateNtimes(iteration)
            self.SystemII.MetropolisUpdateNtimes(iteration)
            flag =  self.CheckingBoundaryCondition()
            if flag: 
                count+=1 
            print(r"Iteration {a} completed.".format(a = i*iteration))
    
        return count
