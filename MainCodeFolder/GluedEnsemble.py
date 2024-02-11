import numpy as np 
from OneDTransverseFieldModel import OneDTFIM

class GluedEnsemble: 
    """
    This class is used to initialize and manage GluedEnsemble
    It will be implement as the 2D Ising Model with 
    dimension SystemSize*(2*ImaginaryTimeAxisLength)
    (2 connected 1D transverse field Ising model in A)
    """

    def __init__(self,SystemSize,Temperature,SubSystemSizeProp):
        """
        This function is used to initialize the GluedEnsemble. 

        @param SystemSize: the size/length of the 1D tranverse field Ising Model 
        @param Temperature: the temperature of the system 
        @param SubSystemSizeProp(A): the proportion of sub-system sizes in the system 

        The following variables are initialized, 
        - SystemSize
        - Temperature
        - ImaginaryTimeAxisLength
        - SubSystemSize
        - ThermalizationIteration (set 10000 times as default)
        - System (glued 2 1D Transverse field Ising Model together)
        
        """

        self.SystemSize = SystemSize
        self.Temperature = Temperature
        self.SubSystemSizeProp = SubSystemSizeProp
        self.ThermalizationIteration = 10000
        self.CountIteration = 10000

        """
        1. Calculate the ImaginaryTimeAxisLength by using 
           - using the relation 1/(Temperature) = SystemSize*ImaginaryTimeAxisLength

           Note: We want set the ImaginaryTimeAxisLength to be 10 
           => set temperature to be 0.025 
        """
        #self.ImaginaryTimeAxisLength = int(1/(SystemSize*Temperature))
        self.ImaginaryTimeAxisLength = 40 # fix it to be 40 first 
    
        """
        2. Calculate the SubSystemSize
        """
        self.SubSystemSize = int(SystemSize*SubSystemSizeProp)

        """
        3. Initialize the 1D tranverse field Ising model with
        dimension SystemSize*(2*ImaginaryTimeAxisLength) 
           
        """
        self.System = OneDTFIM(self.SystemSize,2*self.ImaginaryTimeAxisLength) 

        self.J_x = self.System.J_x 
        self.J_y = self.System.J_y
        self.beta_cl = self.System.beta_cl


    """
    def WolffUpdateGluedEnsemble(self):

        randomX = np.random.randint(0,self.SystemSize)
        randomY = np.random.randint(0,self.ImaginaryTimeAxisLength)

        Pocket = [(randomX,randomY)]

        ## Complete it later 
    """
 

    def MetropolisUpdateGluedEnsemble(self):
        """
        This function is to update the 1D transverse field Ising model(Glued ensemble)
        with Metropolis algorithm. It will change the self.TFIM (np.array) afterwards. 

        Steps 
        1. Randomly choose a Im-L layer and L number 
        2. check neighbour energy (same layer & different layer)
            if it is on the boundary and it is not in the glued part 
            => assign to correct location 
        3. decide whether you will update (deltaE and weight)
        """
        for i in range(self.System.TotalMCsteps): 

            #1. Randomly choose Im-Length layer and SystemSize number 
            SystemSizeIndexRandom = np.random.randint(0,self.SystemSize)
            ImTimeAxisIndexRandom = np.random.randint(0,self.ImaginaryTimeAxisLength)

            #2. Calculate the energy of its neighbour 
            SumOfSpinInSameLayer = 0 
            SumOfSpinInDiffLayer = 0 

            #print(SystemSizeIndexRandom)

            LeftIndex = (SystemSizeIndexRandom-1)%(self.SystemSize)
            RightIndex = (SystemSizeIndexRandom+1)%(self.SystemSize)
            DownIndex = (ImTimeAxisIndexRandom-1)%(self.ImaginaryTimeAxisLength)
            UpIndex = (ImTimeAxisIndexRandom+1)%(self.ImaginaryTimeAxisLength)

            if (ImTimeAxisIndexRandom==(self.ImaginaryTimeAxisLength-1)) and (SystemSizeIndexRandom>=self.SubSystemSize):
                UpIndex = 0

            if (ImTimeAxisIndexRandom==self.ImaginaryTimeAxisLength) and (SystemSizeIndexRandom>=self.SubSystemSize):
                DownIndex = 2*self.ImaginaryTimeAxisLength-1
            
            SumOfSpinInSameLayer = self.System.TFIM[ImTimeAxisIndexRandom][LeftIndex] + self.System.TFIM[ImTimeAxisIndexRandom][RightIndex]
            SumOfSpinInDiffLayer = self.System.TFIM[DownIndex][SystemSizeIndexRandom] + self.System.TFIM[UpIndex][SystemSizeIndexRandom]

            #3. Determine whehter to flip or not  
            DeltaETotal = 2*self.System.TFIM[ImTimeAxisIndexRandom][SystemSizeIndexRandom]* \
                        (self.J_x*SumOfSpinInSameLayer + self.J_y*SumOfSpinInDiffLayer)
            
            weight = np.exp(-DeltaETotal*self.beta_cl)

            #print("test3")

            if (DeltaETotal < 0) or (np.random.rand() < weight):
                self.System.TFIM[ImTimeAxisIndexRandom][SystemSizeIndexRandom] *= -1 
    
    def MetropolisUpdateGluedEnsembleNtimes(self,N):
        """
        This function is just to repeat the method MetropolisUpdateGE for N times 
        @param : no of times of calling MetropolisUpdateGE
        
        """
        for i in range(N):
            self.MetropolisUpdateGluedEnsemble()

    def thermalization(self):
        self.MetropolisUpdateGluedEnsembleNtimes(self.ThermalizationIteration)
    
    def CheckingBoundaryCondition(self):
        """
        This function would check the boundary condition of the glued ensemble 
        If it acts like the independent ensemble, the function will return True 

        @return flag(boolean): true if it acts like independent ensemble
        
        """
       
        LayerIm_LMinus1SubI = self.System.TFIM[self.ImaginaryTimeAxisLength-1][0:self.SubSystemSize] 
        LayerIm_LMinus1SubII = self.System.TFIM[2*self.ImaginaryTimeAxisLength-1][0:self.SubSystemSize] 
        LayerZeroSubI = self.System.TFIM[0][0:self.SubSystemSize]
        LayerZeroSubII= self.System.TFIM[self.ImaginaryTimeAxisLength][0:self.SubSystemSize]

        LayerIm_LMinus1Same = (LayerIm_LMinus1SubI == LayerIm_LMinus1SubII).all()
        LayerZeroSame = (LayerZeroSubI == LayerZeroSubII).all()

        return (LayerIm_LMinus1Same and LayerZeroSame)

    def CountGluedToIndependent(self):
        """
        This function is to return the count of glued ensemble acts like 
        independent ensemble 

        @return count(int) : the number of times that the system acts like
                             independent ensemble  
        """

        """
        1. Thermalization 
           The number of iteration is specified by the instance variable ThermalizationIteration
        """
    
        self.MetropolisUpdateGluedEnsembleNtimes(self.ThermalizationIteration)

        """
        2. Keep on updating the subsystemI and subsystemII. 
           Checking whether the independent ensemble acts like the glued ensemble. 
        """

        count = 0 
        for i in range(0,self.CountIteration):
            self.MetropolisUpdateGluedEnsemble()
            if (i%100==0 and i!=0) and self.CheckingBoundaryCondition():
                count+=1 
            print(r"Iteration {a} completed.".format(a = i*self.CountIteration))
    
        return count

    def AverageCountGluedToIndependent(self,RepeatTimes=10):
        count_array = np.zeros(RepeatTimes)
        for i in range(RepeatTimes):
            print("Trial {no} starts".format(no=i))
            count_array[i] = self.CountGluedToIndependent()
        
        return count_array

        

def TestGluedEnsembleMethod(): 
 
    #1. Test the initialization - See whether the System is initialized correctly 
    Systemsize = 4 
    T = 0.025
    A = 0.5 
    test = GluedEnsemble(Systemsize,T,A)
    # It correctly set up the model with correct dimension 

    #2. Test the MetropolisUpdateGE() method
    # a. Fix the ImTimeAxisIndexRandom = self.ImaginaryTimeAxisLength - 1  and test UpIndex
    # b. Fix the ImTimeAxisIndexRandom = self.ImaginaryTimeAxisLength 1  and test DownIndex
    #test.MetropolisUpdateGluedEnsembleNtimes(10)
    temp_array = test.AverageCountGluedToIndependent(15)
    print(temp_array)
    print(np.mean(temp_array))



    #3. Test the CheckingBoundaryCondition()
    #test.CheckingBoundaryCondition()

    #4. Test the CountGluedToIndependent()
    #iteration = 10000
    #count = test.CountGluedToIndependent(iteration)

    #print(count)


def main():

    TestGluedEnsembleMethod()

    return 0 


if __name__ == '__main__':
    main()
 