"""
1DTransverseFieldModel.py 
This program defines 1DTransverseFieldModel and the updating methods 
Lastest Update on 10 Nov 2023 by Justin Chau 
"""

import numpy as np 
import matplotlib.pyplot as plt 

class OneDTFIM: 
    """
    This class is used to initialize and manage the 1D transverse field Ising Model.
    """

    def __init__(self,RealAxisLength,ImaginaryTimeAxisLength):
        """
        This function is to initialize the 1D tranverse field Ising Model. 
        The 1D TFIM is stored as instance 

        @param RealAxisLength: the size/length of the 1D tranverse field Ising Model 
        @param ImaginaryTimeAxisLength: the imaginary time axis length of 1D tranverse 
                                         field Ising Model 
        """
        
        self.SystemSize = RealAxisLength
        self.ImaginaryTimeAxisLength = ImaginaryTimeAxisLength
        self.TFIM = np.random.choice([1,-1],size=(RealAxisLength,
                                                      ImaginaryTimeAxisLength))
        
        self.J_x = 1 

        # Depends on delta_tau
        self.delta_tau = 0.1  
        # Set h to be 1 as default      
        self.h = 1   
        
        # Other value depends on delta_tau and h 
        self.beta_cl = self.delta_tau
        self.gamma = -0.5*np.log(np.tanh(self.delta_tau*self.h))
        self.J_y = (self.gamma)/(self.beta_cl)
  
        # 4 directions, 2 coordinates (x,y)
        self.neighbour = np.zeros((RealAxisLength,ImaginaryTimeAxisLength, 4, 2),dtype=int) 


        self.TotalMCsteps = RealAxisLength*ImaginaryTimeAxisLength 
        # 1 MC steps = 2D Ising model size

        self.__InitializeNeigbour__()

    def __InitializeNeigbour__(self):
        """
        This function is set to find the coordinate of neighbour of every sites and stored 
        it in neighbour. 

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

        for Re_L in range(0,self.SystemSize):
            for Im_L in range(0,self.ImaginaryTimeAxisLength):

                self.neighbour[Re_L][Im_L][0][0] = (Re_L-1)%self.SystemSize                    #left_x
                self.neighbour[Re_L][Im_L][0][1] = Im_L                                        #left_y
                self.neighbour[Re_L][Im_L][1][0] = (Re_L+1)%self.SystemSize                    #right_x
                self.neighbour[Re_L][Im_L][1][1] = Im_L                                        #right_y
                self.neighbour[Re_L][Im_L][2][0] = Re_L                                        #up_x   
                self.neighbour[Re_L][Im_L][2][1] = (Im_L+1)%self.ImaginaryTimeAxisLength       #up_y   
                self.neighbour[Re_L][Im_L][3][0] = Re_L                                        #down_x
                self.neighbour[Re_L][Im_L][3][1] = (Im_L-1)%self.ImaginaryTimeAxisLength       #down_y

    def SetDetlaTau(self,Dtau):
        """
        This function is to set the delta_tau value of the 1D transverse field Ising 
        model. 

        Since by our convention, we set delta_tau = beta_cl. The function would also change the 
        value of the beta_cl. Other values that depends on beta_cl an delta_tau, will be calculated 
        again. 

        @Dtau: New delta_tau value 
        """

        self.delta_tau = Dtau
        self.beta_cl = Dtau

        self.gamma = -0.5*np.log(np.tanh(Dtau*self.h))
        self.J_y = (self.gamma)/Dtau

        print("The delta tau is reset.")
        print("The delta tau is set to be :", self.delta_tau)
        print("The beta_cl is set to be :", self.beta_cl) 
        print("The gamma is set to be :", self.gamma) 
        print("J_y is set to be :", self.J_y) 
    
    def GetEnergyPerSite(self):
        """
        This function is to return the Energy per site (E/n) value of the 1D transverse 
        field Ising model.

        @return EPerN: Energy per site (E/n) value of the model 
        """
        
        #print("Remember to call __InitializeNeigbour__() before you use GetEnergyPerSite()")

        Energy = 0 

        for Re_L in range(0,self.SystemSize):
            for Im_L in range(0,self.ImaginaryTimeAxisLength):

                SpinAtThisSite = self.TFIM[Re_L][Im_L]

                LeftReLCoor  = self.neighbour[Re_L][Im_L][0][0]
                RightReLCoor = self.neighbour[Re_L][Im_L][1][0]
                UpImLCoor    = self.neighbour[Re_L][Im_L][2][1]
                DownImLCoor  = self.neighbour[Re_L][Im_L][3][1]
            
                SumOfSpinInSameLayer = self.TFIM[LeftReLCoor][Im_L] + \
                                       self.TFIM[RightReLCoor][Im_L]
                SumOfSpinInDiffLayer = self.TFIM[Re_L][UpImLCoor] + \
                                       self.TFIM[Re_L][DownImLCoor]
    
                Energy += -SpinAtThisSite*(self.J_x*SumOfSpinInSameLayer + 
                                           self.J_y*SumOfSpinInDiffLayer)

        return (Energy/(2*self.SystemSize*self.ImaginaryTimeAxisLength))

    def GetAbsolutemagnetization(self):
        """
        This function is to return the |m| value of the 1D transverse field Ising 
        model.

        @return m: Absolute magentization value of the model 
        """

        return np.abs(np.sum(self.TFIM)/(self.SystemSize*self.ImaginaryTimeAxisLength)) 
    
    def GetAbsolutemagnetizationSquare(self):
        """
        This function is to return the m^2 value of the 1D transverse field Ising 
        model.

        @return m: Absolute magentization value of the model 
        """

        return np.square(np.abs(np.sum(self.TFIM)/(self.SystemSize*self.ImaginaryTimeAxisLength))) 
    
    def GetAbsolutemagentizationPowerFour(self):
        """
        This function is to return the m^4 value of the 1D transverse field Ising 
        model.

        @return m: Absolute magentization value of the model 
        """

        return np.power(np.sum(self.TFIM)/(self.SystemSize*self.ImaginaryTimeAxisLength),4)

    def Seth_MagneticFieldStrength(self,H):
        """
        This function is to set the magnetic field strength (h) of the 1D transverse field 
        Ising model. 

        @param h: the magnetic field strength
        """

        self.h = H
        self.gamma = -0.5*np.log(np.tanh(self.delta_tau*H))
        self.J_y = (self.gamma)/(self.beta_cl)

        print("The h is reset.")
        print("h is updated to be :",self.h)
        print("The gamma is set to be :",self.gamma)
        print("The J_y is set to be :",self.J_y)
    
    def Geth_MagneticFieldStrength(self):
        """
        This function is to get the magnetic field strength (h) of the 1D transverse field 
        Ising model. 

        @return h: the magnetic field strength
        """

        return self.h
    
    # Use list to implement Wolff algo 
    def WolffUpdateBFSList(self):
        """
        This function is to update the 1D transverse field Ising model with 
        Wolff algorithm. It will change the self.TFIM (np.array) afterwards. 

        This function will be mainly using list to implement the queue 

        """
        from collections import deque

        # 1. Initialize the queue, visited set, cluster, the randomly choosen site and 
        # calculate the probabilities of adding the site to the queue 
        
        select_ReL = np.random.randint(0,self.SystemSize)
        select_ImL = np.random.randint(0,self.ImaginaryTimeAxisLength)
        
        queue = deque([(select_ReL, select_ImL)])
        cluster = set([(select_ReL,select_ImL)])

        ReL_Prob = 1 - np.exp(-2*self.J_x*self.beta_cl)
        ImL_Prob = 1 - np.exp(-2*np.abs(self.J_y)*self.beta_cl)

        #2. Keep on growing the cluster when the queue is not empty 
        while queue:
            current_ReL,current_ImL= queue.popleft()
            SpinAtThisSite = self.TFIM[current_ReL][current_ImL]
            tempNeighbourList = self.neighbour[current_ReL][current_ImL]

            for tempOrientation in range(0,4):
                NeighbourReL, NeighbourImL = tempNeighbourList[tempOrientation]
                SpinAtNeighbourSite = self.TFIM[NeighbourReL][NeighbourImL]

                if ((NeighbourReL, NeighbourImL) not in cluster) and \
                   (SpinAtNeighbourSite == SpinAtThisSite):
                    
                    if tempOrientation < 2 and (np.random.rand() < ReL_Prob):
                        cluster.add((NeighbourReL, NeighbourImL))
                        queue.append((NeighbourReL, NeighbourImL))
                    
                    elif tempOrientation >= 2 and (np.random.rand() < ImL_Prob):
                        cluster.add((NeighbourReL, NeighbourImL))
                        queue.append((NeighbourReL, NeighbourImL))

        #3. Flip the spin in cluster 
        for FlipReL,FlipImL in cluster:
            self.TFIM[FlipReL][FlipImL]*=-1 

    def WolffUpdateBFSListTest(self):
        """
        This function is to update the 1D transverse field Ising model with 
        Wolff algorithm. It will change the self.TFIM (np.array) afterwards. 

        This function will be mainly using list to implement the queue 

        """
        from collections import deque

        # 1. Initialize the queue, visited set, cluster, the randomly choosen site and 
        # calculate the probabilities of adding the site to the queue 
        
        select_ReL = np.random.randint(0,self.SystemSize)
        select_ImL = np.random.randint(0,self.ImaginaryTimeAxisLength)
        
        queue = deque([(select_ReL,select_ImL)])
        cluster = set([(select_ReL,select_ImL)])

        ReL_Prob = 1 - np.exp(-2*self.J_x*self.beta_cl)
        ImL_Prob = 1 - np.exp(-2*np.abs(self.J_y)*self.beta_cl)
        SpinAtThisSite = self.TFIM[select_ReL][select_ImL]
        notvisited = np.ones((self.SystemSize,self.ImaginaryTimeAxisLength))

        #2. Keep on growing the cluster when the queue is not empty 
        while queue:
            current_ReL,current_ImL= queue.popleft()
            tempNeighbourList = self.neighbour[current_ReL][current_ImL]

            for tempOrientation in range(0,4):
                NeighbourReL, NeighbourImL = tempNeighbourList[tempOrientation]
                SpinAtNeighbourSite = self.TFIM[NeighbourReL][NeighbourImL]

                if (notvisited[NeighbourReL][NeighbourImL]) and \
                   (SpinAtNeighbourSite == SpinAtThisSite):
                    
                    if tempOrientation < 2 and (np.random.rand() < ReL_Prob):
                        cluster.add((NeighbourReL, NeighbourImL))
                        queue.append((NeighbourReL, NeighbourImL))
                        notvisited[NeighbourReL][NeighbourImL] = 0 
                        self.TFIM[NeighbourReL][NeighbourImL]*=-1
                    
                    elif tempOrientation >= 2 and (np.random.rand() < ImL_Prob):
                        cluster.add((NeighbourReL, NeighbourImL))
                        queue.append((NeighbourReL, NeighbourImL))
                        notvisited[NeighbourReL][NeighbourImL]=0
                        self.TFIM[NeighbourReL][NeighbourImL]*=-1
                    
                    #else:
                    #    notvisited[NeighbourReL][NeighbourImL]=0

        #3. Flip the spin in cluster 
        #for FlipReL,FlipImL in cluster:
        #    self.TFIM[FlipReL][FlipImL]*=-1 

    def MetropolisUpdate(self):
        """
        This function is to update the 1D transverse field Ising model with 
        Metropolis algorithm. It will change the self.TFIM (np.array) afterwards. 

        Steps 
        1. Randomly choose a Im-L layer and L number 
        2. check neighbour energy (same layer & different layer)
        3. decide whether you will update
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
        
            SumOfSpinInSameLayer = self.TFIM[LeftReLCoor][ImTimeAxisIndexRandom] + \
                                   self.TFIM[RightReLCoor][ImTimeAxisIndexRandom]
            SumOfSpinInDiffLayer = self.TFIM[SystemSizeIndexRandom][UpImLCoor] + \
                                   self.TFIM[SystemSizeIndexRandom][DownImLCoor]

            #3. Determine whether to flip or not  
            DeltaETotal = 2*SpinAtThisSite*(self.J_x*SumOfSpinInSameLayer + 
                                            self.J_y*SumOfSpinInDiffLayer)
            
            weight = np.exp(-DeltaETotal*self.beta_cl)

            if (DeltaETotal <= 0) or (np.random.rand() < weight):
                self.TFIM[SystemSizeIndexRandom][ImTimeAxisIndexRandom] *= -1 

    """
    Below are the function used for debugging 
    """
    
    def __str__(self):
        """
        This function is to print the 1D TFIM as string 

        @return str(self.TFIM): Return the 1D TFIM as a string 
        """

        return str(self.TFIM)

    def get_all_variables(self):
        """
        This function would print all the instance variables of this class 
        
        """
        variables = vars(self)
        for name, value in variables.items():
            print(name, "is : \n", value)
    
    def get_variable(self,name):
        """
        This function would print all the instance variables of this class 
        @param : name of that variables
        
        """
        variable_value = vars(self).get(name)
        print(f"{name} is set to be: \n",variable_value)
    
    def print_boundary_neigbour(self,ImL_index):
        """
        This function will print the neighbour (only up down) of the configuration in certain
        ImaginaryTimeAxis Level
        @param ImL: The required ImaginaryTimeAxisLevel
        
        """

        for ReL_index in range(0,self.SystemSize):
            print(f"Site ReL = {ReL_index}, ImL = {ImL_index}")
            print(f"Up coordinates are ({self.neighbour[ReL_index][ImL_index][2][0]},{self.neighbour[ReL_index][ImL_index][2][1]})")
            print(f"Down coordinates are ({self.neighbour[ReL_index][ImL_index][3][0]},{self.neighbour[ReL_index][ImL_index][3][1]})")
        print()