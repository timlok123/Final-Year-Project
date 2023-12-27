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
        self.TFIM = np.random.choice([1,-1],size=(RealAxisLength, \
                                                  ImaginaryTimeAxisLength))

        # set h to be 1
        self.h = 1      

        self.delta_tau = 0.1        
        self.beta_cl = self.delta_tau

        self.gamma = -0.5*np.log(np.tanh(self.delta_tau*self.h))
        self.J_x = 1 
        self.J_y = (self.gamma)/(self.beta_cl)

        self.neighbour = np.zeros((RealAxisLength,ImaginaryTimeAxisLength, 4, 2),dtype=int) 
        # 4 directions, 2 coordinates (x,y)

        self.TotalMCsteps = RealAxisLength*ImaginaryTimeAxisLength 
        # 1 MC steps = 2D Ising model size


    def SetDetlaTau(self,Dtau):
        self.delta_tau = Dtau
        self.beta_cl = Dtau
        print("The delta tau is set to be :", self.delta_tau)
        print("The beta_cl is set to be :", self.beta_cl)

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
    
    # Use list to implement Wolff algo 
    def WolffUpdateBFSList(self):
        """
        This function is to update the 1D transverse field Ising model with 
        Wolff algorithm. It will change the self.TFIM (np.array) afterwards. 

        This function will be mainly using list to implement the queue 

        """

        #1. Initialize the queue, visited set, cluster, the randomly choosen site and 
        # calculate the probabilities of adding the site to the queue 
        
        select_ReL = np.random.randint(0,self.SystemSize)
        select_ImL = np.random.randint(0,self.ImaginaryTimeAxisLength)
        
        queue = [(select_ReL,select_ImL)]
        cluster = [(select_ReL,select_ImL)]
        visited = set()
        visited.add((select_ReL,select_ImL))

        ReL_Prob = 1 - np.exp(-2*self.J_x*self.beta_cl)
        ImL_Prob = 1 - np.exp(-2*np.abs(self.J_y)*self.beta_cl)

        #2. Keep on growing the cluster when the queue is not empty 
        while len(queue)!=0:
            current_ReL,current_ImL= queue.pop(0)
            SpinAtThisSite = self.TFIM[current_ReL][current_ImL]

            tempNeighbourList = self.neighbour[current_ReL][current_ImL]

            # left-right 
            for tempOrientation in range(0,2):

                NeighbourReL,NeighbourImL = tempNeighbourList[tempOrientation]
                SpinAtNeighbourSite = self.TFIM[NeighbourReL][NeighbourImL]

                if ((NeighbourReL,NeighbourImL) not in visited) and \
                   (SpinAtNeighbourSite == SpinAtThisSite) and\
                   (np.random.rand() < ReL_Prob):
                    cluster.append((NeighbourReL,NeighbourImL))
                    queue.append((NeighbourReL,NeighbourImL))
                
                visited.add((NeighbourReL,NeighbourImL))
                
            # up-down layer 
            for tempOrientation in range(2,4):

                NeighbourReL,NeighbourImL = tempNeighbourList[tempOrientation]
                SpinAtNeighbourSite = self.TFIM[NeighbourReL][NeighbourImL]

                if ((NeighbourReL,NeighbourImL) not in visited) and \
                   (SpinAtNeighbourSite == SpinAtThisSite) and\
                   (np.random.rand() < ImL_Prob):
                    cluster.append((NeighbourReL,NeighbourImL))
                    queue.append((NeighbourReL,NeighbourImL))
                
                visited.add((NeighbourReL,NeighbourImL))

        #3. Flip the spin in cluster 
        for site in cluster:
            ReL,ImL = site
            self.TFIM[ReL][ImL]*=-1 


    # Too tedious to implement. Give this up 
    def WolffUpdateArray(self):
        """
        This function is to update the 1D transverse field Ising model with 
        Wolff algorithm. It will change the self.TFIM (np.array) afterwards. 

        This function will be mainly using array to implement the code

        """

        select_ReL = np.random.randint(0,self.SystemSize)
        select_ImL = np.random.randint(0,self.ImaginaryTimeAxisLength)

        Pocket = np.zeros((self.ImaginaryTimeAxisLength*self.SystemSize,2),dtype=int)
        Cluster = np.zeros((self.ImaginaryTimeAxisLength*self.SystemSize,2),dtype=int)

        Pocket[0] = [select_ReL,select_ImL]
        PocketPointToTop = 1 
        PocketStart = 0 
        CluesterPointToTop = 0

        # Calculate the prob 
        ReL_Prob = 1 - np.exp(-2*self.J_x*self.beta_cl)
        ImL_Prob = 1 - np.exp(-2*np.abs(self.J_y)*self.beta_cl)

        while (PocketPointToTop - PocketStart) != 0: 
            randomIndex = PocketStart
            RandomReL, RandomImL = Pocket[randomIndex]
            SpinNow = self.TFIM[RandomReL][RandomImL]

            #for NeighbourX,NeighbourY in self.neighbour[RandomReL][RandomImL]:
            for count,NeighbourS in enumerate(self.neighbour[RandomReL][RandomImL]):
                NeighbourReL,NeighbourImL = NeighbourS 
                SpinNeighbour = self.TFIM[NeighbourReL][NeighbourImL]

                if (NeighbourS not in Cluster[0:CluesterPointToTop]) and \
                    (SpinNow == SpinNeighbour):

                    if ((count<2) and (np.random.random_sample()<ReL_Prob)) or \
                       ((count>=2) and (np.random.random_sample()<ImL_Prob)):
                        
                       Pocket[PocketPointToTop] = NeighbourS
                       Pocket[PocketStart] = [0,0]
                       Cluster[CluesterPointToTop] = NeighbourS
                       PocketStart = (PocketStart+1)%len(Pocket)
                       PocketPointToTop = (PocketPointToTop+1)%len(Pocket)  
                       CluesterPointToTop += 1
                    
                    else:
                        Pocket[PocketStart] = [0,0]
                        PocketStart = (PocketStart+1)%len(Pocket)

                else:
                    Pocket[PocketStart] = [0,0]
                    PocketStart = (PocketStart+1)%len(Pocket)
        
        for s in Cluster:
            ReL,ImL = s 
            self.TFIM[ReL][ImL] *= -1

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

            #2. Calculate the energy of its neighbour 
            SumOfSpinInSameLayer = 0 
            SumOfSpinInDiffLayer = 0 

            LeftIndex = (SystemSizeIndexRandom-1)%(self.SystemSize)
            RightIndex = (SystemSizeIndexRandom+1)%(self.SystemSize)
            DownIndex = (ImTimeAxisIndexRandom-1)%(self.ImaginaryTimeAxisLength)
            UpIndex = (ImTimeAxisIndexRandom+1)%(self.ImaginaryTimeAxisLength)
        
            SumOfSpinInSameLayer =  self.TFIM[LeftIndex][ImTimeAxisIndexRandom] + self.TFIM[RightIndex][ImTimeAxisIndexRandom]
            SumOfSpinInDiffLayer = self.TFIM[SystemSizeIndexRandom][DownIndex] + self.TFIM[SystemSizeIndexRandom][UpIndex]

            #3. Determine whether to flip or not  
            DeltaETotal = 2*self.TFIM[SystemSizeIndexRandom][ImTimeAxisIndexRandom]*(self.J_x*SumOfSpinInSameLayer + self.J_y*SumOfSpinInDiffLayer)
            
            weight = np.exp(-DeltaETotal*self.beta_cl)

            if (DeltaETotal < 0) or (np.random.rand() < weight):
                self.TFIM[SystemSizeIndexRandom][ImTimeAxisIndexRandom] *= -1 

    """
    Below are the function used for debugging 
    """
    
    def __str__(self):
        """
        This function is to print the 1D TFIM as string 

        @return TFIM(np.array): Return the 1D TFIM as a string 
        """

        return str(self.TFIM)

    def get_all_variables(self):
        """
        This function would print all the instance variables of this class 
        
        """
        variables = vars(self)
        for name, value in variables.items():
            print(name, "=", value)
    
    def get_variable(self,name):
        """
        This function would print all the instance variables of this class 
        @name : name of that variables
        
        """
        variable_value = vars(self).get(name)
        print(f"{name} is set to be: ",variable_value)
    
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