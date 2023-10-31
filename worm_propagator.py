import math
import numpy
import random

class single_propagator():
    def __init__(self,n_time_slice,starting_x,starting_y,starting_z,starting_time):
        self.trajectory_x = [None for i in range(n_time_slice)]
        self.trajectory_y = [None for i in range(n_time_slice)]
        self.trajectory_z = [None for i in range(n_time_slice)]
        self.trajectory_x[starting_time]=starting_x
        self.trajectory_y[starting_time]=starting_y
        self.trajectory_z[starting_time]=starting_z

class worm_propagator():
    def __init__(self, N, n_time_slice,chemical_potential,epsilon, t):
        self.dist = []
        self.nbosons = 0
        self.N = N
        self.t = t
        self.n_time_slice = n_time_slice
        self.starting_point_x = None
        self.starting_point_y = None
        self.starting_point_z = None
        self.starting_time = None
        self.direction = None 
        self.epsilon = epsilon
        self.chemical_potential = chemical_potential
        # direction: -1 for backward propagator;
        # 1 for forward propagator.
        
        # initialize the distribution of bosons #
    def judge(self,current_time,current_x,current_y,current_z):
        for propagator in self.dist[:-2]:
            if propagator.trajectory_x[current_time] == current_x and \
                propagator.trajectory_y[current_time] == current_y and \
                propagator.trajectory_z[current_time] == current_z :
                return propagator
        return None
    
    def find_and_swap(self, current_x, current_y, current_z):
        for npropagator in len(self.dist[:-1]):
            if (self.dist[npropagator].trajectory_x[self.n_time_slice-1] == current_x) and (self.dist[npropagator].trajectory_y[self.n_time_slice-1] == current_y) and \
                (self.dist[npropagator].trajectory_z[self.n_time_slice-1] == current_z):
                    self.dist[npropagator], self.dist[-1] = self.dist[-1], self.dist[npropagator]
        return 0
    
    def prop_all(self,propagator):
        for i in range(self.n_time_slice):
            if propagator.trajectory_x[i] != None:
                return False
        return True
    
    def prop_any(self,propagator):
        for i in range(self.n_time_slice):
            if propagator.trajectory_x[i] == None:
                return True 
        return False
     
    def renormalize(self):
        first_propagator = None
        for propagator in range(len(self.dist)-1):
            if self.prop_all(self.dist[propagator]):
                del self.dist[propagator]
            elif self.prop_any(self.dist[propagator]) :
                first_propagator = self.dist[propagator]
        if first_propagator != None :
            for time in range(self.n_time_slice):
                if first_propagator.trajectory_x[time] == None :
                    first_propagator.trajectory_x[time] = self.dist[-1].trajectory_x[time]
                    first_propagator.trajectory_y[time] = self.dist[-1].trajectory_y[time]
                    first_propagator.trajectory_z[time] = self.dist[-1].trajectory_z[time]
            del self.dist[-1]
        return 0
    
    def compare(self,time,position_x,position_y,position_z):
        for propagator in self.dist:
                if propagator.trajectory_x[time] == position_x \
                and propagator.trajectory_y[time] == position_y \
                and propagator.trajectory_z[time] == position_z :
                    #Then it overlaps with the current particle#
                    return False
        return True 
    
    def initialization(self):
        st_point = random.randrange(0,self.N**3*self.n_time_slice,1)
        self.starting_time = st_point // (self.N**3) # find the correct time slice
        position = st_point % (self.N**3)
        self.starting_point_z = position // (self.N**2)
        position = position % (self.N**2)
        self.starting_point_y = position // self.N
        self.starting_point_x = position % self.N
        if self.dist != []:
            temp2 = self.judge(self.starting_time, self.starting_point_x, self.starting_point_y, self.starting_point_z)
            if temp2 == None:
                self.dist.append(single_propagator(self.n_time_slice,
                                               self.starting_point_x,self.starting_point_y,self.starting_point_z,self.starting_time))
                if random.random() <= min(1,math.exp(self.epsilon*self.chemical_potential)):
                    self.direction = +1
            else:
                self.direction = -1
                temp2, self.dist[-1] = self.dist[-1], temp2
                #Execute the steps for overlapping. I would rather consider it later. #                                     
        else:
            self.dist.append(single_propagator(self.n_time_slice,self.starting_point_x,
                                               self.starting_point_y,self.starting_point_z,self.starting_time))
            if random.random() <= min(1,math.exp(self.epsilon*self.chemical_potential)):
                self.direction = +1
            else:
                self.direction = -1
        
    def propagator(self, current_time, current_x,current_y,current_z,init=None): 
        # To make sure the path is feasible after MC update.
        # In other words, to make sure we are actually calculating the 
        # diagonal elements of the propagator.
        if (self.direction == 1) :
            if (current_time == self.starting_time) and (current_x == self.starting_point_x) and \
                (current_y == self.starting_point_y) and (current_z == self.starting_point_z) and (init == None):
                return 0
            temp1 = self.judge(current_time, current_x, current_y, current_z)
            if (temp1 != None):
                self.direction = -1
                self.dist[-1].trajectory_x[current_time:-1] = temp1.trajectory_x[current_time:-1]
                self.dist[-1].trajectory_y[current_time:-1] = temp1.trajectory_y[current_time:-1]
                self.dist[-1].trajectory_z[current_time:-1] = temp1.trajectory_z[current_time:-1]
                temp1.trajectory_x[current_time:-1] = None
                temp1.trajectory_y[current_time:-1] = None
                temp1.trajectory_z[current_time:-1] = None
                temp1, self.dist[-1] = self.dist[-1], temp1
                self.propagator(current_time, self.dist[-1].trajectory_x[current_time],
                                self.dist[-1].trajectory_y[current_time],self.dist[-1].trajectory_z[current_time])
                
            elif (current_time == self.n_time_slice-1):
                self.dist.append(single_propagator(self.n_time_slice,
                                current_x,current_y,current_z, 0)) # bond update
                current_time = 0
                if (random.random() <= min(1,math.exp(self.epsilon*self.chemical_potential))):
                    self.direction = +1
                else:
                    self.direction = -1
                    # site update
                    self.propagator(current_time,self.dist[-1].trajectory_x[current_time],
                                self.dist[-1].trajectory_y[current_time],self.dist[-1].trajectory_z[current_time])
            else:
                current_time += 1
                randomnum = random.random()
                if randomnum <= 1-6*self.t*self.epsilon :
                    self.dist[-1].trajectory_x[current_time] = current_x
                    self.dist[-1].trajectory_y[current_time] = current_y
                    self.dist[-1].trajectory_z[current_time] = current_z
                elif randomnum <= 1-5*self.t*self.epsilon:
                    self.dist[-1].trajectory_x[current_time] = (current_x + 1) % self.N
                    self.dist[-1].trajectory_y[current_time] = current_y
                    self.dist[-1].trajectory_z[current_time] = current_z
                elif randomnum <= 1-4*self.t*self.epsilon:
                    self.dist[-1].trajectory_x[current_time] = current_x
                    self.dist[-1].trajectory_y[current_time] = (current_y + 1) % self.N
                    self.dist[-1].trajectory_z[current_time] = current_z
                elif randomnum <= 1-3*self.t*self.epsilon:
                    self.dist[-1].trajectory_x[current_time] = current_x
                    self.dist[-1].trajectory_y[current_time] = current_y
                    self.dist[-1].trajectory_z[current_time] = (current_z + 1) % self.N
                elif randomnum <= 1-2*self.t*self.epsilon:
                    self.dist[-1].trajectory_x[current_time] = (current_x - 1 + self.N) % self.N
                    self.dist[-1].trajectory_y[current_time] = current_y
                    self.dist[-1].trajectory_z[current_time] = current_z
                elif randomnum <= 1-self.t*self.epsilon:
                    self.dist[-1].trajectory_x[current_time] = current_x
                    self.dist[-1].trajectory_y[current_time] = (current_y - 1 + self.N) % self.N
                    self.dist[-1].trajectory_z[current_time] = current_z
                elif randomnum <= 1:
                    self.dist[-1].trajectory_x[current_time] = current_x
                    self.dist[-1].trajectory_y[current_time] = current_y
                    self.dist[-1].trajectory_z[current_time] = (current_z - 1 + self.N) % self.N
                if (random.random() <= min(1,math.exp(self.epsilon*self.chemical_potential))):
                    self.direction = +1
                else:
                    self.direction = -1
                    # site update
                self.propagator(current_time, self.dist[-1].trajectory_x[current_time],
                                self.dist[-1].trajectory_y[current_time],self.dist[-1].trajectory_z[current_time]) # bond update
                    
        if (self.direction == -1):
            self.dist[-1].trajectory_x[current_time] = None
            self.dist[-1].trajectory_y[current_time] = None
            self.dist[-1].trajectory_z[current_time] = None
            if (current_x == self.starting_point_x) and (current_y == self.starting_point_y) and \
                (current_z == self.starting_point_z) and (current_time == self.starting_time):
                    return 0 # Finished #
            if (random.random() <= min(1,math.exp(-self.epsilon*self.chemical_potential))):
                self.direction = -1
            else:
                self.direction = 1
            if current_time > 0 :
                current_time = current_time - 1
                self.propagator(current_time, self.dist[-1].trajectory_x[current_time],
                            self.dist[-1].trajectory_y[current_time],self.dist[-1].trajectory_z[current_time])
            elif current_time == 0:
                self.find_and_swap(current_x, current_y, current_z)
                self.propagator(self.N-1, self.dist[-1].trajectory_x[self.n_time_slice-1],
                                self.dist[-1].trajectory_y[self.n_time_slice-1], self.dist[-1].trajectory_z[self.n_time_slice-1])
        
    def nbosons_sampler(self):
        self.renormalize()
        self.nbosons = len(self.dist)
        return self.nbosons
    
    def energy_sampler(self):
        return 0 # needs work

if __name__ == "__main__":
    
    worm = worm_propagator(2,50,0,0.002,1)
    for i in range(100):
        worm.initialization()
        print(worm.nbosons_sampler())
        worm.propagator(worm.starting_time, worm.starting_point_x,worm.starting_point_y,worm.starting_point_z,init=1)
        print(worm.nbosons_sampler())
    print(worm.dist[0].trajectory_x)
    print(worm.dist[0].trajectory_y)
    print(worm.dist[0].trajectory_z)
