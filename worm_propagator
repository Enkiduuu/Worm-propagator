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
    def __init__(self, N, n_time_slice):
        self.dist = []
        self.nbosons = 0
        self.N = N
        self.n_time_slice = n_time_slice
        self.starting_point_x = None
        self.starting_point_y = None
        self.starting_point_z = None
        self.starting_time = None
        self.direction = None 
        self.epsilon = None
        self.chemical_potential = None 
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
        for npropagator in len(self.dist[:-2]):
            if (self.dist[npropagator].trajectory_x[0] == current_x) and (self.dist[npropagator].trajectory_y[0] == current_y) and \
                (self.dist[npropagator].trajectory_z[0] == current_z):
                    self.dist[npropagator], self.dist[-1] = self.dist[-1], self.dist[npropagator]
        return 0
    
    def renormalize(self):
        first_propagator = None
        second_propagator = None
        for propagator in self.dist[:-1]:
            if all(propagator.trajectory_x == None):
                del propagator
            if any(propagator.trajectory_x == None) and (first_propagator == None):
                first_propagator = propagator
            if any(propagator.trajectory_x == None) and (first_propagator != None):
                second_propagator = propagator
        if (first_propagator != None) and (second_propagator != None):
            for time in range(self.N):
                if first_propagator.trajectory_x[time] == None :
                    first_propagator.trajectory_x[time] = second_propagator.trajectory_x[time]
                    first_propagator.trajectory_y[time] = second_propagator.trajectory_y[time]
                    first_propagator.trajectory_z[time] = second_propagator.trajectory_z[time]
            del second_propagator
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
            if self.compare(self.starting_time, self.starting_point_x, self.starting_point_y, self.starting_point_z):
                self.dist.append(single_propagator(self.n_time_slice,
                                               self.starting_point_x,self.starting_point_y,self.starting_point_z,self.starting_time))
                self.nbosons +=1
                if random.random() <= min(1,math.exp(self.epsilon*self.chemical_potential)):
                    self.direction = +1
            else:
                self.nbosons -= 1
                pass #Execute the steps for overlapping. I would rather consider it later. #                 
                                    
        else:
            self.dist.append(single_propagator(self.n_time_slice,
                                               self.starting_point_x,self.starting_point_z,self.starting_time))
            if random.random() <= min(1,math.exp(self.epsilon*self.chemical_potential)):
                self.direction = +1
            self.nbosons += 1
        
    def propagator(self, current_time, current_x,current_y,current_z): 
        # To make sure the path is feasible after MC update.
        # In other words, to make sure we are actually calculateing the 
        # diagonal elements of the propagator.
        if (self.direction == 1) :
            randomnum = random.random()
            if (current_time == self.starting_time) and (current_x == self.starting_point_x) and \
                (current_y == self.starting_point_y) and (current_z == self.starting_point_z):
                    return 0
            if (current_time == self.N-1):
                self.dist.append(single_propagator(self.n_time_slice,
                                current_x,current_y,current_z, 0)) # bond update
                if (random.random() <= min(1,math.exp(self.epsilon*self.chemical_potential))):
                    self.direction = +1
                else:
                    self.direction = -1
                    # site update
            else:
                current_time = (current_time + 1) % self.N
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
                self.propagator(self, current_time, self.dist[-1].trajectory_x[current_time],
                                self.dist[-1].trajectory_y[current_time],self.dist[-1].trajectory_z[current_time]) # bond update
                if (random.random() <= min(1,math.exp(self.epsilon*self.chemical_potential))):
                    self.direction = +1
                else:
                    self.direction = -1
                    # site update
                    
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
                self.propagator(self, current_time, self.dist[-1].trajectory_x[current_time],
                            self.dist[-1].trajectory_y[current_time],self.dist[-1].trajectory_z[current_time])
            elif current_time == 0:
                self.find_and_swap(self, current_x, current_y, current_z)
                self.propagator(self, current_time, self.dist[-1].trajectory_x[self.N-1],
                                self.dist[-1].trajectory_y[self.N-1], self.dist[-1].trajectory_z[self.N-1])
        
    def nbosons_sampler(self):
        self.nbosons = len(self.dist)
        return self.nbosons
    
    def energy_sampler(self):
        return 0 # needs work
