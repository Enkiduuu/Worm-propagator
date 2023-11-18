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
        self.current_x = None
        self.current_y = None 
        self.current_z = None 
        self.current_time = None
        # direction: -1 for backward propagator;
        # 1 for forward propagator.
        
        # initialize the distribution of bosons #
    def judge(self,current_time,current_x,current_y,current_z):
        for propagator in range(len(self.dist)-1):
            if self.dist[propagator].trajectory_x[current_time] == current_x and \
                self.dist[propagator].trajectory_y[current_time] == current_y and \
                self.dist[propagator].trajectory_z[current_time] == current_z :
                return propagator
        return None
    
    def judge2(self,current_time,current_x,current_y,current_z):
        for propagator in range(len(self.dist)):
            if self.dist[propagator].trajectory_x[current_time] == current_x and \
                self.dist[propagator].trajectory_y[current_time] == current_y and \
                self.dist[propagator].trajectory_z[current_time] == current_z :
                return propagator
        return None
    
    def find_and_swap(self, current_x, current_y, current_z):
        for npropagator in range(len(self.dist[:-1])):
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
        dist2 = []
        for propagator in range(len(self.dist)-1):
            if self.prop_all(self.dist[propagator]):
                pass
            elif self.prop_any(self.dist[propagator]) :
                first_propagator = self.dist[propagator]
            else:
                dist2.append(self.dist[propagator])
        if first_propagator != None :
            #print(first_propagator.trajectory_x, self.dist[-1].trajectory_x)
            for time in range(self.n_time_slice):
                if first_propagator.trajectory_x[time] == None :
                    first_propagator.trajectory_x[time] = self.dist[-1].trajectory_x[time]
                    first_propagator.trajectory_y[time] = self.dist[-1].trajectory_y[time]
                    first_propagator.trajectory_z[time] = self.dist[-1].trajectory_z[time]
            dist2.append(first_propagator)
        elif self.prop_all(self.dist[-1]):
            pass 
        elif not self.prop_any(self.dist[-1]):
            dist2.append(self.dist[-1])
        self.dist = dist2.copy()
        return 0
    
    def initialization(self):
        st_point = random.randrange(0,self.N**3*self.n_time_slice,1)
        self.starting_time = st_point // (self.N**3) # find the correct time slice
        position = st_point % (self.N**3)
        self.starting_point_z = position // (self.N**2)
        position = position % (self.N**2)
        self.starting_point_y = position // self.N
        self.starting_point_x = position % self.N
        self.current_x = self.starting_point_x 
        self.current_y = self.starting_point_y
        self.current_z = self.starting_point_z
        self.current_time = self.starting_time
        if self.dist != []:
            temp2 = self.judge2(self.starting_time, self.starting_point_x, self.starting_point_y, self.starting_point_z)
            if temp2 == None:
                self.dist.append(single_propagator(self.n_time_slice,
                                               self.starting_point_x,self.starting_point_y,self.starting_point_z,self.starting_time))
                if random.random() <= min(1,math.exp(self.epsilon*self.chemical_potential)):
                    self.direction = +1
                else:
                    self.direction = -1
            else:
                self.direction = -1
                self.dist.append(single_propagator(self.n_time_slice,None,None,None,self.starting_time))
                self.dist[-1].trajectory_x[:self.starting_time] = self.dist[temp2].trajectory_x[:self.starting_time]
                self.dist[-1].trajectory_y[:self.starting_time] = self.dist[temp2].trajectory_y[:self.starting_time]
                self.dist[-1].trajectory_z[:self.starting_time] = self.dist[temp2].trajectory_z[:self.starting_time]
                self.dist[temp2].trajectory_x[:self.starting_time] = [None for i in range(self.starting_time)] # Trying
                self.dist[temp2].trajectory_y[:self.starting_time] = [None for i in range(self.starting_time)] # Trying
                self.dist[temp2].trajectory_z[:self.starting_time] = [None for i in range(self.starting_time)] # Trying
                             
        else:
            self.dist.append(single_propagator(self.n_time_slice,self.starting_point_x,
                                               self.starting_point_y,self.starting_point_z,self.starting_time))
            if random.random() <= min(1,math.exp(self.epsilon*self.chemical_potential)):
                self.direction = +1
            else:
                self.direction = -1

        return 0
        
    def bond_update(self):
        if self.direction == 1:
            if self.current_time < self.n_time_slice-1:
                self.current_time += 1
                randomnum = random.random()
                if randomnum <= 1-6*self.t*self.epsilon :
                    self.dist[-1].trajectory_x[self.current_time] = self.current_x
                    self.dist[-1].trajectory_y[self.current_time] = self.current_y
                    self.dist[-1].trajectory_z[self.current_time] = self.current_z
                elif randomnum <= 1-5*self.t*self.epsilon:
                    self.dist[-1].trajectory_x[self.current_time] = (self.current_x + 1) % self.N
                    self.dist[-1].trajectory_y[self.current_time] = self.current_y
                    self.dist[-1].trajectory_z[self.current_time] = self.current_z
                elif randomnum <= 1-4*self.t*self.epsilon:
                    self.dist[-1].trajectory_x[self.current_time] = self.current_x
                    self.dist[-1].trajectory_y[self.current_time] = (self.current_y + 1) % self.N
                    self.dist[-1].trajectory_z[self.current_time] = self.current_z
                elif randomnum <= 1-3*self.t*self.epsilon:
                    self.dist[-1].trajectory_x[self.current_time] = self.current_x
                    self.dist[-1].trajectory_y[self.current_time] = self.current_y
                    self.dist[-1].trajectory_z[self.current_time] = (self.current_z + 1) % self.N
                elif randomnum <= 1-2*self.t*self.epsilon:
                    self.dist[-1].trajectory_x[self.current_time] = (self.current_x - 1 + self.N) % self.N
                    self.dist[-1].trajectory_y[self.current_time] = self.current_y
                    self.dist[-1].trajectory_z[self.current_time] = self.current_z
                elif randomnum <= 1-self.t*self.epsilon:
                    self.dist[-1].trajectory_x[self.current_time] = self.current_x
                    self.dist[-1].trajectory_y[self.current_time] = (self.current_y - 1 + self.N) % self.N
                    self.dist[-1].trajectory_z[self.current_time] = self.current_z
                elif randomnum <= 1:
                    self.dist[-1].trajectory_x[self.current_time] = self.current_x
                    self.dist[-1].trajectory_y[self.current_time] = self.current_y
                    self.dist[-1].trajectory_z[self.current_time] = (self.current_z - 1 + self.N) % self.N
                self.current_x = self.dist[-1].trajectory_x[self.current_time]
                self.current_y = self.dist[-1].trajectory_y[self.current_time]
                self.current_z = self.dist[-1].trajectory_z[self.current_time]
            elif self.current_time == self.n_time_slice -1:
                self.current_time = 0
                self.dist.append(single_propagator(self.n_time_slice,self.current_x,self.current_y,
                                self.current_z, self.current_time)) # n_time_slice,starting_x,starting_y,starting_z,starting_time
        elif self.direction == -1:
            self.dist[-1].trajectory_x[self.current_time] = None 
            self.dist[-1].trajectory_y[self.current_time] = None 
            self.dist[-1].trajectory_z[self.current_time] = None
            if self.current_time > 0:
                self.current_time -= 1
            elif self.current_time == 0 :
                self.current_time = self.n_time_slice-1 
                self.find_and_swap(self.current_x,self.current_y,self.current_z)
            self.current_x = self.dist[-1].trajectory_x[self.current_time]
            self.current_y = self.dist[-1].trajectory_y[self.current_time]
            self.current_z = self.dist[-1].trajectory_z[self.current_time]  
        return 0 
    
    def site_update(self):
        if self.direction == -1:
            if random.random() < min(1,math.exp(-self.epsilon*self.chemical_potential)):
                self.direction = -1
            else:
                self.direction = +1
        elif self.direction == +1:
            if random.random() < min(1,math.exp(self.epsilon*self.chemical_potential)):
                self.direction = +1
            else:
                self.direction = -1
        return 0
    
    def propagator_new(self):
        judgement = True
        #print(self.starting_time, self.starting_point_x,self.starting_point_y,self.starting_point_z)
        while judgement:
            self.bond_update()
            #print(self.current_time,self.current_x,self.current_y,self.current_z)
            if (self.current_x == None):
                print("Error! None of current_x!")
                break
            if (self.direction == 1) and (self.current_time == self.starting_time) and (self.current_x == self.starting_point_x) \
                and (self.current_y == self.starting_point_y) and (self.current_z == self.starting_point_z):
                    judgement = False
            elif self.direction == 1:
                prop = self.judge(self.current_time,self.current_x,self.current_y,self.current_z) #current_time,current_x,current_y,current_z
                if prop == None:
                    self.site_update()
                else:
                    self.dist[-1].trajectory_x[self.current_time:], self.dist[prop].trajectory_x[self.current_time:] = \
                    self.dist[prop].trajectory_x[self.current_time:], self.dist[-1].trajectory_x[self.current_time:]
                    
                    self.dist[-1].trajectory_y[self.current_time:], self.dist[prop].trajectory_y[self.current_time:] = \
                    self.dist[prop].trajectory_y[self.current_time:], self.dist[-1].trajectory_y[self.current_time:]
                    
                    self.dist[-1].trajectory_z[self.current_time:], self.dist[prop].trajectory_z[self.current_time:] = \
                    self.dist[prop].trajectory_z[self.current_time:], self.dist[-1].trajectory_z[self.current_time:]
                    
                    self.dist[prop], self.dist[-1] = self.dist[-1], self.dist[prop]
                    
                    self.direction = -1
                    
            elif self.direction == -1:
                self.site_update()
                if (self.direction == -1) and (self.current_time == self.starting_time) and (self.current_x == self.starting_point_x) \
                and (self.current_y == self.starting_point_y) and (self.current_z == self.starting_point_z):
                    
                    judgement = False
                    self.dist[-1].trajectory_x[self.current_time] = None
                    self.dist[-1].trajectory_y[self.current_time] = None 
                    self.dist[-1].trajectory_z[self.current_time] = None
        return 0
    
    def energy_sampler(self):
        return 0 # needs work

    def evaluate_deviation(self,aveN): 
        aveN = aveN[1:]
        print(aveN)
        ave = sum(aveN)/len(aveN)
        for i in range(len(aveN)):
            aveN[i] -= ave 
            aveN[i] = aveN[i]*aveN[i]
        deviation = math.sqrt(sum(aveN)/len(aveN))
        print(ave, deviation/math.sqrt(len(aveN)))
        return 0
    
    def check1(self):
        for propagator in self.dist:
            for i in range(self.n_time_slice):
                if propagator.trajectory_x[i] == None:
                    return False
        return True
    
    def check2(self):
        for propagator in self.dist:
            judge = False
            for propagator2 in self.dist:
                if (propagator.trajectory_x[0]==propagator2.trajectory_x[self.n_time_slice-1]) and \
                    (propagator.trajectory_y[0]==propagator2.trajectory_y[self.n_time_slice-1]) and \
                    (propagator.trajectory_z[0]==propagator2.trajectory_z[self.n_time_slice-1]):
                    judge = True
            if judge == False:
                return False 
        return True

if __name__ == "__main__":
    Ntimes = int(input())
    aveN = []
    worm = worm_propagator(2,12001,1.4,0.001,1) # N, n_time_slice,chemical_potential,epsilon, t
    for i in range(Ntimes):
        aveone = []
        for j in range(5000):
            worm.initialization()
            worm.propagator_new()
            #print("---------------Before-------------")
            #for k in range(len(worm.dist)):
            #    print(worm.dist[k].trajectory_x)
            #    print(worm.dist[k].trajectory_y)
            #    print(worm.dist[k].trajectory_z)
            worm.renormalize()
            check_worm = worm.check1()
            check_worm_2 = worm.check2()
            print(len(worm.dist),j/5000,i)
            aveone.append(len(worm.dist))
            #print("---------------After-------------")
            #for k in range(len(worm.dist)):
            #    print(worm.dist[k].trajectory_x)
            #    print(worm.dist[k].trajectory_y)
            #    print(worm.dist[k].trajectory_z)
            if not check_worm_2:
                print("Error! Not closed!")
                break
            if not check_worm:
                print("Error! None! ")
                break
        aveN.append(sum(aveone)/len(aveone))
    print(aveN)
    worm.evaluate_deviation(aveN)
