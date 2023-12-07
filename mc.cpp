#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>

class Node {
public:
    int x, y, z, t;
    bool status;
    Node *prev;
    Node *prev_judge;
    Node *next;

    Node(int x_val, int y_val, int z_val, int t_val) {
        x = x_val;
        y = y_val;
        z = z_val;
        t = t_val;
        prev = nullptr;
        prev_judge = nullptr;
        next = nullptr;
        status = false;
    }
};

class WormPropagator {
private:
    int nbosons;
    int N, t, n_time_slice;
    int starting_point_x, starting_point_y, starting_point_z, starting_time;
    int direction;
    double epsilon, chemical_potential;
    int current_x, current_y, current_z, current_time;
    Node *****dist;

public:
    WormPropagator(int N_val, int n_time_slice_val, double chemical_potential_val, double epsilon_val, int t_val) {
        N = N_val;
        t = t_val;
        n_time_slice = n_time_slice_val;
        chemical_potential = chemical_potential_val;
        epsilon = epsilon_val;
        // Initialize dist
        dist = new Node****[n_time_slice];
        for (int i = 0; i < n_time_slice; ++i) {
            dist[i] = new Node***[N];
            for (int j = 0; j < N; ++j) {
                dist[i][j] = new Node**[N];
                for (int k = 0; k < N; ++k) {
                    dist[i][j][k] = new Node*[N];
                    for (int l = 0; l < N; ++l) {
                        dist[i][j][k][l] = new Node(l, k, j, i);
                    }
                }
            }
        }
        nbosons = 0;
        starting_point_x = -1;
        starting_point_y = -1;
        starting_point_z = -1;
        starting_time = -1;
        direction = 0;
        current_x = -1;
        current_y = -1;
        current_z = -1;
        current_time = -1;
    }

    void initialization() {
        int st_point = rand() % (N * N * N * n_time_slice);
        starting_time = st_point / (N * N * N);
        int position = st_point % (N * N * N);
        starting_point_z = position / (N * N);
        position = position % (N * N);
        starting_point_y = position / N;
        starting_point_x = position % N;
        current_x = starting_point_x;
        current_y = starting_point_y;
        current_z = starting_point_z;
        current_time = starting_time;

        if (true) {
            if (!dist[starting_time][starting_point_z][starting_point_y][starting_point_x]->status) {
                dist[starting_time][starting_point_z][starting_point_y][starting_point_x]->status = true;
                if ((double)rand() / (double) RAND_MAX <= std::min(1.0, exp(epsilon * chemical_potential))) {
                    direction = 1;
                } else {
                    direction = -1;
                }
            } else {
                direction = -1;
                dist[starting_time][starting_point_z][starting_point_y][starting_point_x]->status = false;
            }
        } 
    }

    void bond_update() {
    if (direction == 1) {
        if (current_time < n_time_slice - 1) {
            current_time += 1;
            double randomnum = (double)rand() / (double) RAND_MAX;
            if (randomnum <= (double)1 - (double)6 * t * epsilon) {
                dist[current_time - 1][current_z][current_y][current_x]->next = dist[current_time][current_z][current_y][current_x];
                dist[current_time][current_z][current_y][current_x]->prev_judge = dist[current_time - 1][current_z][current_y][current_x];
            } 
            else if (randomnum <= (double)1 - (double)5 * t * epsilon) {
                dist[current_time - 1][current_z][current_y][current_x]->next = dist[current_time][current_z][current_y][(current_x + 1) % N];
                dist[current_time][current_z][current_y][(current_x + 1) % N]->prev_judge = dist[current_time - 1][current_z][current_y][current_x];
            } 
            else if (randomnum <= (double)1 - (double)4 * t * epsilon) {
                dist[current_time-1][current_z][current_y][current_x]->next = dist[current_time][current_z][(current_y+1)% N][current_x];
                dist[current_time][current_z][(current_y+1)% N][current_x]->prev_judge = dist[current_time-1][current_z][current_y][current_x];
            }
            else if (randomnum <= (double)1 - (double)3 * t * epsilon){
                dist[current_time-1][current_z][current_y][current_x]->next = dist[current_time][(current_z+1)% N][current_y][current_x];
                dist[current_time][(current_z+1)% N][current_y][current_x]->prev_judge = dist[current_time-1][current_z][current_y][current_x];
            }
            else if (randomnum <= (double)1 - (double)2 * t * epsilon){
                dist[current_time-1][current_z][current_y][current_x]->next = dist[current_time][current_z][current_y][(current_x-1+N)% N];
                dist[current_time][current_z][current_y][(current_x-1+N)% N]->prev_judge = dist[current_time-1][current_z][current_y][current_x];
            }
            else if (randomnum <= (double)1 - (double)t * epsilon){
                dist[current_time-1][current_z][current_y][current_x]->next = dist[current_time][current_z][(current_y-1+N)% N][current_x];
                dist[current_time][current_z][(current_y-1+N)% N][current_x]->prev_judge = dist[current_time-1][current_z][current_y][current_x];
            }
            else if (randomnum <= (double)1){
                dist[current_time-1][current_z][current_y][current_x]->next = dist[current_time][(current_z-1+N)% N][current_y][current_x];
                dist[current_time][(current_z-1+N)% N][current_y][current_x]->prev_judge = dist[current_time-1][current_z][current_y][current_x];
            }
            // ... More conditions to be translated similarly
            int curr_x = dist[current_time - 1][current_z][current_y][current_x]->next->x;
            int curr_y = dist[current_time - 1][current_z][current_y][current_x]->next->y;
            int curr_z = dist[current_time - 1][current_z][current_y][current_x]->next->z;
            current_x = curr_x;
            current_y = curr_y;
            current_z = curr_z;
        } 
        else if (current_time == n_time_slice - 1) {
            current_time = 0;
            dist[n_time_slice - 1][current_z][current_y][current_x]->next = dist[0][current_z][current_y][current_x];
            dist[0][current_z][current_y][current_x]->prev_judge = dist[n_time_slice - 1][current_z][current_y][current_x];
        }
    } else if (direction == -1) {
        dist[current_time][current_z][current_y][current_x]->status = false;
        int curr_x = dist[current_time][current_z][current_y][current_x]->prev->x;
        int curr_y = dist[current_time][current_z][current_y][current_x]->prev->y;
        int curr_z = dist[current_time][current_z][current_y][current_x]->prev->z;
        dist[current_time][current_z][current_y][current_x]->prev->next = nullptr;
        dist[current_time][current_z][current_y][current_x]->prev = nullptr;
        current_x = curr_x;
        current_y = curr_y;
        current_z = curr_z;
        if (current_time > 0) {
            current_time -= 1;
        } else if (current_time == 0) {
            current_time = n_time_slice - 1;
        }
    }
    return;
    }

    void site_update() {
    if (direction == -1) {
        if ((double)rand() / (double) RAND_MAX < std::min(1.0, exp(-epsilon * chemical_potential))) {
            direction = -1;
        } else {
            direction = 1;
        }
    } else if (direction == 1) {
        if ((double)rand() / (double) RAND_MAX < std::min(1.0, exp(epsilon * chemical_potential))) {
            direction = 1;
        } else {
            direction = -1;
        }
    }
    return;
    }

    void propagator_new() {
    bool judgement = true;
    bool prop;
    int curr_x;
    int curr_y;
    int curr_z;
    while (judgement) {
        bond_update();
        if (direction == 1 && current_time == starting_time && current_x == starting_point_x && current_y == starting_point_y && current_z == starting_point_z) {
            judgement = false;
            dist[current_time][current_z][current_y][current_x]->status = true;
            dist[current_time][current_z][current_y][current_x]->prev = dist[current_time][current_z][current_y][current_x]->prev_judge;
            dist[current_time][current_z][current_y][current_x]->prev_judge = nullptr;
        }
        else if (direction == 1){
            prop = dist[current_time][current_z][current_y][current_x]->status;
            if (prop == false){
                site_update();
                dist[current_time][current_z][current_y][current_x]->status = true;
                dist[current_time][current_z][current_y][current_x]->prev = dist[current_time][current_z][current_y][current_x]->prev_judge;
                dist[current_time][current_z][current_y][current_x]->prev_judge = nullptr;
            }
            else {
                curr_x = dist[current_time][current_z][current_y][current_x]->prev->x;
                curr_y = dist[current_time][current_z][current_y][current_x]->prev->y;
                curr_z = dist[current_time][current_z][current_y][current_x]->prev->z;
                dist[current_time][current_z][current_y][current_x]->prev = dist[current_time][current_z][current_y][current_x]->prev_judge;
                dist[current_time][current_z][current_y][current_x]->prev_judge = nullptr;
                current_x = curr_x;
                current_y = curr_y;
                current_z = curr_z;
                current_time = (current_time - 1 + n_time_slice) % n_time_slice;
                direction = -1;
                site_update();
                if ((current_x == starting_point_x) && (current_y == starting_point_y) && (current_z == starting_point_z) && (current_time == starting_time) && (direction == -1))
                {
                        judgement = false;
                        dist[current_time][current_z][current_y][current_x]->status = false;
                        dist[current_time][current_z][current_y][current_x]->next = nullptr;
                }
            }
        }
        else if (direction == -1)
        {
            site_update();
            if ((direction == -1) && (current_time == starting_time) && (current_x == starting_point_x) && (current_y == starting_point_y) && (current_z == starting_point_z))
                {
                    judgement = false;
                    dist[current_time][current_z][current_y][current_x]->status = false;
                }
        }
    }
    return;
    }
    int particle_sampler(){
        int ncount = 0;
        for (int i = 0; i < N; i++){
            for (int j = 0; j < N; j++){
                for (int k = 0; k < N; k++){
                    if (dist[0][i][j][k]->status){
                        ncount = ncount+1;
                    }
                }
            }
        }
    return ncount;
    }
};

int main() {
    srand(time(NULL));
    int Ntimes;
    std::cin >> Ntimes;
    double aveN[Ntimes];
    double aveE[Ntimes];
    double chemical_potential[1] = {1.0};
    for (int l = 0; l < 1; l++){
    WormPropagator worm((int)24, (int) 111, chemical_potential[l], (double)0.01, (double)1.0); // N, n_time_slice, chemical_potential, epsilon, t

    for (int i = 0; i < Ntimes; ++i) {
        double aveone[5000];
        double aveEone[5000];

        for (int j = 0; j < 5000; ++j) {
            worm.initialization();
            worm.propagator_new();
            // Call other member functions accordingly...
            // Uncomment to print particle_sampler() output
            if (j % 500 == 0) {
                std::cout << worm.particle_sampler() << " " << (j / 5000.0 / Ntimes + i / (double)Ntimes) << std::endl;
            }
            aveone[j] = worm.particle_sampler();
        }

        double sum = 0.0;
        for (int j = 0; j<5000; j++) {
            sum += aveone[j];
        }
        aveN[i] = sum / (double) 5000;
    }

    // Print aveN values
    double sum1 = 0.0;
    for (int i = 1; i < Ntimes; i++) {
        std::cout << aveN[i] << " ";
        sum1 += aveN[i];
    }
    double l1 = sum1/((double)Ntimes-(double)1);
    std::cout << l1 << std::endl;

    double sum2 = 0.0;
    for (int i = 1; i < Ntimes; i++){
        sum2 += (aveN[i]-l1)*(aveN[i]-l1);
    }
    sum2 = sqrt(sum2/(Ntimes-1));
    std::cout << sum2 << std::endl;
    }
    return 0;
}
