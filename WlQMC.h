#ifndef WLQMC_H
#define WLQMC_H

#include <iostream>
#include <vector>
#include <array>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <new>
#include <iterator>
#include <chrono>
#include <fstream>
#include <filesystem>
#include <sstream>
#include <string>


class Worldline_Quantum_Monte_Carlo {
    public:
    int m, N, m2;
    int n_lines;
    double J, Jz; 
    double beta;
    double fs, fd, fn, gs, gd, gn, mass;
    int ns, nd, nn;
    int ind_i, ind_j;
    std::array<std::array<double,3>,5> Probability_Table;

    std::vector<std::vector<int>> Wl;
    std::vector<std::array<int,2>> Wl_zeros;

    std::vector<int> create_vector(int Sz, int N);

    int cycle(int x, int x_max);
    int direction_correct(int dir);
    int get_direction(int tau, int line);
    std::array<int,3> get_direction(std::array<int,3> tau, int line);
    void total_weight(int& ns, int& nd);

    struct Proposal{
        int rand_int, direction;
        std::array<int,2> rand_zero;
    };

    void Generate_Proposal();
    void Generate_Proposal(int seed);

    bool Check_Collision();
    double Get_Weight(int& delta_ns, int& delta_nd, int& delta_nn, std::array<bool,2>& curved);
    bool Check_Rejection(double weight);

    Worldline_Quantum_Monte_Carlo(int m, int N, int Sz, double J, double Jz, double beta, int i, int j);
    Proposal Prop;

    std::array<std::vector<double>*,2> run(int measurements, const std::string& foldername);
    void sweep();
    double energy();
    double spin_spin_corr();
    template<size_t n>
    void saveData(std::vector<std::array<double,n>>& v, const std::string& filename);


    std::vector<double>* Energy_Vec;
    std::vector<double>* SSCorr_Vec;
    void zero_change(int r, std::array<bool,2> curved);
    int zero_search(int r, std::array<int,2> obj, int dir);

    template<typename T>
    double mean(std::vector<T>* v, bool squared);
    template<typename T>
    double mean(std::vector<T>* v, bool squared, int start);
    template<typename T>
    double var(std::vector<T>* v);
    template<typename T>
    double var(std::vector<T>* v, double vmean);
    template<typename T>
    double var(std::vector<T>* v, double vmean, int start);
    template<typename T>
    std::vector<double>* autocorrelation(std::vector<T>* v, double vmean, double vvar, int tmax, double* tau);
    template<typename T>
    std::vector<double>* autocorrelation(std::vector<T>* v, double vmean, double vvar, int tmax, double* tau, int start);

};

#endif
