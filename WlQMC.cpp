#include "WlQMC.h"



Worldline_Quantum_Monte_Carlo::Worldline_Quantum_Monte_Carlo(int m, int N, int Sz, double J, double Jz, double beta, int ind_i, int ind_j){
    if (N % 2){
        throw std::runtime_error("Dimensions of the checker board must be divisible by 2");
    }

    std::vector<int> Wl_start = create_vector(Sz,N);
    n_lines = Wl_start.size();

    if (n_lines >= N){
        throw std::runtime_error("Too many worldlines");
    }
    if (ind_i >= N || ind_j >= N){
        throw std::runtime_error("Indices for Spin-Spin-Correlation out of range");
    }
    srand(time(0));

    this -> N = N;
    this -> m = m;
    this -> J = J;
    this -> Jz = Jz;
    this -> beta = beta;
    this -> ind_i = ind_i;
    this -> ind_j = ind_j;
    m2 = 2*m;

    Wl.reserve(m2);
    Wl_zeros.reserve(m2*n_lines);
    Wl.resize(m2,Wl_start);

    for (int i = 0; i < m2*n_lines; i++){
        Wl_zeros.emplace_back(std::array<int,2>{i % m2, i / m2});
    }


    int ns_helper = 0;

    for (int i = 0; i < n_lines; i++){
        if (Wl_start[i] != Wl_start[cycle(i+1,n_lines)]){
            ns_helper++;
        }
    }

    ns = ns_helper * m2;
    nd = 0;
    nn = N * m - ns;

    double dtau = beta / static_cast<double>(m);
    double deltJ2 = dtau * J / 2;
    //std::cout << "deltJ2: " << deltJ2 << "\n";
    double g_ex = std::exp(dtau * Jz / 4);

        
    gs = g_ex * std::cosh(deltJ2);
    gd = -g_ex * std::sinh(deltJ2);
    gn = 1 / g_ex;
    fs = 1 / (2 * static_cast<double>(m)) * (J * gd / gs - Jz / 2);
    fd = 1 / (2 * static_cast<double>(m)) * (J * gs / gd - Jz / 2 );
    fn = Jz / (4 * static_cast<double>(m));
    std::cout << "fs = " << fs << " , fd = " << fd << ", fn = " << fn << "\n";
    std::cout << "gs = " << gs << " , gd = " << gd << ", gn = " << gn << "\n \n";

    for (int s = 0; s != 5; s++){
        for (int d = 0; d != 3; d++){
            Probability_Table[s][d] = std::pow(gs, 2*s-4) * std::pow(gd, 2*d-2) * std::pow(gn, 6-2*s-2*d);
            std::cout << Probability_Table[s][d] << " ";
        }
        std::cout << "\n";
    }

}

std::vector<int> Worldline_Quantum_Monte_Carlo::create_vector(int Sz, int N){
    std::vector<int> Wls;
    Wls.reserve(N/2 + Sz);
    Wls.resize(N/2);

    for (int i = 0; i < N/2; i++){
        Wls[i] = 2*i;
    }
    for (int i = 0; i < Sz; i++){
        Wls.insert(Wls.begin() + 2*i+1, 2*i+1);
    }
    for (int i = 0; i > Sz; i--){
        Wls.pop_back();
    }
    return Wls;

}

    /*
    int get_position(int tau, int line){
        if(line >= n_lines || tau >= m){
            throw std::runtime_error("Out of bounds");
        }

        int pos = Wl_start[line];

        for (int i = 0; i < tau; i++){
            pos += Wl[i][line];
        }
        return pos;    
    }*/

int Worldline_Quantum_Monte_Carlo::cycle(int x, int x_max){
    if (x >= x_max){
        x -= x_max;
    }
    else if (x < 0){
        x += x_max;
    }
    return x;
}

int Worldline_Quantum_Monte_Carlo::direction_correct(int dir){
    if (dir > 1){
        dir -= N;
    } else if (dir < -1){
        dir += N;
    }
    return dir;
}

int Worldline_Quantum_Monte_Carlo::get_direction(int tau, int line){
    int dir = direction_correct(Wl[cycle(tau + 1,m2)][line] - Wl[tau][line]);
    return dir;
}


std::array<int,3> Worldline_Quantum_Monte_Carlo::get_direction(std::array<int,3> tau, int line){
    std::array<int,3> dir;

    for (int i = 0; i < 3; i++){
        dir[i] = direction_correct(Wl[cycle(tau[i] + 1,m2)][line] - Wl[cycle(tau[i],m2)][line]);
    }

    return dir;
}

void Worldline_Quantum_Monte_Carlo::total_weight(int& ns, int& nd){
    int nsp = 0;
    int nnp = 0;
    for (int i = 0; i < m2; i++){
        for (int j = 0; j < n_lines; j++){
            if (Wl[cycle(i+1,m2)][j] - Wl[i][j] == 0){
                nsp++;
                int plaq_dir = (((Wl[i][j] + i) % 2) *2 -1) * (-1);
                int plaq_line = cycle(Wl[i][j] + plaq_dir, N);
                int plaq_next = Wl[i][cycle(j+plaq_dir,n_lines)];
                if (plaq_line == plaq_next){
                    nnp++;
                }
            }

        }
        nd = m2*n_lines - nsp;
        ns = nsp - nnp;
    }

}


void Worldline_Quantum_Monte_Carlo::Generate_Proposal(){
    Prop.rand_int = rand() % Wl_zeros.size();
    Prop.rand_zero = Wl_zeros[Prop.rand_int];

    int dir_data = 2 * Wl[Prop.rand_zero[0]][Prop.rand_zero[1]] + 2 * Prop.rand_zero[0] + 1;
    Prop.direction = dir_data % 4 - 2;
}

void Worldline_Quantum_Monte_Carlo::Generate_Proposal(int seed){
    Prop.rand_int = seed;
    Prop.rand_zero = Wl_zeros[Prop.rand_int];

    int dir_data = 2 * Wl[Prop.rand_zero[0]][Prop.rand_zero[1]] + 2 * Prop.rand_zero[0] + 1;
    Prop.direction = dir_data % 4 - 2;
}

bool Worldline_Quantum_Monte_Carlo::Check_Collision(){
        
    int new_pos = cycle(Wl[Prop.rand_zero[0]][Prop.rand_zero[1]] + Prop.direction, N);
    int new_pos_1 = cycle(Wl[cycle(Prop.rand_zero[0] + 1, m2)][Prop.rand_zero[1]] + Prop.direction, N);

    int next_pos = Wl[Prop.rand_zero[0]][cycle(Prop.rand_zero[1] + Prop.direction, n_lines)];
    int next_pos_1 = Wl[cycle(Prop.rand_zero[0] + 1, m2)][cycle(Prop.rand_zero[1] + Prop.direction, n_lines)];

    return (new_pos == next_pos || new_pos_1 == next_pos_1);
}




double Worldline_Quantum_Monte_Carlo::Get_Weight(int& delta_ns, int& delta_nd, int& delta_nn, std::array<bool,2>& curved){

    int tau = Prop.rand_zero[0];
    std::array<int,3> tau_arr {tau - 1, tau, tau + 1};
    std::array<int,3> directions = get_direction(tau_arr, Prop.rand_zero[1]);
    int W_i = Wl[tau][Prop.rand_zero[1]];

    curved[0] = !(static_cast<bool>(directions[0]));
    curved[1] = !(static_cast<bool>(directions[2]));

    delta_nn = 0;
    delta_nn -= 2 * static_cast<int>(cycle(Wl[tau][Prop.rand_zero[1]] - Prop.direction,N) == Wl[tau][cycle(Prop.rand_zero[1] - Prop.direction,n_lines)]) - 1;
    delta_nn += 2 * static_cast<int>(cycle(Wl[tau][Prop.rand_zero[1]] + 2*Prop.direction,N) == Wl[tau][cycle(Prop.rand_zero[1] + Prop.direction,n_lines)]) - 1;

    //std::cout << directions[0] << " " << directions[2] << "\n";
    //std::cout << curved[0] << " " << curved[1] << "\n";
    if (curved[0]){
        if (curved[1]){
            delta_ns = -2 - delta_nn;
            delta_nd = 2;
        }
        else{
            delta_ns = -delta_nn;
            delta_nd = 0;
        }
    }
    else{
        if (curved[1]){
            delta_ns = -delta_nn;
            delta_nd = 0;
        }
        else{
            delta_ns = 2 - delta_nn;
            delta_nd = -2;
        }
    }

    //std::cout << "Checking Values: " << delta_ns << " " << delta_nd << " " << delta_nn << "\n";
    //std::cout << Probability_Table[(delta_ns+4)/2][(delta_nd+2)/2] << "\n";

    return Probability_Table[(delta_ns+4)/2][(delta_nd+2)/2];
        
}

bool Worldline_Quantum_Monte_Carlo::Check_Rejection(double weight){
        double r = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
        return (r > weight);
}

double Worldline_Quantum_Monte_Carlo::energy(){
    double s = fs * static_cast<double>(ns);
    double d = fd * static_cast<double>(nd);
    double n = fn * static_cast<double>(nn);
    return(n + s + d);
}

double Worldline_Quantum_Monte_Carlo::spin_spin_corr(){
    int sum = 0;
    int spin_i = -1;
    int spin_j = -1;
    for (auto& i : Wl){
        for (auto& j : i){
            if (j == ind_i){
                spin_i = 1;
            }
            else if(j == ind_j){
                spin_j = 1;
            }
        }
        sum += spin_i * spin_j;
        spin_i = -1;
        spin_j = -1;
    }
    double S = static_cast<double>(sum) / (8 * static_cast<double>(m));
    return S;    
}


std::array<std::vector<double>*,2> Worldline_Quantum_Monte_Carlo::run(int measurements, const std::string& foldername){
    time_t t1 = time(0);
    Energy_Vec = new std::vector<double>;
    SSCorr_Vec = new std::vector<double>;
    Energy_Vec->reserve(measurements+1);
    Energy_Vec->emplace_back(energy());
    SSCorr_Vec->reserve(measurements+1);
    SSCorr_Vec->emplace_back(spin_spin_corr());
    long long cout_count = 1;
    char pm = static_cast<char>(241);
    std::vector<std::array<double,2>> E_data;
    std::array<double,2> E_proto_data;
    std::vector<std::array<double,2>> S_data;
    std::array<double,2> S_proto_data;
    std::filesystem::create_directory(foldername);
    int i_old = 0;


    for (int i = 1; i <= measurements; i++){

        //std::cout << " Probability: " << weight << "\n";

        sweep();

        long long i5 = static_cast<long long>(i) * 100;
        long long c = cout_count * static_cast<long long>(measurements);
        if (i5 >= c){
            E_proto_data[0] = mean(Energy_Vec, false, i_old);
            E_proto_data[1] = std::sqrt(var(Energy_Vec, E_proto_data[0], i_old) / (i-i_old));
            S_proto_data[0] = mean(SSCorr_Vec, false, i_old);
            S_proto_data[1] = std::sqrt(var(SSCorr_Vec, S_proto_data[0], i_old) / (i-i_old));
            i_old = i;
            std::cout << "\n" << "Took Measurement " << i << ": ";
            std::cout << "ns = " << ns << ", nd = " << nd << ", nn = " << nn << "\n" << "Energy: " << E_proto_data[0] << " " << pm << " " << E_proto_data[1] << "\n";
            std::cout << "Spin-Spin-Correlation: " << S_proto_data[0] << " " << pm << " " << S_proto_data[1] << "\n";
            cout_count++;
            E_data.emplace_back(E_proto_data);
            S_data.emplace_back(S_proto_data);
        }
    }


        /*std::cout << "Delta ns: " << delta_ns << ", Delta nd: " << delta_nd << "\n";
        std::cout << "Ns: " << Worldline_Quantum_Monte_Carlo::ns << ", Nd: " << Worldline_Quantum_Monte_Carlo::nd << "\n";
        std::cout << "Curved Bool Array: " << curved[0] << " " << curved[1] << "\n";
        std::cout << P->rand_zero[0] << " " << P->rand_zero[1] << "\n";
        std::cout << "Collisions: " << no << " " << nn << "\n";
        std::cout << "Case Number " << ((int) curved[0]) * 8 + ((int) curved[1]) * 4 + no * 2 +nn << "\n";*/
            

/*
            WlQMC->total_weight(nst,ndt);
            std::cout << "Ns total: " << nst << ", Nd total: " << ndt << "\n";
            for (int i = 0; i < Worldline_Quantum_Monte_Carlo::m; i++){
                for (int j = 0; j < Worldline_Quantum_Monte_Carlo::n_lines; j++){
                    std::cout << WlQMC->Wl[i][j] << " ";
                }
                std::cout << "\n";
            }*/
            
            //system("Pause");
            //P->generate();


    int td = static_cast<int>(difftime(time(0),t1));
    std::cout << "\n" << "Calculation time: " << td / 3600 << " h ";
    std::cout << (td % 3600) / 60 << " min ";
    std::cout << (td % 60) << " sec \n \n";

    std::string E_data_name = foldername + "/E_data";
    saveData(E_data,E_data_name);
    std::string S_data_name = foldername + "/S_data";
    saveData(S_data,S_data_name);


    std::array<std::vector<double>*,2> result {Energy_Vec, SSCorr_Vec};
    return result;
}

void Worldline_Quantum_Monte_Carlo::sweep(){
    long long i = 0;
    long long loops = static_cast<long long>(m2) * static_cast<long long>(N);
    int delta_ns = 0;
    int delta_nd = 0;
    int delta_nn = 0;
    std::array<bool,2> curved;

    while (i < loops)
    {
        Generate_Proposal();
        if (Check_Collision()){
            continue;
        }
        i++;
        double weight = Get_Weight(delta_ns, delta_nd, delta_nn, curved);
        if (Check_Rejection(weight)){
            continue;
        }
        ns += delta_ns;
        nd += delta_nd;
        nn += delta_nn;
        zero_change(Prop.rand_int, curved);
        Wl[Prop.rand_zero[0]][Prop.rand_zero[1]] = cycle(Wl[Prop.rand_zero[0]][Prop.rand_zero[1]] + Prop.direction, N);
        Wl[cycle(Prop.rand_zero[0] + 1, m2)][Prop.rand_zero[1]] = cycle(Wl[cycle(Prop.rand_zero[0] + 1, m2)][Prop.rand_zero[1]] + Prop.direction, N);
        /*
        std::cout << "Delta Ns: " << delta_ns << " Delta Nd: " << delta_nd << "\n";
        std::cout << "Coordinates: " << Prop.rand_zero[0] << " " << Prop.rand_zero[1] << "\n";
        for (int i = 0; i < m2; i++){
            for (auto& w : Wl[i]){
                std::cout << w << " ";
            }
            std::cout << "\n";
        }*/

    }

    Energy_Vec->emplace_back(energy());
    SSCorr_Vec->emplace_back(spin_spin_corr());
    
}

template<size_t n>
void Worldline_Quantum_Monte_Carlo::saveData(std::vector<std::array<double,n>>& v, const std::string& filename){
    std::ofstream file(filename, std::ios::binary);
    int vsize = v[0].size() * 8;

    if (!file) {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }

    for (auto& i : v){
        file.write(reinterpret_cast<const char*>(i.data()), vsize);
    }
}


void Worldline_Quantum_Monte_Carlo::zero_change(int r, std::array<bool,2> curved){
    int ind;
    std::array<int,2> zero_coord = Wl_zeros[r];
    zero_coord[0] = cycle(zero_coord[0] + 1, m2);
    //std::cout << zero_coord[0] << zero_coord[1] << "\n";

    if (curved[1]){
        ind = zero_search(r, zero_coord, 1);
        Wl_zeros.erase(Wl_zeros.begin() + ind);  
        if (ind != r+1){
            r -= 1;
        }
    } else {
        Wl_zeros.insert(Wl_zeros.begin() + r + 1, zero_coord);    
    }
    zero_coord[0] = cycle(zero_coord[0] - 2, m2);
    //std::cout << zero_coord[0] << zero_coord[1] << "\n";


    if (curved[0]){
        ind = zero_search(r, zero_coord, -1);
        //std::cout << ind << "\n";
        Wl_zeros.erase(Wl_zeros.begin() + ind);
    } else{
        Wl_zeros.insert(Wl_zeros.begin() + r, zero_coord);
    }

}
    
int Worldline_Quantum_Monte_Carlo::zero_search(int r, std::array<int,2> obj, int dir){
        
    if (Wl_zeros[cycle(r+dir,Wl_zeros.size())] == obj){
        return cycle(r+dir,Wl_zeros.size());
    }
    //std::cout << "PosInd: " << r+2*dir << "\n";
    int ind;

    for (int i = m2; i > 0; i--){
        ind = cycle(r - i * dir, Wl_zeros.size());
        //std::cout << "PosInd: " << ind << "\n";
        if (Wl_zeros[ind] == obj){
            return ind;
        }
    }

    throw std::runtime_error("Something went wrong in the Storage of Zero Direction Points");
        
}

template<typename T>
double Worldline_Quantum_Monte_Carlo::mean(std::vector<T>* v, bool squared){
    T vsum = 0;
    if (squared){
        for (auto& i : (*v)){
            vsum += i * i;
        }
    } else{
        for (auto& i : (*v)){
            vsum += i;
        }
    }
    double vmean = static_cast<double>(vsum) / static_cast<double>(v->size());
    return vmean;
}

template<typename T>
double Worldline_Quantum_Monte_Carlo::mean(std::vector<T>* v, bool squared, int start){
    T vsum = 0;
    if (squared){
        for (typename std::vector<T>::iterator i = v->begin() + start; i != v->end(); i++){
            vsum += (*i) * (*i);
        }
    } else{
        for (typename std::vector<T>::iterator i = v->begin() + start; i != v->end(); i++){
            vsum += *i;
        }
    }
    double vmean = static_cast<double>(vsum) / static_cast<double>(v->size() - start);
    return vmean;
}

template<typename T>
double Worldline_Quantum_Monte_Carlo::var(std::vector<T>* v){
    double vmean = mean(v,false);
    double vmean2 = mean(v,true);
    double dev = vmean2 - vmean * vmean;
    return dev;
}

template<typename T>
double Worldline_Quantum_Monte_Carlo::var(std::vector<T>* v, double vmean){
    double vmean2 = mean(v,true);
    double dev = vmean2 - vmean * vmean;
    return dev;
}

template<typename T>
double Worldline_Quantum_Monte_Carlo::var(std::vector<T>* v, double vmean, int start){
    double vmean2 = mean(v,true,start);
    double dev = vmean2 - vmean * vmean;
    return dev;
}

template<typename T>
std::vector<double>* Worldline_Quantum_Monte_Carlo::autocorrelation(std::vector<T>* v, double vmean, double vvar, int tmax, double* tau) {
    double vmean2 = vmean * vmean;
    int vsize = v->size();
    std::vector<double>* aut_vec = new std::vector<double>;
    aut_vec->reserve(tmax);
    double aut;

    if (tau){
        *tau = 0.5;
    }

    for (int t = 0; t < tmax; t++) {
        T aut_sum = 0;
        for (int i = 0; i < (vsize - t); i++) {
            aut_sum += (*v)[i] * (*v)[i + t];
        }
        aut = (static_cast<double>(aut_sum) / static_cast<double>(vsize - t) - vmean2) / vvar;
        aut_vec->emplace_back(aut);
        if (tau) {
            *tau += aut * (1 - t / vsize);
            }
    }
    return aut_vec;
}


template<typename T>
std::vector<double>* Worldline_Quantum_Monte_Carlo::autocorrelation(std::vector<T>* v, double vmean, double vvar, int tmax, double* tau, int start) {
    double vmean2 = vmean * vmean;
    int vsize = v->size() - start;
    std::vector<double>* aut_vec = new std::vector<double>;
    aut_vec->reserve(tmax);
    double aut;

    if (tau){
        *tau = 0.5;
    }

    for (int t = 0; t < tmax; t++) {
        T aut_sum = 0;
        for (int i = 0; i < (vsize - t); i++) {
            aut_sum += (*v)[i + start] * (*v)[i + start + t];
        }
        aut = (static_cast<double>(aut_sum) / static_cast<double>(vsize - t) - vmean2) / vvar;
        aut_vec->emplace_back(aut);
        if (tau) {
            *tau += aut * (1 - t / vsize);
            }
    }
    return aut_vec;
}

















