#include "WlQMC.h"
#include "WlQMC.cpp"




template<typename T>
void CheckerBoardInput(T& var, bool even = false, bool spin = false, int Ncont = 0){
    do{
        std::cin >> var;
        std::cout << "\n";

        if (std::cin.fail()){
            std::cout << "Not a valid value \n";
            std::cin.clear();
            continue;
        }

        if (even && static_cast<int>(var) % 2 != 0){
            std::cout << "Number should be divisible by 2 \n";
            continue;
        }
        if (spin && std::abs(var) >= Ncont){
            std::cout << "Electron assembly with this total Spin not possible or pointless \n";
            continue;
        }
        break;
    }
    while (true);
}

void FolderInput(std::string& var){
    std::cin >> var;
    std::cout << "\n";

}





int main(){

    int m;
    int N;
    int Sz;
    int perc;
    int i, j;
    std::vector<int> vector {0,2,4,6,8};
    int J, Jz;
    double beta, mass;
    int measurements;
    std::string foldername;
    char pm = static_cast<char>(241);

    std::cout << "Type in value for m: ";
    CheckerBoardInput(m,true);
    std::cout << "Type in value for N: ";
    CheckerBoardInput(N,true);
    std::cout << "Type in value for total Spin: ";
    CheckerBoardInput(Sz,false,true,N);

    std::cout << "Type in value for J: ";
    CheckerBoardInput(J);
    std::cout << "Type in value for Jz: ";
    CheckerBoardInput(Jz);
    std::cout << "Type in value for beta: ";
    CheckerBoardInput(beta);
    std::cout << "Type in first index for spin spin correlation: ";
    CheckerBoardInput(i);
    std::cout << "Type in second index for spin spin correlation: ";
    CheckerBoardInput(j);
    std::cout << "Type in number of measurements: ";
    CheckerBoardInput(measurements);
    std::cout << "Name the folder for saved data: ";
    FolderInput(foldername);
    


    Worldline_Quantum_Monte_Carlo W(m,N,Sz,J,Jz,beta,i,j);
    std::cout << "Worldline Quantum Monte Carlo Engine Initialized. Press any key to run. \n \n";
    system("Pause");
    try{
        std::array<std::vector<double>*,2> ES = W.run(measurements,foldername);
        std::vector<double>* Energy = ES[0];
        std::vector<double>* SSCor = ES[1];
        std::string command = "python3 main.py " + foldername;
        std::system(command.c_str());
        std::cout << "Type in starting percentage of values to be analysed: ";
        CheckerBoardInput(perc);
        int start = ((measurements+1)*perc)/100;
        double E_m = W.mean(Energy, false, start);
        double E_v = W.var(Energy, E_m, start);
        double* aut_time_e = new double;
        std::vector<double>* aute = W.autocorrelation(Energy, E_m, E_v, 1000, aut_time_e, start);
        double S_m = W.mean(SSCor, false, start);
        double S_v = W.var(SSCor, S_m);
        double* aut_time_s = new double;
        std::vector<double>* auts = W.autocorrelation(SSCor, S_m, S_v, 1000, aut_time_s, start);
        double Epm = std::sqrt(2*(*aut_time_e)*E_v/(measurements-start));
        double Spm = std::sqrt(2*(*aut_time_s)*S_v/(measurements-start));
        std::cout << "Energy: " << E_m << " " << pm << " " << Epm << "\n";
        std::cout << "Spin-Spin-Correlation: " << S_m << " " << pm << " " << Spm << "\n";
    }
    catch(std::exception& e){
        std::cout << e.what() << "\n";
    }
    
    
    system("Pause");
    system("Pause");
    system("Pause");

}




