#include "../include/EquationSolverFactory.h"
#include "json.hpp"

void test_convergence(EquationSolver* solver,
    const HardcodedFunction_7& Func,
    int condition_type,
    const std::vector<int>& Ns,
    double T) {
    std::vector<double> hs;
    std::vector<double> errors;

    for (int N : Ns) {
    double h = T / N;
    hs.push_back(h);

    std::string dummy_file = "temp.txt"; // 或者用 /dev/null
    solver->solve(Func, T, N, condition_type, dummy_file, false);

    std::vector<double> u_start = solver->get_u_start();  //初始解
    std::vector<double> u_end = solver->get_u_end();    //最终解

    double err = solver->error(u_start, u_end);
    errors.push_back(err);
    std::cout << "N = " << N << ", h = " << h << ", Error = " << err << std::endl;
    }

    // 计算收敛阶
    for (int i = 1; i < errors.size(); ++i) {
    double p = log(errors[i - 1] / errors[i]) / log(hs[i - 1] / hs[i]);
    std::cout << "Convergence rate between N = " << Ns[i-1]
    << " and " << Ns[i] << " is approximately p = " << p << std::endl;
    }
}

int main() {
    std::ifstream ifs("../test/main.json");
    nlohmann::json data = nlohmann::json::parse(ifs);
    std::string inputDir = "../test/";

    std::string FuncName = data["FuncName"];
    std::string method = data["method"];
    int condition_type = data["condition_type"];
    int step = data["step"];
    std::string output_file = data["output_file"];
    bool print_accurate = data["print_accurate"];
    bool is_test = data["is_test"];
    
    std::ifstream ifsFunc(inputDir + FuncName + ".txt");
    std::string ofs(inputDir + output_file + ".txt");

    double t1 = 17.06521656015796;
    double t2 = 19.140540691377;
    double T;
    if(condition_type == 1)
        T = t1;
    else if(condition_type == 2)
        T = t2;
    else;
    
    RegisterAllEquationSolvers();
    auto& fac = EquationSolverFactory::getInstance();
    //
    std::unique_ptr<EquationSolver> solver = fac.createEquationSolver(method);
    // DynamicFunction_7 fuc(ifsFunc);
    HardcodedFunction_7 func;
    solver->solve(func, T, step, condition_type, ofs, is_test);
    if(print_accurate){
        std::vector<int> Ns = {10000, 20000, 40000, 80000};
        test_convergence(solver.get(), func, condition_type, Ns, T);
    }

    system("python3 ../src/plot.py ");
    
    return 0;
}