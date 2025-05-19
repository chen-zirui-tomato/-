#include "../include/EquationSolverFactory.h"
#include "json.hpp"
#include <chrono> 

int main() {
    using namespace std::chrono;
    std::ifstream ifs("../test/main.json");
    nlohmann::json data = nlohmann::json::parse(ifs);
    std::string inputDir = "../test/";

    std::string method = data["method"];
    double r = data["r"];
    bool is_exact = data["is_exact"];
    double time = data["time"]; 
    std::string output_file = data["output_file"];

    // std::string ofs(inputDir + output_file + ".txt");

    RegisterAllEquationSolvers();
    auto& fac = EquationSolverFactory::getInstance();
    std::unique_ptr<Solver> solver = fac.createEquationSolver(method);

    auto total_start = high_resolution_clock::now();

    if(is_exact){
        solver->compute_true_solution(time/400);
    }
    else
        solver->solve(r);
    solver->plot(inputDir + output_file + ".txt");

    auto total_end = high_resolution_clock::now();
    auto total_duration = duration_cast<milliseconds>(total_end - total_start);
    std::cout << "Total CPU Time : " << total_duration.count() / 1000.0 << " seconds." << std::endl;

    system("python3 ../src/plot.py ");
    return 0;
}