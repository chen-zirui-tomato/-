#include "../include/Multigrid.h"
#include "json.hpp"

#include <fstream>
#include <string>
#include <iostream>

int main() {
    std::ifstream ifs("../test/test.json");
    // std::ifstream ifs("../test/test2.json");
    // std::ifstream ifs("../test/test3.json");
    nlohmann::json data = nlohmann::json::parse(ifs);

    std::string funcName = data["funcName"].get<std::string>();
    std::string boundaryFuncName  = data["boundaryFuncName"].get<std::string>();
    std::string boundaryFuncNamex = data["boundaryFuncNamex"].get<std::string>();
    std::string boundaryFuncNamey = data["boundaryFuncNamey"].get<std::string>();
    int pixelNum = data["pixelNum"].get<int>();
    std::string boundaryCondition = data["boundaryCondition"].get<std::string>();
    std::string boundaryStructure = data["boundaryStructure"].get<std::string>();
    std::vector<double> center = data["center"].get<std::vector<double>>();
    double radius = data["radius"].get<double>();
    std::string center_Func = data["center_Func"].get<std::string>();
    std::string exact_solution = data["exactSolution"].get<std::string>();

    std::string InteratorFuncName = data["InteratorFuncName"].get<std::string>();
    std::string RestrictionFuncName = data["RestrictionFuncName"].get<std::string>();
    std::string InterpolationFuncName = data["InterpolationFuncName"].get<std::string>();
    int solverTorlerance = data["solverTorlerance"].get<int>();
    int nu1 = data["nu1"].get<int>();
    int nu2 = data["nu2"].get<int>();

    // // 1D 求解器
    // Multigrid<1> mg1d(f, g, Dirichlet);
    // mg1d.Solve(5, full_weighting, linear, V_Cycle, 2, 2, 1e-6);

    // // 2D 求解器
    // Multigrid<2> mg2d(f, g, Neumann);
    // mg2d.Solve(5, injection, quadratic, FMG_Cycle, 3, 3, 1e-6);

    try {
        DynamicFunction func(funcName);
        DynamicFunction boundaryFunc(boundaryFuncName);
        DynamicFunction boundaryFuncx(boundaryFuncNamex);
        DynamicFunction boundaryFuncy(boundaryFuncNamey);
        DynamicFunction exact(exact_solution);
        if(boundaryStructure == "Regular") {
            Multigrid<1> mg(func, boundaryFunc, boundaryCondition, exact);
            mg.Solver(pixelNum, RestrictionFuncName, InterpolationFuncName, 
                      InteratorFuncName, nu1, nu2, solverTorlerance);
            std::cerr << mg.errorAnalysis() << std::endl;
        } else {
            //
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    // system("python3 ../src/plot.py");
    return 0;
}