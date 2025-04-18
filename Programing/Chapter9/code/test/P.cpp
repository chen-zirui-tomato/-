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

    int Dim = data["Dim"].get<int>();
    std::string InteratorFuncName = data["InteratorFuncName"].get<std::string>();
    std::string RestrictionFuncName = data["RestrictionFuncName"].get<std::string>();
    std::string InterpolationFuncName = data["InterpolationFuncName"].get<std::string>();
    double solverTorlerance = data["solverTorlerance"].get<double>();
    int nu1 = data["nu1"].get<int>();
    int nu2 = data["nu2"].get<int>();
    int maxIter = data["maxIter"].get<int>();
    int pixelNum = data["pixelNum"].get<int>();
    std::string boundaryCondition = data["boundaryCondition"].get<std::string>();
    std::string exact_solution = data["exactSolution"].get<std::string>();

    std::string funcName;
    std::string boundaryFuncName;
    std::string boundaryFuncNamex;
    std::string boundaryFuncNamey;
    if(Dim == 2){ 
        funcName= data["funcName2"].get<std::string>();
        boundaryFuncName  = data["boundaryFuncName2"].get<std::string>();
        boundaryFuncNamex = data["boundaryFuncNamex2"].get<std::string>();
        boundaryFuncNamey = data["boundaryFuncNamey2"].get<std::string>();
        try{
            DynamicFunction func(funcName);
            DynamicFunction boundaryFunc(boundaryFuncName);
            DynamicFunction boundaryFuncx(boundaryFuncNamex);
            DynamicFunction boundaryFuncy(boundaryFuncNamey);
            DynamicFunction exact(exact_solution);
            Multigrid<2> mg(func, boundaryFunc, boundaryFuncx, boundaryFuncy, boundaryCondition, exact);
            mg.Solver(pixelNum, RestrictionFuncName, InterpolationFuncName, 
                    InteratorFuncName, nu1, nu2, solverTorlerance, maxIter);
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
        }
    }
    else if(Dim == 1){ 
        funcName = data["funcName1"].get<std::string>();
        boundaryFuncName = data["boundaryFuncName1"].get<std::string>();
        try{
            DynamicFunction func(funcName);
            DynamicFunction boundaryFunc(boundaryFuncName);
            DynamicFunction exact(exact_solution);
            Multigrid<1> mg(func, boundaryFunc, boundaryCondition, exact);
            mg.Solver(pixelNum, RestrictionFuncName, InterpolationFuncName, 
                    InteratorFuncName, nu1, nu2, solverTorlerance, maxIter);
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
        }
    }
    return 0;


    // // 1D 求解器
    // Multigrid<1> mg1d(f, g, Dirichlet);
    // mg1d.Solve(5, full_weighting, linear, V_Cycle, 2, 2, 1e-6);

    // // 2D 求解器
    // Multigrid<2> mg2d(f, g, Neumann);
    // mg2d.Solve(5, injection, quadratic, FMG_Cycle, 3, 3, 1e-6);

//     try {
//         DynamicFunction func(funcName);
//         DynamicFunction boundaryFunc(boundaryFuncName);
//         DynamicFunction boundaryFuncx(boundaryFuncNamex);
//         DynamicFunction boundaryFuncy(boundaryFuncNamey);
//         DynamicFunction exact(exact_solution);
//             Multigrid<1> mg(func, boundaryFunc, boundaryCondition, exact);
//             mg.Solver(pixelNum, RestrictionFuncName, InterpolationFuncName, 
//                       InteratorFuncName, nu1, nu2, solverTorlerance, maxIter);
//     } catch (const std::exception& e) {
//         std::cerr << "Error: " << e.what() << std::endl;
//     }
//     // system("python3 ../src/plot.py");
//     return 0;
}