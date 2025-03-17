#include "../include/FD_Method.h"
#include "../include/I_FD_Method.h"
#include <fstream>
#include <string>
#include "../include/exprtk.hpp"
#include "json.hpp"

int main() {
    std::ifstream ifs("../test/test.json");
    nlohmann::json data = nlohmann::json::parse(ifs);

    std::string funcName = data["funcName"].get<std::string>();
    std::string boundaryFuncNamex = data["boundaryFuncNamex"].get<std::string>();
    std::string boundaryFuncNamey = data["boundaryFuncNamey"].get<std::string>();
    int pixelNum = data["pixelNum"].get<int>();
    std::string boundaryCondition = data["boundaryCondition"].get<std::string>();
    std::string boundaryStructure = data["boundaryStructure"].get<std::string>();
    std::vector<double> center = data["center"].get<std::vector<double>>();
    double radius = data["radius"].get<double>();
    std::string center_Func = data["center_Func"].get<std::string>();
    std::string exact_solution = data["exactSolution"].get<std::string>();
    
    try {
        DynamicFunction func(funcName);
        DynamicFunction boundaryFuncx(boundaryFuncNamex);
        DynamicFunction boundaryFuncy(boundaryFuncNamey);
        DynamicFunction exact(exact_solution);
        if(boundaryStructure == "Regular") {
            FD_Method method(pixelNum, func, boundaryFuncx, boundaryFuncy, boundaryCondition);
            method.construct_equation();
            method.solver();
            method.error_analysis(exact);
        } else {
            DynamicFunction centerFuc(center_Func);
            I_FD_Method method(pixelNum, func, boundaryFuncx, boundaryFuncy, centerFuc, boundaryCondition, center, radius);
            method.I_solver();
            method.error_analysis(exact);
        }
    } catch (const std::exception& e) {
        std::cerr << "錯誤: " << e.what() << std::endl;
    }
    system("python3 ../src/plot.py");
    return 0;
}