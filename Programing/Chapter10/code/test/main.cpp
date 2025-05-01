#include "../include/EquationSolverFactory.h"
#include "json.hpp"

int main() {
    std::ifstream ifs("../test/test.json");
    nlohmann::json data = nlohmann::json::parse(ifs);

    int order = data["order"];
    std::string method = data["method"];
    std::string FuncName = data["FuncName"];
    double step = data["step"];


    RegisterAllEquationSolvers();
    auto& fac = EquationSolverFactory::getInstance();
    //创建四解AB求解器
    std::unique_ptr<EquationSolver> solver = fac.createEquationSolver("AB4");
    DynamicFunction fuc(FuncName);
    solver->solve(fuc, 1, step);

    return 0;
}