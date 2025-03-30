#include "../include/Multigrid.h"
#include <fstream>
#include <string>
#include <Strategy.h>
#include "../include/exprtk.hpp"
#include "json.hpp"

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
            //
        } else {
            //
        }
    } catch (const std::exception& e) {
        std::cerr << "錯誤: " << e.what() << std::endl;
    }
    system("python3 ../src/plot.py");
    return 0;
}

// --------------------------
// 5. 用户输入驱动示例
// --------------------------
int main() {
    // 假设用户输入如下参数（可从配置文件或命令行读取）
    std::string restriction_strategy = "half_weighting";
    std::string interpolation_strategy = "quadratic";
    std::string cycle_strategy = "f_cycle";

    // 创建策略对象
    auto restriction = StrategyFactory::createRestriction(restriction_strategy);
    auto interpolation = StrategyFactory::createInterpolation(interpolation_strategy);
    auto cycle = StrategyFactory::createCycle(cycle_strategy);

    // 配置求解器
    Solver solver;
    solver.setRestriction(std::move(restriction));
    solver.setInterpolation(std::move(interpolation));
    solver.setCycle(std::move(cycle));

    // 执行求解
    solver.solve();

    return 0;
}