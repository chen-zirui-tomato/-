#pragma once
#include <iostream>
#include <string>
#include <map>
#include <memory>
#include "EquationSolver.h"

class EquationSolverFactory {
public:
    using CreateEquationSolverCallback = std::unique_ptr<EquationSolver> (*)(int);
private:
    using CallbackMap = std::map<std::string, CreateEquationSolverCallback>;
public:
    void registerEquationSolver(const std::string& ID, CreateEquationSolverCallback createFn) {
        callbacks_[ID] = createFn;
    }
    std::unique_ptr<EquationSolver> createEquationSolver(const std::string &ID) {
        std::string method;
        int order = 0;
        size_t pos = 0;
    
        // 提取方法名（非数字部分）
        while (pos < ID.size() && !std::isdigit(ID[pos])) {
            method += ID[pos];
            pos++;
        }
    
        // 提取阶数（剩余数字部分）
        try {
            if (pos < ID.size()) {
                order = std::stoi(ID.substr(pos));
            } else {
                std::cerr << "Error: Missing order in solver ID '" << ID << "'." << std::endl;
                return nullptr;
            }
        } catch (const std::invalid_argument& e) {
            std::cerr << "Error: Invalid order format in solver ID '" << ID 
                     << "'. Expected numeric suffix." << std::endl;
            return nullptr;
        } catch (const std::out_of_range& e) {
            std::cerr << "Error: Order value out of range in solver ID '" << ID 
                     << "'." << std::endl;
            return nullptr;
        }
    
        // 检查方法是否存在
        if (!callbacks_.count(method)) {
            std::cerr << "EquationSolver: No such method '" << method << "'." << std::endl;
            return nullptr;
        }
    
        return callbacks_[method](order);
    }
private:
    CallbackMap callbacks_;
private:
    EquationSolverFactory() = default;
    EquationSolverFactory(const EquationSolverFactory&) = default;
    EquationSolverFactory& operator = (const EquationSolverFactory&) = default;
    ~EquationSolverFactory() = default;
public:
    static EquationSolverFactory& getInstance() {
        static EquationSolverFactory instance;
        return instance;
    }
};

inline void RegisterAllEquationSolvers() {
    auto& fac = EquationSolverFactory::getInstance();
    fac.registerEquationSolver("AB", &createAdamsBashforthSolver);
    fac.registerEquationSolver("AM", &createAdamsMoultonSolver);
    fac.registerEquationSolver("BDF", &createBDFSolver);
    
    fac.registerEquationSolver("CRK", &createCRKSolver);
    fac.registerEquationSolver("ESDIRK", &createESDIRKSolver);
    fac.registerEquationSolver("GLRK", &createGLRKSolver);
    fac.registerEquationSolver("F", &createFSolver);
    fac.registerEquationSolver("DP", &createDPSolver);
}
// 注册完成，可以创建实例了