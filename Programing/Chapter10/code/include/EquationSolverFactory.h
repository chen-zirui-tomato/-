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
        //解析ID字符串，获取order，如“AB4” -> method = "AB", order = 4
        std::string method;
        int order = 0;
        size_t pos = 0;
        while (pos < ID.size() && std::isdigit(ID[pos])) {
            method += ID[pos];
            pos++;
        }
        if(pos < ID.size()){
            order = std::stoi(ID.substr(pos));
        }
        if(!callbacks_.count(method)){
            std::cerr << "EquationSolver:: No such Equation Solver called '" << ID << "'." << std::endl;
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
}

// 注册完成，可以创建实例了