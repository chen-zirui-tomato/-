#ifndef DYNAMIC_FUNCTION_7_H
#define DYNAMIC_FUNCTION_7_H

#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <exprtk.hpp>

class Function_with_7_parameters {
public:
    virtual std::vector<double> operator()(const std::vector<double>& x_val, double t_val) const = 0;
    virtual ~Function_with_7_parameters() {}
};

class DynamicFunction_7 : public Function_with_7_parameters {
private:
    std::vector<std::string> raw_expr_list;

public:
    explicit DynamicFunction_7(std::ifstream& expr_file) {
        std::string expr;
        while (std::getline(expr_file, expr, ';')) {
            if (!expr.empty()) {
                raw_expr_list.push_back(expr);
            }
        }

        if (raw_expr_list.size() != 6) {
            throw std::runtime_error("表达式个数不正确，应为6个");
        }
    }

    std::vector<double> operator()(const std::vector<double>& x_val, double t_val) const override {
        if (x_val.size() != 6) {
            throw std::invalid_argument("输入向量必须为6维");
        }

        double u1 = x_val[0], u2 = x_val[1], u3 = x_val[2];
        double u4 = x_val[3], u5 = x_val[4], u6 = x_val[5];
        double t = t_val;

        exprtk::symbol_table<double> table;
        table.add_variable("u1", u1);
        table.add_variable("u2", u2);
        table.add_variable("u3", u3);
        table.add_variable("u4", u4);
        table.add_variable("u5", u5);
        table.add_variable("u6", u6);
        table.add_variable("t", t);
        table.add_constant("PI", 3.14159265358979323846);
        table.add_constant("miu", 0.012277471);

        exprtk::parser<double> parser;
        std::vector<double> result;
        result.reserve(6);

        for (const auto& expr_str : raw_expr_list) {
            exprtk::expression<double> expression;
            expression.register_symbol_table(table);
            if (!parser.compile(expr_str, expression)) {
                throw std::runtime_error("表达式解析失败: " + expr_str);
            }
            result.push_back(expression.value());
        }

        return result;
    }
};

class HardcodedFunction_7 : public Function_with_7_parameters {
    public:
        std::vector<double> operator()(const std::vector<double>& u, double t) const override {
            double u1 = u[0], u2 = u[1], u3 = u[2];
            double u4 = u[3], u5 = u[4], u6 = u[5];
            constexpr double miu = 0.012277471;
    
            double D1 = pow((u1 + miu)*(u1 + miu) + u2*u2 + u3*u3, 1.5);
            double D2 = pow((u1 - 1 + miu)*(u1 - 1 + miu) + u2*u2 + u3*u3, 1.5);
    
            std::vector<double> res(6);
            res[0] = u4;
            res[1] = u5;
            res[2] = u6;
            res[3] = 2*u5 + u1 - (1 - miu)*(u1 + miu)/D1 - miu*(u1 + miu - 1)/D2;
            res[4] = -2*u4 + u2 - (1 - miu)*u2/D1 - miu*u2/D2;
            res[5] = - (1 - miu)*u3/D1 - miu*u3/D2;
            return res;
        }
    };

#endif // DYNAMIC_FUNCTION_7_H
