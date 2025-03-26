#ifndef FUNCTION
#define FUNCTION

#include <iostream>
#include <exprtk.hpp>
enum class ConditionType{
    Dirichlet,
    Mixed,
    Neumann
};

class Function{
public:
    virtual double operator()(double x, double y) const = 0;
    virtual ~Function(){}
};

class DynamicFunction : public Function {
    private:
        exprtk::parser<double> parser;
        exprtk::expression<double> expression;
        exprtk::symbol_table<double> symbol_table;
        mutable double x, y;
    
    public:
        explicit DynamicFunction(const std::string& expr_str) {
            symbol_table.add_variable("x", x);
            symbol_table.add_variable("y", y);
            expression.register_symbol_table(symbol_table);
    
            if (!parser.compile(expr_str, expression)) {
                throw std::runtime_error("表達式解析失敗: " + expr_str);
            }
        }
    
        double operator()(double x_val, double y_val)const override{
            x = x_val;
            y = y_val;
            return expression.value();
        }
    };
    
#endif