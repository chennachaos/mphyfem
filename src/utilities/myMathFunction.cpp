#include "myMathFunction.h"


myMathFunction::myMathFunction()
{
    x = 0.0;
    y = 0.0;
    z = 0.0;
    t = 0.0;
}



myMathFunction::~myMathFunction()
{

    
    
}


myMathFunction::myMathFunction(const string& expr_)
{
    x = 0.0;
    y = 0.0;
    z = 0.0;
    t = 0.0;

    expr = expr_;

    //typedef exprtk::symbol_table<double> symbol_table_t;
    //typedef exprtk::parser<double>       parser_t;
    //symbol_table_t symbol_table;
    //parser_t parser;

    expression.release();

    exprtk::symbol_table<double> symbol_table;
    exprtk::parser<double>       parser;

    symbol_table.add_variable("x",x);
    symbol_table.add_variable("y",y);
    symbol_table.add_variable("z",z);
    symbol_table.add_variable("t",t);
    symbol_table.add_constants();
    expression.register_symbol_table(symbol_table);

    parser.compile(expr, expression);
}




void  myMathFunction::initialise(const string& expr_)
{
    expr = expr_;

    //typedef exprtk::symbol_table<double> symbol_table_t;
    //typedef exprtk::parser<double>       parser_t;
    //symbol_table_t symbol_table;
    //parser_t parser;

    expression.release();

    exprtk::symbol_table<double> symbol_table;
    exprtk::parser<double>       parser;

    symbol_table.add_variable("x",x);
    symbol_table.add_variable("y",y);
    symbol_table.add_variable("z",z);
    symbol_table.add_variable("t",t);
    symbol_table.add_constants();
    expression.register_symbol_table(symbol_table);

    parser.compile(expr, expression);

    return;
}




double  myMathFunction::getValue(double x_, double y_, double z_, double t_)
{
    x = x_; y = y_; z = z_; t=t_;

    return expression.value();
}









