//
// Created by xwellen on 21.10.2021.
//

#ifndef VMLAB6_ODE_H
#define VMLAB6_ODE_H
#include <vector>
#include <functional>
typedef std::vector<std::pair<double, double>> data_t;

/*
[a:b]-интервал дифференцирования
h - шаг
alpha - начальные условия y(a)=alpha
f - функция данного ОДУ
*/
data_t euler(double a, double b, double h, double alpha, const std::function<double(double, double)>& f);
data_t runge(double a, double b, double h, double alpha, const std::function<double(double, double)>& f);
//F - функция точного аналитического решения
void print_table(const data_t &data, std::function<double(double, double)> f, std::function<double(double)> F);

#endif //VMLAB6_ODE_H
