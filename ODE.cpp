//
// Created by xwellen on 21.10.2021.
//

#include "ODE.h"

//Метод Эйлера
data_t euler(double a,
             double b,
             double h,
             double alpha,
             const std::function<double(double, double)>& f) {
    data_t data_euler;
    data_euler.clear();
    double xi = a;
    double yi = alpha;
    int n = (int) (round((b - a) / h));
    data_euler.emplace_back(xi, yi);
    for (int i = 0; i < n; i++) {
        yi += h * f(xi, alpha);
        xi += h;
        data_euler.emplace_back(xi, yi);
        alpha = yi;
    }
    return data_euler;
}

//Метод Рунге-Кутта 4-го порядка
data_t runge(double a,
             double b,
             double h,
             double alpha,
             const std::function<double(double, double)>& f) {
    double n = (b - a) / h;
    double X[(int) n];
    double Y1[(int) n];
    double Y2[(int) n];
    double Y3[(int) n];
    double Y4[(int) n];
    double Y[(int) n];
    //calculate
    X[0] = a;
    Y[0] = alpha;
    for (int i = 1; i <= n; i++) {
        X[i] = a + i * h;
        Y1[i] = h * f(X[i - 1], Y[i - 1]);
        Y2[i] = h * f(X[i - 1] + h / 2.0, Y[i - 1] + Y1[i] / 2.0);
        Y3[i] = h * f(X[i - 1] + h / 2, Y[i - 1] + Y2[i] / 2);
        Y4[i] = h * f(X[i - 1] + h, Y[i - 1] + Y3[i]);
        Y[i] = Y[i - 1] + (Y1[i] + 2 * Y2[i] + 2 * Y3[i] + Y4[i]) / 6;
    }
    data_t data_runge;
    data_runge.clear();
    for (int i = 0; i <= n; i++) {
        data_runge.emplace_back(X[i], Y[i]);
    }
    data_runge[0].first = a;
    return data_runge;
}

void print_table(const data_t &data,
                 std::function<double(double, double)> f,
                 std::function<double(double)> F) {
    const char *title_format = "%3s %10s %10s %10s %30s %24s\n";
    const char *row_format = "%3d % 10.5f % 10.5f % 10.5f % 16.5f % 13.5f\n";
    printf(title_format,
           "i",
           "xi",
           "yi",
           "f()",
           "Точное значение",
           "Погрешность");
    for (int i = 0; i < data.size(); i++) {
        printf(row_format,
               i,
               data[i].first,
               data[i].second,
               f(data[i].first, data[i].second),
               F(data[i].first),
               abs(data[i].second - F(data[i].first)));
    }
}