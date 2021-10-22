//
// Created by xwellen on 21.10.2021.
//

#include "ODE.h"

//Метод Эйлера
data_t euler(double a,
             double b,
             double h1,
             double alpha,
             double eps,
             const std::function<double(double, double)>& f) {
    data_t data_euler;
    double h = h1;
    while (true){
        data_euler.clear();
        double xi = a;
        double yi = alpha;
        int n = (int) (round((b - a) / h));
        data_euler.emplace_back(xi, yi);
        double maximum_error = -1;
        for (int i = 0; i < n; i++) {
            yi += h * f(xi, yi);
            xi += h;
            data_euler.emplace_back(xi, yi);
            if (data_euler.size()>1){
                double s1 = data_euler[data_euler.size()-1].second;
                double s2 = data_euler[data_euler.size()-2].second;
                double error = abs(s1-s2);
                maximum_error = std::max(maximum_error, error);
            }
        }
        if (maximum_error <= eps) return data_euler;
        else h /= 2;
    }
}

//Метод Рунге-Кутта 4-го порядка
data_t runge(double a,
             double b,
             double h1,
             double alpha,
             double eps,
             const std::function<double(double, double)>& f) {
    double h = h1;
    //calculate
    while (true){
        double n = (b - a) / h;
        double X[(int) n];
        double Y1[(int) n];
        double Y2[(int) n];
        double Y3[(int) n];
        double Y4[(int) n];
        double Y[(int) n];
        X[0] = a;
        Y[0] = alpha;
        double maximum_error = -1;
        for (int i = 1; i <= n; i++) {
            X[i] = a + i * h;
            Y1[i] = h * f(X[i - 1], Y[i - 1]);
            Y2[i] = h * f(X[i - 1] + h / 2.0, Y[i - 1] + Y1[i] / 2.0);
            Y3[i] = h * f(X[i - 1] + h / 2, Y[i - 1] + Y2[i] / 2);
            Y4[i] = h * f(X[i - 1] + h, Y[i - 1] + Y3[i]);
            Y[i] = Y[i - 1] + (Y1[i] + 2 * Y2[i] + 2 * Y3[i] + Y4[i]) / 6;
            maximum_error = std::max(maximum_error, abs(Y[i]-Y[i-1]));
        }
        if (maximum_error * (1/15.) <= eps){
            data_t data_runge;
            for (int i = 0; i <= n; i++) {
                data_runge.emplace_back(X[i], Y[i]);
            }
            data_runge[0].first = a;
            return data_runge;
        }
        else{
            h /= 2;
        }
    }
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