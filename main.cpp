#include <iostream>
#include "ODE.h"
#include "gnuplot-iostream.h"

#define EXACT_SMOOTH 16.

using namespace std;

//данное ОДУ y'=f(x, y)
double f(double x, double y) {
    return y + (1 + x) * pow(y, 2);
}

//точное аналитическое решение данного выше уравнения при y(1)=1
double F(double x) {
    return -1 / x;
}

int main() {
    double a = 1, b = 1.5, h = 0.1, alpha = -1, eps = 0.01;

    cout << "Метод Эйлера: " << endl;
    data_t eulerData = euler(a, b, h, alpha, eps, f);
    printf("с шагом h=%f:\n", eulerData[1].first - eulerData[0].first);
    print_table(eulerData, f, F);

    cout << "Метод Рунге-Кутта: " << endl;
    data_t rungeData = runge(a, b, h, alpha, eps, f);
    printf("с шагом h=%f:\n", rungeData[1].first - rungeData[0].first);
    print_table(rungeData, f, F);


    //строим точки для точного решения
    data_t exactData;
    int n = (int) ((b - a) / (h / EXACT_SMOOTH));
    double x = a;
    for (int i = 0; i < n; i++) {
        exactData.emplace_back(x, F(x));
        x += h / EXACT_SMOOTH;
    }

    Gnuplot gp;
    gp << "set title 'Решение ОДУ численными методами'\n";
    gp << "set grid\n";
    gp << "set key bottom right\n";
    gp << "plot '-' with lines title 'Метод Эйлера' linewidth 4 lt 12, "
       << "'-' with lines title 'Метод Рунге-Кутта 4-го порядка 'linewidth 4 lt 19, "
       << "'-' with lines title 'Точное решение' linewidth 1 lt -1\n";
    gp.send1d(eulerData);
    gp.send1d(rungeData);
    gp.send1d(exactData);
    return 0;
}
