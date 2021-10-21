#include <iostream>
#include "ODE.h"
#include "gnuplot-iostream.h"

using namespace std;

//данное ОДУ y'=f(x, y)
/*double f(double x, double y){
    return y*(2-y)/(x+3);
}*/
double f(double x, double y){
    return (5+x)/pow(y,2);
}

//точное аналитическое решение данного выше уравнения
/*
double F(double x){
    return (2*pow((x+3),2))/(x*x+6*x+3);
}
*/
double F(double x){
    return pow(((3*x*x+30*x-17)/2),(1/3.));
}

int main() {
    //double a=1, b=4, h=0.5, alpha=2, eps=0.0001;
    double a,b,h,alpha, eps=0.0001;
    cout << "a=";
    cin >> a;
    cout << "b=";
    cin >> b;
    cout << "h=";
    cin >> h;
    cout << "y("<<a<<")=";
    cin >> alpha;


    cout << "Метод Эйлера: " << endl;
    data_t euler1,euler2;
    double h1 = h;
    do{
        euler1 = euler(a,b,h1,alpha, f);
        cout << "c шагом h=" << h1 << endl;
        print_table(euler1, f, F);

        euler2 = euler(a,b,h1/2,alpha, f);
        cout << "c шагом h=" << h1/2 << endl;
        print_table(euler2, f, F);

        h1 /= 2;
    }while((abs(euler1[euler1.size()-1].second-euler2[euler2.size()-1].second)/pow(2,1/eps)-1)>=eps);



    cout << "Метод Рунге-Кутта: " << endl;
    data_t runge1,runge2;
    double h2 = h;
    do{
        runge1 = runge(a,b,h2,alpha, f);
        cout << "c шагом h=" << h2 << endl;
        print_table(runge1, f, F);

        runge2 = runge(a,b,h2/2,alpha, f);
        cout << "c шагом h=" << h2/2 << endl;
        print_table(runge2, f, F);
        h2 /= 2;
    }while((abs(runge1[runge1.size()-1].second-runge2[runge2.size()-1].second)/pow(2,1/eps)-1)>=eps);

    //строим точки для точного решения
    data_t exact_data;
    int n = (int)((b-a)/(h/4));
    double x=a;
    for (int i=0; i<n; i++){
        exact_data.emplace_back(x,F(x));
        x+=h/4;
    }

    Gnuplot gp;
    gp << "set title 'Решение ОДУ численными методами'\n";
    gp << "set grid\n";
    gp << "plot '-' with lines title 'Метод Эйлера' linewidth 4 lt 12, "
       << "'-' with lines title 'Метод Рунге-Кутта 4-го порядка 'linewidth 4 lt 19, "
       << "'-' with lines title 'Точное решение' linewidth 1 lt -1\n";
    gp.send1d(euler2);
    gp.send1d(runge2);
    gp.send1d(exact_data);
    return 0;
}
