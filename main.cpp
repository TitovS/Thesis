//
//  main.cpp
//  Thesis
//
//  Created by Сергей Титов on 04/05/15.
//  Copyright (c) 2015 Сергей Титов. All rights reserved.
//

#include "Lorenz.h"
#include "iostream"


int main (void) {

    Stat S;

    S.n_cycle = 0;
    S.u_cycle = 0;
    S.l_cycle = new double [100000]; //TODO Оформить динмаичесские массивы
    S.s_cycle = new int [100000];


    Lorenz L;
    L.SetParam(8. / 3., 10, 28); // b, sigma, r

    double a = 0.065;
    double eps = 0.001;

    for (int i = 0; i < 50; ++i) {

        L.GetLine();
        double *x = new double[L.m_dim];
        memcpy(x, L.dot, L.m_dim * sizeof(double));


        double t = 0.;

        for (int k = 0; k < 200; ++k) {

            L.GetTr(x, a, 10, &S, i); // шаг метода, величина решетки, шаг на решетке

            for (int j = 0; j < L.m_dim; ++j) {
                x[j] = L.dot[j] + L.vector[j] * t;
            }
            t += 0.065;
        }

        //L->Save(&S, i);
        //std::cout << "\n";

        L.Reset();
        S.Reset();
        a-=0.001;
    }

    return (0);
}
