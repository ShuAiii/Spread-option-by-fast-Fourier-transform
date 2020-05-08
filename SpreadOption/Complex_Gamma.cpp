//
//  Complex_Gamma.cpp
//  SpreadOption
//
//  Created by Kevin Zhang on 5/6/20.
//  Copyright © 2020 Kevin Zhang. All rights reserved.
//

#include "Complex_Gamma.hpp"

complex <double> round(complex <double> z){
    double real = z.real();
    double imag = z.imag();
    real = floor(real + 0.5);
    imag = floor(imag + 0.5);
    z.real(real);
    z.imag(imag);
    return z;
};

complex <double> cgamma(complex <double> &z){
    const double g = 607.0 / 128.0;
    complex <double> zz = z;
    
    if (z.real()<0){
        z = -1.0 * z;
    }
    vector<double> c;
    
    c.push_back(0.99999999999999709182);
    c.push_back(57.156235665862923517);
    c.push_back(-59.597960355475491248);
    c.push_back(14.136097974741747174);
    c.push_back(-0.49191381609762019978);
    c.push_back(.33994649984811888699e-4);
    c.push_back(.46523628927048575665e-4);
    c.push_back(-.98374475304879564677e-4);
    c.push_back(.15808870322491248884e-3);
    c.push_back(-.21026444172410488319e-3);
    c.push_back(.21743961811521264320e-3);
    c.push_back(-.16431810653676389022e-3);
    c.push_back(.84418223983852743293e-4);
    c.push_back(-.26190838401581408670e-4);
    c.push_back(.36899182659531622704e-5);
    
    z = z - 1.0;
    complex <double> zh = z + 0.5;
    complex <double> zgh = zh + g;
    complex <double> zp = pow(zgh, zh*0.5);
    complex <double> ss = 0;
    for(size_t pp=14;pp>=1;pp--){
        ss = ss + c[pp]/( pp*1.0 + z);
    };
    
    const double sq2pi=  2.5066282746310005024157652848110;
    complex <double> f;
    f = (sq2pi * (c[0] + ss))*((zp * exp(-zgh)) * zp);
    if(z==0.0||z==1.0){
        f= 1.0;
    }
    if (zz.real()<0){
        f = -M_PI / (zz * f * sin(M_PI * zz));
    }
    if(round(zz)==zz && zz.imag()==0 && zz.real()<=0){
        f = INFINITY;
    }
    return f;
};
