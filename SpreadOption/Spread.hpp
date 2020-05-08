//
//  Spread.hpp
//  SpreadOption
//
//  Created by Kevin Zhang on 5/6/20.
//  Copyright © 2020 Kevin Zhang. All rights reserved.
//

#ifndef Spread_hpp
#define Spread_hpp

#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <armadillo>
#include <vector>
#include <complex>
#include <algorithm>
#include <typeinfo>


using namespace std;
using namespace arma;

class FFTSpread{
    private:
        int m_N;
        double m_X1, m_X2, m_K, m_rho, m_tau,m_r;
        mat m_sigma;
        mat m_mu;
        double u_bar = 40.0;
        double e1 = -3.0, e2 = 1.0;
    
        double m_ubar1, m_ubar2;
        
        double m_eta1, m_eta2;
        double m_etastar1, m_etastar2;
    
        vector<double> m_u1, m_u2;
        vector<double> m_x1, m_x2;
    
        double m_p1, m_p2;
    
    public:
        FFTSpread(int N, double S1, double S2, double K, double r, double sigma1, double sigma2, double rho, double tau);
        complex <double> Phi(complex <double> &mu, complex <double> &Sigma);
        complex <double> P_hat(complex <double> &u1, complex <double> &u2);
        int Search(vector<double> &v, double &target);
    
        double Delta1();
        double Delta2();
        vector<double> Greek();
    
};

#endif /* Spread_hpp */
