//
//  Spread.cpp
//  SpreadOption
//
//  Created by Kevin Zhang on 5/6/20.
//  Copyright © 2020 Kevin Zhang. All rights reserved.
//

#include "Spread.hpp"
#include "Complex_Gamma.hpp"

FFTSpread::FFTSpread(int N, double S1, double S2, double K, double r, double sigma1, double sigma2, double rho, double tau){
    m_N = N;
    m_X1 = log(S1/K);
    m_X2 = log(S2/K);
    m_tau = tau;
    m_K = K;
    m_sigma << sigma1 * sigma1 << rho * sigma1 * sigma2 << endr
            << rho * sigma1 * sigma2 << sigma2 * sigma2;
    m_mu << r - 0.5 * sigma1 * sigma1  << endr
         << r - 0.5 * sigma2 * sigma2;
    m_r = r;
    
    int i = 0;
    while(m_ubar1 < u_bar && i<m_N){
        m_ubar1 = M_PI * (i - m_N/2) / m_X1;
        i++;
    };
    
    i = 0;
    while(m_ubar2 < u_bar && i<m_N){
        m_ubar2 = M_PI * (i - m_N/2) / m_X2;
        i++;
    };
    
    m_eta1 = 2 * m_ubar1 / m_N;
    m_eta2 = 2 * m_ubar2 / m_N;
    m_etastar1 = M_PI / m_ubar1;
    m_etastar2 = M_PI / m_ubar2;
    
    vector<double> p1,p2;
    for (i=0;i<m_N;i++){
        m_u1.push_back(-m_ubar1 + m_eta1 * i);
        m_u2.push_back(-m_ubar2 + m_eta2 * i);
        m_x1.push_back(-0.5 * N * m_etastar1 + m_etastar1 * i);
        m_x2.push_back(-0.5 * N * m_etastar2 + m_etastar2 * i);
        p1.push_back(abs((m_X1 + m_etastar1 * m_N * 0.5) / m_etastar1 - i));
        p2.push_back(abs((m_X2 + m_etastar2 * m_N * 0.5) / m_etastar2 - i));
    }
    auto min1 = *min_element(p1.begin(), p1.end());
    auto min2 = *min_element(p2.begin(), p2.end());
    m_p1 = Search(p1,min1);
    m_p2 = Search(p2,min2);
};

int FFTSpread::Search(vector<double> &v, double &target){
    vector<double>::iterator itr = find(v.begin(), v.end(), target);
    return int(distance(v.begin(), itr));
};


complex <double> FFTSpread::P_hat(complex <double> &u1, complex <double> &u2){
    const complex <double> img(0, 1);
    complex <double> c1 = (u1 + u2) * img - 1.0;
    complex <double> c2 = -img * u2;
    complex <double> c3 = img * u1 + 1.0;
    return cgamma(c1) * cgamma(c2) / cgamma(c3);
};

complex <double> FFTSpread::Phi(complex <double> &mu, complex <double> &Sigma){
    const complex <double> img(0, 1);
    return  exp(img * mu * m_tau - Sigma * m_tau * 0.5);
};

double FFTSpread::Delta1(){
    cx_mat u;
    cx_mat *pu = &u;
    cx_mat Sigma(1,1);
    cx_mat mu;
    cx_mat HH1(m_N,m_N);
    cx_mat C(m_N,m_N);
    mat Delta1Mat(m_N,m_N);

    const complex <double> img(0, 1);
    for(int i=0;i<m_N;i++){
        for(int j=0;j<m_N;j++){
            *pu << m_u1[i] + e1 * img << m_u2[j] + e2 * img;
            mu = u * m_mu;
            Sigma = u * m_sigma * (u.t());
            HH1(i,j) = (m_u1[i] * img - e1) * pow(-1, i+j+2) * Phi(mu(0,0), Sigma(0,0)) * P_hat(u(0),u(1));
            C(i,j) = pow(-1, i+j+2) * exp(-e1 * m_x1[i] - e2 * m_x2[j]) * m_eta1 * m_eta2 * pow((m_N/(2*M_PI)), 2);
        }
    }
    Delta1Mat = exp(-m_X1) * exp(-0.05 * m_tau) * real(C % ifft2(HH1));
    return Delta1Mat(m_p1,m_p2);
};

vector<double> FFTSpread::Greek(){
    const complex <double> img(0, 1);
    cx_mat u;
    cx_mat *pu = &u;
    cx_mat Sigma(1,1);
    cx_mat mu;
    complex <double> H;
    cx_mat H1(m_N,m_N);
    cx_mat H2(m_N,m_N);
    cx_mat T(m_N,m_N);
    cx_mat C(m_N,m_N);
    cx_mat HH11(m_N,m_N);
    cx_mat HH12(m_N,m_N);
    cx_mat HH22(m_N,m_N);
    mat Delta1Mat(m_N,m_N);
    mat Delta2Mat(m_N,m_N);
    mat ThetaMat(m_N,m_N);
    mat Gamma11Mat(m_N,m_N);
    mat Gamma12Mat(m_N,m_N);
    mat Gamma22Mat(m_N,m_N);
    vector<double> Greek;

    for(int i=0;i<m_N;i++){
        for(int j=0;j<m_N;j++){
            *pu << m_u1[i] + e1 * img << m_u2[j] + e2 * img;
            mu = u * m_mu;
            Sigma = u * m_sigma * (u.t());
            H = pow(-1, i+j+2) * Phi(mu(0,0), Sigma(0,0)) * P_hat(u(0),u(1));
            H1(i,j) = (m_u1[i] * img - e1) * H;
            H2(i,j) = (m_u2[j] * img - e2) * H;
            T(i,j) = (mu(0,0) * img - 0.5 * Sigma(0,0) - m_r) * H;
            HH11(i,j) = (m_u1[i] * img - e1) * (m_u1[i] * img - e1) * H;
            HH12(i,j) = (m_u1[i] * img - e1) * (m_u2[j] * img - e2) * H;
            HH22(i,j) = (m_u2[j] * img - e2) * (m_u2[j] * img - e2) * H;
            C(i,j) = pow(-1, i+j+2) * exp(-e1 * m_x1[i] - e2 * m_x2[j]) * m_eta1 * m_eta2 * pow((m_N/(2*M_PI)), 2);
        }
    }
    Delta1Mat = exp(-m_X1) * exp(-m_r * m_tau) * real(C % ifft2(H1));
    Delta2Mat = exp(-m_X2) * exp(-m_r * m_tau) * real(C % ifft2(H2));
    ThetaMat = m_K * exp(-m_r * m_tau) * real(C % ifft2(T));
    Gamma11Mat = -(exp(-m_X1) / m_K) * Delta1Mat + (exp(-2.0 * m_X1) / m_K) * exp(-m_r * m_tau) * real(C % ifft2(HH11));
    Gamma12Mat = (exp(-m_X1) * exp(-m_X2) / m_K) * exp(-m_r * m_tau) * real(C % ifft2(HH12));
    Gamma22Mat = -(exp(-m_X2) / m_K) * Delta2Mat + (exp(-2.0 * m_X2) / m_K) * exp(-m_r * m_tau) * real(C % ifft2(HH22));
    
    Greek.push_back(Delta1Mat(m_p1,m_p2));
    Greek.push_back(Delta2Mat(m_p1,m_p2));
    Greek.push_back(ThetaMat(m_p1,m_p2));
    Greek.push_back(Gamma11Mat(m_p1,m_p2));
    Greek.push_back(Gamma12Mat(m_p1,m_p2));
    Greek.push_back(Gamma22Mat(m_p1,m_p2));
    return Greek;
};
