//
//  main.cpp
//  SpreadOption
//
//  Created by Kevin Zhang on 4/29/20.
//  Copyright © 2020 Kevin Zhang. All rights reserved.
//


#include<iostream>
#include<armadillo>

#include "Spread.hpp"


int main() {
    
    clock_t tStart = clock();
    vector<vector <double>> Greeks;
    vector <double> Greek;
    vector <double> delta;
    double delta1;
    for(int i=0;i<1;i++){
        FFTSpread spread(128,50,70,2,0.05,0.2,0.6,0.3,2);
        Greek = spread.Greek();
        Greeks.push_back(Greek);
    }
    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
}
