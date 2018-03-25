//
// Created by qiujiawei on 18-3-25.
//

#include "Predict.h"
#include <cmath>
#include "ARMA.h"
#include <iostream>

std::vector<int> Predict::predictByAIC(std::vector<double> array, int maxP, int maxQ) {
    int p = 0;
    int q = 0;
    std::vector<int> bestModel(2);
    double minAIC = 1.7976931348623157E308D;
    if(maxP == 0 && maxQ == 0){
        std::cout << "p and p cannot be both zero!" << std::endl;
    } else{
        for(int i = 0; i < maxP; i++){
            for(int j = 0; j < maxQ; j++){
                ARMA arma(array, i, j);
                arma.calculateARMACoe();
                double AIC = arma.calculateAIC();
                if( AIC <= minAIC){
                    minAIC = AIC;
                    p = i;
                    q = j;
                }
            }
        }
    }
    bestModel.push_back(p);
    bestModel.push_back(q);
    return bestModel;
}
