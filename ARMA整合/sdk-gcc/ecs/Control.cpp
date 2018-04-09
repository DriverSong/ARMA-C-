//
// Created by qiujiawei on 18-3-25.
//
#include <cmath>
#include "ARMA.h"
#include "MA.h"
#include "AR.h"
#include <iostream>
#include "Control.h"

std::vector<int> predictByAIC(std::vector<double> array, int maxP, int maxQ, std::vector<std::vector<int>> usedPQ , int flag0) {
    int p = 0;
    int q = 0;
    std::vector<int> bestModel(2);
    double minAIC = 1.7976931348623157E308;
    if(maxP == 0 && maxQ == 0){
        std::cout << "p and p cannot be both zero!" << std::endl;
    } else{
        for(int i = 0; i <= maxP; i++){
            for(int j = 0; j <= maxQ; j++){
                int flag = 1;
                if(flag0 != 0){
                    for(int k = 0; k < usedPQ.size(); k++){
                        if(i == usedPQ[k][0] && j == usedPQ[k][1]){
                            flag = 0;
                            break;
                        }
                    }
                }
                if(flag){
                    if(i == 0 && j == 0)continue;
                    double AIC;
                    if(i == 0){
                        MA ma(array,j);
                        ma.calculateCoe();
                        AIC = ma.calculateAIC();
                    } else if(j == 0){
                        AR ar(array, i);
                        ar.calculateCoe();
                        AIC = ar.calculateAIC();
                    } else{
                        ARMA arma(array, i, j);
                        arma.calculateARMACoe();
                        AIC = arma.calculateAIC();
                    }
                    if( AIC <= minAIC){
                        minAIC = AIC;
                        p = i;
                        q = j;
                    }
                }
            }
        }
    }
    bestModel[0] = p;
    bestModel[1] = q;
    return bestModel;
}

std::vector<double> preDiff(std::vector<double> array, int gap){
    int size = array.size();
    int newSize = size - gap;
    std::vector<double> arrayDiff(newSize);
    for(int i = newSize - 1; i >= 0; i--){
        arrayDiff[i] = array[i + gap] - array[i];
    }
    return arrayDiff;
}

std::vector<double> contrastDiff(std::vector<double> array, int gap){
    int size = array.size();
    int newSize = size - gap;
    std::vector<double> arrayReturn(newSize);
    for(int i = 0; i < newSize; i++){
        arrayReturn[i] = array[gap + i] + array[i];
    }
    return arrayReturn;
}

double gaussrand1(){
    static double V1, V2, S;
    static int phase = 0;
    double X;

    if ( phase == 0 ) {
        do {
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;

            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);

        X = V1 * sqrt(-2 * log(S) / S);
    } else
        X = V2 * sqrt(-2 * log(S) / S);

    phase = 1 - phase;

    return X;
}
