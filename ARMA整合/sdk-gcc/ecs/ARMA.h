//
// Created by qiujiawei on 18-3-22.
//

#ifndef MYARMA_ARMA_H
#define MYARMA_ARMA_H

#include <vector>
#include <cmath>
#include "MA.h"
#include "AR.h"
#include "Control.h"
#include <iostream>

class ARMA {
private:
    int p;
    int q;
    std::vector<std::vector<double>> vec;//[c,1,2,...,p;
                                         // var,1,2,...,q]
    std::vector<double> array;


    double gaussrand0(){
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


public:
    ARMA(std::vector<double> array, int p, int q){
        this->p = p;
        this->q = q;
        this->array.assign(array.begin(), array.end());
    }

    std::vector<std::vector<double>> calculateARMACoe(){
        if(this->q == 0){
            AR ar(this->array, this->p);
            this->vec = ar.calculateCoe();
        }else if(this->p == 0){
            MA ma(this->array, this->q);
            this->vec = ma.calculateCoe();
        } else{
            ARMACore arma;
            std::vector<double> armaCoe = arma.ARMACoe(this->array, this->p, this-> q);
            std::vector<double> arCoe(this->p+1);
            std::vector<double> maCoe(this->q+1);
            for(int i = 0; i < armaCoe.size(); i++){
                if(i < arCoe.size()){
                    arCoe[i] = armaCoe[i];
                } else {
                    maCoe[i - this->p - 1] = armaCoe[i];
                }
            }
            this->vec.push_back(arCoe);
            this->vec.push_back(maCoe);
        }
        return this->vec;
    }

    double calculateAIC(){
        ARMACore arma;
        double aic;
        if(this->q == 0){
            aic = arma.calculateAIC(this->array, this->vec, 0);
        } else if(this->p == 0){
            aic = arma.calculateAIC(this->array, this->vec, 1);
        } else{
            aic = arma.calculateAIC(this->array, this->vec, 2);
        }
        return aic;
    }

    std::vector<double> predict(int count){
        std::vector<double> predict(count);
        int p = this->p;
        int q = this->q;
        std::vector<double> array(this->array);

        if(this->q == 0){
            //AR
            //[var,a[1,2,...,p]]
            std::vector<double> ARCoe(p);
            for(int i = 0; i < p; i++){
                ARCoe[i] = this->vec[0][i + 1];
            }

            for(int i = 0; i < count; i++){
                int n = array.size();
                for(int j = 0; j < p; j++){
                    predict[i] += array[n - 1 - j] * ARCoe[j];
                }
                array.push_back(predict[i]);
            }
        } else if(this->p == 0){
            //MA
            //[var, avg, 1(belta[0]), belta[1,2,...,q]]
            double var = this->vec[0][0];
            std::vector<double> MACoe(q);
            for(int i = 0; i < q; i++){
                MACoe[i] = this->vec[0][i + 3];
            }

            std::vector<double> gaussNoise(q);

            int size = array.size();
            for(int i = 0; i < size + count; i++){
                if(i >= size){
                    for(int j = 0; j < q; j++){
                        predict[i - size] += gaussNoise[j] * MACoe[j];
                    }
                }

                for(int j = q - 1; j > 0; j--){
                    gaussNoise[j] = gaussNoise[j - 1];
                }
                gaussNoise[0] = this->gaussrand0() * std::sqrt(var);
            }

        } else{
            //ARMA
            //[c,1,2,...,p;
            // var,1,2,...,q]
            double var = this->vec[1][0];

            std::vector<double> ARCoe(p);
            std::vector<double> MACoe(q);
            for(int i = 0; i < p; i++){
                ARCoe[i] = this->vec[0][i + 1];
            }
            for(int i = 0; i < q; i++){
                MACoe[i] = this->vec[1][i + 1];
            }

            std::vector<double> gaussNoise(q);

            double size = array.size();
            for(int i = 0; i < size + count; i++){
                if(i >= size){
                    for(int j = 0; j < p; j++){
                        predict[i - size] += array[i - 1 - j] * ARCoe[j];
                    }

                    for(int j = 0; j < q; j++){
                        predict[i - size] += gaussNoise[j] * MACoe[j];
                    }
                }
                for(int j = q - 1; j > 0; j--){
                    gaussNoise[j] = gaussNoise[j - 1];
                }
                gaussNoise[0] = this->gaussrand0() * std::sqrt(var);

                array.push_back(predict[i - size]);
            }
        }
        return predict;
    }
};


#endif //MYARMA_ARMA_H
