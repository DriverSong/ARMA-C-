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

    std::vector<double> predictSum(int count){
        //std::cout << "Coe:" << std::endl;
//        for(int i = 0; i < this->vec[0].size(); i++){
//            std::cout << this->vec[0][i] << ";";
//        }
//        std::cout << std::endl;
        std::vector<double> predict(count);
        if(this->q == 0){
            //AR
            //[var,a[1,2,...,p]]
            double var = this->vec[0][0];
            std::vector<double> arrayTemp(this->array);
//            std::cout << "size1" << arrayTemp.size() << std::endl;
//            std::cout << "predict:";
            for(int i = 0; i< count; i++){
                int size = arrayTemp.size();
                predict[i] = std::sqrt(var) * gaussrand();
//                std::cout << predict[i] << ">>>>>";
                for(int j = 0; j < this->p; j++){
                    predict[i] += arrayTemp[size - 1 - j] * this->vec[0][j + 1];
                }
//                std::cout << i << ":" << predict[i] << std::endl;
                arrayTemp.push_back(predict[i]);
            }
        } else if(this->p == 0){
            //MA
            //[var, avg, 1(belta[0]), belta[1,2,...,q]]
            int size = this->array.size();
            int errLength = count + this->q;
            std::vector<double> MACoe = this->vec[0];
            double var = MACoe[0];
            double avg = MACoe[1];
            std::vector<double> errData(errLength);

            for(int i = 0; i < errLength; i++){
                errData[i] = std::sqrt(var) * gaussrand();
            }
            for(int i = 0; i < count; i ++){
                predict[i] = errData[this->q + i] + avg;
                for(int j = 0; j < this->q; j++){
                    predict[i] += errData[this->q + i - 1 - j] * MACoe[3 + j];
                }
            }
        } else{
            //ARMA
            //[c,1,2,...,p;
            // var,1,2,...,q]
            int errLength = count + this->q;
            double c = this->vec[0][0];
            double var = this->vec[1][0];
            std::vector<double> arrayTemp(this->array);
            std::vector<double> errData(errLength);
            std::vector<double> ARCoe(this->p);
            std::vector<double> MACoe(this->q);
            for(int i = 0; i < this->p; i++){
                ARCoe[i] = this->vec[0][i + 1];
            }
            for(int i = 0; i < this->q; i++){
                MACoe[i] = this->vec[1][i + 1];
            }

            for(int i = 0; i < errLength; i++){
                errData[i] = std::sqrt(var) * gaussrand();
            }
            for(int i = 0; i < count; i++){
                int size = arrayTemp.size();
                predict[i] = c + errData[this->q];
                for(int j = 0; j < this->p; j++){
                    predict[i] += ARCoe[j] * arrayTemp[size - 1 - j];
                }
                for(int j = 0; j < this->q; j++){
                    predict[i] += MACoe[j] * errData[this->q + i - 1 -j];
                }
                arrayTemp.push_back(predict[i]);
            }
        }
//        std::cout << "------------" << std::endl;
        return predict;
    }
};


#endif //MYARMA_ARMA_H
