//
// Created by qiujiawei on 18-3-22.
//

#ifndef MYARMA_MA_H
#define MYARMA_MA_H

#include <vector>
#include "ARMACore.h"
#include <cmath>
#include "ARMACore.cpp"

class MA {
private:
    std::vector<double> array;
    int q;
    std::vector<std::vector<double>> vec;//[var, avg, 1(belta[0]), belta[1,2,...,q]]

public:
    MA(std::vector<double> array, int q){
        this->array.assign(array.begin(),array.end());
        this->q = q;
    }

    std::vector<std::vector<double>> calculateCoe(){
        ARMACore ma;
        std::vector<double> maCoe = ma.MACoe(this->array, this->q);
        this->vec.push_back(maCoe);
        return this->vec;
    }

    double calculateAIC(){
        double aic;
        ARMACore ma;
        aic = ma.calculateAIC(this->array,this->vec, 1);
        return aic;
    }

    std::vector<double> predictSum(int count){
        int size = this->array.size();
        std::vector<double> MACoe = this->vec[0];
        std::vector<double> predict(count);
        double var = MACoe[0];
        double avg = MACoe[1];
        int errLength = count + this->q;
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
        return predict;
    }
};


#endif //MYARMA_MA_H
