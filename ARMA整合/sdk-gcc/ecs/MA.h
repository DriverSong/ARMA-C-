//
// Created by qiujiawei on 18-3-22.
//

#ifndef MYARMA_MA_H
#define MYARMA_MA_H

#include <vector>
#include "ARMACore.h"
#include <cmath>
#include "Control.h"

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
};


#endif //MYARMA_MA_H
