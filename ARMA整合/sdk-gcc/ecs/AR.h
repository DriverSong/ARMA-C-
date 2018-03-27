//
// Created by qiujiawei on 18-3-22.
//

#ifndef MYARMA_AR_H
#define MYARMA_AR_H

#include <vector>
#include <cmath>
#include "ARMACore.h"
#include "Control.h"


class AR {
private:
    std::vector<double> array;
    int p;
    std::vector<std::vector<double>> vec;//[var,a[1,2,...,order]]

public:
    AR(std::vector<double> array, int p){
        this->array.assign(array.begin(), array.end());
        this->p = p;
    }

    std::vector<std::vector<double>> calculateCoe(){
        ARMACore ar;
        std::vector<double> arCoe = ar.ARCoe(this->array, this->p);
        this->vec.push_back(arCoe);
        return this->vec;
    }

    double calculateAIC(){
        double aic;
        ARMACore ar;
        aic = ar.calculateAIC(this->array, this->vec, 0);
        return aic;
    }

    std::vector<double> predictSum(int count){
        //[var,a[1,2,...,order]]
        double var = this->vec[0][0];
        std::vector<double> predict(count);
        std::vector<double> arrayTemp(this->array);
        for(int i = 0; i< count; i++){
            int size = arrayTemp.size();
            predict[i] = std::sqrt(var) * gaussrand();
            for(int j = 0; j < this->p; j++){
                predict[i] += arrayTemp[size - 1 - j] * this->vec[0][j];
            }
            arrayTemp.push_back(predict[i]);
        }

        return predict;
    }
};


#endif //MYARMA_AR_H
