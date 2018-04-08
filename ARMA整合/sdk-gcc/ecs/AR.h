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
};


#endif //MYARMA_AR_H
