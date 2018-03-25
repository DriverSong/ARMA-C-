//
// Created by qiujiawei on 18-3-25.
//

#ifndef MYARMA_PREDICT_H
#define MYARMA_PREDICT_H

#include <vector>


class Predict {
public:
    std::vector<int> predictByAIC(std::vector<double> array, int maxP, int maxQ);
};


#endif //MYARMA_PREDICT_H
