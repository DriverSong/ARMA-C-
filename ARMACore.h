//
// Created by qiujiawei on 18-3-21.
//

#ifndef MYARMA_ARMACORE_H
#define MYARMA_ARMACORE_H

#include <vector>


class ARMACore{
public:
    double arraySum(std::vector<double> array);
    double arrayAvg(std::vector<double> array);
    double arrayVar(std::vector<double> array);//方差
    double arrayStd(std::vector<double> array);//标准差

    std::vector<double> autoCor(std::vector<double> array, int order);//自相关
    std::vector<double> autoCov(std::vector<double> array, int order);//自协方差

    std::vector<double> Levison(std::vector<double> autoCov);

    std::vector<double> ARCoe(std::vector<double> array, int p);
    std::vector<double> MACoe(std::vector<double> array, int q);
    std::vector<double> ARMACoe(std::vector<double> array, int p, int q);

    double AIC(std::vector<std::vector<double>> coe, std::vector<double> array, int type);

};

#endif //MYARMA_ARMACORE_H
