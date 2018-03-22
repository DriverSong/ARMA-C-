//
// Created by qiujiawei on 18-3-21.
//

#include <cmath>
#include <cstdlib>
#include <cstring>
#include "ARMACore.h"


double ARMACore::arraySum(std::vector<double> array){
    int size = array.size();
    double sum = 0.0;
    for(int i = 0; i < size; i++){
        sum += array[i];
    }
    return sum;
}

double ARMACore::arrayAvg(std::vector<double> array){
    int size = array.size();
    double sum = 0.0;
    double avg = 0.0;
    sum = this->arraySum(array);
    avg = sum/size;
    return avg;
}

//方差
double ARMACore::arrayVar(std::vector<double> array){
    int size = array.size();
    double avg = this->arrayAvg(array);
    double var = 0.0;
    for(int i = 0; i < size; i++){
        var += (array[i] - avg) * (array[i] - avg);
    }
    var = var/size;
    return var;
}

//标准差
double ARMACore::arrayStd(std::vector<double> array){
    double arrayVar = this->arrayVar(array);
    double arrayStd = std::sqrt(arrayVar);
    return arrayStd;
}

//自相关
std::vector<double> ARMACore::autoCor(std::vector<double> array, int order){
    std::vector<double> arrayAutoCor = this->autoCov(array, order);
    double var = this->arrayVar(array);
    int size = arrayAutoCor.size();
    for(int i = 0; i < size; i++){
        arrayAutoCor[i] = arrayAutoCor[i] / var;
    }
    return arrayAutoCor;
}

//自协方差
//根据order（p || q || p+q）计算需要的自协方差阶数
//协方差由levison递推得到系数（a(1),a(2),...a(order)）以及常数项，共order+1项。
std::vector<double> ARMACore::autoCov(std::vector<double> array, int order){
    std::vector<double> arrayAutoCov(order + 1);//包含常数项与order个系数，所以需+1
    double avg = this->arrayAvg(array);
    for(int i = 0; i < order + 1; i++){
        arrayAutoCov[i] = 0.0;
        for (int j = 0; j < array.size() - i; j++) {
            arrayAutoCov[i] += (array[j] - avg) * (array[j + i] - avg);
        }
        arrayAutoCov[i] = arrayAutoCov[i] / (array.size() - i);
    }
    return arrayAutoCov;//cov[0,1,...,order]
}

//Levison递推公式,计算ARMA的系数
//参考网址：https://wenku.baidu.com/view/48f74ec97d1cfad6195f312b3169a4517723e50e.html
std::vector<double> ARMACore::Levison(std::vector<double> autoCov){
    int order = autoCov.size() - 1;
    std::vector<std::vector<double>> tempCoe;//[var(0),
                                             // var(1), a(11)
                                             // var(2), a(21), a(22)
                                             // ...
                                             // var(size-1), a((size-1)1),...,a((size-1)(size-1))]

    tempCoe.resize(2);
    for(int i = 0; i < tempCoe.size(); i++){
        tempCoe[i].resize(order + 1);
    }


    tempCoe[0][0] = autoCov[0];
    for(int i = 0; i < order; i++){
        if(i == 0){
            tempCoe[1][1] = autoCov[1] / tempCoe[0][0];
            tempCoe[1][0] = tempCoe[0][0]*(1 - tempCoe[1][1]*tempCoe[1][1]);
        }
        else{
            double temp1 = 0.0;
            double temp2 = 0.0;
            for(int j = 0; j < i; j++){
                temp1 += autoCov[i - j]*tempCoe[0][j + 1];
                temp2 += autoCov[j + 1]*tempCoe[0][j + 1];
            }

            tempCoe[1][i + 1] = (autoCov[i + 1] - temp1) / (autoCov[0] - temp2);
            for (int j = 1; j < i + 1; j++) {
                tempCoe[1][j] = tempCoe[0][j] - tempCoe[1][i+1]*tempCoe[0][i+1-j];
            }
            tempCoe[1][0] = tempCoe[0][0] * (1 - tempCoe[i+1][i+1]*tempCoe[i+1][i+1]);
        }
        for(int j = 0; j <= i + 1; j++) {
            tempCoe[0][j] = tempCoe[1][j];
        }
    }
    return tempCoe[0];
}

//计算AR系数
std::vector<double> ARMACore::ARCoe(std::vector<double> array, int p){
    std::vector<double> arrayAutoCov = this->autoCov(array, p);//cov[0,1,..,p]
    std::vector<double> ARCoe = this->Levison(arrayAutoCov);
    return ARCoe;//[var,a[1,2,...,order]]
}

//计算MA系数
//参考https://wenku.baidu.com/view/555776cd2cc58bd63186bd42.html?pn=51NaN
//涉及逆相关函数的计算
std::vector<double> ARMACore::MACoe(std::vector<double> array, int q) {
    int pTemp = (int) std::sqrt(array.size());

    double arrayAvg = this->arrayAvg(array);
//    std::vector arrayZeroMean(array.size());
//    for(int i = 0; i < array.size(); i++){
//        arrayZeroMean[i] = array[i] - arrayAvg;
//    }//由于自协方差计算时会减去均值，所以arrayZeroMean和array的自协方差函数相同
    std::vector<double> arrayAutoCov = this->autoCov(array, pTemp);
    std::vector<double> alpha = this->Levison(arrayAutoCov);
    int alphaSize = alpha.size();
    double alphaVar = alpha[0];
    alpha[0] = -1;
    std::vector<double> alphaCor(q + 1);//样本的逆相关函数alphaCor[0,1,2,...,q]
    for(int i = 0; i < q + 1; i++){
        double sum = 0.0;
        for(int j = 0; j < alphaSize - i; j++){
            sum += alpha[j] * alpha[j + i];
        }
        alphaCor[i] = sum/alphaVar;
    }
    std::vector<double> beltaTemp = this->Levison(alphaCor);//[1/var,-belta[1,2,...,q]]
    std::vector<double> belta(beltaTemp.size() + 2);//[var, avg, 1(belta[0]), belta[1,2,...,q]]
    belta[0] = 1 / beltaTemp[0];
    belta[1] = arrayAvg;
    belta[2] = 1
    for(int i = 1; i < belta.size(); i ++){
        belta[i + 2] = -belta[i];
    }
    return belta;//[var, avg, 1(belta[0]), belta[1,2,...,q]]
}

//计算ARMA系数
//参考https://wenku.baidu.com/view/555776cd2cc58bd63186bd42.html?pn=101NaN
//求AR系数用到了逆矩阵
//求MA系数和MA模型的求法类似
std::vector<double> ARMACore::ARMACoe(std::vector<double> array, int p, int q) {
    std::vector<double> arrayCov = this->autoCov(array, p + q);//Cov[0,1,2,...p+q]
    std::vector<double> ARCov(p);//Cov[q+1,q+2,...,q+p]
    std::vector<std::vector<double>> matrixCov;
    matrixCov.resize(p);

    //ARCov
    for(int i = 0; i < ARCov.size(); i++){
        ARCov[i] = arrayCov[q + i + 1];
    }

    //matrixCov
    for(int i = 0; i < matrixCov.size(); i++){
        matrixCov[i].resize(p);
    }
    for(int i = 0; i < p; i++){
        for(int j = 0; j < p; j++){
            int indexTemp = std::abs(i + q - j);
            matrixCov[i][j] = arrayCov[indexTemp];
        }
    }

    std::vector<double> alpha(p);//=matrixCov^(-1)*ARCov
                                 //alpha[1,2,...,p]

    std::vector<double> alphaAddZero(p+1);//alpha[0,1,2,...,p]
    alphaAddZero[0] = -1;
    for(int i = 0; i < alpha.size(); i++){
        alphaAddZero[i + 1] = alpha[i];
    }

    std::vector<double> MACov(q + 1);
    for(int k = 0; k < q + 1; k++){
        MACov[k] = 0;
        for(int i = 0; i < alphaAddZero.size(); i++){
            for(int j = 0; j < alphaAddZero.size(); j++){
                int indexCov = std::abs(k+j-i);
                MACov[k] += alphaAddZero[i]*alphaAddZero[j]*arrayCov[indexCov];
            }
        }
    }

    

}
