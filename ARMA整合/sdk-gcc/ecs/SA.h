//
// Created by qiujiawei on 18-4-8.
//

#ifndef SA_SIMULATED_ANNEALING_SA_H
#define SA_SIMULATED_ANNEALING_SA_H


#include <string>
#include <vector>

class SA {
private:
    double startT;//初始温度
    double endT;//结束温度
    double r;//降温系数
    int sumCPU;
    int sumMEM;
    int numFla;
    std::string target;
    int* vecFlaCPU;
    int* vecFlaMEM;
    int* vecFlaPre;
    int sumPredictCPU;//预测的总CPU
    int sumPredictMEM;//预测的总MEM
    std::vector<int> flaSeries;
    std::vector<int> des;
    std::vector<std::vector<int>> res;//分配结果
    int score;

    int calculateScore(std::vector<std::vector<int>> res);
    std::vector<int> newFlaSeries();
    std::vector<std::vector<int>> distribution(std::vector<int> flaSeries);
    int calculateCPU(std::vector<int> PHY);
    int calculateMEM(std::vector<int> PHY);

public:
    SA(double startT, double endT, double r, int sumCPU, int sumMEM, int numFla, std::string target,
       int* vecFlaCPU, int* vecFlaMEM, int* vecFlaPre);

    void calculate();


};


#endif //SA_SIMULATED_ANNEALING_SA_H
