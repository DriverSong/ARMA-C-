//
// Created by qiujiawei on 18-4-8.
//

#include "SA.h"
#include <algorithm>

int SA::calculateScore(std::vector<std::vector<int>> res) {
//    int sizePHY = res.size();
//    double score = 0.0;
//    if(this->target == "CPU"){
//        score = ((double)this->sumPredictCPU) / (sizePHY * this->sumCPU);
//    } else if(this->target == "MEM"){
//        score = ((double)this->sumPredictMEM) / (sizePHY * this->sumMEM);
//    }
//    return score;
    return res.size();
}

std::vector<int> SA::newFlaSeries() {
    std::vector<int> newFlaSeries(this->flaSeries);
    std::random_shuffle(this->des.begin(), this->des.end());//define in algorithm
    std::swap(newFlaSeries[this->des[0]], newFlaSeries[this->des[1]]);
    return newFlaSeries;
}

std::vector<std::vector<int>> SA::distribution(std::vector<int> flaSeries){
    std::vector<std::vector<int>> res;
    for(int i = 0; i < flaSeries.size(); i++){
        if (i == 0){//初始添加一个物理机
            std::vector<int> newPHY(this->numFla, 0);
            newPHY[flaSeries[i]]++;
            res.push_back(newPHY);
        } else{
            for(int j = 0; j < res.size(); j++){
                if(this->calculateCPU(res[j]) + this->vecFlaCPU[flaSeries[i]] <= this->sumCPU &&
                        this->calculateMEM(res[j]) + this->vecFlaMEM[flaSeries[i]] <= this->sumMEM){//flaSeries[i]可以存放在第j个物理机中
                    res[j][flaSeries[i]]++;
                    break;//找到存放位置，退出循环
                }
                if(j == res.size() - 1){//遍历所有物理机，均不能存放flaSeries[i]，新增一个物理机
                    std::vector<int> newPHY(this->numFla, 0);
                    newPHY[flaSeries[i]]++;
                    res.push_back(newPHY);
                }
            }
        }
    }
    return res;
}

int SA::calculateCPU(std::vector<int> PHY){
    int totalCPU = 0;
    for(int i = 0; i < PHY.size(); i++){
        totalCPU += PHY[i] * this->vecFlaCPU[i];
    }
    return totalCPU;
}

int SA::calculateMEM(std::vector<int> PHY){
    int totalMEM = 0;
    for(int i = 0; i < PHY.size(); i++){
        totalMEM += PHY[i] * this->vecFlaMEM[i];
    }
    return totalMEM;
}

SA::SA(double startT, double endT, double r, int sumCPU, int sumMEM, int numFla, std::string target, int *vecFlaCPU,
       int *vecFlaMEM, int *vecFlaPre) {
    this->startT = startT;
    this->endT = endT;
    this->r = r;
    this->sumCPU = sumCPU;
    this->sumMEM = sumMEM;
    this->numFla = numFla;
    this->target = target;
    this->vecFlaCPU =vecFlaCPU;
    this->vecFlaMEM = vecFlaMEM;
    this->vecFlaPre = vecFlaPre;

    this->sumPredictCPU = 0;
    this->sumPredictMEM = 0;
    for(int i = 0; i < this->numFla; i++){
        this->sumPredictCPU += vecFlaPre[i] * vecFlaCPU[i];
        this->sumPredictMEM += vecFlaPre[i] * vecFlaMEM[i];
    }

    for(int i = 0; i < this->numFla; i++){
        for(int j = 0; j < this->vecFlaPre[i]; j++){
            this->flaSeries.push_back(i);
        }
    }

    for(int i = 0; i < this->flaSeries.size(); i++){
        this->des.push_back(i);
    }

    this->res = this->distribution(this->flaSeries);

    this->score = this->calculateScore(this->res);
}

void SA::calculate(){
    double T = this->startT;
    while(T > this->endT){
        std::vector<int> newFlaSeries = this->newFlaSeries();
        std::vector<std::vector<int>> tempRES = this->distribution(newFlaSeries);
        int tempScore = this->calculateScore(tempRES);
        if(tempScore < this->score){
            this->flaSeries = newFlaSeries;
            this->res = tempRES;
            this->score = tempScore;
        } else if((std::exp(this->score - tempScore) / T) > (double)std::rand() / RAND_MAX){
            this->flaSeries = newFlaSeries;
            this->res = tempRES;
            this->score = tempScore;
        }
        T = T * this->r;
    }
}
