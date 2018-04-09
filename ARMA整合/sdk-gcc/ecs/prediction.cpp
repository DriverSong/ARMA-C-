#include "prediction.h"
#include "Control.h"
#include "ARMA.h"
#include <cmath>
#include <iostream>

void PredictAll(std::vector<std::vector<double> >& AlldataArray, int period, int numFla, int* predict)
{
    int cnt = 7;//预测一天时选取最优的cnt个模型进行预测，再取平均。
    for(int n = 0; n < numFla; n++){//最大循环，共numfla个服务器
        //第n个服务器的差分预处理
        std::vector<double> dataArray;
        dataArray.assign(AlldataArray[n].begin(), AlldataArray[n].end());
        int diffGap = 7;
        std::vector<double> arrayDiff(preDiff(dataArray, diffGap));

        std::vector<double> predictALL(period);//period天的预测值

        for(int i = 0; i < period; i++){
            //每次预测一天，共预测7次
            std::vector<std::vector<int>> usedPQ;//存放使用过的p、q对
            std::vector<double> bestPredict(cnt);//存放最优的预测值
            double predictAvg = 0.0;

            for(int j = 0; j < cnt; j ++) {
//                int maxP = (int) std::log(arrayDiff.size());
//                int maxP = (int)pow(round(12 * arrayDiff.size() / 100.0), 0.25);
//                if(maxP > arrayDiff.size() - 1){
//                    maxP = arrayDiff.size() -1;
//                }
//                int maxQ = 0;
                int maxP = 10;
                int maxQ = 10;
                int maxPQ = 7;

                std::vector<int> bestModel(predictByAIC(arrayDiff, maxP, maxQ, maxPQ, usedPQ,j));



                int bestP = bestModel[0];
                int bestQ = bestModel[1];
//                int bestP = (int)pow(round(12 * arrayDiff.size() / 100.0), 0.25);
//                if(bestP > arrayDiff.size() - 1){
//                    bestP = arrayDiff.size() -1;
//                }
//
//                int bestQ = 0;

                std::vector<int> PQ(2);
                PQ[0] = bestP;
                PQ[1] = bestQ;
                usedPQ.push_back(PQ);

                //        int bestP = (int)pow(round(12 * arrayDiff.size() / 100.0), 0.25);
                //        if(bestP > arrayDiff.size() - 1){
                //            bestP = arrayDiff.size() -1;
                //        }
                //        int bestQ = 0;

                ARMA arma = ARMA(arrayDiff, bestP, bestQ);
                arma.calculateARMACoe();

                std::vector<double> predictOneDay(arma.predict(1));
                bestPredict[j] = predictOneDay[0];
            }
            for(int k = 0; k < bestPredict.size(); k++){
                predictAvg += bestPredict[k];
            }
            predictAvg = predictAvg/bestPredict.size();

            predictALL[i] = predictAvg;

            arrayDiff.push_back(predictAvg);
        }

        std::vector<double> arrayForDiffReverse;
        for(int i = diffGap - 1; i >= 0; i--){
            arrayForDiffReverse.push_back(dataArray[dataArray.size() - i]);
        }
        for(int i = 0; i < predictALL.size(); i++){
            arrayForDiffReverse.push_back(predictALL[i]);
        }

        std::vector<double> predictDiffReverse(contrastDiff(arrayForDiffReverse, diffGap));


        for(int i = 0; i < predictDiffReverse.size(); i++){
//            std::cout<<predictDiffReverse[i] << ";";
            if(predictDiffReverse[i] > 0){
                predict[n] += predictDiffReverse[i];
            }
        }
        std::cout<< std::cout;
        predict[n] = (int)(predict[n] + 0.5);
    }
}