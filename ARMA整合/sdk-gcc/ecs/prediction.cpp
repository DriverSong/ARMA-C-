#include "prediction.h"
#include "Control.h"
#include "ARMA.h"
#include <cmath>
#include <iostream>

void PredictAll(std::vector<std::vector<double> >& AlldataArray, int period, int numFla, int* predict)
{
    int cnt = 5;//预测一天时选取最优的cnt个模型进行预测，再取平均。
    for(int n = 0; n < numFla; n++){//最大循环，共numfla个服务器
        //第n个服务器的差分预处理
        int weekPredict = 0;
        std::vector<double> dataArray;
        dataArray.assign(AlldataArray[n].begin(), AlldataArray[n].end());
        int diffGap = 7;//差分周期   
        std::vector<double> predictALL(period);//period天的预测值
        std::vector<double> tempdataArray;//去掉原始数据末尾N位
        //设置缺省时间段
        int dataArrayNums = dataArray.size();
        int aheaddays = dataArrayNums - 14;
        //将去掉尾部的数据赋值给tempdataArray
        for (int i = 0; i < aheaddays ; i++)
        {
            tempdataArray.push_back(dataArray[i]);
            std::cout<<dataArray[i]<<std::endl;
        } 
        //预测的天数改为period加上14天
        for(int i = 0; i < period + 14; i++){
            std::vector<std::vector<int>> usedPQ;//存放使用过的p、q对
            std::vector<double> bestPredict(cnt);//存放最优的预测值
            std::vector<double> arrayDiff(preDiff(tempdataArray, diffGap));//生成差分序列
            double predictAvg = 0.0;

            for(int j = 0; j < cnt; j ++) {
                int maxP = (int) std::log(arrayDiff.size());
                int maxQ = 7;

                std::vector<int> bestModel(predictByAIC(arrayDiff, maxP, maxQ, usedPQ,j));

                int bestP = bestModel[0];
                int bestQ = bestModel[1];
                std::vector<int> PQ(2);
                PQ[0] = bestP;
                PQ[1] = bestQ;
                usedPQ.push_back(PQ);

                ARMA arma = ARMA(arrayDiff, bestP, bestQ);
                arma.calculateARMACoe();

                std::vector<double> predictOneDay(arma.predict(1));
                bestPredict[j] = predictOneDay[0];
            }

            for(int k = 0; k < bestPredict.size(); k++){
                predictAvg += bestPredict[k];
            }
            predictAvg = predictAvg/(double)bestPredict.size();
            int aftdeal = (int)std::round(predictAvg + tempdataArray[tempdataArray.size() - 7]);  
            if (aftdeal < 0)
            {
                aftdeal = 0;
            }
            if(i<dataArrayNums-aheaddays)
            {
                if(abs(aftdeal - dataArray[aheaddays+i]) > 22)
                {
                    tempdataArray.push_back(double(aftdeal));

                }
                else
                {
                    tempdataArray.push_back(dataArray[aheaddays+i]);
                }   
            }
            else
            {
                tempdataArray.push_back(double(aftdeal));
            }
            if(tempdataArray.size() > dataArrayNums)
            {
                weekPredict += aftdeal;
            }
            
            //weekPredict += aftdeal;
            //dataArray.push_back(double(aftdeal));      
        }   
        predict[n] = weekPredict;
    }
}