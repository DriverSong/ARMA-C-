#include "prediction.h"
#include "Control.h"
#include "ARMA.h"
#include <cmath> 
#include <iostream>

void PredictAll(std::vector<std::vector<double> >& AlldataArray, int period, int interval, int numFla, int* predict)
{
    int cnt = 5;//预测一天时选取最优的cnt个模型进行预测，再取平均。
    for(int n = 0; n < numFla; n++){//最大循环，共numfla个服务器
        //第n个服务器的差分预处理
        int weekPredict = 0;
        int intervalPredict = 0;
        std::vector<double> dataArray;
        dataArray.assign(AlldataArray[n].begin(), AlldataArray[n].end());
        int diffGap = 7;//差分周期   
        std::vector<double> predictALL(period);//period天的预测值
        //std::vector<double> array1;//去掉原始数据末尾N位
        //for(int i = 0; i < dataArray.size(); i++){
        //    tempdataArray[i] = dataArray[i];
        //}
        //设置缺省时间段
        int dataArrayNums = dataArray.size();
        int aheaddays = dataArrayNums - 0;
        //将去掉尾部的数据赋值给tempdataArray
        //for (int i = 0; i < dataArray.size() / 10 ; i++)
        //{
        //array1.push_back(1);
        
        std::cout<<dataArray[9]<<std::endl;
        //}
        //array1.push_back(1);
        //tempdataArray.assign(dataArray.begin(), dataArray.end());
        //std::cout<<dataArray[9]<<std::endl;
        //预测的天数改为period加上14天
        for(int i = 0; i < period; i++){
            std::vector<std::vector<int>> usedPQ;//存放使用过的p、q对
            std::vector<double> bestPredict(cnt);//存放最优的预测值
            std::vector<double> arrayDiff(preDiff(dataArray, diffGap));//生成差分序列
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
            int aftdeal = (int)std::round(predictAvg + dataArray[dataArray.size() - 7]);  
            if (aftdeal < 0)
            {
                aftdeal = 0;
            }        

            if(i<dataArrayNums-aheaddays)
            {
                if(abs(aftdeal - dataArray[aheaddays+i]) > 22)
                {
                    //tempdataArray.push_back(dataArray[aheaddays+i]);

                }
                else
                {
                    //tempdataArray.push_back(dataArray[aheaddays+i]);
                }   
            }
            else
            {
                if(i<interval)
                {
                    intervalPredict += aftdeal;;
                }
                weekPredict += aftdeal;
                dataArray.push_back(double(aftdeal));
            }
            
            /*if(tempdataArray.size() > dataArrayNums)
            {
                weekPredict += aftdeal;
                std::cout<<aftdeal<<std::endl;

            }
            */
            //weekPredict += aftdeal;
            //dataArray.push_back(double(aftdeal));      
        } 
        //总预测长度减去中间间隔长度
        predict[n] = weekPredict-intervalPredict;
    }
}