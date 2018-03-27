#include "prediction.h"
#include "Control.h"
#include "ARMA.h"
#include <iostream>

void PredictAll(std::vector<std::vector<double> >& AlldataArray, int period, int numFla, int* predict)
{
	// std::vector<int> respredict;
//	std::cout << numFla << std::endl;
	for (int n = 0; n < numFla; n++)
	{
		double weekPredict = 0.0;
		//dataArray:���������ѵ������
		//LastWeek:���������ѵ�����ݵ������λ����
		std::vector<double> dataArray;
		//�����ȡ15̨����������ݵ������ڵ�һά������
		dataArray.assign(AlldataArray[n].begin(), AlldataArray[n].end());

		std::vector<double> arrayDiff(preDiff(dataArray, 7));


		//����arimaģ��
		//ARIMAModel* arima = new ARIMAModel(dataArray);
		int maxP = 7;
		int maxQ = 7;
		std::vector<int> bestModel(predictByAIC(arrayDiff, maxP, maxQ));

		int bestP = bestModel[0];
		int bestQ = bestModel[1];

//		std::cout<< "P" << bestP << std::endl;
//		std::cout<< "Q" << bestQ << std::endl;

//		std::cout << "size0" << arrayDiff.size() << std::endl;
		ARMA arma = ARMA(arrayDiff, bestP, bestQ);

		arma.calculateARMACoe();

		std::vector<double> predict1(arma.predictSum(period));

//		for(int i = 0; i < predict1.size(); i++){
//			std::cout << predict1[i] << ";";
//		}
//		std::cout << std::endl;
//		std::cout << "----------------" << std::endl;


		for(int i = 0; i < predict1.size(); i++){
			arrayDiff.push_back(predict1[i]);
//			std::cout << "size:" << arrayDiff.size() << std::endl;
		}

		std::vector<double> arrayRecover(contrastDiff(arrayDiff, 7));

		std::vector<double> predictRecover(period);

		for(int i = 0; i < period; i++){
			predictRecover[i] = arrayRecover[arrayRecover.size() - period + i];
//			std::cout << predictRecover[i] << ";";
		}
//		for(int i = 0; i < period; i++){
//			std::cout << predictRecover[i] << ";";
//		}
//		std::cout << std::endl;

		for(int i = 0; i < period; i++){
			if(predictRecover[i] > 0){
				weekPredict += predictRecover[i];
			}
		}
//		std::cout << "weekPredict" << weekPredict << std::endl;
		//��ӡ���Ƚ�Ԥ��ֵ����ʵֵ���
		predict[n] = (int)(weekPredict + 0.5);
	}
}