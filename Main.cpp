#include <iostream>
#include <fstream>
#include <vector> 
#include <string> 
#include <Eigen/core>
#include <Eigen/dense>

using namespace Eigen;
using namespace std;

double CalculateMean(vector<vector<double>>::iterator AssetReturns)
{
	double sum = 0;
	int numberOfDataPoints = 0;
	for (vector<double>::iterator aReturn = AssetReturns->begin(); aReturn != AssetReturns->end(); aReturn++)
	{
		sum = sum + *aReturn;
		numberOfDataPoints++;
	}
	return (sum / numberOfDataPoints);
}

double CalculateSingleCovariance(vector<vector<double>>::iterator AssetReturn1, vector<vector<double>>::iterator AssetReturn2)
{
	double covariance = 0;
	int numberOfLines = 0;
	double aMean1 = CalculateMean(AssetReturn1);
	double aMean2 = CalculateMean(AssetReturn2);
	vector<double>::iterator aReturn2 = AssetReturn2->begin();

	for (vector<double>::iterator aReturn1 = AssetReturn1->begin(); aReturn1 != AssetReturn1->end(); aReturn1++)
	{
		covariance = covariance + (*aReturn1 - aMean1) * (*aReturn2 - aMean2);
		aReturn2++;
		numberOfLines++;
	}
	covariance = covariance / (numberOfLines - 1);
	return covariance;
}

void MeanVariance(int numberOfAssets, vector<vector<double>> VectorOfAssetReturns)
{
	// Fill Covariance Matrix 
	int ColumnIndex = 0;
	MatrixXd covMatrix(numberOfAssets, numberOfAssets);
	for (vector<vector<double>>::iterator RcolumnForReturnCov = VectorOfAssetReturns.begin(); RcolumnForReturnCov != VectorOfAssetReturns.end(); RcolumnForReturnCov++)
	{
		int RowIndex = 0;
		for (vector<vector<double>>::iterator ScolumnForReturnCov = VectorOfAssetReturns.begin(); ScolumnForReturnCov != VectorOfAssetReturns.end(); ScolumnForReturnCov++)
		{
			double cov = CalculateSingleCovariance(RcolumnForReturnCov, ScolumnForReturnCov);
			covMatrix(RowIndex, ColumnIndex) = cov;
			RowIndex++;
		}
		ColumnIndex++;
	}

	// Fill vector with means  
	RowVectorXd matMeans(numberOfAssets);
	int i = 0;
	for (vector<vector<double>>::iterator aColumn = VectorOfAssetReturns.begin(); aColumn != VectorOfAssetReturns.end(); aColumn++)
	{
		matMeans(i) = CalculateMean(aColumn);
		i++;
	}

	// Create Ones Vector 
	RowVectorXd onesT(numberOfAssets);
	for (int i = 0; i < numberOfAssets; i++)
	{
		onesT(i) = 1;
	}

	cout << "Enter risk free rate as decimal: ";
	double riskfreeRate;
	cin >> riskfreeRate;

	// Calculate Mean Variance Frontier 
	double A = (onesT).dot(covMatrix.inverse()*matMeans.transpose());
	double B = (matMeans*covMatrix.inverse()).dot(matMeans.transpose());
	double C = (onesT).dot(covMatrix.inverse()*onesT.transpose());
	Matrix2d bigA;
	bigA << B, A, A, C;
	double D = bigA.determinant();

	VectorXd g = (1 / D)*(B*onesT - A*matMeans)*covMatrix.inverse();
	VectorXd h = (1 / D)*(C*matMeans - A*onesT)*covMatrix.inverse();

	VectorXd w = g + h*riskfreeRate;

	//Outputs
	cout << endl;
	cout << "MEAN VARIANCE!!!!" << endl; 
	cout << "Means vector: " << endl << matMeans << endl;
	cout << "Covariance Matrix: " << endl << covMatrix << endl;
	cout << "A is " << endl << A << endl;
	cout << "B is " << endl << B << endl;
	cout << "C is " << endl << C << endl;
	cout << "big A is " << endl << bigA << endl;
	cout << "g is " << endl << g << endl;
	cout << "h is " << endl << h << endl;
	
	cout << "Mean Variance Portfolio Weights are " << endl << w << endl << endl;


}

void BlackLitterMan_Org(int numberOfAssets, vector < vector<double >> VectorOfAssetReturns)
{

	// Fill Covariance Matrix 
	int ColumnIndex = 0;
	MatrixXd covMatrix(numberOfAssets, numberOfAssets);
	for (vector<vector<double>>::iterator RcolumnForReturnCov = VectorOfAssetReturns.begin(); RcolumnForReturnCov != VectorOfAssetReturns.end(); RcolumnForReturnCov++)
	{
		int RowIndex = 0;
		for (vector<vector<double>>::iterator ScolumnForReturnCov = VectorOfAssetReturns.begin(); ScolumnForReturnCov != VectorOfAssetReturns.end(); ScolumnForReturnCov++)
		{
			double cov = CalculateSingleCovariance(RcolumnForReturnCov, ScolumnForReturnCov);
			covMatrix(RowIndex, ColumnIndex) = cov;
			RowIndex++;
		}
		ColumnIndex++;
	}

	// Fill vector with means  
	RowVectorXd matMeans(numberOfAssets);
	int i = 0;
	for (vector<vector<double>>::iterator aColumn = VectorOfAssetReturns.begin(); aColumn != VectorOfAssetReturns.end(); aColumn++)
	{
		matMeans(i) = CalculateMean(aColumn);
		i++;
	}

	// Risk Aversion Parameter
	double delta = 1;

	// Tau - uncertainty of CAPM prior 
	double tau = 0.001;
	//cout << "Enter Tau: ";
	//cin >> tau;
	//cout << endl;

	// Must have multiple views (more than 1 ) 
	int numberOfViews_k;
	cout << "Enter number of views (more than 2): (2) " << endl; 
	cin >> numberOfViews_k; 

	RowVectorXd W_eq_EquilibriumPortfolioWeights(numberOfAssets);  
	for (int i = 0; i < numberOfAssets; i++)
	{
		W_eq_EquilibriumPortfolioWeights(i) = (1.0 / numberOfAssets);
		//W_eq_EquilibriumPortfolioWeights(i) = 1;
	}

	VectorXd Pi_EquilibruimExpectedReturn = delta* covMatrix*W_eq_EquilibriumPortfolioWeights.transpose();

	ifstream Views; 
	Views.open("Views.txt");
	MatrixXd P_Views(numberOfViews_k, numberOfAssets);
	double tempView; 
	while (!Views.eof())
	{
		for (int rows = 0; rows < numberOfViews_k; rows++)
		{
			for (int cols = 0; cols < numberOfAssets; cols++)
			{
				Views >> tempView;
				P_Views(rows,cols)= tempView;
			}
		}

	}
	Views.close(); 

	VectorXd Q_ExcessReturns(numberOfAssets);
	Q_ExcessReturns = P_Views*matMeans.transpose();

	// Create 0s for Omega matrix 
	MatrixXd Omega(numberOfViews_k, numberOfViews_k); // Uncertainity within each view matrix 
	for (int rows = 0; rows < numberOfViews_k; rows++)
	{
		for (int col = 0; col < numberOfViews_k; col++)
		{
			Omega(rows, col) = 0;
		}
	}


	// Put values on diagonal 
	for (int i = 0; i < numberOfViews_k; i++)
	{
		double SinlgeOmegaEntry = P_Views.row(i)*covMatrix*P_Views.row(i).transpose();
		SinlgeOmegaEntry = SinlgeOmegaEntry*tau;
		Omega(i, i) = SinlgeOmegaEntry;
	}
	//cout << endl << Omega << endl<<endl;

	MatrixXd mu1(numberOfAssets, numberOfAssets);
	mu1 = (tau*covMatrix.inverse() + P_Views.transpose()*Omega.inverse()*P_Views).inverse();
	VectorXd mu2(numberOfAssets);
	mu2 = (tau*covMatrix).inverse()*Pi_EquilibruimExpectedReturn + P_Views.transpose()*Omega.inverse()*Q_ExcessReturns;

	RowVectorXd mu(numberOfAssets);
	mu = mu1*mu2;


	MatrixXd pi1_hat(numberOfAssets, numberOfViews_k);
	pi1_hat = (tau*covMatrix*P_Views.transpose())*((P_Views*tau*covMatrix*P_Views.transpose() + Omega).inverse());
	VectorXd pi2_hat(numberOfViews_k);
	pi2_hat = Q_ExcessReturns - P_Views*Pi_EquilibruimExpectedReturn;
	VectorXd pi_hat(numberOfAssets);
	pi_hat = Pi_EquilibruimExpectedReturn + pi1_hat*pi2_hat;


	MatrixXd M_PosteriorVariance(numberOfAssets, numberOfAssets);
	M_PosteriorVariance = tau*covMatrix - (tau*covMatrix*P_Views.transpose())*((P_Views*tau*covMatrix*P_Views.transpose() + Omega).inverse())*(P_Views*tau*covMatrix);

	MatrixXd SigmaP(numberOfAssets, numberOfAssets);
	SigmaP = covMatrix + M_PosteriorVariance;

	VectorXd WeightsW(numberOfAssets);
	WeightsW = pi_hat.transpose()*((delta*SigmaP).inverse());


	// Outputs
	cout << endl << endl << "BLACK LITTERMAN OUTPUT" << endl; 
	cout << "Means vector: " << endl << matMeans << endl;
	cout << "Pi is " << endl << Pi_EquilibruimExpectedReturn << endl;
	cout << "Views Matrix " << endl << P_Views << endl;
	cout << "Excess Returns are:" << endl << Q_ExcessReturns << endl;
	cout << "Mu is " << endl << mu << endl;
	cout << "Pi Hat " << endl << pi_hat << endl;
	cout << "Posterior variance M is " << endl << M_PosteriorVariance << endl;
	cout << "Sigma P is " << endl << SigmaP << endl;
	cout << "Black Litterman Weights are " << endl << WeightsW << endl;


	double sum = 0;
	for (int i = 0; i < numberOfAssets; i++)
	{
		sum = sum + WeightsW(i);
	}
	cout << "The sum of the weights are " << endl << sum << endl;
}

vector<vector<double>> ReadInData(int numberOfAssets, string fileName)
{

	// Read data from file 
	ifstream infile;
	infile.open(fileName);
	infile.ignore(20, '\n'); // Ignore the first line

	vector<vector<double>> VectorOfAssetReturns;
	for (int i = 0; i < numberOfAssets; i++)
	{
		vector<double> Asset;
		VectorOfAssetReturns.push_back(Asset);
	}



	int numberOfLines = 0;
	while (!infile.eof())
	{
		infile.ignore(13, ' '); // Ignore the date 
		numberOfLines++; // Keep track of the number of data points 
		vector<vector<double>>::iterator OneAssetReturns = VectorOfAssetReturns.begin();
		for (int i = 0; i < numberOfAssets; i++)
		{
			double Returns;
			infile >> Returns;
			OneAssetReturns->push_back(Returns); // Add Return to the column of Returns 
			OneAssetReturns++;
		}
		infile.ignore(1, '\n');
	}
	cout << "Number of Data Points " << numberOfLines << endl;

	return VectorOfAssetReturns;
}

vector<vector<double>> CreateReturns(int numberOfAssets, vector<vector<double>> VectorOfPrices)
{
	double aReturn; 
	vector<vector<double>>VectorOfReturns;
	for (int i = 0; i < numberOfAssets; i++)
	{
		vector<double> dummyColumn;
		VectorOfReturns.push_back(dummyColumn);
	}
	vector<vector<double>>::iterator ColumnOfPrices = VectorOfPrices.begin();
	vector<vector<double>>::iterator ColumnOfReturns = VectorOfReturns.begin();
	for (ColumnOfPrices; ColumnOfPrices != VectorOfPrices.end(); ColumnOfPrices++)
	{
		vector<double>::iterator currentPrice = ColumnOfPrices->begin();
		vector<double>::iterator nextPrice = ColumnOfPrices->begin() + 1;
		for (nextPrice; nextPrice != ColumnOfPrices->end(); nextPrice++)
		{
			aReturn = (*nextPrice - *currentPrice) / *currentPrice;
			ColumnOfReturns->push_back(aReturn);
			currentPrice++;
		}
		ColumnOfReturns++;
	}
	return VectorOfReturns;
}

int main()
{
	int numberOfAssets;
	cout << "How many assets? (3) ";
	cin >> numberOfAssets; 

	string fileName; 
	fileName = "DataProyect4.txt"; 

	vector<vector<double>> VectorOfAssetPrices = ReadInData(numberOfAssets, fileName);
	vector<vector<double>> VectorOfAssetReturns = CreateReturns(numberOfAssets, VectorOfAssetPrices);
	MeanVariance(numberOfAssets, VectorOfAssetReturns);
	BlackLitterMan_Org(numberOfAssets, VectorOfAssetReturns);

	return 0;
}