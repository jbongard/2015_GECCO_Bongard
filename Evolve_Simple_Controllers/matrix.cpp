/* ---------------------------------------------------
   FILE:     bodyPlan.cpp
	AUTHOR:   Josh Bongard
	DATE:     October 2, 2000
	FUNCTION: This class contains all information for
				 a single physical segment of an organism
				 in the MathEngine environment.
 -------------------------------------------------- */

#include "math.h"

#ifndef _MATRIX_CPP
#define _MATRIX_CPP

#include "matrix.h"

MATRIX::MATRIX(int ln, int wd, double val) {

	length = ln;
	width = wd;

	vals = new double[length*width];

	for (int i=0;i<length;i++)

		for (int j=0;j<width;j++)

			Set(i,j,val);
}

MATRIX::MATRIX(MATRIX *m) {

	length = m->length;
	width  = m->width;

	vals = new double[length*width];

	for (int i=0;i<length;i++)

		for (int j=0;j<width;j++)

//			Set(i,j,m->Get(i,j));
			vals[i*width+j] = m->vals[i*m->width+j];

}

MATRIX::MATRIX(int ln, int wd, ifstream *inFile, int dummyVariable) {

	length = ln;
	width = wd;

	vals = new double[length*width];

	double temp;

	for (int i=0;i<length;i++)

		for (int j=0;j<width;j++) {

			(*inFile) >> temp;
			Set(i,j,temp);
		}
}

MATRIX::MATRIX(ifstream *inFile) {

	ReadFromFile(inFile);
}

MATRIX::MATRIX(char *fileName) {

	ifstream *inFile = new ifstream(fileName);

	ReadFromFile(inFile);

	inFile->close();
	delete inFile;
	inFile = NULL;
}

MATRIX::~MATRIX(void) {

	delete[] vals;
	vals = NULL;
}

double MATRIX::AbsoluteDifference(MATRIX *m) {

	return( AbsoluteDifference(m,width) );
}

double MATRIX::AbsoluteDifference(MATRIX *m, int startWd, int endWd) {

	if ( (length==m->length) && (width==m->width) ) {

		double sum = 0.0;
		
		for (int i=0;i<length;i++)

			for (int j=startWd;j<=endWd;j++)
				
				sum = sum + fabs( double(Get(i,j)) - double(m->Get(i,j)) );

		return( sum / double(length*(endWd-startWd+1)) );
	}
	else
		return( 0.0 );
}

double MATRIX::AbsoluteDifference(MATRIX *m, int wd) {

	if ( (length==m->length) && (width==m->width) ) {

		double sum = 0.0;
		
		for (int i=0;i<length;i++)

			for (int j=0;j<wd;j++)
				
				sum = sum + fabs( double(Get(i,j)) - double(m->Get(i,j)) );

		return( sum / double(length*wd) );
	}
	else
		return( 0.0 );
}

double MATRIX::AbsoluteDifference(MATRIX *m, int wd, double normFactor) {

	return( AbsoluteDifference(m,wd) * double(length*wd) / normFactor );
}

double MATRIX::AbsoluteDifferenceOfTails(MATRIX *m) {

	if ( length==m->length ) {

		double sum = 0.0;
		
		for (int i=0;i<length;i++)
				
			sum = sum + fabs( double(Get(i,width-1)) - double(m->Get(i,m->width-1)) );

		return( sum / double(length) );
	}
	else
		return( 0.0 );
}

void MATRIX::Add(int i, int j, double val) {

	Set(i,j,Get(i,j)+val);

}

void MATRIX::AddToColumn(int j, double val) {

	for (int i=0;i<length;i++)
		Add(i,j,val);
}

void MATRIX::AverageMatrices(MATRIX **matrices, int numMatrices) {

	double total;
	double mean;

	for (int i=0;i<length;i++)
		
		for (int j=0;j<width;j++) {

			total = 0.0;

			for (int k=0;k<numMatrices;k++)

				total = total + matrices[k]->Get(i,j);

			mean = total / double(numMatrices);
			Set(i,j,mean);
		}
}

int  MATRIX::BinaryToDecimal(int i, int j, int ln) {

	int totalValue = 0;
	int placeValue = 128;
	int temp;

	for (int wd=j;wd<(j+ln);wd++) {

		temp = int(Get(i,wd));

		totalValue = totalValue + placeValue*int(Get(i,wd));
		placeValue = placeValue/2;
	}

	return( totalValue );
}

MATRIX *MATRIX::ColMeans(int lastRow) {

	MATRIX *colMeans = new MATRIX(1,width,0.0);

	for (int j=0;j<width;j++)

		for (int i=0;i<=lastRow;i++)
	
			colMeans->Add(0,j,Get(i,j));

	colMeans->Div(double(lastRow+1));

	return( colMeans );
}

int  MATRIX::ComputeFSM(MATRIX *str, MATRIX *states, int i, int strIndex, int strLength, int q0, MATRIX *F) {

	states->Set(i,strIndex,q0);

	if ( strIndex == strLength ) {

		return( int(F->Get(0,q0)) );
	}
	else {

		int currChar = int(str->Get(i,strIndex));

		return( ComputeFSM(str,states,i,++strIndex,strLength,int(Get(currChar,q0)),F) );
	}
}

void MATRIX::CopyFrom(MATRIX *other) {

	for (int i=0;i<length;i++)

		for (int j=0;j<width;j++)

			vals[i*width + j] = other->vals[i*width + j];
			//Set(i,j,other->Get(i,j));
}

void MATRIX::CopyRow(int sourceRow, int targetRow) {

        for (int j=0;j<width;j++)

                Set(targetRow,j,Get(sourceRow,j));
}

void MATRIX::CopyRow(int myRow, MATRIX *other, int hisRow) {

	for (int j=0;j<width;j++) {

//		double val = other->Get(hisRow,j);
		double val = other->vals[hisRow*other->width+j];

		vals[myRow*width+j] = val;
//		Set(myRow,j,val);
	}
}

MATRIX *MATRIX::CountValues(int i1, int i2, int j1, int j2, int maxVal) {

	MATRIX *valueCount = new MATRIX(1,maxVal,0);

	for (int i=i1;i<i2;i++)

		for (int j=j1;j<j2;j++)

			valueCount->Add(0,int(Get(i,j)),1);

	return( valueCount );
}

void MATRIX::CreateChain(int depth) {

	int i;
	int depthCounter;

	for (int j=0;j<width;j++) {

		depthCounter = depth;
		i = j;

		while ( depthCounter > 0 ) {

			Set(i,j,1);
			i = i + 1;
			if ( i==length )
				i = 0;
			depthCounter--;
		}
	}
}

void MATRIX::CreateIdentity(void) {

	for (int i=0;i<length;i++)
		for (int j=0;j<width;j++)
			if ( i == j )
				Set(i,j,1);
			else
				Set(i,j,0);
}

void MATRIX::CreateParity(void) {

	MATRIX *binaryOfColumn;
	int binarySum;

	for (int j=0;j<width;j++)
		
		for (int i=0;i<length;i++) {
			
			binaryOfColumn = DecimalToBinary(i,length-1);

			binarySum = int(binaryOfColumn->SumOfRow(0));
			
			if ( (binarySum%2) == 0 )
				Set(i,j,0);
			else
				Set(i,j,1);

			delete binaryOfColumn;
			binaryOfColumn = NULL;

		}
}

double MATRIX::DistBetRows(int i1, int i2) {

	double dist = 0.0;

	for (int j=0;j<width;j++)
		dist = dist + pow(Get(i1,j)-Get(i2,j),2.0);

	dist = sqrt(dist);

	return(dist);
}

void MATRIX::Div(double val) {

	for (int i=0;i<length;i++)
		for (int j=0;j<width;j++)
			Div(i,j,val);
}

void MATRIX::Div(int i, int j, double val) {

	Set(i,j,Get(i,j)/val);
}

double MATRIX::EqualColumnVals(int col1, int col2) {

	double sumOfEquals = 0.0;

	for (int i=0;i<length;i++)

		if ( Get(i,col1) == Get(i,col2) )
			sumOfEquals++;

	return( sumOfEquals / double(length) );
}

void MATRIX::FillColumn(int myCol, int hisCol, MATRIX *m) {

	for (int i=0;i<length;i++)
		Set(i,myCol,m->Get(i,hisCol));
}

void MATRIX::FillRow(int myRow, int hisRow, MATRIX *m) {

	for (int j=0;j<width;j++)
		Set(myRow,j,m->Get(hisRow,j));
}

void MATRIX::FindMinMaxInColumn(int j, int i1, int i2, double *min, double *max) {

	(*min) = 1000.0;
	(*max) = -1000.0;

	for (int i=i1;i<=i2;i++) {

		if ( Get(i,j) < (*min) )
			(*min) = Get(i,j);

		if ( Get(i,j) > (*max) )
			(*max) = Get(i,j);
	}
}

void MATRIX::Flip(int i, int j) {

	if ( Get(i,j) )
		Set(i,j,0);
	else
		Set(i,j,1);
}

double MATRIX::Get(int i, int j) {

	return( vals[i*width + j] );
}

double MATRIX::GetDistanceReading(void) {

	double dist = Get(length-1,14);

	if ( dist == 0 )
		dist = 0.00001;

	return( dist );
}

double MATRIX::GetForwardBackTiltReading(void) {

	return( Get(length-1,13) );
}

double MATRIX::GetLeastReadingFromPressureSensors(void) {

	double minVal = 1000.0;

	for (int j=8;j<12;j++) {

		if ( Get(length-1,j) < minVal )
			minVal = Get(length-1,j);
	}

	return( minVal );
}

double MATRIX::GetLeftRightTiltReading(void) {

	return( Get(length-1,12) );
}

double MATRIX::GetMaxReadingFromPressureSensors(void) {

	double maxVal = -1000.0;

	for (int j=8;j<12;j++) {

		if ( Get(length-1,j) > maxVal )
			maxVal = Get(length-1,j);
	}

	return( maxVal );
}

void MATRIX::GetMaxValsInColumn(int j, MATRIX *maxVals) {

	double maxVal = -1000.0;

	int i;
	int k;

	for (i=0;i<length;i++)

		if ( Get(i,j) > maxVal )
			maxVal = Get(i,j);

	maxVals->Set(0,0,maxVal);

	for (k=1;k<maxVals->width;k++) {

		maxVal = -1000.0;

		for (i=0;i<length;i++) {

			if ( (Get(i,j) > maxVal) && !(maxVals->In(Get(i,j))) )
				maxVal = Get(i,j);
		}

		maxVals->Set(0,k,maxVal);
	}
}

void MATRIX::GetMaxValsInColumn(int j, int iFirst, int iLast, MATRIX *maxVals) {

	double maxVal = -1000.0;

	int i;
	int k;

	for (i=iFirst;i<=iLast;i++)

		if ( Get(i,j) > maxVal )
			maxVal = Get(i,j);

	maxVals->Set(0,0,maxVal);

	for (k=1;k<maxVals->width;k++) {

		maxVal = -1000.0;

		for (i=iFirst;i<=iLast;i++) {

			if ( (Get(i,j) > maxVal) && !(maxVals->In(Get(i,j))) )
				maxVal = Get(i,j);
		}

		maxVals->Set(0,k,maxVal);
	}
}

double MATRIX::GetTotalReadingsFromPressureSensors(void) {

	double totalReadings = 0.0;

	for (int j=8;j<12;j++) 
			totalReadings = totalReadings + Get(length-1,j);

	return( totalReadings );
}

double MATRIX::GetTotalSpeed(void) {

	double dist = Get(length-1,15);

	return( dist );
}

int  MATRIX::In(double val) {

	int found = false;

	int i=0;
	int j;

	while ( (i<length) && (!found) ) {

		j = 0;

		while ( (j<width) && (!found) ) {

			if ( Get(i,j) == val )
				found = true;

			j++;
		}

		i++;
	}

	return( found );
}

void MATRIX::Load(ifstream *inFile) {

	for (int i=0;i<length;i++)

		for (int j=0;j<width;j++)

			(*inFile) >> vals[i*width+j];
}

double MATRIX::Max(void) {

	double max = -1000000.0;

	for (int i=0;i<length;i++)

		for (int j=0;j<width;j++)

			if ( Get(i,j) > max )

				max = Get(i,j);

	return( max );
}

double MATRIX::MaxDiff(MATRIX *otherM) {

	double maxDiff = -1.0;

	for (int i=0;i<length;i++)

		for (int j=0;j<width;j++)

			if ( fabs(Get(i,j)-otherM->Get(i,j)) > maxDiff )
				maxDiff = fabs(Get(i,j)-otherM->Get(i,j));

	return( maxDiff );
}

double MATRIX::MaxDistBetRows(MATRIX *otherM) {

	double maxDist = -1000.0;

	for (int i=0;i<length;i++) {

		double dist = 0.0;

		for (int j=0;j<width;j++) {

			dist = dist + pow( Get(i,j) - otherM->Get(i,j) , 2.0);
		}

		dist = sqrt(dist);

		if ( dist > maxDist )
			maxDist = dist;
	}

	return( maxDist );
}

double MATRIX::MaxValInColumn(int j) {

	double maxVal = -1000.0;

	for (int i=0;i<length;i++)

		if ( Get(i,j) > maxVal )

			maxVal = Get(i,j);

	return( maxVal );
}

double MATRIX::MaxValInColumnContingentOn(int j, int k) {

	double maxVal = -1000.0;

	for (int i=0;i<length;i++)

		if ( (Get(i,j) > maxVal) && (Get(i,k) == 0) )

			maxVal = Get(i,j);

	return( maxVal );
}

double MATRIX::Mean(void) {

	double mean = 0.0;

	for (int i=0;i<length;i++)

		for (int j=0;j<width;j++)

			mean = mean + Get(i,j);

	return( mean / double(length*width) );
}

double MATRIX::Mean(int i1, int i2, int j1, int j2) {

	double mean = 0.0;

	for (int i=i1;i<=i2;i++)

		for (int j=j1;j<=j2;j++)

			mean = mean + Get(i,j);

	return( mean / double((i2-i1+1)*(j2-j1+1)) );
}

double MATRIX::MeanDist(MATRIX *otherM) {

	double dist = 0.0;

	for (int i=0;i<length;i++)

		for (int j=0;j<width;j++)

			dist = dist + fabs( Get(i,j) - otherM->Get(i,j) );

	return( dist/double(length*width) );
}

double MATRIX::MeanDist(MATRIX *otherM, int i1, int i2, int j1, int j2) {

	double dist = 0.0;

	for (int i=i1;i<=i2;i++)

		for (int j=j1;j<=j2;j++)

			dist = dist + fabs( Get(i,j) - otherM->Get(i,j) );

	return( dist/double((i2-i1+1)*(j2-j1+1)) );
}

double MATRIX::MeanDistBetRows(void) {

	if ( length==1 )
		return( 1.0 );
	
	double meanDists = 0.0;
	double numPairs = 0;

	for (int i1=0;i1<length-1;i1++) {

		for (int i2=i1+1;i2<length;i2++) {

			double dist = 0.0;

			for (int j=0;j<width;j++)
				dist = dist + pow( Get(i1,j) - Get(i2,j) ,2.0);

			meanDists = meanDists + sqrt(dist);
			numPairs++;
		}
	}

	meanDists = meanDists / double(numPairs);

	return( meanDists );
}

double MATRIX::MeanDistBetRows(MATRIX *otherM) {

	double meanDists = 0.0;

	for (int i=0;i<length;i++) {

		double dist = 0.0;

		for (int j=0;j<width;j++) {

			dist = dist + pow( Get(i,j) - otherM->Get(i,j) , 2.0);
		}
		meanDists = meanDists + sqrt(dist);
	}

	return( meanDists/double(length) );
}

double MATRIX::Min(void) {

	double min = 1000000.0;

	for (int i=0;i<length;i++)

		for (int j=0;j<width;j++)

			if ( Get(i,j) < min )

				min = Get(i,j);

	return( min );
}

double MATRIX::MinDistBetRows(void) {

	if ( length==1 )
		return( 0.0 );
	
	double minDist = 1000.0;

	for (int i1=0;i1<length-1;i1++) {

		for (int i2=i1+1;i2<length;i2++) {

			double dist = 0.0;

			for (int j=0;j<width;j++)
				dist = dist + pow( Get(i1,j) - Get(i2,j) ,2.0);

			dist = sqrt(dist);

			if ( dist < minDist )
				minDist = dist;
		}
	}

	return( minDist );
}

double MATRIX::MinValInColumn(int j) {

	double minVal = 1000.0;

	for (int i=0;i<length;i++)

		if ( Get(i,j) < minVal )

			minVal = Get(i,j);

	return( minVal );
}

int  MATRIX::MostSimilarRow(MATRIX *r, int i1, int i2) {

	int maxSimilarity = 10000;
	int rowDiff;

	for (int i=i1;i<i2;i++) {

		rowDiff = RowDifference(i,0,r);

		if ( rowDiff < maxSimilarity )
			maxSimilarity = rowDiff;
	}

	return( maxSimilarity );
}

double MATRIX::MSE(MATRIX *other) {

	double total = 0.0;

	for (int i=0;i<length;i++)
		for (int j=0;j<width;j++)

			total = total + pow(Get(i,j) - other->Get(i,j),2.0);

	total = total / (length*width);

	return( total );
}

void MATRIX::Mult(int i, int j, double val) {

	Set(i,j,Get(i,j)*val);
}

void MATRIX::MultColumn(int j, double val) {

	for (int i=0;i<length;i++)
		Set(i,j,Get(i,j)*val);
}

void MATRIX::Normalize(void) {

	double sum = Sum();

	for (int i=0;i<length;i++)

		for (int j=0;j<width;j++)

			Set(i,j,Get(i,j)/sum);
}

void MATRIX::Print(void) {

	Print(3);
}

void MATRIX::Print(int decimalPlaces) {

	for (int i=0;i<length;i++) {

		for (int j=0;j<width;j++) {

			if ( vals[i*width+j]>=0 )
				printf(" ");

			if ( decimalPlaces == 0 )

				printf("%1.0f ",vals[i*width+j]);
			else {
				if ( decimalPlaces ==1 )
					printf("%1.1f ",vals[i*width+j]);
				if ( decimalPlaces ==2 )
					printf("%2.2f ",vals[i*width+j]);
				if ( decimalPlaces ==3 )
					printf("%3.3f ",vals[i*width+j]);
				if ( decimalPlaces ==4 )
					printf("%4.4f ",vals[i*width+j]);

			}

		}
		printf("\n");
	}

//	char ch = getchar();
}

void MATRIX::Print(int i1, int i2, int j1, int j2) {

	for (int i=i1;i<=i2;i++) {

		for (int j=j1;j<=j2;j++) {

				if ( vals[i*width+j]>=0 )
					printf(" ");

				printf("%1.0f ",vals[i*width+j]);
		}
		printf("\n");
	}

//	char ch = getchar();
}

void MATRIX::PrintColumnAsRow(int j) {

	for (int i=0;i<length;i++)
		printf("%3.3f ",Get(i,j));

	printf("\n");
}

void MATRIX::PrintColumnMeans(void) {

	for (int j=0;j<width;j++)

		printf("%3.3f ",Mean(0,length-1,j,j));

	printf("\n");
}

void MATRIX::PrintColumnMeans(int i1, int i2) {

	for (int j=0;j<width;j++)

		printf("%3.3f ",Mean(i1,i2,j,j));

	printf("\n");
}

void MATRIX::PrintDiff(MATRIX *other) {

	for (int i=0;i<length;i++) {

		for (int j=0;j<width;j++) {

			printf("%1.0f ",fabs(Get(i,j)-other->Get(i,j)));
		}
		printf("\n");
	}
}


void MATRIX::PrintRow(int i, int decimalPlaces) {

                for (int j=0;j<width;j++) {

                        if ( vals[i*width+j]>=0 )
                                printf(" ");

                        if ( decimalPlaces == 0 )

                                printf("%1.0f ",vals[i*width+j]);
                        else {
                                if ( decimalPlaces ==1 )
                                        printf("%1.1f ",vals[i*width+j]);
                                if ( decimalPlaces ==2 )
                                        printf("%2.2f ",vals[i*width+j]);
                                if ( decimalPlaces ==3 )
                                        printf("%3.3f ",vals[i*width+j]);
                                if ( decimalPlaces ==4 )
                                        printf("%4.4f ",vals[i*width+j]);
                        }
                }
                printf("\n");
}

void MATRIX::PrintRows(int lastRowIndex) {

	for (int i=0;i<=lastRowIndex;i++) {
		for (int j=0;j<width;j++)
			if ( Get(i,j) > 0 )
				printf("+%1.0f ",Get(i,j));
			else
				printf("%1.0f ",Get(i,j));

		printf("\n");
	}
}

void MATRIX::Reduce(int amt) {

	double newLength = int(double(length) / double(amt));

	double *newVals = new double[int(newLength)*width];

	double total;
	double mean;
	int n;

	for (int j=0;j<width;j++) {

		n = 0;

		for (int i=0;i<length;i=i+amt) {

			total = 0.0;

			for (int k=i;k<i+amt;k++) {

				total = total + Get(k,j);
			}
			mean = total / double(amt);
			
			newVals[n*width + j] = mean;

			n++;
		}
	}

	delete vals;
	vals = newVals;

	length = int(newLength);
}

void MATRIX::Replace(MATRIX *m) {

	if ( (length==m->length) && (width==m->width) ) {

		for (int i=0;i<length;i++)

			for (int j=0;j<width;j++)

				m->Set(i,j,Get(i,j));
	}
}

void MATRIX::ReZero(void) {

	for (int i=0;i<length;i++)

		for (int j=0;j<width;j++)

			Set(i,j,0.0);
}

double MATRIX::RollingMean(MATRIX *other, int h, int w) {

	double rm = 0.0;

	for (int j=0;j<width;j++)
		rm = rm + RollingMean(other,j,h,w);

	return( rm / double(width) );
}

void MATRIX::Round(void) {

	for (int i=0;i<length;i++)

		for (int j=0;j<width;j++)

			Set(i,j,int(Get(i,j)));
}

int  MATRIX::RowDifference(int myRow, int hisRow, MATRIX *m) {

	int diff = 0;

	for (int j=0;j<width;j++)
		diff = diff + abs( int(Get(myRow,j)) - int(m->Get(hisRow,j)) );

	return( diff );
}

void MATRIX::Set(int i, int j, double val) {

	vals[i*width + j] = val;
}

void MATRIX::SetAllTo(double val) {

	for (int i=0;i<length;i++)

		for (int j=0;j<width;j++)

			Set(i,j,val);
}

void MATRIX::SetRow(int i, int val) {

	for (int j=0;j<width;j++)
		Set(i,j,val);
}


void MATRIX::Sqrt(void) {
	
	for (int i=0;i<length;i++)

		for (int j=0;j<width;j++)

			Sqrt(i,j);
}

void MATRIX::Sqrt(int i, int j) {

	Set(i,j,sqrt(Get(i,j)));
}

double MATRIX::Sum(void) {

	double sum = 0.0;

	for (int i=0;i<length;i++)

		for (int j=0;j<width;j++)

			sum = sum + Get(i,j);

	return( sum );
}

double MATRIX::Sum(int i1, int i2, int j1, int j2) {

	double sum = 0.0;

	for (int i=i1;i<=i2;i++)
		for (int j=j1;j<=j2;j++)
			sum = sum + Get(i,j);

	return( sum );
}

int MATRIX::SumOfColumn(int j) {

	int sum = 0;

	for (int i=0;i<length;i++)
		sum = sum + int(Get(i,j));

	return( sum );
}

int MATRIX::SumOfColumn(int i1, int i2, int j) {

	int sum = 0;

	for (int i=i1;i<i2;i++)
		sum = sum + int(Get(i,j));

	return( sum );
}

int MATRIX::SumOfIndices(MATRIX *indices, int i, int j1, int j2) {

	int sum = 0;

	for (int j=j1;j<j2;j++)
		sum = sum + int(Get(0,int(indices->Get(i,j))));

	return( sum );
}

double MATRIX::SumOfRow(int i) {

	double sum = 0;

	for (int j=0;j<width;j++)
		sum = sum + Get(i,j);

	return( sum );
}

double MATRIX::SumOfRow(int i, int j1, int j2) {

	double sum = 0;

	for (int j=j1;j<=j2;j++)
		sum = sum + Get(i,j);

	return( sum );
}

void MATRIX::Write(ofstream *outFile) {

	(*outFile) << length << " " << width << "\n";

	WriteWithoutPreamble(outFile);
}

void MATRIX::WriteToFile(char *fileName, int append) {

	ofstream *outFile;
	
	if ( append )
		outFile = new ofstream(fileName,ios::app);
	else
		outFile = new ofstream(fileName);

	for (int i=0;i<length;i++) {
		for (int j=0;j<width;j++)
			(*outFile) << Get(i,j) << " ";
		(*outFile) << "\n";
	}

	outFile->close();
	delete outFile;
	outFile = NULL;
}

void MATRIX::WriteWithoutPreamble(ofstream *outFile) {

	for (int i=0;i<length;i++) {

		for (int j=0;j<width;j++)

			(*outFile) << Get(i,j) << " ";

		(*outFile) << "\n";
	}
}

void MATRIX::ZeroColumn(int j) {

	for (int i=0;i<length;i++)
		Set(i,j,0);
}

// ----------------------------------------------------------------
//                           Private methods
// ----------------------------------------------------------------

int  MATRIX::BinaryToDecimal(int j) {

	int totalValue = 0;
	int placeValue = 1;

	for (int i=0;i<length;i++) {

		totalValue = totalValue + placeValue*int(Get(i,j));
		placeValue = placeValue*2;
	}

	return( totalValue );
}

int MATRIX::Contains(int val) {

	int found = false;
	int i = 0;
	int j;

	while ( (i<length) && (!found) ) {

		j = 0;

		while ( (j<width) && (!found) ) {

			if ( Get(i,j) == val )
				found = true;

			j++;
		}

		i++;
	}

	return( found );
}

MATRIX *MATRIX::DecimalToBinary(int val, int maxValue) {

	int power = 0;

	while ( (maxValue - pow(2,power)) >= 0 ) {

		power = power + 1;
	}

	MATRIX *binaryVal = new MATRIX(1,power,0);

	while ( val > 0 ) {

		power = 0;

		while ( (val - pow(2,power)) >= 0 ) {

			power = power + 1;
		}

		binaryVal->Set(0,binaryVal->width-power,1);
		val = val - int(pow(2,power-1));
	}

	return( binaryVal );
}

int  MATRIX::FindFirstValue(int j) {

	int found = false;
	int index = 0;

	while ( (index<length) && (!found) ) {

		if ( Get(index,j) )
			found = true;
		else
			index++;
	}

	return( index );
}

double MATRIX::Mean(int i, int startJ, int endJ) {

	double mean = 0.0;

	for (int j=startJ;j<=endJ;j++)
		mean = mean + Get(i,j);

	return( mean/double(endJ-startJ+1) );
}

void MATRIX::ReadFromFile(ifstream *inFile) {

	double temp;

	(*inFile) >> length;
	(*inFile) >> width;

	vals = new double[length*width];

	for (int i=0;i<length;i++)

		for (int j=0;j<width;j++) {

			(*inFile) >> temp;
			Set(i,j,temp);
		}
}

double MATRIX::RollingMean(MATRIX *other, int j, int h, int w) {

	double rm = 0.0;
	
	/*
  	double hisSum;
	double mySum;

	int start = int( double(w-1) / 2.0 );
	int end   = start + h;

	int leftOffset  = -start;
	int rightOffset = start;

	for ( int i=start ; i<end ; i++ ) {

		hisSum = 0.0;
		mySum  = 0.0;

		for ( int k=(i+leftOffset) ; k<=(i+rightOffset) ; k++ ) {

			hisSum = hisSum + other->Get(k,j);
			mySum  = mySum  + Get(k,j);
		}

		hisSum = hisSum / double(w);
		mySum  = mySum / double(w);

		rm = rm + fabs( hisSum - mySum );

	}

  	return( rm / double(h) );

	*/

	for (int i=0;i<length;i++) {

		rm = rm + fabs( Get(i,j) - other->Get(i,j) );
	}

	return( rm / double(length) );
}

#endif
