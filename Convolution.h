// header files
#include <iostream>
using namespace std;

// STL used
#include <vector>

// image library used
#include "CImg.h"
using namespace cimg_library;

// class definition
class Convolution{
private:
	int rows,columns;
	int divisor;
	vector< vector<int> > V;
public:
	Convolution(int x,int y,int d,vector< vector<int> > v){
		rows=x;
		columns=y;
		divisor=d;
		V=v;
	}

	// applying convolution matrix on source Image
	void operate(CImg<unsigned int>& sourceImage,CImg<unsigned int>& resultImage){ // 
		for (int i=rows/2;i<sourceImage.width()-rows/2;i++){
			for (int j=columns/2;j<sourceImage.height()-columns/2;j++){
				int sum=0; 
				for (int k=-rows/2;k<=rows/2;k++){
					for (int l=-columns/2;l<=columns/2;l++){
						sum+=sourceImage(i+k,j+l,0,0)*V[k+rows/2][l+columns/2];
					}
				}
				sum/=divisor;
				resultImage(i,j,0,0)=sum;
			}
		}
	}
};
