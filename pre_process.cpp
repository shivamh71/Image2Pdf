#include <iostream>
#include <fstream>
#include <cstdio>
using namespace std;

// STL used
#include <vector>

// Image library used
#include "CImg.h"
using namespace cimg_library;

// Struct for representing standard font templates 
struct Template{
	char label;		// label stores character represented by the template 
	int Image[16][16]; // Image represents binary image of corresponding template as 16x16 matrix
	// mean is average of 16x16 image
	// hmean and vmean are means of horizontal and vertical projections respectively
	float mean,hmean,vmean; 
	int horizontalProjection[16],verticalProjection[16]; // arrays for storing horizontal and vertical projections of template
};

// function to convert a grayscale image into a binary image using Balanced Histogram Thresholding
// all pixel values above a certain threshold are set to 255(white) and below threshhold are set to 0(black) 
void binariseImage(CImg<unsigned int>& sourceImage, CImg<unsigned int>& resultImage){
	vector<int> V(256,0);
	vector<long long int> sum(256,0);
	for(int i=0;i<sourceImage.width();i++){
		for (int j=0;j<sourceImage.height();j++){
			V[sourceImage(i,j,0,0)]++;
		}
	}
	for (int i=0;i<255;i++){
		sum[i+1]=sum[i]+V[i+1]*(i+1);
		V[i+1]+=V[i];
	}
	int t=127; // first threshold is assumed to be mid value that is 127
	int low=V[127],high=V[255]-V[127];
	int lowS=sum[127],highS=sum[255]-sum[127];
	int mean;
	if (low==0 || high==0)
		mean=127;
	else
		mean=(lowS/low + highS/high)/2;
	// mean1 = mean of all pixel values higher than threshold t
	// mean2 = mean of all pixel values lower than threshold t
	// now updated threshold = (mean1 + mean2)/2
	// repeat until new threshold and old threshold differ by less than 0.5
	while (t-mean > 0){
		t=mean;
		low=V[t];
		high=V[255]-V[t];
		lowS=sum[t];
		highS=sum[255]-sum[t];
		if (low==0 || high==0)
			break;
		mean=(lowS/low + highS/high)/2;
	}
	cout<<"threshold\t"<<t<<endl;
	for (int i=0;i<sourceImage.width();i++){
		for (int j=0;j<sourceImage.height();j++){
			if (sourceImage(i,j,0,0)>t)
				resultImage(i,j,0,0)=255;
			else
				resultImage(i,j,0,0)=0;
		}
	}
}

// function for generating a grayscale image 
void grayscale(CImg<unsigned int>& sourceImage, CImg<unsigned int>& grayscaleImage){
	int height=sourceImage.height();
	int width=sourceImage.width();
	for (int i=0;i<width;i++){
		for (int j=0;j<height;j++){
			int r=sourceImage(i,j,0,0);
			int g=sourceImage(i,j,0,1);
			int b=sourceImage(i,j,0,2);
			// formula for calculating grayscale value corresponding to rgb values using wighted sum method
			grayscaleImage(i,j,0,0)=(int)(0.2126*r + 0.7152*g + 0.0722*b);
		}
	}
}

// main file
int main(int argc,char* argv[]){
	// creating a binary file to store information of processed standard font templates
	FILE* fp;
	fp=fopen("char_dict.bin","wb");

	// for each argument of command line reading the image
	// image is converted to grayscale and then binarised 
	// after that extra white pixels are removed from all sides 
	// image is rescaled to 16x16
	// all means are calculated and filename is read from path of the image
	// all this information is writen in binary file
	for (int l=1;argv[l]!=NULL;l++){
		CImg<unsigned int> Image(argv[l]);
		CImg<unsigned int> gImage(Image.width(),Image.height(),1,1,0);
		CImg<unsigned int> Image2(Image.width(),Image.height(),1,1,0);
		if (Image.spectrum()==3)
			grayscale(Image,gImage);
		else
			gImage=Image;
		binariseImage(gImage,Image2);
		int x1=-1,x2=-1,y1=-1,y2=-1;
		for (int i=0;i<Image2.height();i++){
			for (int j=0;j<Image2.width();j++){
				if (Image2(j,i,0,0)==0){
					x1=i;
					break;
				}
			}
			if (x1!=-1)
				break;
		}
		for (int i=Image2.height()-1;i>=0;i--){
			for (int j=0;j<Image2.width();j++){
				if (Image2(j,i,0,0)==0){
					x2=i;
					break;
				}
			}
			if (x2!=-1)
				break;
		}
		for (int i=0;i<Image2.width();i++){
			for (int j=0;j<Image2.height();j++){
				if (Image2(i,j,0,0)==0){
					y1=i;
					break;
				}
			}
			if (y1!=-1)
				break;
		}
		for (int i=Image2.width()-1;i>=0;i--){
			for (int j=0;j<Image2.height();j++){
				if (Image2(i,j,0,0)==0){
					y2=i;
					break;
				}
			}
			if (y2!=-1)
				break;
		}
		Image2.crop(y1,x1,y2,x2);
		Image2=Image2.resize(16,16);
		Template T;
		int sum=0;
		for(int k=0;k<16;k++){
			T.horizontalProjection[k]=0;
			T.verticalProjection[k]=0;
		}
		for (int k=0;k<16;k++){
			for (int j=0;j<16;j++){
				if (Image2(j,k,0,0)==255){
					T.Image[k][j]=1;
					T.horizontalProjection[k]++;
					T.verticalProjection[j]++;
					sum++;
				}
				else{
					T.Image[k][j]=0;
				}
				cout<<T.Image[k][j]<<" ";
			}
			cout<<endl;
		}
		int sum1=0,sum2=0;
		for(int k=0;k<16;k++){
			sum1+=T.horizontalProjection[k];
			sum2+=T.verticalProjection[k];
			cout<<T.horizontalProjection[k]<<"  "<<T.verticalProjection[k]<<endl;
		}
		T.mean=((float)sum/256);
		T.hmean=((float)sum1/16);
		T.vmean=((float)sum2/16);
		int n=0;
		char *ch;
		for(int j=0;argv[l][j]!='j';j++){
			if(argv[l][j]<='9' && argv[l][j]>='0'){
				n=n*10+argv[l][j]-'0';
			}
		}
		T.label=n;
		cout<<T.label<<" "<<T.mean<<" "<<T.hmean<<" "<<T.vmean<<endl;
		fwrite(&T,sizeof(T),1,fp);
	}
	fclose(fp);	
	return 0;
}