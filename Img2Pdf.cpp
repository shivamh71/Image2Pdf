// Headerfiles Included
#include <iostream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <utility>
#include <climits>
#include <cfloat>
#include "Convolution.h"

// STL used
#include <vector>
#include <queue>
#include <list>
using namespace std;

// Graphics Libraries used
#include <gtk/gtk.h>

// Image Library used
#include "CImg.h"
using namespace cimg_library;

// Main Code

// Struct for representing boxes that enclose individual characters
struct charBox{ 
	// (x1,y1) represents co-ordinates of top left corner
	// (x2,y2) represents co-ordinates of bottom right corner
	int x1,y1,x2,y2;  
};

// Struct for representing standard font templates 
struct Template{
	char label;		// label stores character represented by the template 
	int Image[16][16]; // Image represents binary image of corresponding template as 16x16 matrix
	// mean is average of 16x16 image
	// hmean and vmean are means of horizontal and vertical projections respectively
	float mean,hmean,vmean; 
	int horizontalProjection[16],verticalProjection[16]; // arrays for storing horizontal and vertical projections of template
};

// Struct for storing correlation value and corresponding character as a pair
struct Pair{
	float first; 
	char second;
};

vector<Template> charTemplates; // vector of all standard font templates

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

// funtion for taking mod of two images 
// modImage[i][j]= sqrt (Image1[i][j]^2 + Image2[i][j]^2)
void modImage(CImg<unsigned int>& sourceImage1,CImg<unsigned int>& sourceImage2,CImg<unsigned int>& resultImage){
	for (int i=0;i<sourceImage1.width();i++){
		for (int j=0;j<sourceImage1.height();j++){
			resultImage(i,j,0,0)=sqrt(sourceImage1(i,j,0,0)*sourceImage1(i,j,0,0) + sourceImage2(i,j,0,0)*sourceImage2(i,j,0,0));
		}
	}
}

// function for getting an eroded image of input sourceImage
// erosion of image enlarges small patches present in image within connected components
void Erosion(CImg<unsigned int>& sourceImage, CImg<unsigned int>& resultImage){
	vector< vector<int> > se={{1,1,1},{1,1,1},{1,1,1}}; // structuring element used for erosion
	for (int i=1;i<sourceImage.width()-1;i++){
		for (int j=1;j<sourceImage.height()-1;j++){
			int min=sourceImage(i,j,0,0);
			for (int k=-1;k<=1;k++){
				for (int l=-1;l<=1;l++){
					if (i+k>=0 && i+k<sourceImage.width() && j+l>=0 && j+l<sourceImage.height()){
						if (sourceImage(i+k,j+l,0,0)-se[k+1][l+1]<min)
							min=sourceImage(i+k,j+l,0,0)-se[k+1][l+1];
					}	
				}
			}
			if (min<0)
				min=0;
			resultImage(i,j,0,0)=min;
		}
	}
}

// function for dilating input source image
// erosion of image fills up the small patches present in image within connected components
void Dilation(CImg<unsigned int>& sourceImage, CImg<unsigned int>& resultImage){
	vector< vector<int> > se={{1,1,1},{1,1,1},{1,1,1}}; // structuring element used for dilation operation
	for (int i=1;i<sourceImage.width()-1;i++){
		for (int j=1;j<sourceImage.height()-1;j++){
			int max=sourceImage(i,j,0,0);
			for (int k=-1;k<=1;k++){
				for (int l=-1;l<=1;l++){
					if (i+k>=0 && i+k<sourceImage.width() && j+l>=0 && j+l<sourceImage.height()){
						if (sourceImage(i+k,j+l,0,0)+se[k+1][l+1]>max)
							max=sourceImage(i+k,j+l,0,0)+se[k+1][l+1];
					}	
				}
			}
			if (max>255)
				max=255;
			resultImage(i,j,0,0)=max;
		}
	}
}

// function to generate an image that is simply average of two images (all images are grayscale)
// averageImage[i][j]= (Image1[i][j] + Image2[i][j]) / 2
void avgImage(CImg<unsigned int>& sourceImage1,CImg<unsigned int>& sourceImage2,CImg<unsigned int>& resultImage){
	for (int i=0;i<sourceImage1.width();i++){
		for (int j=0;j<sourceImage1.height();j++){
			resultImage(i,j,0,0)=(sourceImage1(i,j,0,0)+sourceImage2(i,j,0,0))/2;
		}
	}
}

// function to convert a grayscale image into a binary image using Balanced Histogram Thresholding
// all pixel values above a certain threshold are set to 255(white) and below threshhold are set to 0(black) 
void binariseImage(CImg<unsigned int>& sourceImage, CImg<unsigned int>& resultImage){
	vector<float> V(256,0);
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
	int lowS=sum[127],highS=sum[255]-sum[127],mean=(lowS/low + highS/high)/2;
	// mean1 = mean of all pixel values higher than threshold t
	// mean2 = mean of all pixel values lower than threshold t
	// now updated threshold = (mean1 + mean2)/2
	// repeat until new threshold and old threshold differ by less than 0.5
	while (t-mean > -0.5 && t-mean < 0.5){
		t=mean;
		low=V[t];
		high=V[255]-V[t];
		lowS=sum[t];
		highS=sum[255]-sum[t];
		mean=(lowS/low + highS/high)/2;
	}
	cout<<t<<endl;
	for (int i=0;i<sourceImage.width();i++){
		for (int j=0;j<sourceImage.height();j++){
			if (sourceImage(i,j,0,0)>t)
				resultImage(i,j,0,0)=255;
			else
				resultImage(i,j,0,0)=0;
		}
	}
}

vector< vector<int> > label; // stores pixel values with its connected component's label
vector<list<int> > charLabel; // stores list containing indices of characters corresponding to a connected component label
vector< vector<int> > label1; // stores pixel values with its text block's label
vector<list<int> > charLabel1; // stores list containing indices of characters corresponding to a text block label

// given an input image it generates vector containing co-ordinates separated character boxes
vector<vector<int> > extractCharacters(CImg<unsigned int> sourceImage,CImg<unsigned int> sourceImage1,CImg<unsigned int> sourceImage2,CImg<unsigned int> Image){
	// initialisations
	label.resize(sourceImage.height()); 
	for (int i=0;i<label.size();i++){
		vector<int> v(sourceImage.width(),0);
		label[i]=v;
	}
	label1.resize(sourceImage.height());
	for (int i=0;i<label.size();i++){
		vector<int> v(sourceImage.width(),0);
		label1[i]=v;
	}
	vector< list<pair<int,int> > > equiLabel;
	int curLabel=1;
	list <pair<int,int> > ls;
	equiLabel.push_back(ls);

	// determining connected components
	for (int i=0;i<sourceImage.height();i++){
		for (int j=0;j<sourceImage.width();j++){
			if (sourceImage(j,i,0,0)==255){
				int min=INT_MAX;
				if (j-1>=0 && sourceImage(j-1,i,0,0)==255){
					int l=label[i][j-1];
					if(min>l && l!=0){
						min=l;
						equiLabel[l].push_back(make_pair(i,j));
						label[i][j]=l;
					}
				}
				if (j-1>=0 && i-1>=0 && sourceImage(j-1,i-1,0,0)==255){
					int l=label[i-1][j-1];
					if (min==INT_MAX){
						min=l;
						equiLabel[l].push_back(make_pair(i,j));
						label[i][j]=l;
					}
					else if (min>l && l!=0){
						for (auto it=equiLabel[min].begin();it!=equiLabel[min].end();it++){
							label[it->first][it->second]=l;
							equiLabel[l].push_back(*it);
						}
						equiLabel[min].clear();
						min=l;
					}
					else if (min<l){
						for (auto it=equiLabel[l].begin();it!=equiLabel[l].end();it++){
							label[it->first][it->second]=min;
							equiLabel[min].push_back(*it);
						}
						equiLabel[l].clear();
					}
				}
				if (i-1>=0 && sourceImage(j,i-1,0,0)==255){

					int l=label[i-1][j];
					if (min==INT_MAX){
						min=l;
						equiLabel[l].push_back(make_pair(i,j));
						label[i][j]=l;
					}
					else if (min>l && l!=0){
						for (auto it=equiLabel[min].begin();it!=equiLabel[min].end();it++){
							label[it->first][it->second]=l;
							equiLabel[l].push_back(*it);
						}
						equiLabel[min].clear();
						min=l;
					}
					else if (min<l){
						for (auto it=equiLabel[l].begin();it!=equiLabel[l].end();it++){
							label[it->first][it->second]=min;
							equiLabel[min].push_back(*it);
						}
						equiLabel[l].clear();
					}
				}
				if (i-1>=0 && j+1<sourceImage.width() && sourceImage(j+1,i-1,0,0)==255){
					int l=label[i-1][j+1];
					if (min==INT_MAX){
						min=l;
						equiLabel[l].push_back(make_pair(i,j));
						label[i][j]=l;
					}
					else if (min>l && l!=0){
						for (auto it=equiLabel[min].begin();it!=equiLabel[min].end();it++){
							label[it->first][it->second]=l;
							equiLabel[l].push_back(*it);
						}
						equiLabel[min].clear();
						min=l;
					}
					else if (min<l){
						for (auto it=equiLabel[l].begin();it!=equiLabel[l].end();it++){
							label[it->first][it->second]=min;
							equiLabel[min].push_back(*it);
						}
						equiLabel[l].clear();
					}
				}
				if(min==INT_MAX){
					list<pair<int,int> > lst;
					lst.push_back(make_pair(i,j));
					label[i][j]=curLabel;
					equiLabel.push_back(lst);
					curLabel++;
				}
			}
			else{
				label[i][j]=0;
			}
		}
	}

	// determing rectangle enclosing each of the connected components
	// boxVector consists of { minX, minY, maxX, maxY }
	// these are co-ordinates of rectangle enclosing single connected component
	vector< vector<int> > boxVector(curLabel,vector<int>(4,0));
	for (int i=0;i<sourceImage.height();i++){
		for (int j=0;j<sourceImage.width();j++){
			int l=label[i][j]; 
			if (i<boxVector[l][0] || boxVector[l][0]==0){
				boxVector[l][0]=i;
			}
			if (j<boxVector[l][1] || boxVector[l][1]==0){
				boxVector[l][1]=j;
			}
			if (i>boxVector[l][2] || boxVector[l][2]==0){
				boxVector[l][2]=i;
			}
			if (j>boxVector[l][3] || boxVector[l][3]==0){
				boxVector[l][3]=j;
			}
		}
	}

	// drawing enclosing boxes using rectage co-ordinates of connected components
	for (int i=1;i<boxVector.size();i++){
		int x1=boxVector[i][0],y1=boxVector[i][1],x2=boxVector[i][2],y2=boxVector[i][3];
		for (int i=x1;i<=x2;i++){
			Image(y1,i,0,0)=0;Image(y1,i,0,1)=0;Image(y1,i,0,2)=255;
			Image(y2,i,0,0)=0;Image(y2,i,0,1)=0;Image(y2,i,0,2)=255;
		}
		for (int i=y1;i<=y2;i++){
			Image(i,x1,0,0)=0;Image(i,x1,0,1)=0;Image(i,x1,0,2)=255;
			Image(i,x2,0,0)=0;Image(i,x2,0,1)=0;Image(i,x2,0,2)=255;
		}
	}
	Image.save("filteredImage12.png");

	// determing rectangle enclosing each of the text blocks
	// finalBoxVector consists of { minX, minY, maxX, maxY }
	// these are co-ordinates of rectangle enclosing single text block
	// each connected component is partitioned if number of black pixels at a certain horizontal 
	// level is less than 20% of total number of pixels in that level
	vector<vector<int> > finalBoxVector;
	int countRed=0;
	for (int i=1;i<boxVector.size();i++){
		int x1=boxVector[i][0],y1=boxVector[i][1],x2=boxVector[i][2],y2=boxVector[i][3];
		for (int x=x1;x<=x2 && x1<x2 && y1<y2;x++){
			vector<int> v(5,0);
			v[4]=i;
			int sum=0;
			for (int y=y1;y<=y2;y++){
				if (label[x][y]==i){
					sum++;
					label1[x][y]=countRed;
				}
			}
			if (5*sum>=(y2-y1)){
				v[0]=x;
				v[1]=y1;
				v[3]=y2;
				x++;
				while (5*sum>=(y2-y1)){
					sum=0;
					for (int y=y1;y<=y2&&x<=x2;y++){
						if (label[x][y]==i){
							sum++;
							label1[x][y]=countRed;
						}
					}
					x++;
				}
				v[2]=x-1;
				if(v[2]-v[0]>5){
					finalBoxVector.push_back(v);
					countRed++;
				}
			}
		}
	}
	
	// drawing enclosing boxes using rectage co-ordinates of text blocks	
	for (int i=0;i<finalBoxVector.size();i++){
		int x1=finalBoxVector[i][0],y1=finalBoxVector[i][1],x2=finalBoxVector[i][2],y2=finalBoxVector[i][3];
		for (int i=x1;i<=x2;i++){
			Image(y1,i,0,0)=255;Image(y1,i,0,1)=0;Image(y1,i,0,2)=0;
			Image(y2,i,0,0)=255;Image(y2,i,0,1)=0;Image(y2,i,0,2)=0;
		}
		for (int i=y1;i<=y2;i++){
			Image(i,x1,0,0)=255;Image(i,x1,0,1)=0;Image(i,x1,0,2)=0;
			Image(i,x2,0,0)=255;Image(i,x2,0,1)=0;Image(i,x2,0,2)=0;
		}
	}
	Image.save("filteredImage13.png");

	// initialising charLabel and charLabel1 
	charLabel.resize(boxVector.size());
	charLabel1.resize(finalBoxVector.size());
	int countChar=0;

	// determing rectangle enclosing single characters
	// finalCharVector consists of { minX, minY, maxX, maxY }
	// these are co-ordinates of rectangle enclosing single character
	// each text block is partitioned if number of black pixels at a certain vertical line is zero
	vector<vector<int> > finalCharVector;
	for (int i=0;i<finalBoxVector.size();i++){
		int x1=finalBoxVector[i][0],y1=finalBoxVector[i][1],x2=finalBoxVector[i][2],y2=finalBoxVector[i][3];
		int lab = finalBoxVector[i][4];
		for(int y=y1;y<=y2;y++){
			int sum=0;
			vector<int> v(6,0);
			v[4]=lab;
			for(int x=x1;x<=x2;x++){
				if(label[x][y]==lab)
					sum+=(sourceImage1(y,x,0,0)==0?1:0);
			}
			if(sum!=0){
				v[0]=x1;
				v[1]=y-1;
				v[2]=x2;
				y++;
				while(sum!=0){
					sum=0;
					for(int x=x1;x<=x2&&y<=y2;x++){
						if(label[x][y]==lab)
							sum+=(sourceImage1(y,x,0,0)==0?1:0);
					}
					y++;
				}
				v[3]=y;
				charLabel[lab].push_back(countChar);
				charLabel1[i].push_back(countChar);
				finalCharVector.push_back(v);
				countChar++;
			}
		}
	}

	// drawing enclosing boxes using rectage co-ordinates of character boxes
	for (int i=0;i<finalCharVector.size();i++){
		int x1=finalCharVector[i][0],y1=finalCharVector[i][1],x2=finalCharVector[i][2],y2=finalCharVector[i][3];
		for (int i=x1;i<=x2;i++){
			Image(y1,i,0,0)=0;Image(y1,i,0,1)=255;Image(y1,i,0,2)=0;
			Image(y2,i,0,0)=0;Image(y2,i,0,1)=255;Image(y2,i,0,2)=0;
			sourceImage1(y1,i,0,0)=0;sourceImage1(y2,i,0,0)=0;
		}
		for (int i=y1;i<=y2;i++){
			Image(i,x1,0,0)=0;Image(i,x1,0,1)=255;Image(i,x1,0,2)=0;
			Image(i,x2,0,0)=0;Image(i,x2,0,1)=255;Image(i,x2,0,2)=0;
			sourceImage1(i,x2,0,0)=0;sourceImage1(i,x1,0,0)=0;
		}
	}
	Image.save("filteredImage14.png");
	sourceImage1.save("filteredImage15.png");
	return finalCharVector;
}

// function to compare two templates by using their normalised 2D correlation and
// correlation between horizontal and vertical projections
float compareTemplates(Template T1,Template T2){
	float correl,sum1=0,sum2=0,sum3=0;
	for (int i=0;i<16;i++){
		for (int j=0;j<16;j++){
			float a=(T1.Image[i][j]-T1.mean);
			float b=(T2.Image[i][j]-T2.mean);
			sum1+=(a*b);
			sum2+=(a*a);
			sum3+=(b*b);
		}
	}
	correl=((sum1*sum1)/(sum2*sum3));
	if(sum1<0)
		correl*=-1;
	float hcorrel,hsum1=0,hsum2=0,hsum3=0;
	float vcorrel,vsum1=0,vsum2=0,vsum3=0;
	for (int i=0;i<16;i++){
		float ha=(T1.horizontalProjection[i]-T1.hmean);
		float hb=(T2.horizontalProjection[i]-T2.hmean);
		float va=(T1.verticalProjection[i]-T1.vmean);
		float vb=(T2.verticalProjection[i]-T2.vmean);
		hsum1+=(ha*hb);
		hsum2+=(ha*ha);
		hsum3+=(hb*hb);
		vsum1+=(va*vb);
		vsum2+=(va*va);
		vsum3+=(vb*vb);
	}
	hcorrel=((hsum1*hsum1)/(hsum2*hsum3));
	if(hsum1<0)
		hcorrel*=-1;
	vcorrel=((vsum1*vsum1)/(vsum2*vsum3));
	if(vsum1<0)
		vcorrel*=-1;
	return (correl);
}

// overloading < operator for Pair struct
bool operator<(Pair p1,Pair p2){
	return (p1.first<p2.first);
}

// checks if given template satisfies basic condition of being character c
bool charMatch(Template T, char c){
	return 1;
}

// takes one character box as input and detects the corresponding character
// using character box a Template object is made and now compares it with all 
// standard templates and adds them to priority queue according to value of 
// correlation between standard template and character template
// checks character properties of top element using charMatch 
// finally returns character found
char detectCharacter(CImg<unsigned int>& sourceImage,charBox B){
	CImg<unsigned int> smallImage(B.y2-B.y1+1,B.x2-B.x1+1,1,1,0);
	for (int i=B.y1;i<=B.y2;i++){
		for (int j=B.x1;j<=B.x2;j++){
			smallImage(i-B.y1,j-B.x1,0,0)=sourceImage(i,j,0,0);
		}
	}
	int x1=-1,x2=-1,y1=-1,y2=-1;
	for (int i=0;i<smallImage.height();i++){
		for (int j=0;j<smallImage.width();j++){
			if (smallImage(j,i,0,0)==0){
				x1=i;
				break;
			}
		}
		if (x1!=-1)
			break;
	}
	for (int i=smallImage.height()-1;i>=0;i--){
		for (int j=0;j<smallImage.width();j++){
			if (smallImage(j,i,0,0)==0){
				x2=i;
				break;
			}
		}
		if (x2!=-1)
			break;
	}
	for (int i=0;i<smallImage.width();i++){
		for (int j=0;j<smallImage.height();j++){
			if (smallImage(i,j,0,0)==0){
				y1=i;
				break;
			}
		}
		if (y1!=-1)
			break;
	}
	for (int i=smallImage.width()-1;i>=0;i--){
		for (int j=0;j<smallImage.height();j++){
			if (smallImage(i,j,0,0)==0){
				y2=i;
				break;
			}
		}
		if (y2!=-1)
			break;
	}
	smallImage.crop(y1,x1,y2,x2);
	smallImage=smallImage.resize(16,16);
	Template T;
	int sum=0;
	for(int k=0;k<16;k++){
		T.horizontalProjection[k]=0;
		T.verticalProjection[k]=0;
	}
	for (int k=0;k<16;k++){
		for (int j=0;j<16;j++){
			if (smallImage(j,k,0,0)==255){
				T.Image[k][j]=1;
				T.horizontalProjection[k]++;
				T.verticalProjection[j]++;
				sum++;
			}
			else{
				T.Image[k][j]=0;
			}
		}
	} 
	int sum1=0,sum2=0;
	for(int k=0;k<16;k++){
		sum1+=T.horizontalProjection[k];
		sum2+=T.verticalProjection[k];
	}
	T.mean=((float)sum/256);
	T.hmean=((float)sum1/16);
	T.vmean=((float)sum2/16);
	float max=0;
	char maxLabel;
	priority_queue<Pair> PQ;
	for (int i=0;i<charTemplates.size();i++){
		float correl=compareTemplates(charTemplates[i],T);
		if (max<correl){
			max=correl;
			maxLabel=charTemplates[i].label;
		}
		if(!(correl!=correl)){
			Pair p = {correl,charTemplates[i].label};
			PQ.push(p);
		}
	}
	cout<<endl;
	for (int k=0;k<16;k++){
		for (int j=0;j<16;j++){
			cout<<T.Image[k][j];
		}
		cout<<endl;
	}
	while(!PQ.empty()&&T.mean!=0){
		Pair p = PQ.top();
		if(charMatch(T,p.second))
			break;
		else
			PQ.pop();
	}
	if(T.mean==0){
		T.label='I';
		return 'I';
	}
	T.label=PQ.top().second;
	cout<<endl;
	return T.label;
}

vector<vector<int> > characterBoxes;
int hi=0,wi=0;

// main function which takes image input and detects text in that image
int main1(char* fileName){
	// initialisations
	label.clear();
	charLabel.clear();
	label1.clear();
	charLabel1.clear();
	characterBoxes.clear();
	vector< vector<int> > structuralElement={{1,1,1},{1,1,1},{1,1,1}}; // structuring element for dilation and erosion
	
	try{
		CImg<unsigned int> Image(fileName); // opening input image
		hi=Image.height();
		wi=Image.width();
		if (Image.width()>800)
			wi=800;
		if (Image.height()>600)
			hi=600;
		Image=Image.resize(wi,hi); // scaling image if its size is greater than 800x600
		Image.save("filteredImage0.png");
		CImg<unsigned int> gImage(Image.width(),Image.height(),1,1,0);	
		int h=Image.height(),w=Image.width();
		CImg<unsigned int> filteredImage1(Image.width(),Image.height(),1,1,0);
		CImg<unsigned int> filteredImage2(Image.width(),Image.height(),1,1,0);
		CImg<unsigned int> filteredImage3(Image.width(),Image.height(),1,1,0);
		CImg<unsigned int> filteredImage4(Image.width(),Image.height(),1,1,0);
		CImg<unsigned int> filteredImage5(Image.width(),Image.height(),1,1,0);
		CImg<unsigned int> filteredImage6(Image.width(),Image.height(),1,1,0);
		CImg<unsigned int> filteredImage7(Image.width(),Image.height(),1,1,0);
		CImg<unsigned int> filteredImage8(Image.width(),Image.height(),1,1,0);
		CImg<unsigned int> filteredImage9(Image.width(),Image.height(),1,1,0);
		CImg<unsigned int> filteredImage10(Image.width(),Image.height(),1,1,0);
		CImg<unsigned int> filteredImage11(Image.width(),Image.height(),1,1,0);

		// matrix for gaussian blur 
		vector< vector<int> > Vector1={{2,4,5,4,2},{4,9,12,9,4},{5,12,15,12,5},{4,9,12,9,4},{2,4,5,4,2}};
		// matrices for sobel operators
		vector< vector<int> > Vector2={{-1,0,1},{-2,0,2},{-1,0,1}};
		vector< vector<int> > Vector3={{1,2,1},{0,0,0},{-1,-2,-1}};
		// Convolution filters
		Convolution Filter1(5,5,159,Vector1); // Gaussian 
		Convolution Filter2(3,3,1,Vector2); // Vertical Sobel
		Convolution Filter3(3,3,1,Vector3); // Horizontal Sobel

		// converting to grayscale
		if (Image.spectrum()==3)
			grayscale(Image,gImage);
		else
			gImage=Image;
		gImage.save("grayscaled.png");

		// applying gaussian filter
		Filter1.operate(gImage,filteredImage1);
		filteredImage1.save("filteredImage1.png");

		// applying sobel operators (horizontal and vertical) and taking modular sum 
		Filter2.operate(filteredImage1,filteredImage2);
		filteredImage2.save("filteredImage2.png");
		Filter3.operate(filteredImage1,filteredImage3);
		filteredImage3.save("filteredImage3.png");
		modImage(filteredImage2,filteredImage3,filteredImage4);
		filteredImage4.save("filteredImage4.png");

		// opening Image
		Erosion(filteredImage4,filteredImage5);
		filteredImage5.save("filteredImage5.png");
		Dilation(filteredImage4,filteredImage6);
		filteredImage6.save("filteredImage6.png");

		// closing Image
		Dilation(filteredImage5,filteredImage7);
		filteredImage7.save("filteredImage7.png");
		Erosion(filteredImage6,filteredImage8);
		filteredImage8.save("filteredImage8.png");

		// averaging opened and closed Images
		avgImage(filteredImage6,filteredImage8,filteredImage9);
		filteredImage9.save("filteredImage9.png");

		// binarising averaged image and grayscaled image
		binariseImage(filteredImage9,filteredImage10);
		filteredImage10.save("filteredImage10.png");
		binariseImage(gImage,filteredImage11);
		filteredImage11.save("filteredImage11.png");

		// calling extract character functions to generate character box vector
		characterBoxes=extractCharacters(filteredImage10,filteredImage11,filteredImage2, Image);
		// for each character box calling detect character function
		for (int i=0;i<characterBoxes.size();i++){
			charBox B={characterBoxes[i][0],characterBoxes[i][1],characterBoxes[i][2],characterBoxes[i][3]};
			characterBoxes[i][5]=detectCharacter(filteredImage11,B);
			cout<<characterBoxes[i][5]<<endl;
		}
	}
	catch(CImgException &e){
		cout<<"Open Image file only\n";
	}
	return 0;
}

// Image object    
GtkWidget *image;

// function to close window
void destroy(GtkWidget *widget, gchar* data){
    cout<<data<<endl;
    gtk_main_quit();
}

// callback function of open button to open file chooser dialog box and then calling main1 function
void openCall(GtkWidget *button, gchar* data){
    GtkWidget *dialog;
    GtkWidget *toplevel = gtk_widget_get_toplevel (button);
    if (gtk_widget_is_toplevel (toplevel)){
       /* Perform action on toplevel. */
        char* filename;
        dialog = gtk_file_chooser_dialog_new ("Choose an Image File to Read",
                              GTK_WINDOW(toplevel),
                              GTK_FILE_CHOOSER_ACTION_OPEN,
                              GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                              GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT,
                              NULL);
        gtk_file_chooser_set_current_folder (GTK_FILE_CHOOSER (dialog),"./");
        if (gtk_dialog_run (GTK_DIALOG (dialog)) == GTK_RESPONSE_ACCEPT){
            filename = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (dialog));
            int i;
            for(i=0;filename[i]!=0;i++){
                data[i]=filename[i];
            }
            data[i]=0;
            cout<<"Opened File"<<data<<endl;
            main1(data);
            gtk_image_set_from_file(GTK_IMAGE(image), "filteredImage0.png");
          }
        gtk_widget_destroy (dialog);
    }
}

// function for writing extracted text to pdf
void toTextCall(GtkWidget *button, gchar* data){
	// opening a latex file first we write the extracted text to that file
	// finally we use system commands to convert that latex file to pdf
	FILE *fp;
	char fname[1000];
	for(int i=0;data[i]!=0;i++){
		fname[i]=data[i];
		if(data[i]=='.'){
			fname[i+1]='t';
			fname[i+2]='e';
			fname[i+3]='x';
			fname[i+4]=0;
			break;
		}
	}
	string command = "pdflatex ";
	command+=string(fname);
	fp = fopen(fname,"w");
	fprintf(fp, "\\documentclass{article}\n\\begin{document}\n");

	// picking each text block one by one and detecting the text contained in it
	// then writing the extracted text to latex file
	for(int j=0;j<charLabel1.size();j++){
		list<int> charList=charLabel1[j];
        list<int>::iterator it;
        char* str = new char[10000];
        int i=0;
        for (it=charList.begin();it!=charList.end();it++){
        	str[i++]=characterBoxes[*it][5];
        }
        str[i]=0;
        fprintf(fp, "%s ", str);
    }
	fprintf(fp, "\\end{document}");
	fclose(fp);
	system(command.c_str());
	system("rm *.tex *.aux *.log");
}

// callback functions of different buttons to change image shown on interface
void sourceCall(GtkWidget *button, gpointer data){
    gtk_image_set_from_file(GTK_IMAGE(image),"filteredImage0.png");
}
void gaussCall(GtkWidget *button, gpointer data){
    gtk_image_set_from_file(GTK_IMAGE(image),"filteredImage1.png");
}
void verSobelCall(GtkWidget *button, gpointer data){
    gtk_image_set_from_file(GTK_IMAGE(image),"filteredImage3.png");
}
void horSobelCall(GtkWidget *button, gpointer data){
    gtk_image_set_from_file(GTK_IMAGE(image),"filteredImage2.png");
}
void edgeImageCall(GtkWidget *button, gpointer data){
    gtk_image_set_from_file(GTK_IMAGE(image),"filteredImage4.png");
}
void openImageCall(GtkWidget *button, gpointer data){
    gtk_image_set_from_file(GTK_IMAGE(image),"filteredImage6.png");
}
void closeImageCall(GtkWidget *button, gpointer data){
    gtk_image_set_from_file(GTK_IMAGE(image),"filteredImage8.png");
}
void avgImageCall(GtkWidget *button, gpointer data){
    gtk_image_set_from_file(GTK_IMAGE(image),"filteredImage9.png");
}
void binImageCall(GtkWidget *button, gpointer data){
    gtk_image_set_from_file(GTK_IMAGE(image),"filteredImage10.png");
}
void connCompCall(GtkWidget *button, gpointer data){
    gtk_image_set_from_file(GTK_IMAGE(image),"filteredImage13.png");
}
void charImageCall(GtkWidget *button, gpointer data){
    gtk_image_set_from_file(GTK_IMAGE(image),"filteredImage14.png");
}

// gets co-ordinates of mouse click on image and finds the text box 
// corresponding to that pixel and prints characters corresponding to that label
static gint clickCall(GtkWidget *widget, GdkEventButton *event, gpointer data){
    cout<<"here";
    if (event->button == 1){
        int x=event->x-405+(wi/2);
        int y=event->y-305+(hi/2);
        cout<<x<<" "<<y<<endl;
        if(x>=0 && y>=0 &&wi>0 && hi>0 && x<=wi && y<=hi){
	        list<int> charList=charLabel1[label1[y][x]];
	        list<int>::iterator it;
	        char* str = new char[10000];
	        int i=0;
	        for (it=charList.begin();it!=charList.end();it++){
	        	cout<<*it<<"  ";
	        	str[i++]=characterBoxes[*it][5];
	        }
	        str[i]=0;
	        cout<<endl<<str[i]<<endl;
	        gtk_label_set_text(GTK_LABEL(data),str);
	    }
    }
}

// main function
int main (int argc, char *argv[]){
	// accessing database to load standard character templates
	FILE* fp;
	fp=fopen("char_dict.bin","rb");
	Template T;
	while (!feof(fp)){
		fread(&T,sizeof(T),1,fp);
		charTemplates.push_back(T);
	}
	fclose(fp);

    GtkWidget *window;
    GtkWidget *hbox1;
    GtkWidget *vbox1, *vbox2;
    GtkWidget *label1;
    GtkWidget *frame1,*frame2,*frame3;
    GtkWidget *open, *toText, *Exit, *source, *gauss, *verSobel, *horSobel, *edgeImage, *openImage, *closeImage, *avgImage, *binImage;
    GtkWidget *connComp, *charImage;
    GtkWidget *clickBox;
    gchar* filename = new gchar[100];

    // Initialize the GTK+ and all of its supporting libraries
    gtk_init (&argc, &argv);

    // Create a new window, give it a title and display it to the user
    window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
    gtk_window_set_position(GTK_WINDOW(window), GTK_WIN_POS_CENTER);
    gtk_window_set_title (GTK_WINDOW (window), "Text Detection");
    gtk_window_set_default_size (GTK_WINDOW(window),910,700);

    // Connect the main window to the destroy 
    g_signal_connect(G_OBJECT(window), "destroy", G_CALLBACK(destroy), NULL);

    // initialising click box, image and putting image inside click box
    clickBox=gtk_event_box_new();
    image = gtk_image_new_from_file("a.png");
    gtk_container_add(GTK_CONTAINER(clickBox), image);
    gtk_widget_set_events(clickBox, GDK_BUTTON_PRESS_MASK);

    // label 1 shows detected text on the interface
    label1 = gtk_label_new("Detected Text");
    // initialising different buttons
    open = gtk_button_new_with_label("Open");
    toText = gtk_button_new_with_label("Write to PDF");
    Exit = gtk_button_new_with_label("Exit");
    source = gtk_button_new_with_label("Source Image");
    gauss = gtk_button_new_with_label("Gaussian Blur");
    verSobel = gtk_button_new_with_label("Sobel Vertical");
    horSobel = gtk_button_new_with_label("Sobel Horizontal");
    edgeImage = gtk_button_new_with_label("Edge Detected");
    openImage = gtk_button_new_with_label("Opened Image");
    closeImage = gtk_button_new_with_label("Closed Image");
    avgImage = gtk_button_new_with_label("Average Image");
    binImage = gtk_button_new_with_label("Binary Image");
    connComp = gtk_button_new_with_label("Connected Parts");
    charImage = gtk_button_new_with_label("Character Partition");

  	// connecting call back functions to buttons
    g_signal_connect(G_OBJECT(open),"clicked", G_CALLBACK(openCall), filename);
    g_signal_connect(G_OBJECT(toText),"clicked", G_CALLBACK(toTextCall), filename);
    g_signal_connect(G_OBJECT(source),"clicked", G_CALLBACK(sourceCall), NULL);
    g_signal_connect(G_OBJECT(gauss),"clicked", G_CALLBACK(gaussCall), NULL);
    g_signal_connect(G_OBJECT(verSobel),"clicked", G_CALLBACK(verSobelCall), NULL);
    g_signal_connect(G_OBJECT(horSobel),"clicked", G_CALLBACK(horSobelCall), NULL);
    g_signal_connect(G_OBJECT(edgeImage),"clicked", G_CALLBACK(edgeImageCall), NULL);
    g_signal_connect(G_OBJECT(openImage),"clicked", G_CALLBACK(openImageCall), NULL);
    g_signal_connect(G_OBJECT(closeImage),"clicked", G_CALLBACK(closeImageCall), NULL);
    g_signal_connect(G_OBJECT(avgImage),"clicked", G_CALLBACK(avgImageCall), NULL);
    g_signal_connect(G_OBJECT(binImage),"clicked", G_CALLBACK(binImageCall), NULL);
    g_signal_connect(G_OBJECT(connComp),"clicked", G_CALLBACK(connCompCall), NULL);
    g_signal_connect(G_OBJECT(charImage),"clicked", G_CALLBACK(charImageCall), NULL);
    g_signal_connect(G_OBJECT(Exit), "clicked", G_CALLBACK(destroy), filename);

    // making boxes for layout of the interface
    // interface consists of a horizontal box (hbox1) inside which there is a 
    // vertical box (vbox1) and a frame (frame3), vbox1 consists of frame1 and frmae2
    // frame1 contains clickbox containing image and frame2 contains label of 
    // detected text, frame3 contains toolbar which is vbox2 which consists of all the buttons 
    hbox1 = gtk_hbox_new(FALSE, 1);
    vbox1 = gtk_vbox_new(FALSE, 1);
    vbox2 = gtk_vbox_new(FALSE, 1);
    gtk_box_pack_start(GTK_BOX(vbox2),open,TRUE,FALSE,2);
    gtk_box_pack_start(GTK_BOX(vbox2),toText,TRUE,FALSE,2);
    gtk_box_pack_start(GTK_BOX(vbox2),source,TRUE,FALSE,2);
    gtk_box_pack_start(GTK_BOX(vbox2),gauss,TRUE,FALSE,2);
    gtk_box_pack_start(GTK_BOX(vbox2),verSobel,TRUE,FALSE,2);
    gtk_box_pack_start(GTK_BOX(vbox2),horSobel,TRUE,FALSE,2);
    gtk_box_pack_start(GTK_BOX(vbox2),edgeImage,TRUE,FALSE,2);
    gtk_box_pack_start(GTK_BOX(vbox2),openImage,TRUE,FALSE,2);
    gtk_box_pack_start(GTK_BOX(vbox2),closeImage,TRUE,FALSE,2);
    gtk_box_pack_start(GTK_BOX(vbox2),avgImage,TRUE,FALSE,2);
    gtk_box_pack_start(GTK_BOX(vbox2),binImage,TRUE,FALSE,2);
    gtk_box_pack_start(GTK_BOX(vbox2),connComp,TRUE,FALSE,2);
    gtk_box_pack_start(GTK_BOX(vbox2),charImage,TRUE,FALSE,2);
    gtk_box_pack_start(GTK_BOX(vbox2),Exit,TRUE,FALSE,2);

    frame1 = gtk_frame_new("Image");
    gtk_frame_set_shadow_type(GTK_FRAME(frame1), GTK_SHADOW_ETCHED_IN);
    gtk_signal_connect (GTK_OBJECT (clickBox),"button_press_event",G_CALLBACK(clickCall), label1);

    frame2 = gtk_frame_new("Text");
    gtk_frame_set_shadow_type(GTK_FRAME(frame2), GTK_SHADOW_ETCHED_IN);
    frame3 = gtk_frame_new("Toolbar");
    gtk_frame_set_shadow_type(GTK_FRAME(frame3), GTK_SHADOW_ETCHED_IN);

    gtk_container_add(GTK_CONTAINER(frame1), clickBox);
    gtk_widget_set_size_request(frame1, 810, 620);
    gtk_container_add(GTK_CONTAINER(frame2), label1);
    gtk_widget_set_size_request(frame2, 810, 50);

    gtk_container_add(GTK_CONTAINER(frame3), vbox2);
    gtk_widget_set_size_request(frame3, 150, 700);

    gtk_box_pack_start(GTK_BOX(vbox1),frame1,TRUE,FALSE,5);
    gtk_box_pack_start(GTK_BOX(vbox1),frame2,TRUE,FALSE,5);

    gtk_box_pack_start(GTK_BOX(hbox1),vbox1,FALSE,FALSE,5);
    gtk_box_pack_start(GTK_BOX(hbox1),frame3,FALSE,FALSE,5);
    
    gtk_window_set_resizable(GTK_WINDOW(window), FALSE);
    gtk_container_add(GTK_CONTAINER(window), hbox1);

    gtk_widget_show_all (window);
    // Hand control over to the main loop
    gtk_main ();
    return 0;
}