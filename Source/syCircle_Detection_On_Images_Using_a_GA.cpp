//
//   Module: "syCircle_Detection_On_Images_Using_a_GA.cpp"
//   Looking for a Circle in an Image (Using a Genetic Algorithm)
//
//   Implementation of the work described in
//
//   Ayala-Ramirez, V., Garcia-Capulin, C. H., Perez-Garcia, A., Sanchez-Yanez, R. E. (2006).
//   Circle detection on images using genetic algorithms. Pattern Recognition Letters, 27(6):652-657.
//   DOI: 10.1016/j.patrec.2005.10.003. (Elsevier)
//   https://www.sciencedirect.com/science/article/abs/pii/S016786550500293X
//
//   Required code modules:
//      syGenetico.cpp          syGenetico.hpp
//      syGeneticoCirculo.cpp   syGeneticoCirculo.hpp
//
//   Raúl E. Sánchez-Yáñez
//   Rev 2022-02-23


#include <iostream>
using namespace std;
#include "opencv2/core.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/opencv.hpp"
using namespace cv;


#include "../Headers/AGScirculo.hpp"


// Function prototypes
int syImCircleDetection_TEST(void);
int syImCircleDetectionUsingGA(Mat src, Mat &Edges, Mat &dst, float &X0, float &Y0, float &R0);

int syImEdges(Mat &src);
int syImNegative(Mat &src);
int syImCountBlackPixels(Mat src, int &nBlackPix);
int syImGetBlackPixelLocation(Mat src,float x[],float y[],int nEdges);

int main(void)
{
   syImCircleDetection_TEST();
   return 0;
}

int syImCircleDetection_TEST(void)
{
   int chosenImage;
   cout << "Looking for a Circle in the Image.\nFor exit, press any key\n";
   do
   {
      cout << "\n\n\nChoose the input image:\n";
      cout << "1. Spiral.\n";
      cout << "2. Moon.\n";
      cout << "3. Cactus.\n";
      cout << "0. Exit.\n";
      cout << "\tYour choice: ? ";
      cin >> chosenImage; cin.ignore();
   } while (chosenImage<0 || chosenImage>3);

   Mat Input;
   switch(chosenImage)
   {
      case 1: Input = imread("circ_espiral.jfif",IMREAD_UNCHANGED); break;
      case 2: Input = imread("circ_luna.jfif",IMREAD_UNCHANGED); break;
      case 3: Input = imread("circ_cactus.jfif",IMREAD_UNCHANGED); break;
      case 0:
      default: cout << "\n\n\nExit...\n";
   }

   if(Input.empty())
      { cout << "\n\n No input file or directory.\n"; return -1; }
   cout << "\n\nThe "<<Input.cols<<"x"<<Input.rows<<" image has "<<Input.rows*Input.cols<<" pixels.\n";

   Mat Edges, Output; // Edge image and an Output with a circle overlaid
   float X0=0.0, Y0=0.0, Radius=0.0; // Coordinates of center (X0,Y0) and radius of circle

   // Call to the main function of this module:
   syImCircleDetectionUsingGA(Input,Edges,Output,X0,Y0,Radius);

   // Draw the outcome
   imshow("Input Image",Input);
   imshow("Edges",Edges);
   imshow("Outcome",Output);

   circle(Input, Point(Y0,X0),Radius, Scalar(0,0,255), 2/*thickness*/, LINE_8, 0/*shift or not*/);
   imshow("Outcome on Input",Input);

   waitKey();
   return 0;
}

int syImCircleDetectionUsingGA(Mat Input, Mat &Edges, Mat &Output, float &X0, float &Y0, float &R0)
{
   // IMAGE PREPROCESSING

   Edges = Input.clone(); // A clone of input, the image to be processed
   syImEdges(Edges); // A 3-channel image with white edges (black background)
   syImNegative(Edges); // Black edges (white background)

   int nEdges;    // For the number of black edges in the image
   syImCountBlackPixels(Edges,nEdges);
   cout << "The image has "<<nEdges<<" edges ("<<100*nEdges/(Edges.rows*Edges.cols)<<"%).\n";

   float r[nEdges],c[nEdges]; // For edge location: k-th edge at (r[k],c[k]) array values
   if ( syImGetBlackPixelLocation(Edges,r,c,nEdges) != nEdges )
      cout<< "\n\nError while retrieving edge location.\n";



   // CALL TO THE GENETIC ALGORITHM
   cout << "\n\n\nCalling to the GA.\n\n";

   Genetico_Circulo_Parametros(r,c,nEdges,X0,Y0,R0);
//   Genetico_Circulo_Parametros();

   cout << "\n\nBest fit is circle with center ("<<(int)X0<<","<<(int)Y0<<") and radius "<< (int)R0 << endl;


   // Draw the outcome
   Output = Edges.clone();
   circle(Output, Point((int)Y0,(int)X0),(int)R0, Scalar(0,0,255), 2/*thickness*/, LINE_8, 0/*shift or not*/);
   return 0;
}


int syImGetBlackPixelLocation(Mat src,float x[],float y[],int nEdges)
{
   int counter = 0; // To verify the number of edge pixels

/*   if(src.channels()==1 && src.type()==CV_8UC1)
      for(int r=0; r < src.rows; r++)
      {
         uchar *pAux = src.ptr<uchar>(r);
         for(int c=0; c < src.cols; c++)
            if(pAux[c]==0)
            {
               if(counter+1>nEdges)
                  { cout << "\n\nERROR processing mono image"; return -1; }
               x[counter] = r;
               y[counter] = c;
               counter++;
            }
      }
   else */
   if(src.channels()==3 && src.type()==CV_8UC3)
      for(int r=0; r < src.rows; r++)
      {
         Vec3b* pAtSrc = src.ptr<Vec3b>(r);
         for(int c = 0; c < src.cols; c++)
            if(!pAtSrc[c][0] && !pAtSrc[c][1] && !pAtSrc[c][2]) // black pixel
            {
               if(counter+1>nEdges)
                  { cout << "\n\nERROR processing 3-channel image"; return -1; }
               x[counter] = r;
               y[counter] = c;
               counter++;
            }
      }
   return counter;
}

int syImCountBlackPixels(Mat src, int &nBlackPix)
{
   if(src.empty())
      { cout << "\n\nSource is empty.\n"; return -1; }

   nBlackPix = 0;
/*   if(src.channels()==1 && src.type()==CV_8UC1)
      for(int r=0; r < src.rows; r++)
      {
         uchar *pAux = src.ptr<uchar>(r);
         for(int c=0; c < src.cols; c++)
            if(pAux[c]==0)
               nBlackPix++;
      }
   else */
   if(src.channels()==3 && src.type()==CV_8UC3)
      for(int r=0; r < src.rows; r++)
      {
         Vec3b* pAtSrc = src.ptr<Vec3b>(r);
         for(int c = 0; c < src.cols; c++)
            if(!pAtSrc[c][0] && !pAtSrc[c][1] && !pAtSrc[c][2]) // black pixel
               nBlackPix++;
      }
   return 0;
}

int syImEdges(Mat &src)
{
   if(src.empty())
      { cout << "\n\nEmpty source.\n"; return -1; }

   if(src.channels()==1)
   {
      GaussianBlur(src,src,Size(7,7),1.5,1.5);
      Canny(src,src,0,100,3);
   }
   else if(src.channels()==3)
   {
      cvtColor(src,src,COLOR_BGR2GRAY); // mono image
      GaussianBlur(src,src,Size(7,7),1.5,1.5);
      Canny(src,src,0,100,3);
      cvtColor(src,src,COLOR_GRAY2BGR);
   }
   return 0;
}

int syImNegative(Mat &src)
{
   if(src.empty())
      { cout << "\n\nEmpty source.\n"; return -1; }

   if(src.channels()==1)
      src = 255-src;
   else if(src.channels()==3)
      src = Scalar(255,255,255)-src;
   return 0;
}
