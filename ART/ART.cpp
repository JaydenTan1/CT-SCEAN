//Author: Zhijie Tan
//Time: 4/3/2025
#include <iostream>
#include "Display.h"
#include <cstring>
#include <math.h>
#include <cmath>
#include <vector>
using namespace std;

const int n = 400;
const int N = 400; // Image Size N x N
int Pixels[N][N];  //The original image

const double PI = 3.14159265358979323846;

const int ROWS = N;
const int COLS = N;
const int IMAGE_SIZE = n;     //The image size is 400×400
const int ANGLE_COUNT = 180;    //Number of angles (0°~179°)
const int DETECTOR_COUNT = 400; // Number of detectors

double pi = 3.14;
int dst[N];

struct Point {
    int x;
    int y;
};

// ART Reconstruction
void ARTReconstruction(
    int measuredProjections[ANGLE_COUNT][DETECTOR_COUNT],
    double image[IMAGE_SIZE][IMAGE_SIZE],
    int iterations,
    double lambda0,
    double alpha
) {
    // Initialize the image as a constant value to prevent the update amount from being too small due to all zeros
    for(int i=0; i<IMAGE_SIZE; i++){
        for(int j=0; j<IMAGE_SIZE; j++){
            image[i][j] = 50.0;
        }
    }

    //The outer iteration loop
    for(int iter = 0; iter < iterations; iter++){
        
        // Calculate the step size lambda_k for the current iteration round
        // Formula:  lambda_k = lambda0 / (1 + alpha * iter)
        double lambda_k = lambda0 / (1.0 + alpha * iter);

        // Iterate over all angles
        for(int angleIndex = 0; angleIndex < ANGLE_COUNT; angleIndex++){
            double angleRad = -angleIndex*(PI/180.0);
            double cosT = cos(angleRad);
            double sinT = sin(angleRad);

            // Iterate over all detectors
            for(int det = 0; det < DETECTOR_COUNT; det++){
                // The offset of the ray in the orthogonal direction
                double s = det - DETECTOR_COUNT/2.0;

                // dot(a_i, x) & ||a_i||^2
                double ai_dot_x = 0.0;
                double ai_norm_sq = 0.0;

                vector<Point> rayPixels;

                //Ray integration: Sample along the (cosT, sinT) direction
                double maxT = IMAGE_SIZE * 1.5;
                double step = 0.2;
                for(double tVal = -maxT; tVal <= maxT; tVal += step){
                    double X = (IMAGE_SIZE/2.0) + s*(-sinT) + tVal*cosT;
                    double Y = (IMAGE_SIZE/2.0) + s*( cosT) + tVal*sinT;
                    int iX = (int)round(X);
                    int iY = (int)round(Y);

                    if(iX >= 0 && iX < IMAGE_SIZE && iY >= 0 && iY < IMAGE_SIZE){
                        // Times the step approximation integral
                        ai_dot_x += image[iX][iY] * step;
                        ai_norm_sq += step;
                        rayPixels.push_back({iX, iY});
                    }
                }

                
                double bi = measuredProjections[angleIndex][det];

                //correction = lambda_k * (b_i - dot)/norm
                if(ai_norm_sq > 0){
                    double diff = bi - ai_dot_x;
                    double correction = lambda_k * (diff / ai_norm_sq);

                    // The same correction is applied to all pixels on the shooting line
                    for(auto &pt : rayPixels){
                        image[pt.x][pt.y] += correction;
                    }
                }
            } // end for det
        } // end for angle
    } // end for iter
}


int main() {
    int sinogram[n][n] = {0};
    ifstream inFile;
     inFile.open("sinogram_400x400.txt");	
	for (int row = 0; row < n; row++)
		for (int col = 0; col < n; col++)
			inFile >> sinogram[row][col];
	
	inFile.close();

    int Pixels[n][n] = {0};
    
    int ReconstructedImage[N][N] = {0};

    static double reconImage[N][N] = {0};

     // ART
     int iterations = 200;
    double lambda0 = 0.2; // Initial step size
    double alpha   = 0.003; // Attenuation coefficient, adjustable
    ARTReconstruction(sinogram, reconImage, iterations, lambda0, alpha);

   //Do visual normalization
    double bgSum = 0.0;
    int bgCount = 0;
    for(int r = 0; r < 400; r++){
        for(int c = 0; c < 400; c++){
            double dx = r - 200.0; // Suppose image center (200,200)
            double dy = c - 200.0;
            double dist = sqrt(dx*dx + dy*dy);

            // Assume a radius of 180~190 as the background
            if(dist > 180 && dist < 190){
                bgSum += reconImage[r][c];
                bgCount++;
            }
        }
    }
    double bgMean = (bgCount>0) ? bgSum / bgCount : 0.0;

    // Subtract bgMean from each pixel of the image, and truncate to 0 if it is less than 0
    for(int r=0; r<400; r++){
        for(int c=0; c<400; c++){
            reconImage[r][c] -= 0.8 * bgMean;
            if(reconImage[r][c] < 0)
                reconImage[r][c] = 0;
        }
    }
    double minVal = 1e9, maxVal = -1e9;
    for(int r = 0; r < N; r++){
        for(int c = 0; c < N; c++){
            if(reconImage[r][c] < minVal) minVal = reconImage[r][c];
            if(reconImage[r][c] > maxVal) maxVal = reconImage[r][c];
        }
    }
    double range = maxVal - minVal;
    if(range < 1e-12) range = 1.0;

    //Linear mapping
    for(int r = 0; r < N; r++){
        for(int c = 0; c < N; c++){
            double val = (reconImage[r][c] - minVal)/range * 255.0;
            reconImage[r][c] = val;
        }
    }
   
    static int finalImage[N][N];
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            finalImage[i][j] = (int) round(reconImage[i][j]);
        }
    }

    Display image1(finalImage, "ReconstructedImage_Final");
    return 0;
}
