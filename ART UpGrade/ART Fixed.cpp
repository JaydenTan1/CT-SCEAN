#include <iostream>
#include "Display.h"
#include <cstring>
#include <math.h>
#include <cmath>
#include <vector>
#include <limits>
#include <fstream>
using namespace std;

const int n = 400;
const int N = 400; // Image Size N x N
const double PI = 3.14159265358979323846;

const int IMAGE_SIZE = 400;     // Image Size N x N
const int ANGLE_COUNT = 180;    // Number of angles (0°~179°)
const int DETECTOR_COUNT = 400; // Number of detectors

struct Point {
    int x;
    int y;
    double weight;
};

double Table2[n][n] = {0};


void ARTReconstruction(
    int measuredProjections[ANGLE_COUNT][DETECTOR_COUNT],
    double image[IMAGE_SIZE][IMAGE_SIZE],
    int iterations,
    double lambda0,
    double alpha
) {
   // 1) Initialize the image to zero
    for(int i=0; i<IMAGE_SIZE; i++){
        for(int j=0; j<IMAGE_SIZE; j++){
            image[i][j] = 0.0;
        }
    }

    // 2) outer iteration loop
    for(int iter = 0; iter < iterations; iter++){
        
        // Calculate the step size lambda_k for the current iteration round
        double lambda_k = lambda0 / (1.0 + alpha * iter);

        // Iterate over all angles
        for(int angleIndex = 0; angleIndex < ANGLE_COUNT; angleIndex++){
            // 修正角度：确保与sinogram的几何一致
            double angleRad = -angleIndex * (PI / 180.0);
            double cosT = cos(angleRad);
            double sinT = sin(angleRad);

            // 遍历所有探测器
            for(int det = 0; det < DETECTOR_COUNT; det++){
                // Fix Angle: make sure it is geometrically consistent with the sinogram
                double s = det - DETECTOR_COUNT/2.0 + 0.5;

                // Calculate ray integrals and weights
                double ai_dot_x = 0.0;
                double ai_norm_sq = 0.0;
                vector<Point> rayPixels;

                // Use finer ray sampling
                double step = 0.3;  //Reduce step size to improve accuracy
                double maxT = IMAGE_SIZE * sqrt(2.0);
                
                for(double tVal = -maxT; tVal <= maxT; tVal += step){
                    double X = (IMAGE_SIZE/2.0 - 0.5) + s*(-sinT) + tVal*cosT;
                    double Y = (IMAGE_SIZE/2.0 - 0.5) + s*( cosT) + tVal*sinT;
                    
                    // Use bilinear interpolation
                    int iX = (int)floor(X);
                    int iY = (int)floor(Y);
                    
                    if(iX >= 0 && iX < IMAGE_SIZE-1 && iY >= 0 && iY < IMAGE_SIZE-1){
                        double fx = X - iX;
                        double fy = Y - iY;
                        
                        // bilinear interpolation weights
                        double w00 = (1-fx) * (1-fy);
                        double w01 = (1-fx) * fy;
                        double w10 = fx * (1-fy);
                        double w11 = fx * fy;
                        
                        // Accumulate ray integrals
                        double interpolated = w00 * image[iX][iY] + 
                                            w01 * image[iX][iY+1] +
                                            w10 * image[iX+1][iY] + 
                                            w11 * image[iX+1][iY+1];
                        
                        ai_dot_x += interpolated * step;
                        ai_norm_sq += step * step;
                        
                        // Record the pixels and weights
                        if(w00 > 1e-6) rayPixels.push_back({iX, iY, w00 * step});
                        if(w01 > 1e-6) rayPixels.push_back({iX, iY+1, w01 * step});
                        if(w10 > 1e-6) rayPixels.push_back({iX+1, iY, w10 * step});
                        if(w11 > 1e-6) rayPixels.push_back({iX+1, iY+1, w11 * step});
                    }
                }

                double bi = measuredProjections[angleIndex][det];

                // ART update
                if(ai_norm_sq > 1e-12){
                    double diff = bi - ai_dot_x;
                    double correction = lambda_k * diff / ai_norm_sq;

                    // Apply a weighted correction to all pixels on the projection line
                    for(auto &pt : rayPixels){
                        image[pt.x][pt.y] += correction * pt.weight;
                        // Non-negative constraints
                        if(image[pt.x][pt.y] < 0) image[pt.x][pt.y] = 0;
                    }
                }
            }
        }
        
       // Calculate the error
        double error = 0.0;
        if(iter > 0){
            for(int i=0; i<N; i++){
                for(int j=0; j<N; j++){
                    error += pow((image[i][j] - Table2[i][j]),2);
                }
            }
        }
        cout << iter <<" times iterations error:" << error << endl;
        
        // Save current image for next error calculation
        for(int i=0; i<N; i++){
            for(int j=0; j<N; j++){
                Table2[i][j] = image[i][j];
            }
        }
    }
}

void fix(double reconImage[N][N], int finalImage[N][N]) {
    // Apply the circular mask
    double centerX = IMAGE_SIZE / 2.0;
    double centerY = IMAGE_SIZE / 2.0;
    double maxRadius = IMAGE_SIZE / 2.0 - 20;  // Shrink radius to avoid border artifacts
    
    for(int r = 0; r < N; r++){
        for(int c = 0; c < N; c++){
            double dx = r - centerX;
            double dy = c - centerY;
            double dist = sqrt(dx*dx + dy*dy);
            
            if(dist > maxRadius){
                reconImage[r][c] = 0;
            }
        }
    }
    
    // Background correction - only in valid area
    double bgSum = 0.0;
    int bgCount = 0;
    for(int r = 0; r < 400; r++){
        for(int c = 0; c < 400; c++){
            double dx = r - 200.0;
            double dy = c - 200.0;
            double dist = sqrt(dx*dx + dy*dy);

            // Calculate the background in the smaller ring area
            if(dist > 160 && dist < 180){
                bgSum += reconImage[r][c];
                bgCount++;
            }
        }
    }
    double bgMean = (bgCount>0) ? bgSum / bgCount : 0.0;

    // Subtract the background
    for(int r=0; r<400; r++){
        for(int c=0; c<400; c++){
            reconImage[r][c] -= 0.5 * bgMean;  // Adjust background minus scale
            if(reconImage[r][c] < 0)
                reconImage[r][c] = 0;
        }
    }
    
    // Find the minimum and maximum valid region
    double minVal = 1e9, maxVal = -1e9;
    for(int r = 0; r < N; r++){
        for(int c = 0; c < N; c++){
            double dx = r - centerX;
            double dy = c - centerY;
            double dist = sqrt(dx*dx + dy*dy);
            
            if(dist <= maxRadius && reconImage[r][c] > 0){
                if(reconImage[r][c] < minVal) minVal = reconImage[r][c];
                if(reconImage[r][c] > maxVal) maxVal = reconImage[r][c];
            }
        }
    }
    
    double range = maxVal - minVal;
    if(range < 1e-12) range = 1.0;

    // Map linearly to [0, 255]
    for(int r = 0; r < N; r++){
        for(int c = 0; c < N; c++){
            if(reconImage[r][c] > 0){
                double val = (reconImage[r][c] - minVal)/range * 255.0;
                reconImage[r][c] = val;
            } else {
                reconImage[r][c] = 0;
            }
        }
    }
   
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            finalImage[i][j] = (int) round(reconImage[i][j]);
        }
    }
}

double computeSSIM(const int img1[][n],const int img2[][n],int width, int height,double L = 255.0){
    int N_pixels = width * height;
    
   // Calculate the mean
    double mu1 = 0.0, mu2 = 0.0;
    for(int i = 0; i < height; ++i) {
        for(int j = 0; j < width; ++j) {
            mu1 += static_cast<double>(img1[i][j]);
            mu2 += static_cast<double>(img2[i][j]);
        }
    }
    mu1 /= N_pixels;
    mu2 /= N_pixels;

   // Calculate the variance and covariance
    double var1 = 0.0, var2 = 0.0, cov = 0.0;
    for(int i = 0; i < height; ++i) {
        for(int j = 0; j < width; ++j) {
            double d1 = static_cast<double>(img1[i][j]) - mu1;
            double d2 = static_cast<double>(img2[i][j]) - mu2;
            var1 += d1 * d1;
            var2 += d2 * d2;
            cov  += d1 * d2;
        }
    }
    var1 /= (N_pixels - 1);
    var2 /= (N_pixels - 1);
    cov  /= (N_pixels - 1);

    // SSIM constant
    const double k1 = 0.01, k2 = 0.03;
    const double C1 = (k1 * L) * (k1 * L);
    const double C2 = (k2 * L) * (k2 * L);

    // SSIM formula
    double numerator   = (2 * mu1 * mu2 + C1) * (2 * cov + C2);
    double denominator = (mu1*mu1 + mu2*mu2 + C1) * (var1 + var2 + C2);

    return numerator / denominator;
}

int main() {
    int sinogram[n][n] = {0};
    static int organize[n][n] = {0};

    ifstream inFile;
    
    // Read sinogram data
    inFile.open("sinogram_400x400.txt");	
    for (int row = 0; row < n; row++)
        for (int col = 0; col < n; col++)
            inFile >> sinogram[row][col];
    inFile.close();
    Display image2(sinogram, "sinogram");

   // Read the original image
    inFile.open("HumanBrain.txt");
    for (int row = 0; row < n; row++)
        for (int col = 0; col < n; col++)
            inFile >> organize[row][col];
    inFile.close();
    Display image10(organize, "HumanBrain");

    static double reconImage[N][N] = {0};
    static int finalImage[N][N];

    // ART reconstruction parameter optimization
    int iterations = 50;      // Increment the iteration count
    double lambda0 = 0.01;    // Reduce the initial step size
    double alpha = 0.0005;    // Decrease the attenuation factor
    
    cout << "Starting ART reconstruction..." << endl;
    ARTReconstruction(sinogram, reconImage, iterations, lambda0, alpha);
    
    cout << "Post-processing..." << endl;
    fix(reconImage, finalImage);
    Display image20(finalImage, "ReconstructedImage_Final_ART");

    // Compute evaluation metrics
    double sum = 0;
    double abs_sum = 0;
    int s = N * N;
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            double ERROR = organize[i][j] - finalImage[i][j];
            sum += ERROR * ERROR;
            abs_sum += abs(ERROR);
        }
    }
    
    double mse  = sum / s;
    double mae  = abs_sum / s;
    double rmse = sqrt(mse);
    const double MAX_I = 255.0;
    double psnr;
    if (mse <= 0.0) { 
        psnr = std::numeric_limits<double>::infinity();
    } else {
        psnr = 10.0 * std::log10((MAX_I * MAX_I) / mse);
    }

    double ssim = computeSSIM(organize, finalImage, N, N, MAX_I);
    
    cout << "\n=== Reconstruction Quality Metrics ===" << endl;
    cout << "MSE:  " << mse << endl;
    cout << "RMSE: " << rmse << endl;
    cout << "MAE:  " << mae << endl;
    cout << "PSNR: " << psnr << " dB" << endl;
    cout << "SSIM: " << ssim << endl;

    return 0;
}