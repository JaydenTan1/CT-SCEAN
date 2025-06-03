//Author: Zhijie Tan
//Time： 6/3/2025
#include <fstream>
#include <iostream>
#include <vector>
#include "Display1.h"
#include <math.h>
#include <iomanip>
#include <algorithm>
#include <cmath>
using namespace std;

const int n = 400;
const int N = 400; 
int newsin[N][N];
int newsin2[N][N];
int result[N][N];  //Sinogram

const double PI = 3.14159265358979323846;


//HannWindow
void createHannWindow(double window[], int size) {
    for (int i = 0; i < size; i++) {
        window[i] = 0.5 * (1 - cos(2 * PI * i / (size - 1)));
    }
}

//Create a ramp filter
void createRampFilter(double ramp[], int n) {
    for (int i = 0; i < n; ++i) {
        //The value of k is calculated based on the distance between the current position and the center, which is positive to the right and negative to the left of the center
        int k = (i <= n / 2) ? i : i - n;
        ramp[i] = std::abs((double)k);
    }

    // Enlarge the entire filter amplitude (2~4 times recommended)
    for (int i = 0; i < n; ++i) {
        ramp[i] *= 2.0;
    }
}


void applyRampFilterToSinogram_RowWise(int input[N][N], double output[N][N]) {
    double ramp[N];
    createRampFilter(ramp, N);

    for (int row = 0; row < N; ++row) {
        double realF[N] = {0.0}, imagF[N] = {0.0};

        // DFT
        for (int k = 0; k < N; ++k) {
            for (int n = 0; n < N; ++n) {
                double angle = -2.0 * PI * k * n / N;
                realF[k] += input[row][n] * cos(angle);
                imagF[k] += input[row][n] * sin(angle);
            }
        }

        // Apply ramp filter
        for (int k = 0; k < N; ++k) {
            realF[k] *= ramp[k];
            imagF[k] *= ramp[k];
        }

        // IDFT
        for (int n = 0; n < N; ++n) {
            double val = 0.0;
            for (int k = 0; k < N; ++k) {
                double angle = 2.0 * PI * k * n / N;
                val += realF[k] * cos(angle) - imagF[k] * sin(angle);
            }
            val /= N;
            output[row][n] = val;  // The negative values are not cropped to preserve the true grayscale structure

        }
    }
}

void convertDisplayToProjections(int filteredSinogram[N][N], double* projections) {
    for (int angle = 0; angle < 180; angle++) {
        for (int t = 0; t < N; t++) {
            // Convert from display value (where 128 is zero) back to mathematical value
            // Scale factor should be adjusted based on your filtering approach
            double scaleFactor = 600.0 / 127.0; // Adjust this based on your filtering strength
            projections[angle * N + t] = (filteredSinogram[angle][t] - 128.0) * scaleFactor;
        }
    }
}


void improvedBackProjection(int filteredSinogram[N][N], int output[N][N]) {
    // 1. Convert filtered sinogram display values back to actual projection values
    double* projections = new double[180 * N];
    convertDisplayToProjections(filteredSinogram, projections);
    
    // 2. Clear output image
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            output[i][j] = 0;
        }
    }
    
    // 3. Perform backprojection with bilinear interpolation
    double centerX = N / 2.0;
    double centerY = N / 2.0;
    
    for (int angle = 0; angle < 180; angle++) {
        double theta = angle * PI / 180.0;
        double cosTheta = cos(theta);
        double sinTheta = sin(theta);
        
        for (int y = 0; y < N; y++) {
            for (int x = 0; x < N; x++) {
                // Calculate exact position in sinogram
                double xCentered = x - centerX;
                double yCentered = y - centerY;
                double t = xCentered * cosTheta + yCentered * sinTheta + centerX;
                
                if (t >= 0 && t < N - 1) {
                    // Bilinear interpolation for better quality
                    int t0 = static_cast<int>(floor(t));
                    int t1 = t0 + 1;
                    double weight = t - t0;
                    
                    double val0 = projections[angle * N + t0];
                    double val1 = projections[angle * N + t1];
                    
                    // Interpolate
                    double contribution = val0 * (1.0 - weight) + val1 * weight;
                    
                    // Add contribution to output
                    output[y][x] += contribution;
                }
            }
        }
    }
    
    // 4. Find Min and Max for normalization
    int minVal = output[0][0];
    int maxVal = output[0][0];
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            minVal = min(minVal, output[i][j]);
            maxVal = max(maxVal, output[i][j]);
        }
    }
    
    // 5. Apply enhanced contrast with gamma correction
    double gamma = 1.8; // Values < 1 enhance contrast in darker regions  1.8
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            // Normalize to [0,1]
            double normalized = static_cast<double>(output[i][j] - minVal) / (maxVal - minVal);
            
            // Apply gamma correction
            normalized = pow(normalized, gamma);
            
            // Convert to [0,255]
            output[i][j] = static_cast<int>(normalized * 255.0);
        }
    }
    
    // Apply edge enhancement for better detail visibility
    int* tempImage = new int[N * N];
    
    // Copy current image
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            tempImage[i * N + j] = output[i][j];
        }
    }
    
    // Apply unsharp mask for edge enhancement
    int radius = 2;
    double amount = 0.5; // Strength of sharpening
    
    for (int i = radius; i < N - radius; i++) {
        for (int j = radius; j < N - radius; j++) {
            int sum = 0;
            for (int di = -radius; di <= radius; di++) {
                for (int dj = -radius; dj <= radius; dj++) {
                    sum += tempImage[(i + di) * N + (j + dj)];
                }
            }
            
            // Calculate average (blurred value)
            int count = (2 * radius + 1) * (2 * radius + 1);
            int blurred = sum / count;
            
            // Apply unsharp mask formula: original + (original - blurred) * amount
            int original = tempImage[i * N + j];
            int sharpened = original + static_cast<int>((original - blurred) * amount);
            
            // Clamp to valid range
            output[i][j] = max(0, min(255, sharpened));
        }
    }
    
    delete[] tempImage;
    delete[] projections;
}

int reconstruction[N][N] = {0};

int main() {
    
	ifstream inFile;
    
	inFile.open("sinogram_400x400.txt");
		
	for (int row = 0; row < n; row++)
		for (int col = 0; col < n; col++)
			inFile >> result[row][col];
	
	inFile.close();
    Display image(result, "Orginal");
    
    double filteredSinogram[N][N];
    applyRampFilterToSinogram_RowWise(result, filteredSinogram);

    double minVal = filteredSinogram[0][0];
    double maxVal = filteredSinogram[0][0];
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            minVal = min(minVal, filteredSinogram[i][j]);
            maxVal = max(maxVal, filteredSinogram[i][j]);
        }

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j){
            newsin2[i][j] = static_cast<int>((filteredSinogram[i][j] - minVal) / (maxVal - minVal) * 255.0 );
        }


    cout << "Filtered Sinogram: min = " << minVal << ", max = " << maxVal << std::endl;
    Display image2(newsin2, "Filtered Sinogram");

    // Reconstructed image
    improvedBackProjection(newsin2, reconstruction);

    Display image3(reconstruction, "Reconstructed Image_f4");
}
