#include <iostream>
#include <fstream>
using namespace std;

const int SIZE = 400;  

class Display
{ 
 public:  
    Display(int A[SIZE][SIZE], string fileName) {  
      fstream inFile("html.dat",ios::in);
      fstream outFile(fileName + ".html",ios::out);
      if (inFile.fail()) cout << "Cannot open file";    
      string line;

      outFile << "<html>" << endl;
      outFile << "<body>" << endl;
      outFile << "<canvas id='CANVAS_ID' width ='" << SIZE << "' height='" << SIZE << "'></canvas>" << endl;
      outFile << "<script>" << endl;
      outFile << "var com     =  4 ; //:4 components( RGBA )" << endl;
      outFile << "var wid     = " << SIZE << ";" << endl;
      outFile << "var hig     = " << SIZE << ";" << endl;

      getline(inFile,line);
      while (inFile) {
        outFile << line << endl;
        getline(inFile,line);
      }

      outFile << "var array = [";
      for (int row = 0; row < SIZE; row++)
          for (int col = 0; col < SIZE; col++)
            outFile << A[row][col] << ",";  

      outFile.seekp(-1,ios::cur);
      outFile << "];" << endl << endl;

      outFile << "	for( pi = 0; pi < array.length; pi++ ){" << endl;
      outFile << "		ci = pi * com;" << endl;
      outFile << "		u8a[ ci + 0 ] = array[pi]; //:Red" << endl;
      outFile << "		u8a[ ci + 1 ] = array[pi]; //:Green" << endl;
      outFile << "		u8a[ ci + 2 ] = array[pi]; //:Blue" << endl;
      outFile << "		u8a[ ci + 3 ] = 255; //:Alpha" << endl;
      outFile << "	}" << endl;
      outFile << "}" << endl;
      outFile << "</script>" << endl;
      outFile << "</body>" << endl;
      outFile << "</html>" << endl;

      inFile.close();
      outFile.close();  
    }
};
