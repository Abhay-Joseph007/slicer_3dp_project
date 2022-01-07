#include<iostream>
#include<fstream>
#include<math.h>

#define MAX_READ_CHUNK 80
#define FACET_READ_SIZE 12

using std::cout;
using std::endl;
using std::cin;
using std::ofstream;
using std::ifstream;
using std::ios;
using std::hex;
using std::dec;


void Facet_BintoAsc(float values[], ofstream &writefile, ifstream &readfile)
{
    readfile.read((char*)values, FACET_READ_SIZE*sizeof(float));
    writefile << "Normal Vector: " << values[0] << ", " << values[1] << ", " << values[2] << endl;
    writefile << "Vertex 1: " << values[3] << ", " << values[4] << ", " << values[5] << endl;
    writefile << "Vertex 2: " << values[6] << ", " << values[7] << ", " << values[8] << endl;
    writefile << "Vertex 3: " << values[9] << ", " << values[10] << ", " << values[11] << endl;
}


int main(int argc, char** argv)
{
    int file_size = 0;                              //to store total file size
    int no_of_triangles = 0;                        //to store no of triangles
    char var[MAX_READ_CHUNK];                       //to read and store header from input file
    float vector_values[FACET_READ_SIZE];             //to store the vector read from input file


    //extract name of stl file from input stream
    std::string output_file_name = argv[argc-1];       //stores name of input file
    output_file_name += "_parsed.txt";                  //name of parsed output file

    //open files
    ofstream FileWrite(output_file_name,ios::out | ios::binary | ios::trunc );
    ifstream FileRead(argv[argc-1], ios::binary | ios::ate);

    //file_size stores size of total stl file read
    file_size = ((int)FileRead.tellg());
    FileRead.seekg(0);
    cout << "Size of input file is: " << file_size <<endl;

    //Print number of triangles in the stl file
    no_of_triangles = (file_size-84)/50;
    cout << "Number of Triangles in the input file is: " << no_of_triangles << endl;


    //Read 80 bytes of header and store in file
    FileRead.read(var, 80);
    FileWrite << "File Header is: "<< endl;
    FileWrite << var << endl;

    //Write to parsed file the number of Triangles
    FileWrite << "No of Triangles in File is: ";
    FileWrite << no_of_triangles << endl;
    FileRead.seekg(84);


    //Reading each facet:
    for(int i = 0; i<no_of_triangles; i++)
    {
        FileWrite << "Facet " << i+1 << ":" << endl;
        Facet_BintoAsc(vector_values, FileWrite, FileRead);
        FileWrite << endl;
        FileRead.read(var,2);
    }
    //int pos = ((int)FileRead.tellg());
    //cout << pos << endl;


    //close file descriptors
    FileWrite.close();
    FileRead.close();

    return 0;
}
