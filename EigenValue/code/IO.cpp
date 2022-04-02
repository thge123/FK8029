#include "HEADER.h"

char *get_filename(char name, int filenumber){

    char *filename;
    filename = new char[10];
    filename[0] = name; filename[1] = '/';
    filename[2] = name; filename[3] = '0';
    filename[4] = '0';  filename[5] = '0';
    filename[6] = '.';  filename[7] = 'd';
    filename[8] = 'a';  filename[9] = 't';
    update_filename(filename,filenumber);

    return filename;
}

void update_filename(char *filename, int filenumber){
    
    char FileNum[4] = {'0','0','0'};
    if (filenumber > 999 || filenumber < 0){
        cout << "Maximum amount of files or invalid filename" << endl;
        exit(1);
    }
    string tmp = to_string(filenumber);
    if (filenumber<10)
        FileNum[2] = tmp[0];
    else
        if (filenumber<100){
            FileNum[1] = tmp[0];
            FileNum[2] = tmp[1];
        } else{
            FileNum[0] = tmp[0];
            FileNum[1] = tmp[1];
            FileNum[2] = tmp[2];
        }
    filename[3] = FileNum[0];
    filename[4] = FileNum[1];
    filename[5] = FileNum[2];
}

struct Files get_files(int filenumber){
    
    struct Files files;
    files.filenameA = get_filename('A',filenumber);
    files.filenameB = get_filename('B',filenumber);
    files.filenameC = get_filename('C',filenumber);
    files.filenameD = get_filename('D',filenumber);
    return files;
}

void update_filenames(struct Files files, int filenumber){

    update_filename(files.filenameA, filenumber);
    update_filename(files.filenameB, filenumber);
    update_filename(files.filenameC, filenumber);
    update_filename(files.filenameD, filenumber);
}

void write_vector(VectorXcd &X, ofstream &OutStream){

    int N = X.rows();
    for (int i=0; i<N; i++) OutStream << X(i) << ';';
    OutStream << endl;
}

void write_vector(VectorXd &X, ofstream &OutStream){

    int N = X.rows();
    for (int i=0; i<N; i++) OutStream << X(i) << ';';
    OutStream << endl;
}

void check_fail(ofstream &OutStream){

    if (OutStream.fail()){
        cout << "Failed to open output file" << endl;
        exit(1);
    }
}
