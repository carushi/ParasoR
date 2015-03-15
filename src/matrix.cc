#include "matrix.hh"

namespace Rfold {

void Matrix::Initialize() {
    try {
        stem = Mat(length+1, Vec(_constraint+2, -INF));
        stemend = Mat(length+1, Vec(_constraint+2, -INF));
        multi = Mat(length+1, Vec(_constraint+2, -INF));
        multibif = Mat(length+1, Vec(_constraint+2, -INF));
        multi1 = Mat(length+1, Vec(_constraint+2, -INF));
        multi2 = Mat(length+1, Vec(_constraint+2, -INF));
        _innerset = true;
    } catch (const std::bad_alloc& e) {
        cerr << "Memory allocation error : please try with more nodes" << endl;
        exit(1);
    }
}

void Matrix::PrintMat(const Mat& mat, const string& str) 
{
    for (LEN i = 0; i < (LEN)str.length(); i++)
        cout << str[i] << " ";
    cout << endl;
    for (Mat::const_iterator it = mat.begin(); it != mat.end(); it++)  
        PrintVec(*it, true);
}

void Matrix::Print(const string& str)
{
    cout << "---stem" << endl;
    PrintMat(stem, str);
    cout << "---stemend" << endl; 
    PrintMat(stemend, str);        
    cout << "---multi" << endl;
    PrintMat(multi, str);        
    cout << "---multiBif" << endl;
    PrintMat(multibif, str);        
    cout << "---multi1" << endl;
    PrintMat(multi1, str);        
    cout << "---multi2" << endl;
    PrintMat(multi2, str);        
    cout << "---douter" << endl;
    PrintMat(douter, str);
    cout << "---outer" << endl;
    PrintVec(outer, true);
}

}