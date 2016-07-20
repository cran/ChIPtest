////////////////////////////
// ChIPtest.h Header File
//
class Xclass {
private:
int XLength;
public:
int GetXLength();
double *XVec;
Xclass(double *InputVec, int LengthInput); // Constructor
~Xclass(); // Destructor;
double *AutoCor();
};
double *ChIPtest(double *InputVec, int LengthInput);
