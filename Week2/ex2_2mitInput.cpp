#include <iostream>
using namespace std;

// void printArray(float arr[][], int N)
// {
// 	for (int i = 0; i < N; ++i) {
// 		for (int j = 0; j < N; ++j) {
// 			cout << arr[i][j]<<" ";
// 		}
// 		cout << endl;
// 	}
// }
void printArray(float** arr, int N)
{
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			cout << arr[i][j]<<" ";
		}
		cout << endl;
	}
}
void printVector(float* arr, int N)
{
	for (int j = 0; j < N; ++j) {
		cout << arr[j]<<" ";
	}
	cout << endl;
}
float* backward_substitution(float** input, float* vector, int N){
    static float out[6]; //geht das so?
    for(int x = N-1; x >= 0; x--){
        float subtraction = 0;
        for(int i = x+1; i < N; i++){
            subtraction += input[x][i]*out[i]; 
        }
        out[x] = (vector[x]-subtraction)/input[x][x];
    }
    return out;
}
void GaussJordan_wo_pivot(float** input, float* vector, int N){
    for(int a=0; a < N-1; a++){
        float factor = input[a][a];
        for(int i=0;i<=N-1; i++){
            input[a][i] = input[a][i] / factor;
        }
        vector[a] = vector[a] / factor;
        for(int c=a+1; c<N; c++){            
            float subtraction_factor = input[c][a];
            for(int i=0;i<N; i++){
                input[c][i] = input[c][i]-subtraction_factor*input[a][i];
            }
            vector[c] = vector[c]-subtraction_factor*vector[a];
        }
    }
}
int index_largest_abs(float** vector, int a){
    int index = a;
    float max = 0;
    for(int i=a;i<5;i++){
        if(vector[i][a] < -max){
            index = i;
            max = -vector[i][a];
        }
        if(vector[i][a] > max){
            index = i;
            max = vector[i][a];
        }
    }
    return index;
}

void GaussJordan_w_pivot(float** input, float* vector, int N){
    for(int a=0; a < N-1; a++){
        // pivoting
        int maxindex = index_largest_abs(input, a);
        // cout<<"maxindex: "<<maxindex<<"\n";
        float tmprow;
        for(int l=0;l<N;l++){
            tmprow = input[a][l];
            input[a][l] = input[maxindex][l];
            input[maxindex][l] = tmprow;
        }
        float tmpvec = vector[a];
        vector[a] = vector[maxindex];
        vector[maxindex] = tmpvec;
        // GaussJordan
        float factor = input[a][a];
        for(int i=0;i<=N-1; i++){ 
            input[a][i] = input[a][i] / factor;
        }
        vector[a] = vector[a] / factor;
        // HIER IRGENDWO KICKT ES MICH
        for(int c=a+1; c<N; c++){            
            float subtraction_factor = input[c][a];
            for(int i=0;i<N; i++){
                input[c][i] = input[c][i]-subtraction_factor*input[a][i];
            }
            vector[c] = vector[c]-subtraction_factor*vector[a];
        }
    }
}

float* matrix_vector_mult(float** input, float* vector, float* out, int N){
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            out[i] += input[i][j]*vector[j];
        }
    }
    return out;
}

float** get_Matrix(){
    float** out = {};
    return out;
}

int main(){
    // versuch 2 mit anderer Größe
    int r = 0;
    cout<<"Select the Dimension N of your Matrix: ";
    cin>>r;
    cout<<"next you have to type in the matrix, value by value: \n";
    const int N = r;
    float eingabe[N][N];
    for(int r = 0;r<N;r++){
        for(int c = 0;c<N;c++){
            cout<<"Matrix["<<r<<"]["<<c<<"]= ";
            cin>>eingabe[r][c];
        }
    }    
    // Das funktioniert nur für 5D... :(
    // Ändern
    float* solution_rows[N] = {eingabe[0], eingabe[1], eingabe[2], eingabe[3], eingabe[4]};
    float** matrixpointer = solution_rows;
    float ursprungsmatrix[N][N];
    for(int i = 0;i<N;i++){
        for(int j = 0;j<N;j++){
            ursprungsmatrix[i][j] = eingabe[i][j];}}
    float* solution_rows2[N] = {ursprungsmatrix[0], ursprungsmatrix[1], ursprungsmatrix[2], ursprungsmatrix[3], ursprungsmatrix[4]};
    float** testpointer = solution_rows2;
    cout<<"Your full matrix is: \n";
    printArray(matrixpointer, N);
    cout<<"choose your Outputvector b0: \n";
    float b0[] = {1,-1,1,-1,1,-1};
    for(int c = 0;c<N;c++){
            cout<<"Vector["<<c<<"]= ";
            cin>>b0[c];
    }
    cout<<"Your Outputvector is: ";
    printVector(b0, N);

    //Ab hier beginnt die eigentliche Rechnung
    GaussJordan_w_pivot(matrixpointer, b0, N);
    cout<<"the Triangular Matrix for the calculation is: \n";
    printArray(matrixpointer, N);
    cout<<"\n";
    float* output = backward_substitution(matrixpointer, b0, N);
    cout<<"der gesuchte Vector x ist : \n";
    printVector(output, N);
    cout<<"\n";
    float newVector[N];
    matrix_vector_mult(testpointer, output, newVector, N);
    cout<<"Test durch Rueckrechnen mit der Eingabematrix: \n";
    printVector(newVector, N);
    
    while(1){
    }
    return 0;
}