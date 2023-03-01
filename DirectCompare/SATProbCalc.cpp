#include <fstream>
#include <iostream>  
#include <bitset>  
#include <iomanip>
#include <math.h>
#include "SATProbCalc.hpp"
#include <vector>
#include <stdlib.h>
#include <iterator>

// no error 
// g++ -c -std=c++11 -O3 -fpic SATProbCalc.cpp SATProbCalc.hpp
// g++ -shared -std=c++11 -O3 -o lib/SATProbCalc.so SATProbCalc.o

//compile using 
//g++ -c -Wall -std=c++11 -Werror -O3 -fpic SATProbCalc.cpp SATProbCalc.hpp
//then 
//g++ -shared -std=c++11 -Wall -Werror -O3 -o lib/SATProbCalc.so SATProbCalc.o

//with parallell
//compile using 
//g++ -c -Wall -Werror -O3 -fopenmp -fpic SATProbCalc.cpp SATProbCalc.hpp
//then 
//g++ -shared -Wall -Werror -O3 -fopenmp -o lib/SATProbCalc.so SATProbCalc.o

using namespace std;  

void printVector(const std::vector<char> &n){

    cout << "State is :" << endl;

    for (int i=0; i<n.size(); i++)
       cout << ' ' << n[i];
    std::cout << '\n';

}

int SolveSATbpp(int n, double * SAT, int c, int i){
    //number of variables, SAT instance, number of clauses, number of iterations
//     srand(time());
    srand(time(NULL));
    vector<char> s(n,'0');
    
//     cout << "capacity : " << s.capacity() << endl;
    
    for(int w = 0; w<n; w++){
        s[w] = '0' + (rand()%2);
    }

    bool solved = true;
//     unsigned long long int toPrint = 100000000;
//     unsigned long long int dP = 100000000;
//     int maxC = 0;
    int numIt = 0;

    for (int it = 0;it < i;it++){
        solved = true;
        for(int t = 0; t < c ; t++){
//             numIt++;
            if( ((s[(abs((int)SAT[3*t])-1)]   == '1' && (SAT[3*t]>0))   || (s[(abs((int)SAT[3*t])-1)]   == '0' && (SAT[3*t]<0)) ) 
            &&  ((s[(abs((int)SAT[3*t+1])-1)] == '1' && (SAT[3*t+1]>0)) || (s[(abs((int)SAT[3*t+1])-1)] == '0' && (SAT[3*t+1]<0)) )
            &&  ((s[(abs((int)SAT[3*t+2])-1)] == '1' && (SAT[3*t+2]>0)) || (s[(abs((int)SAT[3*t+2])-1)] == '0' && (SAT[3*t+2]<0)) ) ){
                s[(abs((int)SAT[3*t])-1)] = '0' + (rand()%2);
                s[(abs((int)SAT[3*t+1])-1)] = '0' + (rand()%2);
                s[(abs((int)SAT[3*t+2])-1)] = '0' + (rand()%2);
                solved = false;
//                 it++;
            }

//             it++;
        }

        if(solved){
//             cout << "Solved! Number of iterations = " << it << endl;
//             printVector(s);
            numIt = it;
            return numIt;
//             break;
        }
    }
//     cout << "exited" << endl;
//     printVector(s);
    return numIt;
}

int Shoning(int n, double * SAT, int c, int i){
    //number of variables, SAT instance, number of clauses, number of iterations
//     srand(time());
    srand(time(NULL));
    vector<char> s(n,'0');
    
//     cout << "capacity : " << s.capacity() << endl;
    
    for(int w = 0; w<n; w++){
        s[w] = '0' + (rand()%2);
    }

    bool solved = true;
//     unsigned long long int toPrint = 100000000;
//     unsigned long long int dP = 100000000;
//     int maxC = 0;
    int numIt = 0;
    int bitToFlip = 0;

    for (int it = 0;it < i;it++){
        solved = true;
        for(int t = 0; t < c ; t++){
//             numIt++;
            if( ((s[(abs((int)SAT[3*t])-1)]   == '1' && (SAT[3*t]>0))   || (s[(abs((int)SAT[3*t])-1)]   == '0' && (SAT[3*t]<0)) ) 
            &&  ((s[(abs((int)SAT[3*t+1])-1)] == '1' && (SAT[3*t+1]>0)) || (s[(abs((int)SAT[3*t+1])-1)] == '0' && (SAT[3*t+1]<0)) )
            &&  ((s[(abs((int)SAT[3*t+2])-1)] == '1' && (SAT[3*t+2]>0)) || (s[(abs((int)SAT[3*t+2])-1)] == '0' && (SAT[3*t+2]<0)) ) ){
                bitToFlip = rand()%3;
                s[(abs((int)SAT[3*t+bitToFlip])-1)] = '0' + (rand()%2);
                solved = false;
//                 it++;
            }

//             it++;
        }

        if(solved){
//             cout << "Solved! Number of iterations = " << it << endl;
//             printVector(s);
            numIt = it;
            return numIt;
//             break;
        }
    }
//     cout << "exited" << endl;
//     printVector(s);
    return numIt;
}

double SolveSAT(int n, double * SAT, int c, int i, int v){//number of variables, SAT instance, number of clauses, number of iterations
    for(int b = 0; b < v; b++){
        cout << SAT[b] << " ";
    }
    cout << endl;
    
    //number of states
    int N = 1<<n;

    //Generate probability array for states
    double p[N];

    // double finalProb[i];

    //counts of how many times state was removed
    int counts[N];

    double sum;
    int count;

    for(int a = 0; a<N; a++){
        p[a] = 1/((float)N);
        counts[a] = 0;
    }

    for(int k = 0; k < i ; k++){
    //loop through all the clauses in the SAT instance
    for(int t = 0; t < c ; t++){
        count = 0;
        sum = 0;

        //get the variables the clause is interested in
        int seive = (1<<(abs((int)SAT[3*t])-1)) + (1<<((int)abs(SAT[3*t + 1])-1)) + (1<<((int)abs(SAT[3*t + 2])-1)); 
        //get the values in the clause to xor
        int gate = 0;
        if(SAT[3*t] > 0)
            gate += 1<<(abs((int)(SAT[3*t]))-1);
        if(SAT[3*t+1] > 0)
            gate += 1<<(abs((int)(SAT[3*t+1]))-1);
        if(SAT[3*t+2] > 0)
            gate += 1<<(abs((int)(SAT[3*t+2]))-1);

        for(int j = 0; j < N; j++){
            int mult = seive&j; 
            int ore = gate^mult;

            //if this clause shoudl be rotated out
            if(ore == 0){
                //get probability from state
                sum += p[j];
                p[j] = 0; 
                count++;
                counts[j]++;
            }
        }
        // cout << "number of states rotated = " << count << endl;
        // cout << "total prob in unsat states = " << sum << endl;
        sum = sum/N;
        // cout << "amount added to each state = " << sum << endl;
        for(int j = 0; j < N; j++){
            p[j] += sum;
        }
    }
    }
    
    
    double total = 0;
    int numSat = 0;

    cout << endl;

    for(int j = 0; j < N; j++){
        if(counts[j] == 0){
            total += p[j];
            numSat++;
        }
    }
    

    cout<< "There was a " << total << " probability of measuring one of the " << numSat << " satisfiable states" << endl;
    // for(int j = 0; j < N; j++){

    //     cout<< p[j]<<", ";
 
    // }
    // cout << endl;

    // finalProb[k] = total;
    // }
    return total;//finalProb;
}

// double SolveSATfast(int n, double * SAT, int c, int i, int v){//number of variables, SAT instance, number of clauses, number of iterations
//     for(int b = 0; b < v; b++){
//         cout << SAT[b] << " ";
//     }
//     cout << endl;

//     double sum;
//     int count;

//     for(int a = 0; a<N; a++){
//         p[a] = 1/((float)N);
//         counts[a] = 0;
//     }

//     bool solved = 1;

//     int start = 0;//start in the 0 state for now

//     for(int k = 0; k < i ; k++){\\loop through the number of iterations
//     //loop through all the clauses in the SAT instance
//     for(int t = 0; t < c ; t++){
//         count = 0;
//         sum = 0;

//         //get the variables the clause is interested in
//         int seive = (1<<(abs((int)SAT[3*t])-1)) + (1<<((int)abs(SAT[3*t + 1])-1)) + (1<<((int)abs(SAT[3*t + 2])-1)); 
//         //get the values in the clause to xor
//         int gate = 0;
//         if(SAT[3*t] > 0)
//             gate += 1<<(abs((int)(SAT[3*t]))-1);
//         if(SAT[3*t+1] > 0)
//             gate += 1<<(abs((int)(SAT[3*t+1]))-1);
//         if(SAT[3*t+2] > 0)
//             gate += 1<<(abs((int)(SAT[3*t+2]))-1);

//         for(int j = 0; j < N; j++){
//             int mult = seive&j; 
//             int ore = gate^mult;

//             //if this clause shoudl be rotated out
//             if(ore == 0){
//                 //get probability from state
//                 sum += p[j];
//                 p[j] = 0; 
//                 count++;
//                 counts[j]++;
//             }
//         }
//     }
//         if(){

//         }
//     }
//     }
    
//     cout<< "There was a " << total << " probability of measuring one of the " << numSat << " satisfiable states" << endl;

//     return total;//finalProb;
// }


double * SolveSATactual(int n, double * SAT, int c, int i, int v){//number of variables, SAT instance, number of clauses, number of iterations
    
    for(int b = 0; b < v; b++){
        cout << SAT[b] << " ";
    }
    cout << endl;
    
    //number of states
    int N = 1<<n;

    //Generate probability array for states
    double * p = new double[N];

    // double finalProb[i];

    //counts of how many times state was removed
    int * counts = new int[N];

    double sum;
    int count;


    unsigned int mult; 
    unsigned int an;

    unsigned int mult2; 
    unsigned int ore;

    //get the variables the clause is interested in
    unsigned int seive;// = (1<<(abs((int)SAT[3*t])-1)) + (1<<((int)abs(SAT[3*t + 1])-1)) + (1<<((int)abs(SAT[3*t + 2])-1)); 
    //get the values in the clause to xor
    unsigned int gate;// = 0;

    double * tots = new double[i];

    double total = 0;
    int numSat = 0;

//    int temp = 0;
    
    for(int a = 0; a<N; a++){
        p[a] = 1/((float)N);
        counts[a] = 0;
    }
    for(int k = 0; k < i ; k++){
    //loop through all the clauses in the SAT instance
    for(int t = 0; t < c ; t++){
        count = 0;
        sum = 0;

//        Shaker mode
//        temp = t;
//        if(t%2 == 0){
//            t = c - t - 1;
//        }

        //get the variables the clause is interested in
        seive = (1<<(abs((int)SAT[3*t])-1)) + (1<<((int)abs(SAT[3*t + 1])-1)) + (1<<((int)abs(SAT[3*t + 2])-1)); 
        //get the values in the clause to xor
        gate = 0;
        if(SAT[3*t] > 0)
            gate += 1<<(abs((int)(SAT[3*t]))-1);
        if(SAT[3*t+1] > 0)
            gate += 1<<(abs((int)(SAT[3*t+1]))-1);
        if(SAT[3*t+2] > 0)
            gate += 1<<(abs((int)(SAT[3*t+2]))-1);

        for(int j = 0; j < N; j++){//loop through all states and if state is not unsat then add 1/8th prob of corresponding unsat state

            mult = j&(~seive); 
            an = mult^gate;

            mult2 = seive&j; 
            ore = gate^mult2;
            // cout << an  << " < an or > " << ore << endl;
            if(ore != 0){
                // cout << p[j] << "before"<< endl;
                p[j] += p[an]/8;
                // cout << p[j] << "after"<< endl;
            }else{

                // cout << "hellooooo"<< p[j] << endl;
            }
        }
        for(int j = 0; j < N; j++){//

            //get porbability from state 
            mult = seive&j; 
            ore = gate^mult;
            //if this clause shoudl be rotated out
            if(ore == 0){
                //get probability from state
                // double ne = p[j]/8;
                p[j] = p[j]/8;//ne;
                // cout << p[j] << endl;
                count++;
                counts[j]++;
            }
        }
//        shaker mode
//        t = temp;
    }
    total = 0;
    numSat = 0;

//    cout << endl;
    for(int j = 0; j < N; j++){
        if(counts[j] == 0){
            // cout << j<< endl;
            total += p[j];
            numSat++;
        }
    }
//    cout<< "There was a " << total << " probability of measuring one of the " << numSat << " satisfiable states for iteration : " <<  k <<  endl;

    tots[k] = total;
//    if(total > 0.5)
//        return tots;
    }
    
    
    // total = 0;
    // numSat = 0;

    // cout << endl;

    // for(int j = 0; j < N; j++){
    //     if(counts[j] == 0){
    //         cout << j<< endl;
    //         total += p[j];
    //         numSat++;
    //     }
    // }
    

    cout<< "There was a " << total << " probability of measuring one of the " << numSat << " satisfiable states" << endl;
    // for(int j = 0; j < N; j++){

    //     cout<< p[j]<<", ";
 
    // }
    // cout << endl;

    // finalProb[k] = total;
    // }
    return tots;//finalProb;
}

int SolveSATiterations(int n, double * SAT, int c, int i, int v){//number of variables, SAT instance, number of clauses, number of iterations
    
    for(int b = 0; b < v; b++){
        cout << SAT[b] << " ";
    }
    cout << endl;
    
    //number of states
    int N = 1<<n;

    //Generate probability array for states
    double * p = new double[N];

    // double finalProb[i];

    //counts of how many times state was removed
    int * counts = new int[N];

    double sum;
    int count;


    unsigned int mult; 
    unsigned int an;

    unsigned int mult2; 
    unsigned int ore;

    //get the variables the clause is interested in
    unsigned int seive;// = (1<<(abs((int)SAT[3*t])-1)) + (1<<((int)abs(SAT[3*t + 1])-1)) + (1<<((int)abs(SAT[3*t + 2])-1)); 
    //get the values in the clause to xor
    unsigned int gate;// = 0;

    double * tots = new double[i];

    double total = 0;
    int numSat = 0;
    #pragma omp parallel
    for(int a = 0; a<N; a++){
        p[a] = 1/((float)N);
        counts[a] = 0;
    }
    for(int k = 0; k < i ; k++){
    //loop through all the clauses in the SAT instance
    for(int t = 0; t < c ; t++){
        count = 0;
        sum = 0;

        //get the variables the clause is interested in
        seive = (1<<(abs((int)SAT[3*t])-1)) + (1<<((int)abs(SAT[3*t + 1])-1)) + (1<<((int)abs(SAT[3*t + 2])-1)); 
        //get the values in the clause to xor
        gate = 0;
        if(SAT[3*t] > 0)
            gate += 1<<(abs((int)(SAT[3*t]))-1);
        if(SAT[3*t+1] > 0)
            gate += 1<<(abs((int)(SAT[3*t+1]))-1);
        if(SAT[3*t+2] > 0)
            gate += 1<<(abs((int)(SAT[3*t+2]))-1);

        for(int j = 0; j < N; j++){//loop through all states and if state is not unsat then add 1/8th prob of corresponding unsat state

            mult = j&(~seive); 
            an = mult^gate;

            mult2 = seive&j; 
            ore = gate^mult2;
            // cout << an  << " < an or > " << ore << endl;
            if(ore != 0){
                // cout << p[j] << "before"<< endl;
                p[j] += p[an]/8;
                // cout << p[j] << "after"<< endl;
            }
        }
        for(int j = 0; j < N; j++){//

            //get porbability from state 
            mult = seive&j; 
            ore = gate^mult;
            //if this clause shoudl be rotated out
            if(ore == 0){
                //get probability from state
                // double ne = p[j]/8;
                p[j] = p[j]/8;//ne;
                // cout << p[j] << endl;
                count++;
                counts[j]++;
            }
        }
        
    }
    total = 0;
    numSat = 0;

    cout << endl;
    for(int j = 0; j < N; j++){
        if(counts[j] == 0){
            // cout << j<< endl;
            total += p[j];
            numSat++;
        }
    }
    cout<< "There was a " << total << " probability of measuring one of the " << numSat << " satisfiable states for iteration : " <<  k <<  endl;

    tots[k] = total;
    if(total > 0.125)
        return k;
    }
    
    
    // total = 0;
    // numSat = 0;

    // cout << endl;

    // for(int j = 0; j < N; j++){
    //     if(counts[j] == 0){
    //         cout << j<< endl;
    //         total += p[j];
    //         numSat++;
    //     }
    // }
    

    cout<< "There was a " << total << " probability of measuring one of the " << numSat << " satisfiable states" << endl;
    // for(int j = 0; j < N; j++){

    //     cout<< p[j]<<", ";
 
    // }
    // cout << endl;

    // finalProb[k] = total;
    // }
    return n*n;//finalProb;
}

double * SolveSATgates(int n, double * SAT, int c, int i, int v){//number of variables, SAT instance, number of clauses, number of iterations

    for(int b = 0; b < v; b++){
        cout << SAT[b] << " ";
    }
    cout << endl;

    //number of states
    int N = 1<<n;

    //Generate probability array for states
    double * p = new double[N];

    // double finalProb[i];

    //counts of how many times state was removed
    int * counts = new int[N];

    double sum;
    int count;


    unsigned int mult;
    unsigned int an;

    unsigned int mult2;
    unsigned int ore;

    //get the variables the clause is interested in
    unsigned int seive;// = (1<<(abs((int)SAT[3*t])-1)) + (1<<((int)abs(SAT[3*t + 1])-1)) + (1<<((int)abs(SAT[3*t + 2])-1));
    //get the values in the clause to xor
    unsigned int gate;// = 0;

//    double * tots = new double[i];
    double * tots = new double[c];

    double total = 0;
    int numSat = 0;

    for(int a = 0; a<N; a++){
        p[a] = 1/((float)N);
        counts[a] = 0;
    }
//    for(int k = 0; k < i ; k++){
    //loop through all the clauses in the SAT instance
    for(int t = 0; t < c ; t++){
    //get the variables the clause is interested in
        seive = (1<<(abs((int)SAT[3*t])-1)) + (1<<((int)abs(SAT[3*t + 1])-1)) + (1<<((int)abs(SAT[3*t + 2])-1));
        //get the values in the clause to xor
        gate = 0;
        if(SAT[3*t] > 0)
            gate += 1<<(abs((int)(SAT[3*t]))-1);
        if(SAT[3*t+1] > 0)
            gate += 1<<(abs((int)(SAT[3*t+1]))-1);
        if(SAT[3*t+2] > 0)
            gate += 1<<(abs((int)(SAT[3*t+2]))-1);

    for(int k = 0; k < i ; k++){
        count = 0;
        sum = 0;


        for(int j = 0; j < N; j++){//loop through all states and if state is not unsat then add 1/8th prob of corresponding unsat state

            mult = j&(~seive);
            an = mult^gate;

            mult2 = seive&j;
            ore = gate^mult2;
            // cout << an  << " < an or > " << ore << endl;
            if(ore != 0){
                // cout << p[j] << "before"<< endl;
                p[j] += p[an]/8;
                // cout << p[j] << "after"<< endl;
            }else{

                // cout << "hellooooo"<< p[j] << endl;
            }
        }
        for(int j = 0; j < N; j++){//

            //get porbability from state
            mult = seive&j;
            ore = gate^mult;
            //if this clause shoudl be rotated out
            if(ore == 0){
                //get probability from state
                // double ne = p[j]/8;
                p[j] = p[j]/8;//ne;
                // cout << p[j] << endl;
                count++;
                counts[j]++;
            }
        }

    }
    total = 0;
    numSat = 0;

//    cout << endl;
    for(int j = 0; j < N; j++){
        if(counts[j] == 0){
            // cout << j<< endl;
            total += p[j];
            numSat++;
        }
    }
    cout<< "There was a " << total << " probability of measuring one of the " << numSat << " satisfiable states after clause : " <<  t <<  endl;

//    tots[k] = total;
    tots[t] = total;
//    if(total > 0.5)
//        return tots;
    }


    // total = 0;
    // numSat = 0;

    // cout << endl;

    // for(int j = 0; j < N; j++){
    //     if(counts[j] == 0){
    //         cout << j<< endl;
    //         total += p[j];
    //         numSat++;
    //     }
    // }


    cout<< "There was a " << total << " probability of measuring one of the " << numSat << " satisfiable states" << endl;
    // for(int j = 0; j < N; j++){

    //     cout<< p[j]<<", ";

    // }
    // cout << endl;

    // finalProb[k] = total;
    // }
    return tots;//finalProb;
}

// int * genSAT(int n){
//     const int v = 3*(7+3*(n-3));
//     int SAT[v];
//     static int startSAT[] = {-1,-2,-3,-1,2,3,1,-2,3,1,2,-3,-1,-2,3,-1,2,-3,1,-2,-3};

//     for(int b = 0; b < 7*3; b++){
//         SAT[b] = startSAT[b];
//     }
//     for(int a = 4; a <= n; a++){
//         SAT[7*3 + 9*(a-4)] = -2;
//         SAT[7*3 + 9*(a-4)+1] = -3;
//         SAT[7*3 + 9*(a-4)+2] = n;
//         SAT[8*3 + 9*(a-4)] = -2;
//         SAT[8*3 + 9*(a-4)+1] = -3;
//         SAT[8*3 + 9*(a-4)+2] = -n;
//         SAT[9*3 + 9*(a-4)] = 2;
//         SAT[9*3 + 9*(a-4)+1] = 3;
//         SAT[9*3 + 9*(a-4)+2] = n;
//     }
//     for(int b = 0; b < v; b++){
//         cout << SAT[b] << " ";
//     }
//     cout << endl;
//     return SAT;
// }

// int main(){
    
//     int n = 4;
//     int numClauses = 7 + 3*(n-3);
//     // numClauses = 2;
//     const int v = 3*(7+3*(n-3));
//     int SAT[v];
//     //{-1,-2,-3,-4,-2,-3,1,5,3,4,2,3};
//     //{-1,-2,-3,-1,2,3,1,-2,3,1,2,-3,-1,-2,3,-1,2,-3,1,-2,-3};
//     //{1,2,3,4,5,6};
//     static int startSAT[] = {-1,-2,-3,-1,2,3,1,-2,3,1,2,-3,-1,-2,3,-1,2,-3,1,-2,-3};

//     for(int b = 0; b < 7*3; b++){
//         SAT[b] = startSAT[b];
//     }
//     for(int a = 4; a <= n; a++){
//         SAT[7*3 + 9*(a-4)] = -2;
//         SAT[7*3 + 9*(a-4)+1] = -3;
//         SAT[7*3 + 9*(a-4)+2] = a;
//         SAT[8*3 + 9*(a-4)] = -2;
//         SAT[8*3 + 9*(a-4)+1] = -3;
//         SAT[8*3 + 9*(a-4)+2] = -a;
//         SAT[9*3 + 9*(a-4)] = 2;
//         SAT[9*3 + 9*(a-4)+1] = 3;
//         SAT[9*3 + 9*(a-4)+2] = a;
//     }

//     // int SAT[34*3] = {-1,-2,-3,-1,2,3,1,-2,-3,1,-2,3,-1,2,-3,1,2,-3,-1,-2,3,-2,-3,4,2,3,4,-2,-3,-4,-2,-3,5,2,3,5,-2,-3,-5,-2,-3,6,2,3,6,-2,-3,-6,-2,-3,7,2,3,7,-2,-3,-7,-2,-3,8,2,3,8,-2,-3,-8,-2,-3,9,2,3,9,-2,-3,-9,-2,-3,10,2,3,10,-2,-3,-10,-2,-3,11,2,3,11,-2,-3,-11,-2,-3,12,2,3,12,-2,-3,-12};
//     for(int i = 1;i<=n*n;i++)
//         float b = SolveSAT(n,SAT,numClauses,i,v);

//     return 0;
// }

