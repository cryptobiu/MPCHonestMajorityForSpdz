
#include <stdlib.h>
#include "Protocol.h"
#include "ZpMersenneLongElement.h"
#include <smmintrin.h>
#include <inttypes.h>
#include <stdio.h>
#include <x86intrin.h>



/**
 * The main structure of our protocol is as follows:
 * 1. Initialization Phase: Initialize some global variables (parties, field, circuit, etc).
 * 2. Preparation Phase: Prepare enough random double-sharings: a random double-sharing is a pair of
 *  two sharings of the same random value, one with degree t, and one with degree 2t. One such double-
 *  sharing is consumed for multiplying two values. We also consume double-sharings for input gates
 *  and for random gates (this is slightly wasteful, but we assume that the number of multiplication
 *  gates is dominating the number of input and random gates).
 * 3. Input Phase: For each input gate, reconstruct one of the random sharings towards the input party.
 *  Then, all input parties broadcast a vector of correction values, namely the differences of the inputs
 *  they actually choose and the random values they got. These correction values are then added on
 *  the random sharings.
 * 4. Computation Phase: Walk through the circuit, and evaluate as many gates as possible in parallel.
 *  Addition gates and random gates can be evaluated locally (random gates consume a random double-
 *  sharing). Multiplication gates are more involved: First, every party computes local product of the
 *  respective shares; these shares form de facto a 2t-sharing of the product. Then, from this sharing,
 *  a degree-2t sharing of a random value is subtracted, the difference is reconstructed and added on
 *  the degree-t sharing of the same random value.
 * 5. Output Phase: The value of each output gate is reconstructed towards the corresponding party.
 * @param argc
 * @param argv[1] = id of parties (1,...,N)
 * @param argv[2] = N: number of parties
 * @param argv[3] = path of inputs file
 * @param argv[4] = path of output file
 * @param argv[5] = path of circuit file
 * @param argv[6] = address
 * @param argv[7] = fieldType
 * @return
 */


int main(int argc, char* argv[])
{

    if(argc != 4)
    {
        cout << "wrong number of arguments";
        return 0;
    }

    int times = 1;

    //string outputTimerFileName = string(argv[5]) + "Times" + string(argv[1]) + ".csv";
//    string outputTimerFileName = string(argv[5]) + "Times" + string(argv[1]) + argv[6] + argv[7] + argv[8] + argv[9] + ".csv";
//    ProtocolTimer p(times, outputTimerFileName);

//    string fieldType(argv[6]);

//
    TemplateField<ZpMersenneLongElement> *field = new TemplateField<ZpMersenneLongElement>(0);

    Protocol<ZpMersenneLongElement> protocol(3, atoi(argv[1]),atoi(argv[2]), atoi(argv[3]), 10,field);

    auto t1 = high_resolution_clock::now();
    for(int i=0; i<times; i++) {
        vector<ZpMersenneLongElement> shareArr;
        vector<ZpMersenneLongElement> shareArr1;
        shareArr.resize(10);
        vector<ZpMersenneLongElement> valueArr(10,2);
        bool flag = protocol.makeShare(0, valueArr, shareArr);

        cout<<"---------input flag is-----------" << flag<<endl;
        flag = protocol.offline();

        cout<<"---------offline flag is-----------" << flag<<endl;

        vector<ZpMersenneLongElement> triples(3);
        flag = protocol.triples(1, triples);

        vector<ZpMersenneLongElement> secrets(2);


        flag = protocol.openShare(2, shareArr, secrets);

        cout<<" x = " <<secrets[0]<< " y = "<< secrets[1] << " xy = " <<secrets[0] * secrets[1]<<endl;


        secrets.resize(3);
        flag = protocol.openShare(3, triples, secrets);

        auto shareA = triples[0];
        auto shareB = triples[1];
        auto shareC = triples[2];



        cout<<"a = " <<secrets[0]<< " b = "<< secrets[1] << " c = " <<secrets[2]<<endl;
        cout<<"a*b = "<< secrets[0] * secrets[1]<<endl;




        vector<ZpMersenneLongElement> sharesToOpen(2);
        sharesToOpen[0] = shareArr[0] - shareA; //[x-a]
        sharesToOpen[1] = shareArr[1] - shareB;//[y-b]



        flag = protocol.openShare(2, sharesToOpen, secrets);

        auto xMinusA = secrets[0];
        auto yMinusB = secrets[0];

        //cout<<"---------secrets----------- x-a = " <<xMinusA<< " y-b = "<< yMinusB <<endl;



        //[z] = (x − a) · [b]            + (y − b) · [a]           + (x − a) · (y − b)       + [c]  = [x · y]

        auto z = secrets[0] * triples[1] + secrets[1] * triples[0] + secrets[0] * secrets[1] + triples[2];

        sharesToOpen[0] = z;

        flag = protocol.openShare(1, sharesToOpen, secrets);

        cout<<"using beaver -- x*y = " <<secrets[0]<<endl;



        //generate one mult of inputs 3 and 4
        vector<ZpMersenneLongElement> shareOfXArr(1);
        vector<ZpMersenneLongElement> shareOfYArr(1);
        vector<ZpMersenneLongElement> shareOfXYArr(1);
        vector<ZpMersenneLongElement> secretOfXYArr(1);
        vector<ZpMersenneLongElement> inputSecrets(10);

        shareOfXArr[0] = shareArr[2];
        shareOfYArr[0] = shareArr[3];

        //mult the two shares
        flag = protocol.multShares(1, shareOfXArr, shareOfYArr,shareOfXYArr);

        cout<<"---------mult flag is-----------" << flag<<endl;

        flag = protocol.openShare(10, shareArr, inputSecrets);

        cout<<"x = " <<inputSecrets[2]<<"y = " <<inputSecrets[3]<<endl;

        flag = protocol.openShare(1, shareOfXYArr, secretOfXYArr);

        cout<<"internal mult -- x*y = " <<secretOfXYArr[0]<<endl;


        shareOfXArr[0] = shareArr[4];
        shareOfYArr[0] = shareArr[5];

        //mult the two shares
        flag = protocol.multShares(1, shareOfXArr, shareOfYArr,shareOfXYArr);
        flag = protocol.addShareAndScalar(shareOfXArr[0], shareOfYArr[0],shareOfXYArr[0]);



        cout<<"x = " <<shareOfXArr[0]<<"s = " << shareOfYArr[0] <<" = " <<shareOfXYArr[0]<< endl;




        cout<<"---------mult flag is-----------" << flag<<endl;


//
////
////        //spdz test for values
////
////        vector<ZpMersenneLongElement> SPDZtriples(3);
////        vector<ZpMersenneLongElement> SPDZshares(2);
////
////        if(atoi(argv[1])==0) {
////            //share a = 869738961550008488; share b = 983104380898725265; share c = 39913426787389368
////            SPDZtriples[0] =869738961550008488;
////            SPDZtriples[1] =983104380898725265;
////            SPDZtriples[2] =39913426787389368;
////
////            //input value 160727542214262586
////            SPDZshares[0] = 160727542214262586;
////            //input value 368007752470225898
////            SPDZshares[1] = 368007752470225898;
////
////        }
////        else if(atoi(argv[1])==1) {
////            //share a = 767820936118086415; share b = 1556751900158381769; share c = 2055043744731651588;
////            SPDZtriples[0] =767820936118086415;
////            SPDZtriples[1] =1556751900158381769;
////            SPDZtriples[2] =2055043744731651588;
////
////            //input value 321455084428525169
////            SPDZshares[0] = 321455084428525169;
////            //input value 736015504940451792
////            SPDZshares[1] = 736015504940451792;
////        }
////        else if(atoi(argv[1])==2) {
////            //share a = 665902910686164342; share b = 2130399419418038273; share c = 1764331053462219857;
////            SPDZtriples[0] =665902910686164342;
////            SPDZtriples[1] =2130399419418038273;
////            SPDZtriples[2] =1764331053462219857;
////
////            //input value 482182626642787752
////            SPDZshares[0] = 482182626642787752;
////            //input value 1104023257410677686
////            SPDZshares[1] = 1104023257410677686;
////
////        }
////
////
////        vector<ZpMersenneLongElement>sharesToOpen(2);
////        vector<ZpMersenneLongElement>secrets(2);
////
////        sharesToOpen[0] = SPDZshares[0] - SPDZtriples[0]; //[x-a]
////        sharesToOpen[1] = SPDZshares[1] - SPDZtriples[1];//[y-b]
////
////
////
////        flag = protocol.openShare(2, sharesToOpen, secrets);
////
////        auto xMinusA = secrets[0];
////        auto yMinusB = secrets[1];
////
////        cout<<"x-a  = " <<secrets[0]<<"y-b  = " <<secrets[1]<<endl;
////
////
////
////
////        //cout<<"---------secrets----------- x-a = " <<xMinusA<< " y-b = "<< yMinusB <<endl;
////
////
////
////        //[z] = (x − a) · [b]            + (y − b) · [a]           + (x − a) · (y − b)       + [c]  = [x · y]
////
////        auto SPDZz = secrets[0] * SPDZtriples[1] + secrets[1] * SPDZtriples[0] + secrets[0] * secrets[1] + SPDZtriples[2];
////
////        cout<<"spdz -- [z] = " <<SPDZz<<endl;
////
////
////        sharesToOpen[0] = SPDZz;
////
////        flag = protocol.openShare(1, sharesToOpen, secrets);
////
////        cout<<"spdz -- xy = " <<secrets[0]<<endl;
//    }
//







//    TemplateField<GF2E> *field = new TemplateField<GF2E>(40);
//
//    Protocol<GF2E> protocol(3, atoi(argv[1]),atoi(argv[2]), atoi(argv[3]), field);
//
//    auto t1 = high_resolution_clock::now();
//    for(int i=0; i<times; i++) {
//        vector<GF2E> shareArr;
//        vector<GF2E> shareArr1;
//        shareArr.resize(10);
//        bool flag = protocol.input(0, shareArr);
//        flag = protocol.input(1, shareArr1);
//
//        cout << "---------input flag is-----------" << flag << endl;
//        auto flag = protocol.offline();
//
//        cout << "---------offline flag is-----------" << flag << endl;
//
//        vector<GF2E> triples(3);
//        flag = protocol.triples(1, triples);
//
//        vector<GF2E> secrets(3);
//
//
//        flag = protocol.openShare(3, shareArr, secrets);
//
//        flag = protocol.multShares(2, shareArr, shareArr,shareArr);
//
//        cout<<" x = " <<secrets[0]<< " y = "<< secrets[1] << " xy = " <<secrets[0] * secrets[1]<< endl;
//
//
//        flag = protocol.openShare(1, shareArr, secrets);
//
//        cout<<" x*x = " <<secrets[0]<< "y*y = " <<secrets[1]<<endl;
//
//        secrets.resize(3);
//        flag = protocol.openShare(3, triples, secrets);
//
//        auto shareA = triples[0];
//        auto shareB = triples[1];
//        auto shareC = triples[2];
//
//
//
//        cout<<"a = " <<secrets[0]<< " b = "<< secrets[1] << " c = " <<secrets[2]<<endl;
//        cout<<"a*b = "<< secrets[0] * secrets[1]<<endl;
//
//
//
//
//        vector<ZpMersenneLongElement> sharesToOpen(2);
//        sharesToOpen[0] = shareArr[0] - shareA; //[x-a]
//        sharesToOpen[1] = shareArr[1] - shareB;//[y-b]
//
//
//
//        flag = protocol.openShare(2, sharesToOpen, secrets);
//
//        auto xMinusA = secrets[0];
//        auto yMinusB = secrets[0];
//
//        //cout<<"---------secrets----------- x-a = " <<xMinusA<< " y-b = "<< yMinusB <<endl;
//
//
//
//        //[z] = (x − a) · [b]            + (y − b) · [a]           + (x − a) · (y − b)       + [c]  = [x · y]
//
//        auto z = secrets[0] * triples[1] + secrets[1] * triples[0] + secrets[0] * secrets[1] + triples[2];
//
//        sharesToOpen[0] = z;
//
//        flag = protocol.openShare(1, sharesToOpen, secrets);
//
//        cout<<"using beaver -- x*y = " <<secrets[0]<<endl;
//
//
//
//        //generate one mult of inputs 3 and 4
//        vector<ZpMersenneLongElement> shareOfXArr(1);
//        vector<ZpMersenneLongElement> shareOfYArr(1);
//        vector<ZpMersenneLongElement> shareOfXYArr(1);
//        vector<ZpMersenneLongElement> secretOfXYArr(1);
//        vector<ZpMersenneLongElement> inputSecrets(10);
//
//
//

    }



    auto t2 = high_resolution_clock::now();

    auto duration = duration_cast<milliseconds>(t2-t1).count();
    cout << "time in milliseconds for " << times << " runs: " << duration << endl;

    delete field;

    //p.writeToFile();

    cout << "end main" << '\n';

    return 0;
}
