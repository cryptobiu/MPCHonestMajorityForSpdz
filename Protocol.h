#ifndef PROTOCOL_H_
#define PROTOCOL_H_

#include <stdlib.h>

#include "VDM.h"
#include "VDMTranspose.h"
#include "HIM.h"
#include "GRRHonestMult.h"
#include "DNHonestMult.h"
#include <vector>
#include <bitset>
#include <iostream>
#include <fstream>
#include <chrono>
#include "TemplateField.h"
#include "ProtocolTimer.h"
#include "MPCCommunication.h"
#include "../../include/infra/Common.hpp"
#include "../../include/primitives/Prg.hpp"
#include<emmintrin.h>
#include <thread>

#define flag_print false
#define flag_print_timings false
#define flag_print_output false



using namespace std;
using namespace std::chrono;

template <class FieldType>
class Protocol {


public:

    ///----------------------------API FOR SPDZ----------------------------------------------------------------------///


    /**
     * This constructor is eqvivalent to the setup procedure of SPDZ.
     * Values are initialized, the communication channels are set and the initialization phase of the protocol is called.
     *
     * @param n the number of parties.
     * @param id the party id.
     * @param numOfOpens the number of opens needed to execute the SPDZ program.
     * @param field the related field, currently restriced to mersenne 61.
     * @param multType the semi honest multiplication method. Possible values are "DN"/"GRR"
     */
    Protocol(int n, int id, int numOfOpens, int numOfMults, int numOfBits, TemplateField<FieldType> *field,
    		 string commFile = "Parties.txt", string multType = "DN");


    /**
     * If the partyID == myPartyId than, the input values are taken from a file. After the values are set the shares
     * for all the parties are created and sent to all the parties. The shareArr is filled locally accordingly.
     * Else if partyID != myPartyId than, we recieve the related shares from the party holding the inputs and
     * the shareArr is filled by those values.
     *
     * NOTE : Due to the SPDZ API this function works for a specific party only and not for the inputs of all the parties together.
     *
     * @param partyID the party that holds the inputs
     * @param shareArr an output array to fill with the shares of the inputs
     */
    bool makeShare(int partyID, vector<FieldType> &valuesArr, vector<FieldType> &shareArr);


    /**
     * This includes the offline part of the protocol.
     * The main steps of the offline protocols are as follows.
     * 1. Generating random shares.
     * 2. Generating 2*numOfOpens beaver triples using the semi-honest multiplication method.
     * 3. Verifying numOfOpens triples using the other numOfOpens triples created in the second step.
     */
    bool offline();


    /**
     * This procedure opens the input shares into their true value and put into the secret vector.
     *
     * @param numOfRandomShares the number of shares to open
     * @param Shares a vector holding all the shares to open
     * @param secrets an output vector to fill with the opened shares, that is, the true values of the related shamir shares
     * from all the parties
     */
    bool openShare(int numOfRandomShares, vector<FieldType> &shares, vector<FieldType> &secrets);



    /**
     * Curently this function does two things:
     * 1. mults every item in the x array by every corresponding item in the y array and returns the result
     *    as the xy shares
     * 2. verifies that the resulting multiplication to detect cheats.
     *
     * NOTE: The second part should be removed later on. The verification should be done for all the multiplication
     *       at once.
     *
     * @param numOfpairs the number of x,y pairs to mult
     * @param xShares an array of all the x's to mult
     * @param yShares an array of all the y's to mult
     * @param outputXYShares the output array containing the shares with the resulting xy shares
     * @return
     */
    bool multShares(int numOfpairs, vector<FieldType> &xShares, vector<FieldType> &yShares, vector<FieldType> &outputXYShares);


    static bool addShareAndScalar(FieldType &xShare, FieldType &scalar, FieldType &outputShare){ outputShare = xShare + scalar; return true;};
    static bool shareSubScalar(FieldType &xShare, FieldType &scalar, FieldType &outputShare){ outputShare = xShare - scalar; return true;};
    static bool scalarSubShare(FieldType &xShare, FieldType &scalar, FieldType &outputShare){ outputShare = scalar - xShare; return true;};



    /**
     * This function currently does nothing, since the verification is embedded in the open procedure.
     */
    bool verify() { return true; }

    bool triples(int numOfTriples, vector<FieldType> &triples);

    bool bits(int numOfBits, vector<FieldType> &bits);

    /**
     * Load the reference shares with input values from the reference integers.
     *
     * @param partyID the party that simulates the input provision
     * @param shareArr an output array to fill with the shares of the inputs
     * @param valueArr the immediate values to be used as inputs
     */
    //bool load_share_immediates(int partyID, vector<FieldType> &shareArr, const vector<std::string> &valueStrArr);


    /**
     * Cleaning up. Currently has nothing to clean, but should still be called by the SPDZ extension.
     */
    ~Protocol();




private:
    /**
     * N - number of parties
     * T - number of malicious
     */
    int N, T, m_partyId;
    VDM<FieldType> matrix_vand;
    vector<FieldType> firstRowVandInverse;
    TemplateField<FieldType> *field;
    vector<shared_ptr<ProtocolPartyData>>  parties;
    vector<FieldType> randomTAnd2TShares;
    int numOfOpens;
    int numOfMults;
    int numOfBits;

    FieldType inv_2;
    vector<FieldType> vecBitShares;


    string genRandomSharesType, multType/*, verifyType*/;
    ProtocolTimer* protocolTimer;

    string s;
    string outputFile;
    vector<FieldType> beta;
    HIM<FieldType> matrix_for_interpolate;
    HIM<FieldType> matrix_for_t;
    HIM<FieldType> matrix_for_2t;


    HIM<FieldType> matrix_him;
    VDMTranspose<FieldType> matrix_vand_transpose;
    HIM<FieldType> m;



    vector<FieldType> randomABShares;//a, b random shares
    vector<FieldType> c;//a vector of a*b shares

    vector<FieldType> randomABSharesForMult;//a, b random shares
    vector<FieldType> cForMult;//a vector of a*b shares
    int multIndex;


    boost::asio::io_service io_service;
    vector<FieldType> alpha; // N distinct non-zero field elements

    vector<FieldType> sharingBufTElements; // prepared T-sharings (my shares)
    vector<FieldType> sharingBuf2TElements; // prepared 2T-sharings (my shares)
    int trippleIndex;
    int bitsIndex;


    HonestMultAbstract<FieldType> *honestMult = nullptr;//


private :


    void roundFunctionSync(vector<vector<byte>> &sendBufs, vector<vector<byte>> &recBufs, int round);
    void exchangeData(vector<vector<byte>> &sendBufs,vector<vector<byte>> &recBufs, int first, int last);
    void roundFunctionSyncBroadcast(vector<byte> &message, vector<vector<byte>> &recBufs);
    void recData(vector<byte> &message, vector<vector<byte>> &recBufs, int first, int last);
    void roundFunctionSyncForP1(vector<byte> &myShare, vector<vector<byte>> &recBufs);
    void recDataToP1(vector<vector<byte>> &recBufs, int first, int last);

    void sendDataFromP1(vector<byte> &sendBuf, int first, int last);
    void sendFromP1(vector<byte> &sendBuf);



    vector<int> firstIndex;
    int counter = 0;

    /**
     * We describe the protocol initialization.
     * In particular, some global variables are declared and initialized.
     */
    void initializationPhase();

    /**
     * Generates the first row of the inverse of the vandermonde matrix and saves it
     * in firstRowVandInverse for future use.
     */
    void initFirstRowInvVDM();

    /**
     * A random double-sharing is a pair of two sharings of the same random value, where the one sharing is
     * of degree t, and the other sharing is of degree 2t. Such random double-sharing are of big help in the
     * multiplication protocol.
     * We use hyper-invertible matrices to generate random double-sharings. The basic idea is as follows:
     * Every party generates one random double-sharing. These n double-sharings are processes through a
     * hyper-invertible matrix. From the resulting n double-sharings, t are checked to be valid (correct degree,
     * same secret), and t are then kept as “good” double-sharings. This is secure due to the diversion property
     * of hyper-invertible matrix: We know that n − t of the input double-sharings are good. So, if there are t
     * valid output double-sharings, then all double-sharings must be valid. Furthermore, the adversary knows
     * his own up to t input double-sharings, and learns t output double sharings. So, n − 2t output double
     * sharings are random and unknown to the adversary.
     * For the sake of efficiency, we do not publicly reconstruct t of the output double-sharings. Rather, we
     * reconstruct 2t output double sharings, each to one dedicated party only. At least t of these parties are
     * honest and correctly validate the reconstructed double-sharing.
     *
     * The goal of this phase is to generate “enough” double-sharings to evaluate the circuit. The double-
     * sharings are stored in a buffer SharingBuf , where alternating a degree-t and a degree-2t sharing (of the same secret)
     * is stored (more precisely, a share of each such corresponding sharings is stored).
     * The creation of double-sharings is:
     *
     * Protocol Generate-Double-Sharings:
     * 1. ∀i: Pi selects random value x-(i) and computes degree-t shares x1-(i) and degree-2t shares x2-(i).
     * 2. ∀i,j: Pi sends the shares x1,j and X2,j to party Pj.
     * 3. ∀j: Pj applies a hyper-invertible matrix M on the received shares, i.e:
     *      (y1,j,..., y1,j) = M(x1,j,...,x1,j)
     *      (y2,j,...,y2,j) = M (x2,j,...,x2,)
     * 4. ∀j, ∀k ≤ 2t: Pj sends y1,j and y2,j to Pk.
     * 5. ∀k ≤ 2t: Pk checks:
     *      • that the received shares (y1,1,...,y1,n) are t-consistent,
     *      • that the received shares (y2,1,...,y2,n) are 2t-consistent, and
     *      • that both sharings interpolate to the same secret.
     *
     * We use this algorithm, but extend it to capture an arbitrary number of double-sharings.
     * This is, as usual, achieved by processing multiple buckets in parallel.
     */
    //bool offline();

    /**
     * The aim of this protocol is to produce random triples that will be later used in the main protocol.
     * The protocol takes as input an integer numOfTriples (the number of triples to output), and proceeds as follows.
     * 1. The parties call a protocol to generate random shares to obtain 2L random
     * 2. The parties run a semi-honest multiplication sub protocol to obtain a sharing of a*b
     */
    void generateBeaverTriples(int numOfTriples, vector<FieldType> &aBShares, vector<FieldType> &cToFill);

    bool generateRandomBits();

    /**
     * This protocol is secure only in the presence of a semi-honest adversary.
     */
    void GRRHonestMultiplication(FieldType *a, FieldType *b, vector<FieldType> &cToFill, int numOfTrupples);

    void DNHonestMultiplication(FieldType *a, FieldType *b, vector<FieldType> &cToFill, int numOfTrupples);

    void offlineDNForMultiplication(int numOfTriples);

    /**
     * The input phase proceeds in two steps:
     * First, for each input gate, the party owning the input creates shares for that input by choosing a random coefficients for the polynomial
     * Then, all the shares are sent to the relevant party
     */
//    void inputPhase();
//    void inputVerification();


    bool inputVer(vector<FieldType> &shareArr);

    void generateRandomShares(int numOfRnadoms, vector<FieldType>& randomElementsToFill);
    void generateRandom2TAndTShares(int numOfRandomPairs, vector<FieldType>& randomElementsToFill);


    /**
     * Check whether given points lie on polynomial of degree d.
     * This check is performed by interpolating x on the first d + 1 positions of α and check the remaining positions.
     */
    bool checkConsistency(vector<FieldType>& x, int d);

    FieldType reconstructShare(vector<FieldType>& x, int d, bool &flag);

    /**
     * The cheap way: Create a HIM from the αi’s onto ZERO (this is actually a row vector), and multiply
     * this HIM with the given x-vector (this is actually a scalar product).
     * The first (and only) element of the output vector is the secret.
     */
    FieldType interpolate(vector<FieldType> x);

    FieldType tinterpolate(vector<FieldType> x);


    /**
     * Walk through the circuit and verify the multiplication gates.
     * We first generate the random elements using a common AES key that was generated by the parties,
     * run the relevant verification algorithm and return accept/reject according to the output
     * of the verification algorithm.
     */
//    void verificationPhase();

    bool verifyTriples();
    bool verifyMults(int numOfpairs, vector<FieldType> &xShares, vector<FieldType> &yShares, vector<FieldType> &outputXYShares);


    vector<byte> generateCommonKey();
    void generatePseudoRandomElements(vector<byte> & aesKey, vector<FieldType> &randomElementsToFill, int numOfRandomElements);

//    bool verificationOfBatchedTriples(FieldType *x, FieldType *y, FieldType *z,
//                                      FieldType *a, FieldType *b, FieldType *c,
//                                      FieldType *randomElements, int numOfTriples);

    bool verificationOfSingleTriples(FieldType *x, FieldType *y, FieldType *z,
                                      FieldType *a, FieldType *b, FieldType *c,
                                     FieldType *randomElements, int numOfTriples);


    friend class DNHonestMult<FieldType>;
    friend class GRRHonestMult<FieldType>;

};


template <class FieldType>
Protocol<FieldType>::Protocol(int n, int id, int numOfOpens, int numOfMults, int numOfBits ,
                              TemplateField<FieldType> *field, string commFile, string multType)
{


    this->multType = multType;
    //this->verifyType = verifyType;
    this->numOfOpens = numOfOpens;
    this->numOfMults = numOfMults;
    this->numOfBits = numOfBits;

    if(multType=="GRR"){
        honestMult = new GRRHonestMult<FieldType>(this);
    }
    if(multType=="DN"){
        //if(verifyType=="Single")
            honestMult = new DNHonestMult<FieldType>(0, this);
//        else if(verifyType=="Batch")
//            honestMult = new DNHonestMult<FieldType>(0, this);
    }

    this->protocolTimer = protocolTimer;
    this->field = field;
    //comm = Communication::getInstance(n, id, address);
    N = n;
    T = n/2 - 1;
    this->outputFile = outputFile;
    if(n%2 > 0)
    {
        T++;
    }

    inv_2 = (field->GetElement(1))/(field->GetElement(2)); // calculate the inverse of 3 in the field


    m_partyId = id;
    s = to_string(m_partyId);
    trippleIndex = 0;
    bitsIndex = 0;
    multIndex = 0;
    counter = 0;

    MPCCommunication comm;

    parties = comm.setCommunication(io_service, m_partyId, N, commFile.c_str(), 0);

    string tmp = "init times";
    //cout<<"before sending any data"<<endl;
    byte tmpBytes[20];
    for (int i=0; i<parties.size(); i++){
        if (parties[i]->getID() < id){
            parties[i]->getChannel()->write(tmp);
            parties[i]->getChannel()->read(tmpBytes, tmp.size());
        } else {
            parties[i]->getChannel()->read(tmpBytes, tmp.size());
            parties[i]->getChannel()->write(tmp);
        }
    }



    auto t1 = high_resolution_clock::now();
    initializationPhase(/*matrix_him, matrix_vand, m*/);

    auto t2 = high_resolution_clock::now();

    auto duration = duration_cast<milliseconds>(t2-t1).count();
    if(flag_print_timings) {
        cout << "time in milliseconds initializationPhase: " << duration << endl;
    }
}

/**
 * the function implements the second step of Input Phase:
 * the party broadcasts for each input gate the difference between
 * the random secret and the actual input value.
 * @param diff
 */
template <class FieldType>
bool Protocol<FieldType>::makeShare(int partyID, vector<FieldType> &valueArr, vector<FieldType> &shareArr)
{
    // the number of random double sharings we need altogether
    vector<FieldType> x1(N),y1(N);
    vector<vector<FieldType>> sendBufsElements(N);
    vector<vector<byte>> sendBufsBytes(N);
    vector<vector<byte>> recBufBytes(N);
    vector<vector<FieldType>> recBufElements(N);

    int index = 0;
    vector<int> sizes(N);


    //cout<< "before input: "<<endl;
    //cout<< "shareArray size is: "<<shareArr.size()<<endl;


    // prepare the shares for the input
    for (int k = 0; k < shareArr.size(); k++)
    {
        sizes[partyID]++;

        if (partyID == m_partyId) {

            index++;
            // the value of a_0 is the input of the party.
            x1[0] = valueArr[k];


            // generate random degree-T polynomial
            for(int i = 1; i < T+1; i++)
            {
                // A random field element, uniform distribution
                x1[i] = field->Random();

            }


            matrix_vand.MatrixMult(x1, y1, T+1); // eval poly at alpha-positions predefined to be alpha_i = i

            // prepare shares to be sent
            for(int i=0; i < N; i++)
            {
                //cout << "y1[ " <<i<< "]" <<y1[i] << endl;
                sendBufsElements[i].push_back(y1[i]);

            }
        }

    }

    int fieldByteSize = field->getElementSizeInBytes();
    for(int i=0; i < N; i++)
    {
        sendBufsBytes[i].resize(sendBufsElements[i].size()*fieldByteSize);
        //cout<< "size of sendBufs1Elements["<<i<<" ].size() is " << sendBufs1Elements[i].size() <<"myID =" <<  m_partyId<<endl;
        recBufBytes[i].resize(sizes[i]*fieldByteSize);
        for(int j=0; j<sendBufsElements[i].size();j++) {
            field->elementToBytes(sendBufsBytes[i].data() + (j * fieldByteSize), sendBufsElements[i][j]);
        }
    }


    roundFunctionSync(sendBufsBytes, recBufBytes,10);


    //turn the bytes to elements
    for(int i=0; i < N; i++)
    {
        recBufElements[i].resize((recBufBytes[i].size()) / fieldByteSize);
        for(int j=0; j<recBufElements[i].size();j++) {
            recBufElements[i][j] = field->bytesToElement(recBufBytes[i].data() + ( j * fieldByteSize));
        }
    }

    for (int k = 0; k < shareArr.size(); k++)
    {

            auto share = recBufElements[partyID][k];

            shareArr[k] = share; // set the share sent from the party owning the input


    }

    return inputVer(shareArr);
}

/*
template <class FieldType>
bool Protocol<FieldType>::load_share_immediates(int partyID, vector<FieldType> &shareArr, const vector<std::string> &valueStrArr)
{
    // the number of random double sharings we need altogether
    vector<FieldType> x1(N),y1(N);
    vector<vector<FieldType> > sendBufsElements(N);
    vector<vector<byte> > sendBufsBytes(N);
    vector<vector<byte> > recBufBytes(N);
    vector<vector<FieldType> > recBufElements(N);

    //int index = 0;
    vector<int> sizes(N);

    // prepare the shares for the input
    for (int k = 0; k < shareArr.size(); k++)
    {
        sizes[partyID]++;

        if (partyID == m_partyId) {

            x1[0] = field->stringToElement(valueStrArr[k]);


            // generate random degree-T polynomial
            for(int i = 1; i < T+1; i++)
            {
                // A random field element, uniform distribution
                x1[i] = field->Random();

            }


            matrix_vand.MatrixMult(x1, y1, T+1); // eval poly at alpha-positions predefined to be alpha_i = i

            // prepare shares to be sent
            for(int i=0; i < N; i++)
            {
                //cout << "y1[ " <<i<< "]" <<y1[i] << endl;
                sendBufsElements[i].push_back(y1[i]);

            }
        }

    }

    int fieldByteSize = field->getElementSizeInBytes();
    for(int i=0; i < N; i++)
    {
        sendBufsBytes[i].resize(sendBufsElements[i].size()*fieldByteSize);
        //cout<< "size of sendBufs1Elements["<<i<<" ].size() is " << sendBufs1Elements[i].size() <<"myID =" <<  m_partyId<<endl;
        recBufBytes[i].resize(sizes[i]*fieldByteSize);
        for(int j=0; j<sendBufsElements[i].size();j++) {
            field->elementToBytes(sendBufsBytes[i].data() + (j * fieldByteSize), sendBufsElements[i][j]);
        }
    }


    roundFunctionSync(sendBufsBytes, recBufBytes,10);


    //turn the bytes to elements
    for(int i=0; i < N; i++)
    {
        recBufElements[i].resize((recBufBytes[i].size()) / fieldByteSize);
        for(int j=0; j<recBufElements[i].size();j++) {
            recBufElements[i][j] = field->bytesToElement(recBufBytes[i].data() + ( j * fieldByteSize));
        }
    }

    for (int k = 0; k < shareArr.size(); k++)
    {

            auto share = recBufElements[partyID][k];

            shareArr[k] = share; // set the share sent from the party owning the input


    }

    return inputVer(shareArr);
}*/


template <class FieldType>
bool Protocol<FieldType>::inputVer(vector<FieldType> &shareArr){

    //first generate the common aes key
    auto key = generateCommonKey();

    //print key
    if (flag_print) {
        for (int i = 0; i < key.size(); i++) {
            cout << "key[" << i << "] for party :" << m_partyId << "is : " << (int) key[i] << endl;
        }
    }

    //calc the number of times we need to run the verification -- ceiling
    //int iterations =   (5 + field->getElementSizeInBytes() - 1) / field->getElementSizeInBytes();

    vector<FieldType> randomElements(shareArr.size());
    generatePseudoRandomElements(key, randomElements, shareArr.size());


    vector<FieldType> r(1);//vector holding the random shares generated
    vector<FieldType> v(1);
    vector<FieldType> secret(1);

    generateRandomShares(1, r);


    for (int i = 0; i < shareArr.size(); i++)
        v[0] += randomElements[i] * shareArr[i];

    v[0] += r[0];


    //if all the the parties share lie on the same polynomial this will not fail
    return openShare(1, v, secret);

}


template <class FieldType>
void Protocol<FieldType>::generateRandomShares(int numOfRandoms, vector<FieldType>& randomElementsToFill){


    int index = 0;
    vector<vector<byte>> recBufsBytes(N);
    int robin = 0;
    int no_random = numOfRandoms;

    vector<FieldType> x1(N),y1(N), x2(N),y2(N), t1(N), r1(N), t2(N), r2(N);;

    vector<vector<FieldType>> sendBufsElements(N);
    vector<vector<byte>> sendBufsBytes(N);

    // the number of buckets (each bucket requires one double-sharing
    // from each party and gives N-2T random double-sharings)
    int no_buckets = (no_random / (N-T))+1;

    //maybe add some elements if a partial bucket is needed
    randomElementsToFill.resize(no_buckets*(N-T));


    for(int i=0; i < N; i++)
    {
        sendBufsElements[i].resize(no_buckets);
        sendBufsBytes[i].resize(no_buckets*field->getElementSizeInBytes());
        recBufsBytes[i].resize(no_buckets*field->getElementSizeInBytes());
    }

    /**
     *  generate random sharings.
     *  first degree t.
     *
     */
    for(int k=0; k < no_buckets; k++)
    {
        // generate random degree-T polynomial
        for(int i = 0; i < T+1; i++)
        {
            // A random field element, uniform distribution, note that x1[0] is the secret which is also random
            x1[i] = field->Random();

        }

        matrix_vand.MatrixMult(x1, y1,T+1); // eval poly at alpha-positions

        // prepare shares to be sent
        for(int i=0; i < N; i++)
        {
            //cout << "y1[ " <<i<< "]" <<y1[i] << endl;
            sendBufsElements[i][k] = y1[i];

        }
    }

    if(flag_print) {
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < sendBufsElements[0].size(); k++) {

                // cout << "before roundfunction4 send to " <<i <<" element: "<< k << " " << sendBufsElements[i][k] << endl;
            }
        }
        cout << "sendBufs" << endl;
        cout << "N" << N << endl;
        cout << "T" << T << endl;
    }

    int fieldByteSize = field->getElementSizeInBytes();
    for(int i=0; i < N; i++)
    {
        for(int j=0; j<sendBufsElements[i].size();j++) {
            field->elementToBytes(sendBufsBytes[i].data() + (j * fieldByteSize), sendBufsElements[i][j]);
        }
    }

    roundFunctionSync(sendBufsBytes, recBufsBytes,4);


    if(flag_print) {
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < sendBufsBytes[0].size(); k++) {

                cout << "roundfunction4 send to " <<i <<" element: "<< k << " " << (int)sendBufsBytes[i][k] << endl;
            }
        }
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < recBufsBytes[0].size(); k++) {
                cout << "roundfunction4 receive from " <<i <<" element: "<< k << " " << (int) recBufsBytes[i][k] << endl;
            }
        }
    }

    for(int k=0; k < no_buckets; k++) {
        for (int i = 0; i < N; i++) {
            t1[i] = field->bytesToElement(recBufsBytes[i].data() + (k * fieldByteSize));

        }
        matrix_vand_transpose.MatrixMult(t1, r1,N-T);

        //copy the resulting vector to the array of randoms
        for (int i = 0; i < N - T; i++) {

            randomElementsToFill[index] = r1[i];
            index++;

        }
    }

}


template <class FieldType>
void Protocol<FieldType>::generateRandom2TAndTShares(int numOfRandomPairs, vector<FieldType>& randomElementsToFill){


    int index = 0;
    vector<vector<byte>> recBufsBytes(N);
    int robin = 0;
    int no_random = numOfRandomPairs;

    vector<FieldType> x1(N),y1(N), x2(N),y2(N), t1(N), r1(N), t2(N), r2(N);;

    vector<vector<FieldType>> sendBufsElements(N);
    vector<vector<byte>> sendBufsBytes(N);

    // the number of buckets (each bucket requires one double-sharing
    // from each party and gives N-2T random double-sharings)
    int no_buckets = (no_random / (N-T))+1;

     //maybe add some elements if a partial bucket is needed
    randomElementsToFill.resize(no_buckets*(N-T)*2);


    for(int i=0; i < N; i++)
    {
        sendBufsElements[i].resize(no_buckets*2);
        sendBufsBytes[i].resize(no_buckets*field->getElementSizeInBytes()*2);
        recBufsBytes[i].resize(no_buckets*field->getElementSizeInBytes()*2);
    }

    /**
     *  generate random sharings.
     *  first degree t.
     *
     */
    for(int k=0; k < no_buckets; k++)
    {
        // generate random degree-T polynomial
        for(int i = 0; i < T+1; i++)
        {
            // A random field element, uniform distribution, note that x1[0] is the secret which is also random
            x1[i] = field->Random();

        }

        matrix_vand.MatrixMult(x1, y1,T+1); // eval poly at alpha-positions

        x2[0] = x1[0];
        // generate random degree-T polynomial
        for(int i = 1; i < 2*T+1; i++)
        {
            // A random field element, uniform distribution, note that x1[0] is the secret which is also random
            x2[i] = field->Random();

        }

        matrix_vand.MatrixMult(x2, y2,2*T+1);

        // prepare shares to be sent
        for(int i=0; i < N; i++)
        {
            //cout << "y1[ " <<i<< "]" <<y1[i] << endl;
            sendBufsElements[i][2*k] = y1[i];
            sendBufsElements[i][2*k + 1] = y2[i];

        }
    }

    if(flag_print) {
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < sendBufsElements[0].size(); k++) {

                // cout << "before roundfunction4 send to " <<i <<" element: "<< k << " " << sendBufsElements[i][k] << endl;
            }
        }
        cout << "sendBufs" << endl;
        cout << "N" << N << endl;
        cout << "T" << T << endl;
    }

    int fieldByteSize = field->getElementSizeInBytes();
    for(int i=0; i < N; i++)
    {
        for(int j=0; j<sendBufsElements[i].size();j++) {
            field->elementToBytes(sendBufsBytes[i].data() + (j * fieldByteSize), sendBufsElements[i][j]);
        }
    }

    roundFunctionSync(sendBufsBytes, recBufsBytes,4);


    if(flag_print) {
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < sendBufsBytes[0].size(); k++) {

                cout << "roundfunction4 send to " <<i <<" element: "<< k << " " << (int)sendBufsBytes[i][k] << endl;
            }
        }
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < recBufsBytes[0].size(); k++) {
                cout << "roundfunction4 receive from " <<i <<" element: "<< k << " " << (int) recBufsBytes[i][k] << endl;
            }
        }
    }

    for(int k=0; k < no_buckets; k++) {
        for (int i = 0; i < N; i++) {
            t1[i] = field->bytesToElement(recBufsBytes[i].data() + (2*k * fieldByteSize));
            t2[i] = field->bytesToElement(recBufsBytes[i].data() + ((2*k +1) * fieldByteSize));

        }
        matrix_vand_transpose.MatrixMult(t1, r1,N-T);
        matrix_vand_transpose.MatrixMult(t2, r2,N-T);

        //copy the resulting vector to the array of randoms
        for (int i = 0; i < (N - T); i++) {

            randomElementsToFill[index*2] = r1[i];
            randomElementsToFill[index*2 +1] = r2[i];
            index++;

        }
    }

}

/**
 * some global variables are initialized
 * @param GateValueArr
 * @param GateShareArr
 * @param matrix_him
 * @param matrix_vand
 * @param alpha
 */
template <class FieldType>
void Protocol<FieldType>::initializationPhase()
{
    beta.resize(1);
//    gateValueArr.resize(M);  // the value of the gate (for my input and output gates)
    //gateShareArr.resize(M - circuit.getNrOfOutputGates()); // my share of the gate (for all gates)
    alpha.resize(N); // N distinct non-zero field elements
    vector<FieldType> alpha1(N-T);
    vector<FieldType> alpha2(T);

    beta[0] = *field->GetZero(); // zero of the field
    matrix_for_interpolate.allocate(1,N, field);


    matrix_him.allocate(N,N,field);
    matrix_vand.allocate(N,N,field);
    matrix_vand_transpose.allocate(N,N,field);
    m.allocate(T, N-T,field);

    // Compute Vandermonde matrix VDM[i,k] = alpha[i]^k
    matrix_vand.InitVDM();
    matrix_vand_transpose.InitVDMTranspose();

    // Prepare an N-by-N hyper-invertible matrix
    matrix_him.InitHIM();

    // N distinct non-zero field elements
    for(int i=0; i<N; i++)
    {
        alpha[i]=field->GetElement(i+1);
    }

    for(int i = 0; i < N-T; i++)
    {
        alpha1[i] = alpha[i];
    }
    for(int i = N-T; i < N; i++)
    {
        alpha2[i - (N-T)] = alpha[i];
    }

    m.InitHIMByVectors(alpha1, alpha2);

    matrix_for_interpolate.InitHIMByVectors(alpha, beta);

    vector<FieldType> alpha_until_t(T + 1);
    vector<FieldType> alpha_from_t(N - 1 - T);

    // Interpolate first d+1 positions of (alpha,x)
    matrix_for_t.allocate(N - 1 - T, T + 1, field); // slices, only positions from 0..d
    //matrix_for_t.InitHIMByVectors(alpha_until_t, alpha_from_t);
    matrix_for_t.InitHIMVectorAndsizes(alpha, T+1, N-T-1);

    vector<FieldType> alpha_until_2t(2*T + 1);
    vector<FieldType> alpha_from_2t(N - 1 - 2*T);

    // Interpolate first d+1 positions of (alpha,x)
    matrix_for_2t.allocate(N - 1 - 2*T, 2*T + 1, field); // slices, only positions from 0..d
    //matrix_for_2t.InitHIMByVectors(alpha_until_2t, alpha_from_2t);
    matrix_for_2t.InitHIMVectorAndsizes(alpha, 2*T + 1, N-(2*T +1));


    //create the first row of the inverse of the nxn vandemonde matrix firstRowVandInverse
    initFirstRowInvVDM();



    if(flag_print){
        cout<< "matrix_for_t : " <<endl;
        matrix_for_t.Print();

        cout<< "matrix_for_2t : " <<endl;
        matrix_for_2t.Print();

    }


}
template <class FieldType>
void Protocol<FieldType>::initFirstRowInvVDM(){
    firstRowVandInverse.resize(N);

    //first calc the multiplication of all the alpha's for the denominator and the diff for the
    FieldType accumMult = *field->GetOne();

    for(int i=0; i<N; i++){
        accumMult*= alpha[i];

    }

    FieldType accum = *field->GetOne();
    for(int j=0; j<N; j++){

        for(int m=0; m<N; m++){
            if(m!=j) {
                accum *= (alpha[m] - alpha[j]);
            }
        }

        firstRowVandInverse[j] = accumMult/alpha[j]/accum;
        accum = *field->GetOne();
    }
}



template <class FieldType>
bool Protocol<FieldType>::offline()
{

    trippleIndex = 0;
    bitsIndex = 0;

    //generate triples for the DN multiplication protocol
    if(multType=="DN") {

        //if(verifyType=="Single")
            offlineDNForMultiplication(2 * numOfOpens + numOfMults + N*numOfBits);
//        else if(verifyType=="Batch") {
//
//            int iterations =   (5 + field->getElementSizeInBytes() - 1) / field->getElementSizeInBytes();
//            offlineDNForMultiplication(
//                    2 * circuit.getNrOfMultiplicationGates() + 4 * iterations* circuit.getNrOfMultiplicationGates());
//        }
    }


    //cout<<"before offline"<<endl;
    honestMult->invokeOffline();//need to pust some

    //cout<<"after offline"<<endl;


    //generate verified triples for the "open" instruction of SPDZ
    generateBeaverTriples(numOfOpens*2, randomABShares,c);

    //cout<<"after generating triples"<<endl;

    bool flag = verifyTriples();

    //generate un-verified triples for the multiplication instruction
    generateBeaverTriples(numOfMults, randomABSharesForMult,cForMult);

    flag = generateRandomBits();

    return flag;



}

template <class FieldType>
void Protocol<FieldType>::generateBeaverTriples(int numOfTriples, vector<FieldType> &aBShares, vector<FieldType> &cToFill){

    //int iterations =   (5 + field->getElementSizeInBytes() - 1) / field->getElementSizeInBytes();

    aBShares.resize(numOfTriples*2);//a, b random shares
    cToFill.resize(numOfTriples);//a vector of a*b shares


    auto t1 = high_resolution_clock::now();

    //first generate 2*numOfTriples random shares
    generateRandomShares(numOfTriples*2 ,aBShares);


    auto t2 = high_resolution_clock::now();

    auto duration = duration_cast<milliseconds>(t2-t1).count();
    if(flag_print_timings) {
        cout << "time in milliseconds generateRandomShares: " << duration << endl;
    }

    t1 = high_resolution_clock::now();
    honestMult->mult(aBShares.data(), aBShares.data()+numOfTriples, cToFill, numOfTriples);

    t2 = high_resolution_clock::now();

    duration = duration_cast<milliseconds>(t2-t1).count();
    if(flag_print_timings) {
        cout << "time in milliseconds honestMult: " << duration << endl;
    }
}



template <class FieldType>
bool Protocol<FieldType>::generateRandomBits(){

    //int iterations =   (5 + field->getElementSizeInBytes() - 1) / field->getElementSizeInBytes();

    vector<vector<FieldType>> vecAllShares(N);//each subvector holds the share sent from the corresponding party
    vector<FieldType> vecValues(numOfBits);//the 1/p-1 values
    vector<FieldType> vecSharesSqaure(numOfBits);
    vector<FieldType> vecValuesSquare(numOfBits);
    vector<FieldType> vecValuesDummy(numOfBits);//a dummy vector needed by the interface but has no


    //generate a random prg and extract bits from it
    PrgFromOpenSSLAES prg((numOfBits/16) + 1);
    auto key =prg.generateKey(128);
    prg.setKey(key);
    vector<byte> randomBytes(numOfBits);
    prg.getPRGBytes(randomBytes, 0, numOfBits);

    //set random -1/1 values
    for(int i=0; i<numOfBits; i++)
    {
        if(randomBytes[i]& 1)
            vecValues[i] = *field->GetOne();
        else
            vecValues[i] = *field->GetZero() - *field->GetOne();
    }


    for(int i=0; i<N ;i++){
        vecAllShares[i].resize(numOfBits);
    }

    //generate lamda_i
    for(int i=0; i<N; i++){

        if(m_partyId==i)
            makeShare(i, vecValues, vecAllShares[i]);
        else
            makeShare(i, vecValuesDummy, vecAllShares[i]);
    }
    auto t2 = high_resolution_clock::now();


    vecBitShares = vecAllShares[0];
    //mult all the shares for each party with the result of the previous round
    //This should be optimized to have only log depth for the case of many parties
    for(int i=1; i<N;i++){

        honestMult->mult(vecBitShares.data(), vecAllShares[i].data(), vecBitShares, numOfBits);
    }

    //square
    honestMult->mult(vecBitShares.data(), vecBitShares.data(), vecSharesSqaure, numOfBits);

    //open the squared shares for checking

    openShare(numOfBits,vecSharesSqaure, vecValuesSquare);


    //now normalize the result to be 0/1

    for(int i=0; i<numOfBits;i++)
    {
        //(lamda_i+1)/2. For performance mult by the inverse of 2
        vecBitShares[i] = (vecBitShares[i]+ *field->GetOne())*inv_2;
    }


    //verify that no one cheated
    for(int i=0; i<numOfBits; i++){

        if(vecValuesSquare[i]!=*field->GetOne())
            return false;
    }

    return true;

}


template <class FieldType>
bool Protocol<FieldType>::triples(int numOfTriples, vector<FieldType> &triples){

    if(triples.size()/3!=numOfTriples) {

        cerr<<"the size of triples array should be a multiplication of 3 and numOfTriples" <<endl;
        return false;
    }

    for(int i=0; i<numOfTriples; i++){

        triples[3*i] = randomABShares[trippleIndex];//set a
        triples[3*i+1] = randomABShares[trippleIndex +2*numOfOpens];//set b
        triples[3*i+2] = c[trippleIndex];//set c

        trippleIndex++;

    }

    return true;
}
template <class FieldType>
bool Protocol<FieldType>::bits(int numOfBits, vector<FieldType> &bits){

    if(bits.size()!=numOfBits) {

        cerr<<"the size of bits array is not numOfBits" <<endl;
        return false;
    }


    for(int i=0; i<numOfBits; i++){

        bits[i] = vecBitShares[bitsIndex];

        bitsIndex++;

    }

    return true;
}


/**
 * Check whether given points lie on polynomial of degree d. This check is performed by interpolating x on
 * the first d + 1 positions of α and check the remaining positions.
 */
template <class FieldType>
bool Protocol<FieldType>::checkConsistency(vector<FieldType>& x, int d)
{
    if(d == T)
    {
        vector<FieldType> y(N - 1 - d); // the result of multiplication
        vector<FieldType> x_until_t(T + 1);

        for (int i = 0; i < T + 1; i++) {
            x_until_t[i] = x[i];
        }


        matrix_for_t.MatrixMult(x_until_t, y);

        // compare that the result is equal to the according positions in x
        for (int i = 0; i < N - d - 1; i++)   // n-d-2 or n-d-1 ??
        {
            if ((y[i]) != (x[d + 1 + i])) {
                return false;
            }
        }
        return true;
    } else if (d == 2*T)
    {
        vector<FieldType> y(N - 1 - d); // the result of multiplication

        vector<FieldType> x_until_2t(2*T + 1);

        for (int i = 0; i < 2*T + 1; i++) {
            x_until_2t[i] = x[i];
        }

        matrix_for_2t.MatrixMult(x_until_2t, y);

        // compare that the result is equal to the according positions in x
        for (int i = 0; i < N - d - 1; i++)   // n-d-2 or n-d-1 ??
        {
            if ((y[i]) != (x[d + 1 + i])) {
                return false;
            }
        }
        return true;

    } else {
        vector<FieldType> alpha_until_d(d + 1);
        vector<FieldType> alpha_from_d(N - 1 - d);
        vector<FieldType> x_until_d(d + 1);
        vector<FieldType> y(N - 1 - d); // the result of multiplication

        for (int i = 0; i < d + 1; i++) {
            alpha_until_d[i] = alpha[i];
            x_until_d[i] = x[i];
        }
        for (int i = d + 1; i < N; i++) {
            alpha_from_d[i - (d + 1)] = alpha[i];
        }
        // Interpolate first d+1 positions of (alpha,x)
        HIM<FieldType> matrix(N - 1 - d, d + 1, field); // slices, only positions from 0..d
        matrix.InitHIMByVectors(alpha_until_d, alpha_from_d);
        matrix.MatrixMult(x_until_d, y);

        // compare that the result is equal to the according positions in x
        for (int i = 0; i < N - d - 1; i++)   // n-d-2 or n-d-1 ??
        {
            //if (field->elementToString(y[i]) != field->elementToString(x[d + 1 + i])) {
            if (y[i] != x[d + 1 + i]) {
                return false;
            }
        }
        return true;
    }
    return true;
}

// Interpolate polynomial at position Zero
template <class FieldType>
FieldType Protocol<FieldType>::interpolate(vector<FieldType> x)
{
    vector<FieldType> y(N); // result of interpolate
    matrix_for_interpolate.MatrixMult(x, y);
    return y[0];
}



template <class FieldType>
FieldType Protocol<FieldType>::reconstructShare(vector<FieldType>& x, int d, bool &flag){



    flag = true;

    if (!checkConsistency(x, d))
    {
        // someone cheated!

        cerr << "cheating!!! checkConsistency fained in reconstruct function" << endl;
        flag = false;

        return FieldType();
    }
    else {
        return interpolate(x);
    }
}



template <class FieldType>
void Protocol<FieldType>::DNHonestMultiplication(FieldType *a, FieldType *b, vector<FieldType> &cToFill, int numOfTrupples) {

    int index = 0;
    //int numOfMultGates = circuit.getNrOfMultiplicationGates();
    int fieldByteSize = field->getElementSizeInBytes();
    vector<FieldType> xyMinusRShares(numOfTrupples);//hold both in the same vector to send in one batch
    vector<byte> xyMinusRSharesBytes(numOfTrupples *fieldByteSize);//hold both in the same vector to send in one batch

    vector<FieldType> xyMinusR;//hold both in the same vector to send in one batch
    vector<byte> xyMinusRBytes;

    vector<vector<byte>> recBufsBytes;

    int offset = 0;


    //generate the shares for x+a and y+b. do it in the same array to send once
    for (int k = 0; k < numOfTrupples; k++)//go over only the logit gates
    {
        //compute the share of xy-r
        xyMinusRShares[k] = a[k]*b[k] - randomTAnd2TShares[offset + 2*k+1];

    }

    //set the acctual number of mult gate proccessed in this layer
    xyMinusRSharesBytes.resize(numOfTrupples*fieldByteSize);
    xyMinusR.resize(numOfTrupples);
    xyMinusRBytes.resize(numOfTrupples*fieldByteSize);

    for(int j=0; j<xyMinusRShares.size();j++) {
        field->elementToBytes(xyMinusRSharesBytes.data() + (j * fieldByteSize), xyMinusRShares[j]);
    }

    if (m_partyId == 0) {



        //just party 1 needs the recbuf
        recBufsBytes.resize(N);

        for (int i = 0; i < N; i++) {
            recBufsBytes[i].resize(numOfTrupples*fieldByteSize);
        }

        //receive the shares from all the other parties
        roundFunctionSyncForP1(xyMinusRSharesBytes, recBufsBytes);

    }
    else {//since I am not party 1 parties[0]->getID()=1

        //send the shares to p1
        parties[0]->getChannel()->write(xyMinusRSharesBytes.data(), xyMinusRSharesBytes.size());

    }

    //reconstruct the shares recieved from the other parties
    if (m_partyId == 0) {

        vector<FieldType> xyMinurAllShares(N), yPlusB(N);

        for (int k = 0;k < numOfTrupples; k++)//go over only the logit gates
        {
            for (int i = 0; i < N; i++) {

                xyMinurAllShares[i] = field->bytesToElement(recBufsBytes[i].data() + (k * fieldByteSize));
            }

            // reconstruct the shares by P0
            xyMinusR[k] = interpolate(xyMinurAllShares);

            //convert to bytes
            field->elementToBytes(xyMinusRBytes.data() + (k * fieldByteSize), xyMinusR[k]);

        }

        //send the reconstructed vector to all the other parties
        sendFromP1(xyMinusRBytes);
    }

    else {//each party get the xy-r reconstruced vector from party 1

        parties[0]->getChannel()->read(xyMinusRBytes.data(), xyMinusRBytes.size());
    }




    //fill the xPlusAAndYPlusB array for all the parties except for party 1 that already have this array filled
    if (m_partyId != 0) {

        for (int i = 0; i < numOfTrupples; i++) {

            xyMinusR[i] = field->bytesToElement(xyMinusRBytes.data() + (i * fieldByteSize));
        }
    }


    index = 0;

    //after the xPlusAAndYPlusB array is filled, we are ready to fill the output of the mult gates
    for (int k = 0; k < numOfTrupples; k++)//go over only the logit gates
    {
        cToFill[k] = randomTAnd2TShares[offset + 2*k] + xyMinusR[k];
    }

}

//
//template <class FieldType>
//int Protocol<FieldType>::processMultGRR() {
//    int index = 0;
//
//    vector<FieldType> x1(N),y1(N);
//
//    vector<vector<FieldType>> sendBufsElements(N);
//    vector<vector<byte>> sendBufsBytes(N);
//    vector<vector<byte>> recBufsBytes(N);
//    vector<FieldType> valBuf(circuit.getLayers()[currentCirciutLayer + 1] - circuit.getLayers()[currentCirciutLayer]); // Buffers for differences
//    FieldType d;
//
//    for(int i=0; i < N; i++)
//    {
//        //sendBufs[i] = "";
//
//        sendBufsElements[i].resize(
//                circuit.getLayers()[currentCirciutLayer + 1] - circuit.getLayers()[currentCirciutLayer]);
//        sendBufsBytes[i].resize((circuit.getLayers()[currentCirciutLayer + 1] - circuit.getLayers()[currentCirciutLayer]) *
//                                field->getElementSizeInBytes());
//        recBufsBytes[i].resize((circuit.getLayers()[currentCirciutLayer + 1] - circuit.getLayers()[currentCirciutLayer]) *
//                               field->getElementSizeInBytes());
//    }
//
//
//    for(int k = circuit.getLayers()[currentCirciutLayer]; k < circuit.getLayers()[currentCirciutLayer + 1] ; k++)//go over only the logit gates
//    {
//        // its a multiplication which not yet processed and ready
//        if(circuit.getGates()[k].gateType == MULT )
//        {
//            //set the secret of the polynomial to be the multiplication of the shares
//            x1[0] = gateShareArr[circuit.getGates()[k].input1] * gateShareArr[circuit.getGates()[k].input2];
//
//            // generate random degree-T polynomial
//            for(int i = 1; i < T + 1; i++)
//            {
//                // A random field element, uniform distribution
//                x1[i] = field->Random();
//
//            }
//
//
//            matrix_vand.MatrixMult(x1, y1, T + 1); // eval poly at alpha-positions
//
//            // prepare shares to be sent
//            for(int i=0; i < N; i++)
//            {
//                //cout << "y1[ " <<i<< "]" <<y1[i] << endl;
//                sendBufsElements[i][index] = y1[i];
//            }
//            index++;
//
//        }
//
//
//    }
//
//    //convert to bytes
//    int fieldByteSize = field->getElementSizeInBytes();
//    for(int i=0; i < N; i++)
//    {
//        for(int j=0; j<sendBufsElements[i].size();j++) {
//            field->elementToBytes(sendBufsBytes[i].data() + (j * fieldByteSize), sendBufsElements[i][j]);
//        }
//    }
//
//    roundFunctionSync(sendBufsBytes, recBufsBytes, 4);
//
//    int fieldBytesSize = field->getElementSizeInBytes();
//    index = 0;
//    for(int k = circuit.getLayers()[currentCirciutLayer]; k < circuit.getLayers()[currentCirciutLayer + 1] ; k++) {
//
//        if(circuit.getGates()[k].gateType == MULT) {
//            // generate random degree-T polynomial
//            for (int i = 0; i < N; i++) {
//                x1[i] = field->bytesToElement(recBufsBytes[i].data() + (index * fieldBytesSize));
//
//
//            }
//
//            FieldType accum = *field->GetZero();
//            for (int i = 0; i < N; i++) {
//
//                accum += firstRowVandInverse[i] * x1[i];
//
//            }
//
//            gateShareArr[circuit.getGates()[k].output] = accum;
//
//            index++;
//        }
//
//    }
//
//    return index;
//}

template <class FieldType>
void Protocol<FieldType>::GRRHonestMultiplication(FieldType *a, FieldType *b, vector<FieldType> &cToFill, int numOfTrupples)
{

    vector<FieldType> x1(N),y1(N);

    vector<vector<FieldType>> sendBufsElements(N);
    vector<vector<byte>> sendBufsBytes(N);
    vector<vector<byte>> recBufsBytes(N);

    FieldType d;

    vector<FieldType> ReconsBuf(numOfTrupples);


    for(int i=0; i < N; i++)
    {
        sendBufsElements[i].resize(numOfTrupples);
        sendBufsBytes[i].resize((numOfTrupples)*field->getElementSizeInBytes());
        recBufsBytes[i].resize((numOfTrupples)*field->getElementSizeInBytes());
    }


    for(int k = 0; k < numOfTrupples ; k++)//go over only the logit gates
    {

        //set the secret of the polynomial to be the multiplication of the shares
        x1[0] = a[k] * b[k];

        // generate random degree-T polynomial
        for(int i = 1; i < T+1; i++)
        {
            // A random field element, uniform distribution
            x1[i] = field->Random();

        }

        matrix_vand.MatrixMult(x1, y1, T+1); // eval poly at alpha-positions

        // prepare shares to be sent
        for(int i=0; i < N; i++)
        {
            //cout << "y1[ " <<i<< "]" <<y1[i] << endl;
            sendBufsElements[i][k] = y1[i];
        }

    }

    //convert to bytes
    int fieldByteSize = field->getElementSizeInBytes();
    for(int i=0; i < N; i++)
    {
        for(int j=0; j<sendBufsElements[i].size();j++) {
            field->elementToBytes(sendBufsBytes[i].data() + (j * fieldByteSize), sendBufsElements[i][j]);
        }
    }

    roundFunctionSync(sendBufsBytes, recBufsBytes,4);

    int fieldBytesSize = field->getElementSizeInBytes();

    for(int k = 0; k < numOfTrupples ; k++) {

       // generate random degree-T polynomial
        for (int i = 0; i < N; i++) {
            x1[i] = field->bytesToElement(recBufsBytes[i].data() + (k * fieldBytesSize));
        }

        FieldType accum = *field->GetZero();
        for (int i = 0; i < N; i++) {

            accum += firstRowVandInverse[i] * x1[i];

        }

        cToFill[k] = accum;
    }

}





template <class FieldType>
void Protocol<FieldType>::offlineDNForMultiplication(int numOfTriples){

   generateRandom2TAndTShares(numOfTriples,randomTAnd2TShares);
}


template <class FieldType>
bool Protocol<FieldType>::verifyTriples() {


    //int iterations =   (5 + field->getElementSizeInBytes() - 1) / field->getElementSizeInBytes();//calc the number of times we need to run the verification -- ceiling


    int numOfRandomelements = numOfOpens;//the number of random elements to create


    //first generate the common aes key
    auto key = generateCommonKey();

    bool flag;

//   if (verifyType == "Batch") {
//
//        vector<FieldType> betaElements(numOfRandomelements);
//        generatePseudoRandomElements(key, betaElements, numOfRandomelements);
//
//
//
//            answer = verificationOfBatchedTriples(randomABShares.data(), randomABShares.data() + 2*numOfOpens , c.data(),
//                                                   randomABShares.data() + numOfOpens, randomABShares.data() + 3*numOfOpens,c.data() + numOfOpens,
//                                                  betaElements.data(), numOfOpens);
//
//    }
//    else if (verifyType == "Single")
//    {

        vector<FieldType> betaElements(numOfRandomelements*2);
        generatePseudoRandomElements(key, betaElements, numOfRandomelements*2);



        flag = verificationOfSingleTriples(randomABShares.data(), randomABShares.data() + 2*numOfOpens , c.data(),
                                    randomABShares.data() + numOfOpens, randomABShares.data() + 3*numOfOpens,c.data() + numOfOpens,
                                    betaElements.data(), numOfOpens);


//    }


    return flag;

    if (flag_print) {
        cout << "answer is:" << flag << endl;
    }
}


template <class FieldType>
bool Protocol<FieldType>::verifyMults(int numOfpairs, vector<FieldType> &xShares, vector<FieldType> &yShares, vector<FieldType> &outputXYShares) {


    //int iterations =   (5 + field->getElementSizeInBytes() - 1) / field->getElementSizeInBytes();//calc the number of times we need to run the verification -- ceiling


    int numOfRandomelements = numOfpairs;//the number of random elements to create


    //first generate the common aes key
    auto key = generateCommonKey();

    bool flag;

    vector<FieldType> betaElements(numOfRandomelements*2);
    generatePseudoRandomElements(key, betaElements, numOfRandomelements*2);



    flag = verificationOfSingleTriples(xShares.data(), yShares.data(), outputXYShares.data(),
                                       randomABSharesForMult.data()+multIndex,
                                       randomABSharesForMult.data() + numOfMults + multIndex,
                                       cForMult.data()+multIndex,
                                       betaElements.data(), numOfpairs);




    if (flag_print) {
        cout << "answer is:" << flag << endl;
    }

    multIndex +=numOfpairs;
    return flag;

}



  template <class FieldType>
  vector<byte> Protocol<FieldType>::generateCommonKey(){

      int fieldByteSize = field->getElementSizeInBytes();

      //calc the number of elements needed for 128 bit AES key
      int numOfRandomShares = 16/field->getElementSizeInBytes() + 1;
      vector<FieldType> randomSharesArray(numOfRandomShares);
      vector<FieldType> aesArray(numOfRandomShares);
      vector<byte> aesKey(numOfRandomShares*fieldByteSize);


      //generate enough random shares for the AES key
      generateRandomShares(numOfRandomShares, randomSharesArray);

      openShare(numOfRandomShares, randomSharesArray, aesArray);


      //turn the aes array into bytes to get the common aes key.
      for(int i=0; i<numOfRandomShares;i++){

          for(int j=0; j<numOfRandomShares;j++) {
              field->elementToBytes(aesKey.data() + (j * fieldByteSize), aesArray[j]);
          }
      }

      //reduce the size of the key to 16 bytes
      aesKey.resize(16);

      return aesKey;

  }

  template <class FieldType>
  bool Protocol<FieldType>::openShare(int numOfRandomShares, vector<FieldType> &Shares, vector<FieldType> &secrets){

      bool flag;

      vector<vector<byte>> sendBufsBytes(N);
      vector<vector<byte>> recBufsBytes(N);

      vector<FieldType> x1(N);
      int fieldByteSize = field->getElementSizeInBytes();

      //calc the number of elements needed for 128 bit AES key

      //resize vectors
      for(int i=0; i < N; i++)
      {
          sendBufsBytes[i].resize(numOfRandomShares*fieldByteSize);
          recBufsBytes[i].resize(numOfRandomShares*fieldByteSize);
      }

      //set the first sending data buffer
      for(int j=0; j<numOfRandomShares;j++) {
          field->elementToBytes(sendBufsBytes[0].data() + (j * fieldByteSize), Shares[j]);
      }

      //copy the same data for all parties
      for(int i=1; i<N; i++){

          sendBufsBytes[i] = sendBufsBytes[0];
      }

      //call the round function to send the shares to all the users and get the other parties share
      roundFunctionSync(sendBufsBytes, recBufsBytes,12);

      //reconstruct each set of shares to get the secret

      for(int k=0; k<numOfRandomShares; k++){

          //get the set of shares for each element
          for(int i=0; i < N; i++) {

              x1[i] = field->bytesToElement(recBufsBytes[i].data() + (k*fieldByteSize));
          }

          secrets[k] = reconstructShare(x1, T, flag);

          if(flag==false)
            return false;

      }

      return true;

  }

template <class FieldType>
bool Protocol<FieldType>::multShares(int numOfpairs, vector<FieldType> &xShares, vector<FieldType> &yShares, vector<FieldType> &outputXYShares){


    //check that the arrays are in the right size
    if(xShares.size()!=numOfpairs || yShares.size()!=numOfpairs || outputXYShares.size()!=numOfpairs){
        return false;
    }
    //call the multiplication protocol

    //mult x shares by y shares
    honestMult->mult(xShares.data(), yShares.data(), outputXYShares, numOfpairs);

    //verify the shares
    return verifyMults(numOfpairs, xShares, yShares, outputXYShares);

}


template <class FieldType>
void Protocol<FieldType>::generatePseudoRandomElements(vector<byte> & aesKey, vector<FieldType> &randomElementsToFill, int numOfRandomElements){


    int fieldSize = field->getElementSizeInBytes();
    int fieldSizeBits = field->getElementSizeInBits();
    bool isLongRandoms;
    int size;

    if(fieldSize<=4 ){

        isLongRandoms = false;
        size = 4;
    }
    else if(fieldSize>4 && fieldSize<=8){
      isLongRandoms = true;
      size = 8;
    }
    else{
        size = 16;

    }


    if (flag_print) {
        cout << "size is" << size << "for party : " << m_partyId;
    }


    PrgFromOpenSSLAES prg((numOfRandomElements*size/16) + 1);
    SecretKey sk(aesKey, "aes");
    prg.setKey(sk);


    vector<byte> vec(fieldSize);

    for(int i=0; i<numOfRandomElements; i++){

        prg.getPRGBytes(vec, 0, fieldSize);

        randomElementsToFill[i] = field->bytesToElement(vec.data());

    }

}

template <class FieldType>
bool Protocol<FieldType>::verificationOfSingleTriples(FieldType *x, FieldType *y, FieldType *z,
                                                       FieldType *a, FieldType *b, FieldType *c,
                                                       FieldType * randomElements, int numOfTriples){


    bool flag;

    vector<FieldType> roAndSigma(numOfTriples*2);//vector holding the random shares generated
    vector<FieldType> secretArr(numOfTriples*2);
    vector<FieldType> secretArrResultsV(numOfTriples);

    vector<FieldType> v(numOfTriples);//vector holding the 4*numOfTriples output of the multiplication

    //prepare the 4k pairs for multiplication
    for(int k=0; k<numOfTriples; k++){

        roAndSigma[2*k] = randomElements[k]*x[k] + a[k];//ro
        roAndSigma[2*k + 1] = y[k] + b[k];//sigma
    }

     //open all the shares at once
    flag = openShare(numOfTriples*2, roAndSigma, secretArr);

    if(flag==false)
        return false;

    vector<FieldType> shareArrV(numOfTriples);


    FieldType VShare(0);
    //compute the output share array
    for(int k=0; k<numOfTriples; k++){
        VShare += (randomElements[k]*z[k] - c[k] + secretArr[2*k+1]*a[k]+ secretArr[2*k]*b[k] - secretArr[2*k] * secretArr[2*k + 1])*randomElements[numOfTriples + k];
    }


    //open [V]
    vector<FieldType> shareArr(1);
    vector<FieldType> secretArrV(1);
    shareArr[0] = VShare;

    flag = openShare(1,shareArr,secretArrV);

    if(flag==false)
        return false;

    //check that V=0
    if(secretArrV[0] != *field->GetZero()) {

        cout<<"-------verification failed-----------"<<endl;
        return false;
    }

    return true;

}


template <class FieldType>
void Protocol<FieldType>::roundFunctionSync(vector<vector<byte>> &sendBufs, vector<vector<byte>> &recBufs, int round) {

    //cout<<"in roundFunctionSync "<< round<< endl;

    int numThreads = 10;//parties.size();
    int numPartiesForEachThread;

    if (parties.size() <= numThreads){
        numThreads = parties.size();
        numPartiesForEachThread = 1;
    } else{
        numPartiesForEachThread = (parties.size() + numThreads - 1)/ numThreads;
    }


    recBufs[m_partyId] = move(sendBufs[m_partyId]);
    //recieve the data using threads
    vector<thread> threads(numThreads);
    for (int t=0; t<numThreads; t++) {
        if ((t + 1) * numPartiesForEachThread <= parties.size()) {
            threads[t] = thread(&Protocol::exchangeData, this, ref(sendBufs), ref(recBufs),
                                t * numPartiesForEachThread, (t + 1) * numPartiesForEachThread);
        } else {
            threads[t] = thread(&Protocol::exchangeData, this, ref(sendBufs), ref(recBufs), t * numPartiesForEachThread, parties.size());
        }
    }
    for (int t=0; t<numThreads; t++){
        threads[t].join();
    }

}


template <class FieldType>
void Protocol<FieldType>::exchangeData(vector<vector<byte>> &sendBufs, vector<vector<byte>> &recBufs, int first, int last){


    //cout<<"in exchangeData";
    for (int i=first; i < last; i++) {

        if ((m_partyId) < parties[i]->getID()) {


            if (sendBufs[parties[i]->getID()].size() > 0) {
                //send shares to my input bits
                parties[i]->getChannel()->write(sendBufs[parties[i]->getID()].data(),
                                                sendBufs[parties[i]->getID()].size());
                //cout<<"write the data:: my Id = " << m_partyId - 1<< "other ID = "<< parties[i]->getID() <<endl;
            }

            if (recBufs[parties[i]->getID()].size() > 0) {
                //receive shares from the other party and set them in the shares array
                parties[i]->getChannel()->read(recBufs[parties[i]->getID()].data(),
                                               recBufs[parties[i]->getID()].size());
                //cout<<"read the data:: my Id = " << m_partyId-1<< "other ID = "<< parties[i]->getID()<<endl;
            }

        } else{

            if (recBufs[parties[i]->getID()].size() > 0) {
                //receive shares from the other party and set them in the shares array
                parties[i]->getChannel()->read(recBufs[parties[i]->getID()].data(),
                                               recBufs[parties[i]->getID()].size());
                //cout<<"read the data:: my Id = " << m_partyId-1<< "other ID = "<< parties[i]->getID()<<endl;
            }

            if (sendBufs[parties[i]->getID()].size() > 0) {

                //send shares to my input bits
                parties[i]->getChannel()->write(sendBufs[parties[i]->getID()].data(),
                                                sendBufs[parties[i]->getID()].size());
                //cout<<"write the data:: my Id = " << m_partyId-1<< "other ID = "<< parties[i]->getID() <<endl;
            }

        }

    }


}


template <class FieldType>
void Protocol<FieldType>::roundFunctionSyncForP1(vector<byte> &myShare, vector<vector<byte>> &recBufs) {

    //cout<<"in roundFunctionSyncBroadcast "<< endl;

    int numThreads = parties.size();
    int numPartiesForEachThread;

    if (parties.size() <= numThreads){
        numThreads = parties.size();
        numPartiesForEachThread = 1;
    } else{
        numPartiesForEachThread = (parties.size() + numThreads - 1)/ numThreads;
    }


    recBufs[m_partyId] = myShare;
    //recieve the data using threads
    vector<thread> threads(numThreads);
    for (int t=0; t<numThreads; t++) {
        if ((t + 1) * numPartiesForEachThread <= parties.size()) {
            threads[t] = thread(&Protocol::recDataToP1, this,  ref(recBufs),
                                t * numPartiesForEachThread, (t + 1) * numPartiesForEachThread);
        } else {
            threads[t] = thread(&Protocol::recDataToP1, this, ref(recBufs), t * numPartiesForEachThread, parties.size());
        }
    }
    for (int t=0; t<numThreads; t++){
        threads[t].join();
    }

}


template <class FieldType>
void Protocol<FieldType>::recDataToP1(vector<vector<byte>> &recBufs, int first, int last){


    //cout<<"in exchangeData";
    for (int i=first; i < last; i++) {

        parties[i]->getChannel()->read(recBufs[parties[i]->getID()].data(), recBufs[parties[i]->getID()].size());
        //cout<<"read the data:: my Id = " << m_partyId-1<< "other ID = "<< parties[i]->getID()<<endl;
    }


}



template <class FieldType>
void Protocol<FieldType>::sendFromP1(vector<byte> &sendBuf) {

    //cout<<"in roundFunctionSyncBroadcast "<< endl;

    int numThreads = parties.size();
    int numPartiesForEachThread;

    if (parties.size() <= numThreads){
        numThreads = parties.size();
        numPartiesForEachThread = 1;
    } else{
        numPartiesForEachThread = (parties.size() + numThreads - 1)/ numThreads;
    }

    //recieve the data using threads
    vector<thread> threads(numThreads);
    for (int t=0; t<numThreads; t++) {
        if ((t + 1) * numPartiesForEachThread <= parties.size()) {
            threads[t] = thread(&Protocol::sendDataFromP1, this,  ref(sendBuf),
                                t * numPartiesForEachThread, (t + 1) * numPartiesForEachThread);
        } else {
            threads[t] = thread(&Protocol::sendDataFromP1, this, ref(sendBuf), t * numPartiesForEachThread, parties.size());
        }
    }
    for (int t=0; t<numThreads; t++){
        threads[t].join();
    }

}

template <class FieldType>
void Protocol<FieldType>::sendDataFromP1(vector<byte> &sendBuf, int first, int last){

    for(int i=first; i < last; i++) {

        parties[i]->getChannel()->write(sendBuf.data(), sendBuf.size());

    }


}

template <class FieldType>
Protocol<FieldType>::~Protocol()
{
}



#endif /* PROTOCOL_H_ */
