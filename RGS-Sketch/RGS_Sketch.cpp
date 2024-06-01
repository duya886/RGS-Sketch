#include <cstdint>
#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <ctime>
#include <iomanip>
#include <vector>
#include "MurmurHash3.h"
#include <fstream>
#include <string>

#define HASH_SEED 51318471
using namespace std;
#define maxRho 25
// int maxRho=min((int)31,(int)33-(int)log2(V_REG_GROUP_NUM*REG_GROUP_SIZE));

const uint32_t R[32]={2725628290,1444492164,1001369998,2291192083,3632694421,3217476757,962572696,3523028760,1564808096,1686744101,2572712615,186525479,1267589045,2944996661,3807606838,2892595386,1879920448,3800713281,4044617030,1088214349,2912533710,4048616015,223249368,1094862426,1603814236,3377932768,121561825,2119485156,1642325476,1361943285,3238649846,992916862};

template<uint32_t BUCKET_NUM,uint32_t BUCKET_SIZE,uint32_t REG_GROUP_NUM,uint32_t REG_GROUP_SIZE,uint32_t V_REG_GROUP_NUM>
class SKETCH{
private:
    vector<vector<pair<uint32_t,uint32_t>>> bucketArray;
    uint8_t* regArray;

    inline bool updateBucket(const pair<uint32_t,uint32_t>& pkt,uint32_t eleHashValue,uint32_t bktIdx,double updateValue){
        //check if the flow has been recorded
        uint32_t minCellIdx=-1,minCellVal=-1;
        for(int cellIdx=0;cellIdx<BUCKET_SIZE;cellIdx++){
            if(bucketArray[bktIdx][cellIdx].second==0){
                uint32_t newVal=floor(updateValue);

                double temp=((double)rand())/RAND_MAX;
                if(temp<(updateValue-newVal)){
                    newVal+=1;
                }

                bucketArray[bktIdx][cellIdx].first=pkt.first;
                bucketArray[bktIdx][cellIdx].second=newVal;
                return true;
            }
            else if(bucketArray[bktIdx][cellIdx].first==pkt.first){
                uint32_t newVal=floor(updateValue);
                double temp=((double)rand())/RAND_MAX;
                if(temp<(updateValue-newVal)){
                    newVal+=1;
                }

                bucketArray[bktIdx][cellIdx].second+=newVal;
                return true;
            }
            else{
                if(minCellIdx==-1 || bucketArray[bktIdx][cellIdx].second<minCellVal){
                    minCellIdx=cellIdx;
                    minCellVal=bucketArray[bktIdx][cellIdx].second;
                }
            }
        }

        uint64_t randVal=((uint64_t)eleHashValue)*(minCellVal+1);
        if(randVal<=4294967295){
            bucketArray[bktIdx][minCellIdx].first=pkt.first;
            bucketArray[bktIdx][minCellIdx].second=minCellVal+1;
            return true;
        }
        return false;
    }

    inline static uint32_t getRegValue(uint32_t groupStartByteIdx,uint32_t regIdxInGroup,uint8_t* regArray){
        uint32_t regStartByteIdx=regIdxInGroup*5/8;
        uint32_t regStartBitIdx=regIdxInGroup*5%8;

        return ((*(uint16_t*)(regArray+groupStartByteIdx+regStartByteIdx))>>regStartBitIdx)&0x1F;
    }
    inline static void setRegValue(uint32_t groupStartByteIdx,uint32_t regIdxInGroup,uint32_t newValue,uint8_t* regArray){
        uint32_t regStartByteIdx=regIdxInGroup*5/8;
        uint32_t regStartBitIdx=regIdxInGroup*5%8;
        uint16_t* value_ptr=(uint16_t*)(regArray+groupStartByteIdx+regStartByteIdx);
        *value_ptr=(*value_ptr)&(~(((uint16_t)0x1F)<<regStartBitIdx));
        *value_ptr=(*value_ptr)|(newValue<<regStartBitIdx);
    }

    static double getCoe(uint32_t regNum){
        double coe=0;
        switch(regNum){
            case 16:
                coe=0.673;
                break;
            case 32:
                coe=0.697;
                break;
            case 64:
                coe=0.709;
            default:
                if(regNum>=128){
                    coe=0.7213/(1.0+1.079/regNum);
                }
                break;
        }
        return coe;
    }

    static double getTotalSpread(uint8_t* regArray){
        uint32_t totalRegNum=REG_GROUP_NUM*REG_GROUP_SIZE;
        double coe= getCoe(totalRegNum);

        double sum=0.0;
        uint32_t numOfZeroRegs=0;
        for(uint32_t groupIdx=0;groupIdx<REG_GROUP_NUM;groupIdx++){
            uint32_t groupStartByteIdx= groupIdx*5*REG_GROUP_SIZE/8;
            for(uint32_t regIdxInGroup=0;regIdxInGroup<REG_GROUP_SIZE;regIdxInGroup++){
                uint32_t regValue= getRegValue(groupStartByteIdx,regIdxInGroup,regArray);
                sum+=1.0/pow(2,regValue);

                if(regValue==0){
                    numOfZeroRegs++;
                }
            }
        }

        double estimatedVal=coe*pow(totalRegNum,2)/sum;
        double totalSpread=estimatedVal;

        if(estimatedVal<2.5*totalRegNum){
            if(numOfZeroRegs>0){
                totalSpread=totalRegNum*log((double)totalRegNum/numOfZeroRegs);
            }
        }else if(estimatedVal>pow(2,32)/30.0){
            totalSpread=-1.0*pow(2,32)*log(1.0-(double)estimatedVal/pow(2,32));
        }

        return totalSpread;
    }

    static double estimateBasedOnRegs(uint32_t key, uint8_t* regArray,double totalSpread){
        uint32_t totalRegNum=REG_GROUP_NUM*REG_GROUP_SIZE;
        uint32_t flowRegNum=V_REG_GROUP_NUM*REG_GROUP_SIZE;
        double coe= getCoe(flowRegNum);

        double sum=0.0;
        uint32_t numOfZeroRegs=0;

        uint32_t hashValue1=0;
        MurmurHash3_x86_32(&key, 4,HASH_SEED, &hashValue1);

        for(uint32_t virGroupIdx=0;virGroupIdx<V_REG_GROUP_NUM;virGroupIdx++){
            uint32_t phyGroupIdx=(R[virGroupIdx]^hashValue1)%REG_GROUP_NUM;

            uint32_t groupStartByteIdx= phyGroupIdx*5*REG_GROUP_SIZE/8;
            for(uint32_t regIdxInGroup=0;regIdxInGroup<REG_GROUP_SIZE;regIdxInGroup++){
                uint32_t regValue= getRegValue(groupStartByteIdx,regIdxInGroup,regArray);
                sum+=1.0/pow(2,regValue);

                if(regValue==0){
                    numOfZeroRegs++;
                }
            }
        }

        double n_s=coe*pow(flowRegNum,2)/sum;
        if(n_s<2.5*flowRegNum){
            if(numOfZeroRegs>0){
                n_s=flowRegNum*log((double)flowRegNum/numOfZeroRegs);
            }
        }else if(n_s>pow(2,32)/30.0){
            n_s=-1.0*pow(2,32)*log(1.0-n_s/pow(2,32));
        }

        double n_f=(1.0*totalRegNum*flowRegNum/(totalRegNum-flowRegNum))*(n_s/flowRegNum-totalSpread/totalRegNum);
        if(n_f<0){
            n_f=0;
        }
        return n_f;
    }
public:
    SKETCH():bucketArray(BUCKET_NUM,vector<pair<uint32_t,uint32_t>>(BUCKET_SIZE)){
        regArray=new uint8_t[int(ceil(REG_GROUP_NUM*REG_GROUP_SIZE*5.0/8.0))]();
    }
    SKETCH(const SKETCH& c):bucketArray(c.bucketArray){
        uint32_t byteNum=int(ceil(REG_GROUP_NUM*REG_GROUP_SIZE*5.0/8.0));
        regArray=new uint8_t[byteNum]();
        for(uint32_t byteIdx=0;byteIdx<byteNum;byteIdx++){
            regArray[byteIdx]=c.regArray[byteIdx];
        }
    }
    void insert(const pair<uint32_t,uint32_t>& pkt){
        uint32_t hashValue1=0;
        MurmurHash3_x86_32(&pkt.first, 4,HASH_SEED, &hashValue1);

        uint32_t hashValue2=0;
        MurmurHash3_x86_32(&pkt.second, 4,HASH_SEED+15417891, &hashValue2);

        uint32_t tempHash1=hashValue1+hashValue2;
        uint32_t virGroupIdx=tempHash1%V_REG_GROUP_NUM;
        uint32_t regIdxInGroup=(tempHash1/V_REG_GROUP_NUM)%REG_GROUP_SIZE;


        //Note, do not use: min(maxRho,(int)__builtin_clz(tempHash1/(V_REG_GROUP_NUM*REG_GROUP_SIZE)+1)
        //since tempHash1/(V_REG_GROUP_NUM*REG_GROUP_SIZE) will be seen as a 32-bit string with the highest bits being 0
        uint32_t rhoValue=min(maxRho,(int)__builtin_clz(tempHash1)+1);

        uint32_t phyGroupIdx=(R[virGroupIdx]^hashValue1)%REG_GROUP_NUM;

        uint32_t groupStartByteIdx=phyGroupIdx*5*REG_GROUP_SIZE/8;
        uint32_t oriRhoValue= getRegValue(groupStartByteIdx,regIdxInGroup,regArray);


        if(rhoValue>oriRhoValue){
            uint64_t value=*(uint64_t*)(regArray+groupStartByteIdx);
            double sum=0.0;
            for(int i=0;i<REG_GROUP_SIZE;i++){
                uint32_t regVal=value&0x1F;
                if(regVal!=maxRho){//when the reg value is maxRho, it cannot be updated again
                    sum+=1<<(maxRho-1-regVal);
                }
                value=value>>5;
            }

            double updateValue= (double)(((uint64_t)REG_GROUP_SIZE)<<(maxRho-1))/sum;

            uint32_t bktIdx=hashValue1%BUCKET_NUM;
            if(updateBucket(pkt,(hashValue1+99985157*hashValue2),bktIdx,updateValue)){
                setRegValue(groupStartByteIdx,regIdxInGroup,rhoValue,regArray);
            }
        }
    }

    void getOnlineEstimatedFlowSpreads(unordered_map<uint32_t, uint32_t>& estimatedFlowSpreads) const{
        for (uint32_t bucketIdx = 0; bucketIdx < BUCKET_NUM; bucketIdx++) {
            for (int cellIdx = 0; cellIdx < BUCKET_SIZE; cellIdx++) {
                uint32_t key=bucketArray[bucketIdx][cellIdx].first;
                uint32_t estSpread=bucketArray[bucketIdx][cellIdx].second;

                estimatedFlowSpreads[key]=estSpread;
            }
        }
    }

    void getVHLLEstimatedFlowSpreads(unordered_map<uint32_t, uint32_t>& estimatedFlowSpreads) const{
        double totalSpread= getTotalSpread(regArray);

        for (uint32_t bucketIdx = 0; bucketIdx < BUCKET_NUM; bucketIdx++) {
            for (int cellIdx = 0; cellIdx < BUCKET_SIZE; cellIdx++) {
                uint32_t key=bucketArray[bucketIdx][cellIdx].first;

                double estSpread= estimateBasedOnRegs(key,regArray,totalSpread);

                estimatedFlowSpreads[key]=round(estSpread);
            }
        }
    }

    static void getMergedEstimatedFlowSpreads(const vector<SKETCH>& sketches,unordered_map<uint32_t, uint32_t>& estimatedFlowSpreads){
        uint8_t* newRegArray=new uint8_t[int(ceil(REG_GROUP_NUM*REG_GROUP_SIZE*5.0/8.0))]();

        vector<uint32_t> flowSpreadsMinVals;
        vector<unordered_map<uint32_t, uint32_t>> estimatedFlowSpreadsOfMultipleSketches;
        unordered_map<uint32_t, uint32_t> estimatedFlowSpreadsMaxVal;//get the max estimated spreads of flows among all measurement points, which are used as the lower bound for merged estimation

        for(uint32_t sketchIdx=0;sketchIdx<sketches.size();sketchIdx++){
            unordered_map<uint32_t, uint32_t> estimatedFlowSpreadsOfASketch;
            sketches[sketchIdx].getOnlineEstimatedFlowSpreads(estimatedFlowSpreadsOfASketch);

            uint32_t minVal=UINT32_MAX;
            for(auto iter=estimatedFlowSpreadsOfASketch.begin();iter!=estimatedFlowSpreadsOfASketch.end();iter++){
                uint32_t key=iter->first;
                uint32_t estVal=iter->second;
                estimatedFlowSpreads[key]=0;

                auto iter2=estimatedFlowSpreadsMaxVal.find(key);
                if(iter2!=estimatedFlowSpreadsMaxVal.end()){
                    if(estVal>iter2->second){
                        iter2->second=estVal;
                    }
                }else{
                    estimatedFlowSpreadsMaxVal[key]=estVal;
                }

                if(estVal<minVal){
                    minVal=estVal;
                }
            }
            flowSpreadsMinVals.push_back(minVal);
            estimatedFlowSpreadsOfMultipleSketches.push_back(std::move(estimatedFlowSpreadsOfASketch));

            for(uint32_t groupIdx=0;groupIdx<REG_GROUP_NUM;groupIdx++){
                uint32_t groupStartByteIdx= groupIdx*5*REG_GROUP_SIZE/8;
                for(uint32_t regIdxInGroup=0;regIdxInGroup<REG_GROUP_SIZE;regIdxInGroup++){
                    uint32_t curRegValue= getRegValue(groupStartByteIdx,regIdxInGroup,newRegArray);
                    uint32_t sketchRegValue= getRegValue(groupStartByteIdx,regIdxInGroup,sketches[sketchIdx].regArray);
                    if(sketchRegValue>curRegValue){
                        setRegValue(groupStartByteIdx,regIdxInGroup,sketchRegValue,newRegArray);
                    }
                }
            }
        }

        unordered_map<uint32_t, uint32_t> estimatedFlowSpreadsSumVal;
        for(auto iter=estimatedFlowSpreads.begin();iter!=estimatedFlowSpreads.end();iter++){
            uint32_t key=iter->first;

            estimatedFlowSpreadsSumVal[key]=0;

            for(uint32_t sketchIdx=0;sketchIdx<sketches.size();sketchIdx++){
                auto iter2=estimatedFlowSpreadsOfMultipleSketches[sketchIdx].find(key);
                if(iter2!=estimatedFlowSpreadsOfMultipleSketches[sketchIdx].end()){
                    estimatedFlowSpreadsSumVal[key]+=estimatedFlowSpreadsOfMultipleSketches[sketchIdx][key];
                }else{
                    estimatedFlowSpreadsSumVal[key]+=flowSpreadsMinVals[sketchIdx];
                }
            }
        }

        double totalSpread= getTotalSpread(newRegArray);


        for(auto iter=estimatedFlowSpreads.begin();iter!=estimatedFlowSpreads.end();iter++) {
            uint32_t key = iter->first;

            double n_f= estimateBasedOnRegs(key,newRegArray,totalSpread);
            if(n_f<estimatedFlowSpreadsMaxVal[key]){
                n_f=estimatedFlowSpreadsMaxVal[key];
            }else if(n_f>estimatedFlowSpreadsSumVal[key]){
                n_f=estimatedFlowSpreadsSumVal[key];
            }
            estimatedFlowSpreads[key]=round(n_f);
        }

        delete []newRegArray;
    }

    ~SKETCH(){
        delete []regArray;
    }
};


bool cmpPairFunc(pair<uint32_t,uint32_t>p1, pair<uint32_t,uint32_t>p2)
{
    return p1.second > p2.second;
}

vector<vector<double>> calculateMetrics(unordered_map<uint32_t, uint32_t> &estFlowSpreads,unordered_map<uint32_t, uint32_t> &actualFlowSpreads,vector<uint32_t> topKs, string outputDirPath,string savedFileName,string info){
    vector<pair<uint32_t,uint32_t>> estFlowSpreadsVec;
    for(auto iter=estFlowSpreads.begin();iter!=estFlowSpreads.end();iter++){
        estFlowSpreadsVec.push_back(make_pair(iter->first,iter->second));
    }
    sort(estFlowSpreadsVec.begin(), estFlowSpreadsVec.end(), cmpPairFunc);

    vector<pair<uint32_t,uint32_t>> actualFlowSpreadsVec;
    for(auto iter=actualFlowSpreads.begin();iter!=actualFlowSpreads.end();iter++){
        actualFlowSpreadsVec.push_back(make_pair(iter->first,iter->second));
    }
    sort(actualFlowSpreadsVec.begin(), actualFlowSpreadsVec.end(), cmpPairFunc);

    string resultFilePath = outputDirPath+savedFileName + ".result";
    ofstream resultFile;
    resultFile.open(resultFilePath,ios::out);

    resultFile<<"flowId\t\tacutalSpread\t\testSpread\t\tsort by est"<<endl;


    for(uint32_t i=0;i<estFlowSpreadsVec.size();i++){
        uint32_t key=estFlowSpreadsVec[i].first;
        uint32_t estSpread=estFlowSpreadsVec[i].second;
        uint32_t actualSpread=1.0;

        auto iter=actualFlowSpreads.find(key);
        if(iter!=actualFlowSpreads.end()){
            actualSpread=actualFlowSpreads[key];
        }

        resultFile<<key<<"\t\t"<<actualSpread<<"\t\t"<<estSpread<<endl;
    }
    resultFile.close();

    string metricsFilePath = outputDirPath+savedFileName + ".metrics";
    ofstream metricsFile;
    metricsFile.open(metricsFilePath,ios::app);

    metricsFile<<"******************"<<endl;
    metricsFile<<"Super Spreader Metrics: "<<info<<endl;
    metricsFile<<"******************"<<endl;
    //trueFlowNum means the true number of super spreaders with that threshold
    //fakeFlowNum means the number of fake flows that returned by reversible calculation operation based on Chinese Remainder Theorem, in other words, these flows do not exist in the packet stream.
    metricsFile<<"Threshold\t\tPR\tRC\tF1\tARE\tAAE\tTP\tFP\tFN\ttrueFlowNum\t\tfakeFlowNum"<<endl;

    vector<double> THs,TPs,FPs,FNs,PRs,RCs,F1s,AREs,AAEs,realFlowNums,fakeFlowNums;
    double fakeFlowNum=0;

    for(uint32_t k:topKs){
        double threshold=0.5*(actualFlowSpreadsVec[k-1].second+actualFlowSpreadsVec[k].second);

        unordered_set<uint32_t> trueSuperSpreads;
        for(uint32_t i=0;i<actualFlowSpreadsVec.size();i++){
            if(actualFlowSpreadsVec[i].second>=threshold){
                trueSuperSpreads.insert(actualFlowSpreadsVec[i].first);
            }
        }

        double TP=0,FP=0,totalRE=0,totalAE=0;

        for(uint32_t i=0;i<estFlowSpreadsVec.size();i++){
            uint32_t key=estFlowSpreadsVec[i].first;
            double estSpread=estFlowSpreadsVec[i].second;

            if(estSpread>=threshold){
                auto iter=trueSuperSpreads.find(key);
                if(iter!=trueSuperSpreads.end()){
                    TP+=1;
                }else{
                    FP+=1;
                }

                double actualSpread=1.0;
                auto iter2=actualFlowSpreads.find(key);
                if(iter2!=actualFlowSpreads.end()){
                    actualSpread=actualFlowSpreads[key];
                }else{
                    fakeFlowNum+=1;
                }

                totalAE+=abs(estSpread-actualSpread);
                totalRE+=abs((estSpread-actualSpread)/actualSpread);
            }else{
                break;
            }
        }

        double FN=trueSuperSpreads.size()-TP;
        double PR=TP/(TP+FP);
        double RC=TP/(TP+FN);
        double F1=2*PR*RC/(PR+RC);
        double ARE=totalRE/(TP+FP);
        double AAE=totalAE/(TP+FP);

        char temp[500]{0};
        sprintf(temp,"%-10.3f\t%-10.3f\t%-10.3f\t%-10.3f\t%-10.3f\t%-10.3f\t%-6d\t%-6d\t%-6d\t%-6d\t\t%-6d",threshold,PR,RC,F1,ARE,AAE,int(TP),int(FP),int(FN),trueSuperSpreads.size(),int(fakeFlowNum));
        metricsFile<<temp<<endl;

        THs.push_back(threshold);
        PRs.push_back(PR);
        RCs.push_back(RC);
        F1s.push_back(F1);
        AREs.push_back(ARE);
        AAEs.push_back(AAE);
        TPs.push_back(TP);
        FPs.push_back(FP);
        FNs.push_back(FN);
        realFlowNums.push_back(trueSuperSpreads.size());
        fakeFlowNums.push_back(fakeFlowNum);
    }

    metricsFile.close();

    vector<vector<double>> metrics={THs,PRs,RCs,F1s,AREs,AAEs,TPs,FPs,FNs,realFlowNums,fakeFlowNums};
    return metrics;
}

//for single measurement point
unsigned int ReadInSingleTraceFile(const char* traceFilePath, vector<pair<uint32_t,uint32_t>>& pkts, unordered_map<uint32_t, uint32_t>& actualFlowSpreads)
{
    unordered_map<uint32_t, unordered_set<uint32_t>> allFlowElements;

    FILE* fin = fopen(traceFilePath, "rb");

    uint64_t aPkt;
    unsigned int count = 0;

    while (fread(&aPkt, 8, 1, fin) == 1) {
        uint32_t src=aPkt&0xFFFFFFFF;
        uint32_t dst=aPkt>>32;

        pkts.push_back(make_pair(src,dst));
        allFlowElements[src].insert(dst);

        count++;
        if (count % 5000000 == 0) {
            printf("Successfully read in %s, %u packets\n", traceFilePath, count);

        }
    }
    printf("Successfully read in %s, %u packets\n", traceFilePath, count);
    fclose(fin);

    for(auto iter=allFlowElements.begin();iter!=allFlowElements.end();iter++){
        uint32_t key=iter->first;
        uint32_t acutalSpread=iter->second.size();
        actualFlowSpreads[key]=acutalSpread;
    }
    printf("Successfully get actual flow spreads\n");

    return count;
}
//for multiple measurement points
unsigned int ReadInSingleTraceFile(const char* traceFilePath, vector<pair<uint32_t,uint32_t>>& pkts, unordered_map<uint32_t, uint32_t>& actualFlowSpreads, unordered_map<uint32_t, unordered_set<uint32_t>>& totalFlowElements)
{
    unordered_map<uint32_t, unordered_set<uint32_t>> allFlowElements;

    FILE* fin = fopen(traceFilePath, "rb");

    uint64_t aPkt;
    unsigned int count = 0;

    while (fread(&aPkt, 8, 1, fin) == 1) {
        uint32_t src=aPkt&0xFFFFFFFF;
        uint32_t dst=aPkt>>32;

        pkts.push_back(make_pair(src,dst));
        allFlowElements[src].insert(dst);
        totalFlowElements[src].insert(dst);

        count++;
        if (count % 5000000 == 0) {
            printf("Successfully read in %s, %u packets\n", traceFilePath, count);
        }
    }
    printf("Successfully read in %s, %u packets\n", traceFilePath, count);
    fclose(fin);

    for(auto iter=allFlowElements.begin();iter!=allFlowElements.end();iter++){
        uint32_t key=iter->first;
        uint32_t acutalSpread=iter->second.size();
        actualFlowSpreads[key]=acutalSpread;
    }
    printf("Successfully get actual flow spreads\n");

    return count;
}

template<uint32_t BUCKET_NUM,uint32_t BUCKET_SIZE,uint32_t REG_GROUP_NUM,uint32_t REG_GROUP_SIZE,uint32_t V_REG_GROUP_NUM>
void processMultipleTracesAndGetAvgMetrics(const char* traceFileDir,uint32_t fileNum,const vector<uint32_t>& topKs, string outputDirPath,string traceInfo){
    uint32_t totalMem=BUCKET_NUM*BUCKET_SIZE*64+REG_GROUP_NUM*5*REG_GROUP_SIZE;
    double bucketMemRatio=(double)BUCKET_NUM*BUCKET_SIZE*64/totalMem;
    cout << "totalMem:" << totalMem/(8*1024.0) << "KB" << endl;
    cout << "*********************" << endl;

    char temp[500]{0};
    sprintf(temp,"%s data %d files src RGS-Sketch v_GNum=%d rGNum=%d rGSz=%d bktNum=%d bktSz=%d bktR=%.2f mem=%.2fKB",traceInfo.c_str(),fileNum,V_REG_GROUP_NUM,REG_GROUP_NUM,REG_GROUP_SIZE,BUCKET_NUM,BUCKET_SIZE,bucketMemRatio,totalMem/(8.0*1024));
    string savedFileName(temp);

    double avgThroughput=0;
    vector<double> THs(topKs.size()),TPs(topKs.size()),FPs(topKs.size()),FNs(topKs.size()),PRs(topKs.size());
    vector<double> RCs(topKs.size()),F1s(topKs.size()),AREs(topKs.size()),AAEs(topKs.size()),realFlowNums(topKs.size()),fakeFlowNums(topKs.size());
    vector<vector<double>> avgMetrics={THs,PRs,RCs,F1s,AREs,AAEs,TPs,FPs,FNs,realFlowNums,fakeFlowNums};

    for(uint32_t i=0;i<fileNum;i++){
        char traceFilePath[256]{0};
        sprintf(traceFilePath,"%s%02d.dat",traceFileDir,i);

        //prepare dataset
        cout << "prepare dataset: " <<traceFilePath<< endl;
        vector<pair<uint32_t,uint32_t>> pkts;
        unordered_map<uint32_t, uint32_t> actualFlowSpreads;
        unsigned int actualItemNum = ReadInSingleTraceFile(traceFilePath, pkts, actualFlowSpreads);

        cout << "number of packets: " << actualItemNum << endl;
        cout << "number of flows: " << actualFlowSpreads.size() << endl;
        cout << "*********************" << endl;


        cout << "prepare algorithm "<< endl;
        SKETCH<BUCKET_NUM,BUCKET_SIZE,REG_GROUP_NUM,REG_GROUP_SIZE,V_REG_GROUP_NUM> sketch;

        clock_t time1 = clock();
        for (unsigned int i = 0; i < pkts.size(); i++) {
            sketch.insert(pkts[i]);
        }
        clock_t time2 = clock();

        double numOfSeconds = (double)(time2 - time1) / CLOCKS_PER_SEC;//the seconds using to insert items
        double throughput = (actualItemNum / 1000000.0) / numOfSeconds;
        cout << "throughput: " << throughput << " Mpps" << ", each insert operation uses " << 1000.0 / throughput << " ns" << endl;
        cout << "*********************" << endl;

        unordered_map<uint32_t, uint32_t> estFlowSpreads;
        sketch.getOnlineEstimatedFlowSpreads(estFlowSpreads);

        vector<vector<double>> metrics=calculateMetrics(estFlowSpreads,actualFlowSpreads,topKs,outputDirPath,savedFileName+"-avg",traceFilePath);

        string metricsFilePath = outputDirPath+savedFileName + "-avg.metrics";
        ofstream metricsFile;
        metricsFile.open(metricsFilePath,ios::app);
        metricsFile<<"throughput: " << throughput << " Mpps" << ", each insert operation uses " << 1000.0 / throughput << " ns" << endl;
        metricsFile.close();

        avgThroughput+=throughput;
        for(uint32_t metricIdx=0;metricIdx<metrics.size();metricIdx++){
            for(uint32_t kIdx=0;kIdx<topKs.size();kIdx++){
                avgMetrics[metricIdx][kIdx]+=metrics[metricIdx][kIdx];
            }
        }
    }
    string metricsFilePath = outputDirPath+savedFileName + "-avg.metrics";
    ofstream metricsFile;
    metricsFile.open(metricsFilePath,ios::app);

    metricsFile<<"******************"<<endl;
    metricsFile<<"Average Super Spreader Metrics"<<endl;
    metricsFile<<"******************"<<endl;
    //trueFlowNum means the true number of super spreaders with that threshold
    //fakeFlowNum means the number of fake flows that returned by reversible calculation operation based on Chinese Remainder Theorem, in other words, these flows do not exist in the packet stream.
    metricsFile<<"Threshold\t\tPR\tRC\tF1\tARE\tAAE\tTP\tFP\tFN\ttrueFlowNum\t\tfakeFlowNum"<<endl;

    //THs,PRs,RCs,F1s,AREs,AAEs,TPs,FPs,FNs,realFlowNums,fakeFlowNums};
    for(uint32_t kIdx=0;kIdx<topKs.size();kIdx++){
        for(uint32_t metricIdx=0;metricIdx<avgMetrics.size();metricIdx++){

            avgMetrics[metricIdx][kIdx]=avgMetrics[metricIdx][kIdx]/fileNum;
            char temp[500]{0};
            sprintf(temp,"%-10.3f\t",avgMetrics[metricIdx][kIdx]);
            metricsFile<<temp;
        }
        metricsFile<<endl;
    }
    avgThroughput=avgThroughput/fileNum;
    metricsFile<<"Average throughput: " << avgThroughput << " Mpps" << endl;
    metricsFile.close();

    cout << "Average throughput: " << avgThroughput << " Mpps" << endl;
}


template<uint32_t BUCKET_NUM,uint32_t BUCKET_SIZE,uint32_t REG_GROUP_NUM,uint32_t REG_GROUP_SIZE,uint32_t V_REG_GROUP_NUM>
void processMultipleTracesAndGetMergedMetrics(const char* traceFileDir,uint32_t fileNum,const vector<uint32_t>& topKs, string outputDirPath,string traceInfo){
    uint32_t totalMem=BUCKET_NUM*BUCKET_SIZE*64+REG_GROUP_NUM*5*REG_GROUP_SIZE;
    double bucketMemRatio=(double)BUCKET_NUM*BUCKET_SIZE*64/totalMem;
    cout << "totalMem:" << totalMem/(8*1024.0) << "KB" << endl;
    cout << "*********************" << endl;

    char temp[500]{0};
    sprintf(temp,"%s data %d files src RGS-Sketch v_GNum=%d rGNum=%d rGSz=%d bktNum=%d bktSz=%d bktR=%.2f mem=%.2fKB",traceInfo.c_str(),fileNum,V_REG_GROUP_NUM,REG_GROUP_NUM,REG_GROUP_SIZE,BUCKET_NUM,BUCKET_SIZE,bucketMemRatio,totalMem/(8.0*1024));
    string savedFileName(temp);

    vector<SKETCH<BUCKET_NUM,BUCKET_SIZE,REG_GROUP_NUM,REG_GROUP_SIZE,V_REG_GROUP_NUM>> sketches;
    unordered_map<uint32_t, unordered_set<uint32_t>> totalFlowElements;

    for(uint32_t i=0;i<fileNum;i++){
        char traceFilePath[256]{0};
        sprintf(traceFilePath,"%s%02d.dat",traceFileDir,i);

        //prepare dataset
        cout << "prepare dataset: " <<traceFilePath<< endl;
        vector<pair<uint32_t,uint32_t>> pkts;
        unordered_map<uint32_t, uint32_t> actualFlowSpreads;
        unsigned int actualItemNum = ReadInSingleTraceFile(traceFilePath, pkts, actualFlowSpreads,totalFlowElements);

        cout << "number of packets: " << actualItemNum << endl;
        cout << "number of flows: " << actualFlowSpreads.size() << endl;
        cout << "*********************" << endl;

        cout << "prepare algorithm "<< endl;
        SKETCH<BUCKET_NUM,BUCKET_SIZE,REG_GROUP_NUM,REG_GROUP_SIZE,V_REG_GROUP_NUM> sketch;

        for (unsigned int i = 0; i < pkts.size(); i++) {
            sketch.insert(pkts[i]);
        }
        cout << "*********************" << endl;

//        unordered_map<uint32_t, uint32_t> estFlowSpreads;
//        sketch.getVHLLEstimatedFlowSpreads(estFlowSpreads);
//        calculateMetrics(estFlowSpreads,actualFlowSpreads,topKs,outputDirPath,savedFileName+"-merge",traceFilePath);

        sketches.push_back(sketch);


        if(i>=1){
            unordered_map<uint32_t, uint32_t> estFlowSpreads;
            SKETCH<BUCKET_NUM,BUCKET_SIZE,REG_GROUP_NUM,REG_GROUP_SIZE,V_REG_GROUP_NUM>::getMergedEstimatedFlowSpreads(sketches,estFlowSpreads);

            unordered_map<uint32_t, uint32_t> totalActualFlowSpreads;
            for(auto iter=totalFlowElements.begin();iter!=totalFlowElements.end();iter++){
                uint32_t key=iter->first;
                uint32_t acutalSpread=iter->second.size();
                totalActualFlowSpreads[key]=acutalSpread;
            }

            char temp[500]{0};
            sprintf(temp,"%s data %d files src RGS-Sketch v_GNum=%d rGNum=%d rGSz=%d bktNum=%d bktSz=%d bktR=%.2f mem=%.2fKB",traceInfo.c_str(),i+1,V_REG_GROUP_NUM,REG_GROUP_NUM,REG_GROUP_SIZE,BUCKET_NUM,BUCKET_SIZE,bucketMemRatio,totalMem/(8.0*1024));
            string tempSavedFileName(temp);

            calculateMetrics(estFlowSpreads,totalActualFlowSpreads,topKs,outputDirPath,tempSavedFileName+"-merge","Merge");
        }
    }
}



int main() {

    vector<uint32_t> topKs={50,100,150,200,250,300,350,400,450,500};
    string outputDirPath="./results/";

    const uint32_t bucketSize=8;
    const uint32_t regGroupSize=8;
    const uint32_t vRegGroupNum=32;

//     params for single measurement point
    double bucketMemRatio=0.35;
    //100kB
    const uint32_t bucketNum_100KB=560;
    const uint32_t regGroupNum_100KB=13312;
    //200kB
    const uint32_t bucketNum_200KB=1120;
    const uint32_t regGroupNum_200KB=26624;
    //300kB
    const uint32_t bucketNum_300KB=1680;
    const uint32_t regGroupNum_300KB=39936;
    //400kB
    const uint32_t bucketNum_400KB=2240;
    const uint32_t regGroupNum_400KB=53248;
    //500kB
    const uint32_t bucketNum_500KB=2800;
    const uint32_t regGroupNum_500KB=66560;


 //    // params for multiple measurement points
    //double bucketMemRatio=0.1;
    //500kB
    const uint32_t bucketNum=800;
    const uint32_t regGroupNum=92160;


    uint32_t fileNum=10;
    const char* traceFileDir2016=R"(../data/2016/)";
    const char* traceFileDir2019=R"(../data/2019/)";
    const char* traceFileDirCos=R"(../data/Cos/)";
    const char* traceFileDirMul=R"(../data/Mul2020.4/)";

// //    single measurement point
    processMultipleTracesAndGetAvgMetrics<bucketNum_100KB,bucketSize,regGroupNum_100KB,regGroupSize,vRegGroupNum> (traceFileDir2016,fileNum,topKs,outputDirPath,"2016");
//    processMultipleTracesAndGetAvgMetrics<bucketNum_200KB,bucketSize,regGroupNum_200KB,regGroupSize,vRegGroupNum> (traceFileDir2016,fileNum,topKs,outputDirPath,"2016");
//    processMultipleTracesAndGetAvgMetrics<bucketNum_300KB,bucketSize,regGroupNum_300KB,regGroupSize,vRegGroupNum> (traceFileDir2016,fileNum,topKs,outputDirPath,"2016");
//    processMultipleTracesAndGetAvgMetrics<bucketNum_400KB,bucketSize,regGroupNum_400KB,regGroupSize,vRegGroupNum> (traceFileDir2016,fileNum,topKs,outputDirPath,"2016");
//    processMultipleTracesAndGetAvgMetrics<bucketNum_500KB,bucketSize,regGroupNum_500KB,regGroupSize,vRegGroupNum> (traceFileDir2016,fileNum,topKs,outputDirPath,"2016");

    processMultipleTracesAndGetAvgMetrics<bucketNum_100KB,bucketSize,regGroupNum_100KB,regGroupSize,vRegGroupNum> (traceFileDir2019,fileNum,topKs,outputDirPath,"2019");
//    processMultipleTracesAndGetAvgMetrics<bucketNum_200KB,bucketSize,regGroupNum_200KB,regGroupSize,vRegGroupNum> (traceFileDir2019,fileNum,topKs,outputDirPath,"2019");
//    processMultipleTracesAndGetAvgMetrics<bucketNum_300KB,bucketSize,regGroupNum_300KB,regGroupSize,vRegGroupNum> (traceFileDir2019,fileNum,topKs,outputDirPath,"2019");
//    processMultipleTracesAndGetAvgMetrics<bucketNum_400KB,bucketSize,regGroupNum_400KB,regGroupSize,vRegGroupNum> (traceFileDir2019,fileNum,topKs,outputDirPath,"2019");
//    processMultipleTracesAndGetAvgMetrics<bucketNum_500KB,bucketSize,regGroupNum_500KB,regGroupSize,vRegGroupNum> (traceFileDir2019,fileNum,topKs,outputDirPath,"2019");

    processMultipleTracesAndGetAvgMetrics<bucketNum_100KB,bucketSize,regGroupNum_100KB,regGroupSize,vRegGroupNum> (traceFileDirCos,fileNum,topKs,outputDirPath,"Cos");
//    processMultipleTracesAndGetAvgMetrics<bucketNum_200KB,bucketSize,regGroupNum_200KB,regGroupSize,vRegGroupNum> (traceFileDirCos,fileNum,topKs,outputDirPath,"Cos");
//    processMultipleTracesAndGetAvgMetrics<bucketNum_300KB,bucketSize,regGroupNum_300KB,regGroupSize,vRegGroupNum> (traceFileDirCos,fileNum,topKs,outputDirPath,"Cos");
//    processMultipleTracesAndGetAvgMetrics<bucketNum_400KB,bucketSize,regGroupNum_400KB,regGroupSize,vRegGroupNum> (traceFileDirCos,fileNum,topKs,outputDirPath,"Cos");
//    processMultipleTracesAndGetAvgMetrics<bucketNum_500KB,bucketSize,regGroupNum_500KB,regGroupSize,vRegGroupNum> (traceFileDirCos,fileNum,topKs,outputDirPath,"Cos");

    processMultipleTracesAndGetAvgMetrics<bucketNum_100KB,bucketSize,regGroupNum_100KB,regGroupSize,vRegGroupNum> (traceFileDirMul,fileNum,topKs,outputDirPath,"Mul 2020.4");
//    processMultipleTracesAndGetAvgMetrics<bucketNum_200KB,bucketSize,regGroupNum_200KB,regGroupSize,vRegGroupNum> (traceFileDirMul,fileNum,topKs,outputDirPath,"Mul 2020.4");
//    processMultipleTracesAndGetAvgMetrics<bucketNum_300KB,bucketSize,regGroupNum_300KB,regGroupSize,vRegGroupNum> (traceFileDirMul,fileNum,topKs,outputDirPath,"Mul 2020.4");
//    processMultipleTracesAndGetAvgMetrics<bucketNum_400KB,bucketSize,regGroupNum_400KB,regGroupSize,vRegGroupNum> (traceFileDirMul,fileNum,topKs,outputDirPath,"Mul 2020.4");
//    processMultipleTracesAndGetAvgMetrics<bucketNum_500KB,bucketSize,regGroupNum_500KB,regGroupSize,vRegGroupNum> (traceFileDirMul,fileNum,topKs,outputDirPath,"Mul 2020.4");

 //     multiple measurement points
    processMultipleTracesAndGetMergedMetrics<bucketNum,bucketSize,regGroupNum,regGroupSize,vRegGroupNum>(traceFileDir2016,fileNum,topKs,outputDirPath,"2016");
//    processMultipleTracesAndGetMergedMetrics<bucketNum,bucketSize,regGroupNum,regGroupSize,vRegGroupNum>(traceFileDir2019,fileNum,topKs,outputDirPath,"2019");
//    processMultipleTracesAndGetMergedMetrics<bucketNum,bucketSize,regGroupNum,regGroupSize,vRegGroupNum>(traceFileDirCos,fileNum,topKs,outputDirPath,"Cos");
//    processMultipleTracesAndGetMergedMetrics<bucketNum,bucketSize,regGroupNum,regGroupSize,vRegGroupNum>(traceFileDirMul,fileNum,topKs,outputDirPath,"Mul 2020.4");

    return 0;
}
