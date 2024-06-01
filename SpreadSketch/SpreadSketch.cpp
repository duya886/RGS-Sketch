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
#include <bitset>

#define HASH_SEED 51318471
using namespace std;


const uint32_t HASH_SEEDS[]={674876105,860938402,3044024791,706198233,1838710314,2602646517};

template<uint32_t BIT_NUM1,uint32_t BIT_NUM2,uint32_t SET_MAX,uint32_t SET_LAST_MAX,uint32_t BITMAP_NUM>
class MultiResolutionBitmap{
private:
    vector<bitset<BIT_NUM1>> baseBitmaps;
    bitset<BIT_NUM2> lastBitmap;
public:
    MultiResolutionBitmap():baseBitmaps(BITMAP_NUM-1){};

    void update(const pair<uint32_t,uint32_t>& pkt,uint32_t mszNum){
        uint32_t hashValue2=0;
        MurmurHash3_x86_32(&pkt, 8,HASH_SEED+15417891, &hashValue2);
        if(mszNum<BITMAP_NUM-1){
            uint32_t bitIdx=hashValue2%BIT_NUM1;
            baseBitmaps[mszNum].set(bitIdx);
        }else{
            uint32_t bitIdx=hashValue2%BIT_NUM2;
            lastBitmap.set(bitIdx);
        }
    }
    uint32_t query() const{
        int32_t tempIdx=BITMAP_NUM-2;
        while(tempIdx>=0 && baseBitmaps[tempIdx].count()<=SET_MAX){
            tempIdx--;
        }
        uint32_t baseIdx=tempIdx+1;

        if(baseIdx==BITMAP_NUM-1 && lastBitmap.count()>SET_LAST_MAX){
            if(lastBitmap.count()>=BIT_NUM2){
                //the last bitmap is full, and cannot get estimated value by linear counting
                //return the estimation when there is 1 zero-bit
                double baseEst=BIT_NUM2*log(BIT_NUM2);
                double factor=pow(2,baseIdx);
                return round(baseEst*factor);
            }else{
                //not accurate
                double baseEst=-1.0*BIT_NUM2*log(((double)BIT_NUM2-lastBitmap.count())/BIT_NUM2);
                double factor=pow(2,baseIdx);
                return round(baseEst*factor);
            }
        }
        double baseEst=0;
        for(uint32_t idx=baseIdx;idx<BITMAP_NUM-1;idx++){
            baseEst-=BIT_NUM1*log(((double)BIT_NUM1-baseBitmaps[idx].count())/BIT_NUM1);
        }
        baseEst-=BIT_NUM2*log(((double)BIT_NUM2-lastBitmap.count())/BIT_NUM2);
        double factor=pow(2,baseIdx);
        return round(baseEst*factor);
    }
    void intersectOp(const MultiResolutionBitmap& MRB){
        for(uint32_t idx=0;idx<BITMAP_NUM-1;idx++){
            for(uint32_t bitIdx=0;bitIdx<BIT_NUM1;bitIdx++){
                if(!MRB.baseBitmaps[idx].test(bitIdx)){
                    baseBitmaps[idx].reset(bitIdx);
                }
            }
        }
        for(uint32_t bitIdx=0;bitIdx<BIT_NUM2;bitIdx++){
            if(!MRB.lastBitmap.test(bitIdx)){
                lastBitmap.reset(bitIdx);
            }
        }
    }
    void unionOp(const MultiResolutionBitmap& MRB){
        for(uint32_t idx=0;idx<BITMAP_NUM-1;idx++){
            for(uint32_t bitIdx=0;bitIdx<BIT_NUM1;bitIdx++){
                if(MRB.baseBitmaps[idx].test(bitIdx)){
                    baseBitmaps[idx].set(bitIdx);
                }
            }
        }
        for(uint32_t bitIdx=0;bitIdx<BIT_NUM2;bitIdx++){
            if(MRB.lastBitmap.test(bitIdx)){
                lastBitmap.set(bitIdx);
            }
        }
    }
};

template<uint32_t BIT_NUM1,uint32_t BIT_NUM2,uint32_t SET_MAX,uint32_t SET_LAST_MAX,uint32_t BITMAP_NUM>
struct BUCKET{
    MultiResolutionBitmap<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM> MRB;
    uint32_t key;
    uint8_t level;
};

template<uint32_t BIT_NUM1,uint32_t BIT_NUM2,uint32_t SET_MAX,uint32_t SET_LAST_MAX,uint32_t BITMAP_NUM,uint32_t ROW_NUM,uint32_t COL_NUM>
class SKETCH{
private:
    vector<vector<BUCKET<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM>>> buckets;

public:
    SKETCH():buckets(ROW_NUM,vector<BUCKET<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM>>(COL_NUM)){}
    void insert(const pair<uint32_t,uint32_t>& pkt) {
        uint32_t hashValue1 = 0;
        MurmurHash3_x86_32(&pkt, 8, HASH_SEED, &hashValue1);
        uint32_t mszNum = __builtin_clz(hashValue1);

        for (uint32_t rowIdx = 0; rowIdx < ROW_NUM; rowIdx++) {
            uint32_t hashValue = 0;
            MurmurHash3_x86_32(&pkt.first, 4, HASH_SEEDS[rowIdx], &hashValue);
            uint32_t colIdx = hashValue % COL_NUM;
            buckets[rowIdx][colIdx].MRB.update(pkt,mszNum);

            if (buckets[rowIdx][colIdx].level <= mszNum) {
                buckets[rowIdx][colIdx].key = pkt.first;
                buckets[rowIdx][colIdx].level = mszNum;
            }
        }
    }
    uint32_t query(uint32_t key)const{
        uint32_t hashValue2 = 0;
        MurmurHash3_x86_32(&key,4,HASH_SEEDS[0], &hashValue2);
        uint32_t colIdx = hashValue2 % COL_NUM;

        MultiResolutionBitmap<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM> newMRB=buckets[0][colIdx].MRB;

        for (uint32_t rowIdx = 1; rowIdx < ROW_NUM; rowIdx++) {
            MurmurHash3_x86_32(&key,4,HASH_SEEDS[rowIdx], &hashValue2);
            colIdx = hashValue2 % COL_NUM;
            newMRB.intersectOp(buckets[rowIdx][colIdx].MRB);
        }
        return newMRB.query();
    }

    void getEstimatedFlowSpreads(unordered_map<uint32_t, uint32_t>& estimatedFlowSpreads) const{
        for (uint32_t rowIdx = 0; rowIdx < ROW_NUM; rowIdx++) {
            for (uint32_t colIdx = 0; colIdx < COL_NUM; colIdx++) {
                uint32_t key=buckets[rowIdx][colIdx].key;

                auto iter2=estimatedFlowSpreads.find(key);
                if(iter2==estimatedFlowSpreads.end()){
                    estimatedFlowSpreads[key]= query(key);
                }
            }
        }
    }
    static void getMergedEstimatedFlowSpreads(const vector<SKETCH>& sketches,unordered_map<uint32_t, uint32_t>& estimatedFlowSpreads){
        SKETCH<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM,ROW_NUM,COL_NUM> mergedSketch;

        for(uint32_t sketchIdx=0;sketchIdx<sketches.size();sketchIdx++){
            for (uint32_t rowIdx = 0; rowIdx < ROW_NUM; rowIdx++) {
                for (uint32_t colIdx = 0; colIdx < COL_NUM; colIdx++) {
                    mergedSketch.buckets[rowIdx][colIdx].MRB.unionOp(sketches[sketchIdx].buckets[rowIdx][colIdx].MRB);
                    if(mergedSketch.buckets[rowIdx][colIdx].level<sketches[sketchIdx].buckets[rowIdx][colIdx].level || sketchIdx==0){
                        mergedSketch.buckets[rowIdx][colIdx].level=sketches[sketchIdx].buckets[rowIdx][colIdx].level;
                        mergedSketch.buckets[rowIdx][colIdx].key=sketches[sketchIdx].buckets[rowIdx][colIdx].key;
                    }
                }
            }
        }

        mergedSketch.getEstimatedFlowSpreads(estimatedFlowSpreads);
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

template<uint32_t BIT_NUM1,uint32_t BIT_NUM2,uint32_t SET_MAX,uint32_t SET_LAST_MAX,uint32_t BITMAP_NUM,uint32_t ROW_NUM,uint32_t COL_NUM>
void processMultipleTracesAndGetAvgMetrics(const char* traceFileDir,uint32_t fileNum,const vector<uint32_t>& topKs, string outputDirPath,string traceInfo){
    uint32_t totalMem=((BIT_NUM1*(BITMAP_NUM-1)+BIT_NUM2)+32+5)*ROW_NUM*COL_NUM;
    cout << "totalMem:" << totalMem/(8*1024.0) << "KB" << endl;
    cout << "*********************" << endl;

    char temp[500]{0};
    sprintf(temp,"%s data %d files src SS b=%d bL=%d sM=%d sML=%d c=%d r=%d w=%d mem=%.2fKB",traceInfo.c_str(),fileNum,BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM,ROW_NUM,COL_NUM,totalMem/(8.0*1024));
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
        SKETCH<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM,ROW_NUM,COL_NUM> sketch;

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
        sketch.getEstimatedFlowSpreads(estFlowSpreads);

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


template<uint32_t BIT_NUM1,uint32_t BIT_NUM2,uint32_t SET_MAX,uint32_t SET_LAST_MAX,uint32_t BITMAP_NUM,uint32_t ROW_NUM,uint32_t COL_NUM>
void processMultipleTracesAndGetMergedMetrics(const char* traceFileDir,uint32_t fileNum,const vector<uint32_t>& topKs, string outputDirPath,string traceInfo){
    uint32_t totalMem=((BIT_NUM1*(BITMAP_NUM-1)+BIT_NUM2)+32+5)*ROW_NUM*COL_NUM;
    cout << "totalMem:" << totalMem/(8*1024.0) << "KB" << endl;
    cout << "*********************" << endl;

    char temp[500]{0};
    sprintf(temp,"%s data %d files src SS b=%d bL=%d sM=%d sML=%d c=%d r=%d w=%d mem=%.2fKB",traceInfo.c_str(),fileNum,BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM,ROW_NUM,COL_NUM,totalMem/(8.0*1024));
    string savedFileName(temp);

    vector<SKETCH<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM,ROW_NUM,COL_NUM>> sketches;
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
        SKETCH<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM,ROW_NUM,COL_NUM> sketch;

        for (unsigned int i = 0; i < pkts.size(); i++) {
            sketch.insert(pkts[i]);
        }
        cout << "*********************" << endl;

//        unordered_map<uint32_t, uint32_t> estFlowSpreads;
//        sketch.getEstimatedFlowSpreads(estFlowSpreads);
//        calculateMetrics(estFlowSpreads,actualFlowSpreads,topKs,outputDirPath,savedFileName+"-merge",traceFilePath);

        sketches.push_back(sketch);


        if(i==fileNum-1){
            unordered_map<uint32_t, uint32_t> estFlowSpreads;
            SKETCH<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM,ROW_NUM,COL_NUM>::getMergedEstimatedFlowSpreads(sketches,estFlowSpreads);

            unordered_map<uint32_t, uint32_t> totalActualFlowSpreads;
            for(auto iter=totalFlowElements.begin();iter!=totalFlowElements.end();iter++){
                uint32_t key=iter->first;
                uint32_t acutalSpread=iter->second.size();
                totalActualFlowSpreads[key]=acutalSpread;
            }

            char temp[500]{0};
            sprintf(temp,"%s data %d files src SS b=%d bL=%d sM=%d sML=%d c=%d r=%d w=%d mem=%.2fKB",traceInfo.c_str(),i+1,BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM,ROW_NUM,COL_NUM,totalMem/(8.0*1024));
            string tempSavedFileName(temp);

            calculateMetrics(estFlowSpreads,totalActualFlowSpreads,topKs,outputDirPath,tempSavedFileName+"-merge","Merge");
        }
    }
}


int main() {

    vector<uint32_t> topKs={50,100,150,200,250,300,350,400,450,500};
    string outputDirPath="./results/";

    const uint32_t BIT_NUM1=64;
    const uint32_t BIT_NUM2=176;
    const uint32_t SET_MAX=60;
    const uint32_t SET_LAST_MAX=162;
    const uint32_t ROW_NUM=4;

    // params for 2016 single measurement point, max flow cardinality 12524
    // 100KB
    const uint32_t BITMAP_NUM_2016_100KB=6;
    const uint32_t COL_NUM_2016_100KB=384;
    // 200KB
    const uint32_t BITMAP_NUM_2016_200KB=6;
    const uint32_t COL_NUM_2016_200KB=768;
    // 300KB
    const uint32_t BITMAP_NUM_2016_300KB=6;
    const uint32_t COL_NUM_2016_300KB=1152;
    // 400KB
    const uint32_t BITMAP_NUM_2016_400KB=6;
    const uint32_t COL_NUM_2016_400KB=1536;
    // 500KB
    const uint32_t BITMAP_NUM_2016_500KB=6;
    const uint32_t COL_NUM_2016_500KB=1921;

    // params for 2019 single measurement point, max flow cardinality 96815
    // 100KB
    const uint32_t BITMAP_NUM_2019_100KB=9;
    const uint32_t COL_NUM_2019_100KB=282;
    // 200KB
    const uint32_t BITMAP_NUM_2019_200KB=9;
    const uint32_t COL_NUM_2019_200KB=564;
    // 300KB
    const uint32_t BITMAP_NUM_2019_300KB=9;
    const uint32_t COL_NUM_2019_300KB=847;
    // 400KB
    const uint32_t BITMAP_NUM_2019_400KB=9;
    const uint32_t COL_NUM_2019_400KB=1129;
    // 500KB
    const uint32_t BITMAP_NUM_2019_500KB=9;
    const uint32_t COL_NUM_2019_500KB=1412;


    // params for Cos single measurement point, max flow cardinality 11350
    // 100KB
    const uint32_t BITMAP_NUM_Cos_100KB=6;
    const uint32_t COL_NUM_Cos_100KB=384;
    // 200KB
    const uint32_t BITMAP_NUM_Cos_200KB=6;
    const uint32_t COL_NUM_Cos_200KB=768;
    // 300KB
    const uint32_t BITMAP_NUM_Cos_300KB=6;
    const uint32_t COL_NUM_Cos_300KB=1152;
    // 400KB
    const uint32_t BITMAP_NUM_Cos_400KB=6;
    const uint32_t COL_NUM_Cos_400KB=1536;
    // 500KB
    const uint32_t BITMAP_NUM_Cos_500KB=6;
    const uint32_t COL_NUM_Cos_500KB=1921;


    //params for Mul 2020.4 single measurement point, max flow cardinality 38570
    //    624bits Bitmap
    //    2644 bis for 4 bucket
    // 100KB
    const uint32_t BITMAP_NUM_Mul_100KB=8;
    const uint32_t COL_NUM_Mul_100KB=309;
    // 200KB
    const uint32_t BITMAP_NUM_Mul_200KB=8;
    const uint32_t COL_NUM_Mul_200KB=619;
    // 300KB
    const uint32_t BITMAP_NUM_Mul_300KB=8;
    const uint32_t COL_NUM_Mul_300KB=929;
    // 400KB
    const uint32_t BITMAP_NUM_Mul_400KB=8;
    const uint32_t COL_NUM_Mul_400KB=1239;
    // 500KB
    const uint32_t BITMAP_NUM_Mul_500KB=8;
    const uint32_t COL_NUM_Mul_500KB=1549;


    // 500KB params for 2016 multiple measurement points
    //max flow cardinality: 2:24203; 3:36107; 4:47826; 5:59361; 6:71777; 7:82874; 8:95056; 9:105997; 10: 113494
    // 2 points
    const uint32_t BITMAP_NUM_2016M_2=7;
    const uint32_t COL_NUM_2016M_2=1715;
    // 3-4 points
    const uint32_t BITMAP_NUM_2016M_3_4=8;
    const uint32_t COL_NUM_2016M_3_4=1549;
    // 5-8 points
    const uint32_t BITMAP_NUM_2016M_5_8=9;
    const uint32_t COL_NUM_2016M_5_8=1412;
    // 9-10 points
    const uint32_t BITMAP_NUM_2016M_9_10=10;
    const uint32_t COL_NUM_2016M_9_10=1297;

    // 500KB params for 2019 multiple measurement points
    //max flow cardinality: 2:183447; 3:265670; 4:340986; 5:412381; 6:476491; 7:536774; 8:591833; 9:642557; 10: 690766
    // 2 points
    const uint32_t BITMAP_NUM_2019M_2=10;
    const uint32_t COL_NUM_2019M_2=1297;
    // 3-5 points
    const uint32_t BITMAP_NUM_2019M_3_5=11;
    const uint32_t COL_NUM_2019M_3_5=1200;
    // 6-10 points
    const uint32_t BITMAP_NUM_2019M_6_10=12;
    const uint32_t COL_NUM_2019M_6_10=1116;

    // 500KB params for Cos multiple measurement points
    //max flow cardinality: 2:5831; 3:11236; 4:18125; 5:26205; 6:35955; 7:43746; 8:51115; 9:56971; 10: 62461
    //2-3 calculated using 11350 (single point max cardinality)
    // 2-3 points
    const uint32_t BITMAP_NUM_CosM_2_3=6;
    const uint32_t COL_NUM_CosM_2_3=1921;
    // 4-5 points
    const uint32_t BITMAP_NUM_CosM_4_5=7;
    const uint32_t COL_NUM_CosM_4_5=1715;
    // 6-8 points
    const uint32_t BITMAP_NUM_CosM_6_8=8;
    const uint32_t COL_NUM_CosM_6_8=1549;
    // 9-10 points
    const uint32_t BITMAP_NUM_CosM_9_10=9;
    const uint32_t COL_NUM_CosM_9_10=1412;

    // 500KB params for Mul2020.4 multiple measurement points
    //max flow cardinality: 2:45486; 3:65339; 4:83443; 5:101751; 6:118450; 7:135596; 8:152409; 9:170053; 10:188666
    // 2 points
    const uint32_t BITMAP_NUM_MulM_2=8;
    const uint32_t COL_NUM_MulM_2=1549;
    // 3-5 points
    const uint32_t BITMAP_NUM_MulM_3_5=9;
    const uint32_t COL_NUM_MulM_3_5=1412;
    // 6-10 points
    const uint32_t BITMAP_NUM_MulM_6_10=10;
    const uint32_t COL_NUM_MulM_6_10=1297;

    const char* traceFileDir2016=R"(../data/2016/)";
    const char* traceFileDir2019=R"(../data/2019/)";
    const char* traceFileDirCos=R"(../data/Cos/)";
    const char* traceFileDirMul=R"(../data/Mul2020.4/)";

    uint32_t fileNum=10;
    
    //single measurement point
//    processMultipleTracesAndGetAvgMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2016_100KB,ROW_NUM,COL_NUM_2016_100KB>(traceFileDir2016,fileNum,topKs,outputDirPath,"2016");
//    processMultipleTracesAndGetAvgMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2016_200KB,ROW_NUM,COL_NUM_2016_200KB>(traceFileDir2016,fileNum,topKs,outputDirPath,"2016");
//    processMultipleTracesAndGetAvgMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2016_300KB,ROW_NUM,COL_NUM_2016_300KB>(traceFileDir2016,fileNum,topKs,outputDirPath,"2016");
//    processMultipleTracesAndGetAvgMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2016_400KB,ROW_NUM,COL_NUM_2016_400KB>(traceFileDir2016,fileNum,topKs,outputDirPath,"2016");
    processMultipleTracesAndGetAvgMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2016_500KB,ROW_NUM,COL_NUM_2016_500KB>(traceFileDir2016,fileNum,topKs,outputDirPath,"2016");
//
//    processMultipleTracesAndGetAvgMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2019_100KB,ROW_NUM,COL_NUM_2019_100KB>(traceFileDir2019,fileNum,topKs,outputDirPath,"2019");
//    processMultipleTracesAndGetAvgMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2019_200KB,ROW_NUM,COL_NUM_2019_200KB>(traceFileDir2019,fileNum,topKs,outputDirPath,"2019");
//    processMultipleTracesAndGetAvgMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2019_300KB,ROW_NUM,COL_NUM_2019_300KB>(traceFileDir2019,fileNum,topKs,outputDirPath,"2019");
//    processMultipleTracesAndGetAvgMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2019_400KB,ROW_NUM,COL_NUM_2019_400KB>(traceFileDir2019,fileNum,topKs,outputDirPath,"2019");
    processMultipleTracesAndGetAvgMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2019_500KB,ROW_NUM,COL_NUM_2019_500KB>(traceFileDir2019,fileNum,topKs,outputDirPath,"2019");

//    processMultipleTracesAndGetAvgMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_Cos_100KB,ROW_NUM,COL_NUM_Cos_100KB>(traceFileDirCos,fileNum,topKs,outputDirPath,"Cos");
//    processMultipleTracesAndGetAvgMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_Cos_200KB,ROW_NUM,COL_NUM_Cos_200KB>(traceFileDirCos,fileNum,topKs,outputDirPath,"Cos");
//    processMultipleTracesAndGetAvgMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_Cos_300KB,ROW_NUM,COL_NUM_Cos_300KB>(traceFileDirCos,fileNum,topKs,outputDirPath,"Cos");
//    processMultipleTracesAndGetAvgMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_Cos_400KB,ROW_NUM,COL_NUM_Cos_400KB>(traceFileDirCos,fileNum,topKs,outputDirPath,"Cos");
    processMultipleTracesAndGetAvgMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_Cos_500KB,ROW_NUM,COL_NUM_Cos_500KB>(traceFileDirCos,fileNum,topKs,outputDirPath,"Cos");
//
//    processMultipleTracesAndGetAvgMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_Mul_100KB,ROW_NUM,COL_NUM_Mul_100KB>(traceFileDirMul,fileNum,topKs,outputDirPath,"Mul 2020.4");
//    processMultipleTracesAndGetAvgMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_Mul_200KB,ROW_NUM,COL_NUM_Mul_200KB>(traceFileDirMul,fileNum,topKs,outputDirPath,"Mul 2020.4");
//    processMultipleTracesAndGetAvgMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_Mul_300KB,ROW_NUM,COL_NUM_Mul_300KB>(traceFileDirMul,fileNum,topKs,outputDirPath,"Mul 2020.4");
//    processMultipleTracesAndGetAvgMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_Mul_400KB,ROW_NUM,COL_NUM_Mul_400KB>(traceFileDirMul,fileNum,topKs,outputDirPath,"Mul 2020.4");
    processMultipleTracesAndGetAvgMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_Mul_500KB,ROW_NUM,COL_NUM_Mul_500KB>(traceFileDirMul,fileNum,topKs,outputDirPath,"Mul 2020.4");



    //multiple measurement points
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2016M_2,ROW_NUM,COL_NUM_2016M_2>(traceFileDir2016,2,topKs,outputDirPath,"2016");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2016M_3_4,ROW_NUM,COL_NUM_2016M_3_4>(traceFileDir2016,3,topKs,outputDirPath,"2016");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2016M_3_4,ROW_NUM,COL_NUM_2016M_3_4>(traceFileDir2016,4,topKs,outputDirPath,"2016");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2016M_5_8,ROW_NUM,COL_NUM_2016M_5_8>(traceFileDir2016,5,topKs,outputDirPath,"2016");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2016M_5_8,ROW_NUM,COL_NUM_2016M_5_8>(traceFileDir2016,6,topKs,outputDirPath,"2016");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2016M_5_8,ROW_NUM,COL_NUM_2016M_5_8>(traceFileDir2016,7,topKs,outputDirPath,"2016");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2016M_5_8,ROW_NUM,COL_NUM_2016M_5_8>(traceFileDir2016,8,topKs,outputDirPath,"2016");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2016M_9_10,ROW_NUM,COL_NUM_2016M_9_10>(traceFileDir2016,9,topKs,outputDirPath,"2016");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2016M_9_10,ROW_NUM,COL_NUM_2016M_9_10>(traceFileDir2016,10,topKs,outputDirPath,"2016");
//
//
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2019M_2,ROW_NUM,COL_NUM_2019M_2>(traceFileDir2019,2,topKs,outputDirPath,"2019");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2019M_3_5,ROW_NUM,COL_NUM_2019M_3_5>(traceFileDir2019,3,topKs,outputDirPath,"2019");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2019M_3_5,ROW_NUM,COL_NUM_2019M_3_5>(traceFileDir2019,4,topKs,outputDirPath,"2019");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2019M_3_5,ROW_NUM,COL_NUM_2019M_3_5>(traceFileDir2019,5,topKs,outputDirPath,"2019");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2019M_6_10,ROW_NUM,COL_NUM_2019M_6_10>(traceFileDir2019,6,topKs,outputDirPath,"2019");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2019M_6_10,ROW_NUM,COL_NUM_2019M_6_10>(traceFileDir2019,7,topKs,outputDirPath,"2019");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2019M_6_10,ROW_NUM,COL_NUM_2019M_6_10>(traceFileDir2019,8,topKs,outputDirPath,"2019");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2019M_6_10,ROW_NUM,COL_NUM_2019M_6_10>(traceFileDir2019,9,topKs,outputDirPath,"2019");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2019M_6_10,ROW_NUM,COL_NUM_2019M_6_10>(traceFileDir2019,10,topKs,outputDirPath,"2019");
//
//
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_CosM_2_3,ROW_NUM,COL_NUM_CosM_2_3>(traceFileDirCos,2,topKs,outputDirPath,"Cos");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_CosM_2_3,ROW_NUM,COL_NUM_CosM_2_3>(traceFileDirCos,3,topKs,outputDirPath,"Cos");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_CosM_4_5,ROW_NUM,COL_NUM_CosM_4_5>(traceFileDirCos,4,topKs,outputDirPath,"Cos");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_CosM_4_5,ROW_NUM,COL_NUM_CosM_4_5>(traceFileDirCos,5,topKs,outputDirPath,"Cos");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_CosM_6_8,ROW_NUM,COL_NUM_CosM_6_8>(traceFileDirCos,6,topKs,outputDirPath,"Cos");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_CosM_6_8,ROW_NUM,COL_NUM_CosM_6_8>(traceFileDirCos,7,topKs,outputDirPath,"Cos");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_CosM_6_8,ROW_NUM,COL_NUM_CosM_6_8>(traceFileDirCos,8,topKs,outputDirPath,"Cos");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_CosM_9_10,ROW_NUM,COL_NUM_CosM_9_10>(traceFileDirCos,9,topKs,outputDirPath,"Cos");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_CosM_9_10,ROW_NUM,COL_NUM_CosM_9_10>(traceFileDirCos,10,topKs,outputDirPath,"Cos");

    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_MulM_2,ROW_NUM,COL_NUM_MulM_2>(traceFileDirMul,2,topKs,outputDirPath,"Mul 2020.4");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_MulM_3_5,ROW_NUM,COL_NUM_MulM_3_5>(traceFileDirMul,3,topKs,outputDirPath,"Mul 2020.4");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_MulM_3_5,ROW_NUM,COL_NUM_MulM_3_5>(traceFileDirMul,4,topKs,outputDirPath,"Mul 2020.4");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_MulM_3_5,ROW_NUM,COL_NUM_MulM_3_5>(traceFileDirMul,5,topKs,outputDirPath,"Mul 2020.4");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_MulM_6_10,ROW_NUM,COL_NUM_MulM_6_10>(traceFileDirMul,6,topKs,outputDirPath,"Mul 2020.4");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_MulM_6_10,ROW_NUM,COL_NUM_MulM_6_10>(traceFileDirMul,7,topKs,outputDirPath,"Mul 2020.4");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_MulM_6_10,ROW_NUM,COL_NUM_MulM_6_10>(traceFileDirMul,8,topKs,outputDirPath,"Mul 2020.4");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_MulM_6_10,ROW_NUM,COL_NUM_MulM_6_10>(traceFileDirMul,9,topKs,outputDirPath,"Mul 2020.4");
//    processMultipleTracesAndGetMergedMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_MulM_6_10,ROW_NUM,COL_NUM_MulM_6_10>(traceFileDirMul,10,topKs,outputDirPath,"Mul 2020.4");

    return 0;
}
