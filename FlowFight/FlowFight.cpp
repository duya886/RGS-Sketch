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

template<uint32_t REG_NUM>
class HyperLogLog{
private:
    uint8_t* regArray;
public:
    HyperLogLog(){
        uint32_t byteNum=int(ceil(REG_NUM*5.0/8.0));
        regArray=new uint8_t[byteNum]();
    }
    HyperLogLog(const HyperLogLog& other){
        uint32_t byteNum=int(ceil(REG_NUM*5.0/8.0));
        regArray=new uint8_t[byteNum]();

        for(uint32_t byteIdx=0;byteIdx<byteNum;byteIdx++){
            regArray[byteIdx]=other.regArray[byteIdx];
        }
    }
    HyperLogLog& operator=(const HyperLogLog& other){
        uint32_t byteNum=int(ceil(REG_NUM*5.0/8.0));

        for(uint32_t byteIdx=0;byteIdx<byteNum;byteIdx++){
            regArray[byteIdx]=other.regArray[byteIdx];
        }
        return *this;
    }
    uint32_t getRegValue(uint32_t regIdx) const{
        uint32_t regStartByteIdx=regIdx*5/8;
        uint32_t regStartBitIdx=regIdx*5%8;
        return ((*(uint16_t*)(regArray+regStartByteIdx))>>regStartBitIdx)&0x1F;
    }
    void setRegValue(uint32_t regIdx,uint32_t newValue){
        uint32_t regStartByteIdx=regIdx*5/8;
        uint32_t regStartBitIdx=regIdx*5%8;

        uint16_t* value_ptr=(uint16_t*)(regArray+regStartByteIdx);
        *value_ptr=(*value_ptr)&(~(((uint16_t)0x1F)<<regStartBitIdx));
        *value_ptr=(*value_ptr)|(newValue<<regStartBitIdx);
    }
    void insert(const pair<uint32_t,uint32_t>& pkt,uint32_t &oldVal, uint32_t &newVal){
        uint32_t hashValue1=0;
        MurmurHash3_x86_32(&pkt,8,HASH_SEED, &hashValue1);

        uint32_t hashValue2=0;
        MurmurHash3_x86_32(&pkt, 8,HASH_SEED+15417891, &hashValue2);

        uint32_t regIdx=hashValue1%REG_NUM;
        newVal=min(31,__builtin_clz(hashValue2)+1);

        oldVal= getRegValue(regIdx);
        if(newVal>oldVal){
            setRegValue(regIdx,newVal);
        }
    }

    void clear(){
        uint32_t byteNum=int(ceil(REG_NUM*5.0/8.0));

        for(uint32_t byteIdx=0;byteIdx<byteNum;byteIdx++){
            regArray[byteIdx]=0;
        }
    }

    uint32_t query() const{
        double coe= getCoe(REG_NUM);
        double sum=0.0;
        uint32_t numOfZeroRegs=0;


        for(uint32_t regIdx=0;regIdx<REG_NUM;regIdx++){
            uint32_t regValue= getRegValue(regIdx);
            sum+=1.0/pow(2,regValue);

            if(regValue==0){
                numOfZeroRegs++;
            }
        }

        double estimatedVal=coe*pow(REG_NUM,2)/sum;
        double spread=estimatedVal;

        if(estimatedVal<2.5*REG_NUM){
            if(numOfZeroRegs>0){
                spread=REG_NUM*log((double)REG_NUM/numOfZeroRegs);
            }
        }else if(estimatedVal>pow(2,32)/30.0){
            spread=-1.0*pow(2,32)*log(1.0-(double)estimatedVal/pow(2,32));
        }

        return round(spread);
    }
    ~HyperLogLog(){
        delete[] regArray;
    }

    static void getIdxVal(const pair<uint32_t,uint32_t>& pkt,uint32_t &regIdx,uint32_t &rhoVal){
        uint32_t hashValue1=0;
        MurmurHash3_x86_32(&pkt,8,HASH_SEED, &hashValue1);

        uint32_t hashValue2=0;
        MurmurHash3_x86_32(&pkt, 8,HASH_SEED+15417891, &hashValue2);

        regIdx=hashValue1%REG_NUM;
        rhoVal=min(31,__builtin_clz(hashValue2)+1);
    }
private:
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
};

class BUCKET{
public:
    //we record key instead of the reverseIdx, since we directly use unordered_map as the hash table
    //however, when compute memory, we suppose it is 16 bis, like a short index
    uint32_t reverseIdx=0;
    uint32_t rough=0;
    uint8_t regCnt=0;//the actual number of non-zero regs is regCnt+1
};

//the K here actually equals the sum of k and the size of intermediate zone
template<uint32_t K,uint32_t PZ_SIZE,uint32_t REG_NUM,uint32_t WEIGHT_TH,uint32_t PKT_MARGIN>
class SKETCH{
private:
    unordered_map<uint32_t,pair<uint32_t,uint32_t>> hashTable;
    vector<HyperLogLog<REG_NUM>> HLLs;
    vector<BUCKET> buckets;
    uint32_t pktCnt=0;
    uint32_t promotionIdx;//promotion zone is implemented as a circular queue, this idx points to the oldest pos of the promotion zone;

    inline void sortOnceFrom(uint32_t bktIdx){
        if(bktIdx>0 && buckets[bktIdx].rough>buckets[bktIdx-1].rough){
            BUCKET tempBkt=buckets[bktIdx-1];
            buckets[bktIdx-1]=buckets[bktIdx];
            buckets[bktIdx]=tempBkt;

            hashTable[buckets[bktIdx-1].reverseIdx].first=bktIdx-1;
            hashTable[buckets[bktIdx].reverseIdx].first=bktIdx;
        }
    }
public:
    SKETCH():HLLs(K+PZ_SIZE),buckets(K+PZ_SIZE){
        promotionIdx=K;
    }
    void insert(const pair<uint32_t,uint32_t>& pkt) {
        pktCnt++;

        uint32_t key=pkt.first;
        if(hashTable.find(key)!=hashTable.end()){
            pair<uint32_t,uint32_t> &idxes=hashTable[key];
            HyperLogLog<REG_NUM>& HLL=HLLs[idxes.second];
            uint32_t oldVal,newVal;
            HLL.insert(pkt,oldVal,newVal);

            if(newVal>oldVal){
                BUCKET& bkt=buckets[idxes.first];
                bkt.rough+=(1<<newVal)-(1<<oldVal);
                if(oldVal==0){
                    bkt.regCnt++;
                }

                if(idxes.first<K){
                    sortOnceFrom(idxes.first);
                }
            }
        }else if(hashTable.size()<K+PZ_SIZE){
            uint32_t bktIdx=hashTable.size();
            uint32_t HLLIdx=hashTable.size();
            hashTable[pkt.first]= make_pair(bktIdx,HLLIdx);

            uint32_t oldVal,newVal;
            HLLs[HLLIdx].insert(pkt,oldVal,newVal);
            buckets[bktIdx].reverseIdx=pkt.first;
            buckets[bktIdx].rough=1<<newVal;

            if(bktIdx<K){
                sortOnceFrom(bktIdx);
            }
        }else{
            uint32_t regIdx,rhoVal;
            HyperLogLog<REG_NUM>::getIdxVal(pkt,regIdx,rhoVal);

            //according to https://github.com/valebru/FlowFight-on-VPP/blob/master/vpp/src/plugins/hll/node.c line 243
            BUCKET& minBkt=buckets[K-1];
            uint32_t nonZeroRegNum=minBkt.regCnt+1;

            uint32_t newFlowNonZeroNum=1;
            if(nonZeroRegNum>WEIGHT_TH){
                newFlowNonZeroNum=1+(nonZeroRegNum>>5);
            }

            uint32_t newFlowFairRoughEst=REG_NUM+((1<<rhoVal)-1)*newFlowNonZeroNum;

            if(newFlowFairRoughEst>minBkt.rough || (nonZeroRegNum<WEIGHT_TH && pktCnt>PKT_MARGIN)){
                //kick old flow
                uint32_t kickedFlowId=minBkt.reverseIdx;
                uint32_t kickedHLLIdx=hashTable[kickedFlowId].second;

                HLLs[kickedHLLIdx].clear();

                auto iter=hashTable.find(kickedFlowId);
                hashTable.erase(iter);

                //move the oldest flow to the position of the kicked flow
                minBkt=buckets[promotionIdx];
                hashTable[buckets[promotionIdx].reverseIdx].first=K-1;

                //add new flow
                hashTable[pkt.first]= make_pair(promotionIdx,kickedHLLIdx);
                buckets[promotionIdx].reverseIdx=pkt.first;
                buckets[promotionIdx].rough=REG_NUM+1<<rhoVal;
                buckets[promotionIdx].regCnt=0;
                HLLs[kickedHLLIdx].setRegValue(regIdx,rhoVal);

                //update promotionIdx
                promotionIdx++;
                if(promotionIdx>=K+PZ_SIZE){
                    promotionIdx=K;
                }

                //update pkt count
                pktCnt=0;
            }

        }
    }
    uint32_t query(uint32_t key)const{
        auto iter=hashTable.find(key);
        if(iter!=hashTable.end()){
            uint32_t HLLIdx=iter->second.second;
            uint32_t val=HLLs[HLLIdx].query();
            return val;
        }else{
            return 0;
        }
    }

    void getEstimatedFlowSpreads(unordered_map<uint32_t, uint32_t>& estimatedFlowSpreads) const{
        for(auto iter=hashTable.begin();iter!=hashTable.end();iter++){
            uint32_t key=iter->first;
            uint32_t val= query(key);
            estimatedFlowSpreads[key]=val;
        }
    }
    void getAllFlowHLLs(unordered_map<uint32_t, HyperLogLog<REG_NUM>>& flowHLLs)const{
        for(auto iter=hashTable.begin();iter!=hashTable.end();iter++){
            uint32_t key=iter->first;
            flowHLLs[key]= HLLs[iter->second.second];
        }
    }

    static void getMergedEstimatedFlowSpreads(const vector<SKETCH>& sketches,unordered_map<uint32_t, uint32_t>& estimatedFlowSpreads){

        unordered_map<uint32_t,HyperLogLog<REG_NUM>> mergedFlowHLLs;
        for(uint32_t sketchIdx=0;sketchIdx<sketches.size();sketchIdx++){
            unordered_map<uint32_t, HyperLogLog<REG_NUM>> flowHLLs;
            sketches[sketchIdx].getAllFlowHLLs(flowHLLs);

            for(auto iter=flowHLLs.begin();iter!=flowHLLs.end();iter++){
                uint32_t key=iter->first;
                HyperLogLog<REG_NUM>& HLL=iter->second;

                if(mergedFlowHLLs.find(key)==mergedFlowHLLs.end()){
                    mergedFlowHLLs[key]=HLL;
                }else{
                    HyperLogLog<REG_NUM>& mergedHLL=mergedFlowHLLs[key];
                    for(uint32_t regIdx=0;regIdx<REG_NUM;regIdx++){
                        uint32_t rhoVal=HLL.getRegValue(regIdx);
                        uint32_t oldVal=mergedHLL.getRegValue(regIdx);
                        if(rhoVal>oldVal){
                            mergedHLL.setRegValue(regIdx,rhoVal);
                        }
                    }
                }
            }
        }

        for(auto iter=mergedFlowHLLs.begin();iter!=mergedFlowHLLs.end();iter++){
            uint32_t key=iter->first;
            uint32_t val= iter->second.query();

            estimatedFlowSpreads[key]=val;
        }
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

template<uint32_t K,uint32_t PZ_SIZE,uint32_t REG_NUM,uint32_t WEIGHT_TH,uint32_t PKT_MARGIN>
void processMultipleTracesAndGetAvgMetrics(const char* traceFileDir,uint32_t fileNum,const vector<uint32_t>& topKs, string outputDirPath,string traceInfo){
    uint32_t totalFlowNum=K+PZ_SIZE;
    uint32_t totalMem=totalFlowNum*(5*REG_NUM+(16+32+8)+(32+16+16));
    cout << "totalMem:" << totalMem/(8*1024.0) << "KB" << endl;
    cout << "*********************" << endl;

    char temp[500]{0};
    sprintf(temp,"%s data %d files src FF k=%d pZSize=%d regNum=%d w_th=%d pkt_mar=%d mem=%.2fKB",traceInfo.c_str(),fileNum,K,PZ_SIZE,REG_NUM,WEIGHT_TH,PKT_MARGIN,totalMem/(8.0*1024));
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
        SKETCH<K,PZ_SIZE,REG_NUM,WEIGHT_TH,PKT_MARGIN> sketch;

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


template<uint32_t K,uint32_t PZ_SIZE,uint32_t REG_NUM,uint32_t WEIGHT_TH,uint32_t PKT_MARGIN>
void processMultipleTracesAndGetMergedMetrics(const char* traceFileDir,uint32_t fileNum,const vector<uint32_t>& topKs, string outputDirPath,string traceInfo){
    uint32_t totalFlowNum=K+PZ_SIZE;
    uint32_t totalMem=totalFlowNum*(5*REG_NUM+(16+32+8)+(32+16+16));
    cout << "totalMem:" << totalMem/(8*1024.0) << "KB" << endl;
    cout << "*********************" << endl;

    char temp[500]{0};
    sprintf(temp,"%s data %d files src FF k=%d pZSize=%d regNum=%d w_th=%d pkt_mar=%d mem=%.2fKB",traceInfo.c_str(),fileNum,K,PZ_SIZE,REG_NUM,WEIGHT_TH,PKT_MARGIN,totalMem/(8.0*1024));
    string savedFileName(temp);

    vector<SKETCH<K,PZ_SIZE,REG_NUM,WEIGHT_TH,PKT_MARGIN>> sketches;
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
        SKETCH<K,PZ_SIZE,REG_NUM,WEIGHT_TH,PKT_MARGIN> sketch;

        for (unsigned int i = 0; i < pkts.size(); i++) {
            sketch.insert(pkts[i]);
        }
        cout << "*********************" << endl;

//        unordered_map<uint32_t, uint32_t> estFlowSpreads;
//        sketch.getEstimatedFlowSpreads(estFlowSpreads);
//        calculateMetrics(estFlowSpreads,actualFlowSpreads,topKs,outputDirPath,savedFileName+"-merge",traceFilePath);

        sketches.push_back(sketch);


        if(i>=1){
            unordered_map<uint32_t, uint32_t> estFlowSpreads;
            SKETCH<K,PZ_SIZE,REG_NUM,WEIGHT_TH,PKT_MARGIN>::getMergedEstimatedFlowSpreads(sketches,estFlowSpreads);

            unordered_map<uint32_t, uint32_t> totalActualFlowSpreads;
            for(auto iter=totalFlowElements.begin();iter!=totalFlowElements.end();iter++){
                uint32_t key=iter->first;
                uint32_t acutalSpread=iter->second.size();
                totalActualFlowSpreads[key]=acutalSpread;
            }

            char temp[500]{0};
            sprintf(temp,"%s data %d files src FF k=%d pZSize=%d regNum=%d w_th=%d pkt_mar=%d mem=%.2fKB",traceInfo.c_str(),i+1,K,PZ_SIZE,REG_NUM,WEIGHT_TH,PKT_MARGIN,totalMem/(8.0*1024));
            string tempSavedFileName(temp);

            calculateMetrics(estFlowSpreads,totalActualFlowSpreads,topKs,outputDirPath,tempSavedFileName+"-merge","Merge");
        }
    }
}


int main() {

    vector<uint32_t> topKs={50,100,150,200,250,300,350,400,450,500};
    string outputDirPath="./results/";

    const uint32_t REG_NUM=256;
    const uint32_t WEIGHT_TH=64;
    const uint32_t PKT_MARGIN=6400;

    // params for single measurement point
//    //100KB
    const uint32_t K_100KB=500+42;
    const uint32_t PZ_SIZE_100KB=43;
//    //200KB
    const uint32_t K_200KB=500+335;
    const uint32_t PZ_SIZE_200KB=335;
//    //200KB
    const uint32_t K_300KB=500+627;
    const uint32_t PZ_SIZE_300KB=628;
//    //400KB
    const uint32_t K_400KB=500+920;
    const uint32_t PZ_SIZE_400KB=920;
//    //500KB
    const uint32_t K_500KB=500+1212;
    const uint32_t PZ_SIZE_500KB=1213;
    
    
    uint32_t fileNum=10;
    const char* traceFileDir2016=R"(../data/2016/)";
    const char* traceFileDir2019=R"(../data/2019/)";
    const char* traceFileDirCos=R"(../data/Cos/)";
    const char* traceFileDirMul=R"(../data/Mul2020.4/)";


    //single measurement point
//    processMultipleTracesAndGetAvgMetrics<K_100KB,PZ_SIZE_100KB,REG_NUM,WEIGHT_TH,PKT_MARGIN>(traceFileDir2016,fileNum,topKs,outputDirPath,"2016");
//    processMultipleTracesAndGetAvgMetrics<K_200KB,PZ_SIZE_200KB,REG_NUM,WEIGHT_TH,PKT_MARGIN>(traceFileDir2016,fileNum,topKs,outputDirPath,"2016");
//    processMultipleTracesAndGetAvgMetrics<K_300KB,PZ_SIZE_300KB,REG_NUM,WEIGHT_TH,PKT_MARGIN>(traceFileDir2016,fileNum,topKs,outputDirPath,"2016");
//    processMultipleTracesAndGetAvgMetrics<K_400KB,PZ_SIZE_400KB,REG_NUM,WEIGHT_TH,PKT_MARGIN>(traceFileDir2016,fileNum,topKs,outputDirPath,"2016");
//    processMultipleTracesAndGetAvgMetrics<K_500KB,PZ_SIZE_500KB,REG_NUM,WEIGHT_TH,PKT_MARGIN>(traceFileDir2016,fileNum,topKs,outputDirPath,"2016");
//
//    processMultipleTracesAndGetAvgMetrics<K_100KB,PZ_SIZE_100KB,REG_NUM,WEIGHT_TH,PKT_MARGIN>(traceFileDir2019,fileNum,topKs,outputDirPath,"2019");
//    processMultipleTracesAndGetAvgMetrics<K_200KB,PZ_SIZE_200KB,REG_NUM,WEIGHT_TH,PKT_MARGIN>(traceFileDir2019,fileNum,topKs,outputDirPath,"2019");
//    processMultipleTracesAndGetAvgMetrics<K_300KB,PZ_SIZE_300KB,REG_NUM,WEIGHT_TH,PKT_MARGIN>(traceFileDir2019,fileNum,topKs,outputDirPath,"2019");
//    processMultipleTracesAndGetAvgMetrics<K_400KB,PZ_SIZE_400KB,REG_NUM,WEIGHT_TH,PKT_MARGIN>(traceFileDir2019,fileNum,topKs,outputDirPath,"2019");
//    processMultipleTracesAndGetAvgMetrics<K_500KB,PZ_SIZE_500KB,REG_NUM,WEIGHT_TH,PKT_MARGIN>(traceFileDir2019,fileNum,topKs,outputDirPath,"2019");
//
//    processMultipleTracesAndGetAvgMetrics<K_100KB,PZ_SIZE_100KB,REG_NUM,WEIGHT_TH,PKT_MARGIN>(traceFileDirCos,fileNum,topKs,outputDirPath,"Cos");
//    processMultipleTracesAndGetAvgMetrics<K_200KB,PZ_SIZE_200KB,REG_NUM,WEIGHT_TH,PKT_MARGIN>(traceFileDirCos,fileNum,topKs,outputDirPath,"Cos");
//    processMultipleTracesAndGetAvgMetrics<K_300KB,PZ_SIZE_300KB,REG_NUM,WEIGHT_TH,PKT_MARGIN>(traceFileDirCos,fileNum,topKs,outputDirPath,"Cos");
//    processMultipleTracesAndGetAvgMetrics<K_400KB,PZ_SIZE_400KB,REG_NUM,WEIGHT_TH,PKT_MARGIN>(traceFileDirCos,fileNum,topKs,outputDirPath,"Cos");
//    processMultipleTracesAndGetAvgMetrics<K_500KB,PZ_SIZE_500KB,REG_NUM,WEIGHT_TH,PKT_MARGIN>(traceFileDirCos,fileNum,topKs,outputDirPath,"Cos");
//
//    processMultipleTracesAndGetAvgMetrics<K_100KB,PZ_SIZE_100KB,REG_NUM,WEIGHT_TH,PKT_MARGIN>(traceFileDirMul,fileNum,topKs,outputDirPath,"Mul 2020.4");
//    processMultipleTracesAndGetAvgMetrics<K_200KB,PZ_SIZE_200KB,REG_NUM,WEIGHT_TH,PKT_MARGIN>(traceFileDirMul,fileNum,topKs,outputDirPath,"Mul 2020.4");
//    processMultipleTracesAndGetAvgMetrics<K_300KB,PZ_SIZE_300KB,REG_NUM,WEIGHT_TH,PKT_MARGIN>(traceFileDirMul,fileNum,topKs,outputDirPath,"Mul 2020.4");
//    processMultipleTracesAndGetAvgMetrics<K_400KB,PZ_SIZE_400KB,REG_NUM,WEIGHT_TH,PKT_MARGIN>(traceFileDirMul,fileNum,topKs,outputDirPath,"Mul 2020.4");
//    processMultipleTracesAndGetAvgMetrics<K_500KB,PZ_SIZE_500KB,REG_NUM,WEIGHT_TH,PKT_MARGIN>(traceFileDirMul,fileNum,topKs,outputDirPath,"Mul 2020.4");


    //multiple measurement points
//    processMultipleTracesAndGetMergedMetrics<K_500KB,PZ_SIZE_500KB,REG_NUM,WEIGHT_TH,PKT_MARGIN>(traceFileDir2016,fileNum,topKs,outputDirPath,"2016");
//    processMultipleTracesAndGetMergedMetrics<K_500KB,PZ_SIZE_500KB,REG_NUM,WEIGHT_TH,PKT_MARGIN>(traceFileDir2019,fileNum,topKs,outputDirPath,"2019");
//    processMultipleTracesAndGetMergedMetrics<K_500KB,PZ_SIZE_500KB,REG_NUM,WEIGHT_TH,PKT_MARGIN>(traceFileDirCos,fileNum,topKs,outputDirPath,"Cos");
    processMultipleTracesAndGetMergedMetrics<K_500KB,PZ_SIZE_500KB,REG_NUM,WEIGHT_TH,PKT_MARGIN>(traceFileDirMul,fileNum,topKs,outputDirPath,"Mul 2020.4");

    return 0;
}
