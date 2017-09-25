void freeCorrespondings() {
    free(corresKIdxs);
    free(corresKIsRev);
    corresKIdxCnt = 0;
}

void correctReadsByKmerRead(int *toBeCorrectedRIds, short int *toBeCorrectedStPos, int toBeCorrectedCnt, int correctiveRId, short int correctivePos, short int *nonRelativeVarPos, short int nonRelativeVarPosLen) {
    int idxInRead = -1;
    uint64_t correctiveFreq = (kmersPtr[kIdxOfMaxFreq].rCnt + kmersPtr[kIdxOfMaxFreq].revRCnt);
    char crrNuc;
    for (int i = 0; i < toBeCorrectedCnt; i++) {
        for (int j = 0; j < nonRelativeVarPosLen; j++) {
            idxInRead = toBeCorrectedStPos[i] + nonRelativeVarPos[j];
            if (readsPtr[toBeCorrectedRIds[i]].crrFreq[idxInRead] < correctiveFreq) {
                readsPtr[toBeCorrectedRIds[i]].crrFreq[idxInRead] = correctiveFreq;
                crrNuc = readsPtr[correctiveRId].seq[correctivePos + nonRelativeVarPos[j]];
                if(checkDel){
                    if ((idxInRead - 1) >= 0 && readsPtr[toBeCorrectedRIds[i]].crr[idxInRead - 1] == crrNuc) { //nuc neigh = correct
                        readsPtr[toBeCorrectedRIds[i]].crr[idxInRead] = CRR_DEL; //remove
                    } else if ((correctivePos + nonRelativeVarPos[j] - 1) >= 0 
                            && readsPtr[toBeCorrectedRIds[i]].crr[idxInRead] == readsPtr[correctiveRId].seq[correctivePos + nonRelativeVarPos[j] - 1]) { //nuc = correct neigh
                        readsPtr[toBeCorrectedRIds[i]].crr[idxInRead] = (char)((int)crrNuc + 100); //insert
                    }  else if ((correctivePos + nonRelativeVarPos[j] + 1) < readsPtr[correctiveRId].len
                            && readsPtr[toBeCorrectedRIds[i]].crr[idxInRead] == readsPtr[correctiveRId].seq[correctivePos + nonRelativeVarPos[j] + 1]) { //nuc = correct neigh
                        readsPtr[toBeCorrectedRIds[i]].crr[idxInRead] = (char)((int)crrNuc + 100); //insert
                    } else {
                        readsPtr[toBeCorrectedRIds[i]].crr[idxInRead] = crrNuc; //subs
                    }
                }else{
                    readsPtr[toBeCorrectedRIds[i]].crr[idxInRead] = crrNuc; //subs
                }
            }
        }
    }
}

void investigateAndCorrectKmer(int toBeCorrectedKIdx, short int *nonRelativeVarPos, short int nonRelativeVarPosLen) {
    char nuc01 = 'N', revNuc01 = 'N', nuc02 = 'N', revNuc02 = 'N', nuc11 = 'N', revNuc11 = 'N', nuc22 = 'N', revNuc22 = 'N';
    short int checkPos1 = nonRelativeVarPos[0] + 1;
    short int checkPos2 = nonRelativeVarPos[nonRelativeVarPosLen - 1] - 1;
    if (kmersPtr[kIdxOfMaxFreq].rCnt > 0) {
        nuc01 = readsPtr[kmersPtr[kIdxOfMaxFreq].rIds[0]].seq[kmersPtr[kIdxOfMaxFreq].stPos[0] + checkPos1];
        nuc11 = readsPtr[kmersPtr[kIdxOfMaxFreq].rIds[0]].seq[kmersPtr[kIdxOfMaxFreq].stPos[0] + checkPos2];
    }
    if (kmersPtr[kIdxOfMaxFreq].revRCnt > 0) {
        revNuc01 = readsPtr[kmersPtr[kIdxOfMaxFreq].revRIds[0]].seq[kmersPtr[kIdxOfMaxFreq].revStPos[0] + checkPos1];
        revNuc11 = readsPtr[kmersPtr[kIdxOfMaxFreq].revRIds[0]].seq[kmersPtr[kIdxOfMaxFreq].revStPos[0] + checkPos2];
    }
    if (kmersPtr[toBeCorrectedKIdx].rCnt > 0) {
        nuc02 = readsPtr[kmersPtr[toBeCorrectedKIdx].rIds[0]].seq[kmersPtr[toBeCorrectedKIdx].stPos[0] + checkPos1];
        nuc22 = readsPtr[kmersPtr[toBeCorrectedKIdx].rIds[0]].seq[kmersPtr[toBeCorrectedKIdx].stPos[0] + checkPos2];
    }
    if (kmersPtr[toBeCorrectedKIdx].revRCnt > 0) {
        revNuc02 = readsPtr[kmersPtr[toBeCorrectedKIdx].revRIds[0]].seq[kmersPtr[toBeCorrectedKIdx].revStPos[0] + checkPos1];
        revNuc22 = readsPtr[kmersPtr[toBeCorrectedKIdx].revRIds[0]].seq[kmersPtr[toBeCorrectedKIdx].revStPos[0] + checkPos2];
    }

    if (nuc01 != 'N' && nuc11 != 'N') {
        if (nuc01 == nuc02 && nuc11 == nuc22) {
            correctReadsByKmerRead(kmersPtr[toBeCorrectedKIdx].rIds, kmersPtr[toBeCorrectedKIdx].stPos, kmersPtr[toBeCorrectedKIdx].rCnt, kmersPtr[kIdxOfMaxFreq].rIds[0], kmersPtr[kIdxOfMaxFreq].stPos[0], nonRelativeVarPos, nonRelativeVarPosLen);
        } else if (nuc01 == revNuc02 && nuc11 == revNuc22) {
            correctReadsByKmerRead(kmersPtr[toBeCorrectedKIdx].revRIds, kmersPtr[toBeCorrectedKIdx].revStPos, kmersPtr[toBeCorrectedKIdx].revRCnt, kmersPtr[kIdxOfMaxFreq].rIds[0], kmersPtr[kIdxOfMaxFreq].stPos[0], nonRelativeVarPos, nonRelativeVarPosLen);
        }
    }

    if (revNuc01 != 'N' && revNuc11 != 'N') {
        if (revNuc01 == nuc02 && revNuc11 == nuc22) {
            correctReadsByKmerRead(kmersPtr[toBeCorrectedKIdx].rIds, kmersPtr[toBeCorrectedKIdx].stPos, kmersPtr[toBeCorrectedKIdx].rCnt, kmersPtr[kIdxOfMaxFreq].revRIds[0], kmersPtr[kIdxOfMaxFreq].revStPos[0], nonRelativeVarPos, nonRelativeVarPosLen);
        } else if (revNuc01 == revNuc02 && revNuc11 == revNuc22) {
            correctReadsByKmerRead(kmersPtr[toBeCorrectedKIdx].revRIds, kmersPtr[toBeCorrectedKIdx].revStPos, kmersPtr[toBeCorrectedKIdx].revRCnt, kmersPtr[kIdxOfMaxFreq].revRIds[0], kmersPtr[kIdxOfMaxFreq].revStPos[0], nonRelativeVarPos, nonRelativeVarPosLen);
        }
    }
}

void correctKmerByMaxFreqKmer(int toBeCorrectedKIdx, bool toBeCorrectedKIsRev, short int *nonRelativeVarPos, short int nonRelativeVarPosLen) {
    bool notCorrectedYet = false;
    if (!kOfMaxFreqIsRev) {
        if (kmersPtr[kIdxOfMaxFreq].rCnt > 0) {
            if (!toBeCorrectedKIsRev) {
                if (kmersPtr[toBeCorrectedKIdx].rCnt > 0) {
                    correctReadsByKmerRead(kmersPtr[toBeCorrectedKIdx].rIds, kmersPtr[toBeCorrectedKIdx].stPos, kmersPtr[toBeCorrectedKIdx].rCnt, kmersPtr[kIdxOfMaxFreq].rIds[0], kmersPtr[kIdxOfMaxFreq].stPos[0], nonRelativeVarPos, nonRelativeVarPosLen);
                } else {
                    notCorrectedYet = true;
                }
            } else {
                if (kmersPtr[toBeCorrectedKIdx].revRCnt > 0) {
                    correctReadsByKmerRead(kmersPtr[toBeCorrectedKIdx].revRIds, kmersPtr[toBeCorrectedKIdx].revStPos, kmersPtr[toBeCorrectedKIdx].revRCnt, kmersPtr[kIdxOfMaxFreq].rIds[0], kmersPtr[kIdxOfMaxFreq].stPos[0], nonRelativeVarPos, nonRelativeVarPosLen);
                } else {
                    notCorrectedYet = true;
                }
            }
        } else {
            notCorrectedYet = true;
        }
    } else {

        if (kmersPtr[kIdxOfMaxFreq].revRCnt > 0) {
            if (!toBeCorrectedKIsRev) {
                if (kmersPtr[toBeCorrectedKIdx].rCnt > 0) {
                    correctReadsByKmerRead(kmersPtr[toBeCorrectedKIdx].rIds, kmersPtr[toBeCorrectedKIdx].stPos, kmersPtr[toBeCorrectedKIdx].rCnt, kmersPtr[kIdxOfMaxFreq].revRIds[0], kmersPtr[kIdxOfMaxFreq].revStPos[0], nonRelativeVarPos, nonRelativeVarPosLen);
                } else {
                    notCorrectedYet = true;
                }
            } else {
                if (kmersPtr[toBeCorrectedKIdx].revRCnt > 0) {
                    correctReadsByKmerRead(kmersPtr[toBeCorrectedKIdx].revRIds, kmersPtr[toBeCorrectedKIdx].revStPos, kmersPtr[toBeCorrectedKIdx].revRCnt, kmersPtr[kIdxOfMaxFreq].revRIds[0], kmersPtr[kIdxOfMaxFreq].revStPos[0], nonRelativeVarPos, nonRelativeVarPosLen);
                } else {
                    notCorrectedYet = true;
                }
            }
        } else {
            notCorrectedYet = true;
        }
    }
    if (notCorrectedYet) {
        investigateAndCorrectKmer(toBeCorrectedKIdx, nonRelativeVarPos, nonRelativeVarPosLen);
    }
}

double getKmerVarsAvgQv(int kIdx, short int *nonRelativeVarPos, short int varPosLen) {
    double avg = 0.0;
    int cnt = 0;
    int i, j;
    for (i = 0; i < kmersPtr[kIdx].rCnt; i++) {
        for (j = 0; j < varPosLen; j++) {
            avg += readsPtr[kmersPtr[kIdx].rIds[i]].qv[kmersPtr[kIdx].stPos[i] + nonRelativeVarPos[j]];
            cnt++;
        }
    }
    for (i = 0; i < kmersPtr[kIdx].revRCnt; i++) {
        for (j = 0; j < varPosLen; j++) {
            avg += readsPtr[kmersPtr[kIdx].revRIds[i]].qv[kmersPtr[kIdx].revStPos[i] + nonRelativeVarPos[j]];
            cnt++;
        }
    }
    return avg / cnt;
}

void generateSimilarKmers(int rIdx, short int stPos, short int *varPos, short int varPosLen, short int k, short int crrVarPos) {

    char varNucs[4] = {'A', 'C', 'G', 'T'};
    //if it's the first iteration, the first set of combinations is the same as the set of characters
    if (corresKStrCnt == 0) {
        corresKStrCnt = 4;
        corresKStrs = (char **) malloc(corresKStrCnt * sizeof (char *));
        if (corresKStrs == 0) {
            printf("ERROR in allocating corresKStrs: Out of memory\n");
            return;
        }
        int j;
        for (int i = 0; i < corresKStrCnt; i++) {
            corresKStrs[i] = (char *) malloc(k * sizeof (char));
            if (corresKStrs[i] == 0) {
                printf("ERROR in allocating corresKStrs[i]: Out of memory\n");
                return;
            }
            for (j = 0; (j + stPos) < varPos[0] && j < k; j++) {
                corresKStrs[i][j] = readsPtr[rIdx].seq[(j + stPos)];
            }
            if ((j + stPos) == varPos[0]) {
                corresKStrs[i][j] = varNucs[i];
            }
        }
    }

    //we're done if we're at size 1
    if (crrVarPos >= varPosLen) {
        for (int i = 0; i < corresKStrCnt; i++) {
            for (int j = varPos[varPosLen - 1] + 1; j < k; j++) {
                corresKStrs[i][j] = readsPtr[rIdx].seq[j];
            }
        }
        return;
    }

    //loop through existing combinations and character set to create strings
    int prevcorresKStrCnt = corresKStrCnt;
    int prevVarPos = crrVarPos - 1;
    corresKStrCnt = 4 * corresKStrCnt;
    corresKStrs = (char **) realloc(corresKStrs, corresKStrCnt * sizeof (char *));
    if (corresKStrs == 0) {
        printf("ERROR in reallocating corresKStrs: Out of memory\n");
        return;
    }
    for (int w = prevcorresKStrCnt; w < corresKStrCnt; w++) {
        corresKStrs[w] = (char *) malloc(k * sizeof (char));
        if (corresKStrs[w] == 0) {
            printf("ERROR in allocating corresKStrs[i]: Out of memory\n");
            return;
        }
    }
    int z, m, y;
    for (int x = 0; x < prevcorresKStrCnt; x++) {
        m = x;
        y = 0;
        while (m < corresKStrCnt && y < 4) {
            for (z = 0; (z + stPos) < varPos[crrVarPos] && z < k; z++) {
                if ((z + stPos) <= varPos[prevVarPos]) {
                    corresKStrs[m][z] = corresKStrs[x][z];
                } else {
                    corresKStrs[m][z] = readsPtr[rIdx].seq[(z + stPos)];
                }
            }
            if ((z + stPos) == varPos[crrVarPos]) {
                corresKStrs[m][z] = varNucs[y];
            }
            m += prevcorresKStrCnt;
            y += 1;
        }
    }

    //call same function again for the next iteration
    generateSimilarKmers(rIdx, stPos, varPos, varPosLen, k, (crrVarPos + 1));
}

void getCorrespondingKIdxs(int rIdx, short int stPos, short int *nonRelativeVarPos, short int *varPos, short int varPosLen, int orgKIdx, bool orgKIsRev, short int k) {

    generateSimilarKmers(rIdx, stPos, varPos, varPosLen, k, 1);


    //        cout << "---------------------------------------------------------------" << endl;
    //        cout << "CorrespondingCnt: " << corresKStrCnt << endl;
    //        for (int i = 0; i < corresKStrCnt; i++) {
    //            cout << i<<"::";
    //            for(int j = 0; j< k; j++){
    //                    cout << corresKStrs[i][j];
    //            }
    //            cout  << endl;
    //        }
    //        cout << "---------------------------------------------------------------" << endl;

    corresKIdxCnt = 1;
    corresKIdxs = (int *) malloc(corresKIdxCnt * sizeof (int));
    corresKIsRev = (bool *) malloc(corresKIdxCnt * sizeof (bool));
    if (corresKIdxs == 0 || corresKIsRev == 0) {
        printf("ERROR in allocating corresKIdxs/corresKIsRev: Out of memory\n");
        return;
    }
    corresKIdxs[0] = orgKIdx;
    corresKIsRev[0] = orgKIsRev;

    int varKIdx = -1;
    bool kIsRev = false;
    kIdxOfMaxFreq = orgKIdx;
    kOfMaxFreqIsRev = orgKIsRev;
    int maxFreq = kmersPtr[kIdxOfMaxFreq].rCnt + kmersPtr[kIdxOfMaxFreq].revRCnt;
    int freq = 0;
    for (int i = 0; i < corresKStrCnt; i++) {
        varKIdx = hashKmer(corresKStrs[i], k, -1, true, k, true, &kIsRev);
        free(corresKStrs[i]);
        if (varKIdx > -1 && varKIdx != orgKIdx) {
            corresKIdxCnt += 1;
            corresKIdxs = (int *) realloc(corresKIdxs, corresKIdxCnt * sizeof (int));
            corresKIsRev = (bool *) realloc(corresKIsRev, corresKIdxCnt * sizeof (bool));
            if (corresKIdxs == 0 || corresKIsRev == 0) {
                printf("ERROR in reallocating corresKIdxs/corresKIsRev: Out of memory\n");
                return;
            }
            corresKIdxs[corresKIdxCnt - 1] = varKIdx;
            corresKIsRev[corresKIdxCnt - 1] = kIsRev;
            freq = kmersPtr[varKIdx].rCnt + kmersPtr[varKIdx].revRCnt;
            if (freq > maxFreq) {
                kIdxOfMaxFreq = varKIdx;
                kOfMaxFreqIsRev = kIsRev;
                maxFreq = kmersPtr[kIdxOfMaxFreq].rCnt + kmersPtr[kIdxOfMaxFreq].revRCnt;
            } else if (freq == maxFreq) {
                if (getKmerVarsAvgQv(varKIdx, nonRelativeVarPos, varPosLen) > getKmerVarsAvgQv(kIdxOfMaxFreq, nonRelativeVarPos, varPosLen)) {
                    kIdxOfMaxFreq = varKIdx;
                    kOfMaxFreqIsRev = kIsRev;
                    maxFreq = kmersPtr[kIdxOfMaxFreq].rCnt + kmersPtr[kIdxOfMaxFreq].revRCnt;
                }
            }
        }
    }
    free(corresKStrs);
    corresKStrCnt = 0;
}

void corrrectWithCorrespondings(int kIdx, int k, short int *nonRelativeVarPos, short int varPosCnt) {
    int stPos;
    int rIdx;
    bool kIsRev;
    if (kmersPtr[kIdx].rCnt > 0) {
        rIdx = kmersPtr[kIdx].rIds[0];
        stPos = kmersPtr[kIdx].stPos[0];
        kIsRev = false;
    } else {
        rIdx = kmersPtr[kIdx].revRIds[0];
        stPos = kmersPtr[kIdx].revStPos[0];
        kIsRev = true;
    }

    short int varPos[varPosCnt];
    for (short int x = 0; x < varPosCnt; x++) {
        varPos[x] = stPos + nonRelativeVarPos[x];
    }

    getCorrespondingKIdxs(rIdx, stPos, nonRelativeVarPos, varPos, varPosCnt, kIdx, kIsRev, k);

    for (int i = 0; i < corresKIdxCnt; i++) {
        if (corresKIdxs[i] != kIdxOfMaxFreq) {
            correctKmerByMaxFreqKmer(corresKIdxs[i], corresKIsRev[i], nonRelativeVarPos, 3);
            kmersPtr[corresKIdxs[i]].validByKIdx = kIdxOfMaxFreq;
        } else if (corresKIdxCnt > 1) {
            kmersPtr[corresKIdxs[i]].validByKIdx = kIdxOfMaxFreq;
        }
    }
}

void correctSolidKmers(int k) {

    short int nonRelativeVarPos[3] = {0, ((k - 1) / 2), k - 1};

    for (int j = 0; j < kmersCnt; j++) {
        if (kmersPtr[j].validByKIdx == -1) {

            corrrectWithCorrespondings(j, k, nonRelativeVarPos, 3);

            freeCorrespondings();
            kIdxOfMaxFreq = -1;
            kOfMaxFreqIsRev = false;
        }
    }
}

void correctMismatchedKmers(short int k) {
    short int varPos[1];
    short int cnt = 3;
    short int lowestPos[cnt];
    double lowestVal[cnt];
    short int j, x;
    double avgQV;
    kIdxOfMaxFreq = -1;
    kOfMaxFreqIsRev = false;
    for (int i = 0; i < kmersCnt; i++) {
        if (kmersPtr[i].validByKIdx == -1) {
            for (x = 0; x < cnt; x++) {
                lowestPos[x] = -1;
                lowestVal[x] = ((double) '~') + 1.0;
            }
            for (j = 0; j < k; j++) {
                varPos[0] = j;
                avgQV = getKmerVarsAvgQv(i, varPos, 1);
                if (avgQV < lowestVal[1]) {
                    if (avgQV < lowestVal[0]) {
                        lowestVal[0] = avgQV;
                        lowestPos[0] = j;
                    } else {
                        lowestVal[1] = avgQV;
                        lowestPos[1] = j;
                    }
                } else if (avgQV < lowestVal[2]) {
                    lowestVal[2] = avgQV;
                    lowestPos[2] = j;
                }
            }
            corrrectWithCorrespondings(i, k, lowestPos, cnt);
            freeCorrespondings();
            kIdxOfMaxFreq = -1;
            kOfMaxFreqIsRev = false;
        }
    }
}

void updateReadsByCorrections() {
    short int oldLen = 0;
    string oldSeq, oldQV;
    short int y;
    for (int i = 0; i < readsCnt; i++) {
        oldLen = readsPtr[i].len;
        oldSeq = readsPtr[i].seq;
        oldQV = readsPtr[i].qv;
        //Setting the corrected sequence
        readsPtr[i].len = 0;
        for (y = 0; y < oldLen; y++) {
            if (readsPtr[i].crr[y] != CRR_DEL) {
                readsPtr[i].seq = (char *) realloc(readsPtr[i].seq, (readsPtr[i].len + 1) * sizeof (char));
                readsPtr[i].qv = (char *) realloc(readsPtr[i].qv, (readsPtr[i].len + 1) * sizeof (char));
                if (readsPtr[i].crr[y] == CRR_NO) {
                    readsPtr[i].seq[readsPtr[i].len] = oldSeq[y];
                    readsPtr[i].qv[readsPtr[i].len] = oldQV[y];
                } else if (readsPtr[i].crr[y] >= CRR_SUBS_A && readsPtr[i].crr[y] <= CRR_SUBS_t) {
                    readsPtr[i].seq[readsPtr[i].len] = readsPtr[i].crr[y];
                    readsPtr[i].qv[readsPtr[i].len] = (oldQV[y] < 'c') ? 'c' : oldQV[y];
                } else {// if ((int)(readsPtr[i].crr[y] - 100) >= CRR_SUBS_A && (int)(readsPtr[i].crr[y] - 100) <= CRR_SUBS_t) {
                    readsPtr[i].seq[readsPtr[i].len] = (char) ((int) readsPtr[i].crr[y] - 100);
                    readsPtr[i].qv[readsPtr[i].len] = (oldQV[y] < 'c') ? 'c' : oldQV[y];
                    readsPtr[i].len += 1;
                    readsPtr[i].seq = (char *) realloc(readsPtr[i].seq, (readsPtr[i].len + 1) * sizeof (char));
                    readsPtr[i].qv = (char *) realloc(readsPtr[i].qv, (readsPtr[i].len + 1) * sizeof (char));
                    readsPtr[i].seq[readsPtr[i].len] = oldSeq[y];
                    readsPtr[i].qv[readsPtr[i].len] = oldQV[y];
                }
                readsPtr[i].len += 1;
            }
        }
        readsPtr[i].seq = (char *) realloc(readsPtr[i].seq, (readsPtr[i].len + 1) * sizeof (char));
        readsPtr[i].qv = (char *) realloc(readsPtr[i].qv, (readsPtr[i].len + 1) * sizeof (char));
        readsPtr[i].seq[readsPtr[i].len] = '\0';
        readsPtr[i].qv[readsPtr[i].len] = '\0';
        //Resetting corrections
        readsPtr[i].crr = (char *) realloc(readsPtr[i].crr, (readsPtr[i].len) * sizeof (char));
        readsPtr[i].crrFreq = (int *) realloc(readsPtr[i].crrFreq, (readsPtr[i].len) * sizeof (int));
        if (readsPtr[i].crr == 0 || readsPtr[i].crrFreq == 0) {
            printf("ERROR in reallocating a readsPtr[i].crr/readsPtr[i].crrFreq: Out of memory\n");
            return;
        }
        for (y = 0; y < readsPtr[i].len; y++) {
            readsPtr[i].crr[y] = CRR_NO;
            readsPtr[i].crrFreq[y] = 0;
        }
    }
}
/*
TP = 16,902,900
FP = 354,564
FN = 9,907,200
TN = 130,157,136
Sensitivity = 63.0468%
Specificity = 99.7283%
Accuracy = 61.7243%
Gain = 61.7243%
        (16,902,900+130,157,136)/(16,902,900+130,157,136+9,907,200+354,564) = 0.93477214219
        
TP = 16,904,268
FP = 355,536
FN = 9,905,832
TN = 130,156,164
Sensitivity = 63.0519%
Specificity = 99.7276%
Accuracy = 61.7257%
Gain = 61.7257%
                (16,904,268+130,156,164)/(16,904,268+130,156,164+9,905,832+355,536) = 0.93477465932	
        
TP = 13,670,784
FP = 542,880
FN = 13,139,316
TN = 129,968,820
Sensitivity = 50.9912%
Specificity = 99.584%
Accuracy = 48.9663%
Gain = 48.9663%
(13,670,784+129,968,820)/(13,670,784+13,139,316+542,880+129,968,820)    = 0.91303051452	    
*/
