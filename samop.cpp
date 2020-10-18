//
// Created by 63175 on 2019/1/15.
//

#include <cstring>
#include <iostream>
#include <algorithm>
#include <cstdio>

#include "samop.h"


SAMOpt_t samGetOpt(string a){
    int stop = 0;
    SAMOpt_t res; res.tag=""; res.type=""; res.value="";
    string tmp = "";
    for(int i=0; i<a.length(); ++i){
        if(a[i]==':' || i==a.length()-1){
            ++stop;
            if(stop == 1) {
                res.tag = tmp;
            } else if(stop == 2) {
                res.type = tmp;
            } else if(stop == 3) {
                res.value = tmp;
            }
            tmp = ""; // 清空赋值的时候，内存会不会变，会不会有新对象产生？
        } else {
        	tmp.append(1, a[i]);
        }
    }
    return res;
}


void samGetNode(char *s, SAMNode_t *a) {
    a->options.clear();
    string tmp;
    int i, stop=0, len=(int)strlen(s);
    for(i=0; i<len; ++i) {
        if(s[i]=='\t' || s[i]=='\n' || i==len-1) {
            ++stop;
            if(stop == 1) { // qname
                a->qname = tmp;
            } else if(stop == 2) { // flag
                a->flag = StringToInt(tmp);
            } else if(stop == 3) { // rname
                a->rname = tmp;
            } else if(stop == 4) { // pos
                a->pos = StringToInt(tmp);
            } else if(stop == 5) { // MapQ
                a->mapq = StringToInt(tmp);
            } else if(stop == 6) { // CIGAR
                a->cigar = tmp;
            } else if(stop == 7) { // ref name of next read
                a->rnext = tmp;
            } else if(stop == 8) { // position of next read
                a->pnext = tmp;
            } else if(stop == 9) { // observed template length
                a->tlen = StringToInt(tmp);
            } else if(stop == 10) { // Seq
                a->seq = tmp;
            } else if(stop == 11) { // Qual
            	a->qual = tmp;
            } else if(stop >= 12) { // Options
                a->options.push_back(samGetOpt(tmp));
            }
            tmp = "";
        } else {
            tmp.append(1, s[i]);
        }
    }
}

void samGetNodes(FILE *f, vector<SAMNode_t> *v) {
    char buf[1024 * 1024];
    v->clear();
    int cntHeader = 0;
    while(fgets(buf, sizeof(buf), f) != NULL) {
    	// TODO: 处理头部信息
        if(buf[0] == '@') {
            ++cntHeader;
            continue;
        }
        SAMNode_t tmp;
        samGetNode(buf, &tmp);
        v->push_back(tmp);
    }
    printf("[%s] get %d SAM header information\n", __func__, cntHeader);
    printf("[%s] get %ld SAM nodes\n", __func__, v->size());
}

void samShowSingle(SAMNode_t *q) {
    cout << "qname:    " << q->qname << endl;
    if((q->flag&samFMulSeg) != 0) {
        cout << "    template having multiple segments in sequencing" << endl;
    } else if((q->flag&samFSegAli) != 0) {
        cout << "    each segment properly aligned according to the aligner" << endl;
    } else if((q->flag&samFUnmap) != 0) {
        cout << "    segment unmapped" << endl;
    } else if((q->flag&samFNxtSeq) != 0) {
        cout << "    next segment in the template unmapped" << endl;
    } else if((q->flag&samFRC) != 0) {
        cout << "    SEQ being reverse complemented" << endl;
    } else if((q->flag&samFNRC) != 0) {
        cout << "    SEQ of the next segment in the template being reverse complemented" << endl;
    } else if((q->flag&samFFst) != 0) {
        cout << "    the first segment in the template" << endl;
    } else if((q->flag&samFLst) != 0) {
        cout << "    the last segment in the template" << endl;
    } else if((q->flag&samFSecAli) != 0) {
        cout << "    secondary alignment" << endl;
    } else if((q->flag&samFNPF) != 0) {
        cout << "    not passing filters, such as platform/vendor quality controls" << endl;
    } else if((q->flag&samFPCR) != 0) {
        cout << "    PCR or optical duplicate" << endl;
    } else if((q->flag&samFSupAli) != 0) {
        cout << "    supplementary alignment" << endl;
    }
    cout << "rname:    " << q->rname << endl;
    cout << "pos:      " << q->pos << endl;
    cout << "mapq:     " << q->mapq << endl;
    cout << "CIGAR:    " << q->cigar << endl;
    cout << "rnext:    " << q->rnext << endl;
    cout << "pnext:    " << q->pnext << endl;
    cout << "tlen:     " << q->tlen << endl;
    cout << "seq:      " << q->seq << endl;
    cout << "qual:     " << q->qual << endl;
    vector<SAMOpt_t> *p = &q->options;
    for(int i=0; i<p->size(); ++i) {
        SAMOpt_t *opt = &p->at(i);
        cout << "    " << opt->tag << ":" << opt->type << ":" << opt->value << endl;
    }
    cout << endl << endl;
}

void samDedup(vector<SAMNode_t> *v) {
	// TODO: SAMTools是否支持对SAM的去冗余
    printf("[%s] get %ld SAM nodes\n", __func__, v->size());
}
