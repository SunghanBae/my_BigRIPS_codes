#include "TFile.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TTree.h"
#include "TString.h"
#include "TRandom.h"
#include "TCut.h"
#include "TCanvas.h"
#include "Riostream.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <utility>
#include <sstream>
#include <cmath>
#include <TProfile.h>
#include <TXMLNode.h>
#include <signal.h>
#include <time.h>

#include "BigRIPS.h"

class reco {

		public:
				time_t tm;
				int time0;
				
				//	char outfile[128];
				int i;
				Long64_t timestamp,evt;
				int up[3],down[2];
				
				int index;
				double aoq_0,daoq[4];

				double brho[2], delta[2],beta[2],gam[2],aoq[2];
				double xa35, xx35, xd35, xa511, xx511, xd511;
				double de_v, ionpair;
				double a1,alpha,tof,zet;
				double iccoef[2];
				//From 2014 PDG booklet
				const double c=299.792458;
				const double mu=931.494061;
				//matrix elements for delta calculation

				double l35,l511;

				double PPAC3X[2],PPAC5X[2],PPAC11X[2];
				double PPAC3Xdiff[2],PPAC5Xdiff[2],PPAC11Xdiff[2];
				double PPACTSumX[71],PPACTSumY[71];
				//double* PPACTSumX;
				//double* PPACTSumY;

				double reco_PPAC3X,reco_PPAC5X,reco_PPAC11X,reco_PPAC3A,reco_PPAC5A,reco_PPAC11A;
				double reco_PPAC5X_d,reco_PPAC11X_d;

				int array_map[12]={4,5,6,7,9,10,11,12,32,33,34,35};
				double c_cen[12],c_wid[12];

				bool x_flag[4];

				TCut pcut1[1];
				TCut pcut2[1];								//cut for sort out bad data

//				TGraphErrors* g_tofslope;

				Long64_t pre_ts,ts_now;
				
				reco();
				virtual ~reco();
				//int upstXA(double x1, double x2, double x3, double x4, double* x, double* d);
				//int downstX(double x1, double x2, double x3, double x4, double* x);
				int upstXA(double *xin, double *tin, int fl, double* x, double* d);
				int downstX(double *xin, double *tin, int fl, double* x);
				bool TSumCutSet(char* filename);
				void PID_draw(TString name);
				void recoPID(TTree* newt,BigRIPS &beam, double tofoff);
				void toff_slope(TTree* newt,BigRIPS &beam, double tofoff);
				void toff_fit(TTree* t, TGraphErrors* g);
				bool testx(double *xin,double *tin, int fl);


};

