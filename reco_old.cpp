#include "reco.h"

reco::reco(){
}

reco::~reco(){  
}

int reco::TSumCutSet(char* filename){

	std::ifstream TSumCutFile(filename.data());
	std::string line;
	i=0;
	while(TSumCutFile.good()){
		getline(TSumCutFile,line);
		auto commentline= line.find("#");
		auto newline = line.substr(0,commentline);
		if(newline.size>0){
				cin>>line>>c_cen[i]>>c_wid[i];
				i++;
		}
	}
}

int reco::upstXA(double* xin, double* tin, int fl, double* x, double* d){
//int reco::upstXA(double x1, double x2, double x3, double x4, double* x, double* d){
		if( (x1<-500 && x2<-500) || (x3<-500 && x4<-500) ){		//this gives no A
				*x=-9999;
				*d=-9999;
				return -1;
		}else if( x1>-500 )
		{								//x1 alive
				if(x2>-500){						//x2 alive
						if(x3>-500){
								if(x4>-500){
										*x=0.25*(x1+x2+x3+x4);
										*d=0.5*(x3+x4-x1-x2);
										return 1;
								}else{
										*x=(x1+x2+x3)/3.;
										*d=x3-0.5*(x1+x2);
										return 1;
								}
						}else {
								*x=(x1+x2+x4)/3.;
								*d=x4-0.5*(x1+x2);
								return 1;
						}
				}else if(x3>-500){					//x2 dead
						if(x4>-500){
								*x=(x1+x3+x4)/3.;
								*d=0.5*(x3+x4)-x1;
								return 1;
						}else {
								*x=(x1+x3)*0.5;
								*d=x3-x1;
								return 1;
						}
				}else if(x4>-500){
						*x=(x1+x4)*0.5;
						*d=x4-x1;
						return 1;
				}
		}else{								// x1 dead, x2 alive
				if(x3>-500){					
						if(x4>-500){
								*x=(x2+x3+x4)/3.;
								*d=0.5*(x3+x4)-x2;
								return 1;
						}else{
								*x=(x2+x3)*0.5;
								*d=x3-x2;
								return 1;
						}
				}else {
						*x=(x2+x4)*0.5;
						*d=x4-x2;
						return 1;
				}
		}
}


int reco::downstX(double x1, double x2, double x3, double x4, double t1, double t2, double t3, double t4, double* x){
//int reco::downstX(double x1, double x2, double x3, double x4, double* x){
		int nnn=0;
		*x=0;
		if(x1>-500) {*x+=x1;nnn+=1;}
		if(x2>-500) {*x+=x2;nnn+=1;}
		if(x3>-500) {*x+=x3;nnn+=1;}
		if(x4>-500) {*x+=x4;nnn+=1;}
		if(nnn>0){
				*x=1.0*(*x)/nnn;
				return 1;
		}else{
				*x=-9999;
				return -1;
		}
}


void reco::recoPID(TTree* newt,BigRIPS &beam, double tofoff, TString nameoption =""){


		//run 1179 setting
		//iccoef[0]=1.459048;
		//iccoef[1]=1.628571;
		//run 1190 setting
		//iccoef[0]=1.459048*0.9794;
		//iccoef[1]=1.628571*0.9794-0.1346;
		//run 1181 setting
		//iccoef[0]=1.459048*0.9852;
		//iccoef[1]=1.628571*0.9852+0.03037;
		//run 1193 setting
		iccoef[0]=1.459048*0.9697;
		iccoef[1]=1.628571*0.9697+0.2369;

		ionpair= 4886.00; //From Kathrin BigRIPSIC xml, F11IC
		//	tofoff = 449.18;
		//	tofoff = 446.98; //for run1179
		//	tofoff = 444.38+0.6032; //for run1190
		//	tofoff = 443.99+3.3367; //for run1181
		//	tofoff = 443.56+1.413242; //for run1193
		l35  = 54917-31633 ;
		l511 = 125984-54917;


		//	TTree* btree = beam.fChain;

		newt->Branch("timestamp",&timestamp,"timestamp/L");
		newt->Branch("evt",&evt,"evt/L");
		newt->Branch("reco_PPAC3X",&reco_PPAC3X,"reco_PPAC3X/D");
		newt->Branch("reco_PPAC3A",&reco_PPAC3A,"reco_PPAC3A/D");
		newt->Branch("reco_PPAC5X",&reco_PPAC5X,"reco_PPAC5X/D");
		newt->Branch("reco_PPAC5A",&reco_PPAC5A,"reco_PPAC5A/D");
		newt->Branch("reco_PPAC11X",&reco_PPAC11X,"reco_PPAC11X/D");
		newt->Branch("reco_PPAC11A",&reco_PPAC11A,"reco_PPAC11A/D");
		newt->Branch("reco_PPAC5X_d",&reco_PPAC5X_d,"reco_PPAC5X_d/D");
		newt->Branch("reco_PPAC11X_d",&reco_PPAC11X_d,"reco_PPAC11X_d/D");

		newt->Branch("PPAC3X",PPAC3X,"PPAC3X[2]/D");
		newt->Branch("PPAC5X",PPAC5X,"PPAC5X[2]/D");
		newt->Branch("PPAC11X",PPAC11X,"PPAC11X[2]/D");
		newt->Branch("PPAC3Xdiff",PPAC3Xdiff,"PPAC3Xdiff[2]/D");
		newt->Branch("PPAC5Xdiff",PPAC5Xdiff,"PPAC5Xdiff[2]/D");
		newt->Branch("PPAC11Xdiff",PPAC11Xdiff,"PPAC11Xdiff[2]/D");
		newt->Branch("PPACTSumX",PPACTSumX,"PPACTSumX[71]/D");
		newt->Branch("PPACTSumY",PPACTSumY,"PPACTSumY[71]/D");

		newt->Branch("beta",beta,"beta[2]/D");
		newt->Branch("brho",brho,"brho[2]/D");
		//	newt->Branch("beta2",&beta2,"beta2/D");
		newt->Branch("gam",gam,"gam[2]/D");
		//	newt->Branch("gam2",&gam2,"gam2/D");
		newt->Branch("alpha",&alpha,"alpha/D");
		newt->Branch("delta",delta,"delta[2]/D");
		//	newt->Branch("delta2",&delta2,"delta2/D");
		newt->Branch("aoq",aoq,"aoq[2]/D");
		//	newt->Branch("aoq2",&aoq2,"aoq2/D");
		newt->Branch("zet",&zet,"zet/D");
		newt->Branch("tof",&tof,"tof/D");
		newt->Branch("flag_upPPAC",&up,"flag_upPPAC[3]/I");
		newt->Branch("flag_downPPAC",&down,"flag_downPPAC[2]/I");



		//	TString ppaccut;
		pcut1[0] = "PPAC3X[0]-0.5*PPAC3Xdiff[0] >-400 && PPAC3X[0]+0.5*PPAC3Xdiff[0]>-400 && PPAC5X[0]-0.5*PPAC5Xdiff[0] >-400 && PPAC5X[0]+0.5*PPAC5Xdiff[0]>-400 && PPAC11X[0]-0.5*PPAC11Xdiff[0] >-400 && PPAC11X[0]+0.5*PPAC11Xdiff[0]>-400 && PPAC3X[1]-0.5*PPAC3Xdiff[1] >-400 && PPAC3X[1]+0.5*PPAC3Xdiff[1]>-400 && PPAC5X[1]-0.5*PPAC5Xdiff[1] >-400 && PPAC5X[1]+0.5*PPAC5Xdiff[1]>-400 && PPAC11X[1]-0.5*PPAC11Xdiff[1] >-400 && PPAC11X[1]+0.5*PPAC11Xdiff[1]>-400";

		pcut2[0] = "PPAC3Xdiff[0] > -1 && PPAC3Xdiff[0]<1 && PPAC3Xdiff[1] >-1 && PPAC3Xdiff[1]<1 && PPAC5Xdiff[0] > -1 && PPAC5Xdiff[0]<1 && PPAC5Xdiff[1] >-1 && PPAC5Xdiff[1]<1 && PPAC11Xdiff[0] > -1 && PPAC11Xdiff[0]<1 && PPAC11Xdiff[1] >-1 && PPAC11Xdiff[1]<1";
		//	cut[3] = "BigRIPSPPAC.fTSumX[4]


		//	btree->Draw(">>li",cut[1]);
		//	TEventList* list = (TEventList*)gDirectory->Get("li");

		//	int nentries = list->GetN();

		int nentries = beam.GetEntriesFast();
		std::cout<<nentries<<" entries"<<std::endl;
		xx35=0.970669	;			xx511=1.06073	;
		xa35=0.0666164	;			xa511=0.073679	;
		xd35=64.2924	;			xd511=-68.1397	;

		//	tmin=btree->GetMinimum("EventInfo.timestamp");

		time0=time(&tm);
		for(i=0; i<nentries; i++){
				up[0]=0;
				up[1]=0;
				up[2]=0;
				down[0]=0;
				down[1]=0;

				if(i%100000 ==0) {std::cout<<"\rreco PID processing... "<<i<<"entires done. "<<time(&tm)-time0<<"s passed"/*<<" ADCSq: "<<beam.BigRIPSIC_fRawADCSqSum[2]<<" zet: "<<zet<<" aoq: "<<aoq2*/; std::cout.flush();}
				//	beam.GetEntry(list->GetEntry(i));
				beam.GetEntry(i);

				timestamp = beam.EventInfo_timestamp[0];
				evt=beam.EventInfo_eventnumber[0];
				ts_now=timestamp;

				tof = beam.BigRIPSPlastic_fTime[9]-beam.BigRIPSPlastic_fTime[1]+tofoff;

				PPAC3X[0] = 0.5*(beam.BigRIPSPPAC_fX[4]+beam.BigRIPSPPAC_fX[5]);
				PPAC3X[1] = 0.5*(beam.BigRIPSPPAC_fX[6]+beam.BigRIPSPPAC_fX[7]);
				PPAC5X[0] = 0.5*(beam.BigRIPSPPAC_fX[9]+beam.BigRIPSPPAC_fX[10]);
				PPAC5X[1] = 0.5*(beam.BigRIPSPPAC_fX[11]+beam.BigRIPSPPAC_fX[12]);
				PPAC11X[0] = 0.5*(beam.BigRIPSPPAC_fX[32]+beam.BigRIPSPPAC_fX[33]);
				PPAC11X[1] = 0.5*(beam.BigRIPSPPAC_fX[34]+beam.BigRIPSPPAC_fX[35]);

				PPAC3Xdiff[0]  = (-beam.BigRIPSPPAC_fX[4]+beam.BigRIPSPPAC_fX[5]);
				PPAC3Xdiff[1]  = (-beam.BigRIPSPPAC_fX[6]+beam.BigRIPSPPAC_fX[7]);
				PPAC5Xdiff[0]  = (-beam.BigRIPSPPAC_fX[9]+beam.BigRIPSPPAC_fX[10]);
				PPAC5Xdiff[1]  = (-beam.BigRIPSPPAC_fX[11]+beam.BigRIPSPPAC_fX[12]);
				PPAC11Xdiff[0] = (-beam.BigRIPSPPAC_fX[32]+beam.BigRIPSPPAC_fX[33]);
				PPAC11Xdiff[1] = (-beam.BigRIPSPPAC_fX[34]+beam.BigRIPSPPAC_fX[35]);
				//	if( PPAC3X[0]-0.5*PPAC3Xdiff[0] >-400 && PPAC3X[0]+0.5*PPAC3Xdiff[0]>-400 && PPAC5X[0]-0.5*PPAC5Xdiff[0] >-400 && PPAC5X[0]+0.5*PPAC5Xdiff[0]>-400 && PPAC11X[0]-0.5*PPAC11Xdiff[0] >-400 && PPAC11X[0]+0.5*PPAC11Xdiff[0]>-400 && PPAC3X[1]-0.5*PPAC3Xdiff[1] >-400 && PPAC3X[1]+0.5*PPAC3Xdiff[1]>-400 && PPAC5X[1]-0.5*PPAC5Xdiff[1] >-400 && PPAC5X[1]+0.5*PPAC5Xdiff[1]>-400 && PPAC11X[1]-0.5*PPAC11Xdiff[1] >-400 && PPAC11X[1]+0.5*PPAC11Xdiff[1]>-400 ){

				memcpy(PPACTSumX, beam.BigRIPSPPAC_fTSumX,sizeof(beam.BigRIPSPPAC_fTSumX));
				memcpy(PPACTSumY, beam.BigRIPSPPAC_fTSumY,sizeof(beam.BigRIPSPPAC_fTSumY));

			/*	up[0]=upstXA(beam.BigRIPSPPAC_fX[4],beam.BigRIPSPPAC_fX[5],beam.BigRIPSPPAC_fX[6],beam.BigRIPSPPAC_fX[7],&reco_PPAC3X,&reco_PPAC3A);
				up[1]=upstXA(beam.BigRIPSPPAC_fX[9],beam.BigRIPSPPAC_fX[10],beam.BigRIPSPPAC_fX[11],beam.BigRIPSPPAC_fX[12],&reco_PPAC5X,&reco_PPAC5A);
				up[2]=upstXA(beam.BigRIPSPPAC_fX[32],beam.BigRIPSPPAC_fX[33],beam.BigRIPSPPAC_fX[34],beam.BigRIPSPPAC_fX[35],&reco_PPAC11X,&reco_PPAC11A);

				down[0]=downstX(beam.BigRIPSPPAC_fX[9],beam.BigRIPSPPAC_fX[10],beam.BigRIPSPPAC_fX[11],beam.BigRIPSPPAC_fX[12],&reco_PPAC5X_d);
				down[1]=downstX(beam.BigRIPSPPAC_fX[32],beam.BigRIPSPPAC_fX[33],beam.BigRIPSPPAC_fX[34],beam.BigRIPSPPAC_fX[35],&reco_PPAC11X_d);*/
				
				up[0]=upstXA(beam.BigRIPSPPAC_fX,PPACTSumX,0,&reco_PPAC3X,&reco_PPAC3A);
				up[1]=upstXA(beam.BigRIPSPPAC_fX,PPACTSumX,1,&reco_PPAC5X,&reco_PPAC5A);
				up[2]=upstXA(beam.BigRIPSPPAC_fX,PPACTSumX,2,&reco_PPAC11X,&reco_PPAC11A);

				down[0]=downstX(beam.BigRIPSPPAC_fX,PPACTSumX,0,&reco_PPAC5X_d);
				down[1]=downstX(beam.BigRIPSPPAC_fX,PPACTSumX,1,&reco_PPAC11X_d);

				delta[0] = ( reco_PPAC5X_d
								-xa35*reco_PPAC3A/890.*1000
								-xx35*reco_PPAC3X
						   )/xd35;

				delta[1] = ( reco_PPAC11X_d
								-xa511*reco_PPAC5A/650.*1000
								-xx511*reco_PPAC5X
						   )/xd511;
				/*
				   delta[0] = ( 0.5*(PPAC5X[0]+PPAC5X[1])
				   -xa35*(PPAC3X[1]-PPAC3X[0])/890.*1000
				   -xx35*0.5*(PPAC3X[0]+PPAC3X[1])
				   )/xd35;

				   delta[1] = ( 0.5*(PPAC11X[0]+PPAC11X[1])
				   -xa511*(PPAC5X[1]-PPAC5X[0])/650.*1000
				   -xx511*0.5*(PPAC5X[0]+PPAC5X[1])
				   )/xd511;
				 */
				brho[0]  = 7.7565*(1+delta[0]*0.01);
				brho[1]  = 7.7490*(1+delta[1]*0.01);

				//		brho1  = 7.7565;
				//		brho2  = 7.7490;

				alpha  = brho[1]/brho[0];
				a1     = TMath::Sqrt(alpha*alpha*c*c*tof*tof
								+(TMath::Power(alpha,4)-alpha*alpha)*l35*l35+
								+(1-alpha*alpha)*l511*l511);

				beta[0]  = (l35*tof*c*alpha*alpha+l511*a1)/(alpha*alpha*c*c*tof*tof-l511*l511*(alpha*alpha-1));
				gam[0]=1./sqrt(1-pow(beta[0],2));

				beta[1]  = (l35*a1+l511*c*tof)/(tof*tof*c*c+l35*l35*(alpha*alpha-1));
				gam[1]=1./sqrt(1-pow(beta[1],2));

				aoq[0] = brho[0]*c/mu/beta[0]/gam[0];
				aoq[1] = brho[1]*c/mu/beta[1]/gam[1];

				de_v = TMath::Log(ionpair*beta[1]*beta[1])-TMath::Log((1-beta[1]*beta[1]))-beta[1]*beta[1];
				zet  = iccoef[0] *TMath::Sqrt(beam.BigRIPSIC_fRawADCSqSum[2]/de_v)*beta[1]+iccoef[1];

				newt->Fill();

				//		if(pre_ts > ts_now || ts_now==0) {
				//				newt->Reset();				//run1193 setting
				//				break;						//run1201 setting
				//		}
				//		pre_ts=ts_now;
				//}
		}
		//ppaccut.SetName("ppaccut");
		std::cout<<std::endl<<"job finished"<<std::endl;
		pcut1->SetName("pcut1");
		pcut2->SetName("pcut2");

}


void reco::PID_draw(TString name){
		TFile* f = (TFile*)gROOT->GetListOfFiles()->FindObject(name);
		if(!f || !f->IsOpen()){
				f =new TFile(name);
		}

		TTree* t = (TTree*)f->Get("tree");
		TCanvas* c =new TCanvas("c","",1200,600);
		//TCut p1=(TCut)f->Get("pcut1");
		t->Draw("zet:aoq[1]>>h(1000,2.52,2.74,1000,28,42)","","colz");
}
