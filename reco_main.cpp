#include "reco.cpp"

using namespace std;

int main(int argv, char* argc[]){
	
	if(argv<3 || argv>4){
	std::cout<<"Usage: ./reco runnumber tofoff (nameoption)"<<std::endl;
	return -1;
	}

	reco reco_manager;

	int runN = atoi(argc[1]);
	double tofoff= atof(argc[2]);
	char beamfile[128],TSumXCutFile[128];
		
	TString outfile= Form("$NP1306DIR/BigRIPS/run%04d_recoPID",runN);
	sprintf(beamfile,"$NP1306DIR/BigRIPS/run%04d.root",runN);
	sprintf(TSumXCutFile,"./run%04d_TSumXCut.dat",runN);

	outfile+=argc[3];
	outfile+=".root";
	TFile* f =new TFile(outfile,"RECREATE");
	TTree* t = new TTree("tree","tr");
	
	if(!reco_manager.TSumCutSet(TSumXCutFile)){
			cout<<"TSumX Cut file "<<TSumXCutFile<<" is missing! Stop!"<<endl;
			return -1;
	}

	BigRIPS beam;
	beam.GetTree(beamfile);
	
	reco_manager.recoPID(t,beam,tofoff);

	f->WriteTObject(t);
	std::cout<<outfile<<" is made"<<std::endl;
	f->WriteTObject(reco_manager.pcut1);
	f->WriteTObject(reco_manager.pcut2);
    
//	beam.CloseFile();
	f->Close();

	return 1;

}
