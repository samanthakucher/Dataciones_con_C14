void CousinsHighlandAMS2(){

#include <iostream>

int Nm=22;
double_t RmObs=3.045E-13;
double_t ErrRmObs=0.11/pow(Nm,0.5); // Error relativo

double_t Edadmax=14000;

int Ns=8;
double_t RsObs=6.95E-13;
double_t ErrRsObs=0.1022/pow(Ns,0.5); // Error relativo

int Nf=9;
double_t RfObs=1.068E-13;
double_t ErrRfObs=0.38/pow(Nf,0.5);  // Error relativo

//std::cout << std::setprecision(5);
//std::cout << "scientific:\n" << std::scientific;
std::cout << "default:\n" << std::scientific;

cout<<"This calculations were perfomed for:"<<endl;
cout<<"Sample Ratio =    "<<RmObs <<"+/-"<<ErrRmObs*RmObs<<endl;
cout<<"Standar Ratio =   "<<RsObs <<"+/-"<<ErrRsObs*RsObs<<endl;
cout<<"Background Ratio ="<<RfObs <<"+/-"<<ErrRfObs*RfObs<<endl<<endl;
        
//This scripts uses hybrid toy MC method: bayesian for systematics, frequentist for signal 

int steps = 50;
double T = 5700/0.693147; // 14C lifetime in years

double Edad[steps];
double alpha[steps];
double beta[steps];
//double proba1[steps]; 
//double proba2[steps]; 

// each step corresponds to a different value for the signal s
// and therefore, to a different value of  u =(efiObs*s[k]+bObs)*T
for (int k=0;k<steps;++k){ 

  Edad[k]=k*1./(steps-1)*Edadmax;
  //std::cout << "default:\n" << std::scientific;
  //cout<<"Edad en este paso: "<<Edad[k]<<endl;
  
  // mu_hat corresponds to the case where Rs and Rf has not uncertainty.
  double mu_hat = RfObs+(RsObs-RfObs)*exp(-Edad[k]/T);
  //std::cout << "default:\n" << std::scientific;
  //cout<<"mu calculado sin errores en Rf y Rs: "<<mu_hat<<endl;
		
	  /*
	  // Probability of obtaining the observed n or less 
	  for (int m=0;m<nObs+1;++m) {
		  proba1[k] += exp(m*log(mu_hat)-mu_hat)/TMath::Factorial(m);
	  }
	  proba1[k]=proba1[k]*100;
	  proba2[k]=100-proba1[k];
	  */
	  
  // Here starts the loop for the Monte Carlo simulation
  // or computational "integration" over systematic errors	   
  TRandom *rnd=new TRandom3();
  int nToys = 10000;  // number of toys
  double nexp1 = 0;  // counter for Rm<=RmObs events 
  //double nexp2 = 0;  // counter for Rm>=RmObs events 
 
	  for(int iToy=0; iToy<nToys; iToy++) {

		// get Ratio in the standar from truncated Gaussian
		double Rs=0;
		if(ErrRsObs==0.0) 
		   Rs=RsObs;
		else
		   while (Rs<=0) Rs=rnd->Gaus(RsObs,ErrRsObs*RsObs);

		// get Ratio in the background from truncated Gaussian 
		double Rf=0;
		if(ErrRfObs==0.0) 
		   Rf=RfObs;
		else
		   while (Rf<=0) Rf=rnd->Gaus(RfObs,ErrRfObs*RfObs);

		// Generate number of events from Gauss(mu, ErrRmObs);
		// This Gaussian represents the likelihood in this case
		
		double mu = Rf+(Rs-Rf)*exp(-1*Edad[k]/T);
		//std::cout << "default:\n" << std::scientific;
		//cout<<"mu del Montecarlo: "<<mu<<endl;
		
		double R = rnd->Gaus(mu, ErrRmObs*mu);
		//std::cout << "default:\n" << std::scientific;
		//cout<<"R del Montecarlo:            "<<R<<endl;
		//std::cout << "default:\n" << std::scientific;
		//cout<<"RmObs contra el que compara: "<<RmObs<<endl;
		
		// Count N<=nObs events  
		if(R<=RmObs) nexp1 += 1;
		
		// Count N>=nObs events  
		//if(R>=RmObs) nexp2 += 1;
	  }

        alpha[k] = nexp1 / nToys*100; // That is (1- CL)/2
        //beta[k] = nexp2 / nToys*100; // That is (1- CL)/2
        beta[k]=100-alpha[k];

        cout.precision(3);      // Print using four decimals precision
        cout<<fixed;            // ...exactly four decimals

        cout<<"For a age equal to "<<Edad[k]<<" (1-CL)/2 is "<<alpha[k]<<" %"<<" (1-CL)/2 is "<<beta[k]<<endl;

}
                
/////////////////////////////////////////////////////////////////////////////////////
    TCanvas *c1;
    c1 = new TCanvas("canvas1","canvas1",800,800);
    c1->Clear();
    c1->SetGrid();

    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle("Integral a izquierda (ROJO) y a derecha (AZUL) hasta el valor observado; Valor teorico de la Edad; Probabilidad acumulada %");
        
            TGraph *gr1 = new TGraph(steps,Edad,alpha);
            //gr1->SetTitle("Confidence level as a function of s; 1-CL ; Signal");
			//gr1->SetMarkerColor(kRed);
			//gr1->GetXaxis()->SetTitle("UL for the signal, cosmic UHE neutrino flux");
			//gr1->GetYaxis()->SetTitle("(1 - Confidence Level)/2 %");
			gr1->GetXaxis()->SetLimits(0.0,Edadmax);
			gr1->GetYaxis()->SetRangeUser(0.0,100);
			gr1->SetLineColor(2);
			gr1->SetLineWidth(2);
			gr1->Draw("AP*");
			
			TGraph *gr2 = new TGraph(steps,Edad,beta);
            //gr2->SetTitle("Confidence level as a function of s; 1-CL ; Signal");
			//gr2->SetMarkerColor(kBlue);
			//gr2->GetXaxis()->SetTitle("UL for the signal, cosmic UHE neutrino flux");
			//gr2->GetYaxis()->SetTitle("(1 - Confidence Level)/2 %");
			gr2->GetXaxis()->SetLimits(0.0,Edadmax);
			gr2->GetYaxis()->SetRangeUser(0.0,100);
			gr2->SetLineColor(4);
			gr2->SetLineWidth(2);
			gr2->Draw("AP*");
			
			//TGraph *gr3 = new TGraph(steps,Edad,proba1);
			//TGraph *gr4 = new TGraph(steps,Edad,proba2);
					
			auto l1 = new TLine(0,16,Edadmax,16);
			l1->SetLineStyle(kDashed);
			
			auto l2 = new TLine(0,84,Edadmax,84);
			l2->SetLineStyle(kDashed);
					
			mg->Add(gr1); 
			mg->Add(gr2); 
			//mg->Add(gr3); 
			//mg->Add(gr4); 
			
			mg->Draw("AL*");
			l1->Draw("same");
			l2->Draw("same");
}
