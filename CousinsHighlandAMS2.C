void CousinsHighlandAMS2(){

#include <iostream>

/*
Este codigo calcula la incerteza en la edad de una muestra medida por AMS 
usando tres diferentes metodos basados en la verosimilitud.

- Bayesiano. Aqui los nuisance parameters se eliminan por integracion hecha por Monte Carlo.
  Luego, se invierte la likelihood para finalmente buscar el rango que integra una probabilidad del 68 % centrada

- Frecuentista - Hace profiling para eliminar los nuisance parameters
  Para encontrar las relaciones isotopicas de fondo y estandar que maximizan la verosimilitud
  usa los ensayos hechos Montecarlo del caso anterior realizados para integrar los nuisance parameters
  
- Hibrido: Usa la likelihood que quedo luego de eliminar los nuisance parameters por integracion bayesiana
  y el metodo del cinturon frecuentista para calcular el itervalo
   
*/

// A continuacion se indican los resultados experimentales
// Relacion isotopica medida en la muestra con su error relativo
int Nm=22;
double_t RmObs=3.045E-13;
double_t ErrRmObs=0.11/pow(Nm,0.5); 

// Relacion isotopica medida en el estandar con su error relativo
int Ns=8;
double_t RsObs=6.95E-13;
double_t ErrRsObs=0.1022/pow(Ns,0.5);

// Relacion isotopica medida en el fondo con su error relativo
int Nf=9;
double_t RfObs=1.068E-13;
double_t ErrRfObs=0.38/pow(Nf,0.5); 

double sigmaRm = ErrRmObs*RmObs;
double sigmaRs = ErrRsObs*RsObs;
double sigmaRf = ErrRfObs*RfObs;
		
int steps = 100;
double T = 5700/0.693147; // vida media del 14C
int nbines = 201;

std::cout << std::scientific;

cout<<endl;
cout<<"Este calculo se realizo para:"<<endl;
cout<<"Relacion isotopica en la muestra          = "<<RmObs <<"+/-"<<ErrRmObs*RmObs<<endl;
cout<<"Relacion isotopica en el estandar         = "<<RsObs <<"+/-"<<ErrRsObs*RsObs<<endl;
cout<<"Relacion isotopica en el fondo            = "<<RfObs <<"+/-"<<ErrRfObs*RfObs<<endl;
double fecha = -T*log((RmObs-RfObs)/(RsObs-RfObs));
cout<<"La edad calculada con estos parametros es = "<< fecha  <<" anios"<<endl<<endl;

// Edades maximas y minimas a considerar al hacer el barrido en este parametro
double_t Edadmin=fecha-1000;
double_t Edadmax=fecha+1000;

double Edad[steps];
double alpha[steps];
double beta[steps];
double LL[steps];
double R_bin[nbines][steps];

double Proba = 0;
double integral = 0;
double paso = 0;

bool entre = true;
bool entre1 = true;
bool entre2 = true;
bool lodije = true;
bool lodijeh = true;
bool pentre1 = true;
bool plodije = true;

TCanvas *c2;
c2 = new TCanvas("canvas2","canvas2",800,800);
TH1D * verosimilitud  =  new TH1D("Verosimilitud", "Verosimilitud", nbines, RmObs-10*ErrRmObs*RmObs, RmObs+10*ErrRmObs*RmObs);

cout<<endl;
cout<<"Los resultados para el calculo Bayesiano-Frecuentista son los siguientes:"<<endl;
		
// cada paso corresponde a un valor diferente de la Edad
// y por lo tanto, a un valor diferente de la relacion isotopica esperada en la muestra		
for (int k=0;k<steps;++k){ 
  
  paso = 1./(steps-1)*(Edadmax-Edadmin);
    
  Edad[k]=Edadmin+k*paso;
  //std::cout << std::scientific;
  //cout<<"Edad en este paso: "<<Edad[k]<<endl;
  
  // mu_hat corresponde al caso donde Rs y Rf no tienen incerteza.
  double mu_hat = RfObs+(RsObs-RfObs)*exp(-Edad[k]/T);
  //std::cout << std::scientific;
  //cout<<"mu calculado sin errores en Rf y Rs: "<<mu_hat<<endl;
	  
  // Aca empieza el loop para la integracion Bayesiana por MonteCarlo 
  // realizada sobre los nuisance parameters o errores sistematicos
  TRandom *rnd=new TRandom3();
  int nToys = 100000;// numero de toys
  double nexp = 0;  // contara eventos donde Rm<=RmObs 

  LL[k]=0; //Inicializo la Likelihood que usare para profiling
  
	  for(int iToy=0; iToy<nToys; iToy++) {

		// obtiene la relacion isotopica en el estandar a partir de una gaussiana truncada
		double Rs=0;
		if(ErrRsObs==0.0) 
		   Rs=RsObs;
		else
		   while (Rs<=0) Rs=rnd->Gaus(RsObs,ErrRsObs*RsObs);

		// obtiene la relacion isotopica en el fondo a partir de una gaussiana truncada
		double Rf=0;
		if(ErrRfObs==0.0) 
		   Rf=RfObs;
		else
		   while (Rf<=0) Rf=rnd->Gaus(RfObs,ErrRfObs*RfObs);
	
		//Calcula el valor esperado para la relacion isotopica de la muestra
		//usando el par de valores obtenidos para Rs y Rf
		double mu = Rf+(Rs-Rf)*exp(-Edad[k]/T);
		//std::cout << std::scientific;
		//cout<<"Rm del Montecarlo: "<<mu<<endl;
		
		// Ahora se genera un valor R para la relacion isotopica medida Rm
		// Dicho valor se calcula usando la likelihood como distribucion de probabilidad
		// Esa likelihood se va considerar una Gaussiana de la forma Gaus(mu, ErrRmObs);
		
		// Uno podria preguntarte donde estan las otras dos distribuciones, la de Rf y Rs
		// que analiticamente hay que multiplicar por la de Rm para obtener la verosimilitud,
		// las mismas fueron incorporadas al darle peso al Rf y al Rs con que se calculo el mu
		// de la distribucion de la que proviene ahora R.
		
		double sigma = ErrRmObs*mu;
		
		double R = rnd->Gaus(mu, sigma);
		//std::cout << std::scientific;
		//cout<<"R del Montecarlo:            "<<R<<endl;
		//cout<<"RmObs contra el que compara: "<<RmObs<<endl;
		
		// Para aplicar ahora la forma frecuentista de hallar el intervalo
		// cuento cuantas veces, para este valor fijo de la edad, los valores
		// simulados de la likelihood caen a la izquierda del valor observado.
		if(R<=RmObs) nexp += 1;
		
		// Profiling ///////////////////////////////////////////////////
		// Calculo el valor que toma la likelihood con este R y si es mas grande que los anteriores lo actualizo
		// esto permite encontrar el maximo likelihood usando el mismo MonteCarlo
		double pRm=exp(-1/2*pow((mu-RmObs)/sigma,2))/(pow(2*3.141592654,0.5)*sigma);				
		double pRs=exp(-1/2*pow((Rs-RsObs)/sigmaRs,2))/(pow(2*3.141592654,0.5)*sigmaRs);				
		double pRf=exp(-1/2*pow((Rf-RfObs)/sigmaRf,2))/(pow(2*3.141592654,0.5)*sigmaRf);				
		double likelihood = pRm*pRs*pRf;
		if (LL[k]<likelihood){LL[k]=likelihood;}
			
		// Lleno un histograma de R que luego representara su verosimilitud L(R/Edad)
		// Con esos histogramas se hara una inversion que permitira hallar  L(Edad/R)
		// Esto permitira obtener el intervalo por el enfoque completamente Bayesiano
		verosimilitud -> Fill(R);
	  
	  } // Finaliza Loop MonteCarlo
	  
	  
	  // Normalizo el histograma de los R obtenidos para esta Edad[k]
		double norm=1; //Normalizacion    
		Double_t scale = norm/verosimilitud->Integral();
		verosimilitud->Scale(scale); //L(R/Edad) normalizada
		

        alpha[k] = (nexp/nToys)*100; // esto es la significancia, o bien(1- CL)/2
        beta[k]=100-alpha[k];

        //cout.precision(3);      // Usa precision decimal
        //cout<<fixed;            // con exactamente cuatro decimales
        //cout<<"Para una edad de "<<Edad[k]<<" anios  "<<" (1-CL)/2 es "<<alpha[k]<<" % y"<<" (1-CL)/2 es "<<beta[k]<<" %"<<endl;
   
        // Cuando alpha = 16%, voy a haber encontrado la edad maxima del intervalo
        // la que deja un 16 % de cola  a la izquierda del RmObs
		// Cuando alpha = 84 % voy a haber encontrado la edad minima del intervalo
			
		  if (alpha[k] > 15.9 && entre1){
			cout<<"El limite inferior del intervalo es:      "<<Edad[k]<<" anios"<<endl;
			entre1 = false;
		  }
		  
		  if (alpha[k] > 49.9 && lodijeh){
			cout<<"La edad de la muestra por el metodo es:   "<<Edad[k]<<" anios"<<endl;
			lodijeh = false;
		  }
		  
		  if (alpha[k] > 83.9 && entre2){
			cout<<"El limite superior del intervalo es:      "<<Edad[k]<<" anios"<<endl;
			entre2 = false;
		  }				
			
        // Guardo la probabilidad de caer en cada uno de los i bines del histograma verosimilitud
        // obtenidos para cada valor de la Edad[k]
		for(int i=1; i<nbines+1; i++) {
			R_bin[i][k]=verosimilitud -> GetBinContent(i); 
        }

       
   	c2->Clear();
	//c2->SetGrid();
	//gStyle->SetOptStat(0);
	//gPad->DrawFrame(0.5,0.85,3.5,1.1,"Summary of fit results;;Fit Bias");
					   
			verosimilitud -> SetLineColor(kRed +3);         
			verosimilitud -> SetMarkerColor(kRed +3);       
			verosimilitud -> SetMarkerSize(1.0);            
			verosimilitud -> SetMarkerStyle(24);            
			//verosimilitud->SetMinimum(0.000001);
			//verosimilitud->SetMaximum(yMaxim);
			verosimilitud ->Draw("HIST E1");
			//clusters->SaveAs( "./figures/Dist_n_exp.png"); 
			//getchar();

}

////////////////////////////////////////////////////////////////////////
cout<<endl;
cout<<"Los resultados para el calculo Frecuentista puro (profiling) son los siguientes:"<<endl;

double area = 0;	  

	for (int k=0;k<steps;++k){area +=LL[k]*paso;}
			
		cout<<"El area de LL dio: "<<area<<endl;
			
		for (int k=0;k<steps;++k){
			integral+=LL[k]*paso/area;
						
			  if (integral > 0.159 && pentre1){
				cout<<"El limite inferior del intervalo es:      "<<Edad[k]<<" anios"<<endl;
				pentre1 = false;
			  }
			  
			  if (integral > 0.499 && plodije){
				cout<<"La edad de la muestra por este metodo es: "<<Edad[k]<<" anios"<<endl;
				plodije = false;
			  }
			  
			  if (integral > 0.839){
				cout<<"El limite superior del intervalo es:      "<<Edad[k]<<" anios"<<endl;
				break;
			  }	
	}
	  
////////////////////////////////////////////////////////////////////////
cout<<endl;
cout<<"Los resultados para el calculo Bayesiano puro son los siguientes:"<<endl;

double Normalizacion = 0;

for (int k=0;k<steps;++k){
	Normalizacion+=R_bin[1+(nbines-1)/2][k];
}

//cout<< " La normalizacion resulto igual a: "<<Normalizacion<<endl;

		for (int k=0;k<steps;++k){ 

			  //1+(nbines-1)/2 corresponde al bin central, que deberia contener a Robs por construccion  
			  Proba+=R_bin[1+(nbines-1)/2][k]/Normalizacion;
			  
			  if (Proba > 0.159 && entre){
				cout<<"El limite inferior del intervalo es:      "<<Edad[k]<<" anios"<<endl;
				entre = false;
			  }
			  
			  if (Proba > 0.499 && lodije){
				cout<<"La edad de la muestra por este metodo es: "<<Edad[k]<<" anios"<<endl;
				lodije = false;
			  }
			  
			  if (Proba > 0.839){
				cout<<"El limite superior del intervalo es:      "<<Edad[k]<<" anios"<<endl;
				break;
			  }
		}

// Plot histogramas /////////////////////////////////////////////////////

    TCanvas *c1;
    c1 = new TCanvas("canvas1","canvas1",800,800);
    c1->Clear();
    c1->SetGrid();

    TMultiGraph *mg = new TMultiGraph();
    //mg->SetTitle("Integral a izquierda (ROJO) y a derecha (AZUL) hasta el valor observado; Valor teorico de la Edad; Significancia %");
        
            TGraph *gr1 = new TGraph(steps,Edad,alpha);
            //gr1->SetTitle("Confidence level as a function of s; 1-CL ; Signal");
			//gr1->SetMarkerColor(kRed);
			//gr1->GetXaxis()->SetTitle("UL for the signal, cosmic UHE neutrino flux");
			//gr1->GetYaxis()->SetTitle("(1 - Confidence Level)/2 %");
			gr1->GetXaxis()->SetLimits(Edadmin,Edadmax);
			gr1->GetYaxis()->SetRangeUser(0.0,100);
			gr1->SetLineColor(2);
			gr1->SetLineWidth(2);
			//gr1->Draw("AP*");
			
			TGraph *gr2 = new TGraph(steps,Edad,beta);
            //gr2->SetTitle("Confidence level as a function of s; 1-CL ; Signal");
			//gr2->SetMarkerColor(kBlue);
			//gr2->GetXaxis()->SetTitle("UL for the signal, cosmic UHE neutrino flux");
			//gr2->GetYaxis()->SetTitle("(1 - Confidence Level)/2 %");
			gr2->GetXaxis()->SetLimits(Edadmin,Edadmax);
			gr2->GetYaxis()->SetRangeUser(0.0,100);
			gr2->SetLineColor(4);
			gr2->SetLineWidth(2);
			//gr2->Draw("AP*");
								
			auto l1 = new TLine(0,16,Edadmax,16);
			l1->SetLineStyle(kDashed);
			
			auto l2 = new TLine(0,84,Edadmax,84);
			l2->SetLineStyle(kDashed);
					
			//mg->Add(gr1); No encontre porque grafica mal alpha ...
			mg->Add(gr2); 
			mg->Draw("AL*");
			l1->Draw("same");
			l2->Draw("same");

}
