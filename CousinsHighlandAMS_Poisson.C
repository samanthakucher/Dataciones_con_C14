void CousinsHighlandAMS_Poisson(){

#include <iostream>

/*
Este codigo no tiene en cuenta el error en 12C

Considera que los conteos de 14C son poissonianos, tanto muestra, fondo como estandar.

El codigo calcula la incerteza en la edad de una muestra medida por AMS 
usando cuatro metodos basados en la verosimilitud.

- Bayesiano. Aqui los nuisance parameters se eliminan por integracion hecha por Monte Carlo.
  Luego, se invierte la likelihood para finalmente buscar el rango que integra una probabilidad del 68 % centrada

- Frecuentista - Hace profiling para eliminar los nuisance parameters
  Para encontrar las relaciones isotopicas de fondo y estandar que maximizan la verosimilitud
  usa los ensayos hechos Montecarlo del caso anterior realizados para integrar los nuisance parameters
  
- Hibrido1: Usa la likelihood que quedo luego de eliminar los nuisance parameters por integracion bayesiana
  y el metodo del cinturon frecuentista para calcular el itervalo

- Hibrido2: Usa la los conteos de fondo y estandar obtenidos del profiling
  y luego cuenta cuantos eventos caen en el bin central para invertir y obtener la verosimilitud de la edad
  Luego integra esta verosimilitud buscando el intervalo bayesiano.
   
*/

// A continuacion se indican los resultados experimentales
// Numero de eventos TOTALES de 14C y 12C
int Nm=22;
double_t M14C=4673;
double_t M12C=1.52E16;

int Ns=8;
double_t S14C=1424;
double_t S12C=2.05E15;

int Nf=9;
double_t F14C=330;
double_t F12C=3.7E15; 

// Relaciones isotopicas
double RmObs = M14C/M12C;
double RsObs = S14C/S12C;
double RfObs = F14C/F12C;

// Errores relativos
double ErrRmObs = pow(M14C,0.5)/M14C;
double ErrRsObs = pow(S14C,0.5)/S14C;
double ErrRfObs = pow(F14C,0.5)/F14C;

//Errores absolutos
double sigmaRm = pow(M14C,0.5)/M12C;;
double sigmaRs = pow(S14C,0.5)/S12C;
double sigmaRf = pow(F14C,0.5)/F12C;;
		
std::cout << std::scientific;

cout<<endl;
cout<<"Este calculo se realizo para:"<<endl;
cout<<"Relacion isotopica en la muestra          = "<<RmObs <<"+/-"<<sigmaRm<<endl;
cout<<"Relacion isotopica en el estandar         = "<<RsObs <<"+/-"<<sigmaRs<<endl;
cout<<"Relacion isotopica en el fondo            = "<<RfObs <<"+/-"<<sigmaRf<<endl;

double T = 5730/0.693147; // vida media del 14C
double fecha = -T*log((RmObs-RfObs)/(RsObs-RfObs));

cout.precision(0);      // Usa precision decimal
cout<<fixed;            // con exactamente cero decimales
cout<<"La edad calculada con estos parametros es = "<< fecha  <<" anios"<<endl<<endl;

// Edades maximas y minimas a considerar al hacer el barrido en este parametro
int rango = 4000;
// El resultado depende del rango y no deberia ...
double_t Edadmin=fecha-rango;
double_t Edadmax=fecha+rango;

int steps = rango/5;
int nbines = rango/30+1;
int nToys = 10000;// numero de toys

double integral = 0;
double paso = 0;

double R = 0;
double RM14C = 0;
double RS14C = 0;
double RF14C = 0;

double Edad[steps];
double alpha[steps];
double beta[steps];
double LL[steps];
double LBayesiana[steps];
double entradas[steps];
double RS14C_hat[steps];
double RF14C_hat[steps];

for (int k=0;k<steps;++k){
	entradas[k]=0;
	Edad[k]=0;
	alpha[k]=0;
	beta[k]=0;
	LL[k]=0;
	LBayesiana[k]=0;
	RS14C_hat[k]=0;
	RF14C_hat[k]=0;
}

double R_bin[nbines][steps];
double K_bin[nbines][steps];

double Norma=0;
double Proba=0;

bool entre = true;
bool entre1 = true;
bool entre2 = true;
bool lodije = true;
bool lodijeh = true;
bool pentre1 = true;
bool plodije = true;
bool limiteinferior = true;
bool central = true;

//TCanvas *c2;
//c2 = new TCanvas("canvas2","canvas2",800,800);

TH1D * verosimilitud  =  new TH1D("Verosimilitud", "Verosimilitud", nbines, RmObs-20*sigmaRm, RmObs+20*sigmaRm);
TH1D * verosimilitudp  =  new TH1D("Verosimilitud", "Verosimilitud", nbines, RmObs-20*sigmaRm, RmObs+20*sigmaRm);

cout<<endl;
cout<<"Los resultados para el calculo Bayesiano-Frecuentista son los siguientes:"<<endl;
		
// cada paso corresponde a un valor diferente de la Edad, y por lo tanto,
// a un valor diferente de la relacion isotopica esperada en la muestra		
for (int k=0;k<steps;++k){ 
  
  paso = 1./(steps-1)*(Edadmax-Edadmin);
    
  Edad[k]=Edadmin+k*paso;
  //std::cout << std::scientific;
  //cout<<"Edad en este paso: "<<Edad[k]<<endl;
  
  // Rmesp_hat corresponde al caso donde Rs y Rf no tienen incerteza.
  double Rmesp_hat = RfObs+(RsObs-RfObs)*exp(-Edad[k]/T);
  //std::cout << std::scientific;
  //cout<<"Rmesp calculado sin errores en Rf y Rs: "<<Rmesp_hat<<endl;
	  
  // Aca empieza el loop para la integracion Bayesiana por MonteCarlo 
  // realizada sobre los nuisance parameters o errores sistematicos
  
  double nexp = 0;  // contara eventos donde Rm<=RmObs 
  
  TRandom *rnd=new TRandom3();
  
  LL[k]=0; //Inicializo la Likelihood que usare para hacer profiling
  
	    for(int iToy=0; iToy<nToys; iToy++) {

					double Rs=0;
					   RS14C=rnd->Poisson(S14C);
					   Rs=RS14C/S12C;

					double Rf=0;
					   RF14C=rnd->Poisson(F14C);
					   Rf=RF14C/F12C;
				
					//Calcula la relacion isotopica esperada de la muestra
					//usando el par de valores obtenidos para Rs y Rf
					double Rmesp = Rf+(Rs-Rf)*exp(-Edad[k]/T);
					//std::cout << std::scientific;
					//cout<<"Rm del Montecarlo: "<<Rmesp<<endl;
					
					// Ahora se genera un valor Rmesp para la relacion isotopica medida Rm
					// Dicho valor se calcula usando la likelihood como distribucion de probabilidad
					// Esa likelihood sera una Poissoniana de la forma Poisson(Rmesp*M12C);
					
					// Uno podria preguntarse donde estan las otras dos distribuciones, la de Rf y Rs
					// que analiticamente hay que multiplicar por la de Rm para obtener la verosimilitud,
					// las mismas fueron incorporadas al darle peso al Rf y al Rs con que se calculo el Rmesp
					
					double sigma = ErrRmObs*Rmesp;
					// Para esta Edad y con los valores simulados para RF14C y RS14C se obtiene el siguiente mu
					double mu = Rmesp*M12C; // Valor esperado para M14C
							
					RM14C = rnd->Poisson(mu);
					R=RM14C/M12C;
					
					//std::cout << std::scientific;
					//cout<<"R del Montecarlo:            "<<R<<endl;
					//cout<<"RmObs contra el que compara: "<<RmObs<<endl;
					
					// Para aplicar ahora la forma frecuentista exacta (cinturon) de hallar el
					// intervalo, cuento cuantas veces, para este valor fijo de la edad, los valores
					// simulados de la likelihood caen a la izquierda del valor observado.
					if(R<=RmObs) nexp += 1;
					
					
					// Cuento cuantas veces aparece el RObs para este valor de Edad[k]
					// para luego armar su Distribucion de probabilidad Bayesiana 
					if ((RmObs - sigmaRm/2) < R && R <(RmObs + sigmaRm/2)) {
						entradas[k] += 1;
					// El resultado depende debilmente del rango ...
					}
					
					// Profiling ///////////////////////////////////////////////////
					// Calculo el valor que toma la likelihood LL con este mu.
					// Si LL es mas grande que los anteriores lo actualizo.
					// Asi encuentro el maximo likelihood usando el mismo MonteCarlo
					
					
					//double Rmesp = Rf+(Rs-Rf)*exp(-Edad[k]/T);
					//double mu = Rmesp*M12C;
					
					// Porque LL debe optimizarse evaluada en los valores observados
					// Entonces, la likelihood debe estar valuada en el M14C porque es el valor obervado
					// Pero los valores que van donde la k de Poisson estan fijos en los experimentales.
					
					// Las R delante de RS14C y RF14C significan RANDOM no RATIO.
					
					double pRm = exp(M14C*log(mu)-mu-TMath::LnGamma(M14C+1));
					double pRs = exp(S14C*log(RS14C)-RS14C-TMath::LnGamma(S14C+1));
					double pRf = exp(F14C*log(RF14C)-RF14C-TMath::LnGamma(F14C+1));
					
					double likelihood = pRm*pRs*pRf;
						if (LL[k]<likelihood) {
							LL[k]=likelihood;
							RS14C_hat[k]=RS14C;
							RF14C_hat[k]=RF14C;
						}
						
					// Lleno un histograma de R que normalizado representara su verosimilitud L(R/Edad)
					// Con esos histogramas se hara una inversion para obtener L(Edad/R)
					// Esto permitira hallar el intervalo de confianza con el enfoque Bayesiano puro
					verosimilitud -> Fill(R);
		
        }// Finaliza Loop MonteCarlo ////////////////////////////////////
        

			  /* 
			  Normalizo el histograma de los R obtenidos para esta Edad[k]
			  double norm=1; //Normalizacion    
			  Double_t scale = norm/verosimilitud->Integral();
			  verosimilitud->Scale(scale); //verosimilitud = L(R/Edad)
			  */

    		  //Guardo el numero de eventos en cada uno de los nbines del histograma verosimilitud
			  //obtenidos para este valor particular de la Edad[k]
			  for(int i=1; i<nbines+1; i++) {
					R_bin[i][k]=verosimilitud -> GetBinContent(i); 
			   } 
	  	
        alpha[k] = (nexp/nToys)*100; // esto es la significancia, o bien (1- CL)/2
        // Notar que conforme la Edad aumenta R disminuye y por lo tanto alpha aumenta
        beta[k]=100-alpha[k];

        //cout<<"Para una edad de "<<Edad[k]<<" anios  "<<" (1-CL)/2 es "<<alpha[k]<<" % y"<<" (1-CL)/2 es "<<beta[k]<<" %"<<endl;
   
        // Cuando alpha = 16%, voy a haber encontrado la edad maxima del intervalo
        // la que deja un 16 % de cola  a la izquierda del RmObs
		// Cuando alpha = 84 % voy a haber encontrado la edad minima del intervalo
		
		//cout<<"Edad "<<Edad[k]<<"  Alpha  "<<alpha[k]<<endl;
						
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
			
       
			//c2->Clear();				   
			//verosimilitud -> SetLineColor(kRed +3);         
			//verosimilitud -> SetMarkerColor(kRed +3);       
			//verosimilitud -> SetMarkerSize(1.0);            
			//verosimilitud -> SetMarkerStyle(24);            
			//verosimilitud -> SetMinimum(0.000001);
			//verosimilitud -> SetMaximum(yMaxim);
			//verosimilitud -> Draw("HIST E1");
			
			//cout<<"k: "<<k<<"  Entradas: "<<entradas[k]<<endl;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
cout<<endl;
cout<<"Los resultados para el calculo Frecuentista (profiling) con intervalo bayesiano son:"<<endl;

double area = 0;	  

	for (int k=0;k<steps;++k) area +=LL[k];
			
		//cout<<"El area de LL dio: "<<area<<endl;
			
		for (int k=0;k<steps;++k){
			integral+=LL[k]/area; //Normalizo por la integral sobre las edades.
						
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
////////////////////////////////////////////////////////////////////////
cout<<endl;
cout<<"Los resultados para el calculo Frecuentista puro (profiling + cinturon) son:"<<endl;

for (int k=0;k<steps;++k){ //Barre sobre sobre Edad[k] para armar las L(M14C/Edad[k])

	  // Para cada edad calculo el valor medio esperado para RM14C considerando los hat de fondo y estandar
      double mu_RM14C = (RF14C_hat[k]/F12C+(RS14C_hat[k]/S12C-RF14C_hat[k]/F12C)*exp(-Edad[k]/T))*M12C;

	  //Ahora lleno un histograma usando una Poisson centrada en este mu_RM14C y divido los resultados por M12C
	  TRandom *rnd=new TRandom3();
	  
	  for(int iToy=0; iToy<nToys; iToy++) {
				// Simulo conteos de 14C
				double KM14C = rnd->Poisson(mu_RM14C);
				double KR = KM14C/M12C; //Relacion isotopica
			    verosimilitudp -> Fill(KR);
       }
         
	  //Normalizo
	  double scale = 1/verosimilitudp->Integral();
	  verosimilitudp->Scale(scale); //verosimilitud normalizada
	  
	  double Proba2 = verosimilitudp->Integral(1,1+(nbines-1)/2);
	  
	  //cout<<"Proba2: "<<Proba2<<endl;
	  
	  		  // Ahora miro si Proba2 alcanzo al 16%			  
			  if (Proba2 > 0.159 && limiteinferior){
				cout<<"El limite inferior del intervalo es:      "<<Edad[k]<<" anios"<<endl;
				limiteinferior = false;
			  }
			  
			  // Ahora miro si Proba2 alcanzo al 50%
			  if (Proba2 > 0.499 && central){
				cout<<"La edad de la muestra por este metodo es: "<<Edad[k]<<" anios"<<endl;
				central = false;
			  }
			  
			  // Ahora miro si Proba2 alcanzo al 84%
			  if (Proba2 > 0.839){
				cout<<"El limite superior del intervalo es:      "<<Edad[k]<<" anios"<<endl;
				break;
			  }
}

////////////////////////////////////////////////////////////////////////	  
////////////////////////////////////////////////////////////////////////
cout<<endl;
cout<<"Los resultados para el calculo Bayesiano puro son los siguientes:"<<endl;

double Normalizacion = 0;
for (int k=0;k<steps;++k){ // Suma sobre todas las edades
	Normalizacion+=entradas[k];
}
		
		for (int k=0;k<steps;++k){ 
			  
			  Proba+=entradas[k]/Normalizacion;
			  			  
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
		
////////////////////////////////////////////////////////////////////////
// Plot histogramas /////////////////////////////////////////////////////
/*
	TCanvas *c3;
    c3 = new TCanvas("canvas3","canvas3",800,800);
    c3->Clear();
    c3->SetGrid();

	TGraph *gr5 = new TGraph(steps,Edad,entradas);
	gr5->Draw("AP*");
*/
////////////////////////////////////////////////////////////////////////

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
