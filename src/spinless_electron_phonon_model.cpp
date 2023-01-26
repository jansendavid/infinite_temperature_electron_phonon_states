#include "itensor/all.h"
#include "Holstein_spinless.hpp"
#include "makebasis.hpp"


using namespace itensor;


int
main(int argc, char *argv[])
    {
      int L=6;
      int M=10;
      int N=2*L;
      int Maxd=150;    
      int Mind=0;    
      int n_e=int(L/2);
      double cutoff=1e-09;
      
        auto sites = Holstein_purified(N,{"ConserveNf=",true,
                              "ConserveNb=",false,
       				"MaxOcc",M});

        auto psi=makeBasis_phonons(sites, M);
	print_bonddims(psi);
	auto A=psi;

	for(int j=1; j<=n_e;j++)
	  {	std::vector<MPS> states;
	for(int i=1; i<=N; i+=2){
	  auto psi2=A;
   	  applyCdag(psi2, i, sites);
  applyCdag(psi2, i+1, sites);
  states.push_back(psi2);

}
	A = sum(states,{"MaxDim",Maxd,"MinDim", Mind,"Cutoff",cutoff,"ShowEigs", false});
    }

	psi=A;
       print_bonddims(psi);
	//       auto ampo=makeHolstHam_FT_halffill(sites, t0, gamma, omega);

	//auto H=toMPO(ampo);
 




   
   for(int b = 1; b <= N; b+=2)
        {
   	  psi.position(b);

   	     auto bb = AutoMPO(sites);
   	    
  
   	     bb += 1,"X",b;
   	 auto gg = AutoMPO(sites);
  
   		        gg += 1,"Nph",b ;
   		      auto GG = toMPO(gg);
   		      auto en2 = innerC(psi,GG,psi)/innerC(psi, psi);
   		         printfln("Nph %d %.12f",b,en2);
   			 auto gg2 = AutoMPO(sites);
  
   		        gg2 += 1,"N",b ;
   		      auto GG2 = toMPO(gg2);
   		      auto en22 = innerC(psi,GG2,psi)/innerC(psi, psi);
   		         printfln("Ne %d %.12f",b,en22);
		        

   
   }  //
		
		 
    return 0;
    
    }

