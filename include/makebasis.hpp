
#ifndef MAKEBASIS_H
#define MAKEBASIS_H
#include "itensor/all.h"
#include <string>

namespace itensor{
  // helper function to see bond dimensios
  void print_bonddims(MPS psi)
  {  for(int j=1; j<length(psi); j++)
    {
      auto ind_K=rightLinkIndex(psi,j);



      std::cout<< "j " << j << " bond spdim "<<static_cast<double>(dim(ind_K))<<std::endl;

    }

    return;}

  // helper functions to apply fermionic gates  
template<typename SiteSet>
  void applyCdag(MPS& psi, int i, SiteSet& sites)
  {
    psi.position(i);
    psi.Aref(i)*= op(sites,"Adag",i);
    psi.Aref(i).noPrime("Site");
    for(int k = i-1; k >0; --k)
      {
	psi.position(k);
	psi.Aref(k) *= op(sites,"F",k);
	psi.Aref(k).noPrime("Site");
      }
  
  }

 
  template<typename SiteSet>
  void applyC(MPS& psi, int i, SiteSet& sites)
  {
    for(int k = 1; k <=i-1; ++k)
      {
	psi.position(k);
	psi.Aref(k) *= op(sites,"F",k);
	psi.Aref(k).noPrime("Site");
      }
    psi.position(i);
    psi.Aref(i) *= op(sites,"A",i);
    psi.Aref(i).noPrime("Site");
  
  }

  template<typename Sites, typename MPS, typename Index>
  void  make_loc(const Sites& sites, MPS& psi, Index& I, int n, int offset=0)
    {
      auto s1 = sites(n);
    	 auto J = Index(QN(qn(s1, 1)),1, "U,Link");
         auto s2 =(sites(n+1));
         auto wf =ITensor( s1,s2,  J, dag(I));
	 auto d=blocksize(s1, 1);
    	for(int i=1; i<=d; i++)
    	  {
	    wf.set(s1(i+offset),s2(i+offset),J(1), dag(I)(1),  1);
      }
	psi.Aref(n) = ITensor(s1, dag(I));
    	psi.Aref(n+1) = ITensor(s2, J);
	ITensor D;
    	svd(wf,psi.Aref(n),D,psi.Aref(n+1));
    	psi.Aref(n) *= D;
    	 psi.Aref(n+1)=noPrime(psi.Aref(n+1));
	 I=J;
    }
    template<typename Sites, typename MPS, typename Index>
    void  make_loc_last(const Sites& sites, MPS& psi, Index& I,int offset=0)
    {
      	 auto N=length(sites);
     	  auto s1 = sites(N-1);
         auto s2 =(sites(N));
	 auto wf = ITensor( s1,s2,dag(I));
auto d=blocksize(s1, 1);

	
    	for(int i=1; i<=d; i++)
    	  {
    	    wf.set(s1(i+offset),s2(i+offset), dag(I)(1), 1);
    	  }
    	psi.Aref(N-1) = ITensor(s1, dag(I));
    	psi.Aref(N) = ITensor(s2);
		ITensor D;
    	svd(wf,psi.Aref(N-1),D,psi.Aref(N));
    	psi.Aref(N-1) *= D;
    	 psi.Aref(N)=noPrime(psi.Aref(N));
    }
   template<typename Sites, typename MPS, typename Index>
   void  make_loc_first(const Sites& sites, MPS& psi, Index& I, int offset=0)
    {
      	 auto s1 = sites(1);
     auto s2 = (sites(2));
     auto wf = ITensor(s1,s2, I);
	auto d=blocksize(s1, 1);
	
	for(int i=1; i<=d; i++)
	  {
	     wf.set(s1(i+offset),s2(i+ offset),I(1), 1);


	  }
         ITensor D;

 
    	psi.Aref(1) = ITensor(s1);
    	psi.Aref(2) = ITensor(s2, I);
    	svd(wf,psi.Aref(1),D,psi.Aref(2));
    	psi.Aref(1) *= D;

   
     	 psi.Aref(2)=noPrime(psi.Aref(2));
    }
  /* function to make the bare infinite temperature phonon state*/

  template<typename Sites>
auto  makeBasis_phonons(const Sites& sites, int M)
 ->MPS  {
    
  int N=length(sites);
     auto psi = MPS(sites);
     auto s1 = sites(1);
     auto s2 = (sites(2));
    
     auto I = Index(QN(qn(s1, 1)),1, "U,Link");
 
     make_loc_first( sites, psi,I);
    for(int n = 3; n <= N-2; n += 2)
        {

	  make_loc( sites, psi,I, n);
	 
        }
    make_loc_last( sites, psi,I);
	

  return psi;
	   
  }
  /* function to make the infinite temperature polaron state phonons + pne electron*/
  template<typename Sites>
auto  makeBasis_one_electron(const Sites& sites, int M)
 ->MPS {
  int N=length(sites);
  // make site 1
     auto psi = MPS(sites);
     auto s1 = sites(1);
     auto s2 = (sites(2));    
     auto I = Index(QN(qn(s1, 1)),1, "U,Link");
     auto wf = ITensor(s1,s2, I);
     auto d=blocksize(s1, 1);
     std::vector<MPS> MPSs;

     make_loc_first(sites, psi, I, d);
     	 psi.Aref(2)=noPrime(psi.Aref(2));
    for(int n = 3; n <= N-2; n += 2)
        {
	  
    	  make_loc( sites, psi,I, n);	 
	 
        }
    	 
	  make_loc_last(sites, psi, I);
	  MPSs.push_back(psi);
	 for(int j=3; j<N-2; j+=2)
	   {
	       auto I = Index(QN(qn(s1, 1)),1, "U,Link");
	     auto s1 = sites(1);    
   
     auto psi2 = MPS(sites);
     make_loc_first(sites, psi2, I);
	      for(int n=3; n<N-2; n+=2)
	    {
	     if(n!=j)
	       {
	      make_loc( sites, psi2,I, n);
	       }
	     if(j==n)
	       {
		 make_loc( sites, psi2,I, n, d);
	       }
	    }
	     make_loc_last( sites, psi2,I);
	     MPSs.push_back(psi2);

	   }
	 // last site
	 I = Index(QN(qn(s1, 1)),1, "U,Link");
	 s1 = sites(1);
     auto psi3 = MPS(sites);
     make_loc_first(sites, psi3, I);
        for(int n = 3; n <= N-2; n += 2)
        {
	  
    	  make_loc( sites, psi3,I, n);	 
	 
        }
     	  
	make_loc_last( sites, psi3,I, d);
	  MPSs.push_back(psi3);

	  return sum(MPSs,{"MaxDim",500,"Cutoff",1E-15, "Normalize", false} );
	   
  }

}



#endif /* MAKEBASIS_H */
