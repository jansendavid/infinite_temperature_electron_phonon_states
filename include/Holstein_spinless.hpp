
//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#pragma once
#include "itensor/mps/siteset.h"
#include "itensor/all.h"
#include"extra.hpp"
namespace itensor {





class Holstein_spinlessSite_up;

using Holstein_spinless_up = BasicSiteSet<Holstein_spinlessSite_up>;

class Holstein_spinlessSite_down;

using Holstein_down = BasicSiteSet<Holstein_spinlessSite_down>;

using Holstein_spinless = BasicSiteSet<Holstein_spinlessSite>;  


  // holstein up
class Holstein_spinlessSite_up
    {
      Index s;


      
      std::vector<int> mvec;
      std::vector<Index> indexvec;
    public:

      Holstein_spinlessSite_up(Index I): s(I) { }

      Holstein_spinlessSite_up( Args const& args = Args::global() )
        {
        auto conserveQNs = args.getBool("ConserveQNs",true);
        // auto conserveNb = args.getBool("ConserveNb",conserveQNs);
        auto conserve_Nf = args.getBool("ConserveNf",conserveQNs);
	        auto oddevenupdown = args.getBool("OddEvenUpDown",false);
	auto tags = TagSet("Site,Hol");
        auto n = 1;
        if(args.defined("SiteNumber") )
            {
            n = args.getInt("SiteNumber");
            tags.addTags("n="+str(n));
            }
	auto diffMaxOcc = args.getBool("DiffMaxOcc",false);
	auto maxOcc{0};
	// To use this feature, minor adjustments are needed in the args file to accept a vector with ints
	// if(diffMaxOcc)
	//   {
	//     auto occVec=args.getVecInt("MaxOccVec");
	//     maxOcc= occVec[n-1];
	//   }
	else{
        maxOcc = args.getInt("MaxOcc",1);
	}
	        if(not conserveQNs)
            {
	      s = Index(2*(maxOcc+1),tags);
            }
        else if(not oddevenupdown) //usual case
            {

	    	      if( conserve_Nf) //usual case
{
		auto q_emp = QN({"Nf",0, -1},{"Sz",0});
		auto q_occ = QN({"Nf",1, -1},{"Sz",-1});
            s = Index(q_emp,(maxOcc+1),
                      q_occ,(maxOcc+1),Out,tags);
	      }
	      else{
		auto q_emp = QN({"Pf",0,-2});
		auto q_occ = QN({"Pf",1,-2});
		s = Index(q_emp,(maxOcc+1),
		q_occ,(maxOcc+1),Out,tags);
}
}
        else
            {
            QN q_occ;
	    auto q_emp = QN({"Sz",0},{"Nf",0,-1});
            if(n%2==1) q_occ = QN({"Sz",+1},{"Nf",1,-1});
            else       q_occ = QN({"Sz",-1},{"Nf",1,-1});
            s = Index(q_emp,(maxOcc+1),
                      q_occ,(maxOcc+1),Out,tags);

            }
	}
      // not applicabale yet
        // if(conserveQNs)
        //     {
        //     if(conserveNb)
        //         {
        //         auto qints = Index::qnstorage(1+maxOcc);
        //         for(int n : range(1+maxOcc)) 
        //             {
        //             qints[n] = QNInt(QN({"Nb",n}),1);
        //             }
        //         s = Index(std::move(qints),tags);
        //         }
        //     else
        //         {
        //         s = Index(QN(),1+maxOcc,tags);
        //         }
        //     }
        // else
        //     {
        //     if(conserveNb) Error("ConserveNb cannot be true when ConserveQNs=false");
        //     s = Index(1+maxOcc,tags);
        //     }
        // }


    Index
    index() const { return s; }

    IndexVal
    state(std::string const& state)
        {
	  auto d = static_cast<int>(dim(s)/2);
        if(state == "Emp" || state == "0") 
            {
            return s(1);
            }
        else 
        if(state == "Occ" || state == "1") 
            {
            return s(d+1);
            }
	 if(state == "EmpPh") 
            {
            return s(2);
            }
        else 
        if(state == "OccPh") 
            {
            return s(d+2);
            }
	 	 if(state == "EmpPh2") 
            {
            return s(3);
            }
        else 
        if(state == "OccPh2") 
            {
            return s(d+3);
            }
        else
            {
            Error("State " + state + " not recognized");
            }
        return IndexVal{};
        }

	ITensor
	op(std::string const& opname,
	   Args const& args) const
        {
        auto sP = prime(s);
	auto d = static_cast<int>(dim(s)/2);
	auto Op = ITensor(dag(s),sP);
	

        if(opname == "N" || opname == "n")
            {

	      for (auto n :range1(d))
		{
		
		  Op.set(s=(d+n),sP=(d+n),1);
	      }
            }
	else
	  if(opname == "Nph")
	    {


	       
	      for(auto n :range1(d))
		{
	
		 Op.set(s=n, sP=n, n-1);
		 Op.set(s=(n+d), sP=(n+d), n-1);
 	      }
	    }
	  else
	     if(opname == "X" )
            {


            for(auto n : range1(d-1))
                {
                Op.set(s=n,sP=1+n,std::sqrt(n));
		Op.set(s=(d+n),sP=(d+1+n),std::sqrt(n));
		Op.set(s=1+n,sP=n,std::sqrt(n));
		Op.set(s=(d+n+1),sP=(d+n),std::sqrt(n));
                }
	      if(d-1==0)
	       {
		 Op.set(s=(1),sP=(1),0);
	       }

            }
	else
	  if(opname == "Bdag" )
            {


            for(auto n : range1(d-1))
                {
                Op.set(s=n,sP=1+n,std::sqrt(n));
		Op.set(s=(d+n),sP=(d+1+n),std::sqrt(n));
                }
	      if(d-1==0)
	       {
		 Op.set(s=(1),sP=(1),0);
	       }

            }
	else
	  if(opname == "B" )
            {

            for(auto n : range1(d-1))
                {
                Op.set(s=1+n,sP=n,std::sqrt(n));
		Op.set(s=(d+n+1),sP=(d+n),std::sqrt(n));
                }
	        if(d-1==0)
	       {
		 Op.set(s=(1),sP=(1),0);
	       }
 	      
            }
	   else
	  if(opname == "Bdagnorm" )
            {

	      // normalized operator to implement projected purification PP
            for(auto n : range1(d-1))
                {
                Op.set(s=n,sP=1+n,1);
		Op.set(s=(d+n),sP=(d+1+n),1);
                }
	      if(d-1==0)
	       {
		 Op.set(s=(1),sP=(1),0);
	       }

            }
	else
	  if(opname == "Bnorm" )
            {
	      // normalized operator to implement projected purification PP
            for(auto n : range1(d-1))
                {
                Op.set(s=1+n,sP=n,1);
		Op.set(s=(d+n+1),sP=(d+n),1);
                }
	      if(d-1==0)
	       {
		 Op.set(s=(1),sP=(1),0);
	       }
	      
 	      
            }
	  else
	    	  if(opname == "NX" )
            {


            for(auto n : range1(d-1))
                {
		  
		Op.set(s=(d+n),sP=(d+1+n),std::sqrt(n));
		Op.set(s=(d+n+1),sP=(d+n),std::sqrt(n));
		
                }
	     if(d-1==0)
	       {

		 Op.set(s=(1),sP=(1),0);
	       }

            }
	else
		  if(opname == "NBdag" )
            {


            for(auto n : range1(d-1))
                {
		  
		Op.set(s=(d+n),sP=(d+1+n),std::sqrt(n));
                }
	     if(d-1==0)
	       {

		 Op.set(s=(1),sP=(1),0);
	       }

            }
	else
	  if(opname == "NB" )
            {

            for(auto n : range1(d-1))
                {
                
		Op.set(s=(d+n+1),sP=(d+n),std::sqrt(n));
                }
	     if(d-1==0)
	       {
		 Op.set(s=(1),sP=(1),0);
	       }
 	      
            }
        else
        if(opname == "C")
            {
	      
 
		for(auto n : range1(d))
		  {
 
		 Op.set(s=(d+n), sP=(n), 1);
	      }

            }
        else
        if(opname == "Cdag")
            {

		for(auto n : range1(d))
                {

		Op.set( s=n, sP=(d+n), 1);
	      }

            }
        else
        if(opname == "A")
	  {

		for(auto n : range1(d))
                {
	
		 Op.set(s=(d+n), sP=(n), 1);
	      }

	    
            }
        else
        if(opname == "Adag")
            {

		for(auto n : range1(d))
                {
	
		Op.set( s=n, sP=(d+n), 1);
	      }	      

            }
        else
        if(opname == "F" || opname == "FermiPhase")
            {

	      		for(auto n : range1(d))
                {

            Op.set(s=n,sP=n,1);
            Op.set(s=(d+n), sP=(d+n),-1);
	      }

            }
        else
        if(opname == "projEmp")
            {
	      	for(auto n : range1(d))
		  {
            Op.set(s=n,sP=n,1);
		  }
            }
        else
        if(opname == "projOcc")
            {
	      	for(auto n : range1(d))
		  {
		    Op.set(s=(d+n),sP=(d+n),1); 
            }
	    }
        else
	  if(opname == "maxPhProj")
            {
	      	// for(auto n : range1(d))
		//   {
            Op.set(s=d,sP=d,1);
	    Op.set(s=2*d,sP=2*d,1);
	    //  }
            }
        else
            {
            Error("Operator \"" + opname + "\" name not recognized");
            }

        return Op;
        }

      HolsteinSite_up(int n, Args const& args = Args::global())
      {
        *this = HolsteinSite_up({args,"SiteNumber=",n});
      }
};



  // down
  // holstein 
class Holstein_spinlessSite_down
    {
      Index s;


      
      std::vector<int> mvec;
      std::vector<Index> indexvec;
    public:

      Holstein_spinlessSite_down(Index I): s(I) { }

      //HolsteinSite(IQIndex I, int Mm=4) : M(Mm), s(I)  { } // not quite safe, how to control f the index contains right phonon number

      Holstein_spinlessSite_down( Args const& args = Args::global() )
        {
        auto conserveQNs = args.getBool("ConserveQNs",true);
        // auto conserveNb = args.getBool("ConserveNb",conserveQNs);
        auto conserve_Nf = args.getBool("ConserveNf",conserveQNs);
	auto oddevenupdown = args.getBool("OddEvenUpDown",false);
	auto tags = TagSet("Site,Hol");
        auto n = 1;
        if(args.defined("SiteNumber") )
            {
            n = args.getInt("SiteNumber");
            tags.addTags("n="+str(n));
            }
	auto diffMaxOcc = args.getBool("DiffMaxOcc",false);
	auto maxOcc{0};
	// To use this feature, minor adjustments are needed in the args file to accept a vector with ints
	// if(diffMaxOcc)
	//   {
	//     auto occVec=args.getVecInt("MaxOccVec");
	//     maxOcc= occVec[n-1];
	//   }
	else{
        maxOcc = args.getInt("MaxOcc",1);
	}
	        if(not conserveQNs)
            {
	      s = Index(2*(maxOcc+1),tags);
            }
        else if(not oddevenupdown) //usual case
            {

	      if( conserve_Nf) //usual case
{
		auto q_emp = QN({"Nf",0, -1},{"Sz",0});
		auto q_occ = QN({"Nf",1, -1},{"Sz",-1});
            s = Index(q_emp,(maxOcc+1),
                      q_occ,(maxOcc+1),Out,tags);
	      }
	      else{
		auto q_emp = QN({"Pf",0,-2});
		auto q_occ = QN({"Pf",1,-2});
s = Index(q_emp,(maxOcc+1),
                      q_occ,(maxOcc+1),Out,tags);
            }
        else
            {
            QN q_occ;
	    auto q_emp = QN({"Sz",0},{"Nf",0,-1});
            if(n%2==1) q_occ = QN({"Sz",+1},{"Nf",1,-1});
            else       q_occ = QN({"Sz",-1},{"Nf",1,-1});
            s = Index(q_emp,(maxOcc+1),
                      q_occ,(maxOcc+1),Out,tags);

            }
	}
        // if(conserveQNs)
        //     {
        //     if(conserveNb)
        //         {
        //         auto qints = Index::qnstorage(1+maxOcc);
        //         for(int n : range(1+maxOcc)) 
        //             {
        //             qints[n] = QNInt(QN({"Nb",n}),1);
        //             }
        //         s = Index(std::move(qints),tags);
        //         }
        //     else
        //         {
        //         s = Index(QN(),1+maxOcc,tags);
        //         }
        //     }
        // else
        //     {
        //     if(conserveNb) Error("ConserveNb cannot be true when ConserveQNs=false");
        //     s = Index(1+maxOcc,tags);
        //     }
        // }


    Index
    index() const { return s; }

    IndexVal
    state(std::string const& state)
        {
	  auto d = static_cast<int>(dim(s)/2);
        if(state == "Emp" || state == "0") 
            {
            return s(1);
            }
        else 
        if(state == "Occ" || state == "1") 
            {
            return s(d+1);
            }
	 if(state == "EmpPh") 
            {
            return s(2);
            }
        else 
        if(state == "OccPh") 
            {
            return s(d+2);
            }
	 	 if(state == "EmpPh2") 
            {
            return s(3);
            }
        else 
        if(state == "OccPh2") 
            {
            return s(d+3);
            }
        else
            {
            Error("State " + state + " not recognized");
            }
        return IndexVal{};
        }

	ITensor
	op(std::string const& opname,
	   Args const& args) const
        {
        auto sP = prime(s);
	auto d = static_cast<int>(dim(s)/2);
	auto Op = ITensor(dag(s),sP);
	

        if(opname == "N" || opname == "n")
            {

	      for (auto n :range1(d))
		{
		
		  Op.set(s=(d+n),sP=(d+n),1);
	      }
            }
	else
	  if(opname == "Nph")
	    {


	       
	      for(auto n :range1(d))
		{
	
		 Op.set(s=n, sP=n, n-1);
		 Op.set(s=(n+d), sP=(n+d), n-1);
 	      }
	    }
	  else
	     if(opname == "X" )
            {


            for(auto n : range1(d-1))
                {
                Op.set(s=n,sP=1+n,std::sqrt(n));
		Op.set(s=(d+n),sP=(d+1+n),std::sqrt(n));
		Op.set(s=1+n,sP=n,std::sqrt(n));
		Op.set(s=(d+n+1),sP=(d+n),std::sqrt(n));
                }
	      if(d-1==0)
	       {
		 Op.set(s=(1),sP=(1),0);
	       }

            }
	else
	  if(opname == "Bdag" )
            {


            for(auto n : range1(d-1))
                {
                Op.set(s=n,sP=1+n,std::sqrt(n));
		Op.set(s=(d+n),sP=(d+1+n),std::sqrt(n));
                }
	      if(d-1==0)
	       {
		 Op.set(s=(1),sP=(1),0);
	       }

            }
	else
	  if(opname == "B" )
            {

            for(auto n : range1(d-1))
                {
                Op.set(s=1+n,sP=n,std::sqrt(n));
		Op.set(s=(d+n+1),sP=(d+n),std::sqrt(n));
                }
	        if(d-1==0)
	       {
		 Op.set(s=(1),sP=(1),0);
	       }
 	      
            }
	   else
	  if(opname == "Bdagnorm" )
            {


            for(auto n : range1(d-1))
                {
                Op.set(s=n,sP=1+n,1);
		Op.set(s=(d+n),sP=(d+1+n),1);
                }
	      if(d-1==0)
	       {
		 Op.set(s=(1),sP=(1),0);
	       }

            }
	else
	  if(opname == "Bnorm" )
            {

            for(auto n : range1(d-1))
                {
                Op.set(s=1+n,sP=n,1);
		Op.set(s=(d+n+1),sP=(d+n),1);
                }
	      if(d-1==0)
	       {
		 Op.set(s=(1),sP=(1),0);
	       }
	      
 	      
            }
	  else
	    	  if(opname == "NX" )
            {


            for(auto n : range1(d-1))
                {
		  
		Op.set(s=(d+n),sP=(d+1+n),std::sqrt(n));
		Op.set(s=(d+n+1),sP=(d+n),std::sqrt(n));
		
                }
	     if(d-1==0)
	       {

		 Op.set(s=(1),sP=(1),0);
	       }

            }
	else
		  if(opname == "NBdag" )
            {


            for(auto n : range1(d-1))
                {
		  
		Op.set(s=(d+n),sP=(d+1+n),std::sqrt(n));
                }
	     if(d-1==0)
	       {

		 Op.set(s=(1),sP=(1),0);
	       }

            }
	else
	  if(opname == "NB" )
            {

            for(auto n : range1(d-1))
                {
                
		Op.set(s=(d+n+1),sP=(d+n),std::sqrt(n));
                }
	     if(d-1==0)
	       {
		 Op.set(s=(1),sP=(1),0);
	       }
 	      
            }
        else
        if(opname == "C")
            {
	      
 
		for(auto n : range1(d))
		  {
 
		 Op.set(s=(d+n), sP=(n), 1);
	      }

            }
        else
        if(opname == "Cdag")
            {

		for(auto n : range1(d))
                {

		Op.set( s=n, sP=(d+n), 1);
	      }

            }
        else
        if(opname == "A")
	  {

		for(auto n : range1(d))
                {
	
		 Op.set(s=(d+n), sP=(n), 1);
	      }

	    
            }
        else
        if(opname == "Adag")
            {

		for(auto n : range1(d))
                {
	
		Op.set( s=n, sP=(d+n), 1);
	      }	      

            }
        else
        if(opname == "F" || opname == "FermiPhase")
            {

	      		for(auto n : range1(d))
                {

            Op.set(s=n,sP=n,1);
            Op.set(s=(d+n), sP=(d+n),-1);
	      }

            }
        else
        if(opname == "projEmp")
            {
	      	for(auto n : range1(d))
		  {
            Op.set(s=n,sP=n,1);
		  }
            }
        else
        if(opname == "projOcc")
            {
	      	for(auto n : range1(d))
		  {
		    Op.set(s=(d+n),sP=(d+n),1); 
            }
	    }
        else
	  if(opname == "maxPhProj")
            {
	      	// for(auto n : range1(d))
		//   {
            Op.set(s=d,sP=d,1);
	    Op.set(s=2*d,sP=2*d,1);
	    //  }
            }
        else
            {
            Error("Operator \"" + opname + "\" name not recognized");
            }

        return Op;
        }

      HolsteinSite_down(int n, Args const& args = Args::global())
      {
        *this = HolsteinSite_down({args,"SiteNumber=",n});
      }
};





 using  Holstein_exp = MixedSiteSet<Holstein_spinlessSite_down,Holstein_spinlessSite_up>;






} //namespace itensor

