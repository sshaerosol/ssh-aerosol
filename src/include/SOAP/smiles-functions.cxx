//!!-----------------------------------------------------------------------
//!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
//!!     SSH-aerosol is distributed under the GNU General Public License v3
//!!-----------------------------------------------------------------------

#include <fstream>
using namespace ssh_soap;

void get_branch(int &debug_mode, model_config &config, vector<species>& surrogate,
		Array<int, 1> &carbon_alcool, Array<int, 1> &carbon_near_alcool, Array<int, 1> &carbon_tail, Array<int, 1> &carbon_hydrogen, Array<int, 1> &carbon_nh, Array<int, 1> &carbon_taken,
		Array<int, 1> &carbon_arom, Array<int, 1> &linked_to_arom, Array<int, 1> &carbon_near_alcool_tmp, Array<int, 1> &carbon_near_alcool_save, Array<int, 1> &carbon_near_alcool_bef,
		Array<int, 1> &carbon_anhydre, Array<int, 1> &carbon_ester,
		Array<int, 2> &carbon_cycle, Array<int, 2> carbon_ether, Array<int, 2> carbon_ketone, Array<int, 2> &carbon_peroxide,
		int &arom, int &nitrate, int &nitrite, int &ester, int &hydroxyperoxide,
		int &iket, int &ipero, int &iether,
		int i, int &ngr, int sum_len_group, int ipar)
{
  int j,icycle,igr;

  if (surrogate[i].smile.substr(ipar,4)=="(OO)")
    hydroxyperoxide++;			
  if (surrogate[i].smile.substr(ipar,10)=="(ON(=O)=O)" or surrogate[i].smile.substr(ipar,14)=="(O[N+](=O)[O-]")
    nitrate++;
  if (surrogate[i].smile.substr(ipar,9)=="(N(=O)=O)" or surrogate[i].smile.substr(ipar,13)=="([N+](=O)[O-]")
    nitrite++;
			  
  if (surrogate[i].smile.substr(ipar,2)=="(G" and surrogate[i].smile.substr(ipar,3)!="(GK" and surrogate[i].smile.substr(ipar,3)!="(GL" and surrogate[i].smile.substr(ipar,5)!="(GC=O")
    {
      if (debug_mode==1) cout << "ester loc2" << endl;
      ester++;
    }

  if (surrogate[i].smile.substr(ipar,2)=="C(" and ipar==0)
    {
      if (arom==1 and carbon_arom(sum_len_group+2)==0)
	{
	  linked_to_arom(sum_len_group+2)=1;
	  linked_to_arom(sum_len_group)=1;
	}
      else if (arom==0 and carbon_arom(sum_len_group+2)==1)
	{
	  linked_to_arom(sum_len_group+2)=1;
	  linked_to_arom(sum_len_group)=1;
	}
    }

  if (surrogate[i].smile.substr(ipar,1)!="(")
    {
      if (arom==1 and carbon_arom(sum_len_group+ipar+1)==0)
	{
	  linked_to_arom(sum_len_group)=1;
	  linked_to_arom(sum_len_group+ipar+1)=1;
	}
      else if (arom==0 and carbon_arom(sum_len_group+ipar+1)==1)
	{
	  linked_to_arom(sum_len_group+ipar+1)=1;			
	  linked_to_arom(sum_len_group)=1;		    
	}

      if (surrogate[i].smile.substr(ipar+1,1)=="O" and ipar+1==int(surrogate[i].smile.length())-1)
	{
	  if (debug_mode==1) cout << "alcool1 " << sum_len_group+ipar+1 << " " << arom << endl;
	  if (arom==1)
	    {
	      carbon_arom(sum_len_group)=0;
	      surrogate[i].groups[28]+=1;
	    }
	  else
	    carbon_alcool(sum_len_group)++;
	}
      if (surrogate[i].smile.substr(ipar+1,2)=="O)")
	{
	  if (debug_mode==1) cout << "alcool2 " << sum_len_group+ipar+1 << " " << arom << endl;
	  if (arom==1)
	    {
	      carbon_arom(sum_len_group)=0;
	      surrogate[i].groups[28]+=1;
	    }
	  else
	    carbon_alcool(sum_len_group)++;
	}
    }

  if (surrogate[i].smile.substr(ipar+1,1)=="e" and carbon_anhydre(sum_len_group+ipar+1)<=0)
    {
      if (debug_mode==1) cout << "ester ipar+1" << endl;
      //carbon_ester(sum_len_group+ipar+1)++;
      //carbon_taken(sum_len_group)++;
      carbon_taken(sum_len_group+ipar+1)++;
      ester++;
    }
  
  if (surrogate[i].smile.substr(ipar+1,2)=="OC")
    {		      
      carbon_ether(sum_len_group,iether)=1;
      carbon_ether(sum_len_group+ipar+2,iether)=1;
      iether++;
    }
  if (surrogate[i].smile.substr(ipar+1,2)=="OA")
    {		      
      carbon_ether(sum_len_group,iether)=1;
      carbon_ether(sum_len_group+ipar+2,iether)=1;
      iether++;
    }
  if (surrogate[i].smile.substr(ipar+1,2)=="=O")
    {		      			
      carbon_ketone(sum_len_group,iket)=1;
      carbon_ketone(sum_len_group+2,iket)=1;
      carbon_taken(sum_len_group)=1;
      iket++;		      
    }
  
  if (surrogate[i].smile.substr(ipar+1,1)=="X" and carbon_taken(sum_len_group+ipar+1)==0)
    {		      
      carbon_ether(sum_len_group,iether)=1;
      for (icycle=0;icycle<5;icycle++)
	if (carbon_cycle(sum_len_group+ipar+1,icycle)>0)
	  for (igr=0;igr<ngr;igr++)
	    if (igr!=sum_len_group+ipar+1 and carbon_cycle(igr,icycle)>0)			 
	      carbon_ether(igr,iether)=1;
      iether++;
    }

  if ((surrogate[i].smile.substr(ipar+1,2)=="KF" and surrogate[i].smile.substr(ipar+1,3)!="KFK" and surrogate[i].smile.substr(ipar+1,3)!="KFL") or surrogate[i].smile.substr(ipar+1,1)=="H")
    {
      if (debug_mode==1) cout << "in: " << endl;
      carbon_ketone(sum_len_group+ipar,iket)=1;
      carbon_ketone(sum_len_group+2+ipar,iket)=1;
      carbon_taken(sum_len_group+1+ipar)=2;		    
      iket++;
      ester++;
    }
  else if (surrogate[i].smile.substr(ipar+1,1)=="K")
    {
      if (debug_mode==1) cout << "laK" << endl;
      carbon_ketone(sum_len_group,iket)=1;
      carbon_ketone(sum_len_group+ipar+2,iket)=1;
      carbon_taken(sum_len_group+ipar+1)++;
      iket++;
      /*
      carbon_ketone(sum_len_group+ipar+1,iket)=1;
      carbon_ketone(sum_len_group+ipar+3,iket)=1;
      carbon_taken(sum_len_group+ipar+2)=1;
      iket++;
      carbon_ketone(sum_len_group+ipar+2,iket)=1;
      carbon_ketone(sum_len_group+ipar+4,iket)=1;
      carbon_taken(sum_len_group+ipar+3)=1;
      iket++;*/
      int igrloc=ipar+2;
      int igrl=-1;
      int igrn=-1;
      if (surrogate[i].smile.substr(igrloc,1)=="L")
	igrl=igrloc;
      if (surrogate[i].smile.substr(igrloc,2)=="GK" or surrogate[i].smile.substr(igrloc,2)=="GL" or surrogate[i].smile.substr(igrloc,4)=="GC=O")
	{
	  //surrogate[i].groups[56]+=1;
	  //carbon_taken(igrloc)++;
	  //carbon_taken(igrloc+1)++;
	}        
      else if (surrogate[i].smile.substr(igrloc,1)=="G")
	{
	  if (debug_mode==1) cout << "Flaaa" << endl;
	  surrogate[i].groups[33]+=1;
	  carbon_taken(sum_len_group+igrloc)++;
	  carbon_taken(sum_len_group+igrloc-1)++;
	  if (debug_mode==1) cout << surrogate[i].smile.substr(igrloc,1) << " " <<  carbon_taken(sum_len_group+igrloc) << endl;
	}
      else if (surrogate[i].smile.substr(igrloc,1)=="e")
	{
	  if (carbon_taken(sum_len_group+igrloc)<=0)
	    {
	      if (debug_mode==1) cout << "KKKe" << endl;
	      surrogate[i].groups[33]+=1;
	      carbon_taken(sum_len_group+igrloc)++;
	      carbon_taken(sum_len_group+igrloc-1)++;
	    }
	}

      while (igrloc<int(surrogate[i].smile.length()))
	if (surrogate[i].smile.substr(igrloc,1)=="K")
	  {
	    carbon_ketone(sum_len_group+igrloc-1,iket)=1;
	    carbon_ketone(sum_len_group+igrloc+1,iket)=1;
	    carbon_taken(sum_len_group+igrloc)=1;
	    if (surrogate[i].smile.substr(igrloc+1,1)=="L")
	      igrl=igrloc+1;
	    if (surrogate[i].smile.substr(igrloc+1,2)=="GK" or surrogate[i].smile.substr(igrloc+1,2)=="GL" or surrogate[i].smile.substr(igrloc+1,4)=="GC=O")
	      {
		//surrogate[i].groups[56]+=1;
		//carbon_taken(igrloc+1)++;
		//carbon_taken(igrloc+2)++;
		//igrloc=igrloc+3;
	      }
	    else if (surrogate[i].smile.substr(igrloc+1,1)=="G")
	      {
		if (debug_mode==1) cout << "Flaaa2" << endl;
		surrogate[i].groups[33]+=1;
		carbon_taken(sum_len_group+igrloc+1)++;
		carbon_taken(sum_len_group+igrloc)++;
	      }
	    else if (surrogate[i].smile.substr(igrloc+1,1)=="e")
	      {
		if (carbon_taken(sum_len_group+igrloc+1)<=0)
		  {
		    if (debug_mode==1) cout << "KKKe2" << endl;
		    surrogate[i].groups[33]+=1;
		    carbon_taken(sum_len_group+igrloc+1)++;
		    carbon_taken(sum_len_group+igrloc)++;
		  }
	      }
				
	    iket++;
	    igrloc++;
	  }
	else
	  {
	    igrloc=surrogate[i].smile.length()+1;
	  }
		    
      if (igrl>0)
	{
	  if (debug_mode==1) cout << "find La " << surrogate[i].smile.substr(igrl,1) << " " << sum_len_group+igrl << " " << carbon_cycle(sum_len_group+igrl,0) << endl;
	  //cout << carbon_cycle << endl;
	  
	  int igr2=-1;
	  for (icycle=0;icycle<5;icycle++)
	    if (carbon_cycle(sum_len_group+igrl,icycle)>0)
	      {
		if (debug_mode==1) cout << "Lcycle " << carbon_cycle(sum_len_group+igrl,icycle) << " " << icycle << endl;
		for (igr=0;igr<ngr;igr++)
		  if (igr!=sum_len_group+igrl and carbon_cycle(igr,icycle)>0)
		    {
		      if (debug_mode==1) cout << "lallL " << igr << endl;
		      igr2=igr;			     
		    }
	      }

	  if (carbon_ketone(sum_len_group+igrl-1,iket)<0 and carbon_ketone(igr2,iket)<0)
	    {
	      if (debug_mode==1) cout << "add ketone " << endl;
	      carbon_ketone(sum_len_group+igrl-1,iket)=1;
	      carbon_ketone(igr2,iket)=1;
	      carbon_taken(sum_len_group+igrl)++;
	    }
		    
	  iket++;			

	}
			
    }
  /*
  else if (surrogate[i].smile.substr(ipar+1,3)=="KKK")
    {
      if (debug_mode==1) cout << "laKKK" << endl;
      carbon_ketone(sum_len_group,iket)=1;
      carbon_ketone(sum_len_group+ipar+2,iket)=1;
      carbon_taken(sum_len_group+ipar+1)=1;
      iket++;
      carbon_ketone(sum_len_group+ipar+1,iket)=1;
      carbon_ketone(sum_len_group+ipar+3,iket)=1;
      carbon_taken(sum_len_group+ipar+2)=1;
      iket++;
      carbon_ketone(sum_len_group+ipar+2,iket)=1;
      carbon_ketone(sum_len_group+ipar+4,iket)=1;
      carbon_taken(sum_len_group+ipar+3)=1;
      iket++;
      int igrloc=ipar+4;
      int igrl=-1;
      if (surrogate[i].smile.substr(igrloc,1)=="L")
	igrl=igrloc;
      if (surrogate[i].smile.substr(igrloc,2)=="GK")
	{
	  //surrogate[i].groups[56]+=1;
	  //carbon_taken(igrloc)++;
	  //carbon_taken(igrloc+1)++;
	}
      else if (surrogate[i].smile.substr(igrloc,1)=="F")
	{
	  if (debug_mode==1) cout << "Flaaa" << endl;
	  surrogate[i].groups[33]+=1;
	  carbon_taken(sum_len_group+igrloc)++;
	  carbon_taken(sum_len_group+igrloc-1)++;
	  if (debug_mode==1) cout << surrogate[i].smile.substr(igrloc,1) << " " <<  carbon_taken(sum_len_group+igrloc) << endl;
	}
      else if (surrogate[i].smile.substr(igrloc,1)=="e")
	{
	  if (carbon_taken(sum_len_group+igrloc)<=0)
	    {
	      if (debug_mode==1) cout << "KKKe" << endl;
	      surrogate[i].groups[33]+=1;
	      carbon_taken(sum_len_group+igrloc)++;
	      carbon_taken(sum_len_group+igrloc-1)++;
	    }
	}

      while (igrloc<int(surrogate[i].smile.length()))
	if (surrogate[i].smile.substr(igrloc,1)=="K")
	  {
	    carbon_ketone(sum_len_group+igrloc-1,iket)=1;
	    carbon_ketone(sum_len_group+igrloc+1,iket)=1;
	    carbon_taken(sum_len_group+igrloc)=1;
	    if (surrogate[i].smile.substr(igrloc+1,1)=="L")
	      igrl=igrloc+1;
	    if (surrogate[i].smile.substr(igrloc+1,2)=="FK")
	      {
		//surrogate[i].groups[56]+=1;
		//carbon_taken(igrloc+1)++;
		//carbon_taken(igrloc+2)++;
		//igrloc=igrloc+3;
	      }
	    else if (surrogate[i].smile.substr(igrloc+1,1)=="F")
	      {
		if (debug_mode==1) cout << "Flaaa2" << endl;
		surrogate[i].groups[33]+=1;
		carbon_taken(sum_len_group+igrloc+1)++;
		carbon_taken(sum_len_group+igrloc)++;
	      }
	    else if (surrogate[i].smile.substr(igrloc+1,1)=="e")
	      {
		if (carbon_taken(sum_len_group+igrloc+1)<=0)
		  {
		    if (debug_mode==1) cout << "KKKe2" << endl;
		    surrogate[i].groups[33]+=1;
		    carbon_taken(sum_len_group+igrloc+1)++;
		    carbon_taken(sum_len_group+igrloc)++;
		  }
	      }
				
	    iket++;
	    iket++;
	    igrloc++;
	  }
	else
	  {
	    igrloc=surrogate[i].smile.length()+1;
	  }
		    
      if (igrl>0)
	{
	  if (debug_mode==1) cout << "find L" << endl;
	  int igr2=-1;
	  for (icycle=0;icycle<5;icycle++)
	    if (carbon_cycle(sum_len_group+igrl,icycle)>0)
	      for (igr=0;igr<ngr;igr++)
		if (igr!=sum_len_group+igrl-1 and carbon_cycle(igr,icycle)>0)
		  {
		    //if (debug_mode==1) cout << "lall " << smile2 << " " << smile2.substr(igr,1) << endl;
		    igr2=igr;			     
		  }

	  if (carbon_ketone(sum_len_group+igrl-1,iket)<0 and carbon_ketone(igr2,iket)<0)
	    {
	      carbon_ketone(sum_len_group+igrl-1,iket)=1;
	      carbon_ketone(igr2,iket)=1;
	      carbon_taken(sum_len_group+igrl)++;
	    }
		    
	  iket++;			

	}
			
    }
  else if (surrogate[i].smile.substr(ipar+1,2)=="KL")
    {
      if (debug_mode==1) cout << "C2 " << endl;
      carbon_ketone(sum_len_group,iket)=1;
      carbon_ketone(sum_len_group+ipar+2,iket)=1;
      carbon_taken(sum_len_group+ipar+1)=1;
      iket++;
      //carbon_ketone(sum_len_group+1,iket)=1;
      //carbon_ketone(sum_len_group+3,iket)=1;
      //carbon_taken(sum_len_group+2)=1;

      int igr2=-1;
      for (icycle=0;icycle<5;icycle++)
	if (carbon_cycle(sum_len_group+ipar+2,icycle)>0)
	  for (igr=0;igr<ngr;igr++)
	    if (igr!=sum_len_group+ipar+1 and carbon_cycle(igr,icycle)>0)
	      {
		//if (debug_mode==1) cout << "lall " << smile2 << " " << smile2.substr(igr,1) << endl;
		igr2=igr;			     
	      }

      if (carbon_ketone(sum_len_group+ipar+1,iket)<0 and carbon_ketone(igr2,iket)<0)
	{
	  carbon_ketone(sum_len_group+ipar+1,iket)=1;
	  carbon_ketone(igr2,iket)=1;
	  carbon_taken(sum_len_group+ipar+2)=1;
	}
		    
      iket++;		    
    }
  else if (surrogate[i].smile.substr(ipar+1,2)=="KK")
    {
      if (debug_mode==1) cout << "la par+1" << endl;
      carbon_ketone(sum_len_group,iket)=1;
      carbon_ketone(sum_len_group+ipar+2,iket)=1;
      carbon_taken(sum_len_group+ipar+1)=1;
      iket++;
      carbon_ketone(sum_len_group+ipar+1,iket)=1;
      carbon_ketone(sum_len_group+ipar+3,iket)=1;
      carbon_taken(sum_len_group+ipar+2)=1;
      iket++;
      if (surrogate[i].smile.substr(ipar+3,1)=="F" and surrogate[i].smile.substr(ipar+3,2)!="FK" and surrogate[i].smile.substr(ipar+3,2)!="FL")
	{
	  //ester++;
	  //if (debug_mode==1) cout << "ester1 " << sum_len_group+igrloc << endl;
	  if (debug_mode==1) cout << "Flaes" << endl;
	  surrogate[i].groups[33]+=1;
	  carbon_taken(sum_len_group+ipar+1)++;
	  carbon_taken(sum_len_group+ipar+1)++;
	}
			
      int igrl=-1;
      if (surrogate[i].smile.substr(ipar+3,1)=="L")
	{
	  igrl=ipar+3;
	}
			
      if (igrl>0)
	{
	  if (debug_mode==1) cout << "find L" << endl;
	  int igr2=-1;
	  for (icycle=0;icycle<5;icycle++)
	    if (carbon_cycle(sum_len_group+igrl,icycle)>0)
	      for (igr=0;igr<ngr;igr++)
		if (igr!=sum_len_group+igrl-1 and carbon_cycle(igr,icycle)>0)
		  {
		    //if (debug_mode==1) cout << "lall " << smile2 << " " << smile2.substr(igr,1) << endl;
		    igr2=igr;			     
		  }
			    
	  if (carbon_ketone(sum_len_group+igrl-1,iket)<0 and carbon_ketone(igr2,iket)<0)
	    {
	      carbon_ketone(sum_len_group+igrl-1,iket)=1;
	      carbon_ketone(igr2,iket)=1;
	      carbon_taken(sum_len_group+igrl)++;
	    }
		    
	  iket++;			
			    
	}
    }
  else if (surrogate[i].smile.substr(ipar+1,1)=="K")
    {
      if (debug_mode==1) cout << "la2" << surrogate[i].smile.substr(ipar+1,2) << endl;
      carbon_ketone(sum_len_group,iket)=1;
      carbon_ketone(sum_len_group+ipar+2,iket)=1;
      carbon_taken(sum_len_group+ipar+1)=1;
      iket++;
    }*/
  if (surrogate[i].smile.substr(ipar+1,1)=="L")
    {
      if (debug_mode==1) cout << "ici4" << endl;
      //carbon_ketone(sum_len_group,iket)=1;
      int igr2=-1;
      for (icycle=0;icycle<5;icycle++)
	if (carbon_cycle(sum_len_group+ipar+1,icycle)>0)
	  for (igr=0;igr<ngr;igr++)
	    if (igr!=sum_len_group and carbon_cycle(igr,icycle)>0)
	      {
		//if (debug_mode==1) cout << "lall " << smile2 << " " << smile2.substr(igr,1) << endl;
		igr2=igr;			     
	      }
			
      if (carbon_ketone(sum_len_group,iket)<0 and carbon_ketone(igr2,iket)<0)
	{
	  if (debug_mode==1) cout << "vaaa" << endl;
	  carbon_ketone(sum_len_group,iket)=1;
	  carbon_ketone(igr2,iket)=1;
	  carbon_taken(sum_len_group+ipar+1)++;
	}
      iket++;
    }
  if (surrogate[i].smile.substr(ipar+1,3)=="OOC" or surrogate[i].smile.substr(ipar+1,3)=="OOK" or surrogate[i].smile.substr(ipar+1,3)=="OOA" or surrogate[i].smile.substr(ipar+1,3)=="OOL" or surrogate[i].smile.substr(ipar+1,3)=="OOG")
    {		      
      carbon_peroxide(sum_len_group,ipero)=1;
      carbon_peroxide(sum_len_group+ipar+3,ipero)=1;		    
      ipero++;
      carbon_taken(sum_len_group+ipar+1)=1;
      carbon_taken(sum_len_group+ipar+2)=1;
    }
  
  if (surrogate[i].smile.substr(ipar+1,3)=="OEC" or surrogate[i].smile.substr(ipar+1,3)=="OEK" or surrogate[i].smile.substr(ipar+1,3)=="OEL")
    {
      if (debug_mode==1) cout << "OOC" << endl;
      carbon_peroxide(sum_len_group,ipero)=1;
      carbon_peroxide(sum_len_group+ipar+2,ipero)=1;
      //carbon_taken(sum_len_group)++;
      //carbon_taken(sum_len_group+ipar+2)++;
      carbon_ketone(sum_len_group+ipar+1,iket)=1;
      carbon_ketone(sum_len_group+ipar+2,iket)=1;
      carbon_taken(sum_len_group+ipar+1)++;
      carbon_taken(sum_len_group+ipar+2)++;
      ipero++;
      iket++;
    }
  
  /*
  if (surrogate[i].smile.substr(ipar+1,1)=="Y")
    {		      
      carbon_peroxide(sum_len_group,ipero)=1;
      int igr3=-1;
      for (icycle=0;icycle<5;icycle++)
	if (carbon_cycle(sum_len_group+ipar+1,icycle)>0)
	  for (igr2=0;igr2<ngr;igr2++)
	    if (igr2!=igr and carbon_cycle(igr2,icycle)>0)
	      igr3=igr2;
      carbon_peroxide(sum_len_group+ipar+1,ipero)=1;		    
      ipero++;
    }*/
  if (surrogate[i].smile.substr(ipar+1,2)=="OX")
    {
      //cout << "OX" << endl;
      carbon_peroxide(sum_len_group,ipero)=1;
      for (icycle=0;icycle<5;icycle++)
	if (carbon_cycle(sum_len_group+ipar+2,icycle)>0)
	  for (igr=0;igr<ngr;igr++)
	    if (igr!=sum_len_group+ipar+2 and carbon_cycle(igr,icycle)>0)
	      {
		carbon_peroxide(igr,ipero)=1;
		//cout << "cycle " << icycle << " " << igr << " " << carbon_ketone(igr,0) << endl;
	      }
      carbon_taken(sum_len_group+ipar+1)++;
      carbon_taken(sum_len_group+ipar+2)++;
      ipero++;
    }

  

  /*
    if (surrogate[i].smile.substr(ipar+1,1)=="f")
    {
    if (debug_mode==1) cout << "esterlaf3" << endl;
    ester++;
    }*/

  if (surrogate[i].smile.substr(ipar,1)!="(")
    {
      		    
      if (surrogate[i].smile.substr(ipar+1,1)=="G" and surrogate[i].smile.substr(ipar+1,2)!="GK" and surrogate[i].smile.substr(ipar+1,2)!="GL" and surrogate[i].smile.substr(ipar+1,4)!="GC=O")
	{
	  if (debug_mode==1) cout << "esterA" << endl;
	  ester++;
	}
      if (surrogate[i].smile.substr(ipar+1,3)=="(O)")
	{
	  if (debug_mode==1) cout << "alcool3 " << sum_len_group+ipar+2 << endl;
	  carbon_alcool(sum_len_group)++;
	}
      if (surrogate[i].smile.substr(ipar+1,surrogate[i].smile.length()-1-ipar)=="OO")
	hydroxyperoxide++;
      if (surrogate[i].smile.substr(ipar+1,3)=="OO)")
	hydroxyperoxide++;		    
      //if (debug_mode==1) cout << "la " << surrogate[i].smile.substr(ipar+1,surrogate[i].smile.length()-1-ipar) << endl;
      if (surrogate[i].smile.substr(ipar+1,8)=="ON(=O)=O" or surrogate[i].smile.substr(ipar+1,13)=="O[N+](=O)[O-]")
	nitrate++;
      if (surrogate[i].smile.substr(ipar+1,7)=="N(=O)=O" or surrogate[i].smile.substr(ipar+1,12)=="[N+](=O)[O-]")
	nitrite++;
    }

  if (ngr>10)
    if (debug_mode==1) cout << carbon_taken(8) << endl;

}

