//!!-----------------------------------------------------------------------
//!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
//!!     SSH-aerosol is distributed under the GNU General Public License v3
//!!-----------------------------------------------------------------------

#include <fstream>
#include "smiles-functions.cxx"
#include <sstream>
using namespace ssh_soap;

void get_smiles(model_config &config, vector<species>& surrogate)
{
  int i,j;
  int n=surrogate.size();

  //Debug_mode activate prints to help find the issues in SMILES decomposition (for developers)
  int debug_mode=0;

  // output flux: 1=screen, 2=file
  int outFlux;
  if (config.SOAPlog==3)
    outFlux=1;
  else
    outFlux=config.SOAPlog;
  
  std::ofstream fileFlux;
  
  if (outFlux==2)
    fileFlux.open("smile2UNIFAC.decomp", ios::out); 

  Array<string,1> name_group;
  Array<int, 1> nc_group,no_group;
  name_group.resize(60);
  nc_group.resize(60);
  no_group.resize(60);
  nc_group=1;
  no_group=0;
  name_group(0)="CH3: ";
  name_group(1)="CH2: ";
  name_group(2)="CH: ";
  name_group(3)="C: ";
  name_group(4)="CH3 linked to alcohol: ";
  name_group(5)="CH2 linked to alcohol: ";
  name_group(6)="CH linked to alcohol: ";  
  name_group(7)="C linked to alcohol: ";
  name_group(8)="CH3 between several alcohol groups: ";
  name_group(9)="CH2 between several alcohol groups: ";
  name_group(10)="CH between several alcohol groups: ";
  name_group(11)="C between several alcohol groups: ";
  name_group(12)="CH3 in tails of alcohol molecule: ";
  name_group(13)="CH2 in tails of alcohol molecule: ";
  name_group(14)="CH in tails of alcohol molecule: ";
  name_group(15)="C in tails of alcohol molecule: ";
  name_group(16)="CH2=CH double bound: ";
  nc_group(16)=2;
  name_group(17)="CH=CH double bound: ";
  nc_group(17)=2;
  name_group(18)="CH2=C double bound: ";
  nc_group(18)=2;
  name_group(19)="CH=C double bound: ";
  nc_group(19)=2;
  name_group(20)="C=C double bound: ";
  nc_group(20)=2;
  name_group(21)="aromatic carbon with 1 hydrogen: ";
  name_group(22)="aromatic carbon with 0 hydrogen: ";
  name_group(23)="aromatic carbon linked to CH3: ";
  nc_group(23)=2;
  name_group(24)="aromatic carbon linked to CH2: ";
  nc_group(24)=2;
  name_group(25)="aromatic carbon linked to CH: ";
  nc_group(25)=2;
  name_group(26)="alcohol group: ";
  nc_group(26)=0;
  no_group(26)=1;
  name_group(27)="water: ";
  name_group(28)="phenol: ";
  no_group(28)=1;
  name_group(29)="CH3CO ketone: ";
  nc_group(29)=2;
  no_group(29)=1;
  name_group(30)="CH2CO ketone: ";
  nc_group(30)=2;
  no_group(30)=1;
  name_group(31)="aldehyde: ";
  no_group(31)=1;
  name_group(32)="CH3COO ester: "; 
  nc_group(32)=2;
  no_group(32)=2; 
  name_group(33)="CH2COO ester: ";
  nc_group(33)=2; 
  no_group(33)=2; 
  name_group(34)="CH3O ether: ";
  no_group(34)=1;
  name_group(35)="CH2O ether: ";
  no_group(35)=1;
  name_group(36)="CHO ether: ";
  no_group(36)=1;
  name_group(37)="acid: ";
  no_group(37)=2;
  name_group(38)="Nitro aromatic: ";
  no_group(38)=2;
  name_group(39)="CH2ONO2 nitrate: ";
  no_group(39)=3;
  name_group(40)="CHONO2 nitrate: ";
  no_group(40)=3;
  name_group(41)="CONO2 nitrate: ";
  no_group(41)=3;
  name_group(42)="CH2OOH hydroxyperoxide: ";
  no_group(42)=2;
  name_group(43)="CHOOH hydroxyperoxide: ";
  no_group(43)=2;
  name_group(44)="COOH hydroxyperoxide: ";
  no_group(44)=2;
  name_group(45)="CH3OOCH2 peroxide: ";
  no_group(45)=2;
  nc_group(45)=2;
  name_group(46)="CH3OOCH peroxide: ";
  no_group(46)=2;
  nc_group(46)=2;
  name_group(47)="CH3OOC peroxide: ";
  no_group(47)=2;
  nc_group(47)=2;
  name_group(48)="CH2OOCH2 peroxide: ";
  no_group(48)=2;
  nc_group(48)=2;
  name_group(49)="CH2OOCH peroxide: ";
  no_group(49)=2;
  nc_group(49)=2;
  name_group(50)="CH2OOC peroxide: ";
  no_group(50)=2;
  nc_group(50)=2;
  name_group(51)="CHOOCH peroxide: ";
  no_group(51)=2;
  nc_group(51)=2;
  name_group(52)="CHOOC peroxide: ";
  no_group(52)=2;
  nc_group(52)=2;
  name_group(53)="COOC peroxide: ";
  no_group(53)=2;
  nc_group(53)=2;
  name_group(54)="PAN: ";
  no_group(54)=5;
  name_group(55)="Peroxyacetyl acid: ";
  no_group(55)=3;
  name_group(56)="O=COC=O group: ";
  no_group(56)=3;
  nc_group(56)=2;
  name_group(57)="CH3NO2 group: ";
  no_group(57)=2;
  name_group(58)="CH2NO2 group: ";
  no_group(58)=2;
  name_group(59)="CHNO2 group: ";
  no_group(59)=2;
  
  for (i=0;i<n;i++)
    if (surrogate[i].smile!="" and surrogate[i].smile[0]!='&' and surrogate[i].smile[0]!='$' and surrogate[i].smile!="-")
      {
	if (outFlux==1)
	  {
	    cout <<"===="<< i <<"===="<<endl;
	    cout << surrogate[i].name << " is constructed from smiles: " << surrogate[i].smile << endl;
	  }
	else if (outFlux==2)
	  {
	    fileFlux <<"===="<< i <<"===="<<endl;
	    fileFlux << surrogate[i].name << " is constructed from smiles: " << surrogate[i].smile << endl;
	  }
		
	for (j=0; j<60; j++)
	  surrogate[i].groups[j]=0.;

	int ipos=0;
	string previous="";
	string smile_save=surrogate[i].smile;
	int nc_double_save=-1;
	int will_be_double=0;
	int will_be_ester=0;

	int ngr_max=10;

	Array<int, 1> carbon_alcool,carbon_near_alcool,carbon_tail,carbon_hydrogen,carbon_nh,carbon_taken,carbon_arom;
	Array<int, 1> linked_to_arom,carbon_near_alcool_tmp,carbon_near_alcool_save,carbon_near_alcool_bef,carbon_anhydre,carbon_ester;
	Array<int, 2> carbon_cycle,carbon_ether,carbon_ketone,carbon_peroxide;	
	int igr,igr2;
	int ngr=surrogate[i].smile.length();
        int total_c=0;
        int total_o=0;
	int icycle;
	int icycle2;
        for (igr=0;igr<ngr;igr++)
           {
	     if (surrogate[i].smile.substr(igr,1)=="O")
	       total_o++;
	     if (surrogate[i].smile.substr(igr,1)=="o")
	       total_o++;
	     if (surrogate[i].smile.substr(igr,1)=="c")
	       total_c++;
	     if (surrogate[i].smile.substr(igr,1)=="C")
	       total_c++;
	   }

	int found_arom=1;
	while (found_arom==1)
	  {
	    found_arom=0;
	    for (icycle=0;icycle<5;icycle++)
	      {
	    
		std::stringstream ss;
		ss << icycle+1;
		std::string s = ss.str();
		//cout << "C"+s << endl;
		int found=0;
		int igr1=-1;
		int igr2=-1;
		int check_consecutive_double=-1;
		//cout << "Search: C"+s << endl;
		for (igr=0;igr<ngr;igr++)
		  {
		    if ((surrogate[i].smile.substr(igr,2)=="C"+s or surrogate[i].smile.substr(igr,2)=="c"+s) and found==0)
		      {
			//cout << icycle << " " << igr << " " << "C"+s << endl;
			found++;
			igr1=igr;
		      }
		    else if ((surrogate[i].smile.substr(igr,2)=="C"+s or surrogate[i].smile.substr(igr,2)=="C"+s) and found==1)
		      {
			//cout << icycle << " " << igr << " " << "C"+s << endl;
			found++;
			igr2=igr;
		      }
		    else
		      {
			for (icycle2=0;icycle2<5;icycle2++)
			  if (icycle!=icycle2)
			    {
			      std::stringstream ss2;
			      ss2 << icycle2+1;
			      std::string s2 = ss2.str();
			      //cout << "C"+s2+s << endl;
			      if ((surrogate[i].smile.substr(igr,3)=="C"+s+s2 or surrogate[i].smile.substr(igr,3)=="C"+s2+s) and found==0)
				{
				  //cout << icycle << " " << igr << " " << "C"+s << endl;
				  found++;
				  igr1=igr;
				}
			      else if ((surrogate[i].smile.substr(igr,3)=="C"+s+s2 or surrogate[i].smile.substr(igr,3)=="C"+s2+s) and found==1)
				{
				  //cout << icycle << " " << igr << " " << "C"+s << endl;
				  found++;
				  igr2=igr;
				}
			      else if ((surrogate[i].smile.substr(igr,3)=="c"+s+s2 or surrogate[i].smile.substr(igr,3)=="c"+s2+s) and found==0)
				{
				  //cout << icycle << " " << igr << " " << "C"+s << endl;
				  found++;
				  igr1=igr;
				}
			      else if ((surrogate[i].smile.substr(igr,3)=="c"+s+s2 or surrogate[i].smile.substr(igr,3)=="c"+s2+s) and found==1)
				{
				  //cout << icycle << " " << igr << " " << "C"+s << endl;
				  found++;
				  igr2=igr;
				}
			    }
		      }
		  }


	    
		if (igr1>-1 and igr2>-1)
		  {
		    vector <int> Carom_position,Carom_dbound;
		    //Carom_position.push_back(igr2);
		    double ndouble=0;
		    int nC=1;
		    igr=igr2+1;
		    found=0;
		    for (icycle2=0;icycle2<5;icycle2++)
		      {
			std::stringstream ss2;
			ss2 << icycle2+1;
			std::string s2 = ss2.str();
			//cout << "C"+s2+s << endl;
			if (surrogate[i].smile.substr(igr+1,1)==s2)
			  {
			    found++;
			  }
		      }
		    if (found==1)
		      igr++;

		    
		    int igrn;
		    if (surrogate[i].smile.substr(igr1,1)=="c")
		      ndouble+=0.5;
		    while (igr>igr1)
		      {
			int igrf=igr;
			//cout << igr << " " << igr1 << " " << surrogate[i].smile.substr(igr,1) << endl;
			for (icycle2=0;icycle2<5;icycle2++)
			  {
			    std::stringstream ss;
			    ss << icycle2+1;
			    std::string s = ss.str();
			    if (surrogate[i].smile.substr(igr,1)==s)
			      {
				if (icycle2==icycle)
				  igrf--;
				else
				  {
				    int igrm=-1;
				    for (igrn=igr-2;igrn>igr1;igrn--)
				      {
					if (surrogate[i].smile.substr(igrn,2)=="C"+s or surrogate[i].smile.substr(igrn,2)=="c"+s)
					  {
					    igrm=igrn;
					    //if (igrm!=igr-1 and igrm!=
					    //cout << "Jump loop to " << igrm << " from " << igr << endl;
					    /*
					      if (surrogate[i].smile.substr(igr-2,2)=="=C")
					      {
					      //ndouble+=1;
					      nC+=1;
					      }*/
					    found=0;
					    int igrb=igr-1;
					    for (icycle2=0;icycle2<5;icycle2++)
					      {
						std::stringstream ss2;
						ss2 << icycle2+1;
						std::string s2 = ss2.str();
						//cout << "C"+s2+s << endl;
						if (surrogate[i].smile.substr(igrb,1)==s2)
						  {
						    found++;
						  }
					      }
					    if (found==1)
					      igrb--;
					    
					    if (surrogate[i].smile.substr(igrb,1)=="C")
					      {
						check_consecutive_double--;
						//cout << "C in " << igrb << endl;
						nC+=1;
					      }
					    if (surrogate[i].smile.substr(igrb,1)=="c")
					      {
						check_consecutive_double--;
						nC+=1;
						ndouble+=0.5;
					      }					
					  }
				      }

				    if (igrm<0 or igrm<=igr1)
				      igrf--;
				    else
				      igrf=igrm;
				  }
			      }
			  }
			if (surrogate[i].smile.substr(igr,1)=="O" or surrogate[i].smile.substr(igr,1)=="N")
			  {
			    //cout << "O or N found in " << igr << endl;
			    nC=0;
			    igrf=0;
			  }
			if (surrogate[i].smile.substr(igr,1)=="C")
			  {
			    //cout << "C in " << igr << endl;
			    nC++;
			    igrf=igr-1;
			    Carom_position.push_back(igr);
			    check_consecutive_double--;
			  }
			if (surrogate[i].smile.substr(igr,1)=="c")
			  {
			    //cout << "c in " << igr << endl;
			    nC++;
			    ndouble=ndouble+0.5;
			    igrf=igr-1;
			    check_consecutive_double--;
			  }		    
			if (surrogate[i].smile.substr(igr,1)=="=")
			  {
			    if (check_consecutive_double>0)
			      {
				//cout << "warning consecutive" << endl;
				igr=0;
				ndouble=0;
				nC=0;
			      }
			    //cout << "= in " << igr << endl;
			    ndouble++;
			    igrf=igr-1;
			    Carom_dbound.push_back(igr);
			    check_consecutive_double=2;
			  }
			if (surrogate[i].smile.substr(igr,1)=="(")
			  {
			    igrf=igr-1;
			  }
			if (surrogate[i].smile.substr(igr,1)==")")
			  {
			    //cout << "ramification" << endl;
			    int ipar=0;
			    for (igrn=igr-1;igrn>igr1;igrn--)
			      {			
				if (ipar==0 and surrogate[i].smile.substr(igrn,1)=="(")
				  {
				    ipar=igrn;
				    igrn=-1;
				  }
				else if (ipar<0 and surrogate[i].smile.substr(igrn,1)=="(")			
				  ipar=ipar+1;
				else if (ipar<=0 and surrogate[i].smile.substr(igrn,1)==")")
				  {
				    //cout << "une ramification s'est ouverte" << endl;
				    ipar=ipar-1;
				  }
			      }
			    igrf=ipar;
			  }
			igr=igrf;
		      }

		    //cout << "carbon and double bounds in loop: " << nC << " " << ndouble << endl;
		    if (nC==6 and ndouble>2.2 and int(Carom_dbound.size())>0)
		      {
			//cout << nC << " " << ndouble << endl;
			//cout << "aromatic ring detected" << endl;
			found_arom=1;
			//cout << "Carom: " << endl;
			Carom_position.push_back(igr1);
			/*for (igr=0;igr<int(Carom_position.size());igr++)
			  cout << "Arom " << igr << " " << Carom_position[igr] << endl;

			for (igr=0;igr<int(Carom_dbound.size());igr++)
			  cout << "Double to remove " << igr << " " << Carom_dbound[igr] << endl;

			  cout << int(Carom_dbound.size()) << endl;*/
			surrogate[i].smile=smile_save.substr(0,igr1);
			nC=0;
			ndouble=0;
			//cout << "next arom: " << Carom_position[int(Carom_position.size())-1-nC] << endl;
			//cout << "next double: " << Carom_dbound[int(Carom_dbound.size())-1-ndouble] << endl;
			for (igr=igr1;igr<=igr2;igr++)
			  {
			    if (Carom_position[int(Carom_position.size())-1-nC]==igr)
			      {

				surrogate[i].smile+="c";
				if (nC<int(Carom_position.size())-1)
				  {
				    nC++;
				    //cout << "next arom: " << Carom_position[int(Carom_position.size())-1-nC] << endl;
				  }
			      }
			    else if (Carom_dbound[int(Carom_dbound.size())-1-ndouble]==igr)
			      {
				if (ndouble<int(Carom_dbound.size())-1)
				  {
				    ndouble++;
				    //cout << "next double: " << Carom_dbound[int(Carom_dbound.size())-1-ndouble] << endl;
				  }
			      }
			    else
			      surrogate[i].smile+=smile_save.substr(igr,1);
			  }
			surrogate[i].smile+=smile_save.substr(igr2+1,ngr-igr2);
			//cout << "update: " << surrogate[i].smile << " " << smile_save << endl;
			//exit(0);
			smile_save=surrogate[i].smile;
			ngr=surrogate[i].smile.length();
		      }
		  }
	   
	      }
	  }
	//exit(0);

	//cout << surrogate[i].smile << endl;
	//Reposition loop number (for example C32 is replace by C23)
	for (icycle=0;icycle<5;icycle++)
	  {
	    
	    std::stringstream ss;
	    ss << icycle+1;
	    std::string s = ss.str();
	    for (icycle2=icycle+1;icycle2<5;icycle2++)
	      {
		std::stringstream ss2;
		ss2 << icycle2+1;
		std::string s2 = ss2.str();
		//cout << s2+s << endl;
		for (igr=0;igr<ngr;igr++)
		  {
		    
		    if (surrogate[i].smile.substr(igr,2)==s2+s)
		      {
			//cout << "found " << endl;
			surrogate[i].smile=smile_save.substr(0,igr)+s+s2+smile_save.substr(igr+2,ngr-igr-3);
		      }
		  }
	      }
	  }
	//cout << surrogate[i].smile << endl;
	smile_save=surrogate[i].smile;
	ngr=surrogate[i].smile.length();

	//exit(0);

	carbon_arom.resize(ngr+3);
	carbon_alcool.resize(ngr+3);
	carbon_tail.resize(ngr+3);
	carbon_near_alcool.resize(ngr+3);
	carbon_near_alcool_tmp.resize(ngr+3);
	carbon_near_alcool_save.resize(ngr+3);
	carbon_near_alcool_bef.resize(ngr+3);
	carbon_hydrogen.resize(ngr+3);
	carbon_cycle.resize(ngr+3,5);
	carbon_ether.resize(ngr+3,ngr_max);
	carbon_ester.resize(ngr+3);
	carbon_ketone.resize(ngr+3,ngr_max);
	carbon_peroxide.resize(ngr+3,ngr_max);
	carbon_anhydre.resize(ngr+3);
	carbon_nh.resize(ngr+3);
	carbon_taken.resize(ngr+3);
	linked_to_arom.resize(ngr+3);
	carbon_alcool=0;
	carbon_near_alcool=0;
	carbon_near_alcool_tmp=0;
	carbon_near_alcool_save=0;
	carbon_tail=0;
	carbon_hydrogen=-1;
	carbon_cycle=0;
	carbon_ether=-1;
	carbon_ketone=-1;
	carbon_taken=0;
	carbon_arom=0;
	carbon_nh=-1;
	carbon_ester=-1;
	linked_to_arom=0;
	carbon_near_alcool_bef=0;
	carbon_peroxide=-1;
	carbon_anhydre=-1;
	int iether=0;
	int iket=0;
	int ipero=0;
	
	string smile2="";   
	int last_pos=0;      
	igr=0;
	
	while (igr<ngr)
	  {
	    //cout << "ipos: " << igr << " " << surrogate[i].smile.substr(igr,1) << endl;
	    if (surrogate[i].smile.substr(igr,6)=="O=C1OO")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"LOO";
		carbon_cycle(ipos,0)=1;
		last_pos=igr+6;
		ipos=ipos+2;
	      }
	    else if (surrogate[i].smile.substr(igr,6)=="O=C2OO")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"LOO";
		carbon_cycle(ipos,1)=1;
		last_pos=igr+6;
		ipos=ipos+2;
	      }
	    else if (surrogate[i].smile.substr(igr,6)=="O=C3OO")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"LOO";
		carbon_cycle(ipos,2)=1;
		last_pos=igr+6;
		ipos=ipos+2;
	      }
	    else if (surrogate[i].smile.substr(igr,6)=="O=C4OO")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"LOO";
		carbon_cycle(ipos,3)=1;
		last_pos=igr+6;
		ipos=ipos+2;
	      }
	    else if (surrogate[i].smile.substr(igr,6)=="O=C5OO")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"LOO";
		carbon_cycle(ipos,4)=1;
		last_pos=igr+6;
		ipos=ipos+2;
	      }	
	    else if (surrogate[i].smile.substr(igr,5)=="O=C1O")
	      {
		if (debug_mode==1) cout << "f ester la beg" << endl;
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"f";
		carbon_cycle(ipos,0)=1;
		last_pos=igr+5;
	      }
	    else if (surrogate[i].smile.substr(igr,5)=="O=C2O")
	      {
		if (debug_mode==1) cout << "f ester la beg" << endl;
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"f";
		carbon_cycle(ipos,1)=1;
		last_pos=igr+5;
	      }
	    else if (surrogate[i].smile.substr(igr,5)=="O=C3O")
	      {
		if (debug_mode==1) cout << "f ester la beg" << endl;
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"f";
		carbon_cycle(ipos,2)=1;
		last_pos=igr+5;
	      }
	    else if (surrogate[i].smile.substr(igr,5)=="O=C4O")
	      {
		if (debug_mode==1) cout << "f ester la beg" << endl;
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"f";
		carbon_cycle(ipos,3)=1;
		last_pos=igr+5;
	      }
	    else if (surrogate[i].smile.substr(igr,5)=="O=C5O")
	      {
		if (debug_mode==1) cout << "f ester la beg" << endl;
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"f";
		carbon_cycle(ipos,4)=1;
		last_pos=igr+5;
	      }	
	    else if (surrogate[i].smile.substr(igr,4)=="O=C1")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"L";
		carbon_cycle(ipos,0)=1;	       
		last_pos=igr+4;
	      }	    
	    else if (surrogate[i].smile.substr(igr,4)=="O=C2")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"L";
		carbon_cycle(ipos,1)=1;	       
		last_pos=igr+4;
	      }
	    else if (surrogate[i].smile.substr(igr,4)=="O=C3")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"L";
		carbon_cycle(ipos,2)=1;	       
		last_pos=igr+4;
	      }
	    else if (surrogate[i].smile.substr(igr,4)=="O=C4")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"L";
		carbon_cycle(ipos,3)=1;	       
		last_pos=igr+4;
	      }
	    else if (surrogate[i].smile.substr(igr,4)=="O=C5")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"L";
		carbon_cycle(ipos,4)=1;	       
		last_pos=igr+4;
	      }
	    else if (surrogate[i].smile.substr(igr,4)=="C1=O")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"L";
		//cout << "Linit cycle " << ipos << endl;
		carbon_cycle(ipos,0)=1;	       
		last_pos=igr+4;
	      }	
	    else if (surrogate[i].smile.substr(igr,4)=="C2=O")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"L";
		carbon_cycle(ipos,1)=1;	       
		last_pos=igr+4;
	      }
	    else if (surrogate[i].smile.substr(igr,4)=="C3=O")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"L";
		carbon_cycle(ipos,2)=1;	       
		last_pos=igr+4;
	      }
	    else if (surrogate[i].smile.substr(igr,4)=="C4=O")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"L";
		carbon_cycle(ipos,3)=1;	       
		last_pos=igr+4;
	      }
	    else if (surrogate[i].smile.substr(igr,4)=="C5=O")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"L";
		carbon_cycle(ipos,4)=1;	       
		last_pos=igr+4;
	      }
	    else if (surrogate[i].smile.substr(igr,6)=="OOC1=O")
	      {
		if (debug_mode==1) cout << "OO ester la " << endl;
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"OO";		
		last_pos=igr+2;
		ipos=ipos+1;
	      }
	    else if (surrogate[i].smile.substr(igr,6)=="OOC2=O")
	      {
		if (debug_mode==1) cout << "OO ester la " << endl;
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"OO";		
		last_pos=igr+2;
		ipos=ipos+1;
	      }
	    else if (surrogate[i].smile.substr(igr,6)=="OOC3=O")
	      {
		if (debug_mode==1) cout << "OO ester la " << endl;
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"OO";		
		last_pos=igr+2;
		ipos=ipos+1;
	      }
	    else if (surrogate[i].smile.substr(igr,6)=="OOC4=O")
	      {
		if (debug_mode==1) cout << "OO ester la " << endl;
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"OO";		
		last_pos=igr+2;
		ipos=ipos+1;
	      }
	    else if (surrogate[i].smile.substr(igr,6)=="OOC5=O")
	      {
		if (debug_mode==1) cout << "OO ester la " << endl;
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"OO";		
		last_pos=igr+2;
		ipos=ipos+1;
	      }
	    else if (surrogate[i].smile.substr(igr,5)=="OC1=O")
	      {
		if (debug_mode==1) cout << "e ester la " << endl;
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"f";
		carbon_cycle(ipos,0)=1;
		/*
		int igr2;
		for (igr2=0;igr2<ipos;igr2++)
		  if (carbon_cycle(igr2,0)>0)
		    {
		      carbon_ester(igr2)=1;
		      if (debug_mode==1) cout << "found ester1 " << igr2 << endl;
		    }*/
		
		last_pos=igr+5;
	      }	
	    else if (surrogate[i].smile.substr(igr,5)=="OC2=O")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"f";
		carbon_cycle(ipos,1)=1;
		/*
		int igr2;
		for (igr2=0;igr2<ipos;igr2++)
		  if (carbon_cycle(igr2,1)>0)
		  carbon_ester(igr2)=1;*/
		last_pos=igr+5;
	      }
	    else if (surrogate[i].smile.substr(igr,5)=="OC3=O")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"f";
		carbon_cycle(ipos,2)=1;
		/*
		int igr2;
		for (igr2=0;igr2<ipos;igr2++)
		  if (carbon_cycle(igr2,2)>0)
		    carbon_ester(igr2)=1;*/
		last_pos=igr+5;
	      }
	    else if (surrogate[i].smile.substr(igr,5)=="OC4=O")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"f";
		carbon_cycle(ipos,3)=1;
		/*
		int igr2;
		for (igr2=0;igr2<ipos;igr2++)
		  if (carbon_cycle(igr2,3)>0)
		  carbon_ester(igr2)=1;*/
		last_pos=igr+5;
	      }
	    else if (surrogate[i].smile.substr(igr,5)=="OC5=O")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"f";
		carbon_cycle(ipos,4)=1;
		/*
		int igr2;
		for (igr2=0;igr2<ipos;igr2++)
		  if (carbon_cycle(igr2,4)>0)
		  carbon_ester(igr2)=1;*/
		last_pos=igr+5;
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="C12")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos+1);
		carbon_cycle(ipos,0)=1;
		carbon_cycle(ipos,1)=1;
		last_pos=igr+3;
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="C13")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos+1);
		carbon_cycle(ipos,0)=1;
		carbon_cycle(ipos,2)=1;
		last_pos=igr+3;
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="C14")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos+1);
		carbon_cycle(ipos,0)=1;
		carbon_cycle(ipos,3)=1;
		last_pos=igr+3;
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="C15")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos+1);
		carbon_cycle(ipos,0)=1;
		carbon_cycle(ipos,4)=1;
		last_pos=igr+3;
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="C23")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos+1);
		carbon_cycle(ipos,1)=1;
		carbon_cycle(ipos,2)=1;
		last_pos=igr+3;
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="C24")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos+1);
		carbon_cycle(ipos,1)=1;
		carbon_cycle(ipos,3)=1;
		last_pos=igr+3;
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="C25")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos+1);
		carbon_cycle(ipos,1)=1;
		carbon_cycle(ipos,4)=1;
		last_pos=igr+3;
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="C34")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos+1);
		carbon_cycle(ipos,2)=1;
		carbon_cycle(ipos,3)=1;
		last_pos=igr+3;
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="C35")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos+1);
		carbon_cycle(ipos,2)=1;
		carbon_cycle(ipos,4)=1;
		last_pos=igr+3;
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="C45")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos+1);
		carbon_cycle(ipos,3)=1;
		carbon_cycle(ipos,4)=1;
		last_pos=igr+3;
	      }
	    else if (surrogate[i].smile.substr(igr,2)=="C1")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos+1);
		carbon_cycle(ipos,0)=1;	       
		last_pos=igr+2;
	      }
	    else if (surrogate[i].smile.substr(igr,2)=="C2")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos+1);
		carbon_cycle(ipos,1)=1;	       
		last_pos=igr+2;
	      }
	    else if (surrogate[i].smile.substr(igr,2)=="C3")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos+1);
		carbon_cycle(ipos,2)=1;	       
		last_pos=igr+2;
	      }
	    else if (surrogate[i].smile.substr(igr,2)=="C4")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos+1);
		carbon_cycle(ipos,3)=1;	       
		last_pos=igr+2;
	      }
	    else if (surrogate[i].smile.substr(igr,2)=="C5")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos+1);
		carbon_cycle(ipos,4)=1;	       
		last_pos=igr+2;
	      }
	    /*
	    else if (surrogate[i].smile.substr(igr,3)=="OO1")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"Y";
		carbon_cycle(ipos,0)=1;	       
		last_pos=igr+3;		
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="OO2")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"Y";
		carbon_cycle(ipos,1)=1;	       
		last_pos=igr+3;
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="OO3")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"Y";
		carbon_cycle(ipos,2)=1;	       
		last_pos=igr+3;
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="OO4")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"Y";
		carbon_cycle(ipos,3)=1;	       
		last_pos=igr+3;
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="OO5")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"Y";
		carbon_cycle(ipos,4)=1;	       
		last_pos=igr+3;
		}*/
	    else if (surrogate[i].smile.substr(igr,2)=="O1")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"X";
		carbon_cycle(ipos,0)=1;	       
		last_pos=igr+2;		
	      }
	    else if (surrogate[i].smile.substr(igr,2)=="O2")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"X";
		carbon_cycle(ipos,1)=1;	       
		last_pos=igr+2;
	      }
	    else if (surrogate[i].smile.substr(igr,2)=="O3")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"X";
		carbon_cycle(ipos,2)=1;	       
		last_pos=igr+2;
	      }
	    else if (surrogate[i].smile.substr(igr,2)=="O4")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"X";
		carbon_cycle(ipos,3)=1;	       
		last_pos=igr+2;
	      }
	    else if (surrogate[i].smile.substr(igr,2)=="O5")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"X";
		carbon_cycle(ipos,4)=1;	       
		last_pos=igr+2;
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="Oc1" and igr==0)
	      {
		//if (debug_mode==1) cout << "found" << endl;
		//carbon_alcool(1)=1;
		smile2+="C(O)";
		last_pos=igr+3;
		carbon_arom(0)=1;
		ipos=ipos+3;
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="c12")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"C";
		carbon_cycle(ipos,0)=1;
		carbon_cycle(ipos,1)=1;		
		carbon_arom(ipos)=1;
		last_pos=igr+3;
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="c13")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"C";
		carbon_cycle(ipos,0)=1;
		carbon_cycle(ipos,2)=1;		
		carbon_arom(ipos)=1;
		last_pos=igr+3;
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="c14")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"C";
		carbon_cycle(ipos,0)=1;
		carbon_cycle(ipos,3)=1;		
		carbon_arom(ipos)=1;
		last_pos=igr+3;
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="c15")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"C";
		carbon_cycle(ipos,0)=1;
		carbon_cycle(ipos,4)=1;		
		carbon_arom(ipos)=1;
		last_pos=igr+3;
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="c23")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"C";
		carbon_cycle(ipos,1)=1;
		carbon_cycle(ipos,2)=1;		
		carbon_arom(ipos)=1;
		last_pos=igr+3;
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="c24")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"C";
		carbon_cycle(ipos,1)=1;
		carbon_cycle(ipos,3)=1;		
		carbon_arom(ipos)=1;
		last_pos=igr+3;
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="c25")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"C";
		carbon_cycle(ipos,1)=1;
		carbon_cycle(ipos,4)=1;		
		carbon_arom(ipos)=1;
		last_pos=igr+3;
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="c34")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"C";
		carbon_cycle(ipos,2)=1;
		carbon_cycle(ipos,3)=1;		
		carbon_arom(ipos)=1;
		last_pos=igr+3;
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="c35")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"C";
		carbon_cycle(ipos,2)=1;
		carbon_cycle(ipos,4)=1;		
		carbon_arom(ipos)=1;
		last_pos=igr+3;
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="c45")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"C";
		carbon_cycle(ipos,3)=1;
		carbon_cycle(ipos,4)=1;		
		carbon_arom(ipos)=1;
		last_pos=igr+3;
	      }
	    else if (surrogate[i].smile.substr(igr,2)=="c1")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"C";
		carbon_cycle(ipos,0)=1;
		carbon_arom(ipos)=1;
		last_pos=igr+2;
	      }
	    else if (surrogate[i].smile.substr(igr,2)=="c2")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"C";
		carbon_cycle(ipos,1)=1;
		carbon_arom(ipos)=1;
		last_pos=igr+2;
	      }
	    else if (surrogate[i].smile.substr(igr,2)=="c3")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"C";
		carbon_cycle(ipos,2)=1;
		carbon_arom(ipos)=1;
		last_pos=igr+2;
	      }
	    else if (surrogate[i].smile.substr(igr,2)=="c4")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"C";
		carbon_cycle(ipos,3)=1;
		carbon_arom(ipos)=1;
		last_pos=igr+2;
	      }
	    else if (surrogate[i].smile.substr(igr,2)=="c5")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"C";
		carbon_cycle(ipos,4)=1;
		carbon_arom(ipos)=1;
		last_pos=igr+2;
	      }
	    else if (surrogate[i].smile.substr(igr,1)=="c")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"C";
		carbon_arom(ipos)=1;
		last_pos=igr+1;
	      }
	    else if (surrogate[i].smile.substr(igr,2)=="o1")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"X";
		carbon_cycle(ipos,0)=1;	       
		last_pos=igr+2;		
	      }
	    else if (surrogate[i].smile.substr(igr,2)=="o2")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"X";
		carbon_cycle(ipos,1)=1;	       
		last_pos=igr+2;
	      }
	    else if (surrogate[i].smile.substr(igr,2)=="o3")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"X";
		carbon_cycle(ipos,2)=1;	       
		last_pos=igr+2;
	      }
	    else if (surrogate[i].smile.substr(igr,2)=="o4")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"X";
		carbon_cycle(ipos,3)=1;	       
		last_pos=igr+2;
	      }
	    else if (surrogate[i].smile.substr(igr,2)=="o5")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"X";
		carbon_cycle(ipos,4)=1;	       
		last_pos=igr+2;
	      }
	    else if (surrogate[i].smile.substr(igr,1)=="o")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"O";
		//carbon_arom(ipos)=1;
		last_pos=igr+1;
	      }
	 

	    if (surrogate[i].smile.substr(igr,14)=="O=C(OON(=O)=O)" and igr==0)
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"P";	      
		last_pos=igr+14;
	      }
	    else if (surrogate[i].smile.substr(igr,19)=="C(=O)OO[N+](=O)[O-]")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"P";	      
		last_pos=igr+19;		
	      }
	    else if (surrogate[i].smile.substr(igr,19)=="[O-][N+](=O)OOC(=O)")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"P";	
		last_pos=igr+19;		
	      }
	    //VIC PAN GECKO//
	    else if (surrogate[i].smile.substr(igr,16)=="C(=O)(OON(=O)=O)")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"P";	
		last_pos=igr+16;		
	      }
	    ///////////////
	    else if (surrogate[i].smile.substr(igr,14)=="O=N(=O)OOC(=O)")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"P";	
		last_pos=igr+14;		
	      }
	    else if (surrogate[i].smile.substr(igr,14)=="C(=O)OON(=O)=O")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"P";	      
		last_pos=igr+14;		
	      }
	    else if (surrogate[i].smile.substr(igr,10)=="OOC(=O)C=O" and igr!=0)
	      {	       
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"OOK";	      
		last_pos=igr+7;
		ipos=ipos+2;
	      }
	    else if (surrogate[i].smile.substr(igr,13)=="OOC(=O)C(=O)O" and igr!=0)
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"OOK";	      
		last_pos=igr+7;
		ipos=ipos+2;
	      }
	    else if (surrogate[i].smile.substr(igr,9)=="OC(=O)C=O" and igr!=0)
	      {	       
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"E";	      
		last_pos=igr+6;		
	      }
	    else if (surrogate[i].smile.substr(igr,12)=="OC(=O)C(=O)O" and igr!=0)
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"E";	      
		last_pos=igr+6;		
	      }
	    /*
	    else if (surrogate[i].smile.substr(igr,8)=="OOC(=O)C" and igr!=0)
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"E";	      
		last_pos=igr+7;
		if (outFlux == 1)
		  {
		    cout << "WARNING: group OOC(=O) not available" << endl;
		    cout << "Group OC(=O) used instead" << endl;
		  }
		else if (outFlux == 2)
		  {
		    fileFlux << "WARNING: group OOC(=O) not available" << endl;
		    fileFlux << "Group OC(=O) used instead" << endl;
		  }
		total_o--;
	      }*/
	    else if (surrogate[i].smile.substr(igr,7)=="OC(=O)C" and igr!=0)
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"E";	      
		last_pos=igr+6;		
	      }
	    else if (surrogate[i].smile.substr(igr,15)=="OO[N+](=O)[O-])")
	      {
		if (outFlux == 1)
		  {
		    cout << "WARNING: group OO[N+](=O)[O-]) not available" << endl;
		    cout << "Group O[N+](=O)[O-] used instead" << endl;
		  }
		else if (outFlux == 2)
		  {
		    fileFlux << "WARNING: group OO[N+](=O)[O-]) not available" << endl;
		    fileFlux << "Group O[N+](=O)[O-] used instead" << endl;
		  }
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"O[N+](=O)[O-])";	
		last_pos=igr+15;
		ipos=ipos+13;
		total_o--;
	      }
	    else if (surrogate[i].smile.substr(igr,9)=="OON(=O)=O")
	      {
		if (outFlux == 1)
		  {
		    cout << "WARNING: group OO[N+](=O)[O-]) not available" << endl;
		    cout << "Group O[N+](=O)[O-]) used instead" << endl;
		  }
		else if (outFlux == 2)
		  {
		    fileFlux << "WARNING: group OO[N+](=O)[O-]) not available" << endl;
		    fileFlux << "Group O[N+](=O)[O-]) used instead" << endl;
		  }
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"ON(=O)=O";	
		last_pos=igr+9;
		ipos=ipos+7;
		total_o--;
	      }
	    else if (surrogate[i].smile.substr(igr,7)=="C(=O)O1")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"e";
		if (debug_mode==1) cout << "here esterrr" << endl;
		//carbon_ester(last_pos-1)=1;
		last_pos=igr+7;	 		
		carbon_cycle(ipos,0)=1;						
	      }
	    else if (surrogate[i].smile.substr(igr,7)=="C(=O)O2")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"e";
		if (debug_mode==1) cout << "here esterrr" << endl;
		//carbon_ester(last_pos-1)=1;
		last_pos=igr+7;	 		
		carbon_cycle(ipos,1)=1;						
	      }
	    else if (surrogate[i].smile.substr(igr,7)=="C(=O)O3")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"e";
		if (debug_mode==1) cout << "here esterrr" << endl;
		//carbon_ester(last_pos-1)=1;
		last_pos=igr+7;	 		
		carbon_cycle(ipos,2)=1;						
	      }
	    else if (surrogate[i].smile.substr(igr,7)=="C(=O)O4")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"e";
		if (debug_mode==1) cout << "here esterrr" << endl;
		//carbon_ester(last_pos-1)=1;
		last_pos=igr+7;	 		
		carbon_cycle(ipos,3)=1;						
	      }
	    else if (surrogate[i].smile.substr(igr,7)=="C(=O)O5")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"e";
		if (debug_mode==1) cout << "here esterrr" << endl;
		//carbon_ester(last_pos-1)=1;
		last_pos=igr+7;	 		
		carbon_cycle(ipos,4)=1;						
	      }
	    else if (surrogate[i].smile.substr(igr,6)=="OC(=O)" and igr==0)
	      {
		smile2+="A";	      
		last_pos=igr+6;		
	      }
	     else if (surrogate[i].smile=="OC=O" and igr==0)
	      {
		smile2+="A";	      
		last_pos=igr+4;		
	      }
	    else if (surrogate[i].smile.substr(igr,7)=="OOC(=O)" and igr==0)
	      {
		smile2+="B";	      
		last_pos=igr+7;		
	      }
	    else if (surrogate[i].smile.substr(igr,4)=="OC12" and igr==0)
	      {
		//carbon_alcool(1)=1;
		carbon_cycle(ipos,0)=1;
		carbon_cycle(ipos,1)=1;
		smile2+="C(O)";
		last_pos=igr+4;
		ipos=ipos+3;
	      }
	    else if (surrogate[i].smile.substr(igr,4)=="OC13" and igr==0)
	      {
		//carbon_alcool(1)=1;
		carbon_cycle(ipos,0)=1;
		carbon_cycle(ipos,2)=1;
		smile2+="C(O)";
		last_pos=igr+4;
		ipos=ipos+3;
	      }
	    else if (surrogate[i].smile.substr(igr,4)=="OC14" and igr==0)
	      {
		//carbon_alcool(1)=1;
		carbon_cycle(ipos,0)=1;
		carbon_cycle(ipos,3)=1;
		smile2+="C(O)";
		last_pos=igr+4;
		ipos=ipos+3;
	      }
	    else if (surrogate[i].smile.substr(igr,4)=="OC15" and igr==0)
	      {
		//carbon_alcool(1)=1;
		carbon_cycle(ipos,0)=1;
		carbon_cycle(ipos,4)=1;
		smile2+="C(O)";
		last_pos=igr+4;
		ipos=ipos+3;
	      }
	    else if (surrogate[i].smile.substr(igr,4)=="OC23" and igr==0)
	      {
		//carbon_alcool(1)=1;
		carbon_cycle(ipos,1)=1;
		carbon_cycle(ipos,2)=1;
		smile2+="C(O)";
		last_pos=igr+4;
		ipos=ipos+3;
	      }
	    else if (surrogate[i].smile.substr(igr,4)=="OC24" and igr==0)
	      {
		//carbon_alcool(1)=1;
		carbon_cycle(ipos,1)=1;
		carbon_cycle(ipos,3)=1;
		smile2+="C(O)";
		last_pos=igr+4;
		ipos=ipos+3;
	      }
	    else if (surrogate[i].smile.substr(igr,4)=="OC25" and igr==0)
	      {
		//carbon_alcool(1)=1;
		carbon_cycle(ipos,1)=1;
		carbon_cycle(ipos,4)=1;
		smile2+="C(O)";
		last_pos=igr+4;
		ipos=ipos+3;
	      }
	    else if (surrogate[i].smile.substr(igr,4)=="OC34" and igr==0)
	      {
		//carbon_alcool(1)=1;
		carbon_cycle(ipos,2)=1;
		carbon_cycle(ipos,3)=1;
		smile2+="C(O)";
		last_pos=igr+4;
		ipos=ipos+3;
	      }
	    else if (surrogate[i].smile.substr(igr,4)=="OC35" and igr==0)
	      {
		//carbon_alcool(1)=1;
		carbon_cycle(ipos,2)=1;
		carbon_cycle(ipos,4)=1;
		smile2+="C(O)";
		last_pos=igr+4;
		ipos=ipos+3;
	      }
	    else if (surrogate[i].smile.substr(igr,4)=="OC45" and igr==0)
	      {
		//carbon_alcool(1)=1;
		carbon_cycle(ipos,3)=1;
		carbon_cycle(ipos,4)=1;
		smile2+="C(O)";
		last_pos=igr+4;
		ipos=ipos+3;
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="OC1" and igr==0)
	      {
		//carbon_alcool(1)=1;
		carbon_cycle(ipos,0)=1;
		smile2+="C(O)";
		last_pos=igr+3;
		ipos=ipos+3;
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="OC2" and igr==0)
	      {
		//carbon_alcool(1)=1;
		carbon_cycle(ipos,1)=1;
		smile2+="C(O)";
		last_pos=igr+3;
		ipos=ipos+3;
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="OC3" and igr==0)
	      {
		//carbon_alcool(1)=1;
		carbon_cycle(ipos,2)=1;
		smile2+="C(O)";
		last_pos=igr+3;
		ipos=ipos+3;
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="OC4" and igr==0)
	      {
		//carbon_alcool(1)=1;
		carbon_cycle(ipos,3)=1;
		smile2+="C(O)";
		last_pos=igr+3;
		ipos=ipos+3;
	      }
	    else if (surrogate[i].smile.substr(igr,3)=="OC5" and igr==0)
	      {
		//carbon_alcool(1)=1;
		carbon_cycle(ipos,4)=1;
		smile2+="C(O)";
		last_pos=igr+3;
		ipos=ipos+3;
	      }
	    
	    else if (surrogate[i].smile.substr(igr,2)=="OC" and igr==0)
	      {
		//carbon_alcool(1)=1;
		smile2+="C(O)";
		last_pos=igr+2;
		ipos=ipos+3;
	      }
	    /*
	    else if (surrogate[i].smile.substr(igr,12)=="C(=O)C(=O)OC")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"H";	      
		last_pos=igr+11;
	      }*/	    
	    else if (surrogate[i].smile.substr(igr,7)=="C(=O)OC")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"G";	      
		last_pos=igr+6;
	      }
	    else if (surrogate[i].smile.substr(igr,7)=="C(=O)OO" and surrogate[i].smile.substr(igr,8)!="C(=O)OOC" and surrogate[i].smile.substr(igr,8)!="C(=O)OON")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"B";	      
		last_pos=igr+7;	     
	      }
	    //VIC Acid GECKO//
	    else if (surrogate[i].smile.substr(igr,9)=="C(=O)(OO)")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"B";	
		last_pos=igr+9;		
	      }
	    ///////////////	    
	    else if (surrogate[i].smile.substr(igr,6)=="C(=O)O" and surrogate[i].smile.substr(igr,8)!="C(=O)OOC" and surrogate[i].smile.substr(igr,7)!="C(=O)OC" and surrogate[i].smile.substr(igr,7)!="C(=O)ON")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"A";	      
		last_pos=igr+6;	     
	      }
            //VIC Acid GECKO//
	    else if (surrogate[i].smile.substr(igr,8)=="C(=O)(O)")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"A";	
		last_pos=igr+8;		
	      }
	    ///////////////
	    else if (surrogate[i].smile.substr(igr,5)=="C(=O)" and igr==0)
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"O=C";	      
		last_pos=igr+5;	     
	      }
	    else if (surrogate[i].smile.substr(igr,5)=="C(=O)")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"K";	      
		last_pos=igr+5;	     
	      }

	    if (surrogate[i].smile.substr(igr,1)=="1" or surrogate[i].smile.substr(igr,1)=="2" or surrogate[i].smile.substr(igr,1)=="3" or surrogate[i].smile.substr(igr,1)=="4" or surrogate[i].smile.substr(igr,1)=="5")
	      ipos=ipos;
	    else
	      ipos++;
	    igr=max(last_pos,igr+1);
	  }

	if (last_pos<ngr)
	  smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos);	

	ngr=smile2.length();
	/*
	if (ngr>=2)
	  if (smile2.substr(ngr-2,2)=="CO" and carbon_arom(ngr-2)==0)
	    {
	      if (debug_mode==1) cout << "alcool0 " << endl;
	      carbon_alcool(ngr-2)=1;	      
	      }*/


	//Specific Case: aldehyde followed by an ether
	if (smile2.substr(0,5)=="O=COC")
	  {
	    carbon_taken(2)=1;
	    carbon_ether(2,iether)=1;
	    carbon_ether(4,iether)=1;
	    iether++;
	    previous="C";
	  }
	
	/*
	  if (sum(carbon_alcool)>0)	  
	  carbon_tail=1;			*/

	int will_be_hydroperoxide=0;
	int will_be_nitrate=0;
	int will_be_nitrite=0;
	int sum_len_group=0;
	surrogate[i].smile=smile2;
	smile_save=smile2;
	if (debug_mode==1) cout << smile2 << endl;
	if (debug_mode==1) cout << "cycle: " << carbon_cycle << endl;
	
	//Detect anhydre (checl if e is link to L and remove the carbon_ester that was created previously 
	for (igr=0;igr<ngr;igr++)
	  if (surrogate[i].smile.substr(igr,1)=="e")
	    {
	      if (debug_mode==1) cout << "la e " << endl;

	      for (icycle=0;icycle<5;icycle++)
		if (carbon_cycle(igr,icycle)>0)
		  for (igr2=0;igr2<ngr;igr2++)
		    if (igr2!=igr and carbon_cycle(igr2,icycle)>0)
		      {
			if (debug_mode==1) cout << igr2 << " " << icycle << " " << surrogate[i].smile.substr(igr2,1) << endl; 
			if (surrogate[i].smile.substr(igr2,1)=="L")
			  {
			    if (debug_mode==1) cout << "eL" << endl;
			    int igr3=igr-1;
			    while (igr3>=0)
			      {
				if (surrogate[i].smile.substr(igr3,1)==")")
				  {
				    int igr4=igr3-1;
				    int ramification=0;
				    while (igr4>=0)
				      if (surrogate[i].smile.substr(igr4,1)=="(" and ramification>0)
					{
					  ramification--;
					  igr4--;
					}
				      else if (surrogate[i].smile.substr(igr4,1)==")")
					{
					  ramification++;
					  igr4--;
					}
				      else if (surrogate[i].smile.substr(igr4,1)=="(" and ramification==0)
					{
					  igr3=igr4;
					  igr4=-1;
					}
				      else
					{
					  igr4--;
					}
				  }
				else
				  if (surrogate[i].smile.substr(igr3,1)=="C" and carbon_ester(igr3)>0)
				    {
				      carbon_ester(igr3)=-1;
				      igr3=-1;
				    }
				  else
				    igr3--;
			      }
			    carbon_ester(igr-1)=-1;
			    carbon_anhydre(igr)=1;
			    carbon_anhydre(igr2)=1;
			    carbon_taken(igr)=1;
			    carbon_taken(igr2)=1;
			  }
			/*
			else if (surrogate[i].smile.substr(igr2,1)=="C")
			  {
			    if (debug_mode==1) cout << "loop ester" << endl;
			    carbon_ester(igr2)=1;
			  }*/
		      }
	    }

	for (igr=0;igr<ngr;igr++)
	  if (surrogate[i].smile.substr(igr,1)=="f" and surrogate[i].smile.substr(igr,2)!="fK")
	    {
	      carbon_taken(igr)++;
	      if (debug_mode==1) cout << "la f " << endl;
	     
	      for (icycle=0;icycle<5;icycle++)
		if (carbon_cycle(igr,icycle)>0)
		  for (igr2=0;igr2<ngr;igr2++)
		    if (igr2!=igr and carbon_cycle(igr2,icycle)>0)
		      if (surrogate[i].smile.substr(igr2,1)=="L")
			{
			  if (debug_mode==1) cout << igr2 << " " << icycle << " " << surrogate[i].smile.substr(igr2,1) << endl; 
			  if (debug_mode==1) cout << "loop ester connecter to carbonyl" << endl;
			  carbon_taken(igr2)=1;
			  carbon_taken(igr)=1;
			  surrogate[i].groups[33]++;
			}
		      else
			{
			  if (debug_mode==1) cout << igr2 << " " << icycle << " " << surrogate[i].smile.substr(igr2,1) << endl; 
			  if (debug_mode==1) cout << "loop ester" << endl;
			  carbon_ester(igr2)=1;
			  
			}
	    }
	/*
	for (igr=0;igr<ngr;igr++)
	  if (surrogate[i].smile.substr(igr,1)=="L")
	    {
	      for (icycle=0;icycle<5;icycle++)
		if (carbon_cycle(igr,icycle)>0)
		  for (igr2=0;igr2<ngr;igr2++)
		    if (igr2!=igr and carbon_cycle(igr2,icycle)>0)
		      if (surrogate[i].smile.substr(igr2,1)=="O")
			{
			  if (debug_mode==1) cout << igr2 << " " << icycle << " " << surrogate[i].smile.substr(igr2,1) << endl; 
			  if (debug_mode==1) cout << "loop ester connecter to carbonyl" << endl;
			  carbon_taken(igr2)=1;
			  carbon_taken(igr)=1;
			  surrogate[i].groups[33]++;
			}
	    }*/
	  
	if (debug_mode==1) cout << carbon_ester << endl;	       
	while (surrogate[i].smile!="" and surrogate[i].smile[0]!='&')
	  {
	    //if (debug_mode==1) cout << surrogate[i].smile << endl;
	    int len_group=0; 
	    if (surrogate[i].smile.substr(0,2)=="AK" or surrogate[i].smile.substr(0,2)=="AG")
	      {	
		//if (debug_mode==1) cout << "is acid" << endl;
		surrogate[i].groups[37]+=1;
		//carbon_ketone(sum_len_group+2,iket)=1;
		//carbon_ketone(sum_len_group,iket)=1;
		carbon_taken(sum_len_group)=1;
		//iket++;
		len_group=-1;
		previous="C";
	      }
	    else if (surrogate[i].smile.substr(0,2)=="LK" and sum_len_group==0)
	      {
		if (debug_mode==1) cout << "C0 " << endl;
        
		/*
		carbon_ketone(sum_len_group,iket)=1;
		carbon_ketone(sum_len_group+2,iket)=1;
		carbon_taken(sum_len_group+1)=1;
		iket++;*/

		int igr2=-1;
		for (icycle=0;icycle<5;icycle++)
		  if (carbon_cycle(sum_len_group,icycle)>0)
		    for (igr=sum_len_group+1;igr<ngr;igr++)
		      if (carbon_cycle(igr,icycle)>0)
			{
			  if (debug_mode==1) cout << "lall " << smile2 << " " << smile2.substr(igr,1) << endl;
			  igr2=igr;			     
			}

		if (carbon_ketone(sum_len_group+1,iket)<0 and carbon_ketone(igr2,iket)<0)
		  {
		    carbon_ketone(sum_len_group+1,iket)=1;
		    carbon_ketone(igr2,iket)=1;
		    carbon_taken(sum_len_group)++;
		  }
		    
		iket++;
		len_group=-1;
		previous="O=";
	      }
	    else if (surrogate[i].smile.substr(0,3)=="AOC" or surrogate[i].smile.substr(0,4)=="AOOC" or surrogate[i].smile.substr(0,3)=="AOA")
	      {	
		//if (debug_mode==1) cout << "is acid" << endl;
		surrogate[i].groups[37]+=1;
		carbon_taken(sum_len_group)=1;
		//carbon_ether(sum_len_group,iether)=1;
		//carbon_ether(sum_len_group+2,iether)=1;
		//iether++;
		len_group=-1;
		previous="C";
	      }	    
	    else if (surrogate[i].smile.substr(0,1)=="A")
	      {	
		if (debug_mode==1) cout << "is acid" << endl;
		surrogate[i].groups[37]+=1;
		carbon_taken(sum_len_group)++;
		len_group=1;	      
	      }
	    /*
	    else if (surrogate[i].smile.substr(0,2)=="EK")
	      {	
		if (debug_mode==1) cout << "is EK" << endl;
		surrogate[i].groups[56]+=1;		
		len_group=2;
		carbon_taken(sum_len_group)=1;
		carbon_taken(sum_len_group+1)=1;
	      }*/
	    else if (surrogate[i].smile.substr(0,2)=="fK")
	      {	
		if (debug_mode==1) cout << "is fK" << endl;
		surrogate[i].groups[56]+=1;		
		len_group=-2;
		carbon_taken(sum_len_group)=1;
		carbon_taken(sum_len_group+1)=1;
		previous="C";
		//if (debug_mode==1) cout << "len: " << len_group << endl;
	      }
	    else if (surrogate[i].smile.substr(0,5)=="O=COK")
	      {	
		if (debug_mode==1) cout << "is O=COK" << endl;
		surrogate[i].groups[56]+=1;		
		len_group=-5;
		carbon_taken(sum_len_group)=1;
		carbon_taken(sum_len_group+4)=1;
		//if (debug_mode==1) cout << "len: " << len_group << endl;
	      }
	    else if (surrogate[i].smile.substr(0,4)=="GC=O")
	      {	
		if (debug_mode==1) cout << "is GC=O" << endl;
		surrogate[i].groups[56]+=1;		
		len_group=4;
		carbon_taken(sum_len_group)=1;
		carbon_taken(sum_len_group+1)=1;
		//if (debug_mode==1) cout << "len: " << len_group << endl;
	      }
	    else if (surrogate[i].smile.substr(0,2)=="GL")
	      {	
		if (debug_mode==1) cout << "is GL" << endl;
		surrogate[i].groups[56]+=1;		
		len_group=2;
		carbon_taken(sum_len_group)=1;
		carbon_taken(sum_len_group+1)=1;
		//if (debug_mode==1) cout << "len: " << len_group << endl;
	      }
	    else if (surrogate[i].smile.substr(0,4)=="GKOE")
	      {
		if (debug_mode==1) cout << "is GKOE" << endl;
		carbon_ester(sum_len_group)=-1;
		surrogate[i].groups[56]+=1;		
		len_group=-2;

		carbon_taken(sum_len_group)=1;
		carbon_taken(sum_len_group+1)=1;
		
	      }	    
	    else if (surrogate[i].smile.substr(0,4)=="GKOO")
	      {
		if (debug_mode==1) cout << "is GKOO" << endl;
		carbon_ester(sum_len_group)=-1;
		surrogate[i].groups[56]+=1;		
		len_group=-2;

		carbon_taken(sum_len_group)=1;
		carbon_taken(sum_len_group+1)=1;
		
	      }
	    else if (surrogate[i].smile.substr(0,2)=="GK")
	      {	
		if (debug_mode==1) cout << "is FK" << endl;
		carbon_ester(sum_len_group)=-1;
		surrogate[i].groups[56]+=1;		
		len_group=2;

		carbon_taken(sum_len_group)=1;
		carbon_taken(sum_len_group+1)=1;
		int igrloc=2;
		int igrl=-1;
		
		if (surrogate[i].smile.substr(igrloc,1)=="G" and surrogate[i].smile.substr(igrloc,4)!="GC=O")
		  {
		    if (debug_mode==1) cout << "is GKG" << endl;
		    carbon_taken(sum_len_group+1)++;
		    carbon_taken(sum_len_group+2)++;
		    surrogate[i].groups[33]++;
		  }
		
		if (surrogate[i].smile.substr(igrloc,1)=="L")
		  igrl=igrloc;
		
		while (igrloc<int(surrogate[i].smile.length()))
		  if (surrogate[i].smile.substr(igrloc,1)=="K")
		    {
		      if (debug_mode==1) cout << "found K" << endl;
		      carbon_ketone(sum_len_group+igrloc-1,iket)=1;
		      carbon_ketone(sum_len_group+igrloc+1,iket)=1;
		      carbon_taken(sum_len_group+igrloc)=1;
		      if (surrogate[i].smile.substr(igrloc+1,1)=="L")
			{
			  if (debug_mode==1) cout << "here L" << endl;
			  igrl=igrloc+1;
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
	    /*
	    else if (surrogate[i].smile.substr(0,2)=="EL")
	      {
		//if (debug_mode==1) cout << "is acid" << endl;
		surrogate[i].groups[56]+=1;		
		len_group=2;
		carbon_taken(sum_len_group)=1;
		carbon_taken(sum_len_group+1)=1;
		}*/
	    else if (surrogate[i].smile.substr(0,1)=="f")
	      {	
		if (debug_mode==1) cout << "final position loop f ester" << endl;
		//carbon_taken(sum_len_group)++;
		//if (debug_mode==1) cout << "taken: " << carbon_taken << endl;
		//surrogate[i].groups[33]+=1;		
		len_group=1;	      
	      }
	    else if (surrogate[i].smile.substr(0,1)=="e")
	      {	
		if (debug_mode==1) cout << "final position loop ester" << endl;
		//carbon_taken(sum_len_group)++;
		//if (debug_mode==1) cout << "taken: " << carbon_taken << endl;
		//surrogate[i].groups[33]+=1;		
		len_group=1;	      
	      }
	    /*
	    else if (surrogate[i].smile.substr(0,3)=="BKK")
	      {		
		//if (debug_mode==1) cout << "is acid" << endl;
		carbon_ketone(sum_len_group,iket)=1;
		carbon_ketone(sum_len_group+2,iket)=1;
		carbon_taken(sum_len_group)=1;
		carbon_taken(sum_len_group+1)=1;
		iket++;
		carbon_ketone(sum_len_group+1,iket)=1;
		carbon_ketone(sum_len_group+3,iket)=1;
		carbon_taken(sum_len_group+2)=1;
		iket++;
		surrogate[i].groups[55]+=1;		
		len_group=1;
		int igrloc=3;
		while (igrloc<int(surrogate[i].smile.length()))
		  if (surrogate[i].smile.substr(igrloc,1)=="K")
		    {
		      carbon_ketone(sum_len_group+igrloc-1,iket)=1;
		      carbon_ketone(sum_len_group+igrloc+1,iket)=1;
		      carbon_taken(sum_len_group+igrloc)=1;
		      iket++;
		      igrloc++;
		    }
		  else
		    {
		      igrloc=surrogate[i].smile.length()+1;
		    }		
	      }*/
	    else if (surrogate[i].smile.substr(0,2)=="BK" or surrogate[i].smile.substr(0,2)=="BG")
	      {		
		//if (debug_mode==1) cout << "is acid" << endl;
		//carbon_ketone(sum_len_group,iket)=1;
		//carbon_ketone(sum_len_group+2,iket)=1;
		//carbon_taken(sum_len_group+1)=1;
		carbon_taken(sum_len_group)=1;
		//iket++;
		surrogate[i].groups[55]+=1;		
		len_group=-1;
		previous="C";
	      }
	    else if (surrogate[i].smile.substr(0,1)=="B")
	      {		
		//if (debug_mode==1) cout << "is acid" << endl;
		surrogate[i].groups[55]+=1;
		carbon_taken(sum_len_group)=1;
		len_group=1;
	      }
	    /*
	    else if (surrogate[i].smile.substr(0,3)=="PKK")
	      {		
		//if (debug_mode==1) cout << "is acid" << endl;
		carbon_ketone(sum_len_group,iket)=1;
		carbon_ketone(sum_len_group+2,iket)=1;
		carbon_ketone(sum_len_group+1,iket+1)=1;
		carbon_ketone(sum_len_group+3,iket+1)=1;
		carbon_taken(sum_len_group)=1;
		carbon_taken(sum_len_group+1)=1;
		carbon_taken(sum_len_group+2)=1;
		iket=iket+2;
		surrogate[i].groups[54]+=1;
		len_group=2;
		int igrloc=3;
		while (igrloc<int(surrogate[i].smile.length()))
		  if (surrogate[i].smile.substr(igrloc,1)=="K")
		    {
		      carbon_ketone(sum_len_group+igrloc-1,iket)=1;
		      carbon_ketone(sum_len_group+igrloc+1,iket)=1;
		      carbon_taken(sum_len_group+igrloc)=1;
		      iket++;
		      igrloc++;
		    }
		  else
		    {
		      igrloc=surrogate[i].smile.length()+1;
		    }
	      }	   */ 
	    else if (surrogate[i].smile.substr(0,2)=="PK")
	      {		
		//if (debug_mode==1) cout << "is acid" << endl;
		//carbon_ketone(sum_len_group,iket)=1;
		//carbon_ketone(sum_len_group+2,iket)=1;
		carbon_taken(sum_len_group)=1;
		//carbon_taken(sum_len_group+1)=1;
		//iket++;
		surrogate[i].groups[54]+=1;
		len_group=-1;
		previous="C";
	      }
	    else if (surrogate[i].smile.substr(0,2)=="PG")
	      {		
		//if (debug_mode==1) cout << "is acid" << endl;
		//carbon_ketone(sum_len_group,iket)=1;
		//carbon_ketone(sum_len_group+2,iket)=1;
		carbon_taken(sum_len_group)=1;
		//carbon_taken(sum_len_group+1)=1;
		//iket++;
		surrogate[i].groups[54]+=1;
		len_group=-1;
		previous="C";
	      }	    
	    else if (surrogate[i].smile.substr(0,1)=="P")
	      {		
		//if (debug_mode==1) cout << "is acid" << endl;
		surrogate[i].groups[54]+=1;
		carbon_taken(sum_len_group)=1;
		len_group=1;
	      }
	    /*
	    else if (surrogate[i].smile.substr(0,2)=="LF")
	      {		
		if (debug_mode==1) cout << "LF" << endl;
		surrogate[i].groups[56]+=1;
		len_group=-2;
		carbon_taken(sum_len_group)=1;
		carbon_taken(sum_len_group+1)=1;
		previous="O=";
		}*/
	    else if (surrogate[i].smile.substr(0,3)=="LEK")
	      {		
		if (debug_mode==1) cout << "LEK" << endl;
		surrogate[i].groups[56]+=1;
		len_group=-2;
		carbon_taken(sum_len_group)=1;
		carbon_taken(sum_len_group+1)=1;
		previous="C";
	      }
	    else if (surrogate[i].smile.substr(0,2)=="LE")
	      {		
		//if (debug_mode==1) cout << "is acid" << endl;
		surrogate[i].groups[56]+=1;
		len_group=2;
		carbon_taken(sum_len_group)=1;
		carbon_taken(sum_len_group+1)=1;
	      }
	    
	    else if (surrogate[i].smile.substr(0,1)=="L" and (surrogate[i].smile.substr(0,2)!="LO" or surrogate[i].smile.substr(0,3)=="LOO") and surrogate[i].smile.substr(0,2)!="LF" and surrogate[i].smile.substr(0,2)!="LE")
	      {
		len_group=1;
		previous="C";
		if (debug_mode==1) cout << "Lbegin" << endl;	       
		
		if (carbon_anhydre(sum_len_group)>0)
		  {
		    surrogate[i].groups[56]+=1;
		    carbon_taken(sum_len_group)=1;
		  }
		else
		  {
		    if (sum_len_group==0)
		      {
			if (debug_mode==1) cout << "ici L beg" << endl;
			int iketloc=-1;
			int iket2;
			for (iket2=0;iket2<5;iket2++)
			  if (carbon_ketone(sum_len_group)==1)
			    {
			      iketloc=iket2;
			      carbon_ketone(sum_len_group)=-1;
			    }

			int igr2=-1;
			for (icycle=0;icycle<5;icycle++)
			  if (carbon_cycle(sum_len_group,icycle)>0)
			    for (igr=0;igr<ngr;igr++)
			      if (igr!=sum_len_group and carbon_cycle(igr,icycle)>0)
				{
				  if (debug_mode==1) cout << "lall " << smile2 << " " << smile2.substr(igr,1) << " " << igr2 << endl;
				  igr2=igr;
				}

			if (smile2.substr(igr2,1)=="X" and smile2.substr(max(igr2-1,0),2)!="OX")
			  {
			    if (debug_mode==1) cout << "X is found, will be ester" << " " << previous << endl;
			    will_be_ester=1;
			    carbon_taken(sum_len_group)=1;
			    carbon_taken(igr2)=1;
			  }
			else
			  if (iketloc==-1)
			    {
			      if (debug_mode==1) cout << "iketloc " << endl;
			      
			      iketloc=iket;
			      if (sum_len_group>0)
				carbon_ketone(sum_len_group-1,iketloc)=1;
			      else
				carbon_ketone(sum_len_group+1,iketloc)=1;
			      carbon_taken(sum_len_group)=1;
			      carbon_ketone(igr2,iketloc)=1;
			      iket2=-1;
			    }


		    
			if (iket2==-1)
			  iket++;
		      }
		  }

		if (surrogate[i].smile.substr(0,2)=="LG" and surrogate[i].smile.substr(0,3)!="LGK")
		  {
		    len_group=-1;
		    carbon_taken(sum_len_group+1)++;
		    if (debug_mode==1) cout << "LG is taken " << carbon_taken(sum_len_group) << endl;
		  }

		/*
		if (surrogate[i].smile.substr(0,3)=="LC=")
		  {
		    if (debug_mode==1) cout << "LC=" << endl;
		    nc_double_save=1;
		    carbon_taken(0)=2;
		    }*/
	      }    
    
	    else if (surrogate[i].smile.substr(0,15)=="C([N+](=O)[O-])" and carbon_arom(sum_len_group)==1)
	      {		
		//if (debug_mode==1) cout << "is nitro aromatique" << endl;
		surrogate[i].groups[38]+=1;
		len_group=15;
		carbon_arom(sum_len_group)=0;
	      }
	    else if (surrogate[i].smile.substr(0,4)=="C(O)" and carbon_arom(sum_len_group)==1)
	      {			       
		surrogate[i].groups[28]+=1;
		len_group=4;
		carbon_arom(sum_len_group)=0;
	      }
	    else if (surrogate[i].smile.substr(0,4)=="CO" and carbon_arom(sum_len_group)==1)
	      {			       
		surrogate[i].groups[28]+=1;
		len_group=2;
		carbon_arom(sum_len_group)=0;
	      }
	    else if (surrogate[i].smile.substr(0,3)=="C=O" or surrogate[i].smile.substr(0,3)=="O=C")
	      {
		int aldester=-1;
		int aldanhyd=-1;
		int aldOE=-1;
		int aldoo=-1;
		int aldnitrite=-1;
		if (debug_mode==1) cout << "begin with carbonyl " << endl;
		if (will_be_double>0 and nc_double_save>-1)
		  {
		    if (debug_mode==1) cout << "ald: " << will_be_double << " " << previous << endl;
		
		    //if (debug_mode==1) cout << "double bounds between " << nc_double_save <<  "hydrogen carbon and " << nC << " hydrogen carbon" << endl;

		    int nC1=nc_double_save;
		    int nC2=0;		    
		
		    if (nC1==2 and nC2==1)
		      surrogate[i].groups[16]+=1;
		    else if (nC1==1 and nC2==1)
		      surrogate[i].groups[17]+=1;
		    else if (nC1==2 and nC2==0)
		      surrogate[i].groups[18]+=1;
		    else if (nC1==1 and nC2==0)
		      surrogate[i].groups[19]+=1;
		    else if (nC1==0 and nC2==0)
		      surrogate[i].groups[20]+=1;
		    // VIC C+C GECKO //
		    else if  (nC1==2 and nC2==2)
		      {
			surrogate[i].groups[16]+=1;
			if (outFlux==1)
			  {
			    cout << "WARNING: CH2=CH2 not available" << endl;
			    cout << "Group CH2=CH used instead" << endl;
			  }
			else if (outFlux==2)
			  {
			    fileFlux <<"WARNING: CH2=CH2 not available" << endl;
			    fileFlux << "Group CH2=CH used instead" << endl;
			  }			   
		      }
		    ///////////////////  
		    else
		      {			
			cout << "error double bounds type " << nC1 << " " << nC2 << endl;
			cout << "problem in smiles.cxx. Try manual group decomposition" << endl;	
			if (outFlux==2)
			  {
			    fileFlux << "error double bounds type " << nC1 << " " << nC2 << endl;
			    fileFlux << "problem in smiles.cxx. Try manual group decomposition" << endl;
			  }
			exit(0);
		      }
		    will_be_double=0;
		    nc_double_save=-1;
		    previous="C=C";
		    carbon_taken(sum_len_group)++;
		  }
		int is_pan=0;
		//int double_ketone=0;
		int ketone=0;


		
		if (surrogate[i].smile.substr(0,4)=="O=C(")
		  {
		    previous="O=C";
		    int ipar=0;
		    for (igr=4;igr<int(surrogate[i].smile.length());igr++)
		      {			
			if (ipar==0 and surrogate[i].smile.substr(igr,1)==")")			
			  ipar=igr;
			else if (ipar<0 and surrogate[i].smile.substr(igr,1)==")")			
			  ipar=ipar+1;
			else if (ipar<=0 and surrogate[i].smile.substr(igr,1)=="(" and igr>2)
			  {
			    //cout << "une ramification s'est ouverte" << endl;
			    ipar=ipar-1;
			  }
		      }
		    if (surrogate[i].smile.length()-ipar-8>0)
		      if (surrogate[i].smile.substr(ipar+1,8)=="ON(=O)=O")
			is_pan=1;

		    if (surrogate[i].smile.length()-ipar-1>0)
		      if (surrogate[i].smile.substr(ipar+1,1)=="K")
			{
			  if (debug_mode==1) cout << "here " << ipar+1+sum_len_group << endl;
			  if (debug_mode==1) cout << iket << " ketones before" << endl;
			  ketone=1;
			  carbon_ketone(sum_len_group+ipar+1,iket)=1;
			  carbon_ketone(sum_len_group+4,iket)=1;
			  carbon_taken(sum_len_group+2)=1;
			  iket++;
			  //carbon_ketone(sum_len_group+4,iket)=1;
			  carbon_ketone(sum_len_group+2,iket)=1;
			  carbon_ketone(sum_len_group+ipar+2,iket)=1;
			  carbon_taken(sum_len_group+ipar+1)=1;
			  
			  iket++;
			  int igrloc=ipar+2;
			  while (igrloc<int(surrogate[i].smile.length()))
			    if (surrogate[i].smile.substr(igrloc,1)=="K")
			      {
				if (debug_mode==1) cout << "find K" << endl;
				carbon_ketone(sum_len_group+igrloc-1,iket)=1;
				carbon_ketone(sum_len_group+igrloc+1,iket)=1;
				carbon_taken(sum_len_group+igrloc)=1;
				iket++;
				igrloc++;
			      }
			    else
			      {
				igrloc=surrogate[i].smile.length()+1;
			      }
			  if (debug_mode==1) cout << iket << " ketones after" << endl;
			}
		    

		    if (surrogate[i].smile.length()-ipar-1>0)
		      if (surrogate[i].smile.substr(ipar+1,2)=="OC" or surrogate[i].smile.substr(ipar+1,2)=="OA")
			{
			  if (debug_mode==1) cout << "aldester" << endl;
			  aldester=1;
			  carbon_taken(2)=1;
			  //carbon_ketone(ipar+1,iket)=1;
			  //carbon_ketone(2,iket)=1;
			  //iket++;
			  //carbon_ester(4)=1;			  
			  surrogate[i].groups[33]++;
			  carbon_taken(ipar+1)=1;
			  carbon_taken(4)++;
			}		    
		    
		    if (surrogate[i].smile.length()-ipar-1>0)
		      if (surrogate[i].smile.substr(ipar+1,1)=="E")
			{
			  if (debug_mode==1) cout << "aldanhyd" << endl;
			  aldanhyd=1;
			  carbon_taken(ipar+1)++;
			  carbon_taken(2)=1;
			  carbon_anhydre(ipar+1)=1;
			}

		    if (surrogate[i].smile.length()-ipar-1>0)
		      if (surrogate[i].smile.substr(ipar+1,2)=="OE")
			{
			  if (debug_mode==1) cout << "aldOE" << endl;
			  aldOE=1;
			  //carbon_taken(ipar+1)++;
			  carbon_taken(2)=1;
			  carbon_ketone(ipar+1,iket)=1;
			  carbon_ketone(4,iket)=1;
			  iket++;
			}
		    if (surrogate[i].smile.length()-ipar-1>0)
		      if (surrogate[i].smile.substr(ipar+1,2)=="OO")
			{
			  if (debug_mode==1) cout << "aldOO" << endl;
			  aldoo=1;
			  //carbon_taken(ipar+1)++;
			  carbon_taken(2)=1;
			  carbon_ketone(ipar+1,iket)=1;
			  carbon_ketone(4,iket)=1;
			  iket++;
			}

		    if (surrogate[i].smile.length()-ipar-1>0)
		      if (surrogate[i].smile.substr(ipar+1,7)=="N(=O)=O")
			{
			  if (debug_mode==1) cout << "aldnitrite" << endl;
			  aldnitrite=1;
			  carbon_taken(ipar+1)++;
			  carbon_taken(2)=1;
			  carbon_ketone(ipar+1,iket)=1;
			  carbon_ketone(4,iket)=1;
			  iket++;
			}

		    if (surrogate[i].smile.substr(0,6)=="O=C(OO")
		      {
			if (debug_mode==1) cout << "aldanhyd O=C(OO" << endl;
			aldoo=1;
			carbon_taken(2)++;
			carbon_ketone(ipar+1,iket)=1;
			carbon_ketone(4,iket)=1;
			iket++;
		      }

		    	
		    if (surrogate[i].smile.substr(0,11)=="O=C(N(=O)=O")
		      {
			if (debug_mode==1) cout << "aldnitrite" << endl;
			aldnitrite=1;
			carbon_taken(4)++;
			carbon_taken(2)=1;
			carbon_ketone(ipar+1,iket)=1;
			carbon_ketone(4,iket)=1;
			iket++;
		      }

		    if (surrogate[i].smile.substr(0,6)=="O=C(OE")
		      {
			if (debug_mode==1) cout << "aldOE" << endl;
			aldOE=1;
			//carbon_taken(ipar+1)++;
			carbon_taken(2)=1;
			carbon_ketone(ipar+1,iket)=1;
			carbon_ketone(4,iket)=1;
			iket++;
		      }

		    if (surrogate[i].smile.substr(0,5)=="O=C(E")
		      {
			if (debug_mode==1) cout << "aldanhyd O=C(E" << endl;
			aldanhyd=1;
			carbon_taken(4)++;
			carbon_taken(2)=1;
			carbon_anhydre(4)=1;
		      }
		    if (surrogate[i].smile.substr(0,5)=="O=C(F")
		      {
			if (debug_mode==1) cout << "aldanhyd O=C(F" << endl;
			aldanhyd=1;
			carbon_taken(4)++;
			carbon_taken(2)=1;
			carbon_anhydre(4)=1;
		      }
		    if (surrogate[i].smile.substr(0,6)=="O=C(OC" or surrogate[i].smile.substr(0,6)=="O=C(OA")
		      {
			if (debug_mode==1) cout << "aldester" << endl;
			aldester=1;		     
			//carbon_ketone(ipar+1,iket)=1;
			//carbon_ketone(2,iket)=1;
			//iket++;
			//carbon_ester(ipar+1)=1;			  
			surrogate[i].groups[33]++;
			carbon_taken(2)=1;
			carbon_taken(4)++;
			carbon_taken(ipar+1)++;
		      }


		    if (aldanhyd<=0 and aldester<=0 and aldoo<=0 and aldOE<=0 and aldnitrite<=0)
		      if (surrogate[i].smile.length()-ipar-1>0)
			if (surrogate[i].smile.substr(ipar+1,1)=="C")
			  {
			    if (debug_mode==1) cout << "kethere" << endl;
			    ketone=1;
			    carbon_ketone(sum_len_group+4,iket)=1;
			    carbon_ketone(ipar+1,iket)=1;
			    carbon_taken(sum_len_group+2)=1;
			    iket++;
			  }

		    if (surrogate[i].smile.substr(0,5)=="O=C(K" and aldoo==-1 and aldOE==-1 and aldnitrite==-1)
		      {
			if (debug_mode==1) cout << "found O=C(K" << endl;
			carbon_ketone(sum_len_group+2,iket)=1;
			carbon_ketone(sum_len_group+5,iket)=1;
			carbon_taken(sum_len_group+4)++;
			iket++;
			/*
			carbon_ketone(sum_len_group+4,iket)=1;
			carbon_ketone(sum_len_group+6,iket)=1;
			carbon_taken(sum_len_group+5)=1;
			iket++;*/
			int igrloc=5;
			while (igrloc<int(surrogate[i].smile.length()))
			  if (surrogate[i].smile.substr(igrloc,1)=="K")
			    {
			      carbon_ketone(sum_len_group+igrloc-1,iket)=1;
			      carbon_ketone(sum_len_group+igrloc+1,iket)=1;
			      carbon_taken(sum_len_group+igrloc)++;
			      iket++;
			      igrloc++;
			    }
			  else
			    {
			      igrloc=surrogate[i].smile.length()+1;
			    }
		      }

		    //if (surrogate[i].smile.substr(0,4)=="O=C=")
		    //  nc_double_save=0;
		    /*
		    else if (surrogate[i].smile.substr(0,5)=="O=C(F")
		      {
			aldanhyd=1;
			carbon_taken(sum_len_group+4)=1;
			carbon_taken(sum_len_group+2)=1;
			}*/

		    		    
		    //if (debug_mode==1) cout << surrogate[i].smile.substr(ipar+1,6) << endl;
		    //ON(=O)=O
		  }

		if (debug_mode==1) cout << "aldanhyd: " << aldanhyd << " " << ketone << endl;
		if (is_pan==1)
		  {
		    total_o++;
		    surrogate[i].groups[54]+=1;
		    len_group=3;
		    carbon_taken(sum_len_group+2)++;
		  }
		else if (aldanhyd>0)
		  {
		    if (debug_mode==1) cout << "aldanhyd: " << endl;
		    surrogate[i].groups[56]+=1;
		    len_group=3;
		  }
		else if (aldester>0)
		  {
		    if (debug_mode==1) cout << "aldester: " << endl;
		    len_group=3;
		    previous="reset";
		  }
		else if (aldOE>0)
		  {
		    if (debug_mode==1) cout << "aldOE: " << endl;
		    len_group=-3;
		    previous="reset";
		  }
		else if (aldoo>0)
		  {
		    if (debug_mode==1) cout << "aldoo: " << endl;
		    len_group=-3;
		    previous="reset";
		  }	
		else if (aldnitrite>0)
		  {
		    if (debug_mode==1) cout << "aldnitrite: " << endl;
		    len_group=-3;
		    previous="reset";
		  }	
		else if (ketone==1)
		  {
		    if (debug_mode==1) cout << "aldketone " << endl;
		    len_group=3;
		  }
		else if (surrogate[i].smile.substr(0,4)=="O=CE")
		  {
		    if (debug_mode==1) cout << "O=CE" << endl;
		    carbon_taken(sum_len_group+3)=1;
		    surrogate[i].groups[56]+=1;
		    len_group=-4;
		    previous="O";
		  }		
		else 
		  {
		    previous="O=C";
		    if (debug_mode==1) cout << "is aldehyde" << " " << previous << endl;
		    surrogate[i].groups[31]+=1;
		    len_group=3;
		    /*
		    if (surrogate[i].smile.substr(0,6)=="O=CKKK")
		      {
			//if (debug_mode==1) cout << "ka" << endl;
			carbon_taken(sum_len_group+2)=1;
			carbon_ketone(sum_len_group+2,iket)=1;
			carbon_ketone(sum_len_group+4,iket)=1;			
			carbon_taken(sum_len_group+3)=1;
			iket++;
			carbon_ketone(sum_len_group+3,iket)=1;
			carbon_ketone(sum_len_group+5,iket)=1;
			carbon_taken(sum_len_group+4)=1;
			iket++;
			carbon_ketone(sum_len_group+4,iket)=1;
			carbon_ketone(sum_len_group+6,iket)=1;
			carbon_taken(sum_len_group+5)=1;
			iket++;
			int igrloc=6;
			while (igrloc<int(surrogate[i].smile.length()))
			  if (surrogate[i].smile.substr(igrloc,1)=="K")
			    {
			      carbon_ketone(sum_len_group+igrloc-1,iket)=1;
			      carbon_ketone(sum_len_group+igrloc+1,iket)=1;
			      carbon_taken(sum_len_group+igrloc)=1;
			      iket++;
			      igrloc++;
			    }
			  else
			    {
			      igrloc=surrogate[i].smile.length()+1;
			    }			
		      }
		    else if (surrogate[i].smile.substr(0,5)=="O=CKK")
		      {
			//if (debug_mode==1) cout << "ka2" << endl;
			carbon_ketone(sum_len_group+2,iket)=1;
			carbon_ketone(sum_len_group+4,iket)=1;
			carbon_taken(sum_len_group+2)=1;
			carbon_taken(sum_len_group+3)=1;
			iket++;
			carbon_ketone(sum_len_group+3,iket)=1;
			carbon_ketone(sum_len_group+5,iket)=1;
			carbon_taken(sum_len_group+4)=1;
			iket++;
		      }
		    else if (surrogate[i].smile.substr(0,5)=="O=CKF")
		      {
			//if (debug_mode==1) cout << "ka3" << endl;
			carbon_ketone(sum_len_group+2,iket)=1;
			carbon_ketone(sum_len_group+4,iket)=1;
			surrogate[i].groups[33]+=1;
			carbon_taken(sum_len_group+2)=1;
			carbon_taken(sum_len_group+3)=2;
			carbon_taken(sum_len_group+4)=1;
			len_group=5;
			iket++;
		      }		    
		    else if (surrogate[i].smile.substr(0,4)=="O=CK")
		      {
			//if (debug_mode==1) cout << "ka3" << endl;
			carbon_ketone(sum_len_group+2,iket)=1;
			carbon_ketone(sum_len_group+4,iket)=1;
			carbon_taken(sum_len_group+2)=1;
			carbon_taken(sum_len_group+3)=1;
			iket++;
			}*/
		    if (surrogate[i].smile.substr(0,4)=="O=CK")
		      {
			if (debug_mode==1) cout << "O=CK" << endl;
			carbon_taken(sum_len_group+2)=1;
			len_group=-3;
			previous="O=";
		      }

		    if (surrogate[i].smile.substr(0,5)=="O=COO")
		      {
			if (debug_mode==1) cout << "O=COO" << endl;
			carbon_taken(sum_len_group+2)=1;
			len_group=-3;
			previous="reset";
		      }		    		      
		    if (surrogate[i].smile.substr(0,5)=="O=COE")
		      {
			if (debug_mode==1) cout << "O=COE" << endl;
			carbon_taken(sum_len_group+2)=1;
			len_group=-3;
			previous="reset";
		      }
		    else if (surrogate[i].smile.substr(0,4)=="O=CG")
		      {
			//carbon_ester(sum_len_group+2)=1;
			//surrogate[i].groups[33]+=1;
			//carbon_ketone(sum_len_group+4,iket)=1;
			carbon_taken(sum_len_group+2)=1;
			//iket++;
			len_group=-3;
			previous="O=";
			if (debug_mode==1) cout << "ici 1 " << endl; 
		      }
		    /*
		    else if (surrogate[i].smile.substr(0,4)=="O=CG")
		      {
			if (debug_mode==1) cout << "aldesterG" << endl;
			//aldester=1;		     
			//carbon_ketone(ipar+1,iket)=1;
			//carbon_ketone(2,iket)=1;
			//iket++;
			//carbon_ester(ipar+1)=1;			  
			surrogate[i].groups[33]++;
			//surrogate[i].groups[31]++;
			carbon_taken(2)=2;
			carbon_taken(3)++;		       
			}*/
		    /*
		    else if (surrogate[i].smile.substr(0,5)=="O=CFC")
		      {
			//carbon_ester(sum_len_group+2)=1;
			surrogate[i].groups[33]+=1;
			//carbon_ketone(sum_len_group+4,iket)=1;
			carbon_taken(sum_len_group+2)=2;
			//iket++;
			len_group=4;
			if (debug_mode==1) cout << "ici 1 " << endl; 
		      }*/
		    else
		      carbon_taken(sum_len_group+2)=1;
		      

		    
		    if (surrogate[i].smile.substr(0,3)=="C=O")
		      carbon_taken(sum_len_group)++;
		  }
		if (surrogate[i].smile.substr(0,4)=="O=C=")
		  {
		    nc_double_save=0;
		    carbon_taken(2)=2;
		  }
	      }
	    else if (surrogate[i].smile.substr(0,3)=="OOC" and sum_len_group==0)	      
	      {
		//if (debug_mode==1) cout << "is hydroxyperoxide" << endl;
		len_group=2;
		will_be_hydroperoxide=1;
		/*
		if (surrogate[i].smile.substr(3,1)=="" or previous=="" or surrogate[i].smile.substr(0,4)=="OOC)")
		  {
		    //if (debug_mode==1) cout << "is terminal" <<endl;
		    surrogate[i].groups[42]+=1;
		  }
		else if  (surrogate[i].smile.substr(0,4)=="OOCC")
		  {
		    //if (debug_mode==1) cout << "has one oxygen" <<endl;
		    surrogate[i].groups[43]+=1;
		  }
		else
		  {
		    //if (debug_mode==1) cout << "has no oxygen" <<endl;
		    surrogate[i].groups[44]+=1;
		    }*/	      		
	      }
	    
	    else if (surrogate[i].smile.substr(0,14)=="[O-][N+](=O)OC" and sum_len_group==0)	      
	      {
		len_group=13;
		will_be_nitrate=1;      		
	      }
	    else if (surrogate[i].smile.substr(0,9)=="O=N(=O)OC" and sum_len_group==0)	      
	      {
		len_group=8;
		will_be_nitrate=1;      		
	      }
	    else if (surrogate[i].smile.substr(0,14)=="[O-][N+](=O)C" and sum_len_group==0)	      
	      {
		len_group=12;
		will_be_nitrite=1;      		
	      }
	    else if (surrogate[i].smile.substr(0,8)=="O=N(=O)K" and sum_len_group==0)	      
	      {
		len_group=-7;
		carbon_taken(6)++;
		carbon_taken(7)++;
		surrogate[i].groups[59]++;
	      }
	    else if (surrogate[i].smile.substr(0,8)=="O=N(=O)C" and sum_len_group==0)	      
	      {
		len_group=7;
		will_be_nitrite=1;	      
	      }
	    /*
	    else if (surrogate[i].smile.substr(0,9)=="O=N(=O)EK" and sum_len_group==0)	      
	      {
		len_group=-8;
		if (surrogate[i].smile.length()>8)
		  surrogate[i].groups[41]++;
		else
		  surrogate[i].groups[40]++;
		carbon_ketone(6,iket)=1;
		carbon_ketone(8,iket)=1;
		previous="C";
		carbon_taken(6)++;
		carbon_taken(7)=carbon_taken(7)+2;
		}	 */  
	    else if (surrogate[i].smile.substr(0,8)=="O=N(=O)E" and sum_len_group==0)	      
	      {
		len_group=-8;
		if (surrogate[i].smile.length()>8)
		  surrogate[i].groups[41]++;
		else
		  surrogate[i].groups[40]++;
		carbon_ketone(6,iket)=1;
		carbon_ketone(8,iket)=1;
		iket++;
		previous="reset";
		carbon_taken(6)++;
		carbon_taken(7)=carbon_taken(7)+2;
	      }	    
	    

	    /*
	      else if (surrogate[i].smile.substr(0,5)=="C(OO)")	      
	      {
	      if (debug_mode==1) cout << "is hydroperoxide" << endl;
	      len_group=5;
	      int nC=3;
	      if (previous!="")
	      nC=nC-1;

	      nC=nC-(carbon_cycle(sum_len_group,0)+carbon_cycle(sum_len_group,1)+carbon_cycle(sum_len_group,2)+carbon_cycle(sum_len_group,3)+carbon_cycle(sum_len_group,4));
	   
	      if (surrogate[i].smile.substr(0,6)=="C(OO)C")
	      nC=nC-1;
	      if (surrogate[i].smile.substr(0,6)=="C(OO)(")
	      nC=nC-2;
	      if (nC>=2)
	      {
	      if (debug_mode==1) cout << "hydroperoxyde with tho hydrogens" <<endl;
	      surrogate[i].groups[42]+=1;
	      }
	      else if  (nC==1)
	      {
	      if (debug_mode==1) cout << "hydroperoxyde with one hydrogen" <<endl;
	      surrogate[i].groups[43]+=1;
	      }
	      else
	      {
	      if (debug_mode==1) cout << "hydroperoxyde with no hydrogen" <<endl;
	      surrogate[i].groups[44]+=1;
	      }	      		
	      }*/
	    else if (surrogate[i].smile.substr(0,1)=="C" or surrogate[i].smile.substr(0,1)=="c") // or (surrogate[i].smile.substr(0,2)=="LK" and sum_len_group==0))
	      {
		
		int hydroxyperoxide=0;
		int nitrate=0;
		int nitrite=0;
		int ester=0;
		int arom=0;		
		if (will_be_ester==1)
		  {
		    will_be_ester=0;
		    if (debug_mode==1) cout << "ester from will be " << endl;
		    ester++;		   
		  }
		if (will_be_hydroperoxide==1)
		  {
		    will_be_hydroperoxide=0;
		    hydroxyperoxide++;
		  }
		if (will_be_nitrate==1)
		  {
		    will_be_nitrate=0;
		    nitrate++;
		  }
		if (will_be_nitrite==1)
		  {
		    will_be_nitrite=0;
		    nitrite++;
		  }
		if (carbon_ester(sum_len_group)==1)
		  {
		    if (debug_mode==1) cout << "ester detected " << endl;
		    ester++;
		  }

		if (carbon_arom(sum_len_group)==1)
		  arom=1;
		
		//if (debug_mode==1) cout << "is carbon" << endl;
		int nC;
		if (arom==0)
		  {
		    nC=4;	      
		    nC=nC-(carbon_cycle(sum_len_group,0)+carbon_cycle(sum_len_group,1)+carbon_cycle(sum_len_group,2)+carbon_cycle(sum_len_group,3)+carbon_cycle(sum_len_group,4));				    
		    if (previous!="")
		      nC=nC-1;

		    if (previous=="reset")
		      {
			if (debug_mode==1) cout << "reset" << endl;
			nC=nC-2;
		      }
		  }
		else
		  {
		    //cout << "Ncycle on carbon " << carbon_cycle(sum_len_group,0)+carbon_cycle(sum_len_group,1)+carbon_cycle(sum_len_group,2)+carbon_cycle(sum_len_group,3)+carbon_cycle(sum_len_group,4) << endl;
		    nC=2-(carbon_cycle(sum_len_group,0)+carbon_cycle(sum_len_group,1)+carbon_cycle(sum_len_group,2)+carbon_cycle(sum_len_group,3)+carbon_cycle(sum_len_group,4));
		    if (previous!="")
		      nC=nC-1;
		    if (carbon_arom(sum_len_group+1)==1)
		      nC=nC+1;
		    if (surrogate[i].smile.substr(0,3)=="C(C" and carbon_arom(sum_len_group+2)==1)
		      nC=nC+1;
		    if (surrogate[i].smile.substr(1,1)==")")
		      nC=nC+1;		    
		    if (linked_to_arom(sum_len_group)==1)
		      {
			//if (debug_mode==1) cout << "found" << endl;
			nC=0;
		      }

		    
		      
		    //nC=nC-(carbon_cycle(sum_len_group,0)+carbon_cycle(sum_len_group,1)+carbon_cycle(sum_len_group,2)+carbon_cycle(sum_len_group,3)+carbon_cycle(sum_len_group,4));
		    /*
		      if (previous!="")
		      if (carbon_arom(sum_len_group-1)<1)
		      nC=nC-1;*/
		  }
		if (debug_mode==1) cout << sum_len_group << " is carbon " << nC << " " << previous << endl;
		
		int ipar;
		int ipar2;
		int igr;
		int is_double=0;

		get_branch(debug_mode, config, surrogate, carbon_alcool, carbon_near_alcool, carbon_tail, carbon_hydrogen, carbon_nh, carbon_taken,
			   carbon_arom, linked_to_arom, carbon_near_alcool_tmp, carbon_near_alcool_save, carbon_near_alcool_bef,
			   carbon_anhydre, carbon_ester,
			   carbon_cycle, carbon_ether, carbon_ketone, carbon_peroxide,
			   arom, nitrate, nitrite, ester, hydroxyperoxide,
			   iket, ipero, iether,
			   i, ngr, sum_len_group, 0);
		
		if (surrogate[i].smile.substr(0,2)=="C(")
		  {
		    get_branch(debug_mode,config, surrogate, carbon_alcool, carbon_near_alcool, carbon_tail, carbon_hydrogen, carbon_nh, carbon_taken,
			       carbon_arom, linked_to_arom, carbon_near_alcool_tmp, carbon_near_alcool_save, carbon_near_alcool_bef,
			       carbon_anhydre, carbon_ester,
			       carbon_cycle, carbon_ether, carbon_ketone, carbon_peroxide,
			       arom, nitrate, nitrite, ester, hydroxyperoxide,
			       iket, ipero, iether,
			       i, ngr, sum_len_group, 1);
		    
		    if (surrogate[i].smile.substr(0,3)=="C(=")
		      {
			is_double=1;
			nC=nC+1; //To avoid double count a carbon bound
		      }
		    
		    nC=nC-1;
		    ipar=0;
		    for (igr=1;igr<int(surrogate[i].smile.length());igr++)
		      {			
			if (ipar==0 and surrogate[i].smile.substr(igr,1)==")")			
			  ipar=igr;
			else if (ipar<0 and surrogate[i].smile.substr(igr,1)==")")			
			  ipar=ipar+1;
			else if (ipar<=0 and surrogate[i].smile.substr(igr,1)=="(" and igr>2)
			  {
			    //if (debug_mode==1) cout << "une ramification s'est ouverte" << endl;
			    ipar=ipar-1;
			  }
		      }

		    get_branch(debug_mode,config, surrogate, carbon_alcool, carbon_near_alcool, carbon_tail, carbon_hydrogen, carbon_nh, carbon_taken,
			       carbon_arom, linked_to_arom, carbon_near_alcool_tmp, carbon_near_alcool_save, carbon_near_alcool_bef,
			       carbon_anhydre, carbon_ester,
			       carbon_cycle, carbon_ether, carbon_ketone, carbon_peroxide,
			       arom, nitrate, nitrite, ester, hydroxyperoxide,
			       iket, ipero, iether,
			       i, ngr, sum_len_group, ipar);
		    
		    if (surrogate[i].smile.substr(ipar+1,1)=="(")
		      {
			get_branch(debug_mode,config, surrogate, carbon_alcool, carbon_near_alcool, carbon_tail, carbon_hydrogen, carbon_nh, carbon_taken,
				   carbon_arom, linked_to_arom, carbon_near_alcool_tmp, carbon_near_alcool_save, carbon_near_alcool_bef,
				   carbon_anhydre, carbon_ester,
				   carbon_cycle, carbon_ether, carbon_ketone, carbon_peroxide,
				   arom, nitrate, nitrite, ester, hydroxyperoxide,
				   iket, ipero, iether,
				   i, ngr, sum_len_group, ipar+1);
			
			nC=nC-1;
			ipar2=0;
			for (igr=ipar+1;igr<int(surrogate[i].smile.length());igr++)
			  if (ipar2==0 and surrogate[i].smile.substr(igr,1)==")")			
			    ipar2=igr;
			  else if (ipar2<0 and surrogate[i].smile.substr(igr,1)==")")			
			    ipar2=ipar2+1;
			  else if (ipar2<=0 and surrogate[i].smile.substr(igr,1)=="(" and igr>ipar+1)
			    {
			      //if (debug_mode==1) cout << "une ramification s'est ouverte" << endl;
			      ipar2=ipar2-1;
			    }

			//if (debug_mode==1) cout << "I: " << ipar2 << " " << surrogate[i].smile.substr(ipar2+1,1) << " " << carbon_cycle(sum_len_group+ipar2+1,icycle) << endl;
			get_branch(debug_mode,config, surrogate, carbon_alcool, carbon_near_alcool, carbon_tail, carbon_hydrogen, carbon_nh, carbon_taken,
			       carbon_arom, linked_to_arom, carbon_near_alcool_tmp, carbon_near_alcool_save, carbon_near_alcool_bef,
			       carbon_anhydre, carbon_ester,
			       carbon_cycle, carbon_ether, carbon_ketone, carbon_peroxide,
			       arom, nitrate, nitrite, ester, hydroxyperoxide,
			       iket, ipero, iether,
			       i, ngr, sum_len_group, ipar2);
		    

			if (surrogate[i].smile.substr(ipar2+1,1)=="(")
			  {
			    get_branch(debug_mode,config, surrogate, carbon_alcool, carbon_near_alcool, carbon_tail, carbon_hydrogen, carbon_nh, carbon_taken,
				       carbon_arom, linked_to_arom, carbon_near_alcool_tmp, carbon_near_alcool_save, carbon_near_alcool_bef,
				       carbon_anhydre, carbon_ester,
				       carbon_cycle, carbon_ether, carbon_ketone, carbon_peroxide,
				       arom, nitrate, nitrite, ester, hydroxyperoxide,
				       iket, ipero, iether,
				       i, ngr, sum_len_group, ipar2+1);

			    nC=nC-1;
			    int ipar3=0;
			    for (igr=ipar2+1;igr<int(surrogate[i].smile.length());igr++)
			      if (ipar3==0 and surrogate[i].smile.substr(igr,1)==")")			
				ipar3=igr;
			      else if (ipar3<0 and surrogate[i].smile.substr(igr,1)==")")			
				ipar3=ipar3+1;
			      else if (ipar3<=0 and surrogate[i].smile.substr(igr,1)=="(" and igr>ipar2+1)
				{
				  //if (debug_mode==1) cout << "une ramification s'est ouverte" << endl;
				  ipar3=ipar3-1;
				}
			    
			    get_branch(debug_mode,config, surrogate, carbon_alcool, carbon_near_alcool, carbon_tail, carbon_hydrogen, carbon_nh, carbon_taken,
				       carbon_arom, linked_to_arom, carbon_near_alcool_tmp, carbon_near_alcool_save, carbon_near_alcool_bef,
				       carbon_anhydre, carbon_ester,
				       carbon_cycle, carbon_ether, carbon_ketone, carbon_peroxide,
				       arom, nitrate, nitrite, ester, hydroxyperoxide,
				       iket, ipero, iether,
				       i, ngr, sum_len_group, ipar3);

			    if (surrogate[i].smile.substr(ipar3+1,1)=="(")
			      { 
				if (outFlux==1 or outFlux==0)
				  {
				    cout << "Too much branches " << sum_len_group << " " << sum_len_group+ipar2+1 << endl;
				    cout << "problem in smiles.cxx. Try manual group decomposition" << endl;
				  }
				else if (outFlux==2)
				  {
				    cout << "Too much branches" << endl;
				    cout << "problem in smiles.cxx. Try manual group decomposition" << endl;
				    fileFlux << "encore un groupe" << endl;
				    fileFlux << "problem in smiles.cxx. Try manual group decomposition" << endl;				
				  }
				exit(0);
			      }
			  }

			if (surrogate[i].smile.substr(ipar2+1,1)=="=")
			  {
			    is_double=1;
			  }
			if (surrogate[i].smile.substr(ipar2+1,1)!="=" and surrogate[i].smile.substr(ipar2+1,1)!="" and surrogate[i].smile.substr(ipar2+1,1)!=")")
			  {
			    nC=nC-1;
			  }
		      
		      }
		    else
		      {
			if (surrogate[i].smile.substr(ipar+1,1)=="=")
			  {
			    if (debug_mode==1) cout << "double B" << endl;
			    is_double=1;
			  }
			if (surrogate[i].smile.substr(ipar+1,1)!="=" and surrogate[i].smile.substr(ipar+1,1)!="" and surrogate[i].smile.substr(ipar+1,1)!=")")
			  {
			    nC=nC-1;
			  }
		      }
		       
		  }
		else
		  {
		    if (surrogate[i].smile.substr(1,1)=="=")
		      {
			if (debug_mode==1) cout << "double A " << endl;
			is_double=1;
		      }
		    
		    if (surrogate[i].smile.substr(1,1)!="=" and surrogate[i].smile.substr(1,1)!="" and surrogate[i].smile.substr(1,1)!=")")
		      {
			nC=nC-1;
		      }

		  }

		//cout << nC << endl;
		if (is_double>0)
		  nC=nC-2;
		else if (will_be_double==1)
		  nC=nC-1;
		carbon_nh(sum_len_group)=nC;
		//cout << nC << endl;
		//exit(0);
		
		if (ester>0)
		  {
		    if (nC==3)
		      {
			if (debug_mode==1) cout << "ester with three hydrogens" <<endl;
			surrogate[i].groups[32]+=ester;
		      }
		    else 
		      {			
			if (debug_mode==1) cout << "ester with two or less hydrogens" <<endl;
			surrogate[i].groups[33]+=ester;
		      }
		    carbon_taken(sum_len_group)+=ester;
		  }		
		if (hydroxyperoxide>0)
		  {
		    if (nC>=2)
		      {
			//if (debug_mode==1) cout << "hydroperoxyde with tho hydrogens" <<endl;
			surrogate[i].groups[42]+=hydroxyperoxide;
		      }
		    else if  (nC==1)
		      {
			//if (debug_mode==1) cout << "hydroperoxyde with one hydrogen" <<endl;
			surrogate[i].groups[43]+=hydroxyperoxide;
		      }
		    else
		      {
			//if (debug_mode==1) cout << "hydroperoxyde with no hydrogen" <<endl;
			surrogate[i].groups[44]+=hydroxyperoxide;
		      }
		    
		    if (carbon_taken(sum_len_group)>0 or hydroxyperoxide>1)
		      {
			if (outFlux==1)
			  {
			    cout << "Waring: already occupied " <<  endl;
			  }
			else if (outFlux==2)
			  {
			    fileFlux << "Waring: already occupied" << endl;
			  }
		      }
		    
		    carbon_taken(sum_len_group)+=hydroxyperoxide;
		    
		  }
		if (nitrate>0)
		  {
		    //if (debug_mode==1) cout << "nitrate " << surrogate[i].smile << " " << nitrate << endl;
		    if (nC>=2)
		      {
			//if (debug_mode==1) cout << "nitrate with tho hydrogens" <<endl;
			surrogate[i].groups[39]+=nitrate;
		      }
		    else if  (nC==1)
		      {
			//if (debug_mode==1) cout << "nitrate with one hydrogen" <<endl;
			surrogate[i].groups[40]+=nitrate;
		      }
		    else
		      {
			//if (debug_mode==1) cout << "nitrate with no hydrogen" <<endl;
			surrogate[i].groups[41]+=nitrate;
		      }

		    if (carbon_taken(sum_len_group)>0 or nitrate>1)
		      {
			if (outFlux==1)
			  {
			    cout << "Waring: already occupied" << endl;
			  }
			else if (outFlux==2)
			  {
			    fileFlux << "Waring: already occupied" << endl;
			  }
		      }
		    
		    carbon_taken(sum_len_group)+=nitrate;
		    
		  }
		 if (nitrite>0)
		  {
		    if (carbon_arom(sum_len_group)>0)
		      {
			surrogate[i].groups[38]+=nitrite;
			carbon_arom(sum_len_group)=0;
		      }
		    else
		      {
			if (nC==3)
			  {
			    if (outFlux==1)			  
			      cout << "nitrite with three hydrogens" <<endl;
			    else if (outFlux==2)
			      fileFlux << "nitrite with three hydrogens" <<endl;
			    surrogate[i].groups[57]+=nitrite;
			  }
			else if  (nC==2)
			  {
			    if (outFlux==1)			  
			      cout << "nitrite with two hydrogens" <<endl;
			    else if (outFlux==2)
			      fileFlux << "nitrite with two hydrogens" <<endl;
			    surrogate[i].groups[58]+=nitrite;
			  }
			else
			  {
			    if (outFlux==1)
			      cout << "nitrite with one or zero hydrogen" <<endl;
			    else if (outFlux==2)
			      fileFlux << "nitrite with one or zero hydrogen" <<endl;
			    surrogate[i].groups[59]+=nitrite;
			  }

			if (carbon_taken(sum_len_group)>0 or nitrite>1)
			  {
			    if (outFlux==1)
			      cout << "Waring: already occupied" << endl;
			    else if (outFlux==2)
			      fileFlux << "Waring: already occupied" << endl;			  
			  }
		    
			carbon_taken(sum_len_group)+=nitrite;
		      }
		  }


		if (is_double==1 and nc_double_save>=0)
		  {
		    if (outFlux==1)
		      cout << "another double bound is currently treated" << endl;
		    else if (outFlux==2)
		      fileFlux << "another double bound is currently treated" << endl;
		    //break;
		  }
		else if (is_double==1)
		  {
		    nc_double_save=nC;
		    if (debug_mode==1) cout << "beg double with " << nC-2 << endl;
		    carbon_taken(sum_len_group)++;
		  }
		else if (will_be_double==1)
		  {
		    if (debug_mode==1) cout << "double bounds between " << nc_double_save <<  "hydrogen carbon and " << nC << " hydrogen carbon" << endl;

		    //nc_double_save=max(nc_double_save,0);
		    int nC1=max(nC,nc_double_save);
		    int nC2=min(nC,nc_double_save);		    

		    if (nC1==2 and nC2==1)
		      surrogate[i].groups[16]+=1;
		    else if (nC1==1 and nC2==1)
		      surrogate[i].groups[17]+=1;
		    else if (nC1==2 and nC2==0)
		      surrogate[i].groups[18]+=1;
		    else if (nC1==1 and nC2==0)
		      surrogate[i].groups[19]+=1;
		    else if (nC1==0 and nC2==0)
		      surrogate[i].groups[20]+=1;
		    // VIC C+C GECKO //
		    else if  (nC1==2 and nC2==2)
		      {
		        surrogate[i].groups[16]+=1;
			if (outFlux==1)
			  {
			    cout << "WARNING: CH2=CH2 not available" << endl;
			    cout << "Group CH2=CH used instead" << endl;
			  }
			else if (outFlux==2)
			  {
			    fileFlux <<"WARNING: CH2=CH2 not available" << endl;
			    fileFlux << "Group CH2=CH used instead" << endl;
			  }			   
                      }
		    ///////////////////  
		    else
		      {			
			cout << "error double bounds type " << nC1 << " " << nC2 << endl;
			cout << "problem in smiles.cxx. Try manual group decomposition" << endl;	
			if (outFlux==2)
			  {
			    fileFlux << "error double bounds type " << nC1 << " " << nC2 << endl;
			    fileFlux << "problem in smiles.cxx. Try manual group decomposition" << endl;
			  }
			exit(0);
		      }
		    will_be_double=0;
		    nc_double_save=-1;
		    previous="C=C";
		    carbon_taken(sum_len_group)++;
		  }
		else if (arom==1)
		  {
		    //if (debug_mode==1) cout << "is aromatic carbon with " << nC << " hydrogens" << endl;
		    previous="c";
		  }
		else if (hydroxyperoxide==0 and nitrate==0 and ester==0 and nitrite==0 and carbon_taken(sum_len_group)==0)
		  {
		    //if (debug_mode==1) cout << "is carbon with " << nC << " hydrogens" << endl;
		    previous="C";
		    carbon_hydrogen(sum_len_group)=nC;
		  }
		if (debug_mode==1) cout << "vla " << endl;
		len_group=1;
	      }
	    else if (surrogate[i].smile.substr(0,1)=="=")
	      {
		len_group=1;
		will_be_double=1;
	      }
	    else if (surrogate[i].smile.substr(0,4)=="(OO)")
	      {
		len_group=4;		
	      }
	    else if (surrogate[i].smile.substr(0,3)=="(O)")
	      {
		len_group=3;
		previous="O";
	      }
	    else if (surrogate[i].smile.substr(0,1)=="1" or surrogate[i].smile.substr(0,1)=="1" or surrogate[i].smile.substr(0,1)=="(" or surrogate[i].smile.substr(0,1)==")")
	      {
		len_group=1;
	      }
	    else if (surrogate[i].smile.substr(0,2)=="O")
	      {
		len_group=1;
		previous="O";
	      }
	    else if (surrogate[i].smile.substr(0,2)=="O)")
	      {
		len_group=2;
		previous="O";
	      }
	    else if (surrogate[i].smile.substr(0,3)=="OOK" or surrogate[i].smile.substr(0,3)=="OOA" or surrogate[i].smile.substr(0,3)=="OOL" or surrogate[i].smile.substr(0,3)=="OOG")
	      {
		len_group=-2;
		previous="reset";
	      }	  
	    else if (surrogate[i].smile.substr(0,2)=="OO")
	      {
		len_group=2;
		previous="OO";
	      }
	    else if (surrogate[i].smile.substr(0,3)=="OEK" or surrogate[i].smile.substr(0,3)=="OEL")
	      {
		len_group=-2;
		previous="OE";
	      }	  
	    else if (surrogate[i].smile.substr(0,2)=="OE")
	      {
		len_group=2;
		previous="OE";
	      }	  
	    else if (surrogate[i].smile.substr(0,8)=="ON(=O)=O")
	      {
		len_group=8;
	      }
	    else if (surrogate[i].smile.substr(0,7)=="N(=O)=O")
	      {
		len_group=7;
	      }
	    else if (surrogate[i].smile.substr(0,2)=="OC")
	      {
		len_group=1;
		previous="O";
	      }	    
	    else if (surrogate[i].smile.substr(0,6)=="KKKKON")
	      {			
		len_group=-4;
		previous="reset";
		carbon_taken(sum_len_group+1)++;
	      }
	    else if (surrogate[i].smile.substr(0,5)=="KKKON")
	      {			
		len_group=-3;
		previous="reset";
		carbon_taken(sum_len_group+1)++;
	      }
	    else if (surrogate[i].smile.substr(0,5)=="KKKKN")
	      {			
		len_group=-4;
		previous="reset";
		carbon_taken(sum_len_group+4)++;
	      }	    
	    else if (surrogate[i].smile.substr(0,4)=="KKKN")
	      {			
		len_group=-3;
		previous="reset";
		carbon_taken(sum_len_group+3)++;
	      }
	    else if (surrogate[i].smile.substr(0,3)=="KKN")
	      {			
		len_group=-2;
		previous="reset";
		carbon_taken(sum_len_group+2)++;
	      }
	    else if (surrogate[i].smile.substr(0,2)=="KN")
	      {			
		len_group=-1;
		previous="reset";
		carbon_taken(sum_len_group+1)++;
	      }
	    else if (surrogate[i].smile.substr(0,4)=="KKON")
	      {			
		len_group=-2;
		previous="reset";
		carbon_taken(sum_len_group+1)++;
	      }
	    else if (surrogate[i].smile.substr(0,3)=="KON")
	      {			
		len_group=-1;
		previous="reset";
		carbon_taken(sum_len_group+1)++;
	      }
	    else if (surrogate[i].smile.substr(0,6)=="KKKKOO" or surrogate[i].smile.substr(0,6)=="KKKKOE")
	      {			
		len_group=-4;
		previous="reset";
		
		//carbon_taken(sum_len_group+1)++;
	      }	
	    else if (surrogate[i].smile.substr(0,5)=="KKKOO" or surrogate[i].smile.substr(0,5)=="KKKOE")
	      {			
		len_group=-3;
		previous="reset";
		
		//carbon_taken(sum_len_group+1)++;
	      }	
	    else if (surrogate[i].smile.substr(0,4)=="KKOO" or surrogate[i].smile.substr(0,4)=="KKOE")
	      {			
		len_group=-2;
		previous="reset";
		
		//carbon_taken(sum_len_group+1)++;
	      }		    
	    else if (surrogate[i].smile.substr(0,3)=="KOO" or surrogate[i].smile.substr(0,3)=="KOE")
	      {			
		len_group=-1;
		previous="reset";
		
		//carbon_taken(sum_len_group+1)++;
	      }	    
	    else if (surrogate[i].smile.substr(0,4)=="KKKK")
	      {
		/*
		carbon_ketone(sum_len_group+1,iket)=1;
		carbon_ketone(sum_len_group+2,iket)=1;
		carbon_taken(sum_len_group+1)++;
		iket++;

		carbon_ketone(sum_len_group+2,iket)=1;
		carbon_ketone(sum_len_group+3,iket)=1;
		carbon_taken(sum_len_group+2)++;
		iket++;

		carbon_ketone(sum_len_group+3,iket)=1;
		carbon_ketone(sum_len_group+4,iket)=1;
		carbon_taken(sum_len_group+3)++;
		iket++;		*/		
		
		len_group=4;
		previous="C";
	      }
	    else if (surrogate[i].smile.substr(0,3)=="KKK")
	      {
		/*
		cout << "non..." << endl;
		carbon_ketone(sum_len_group+1,iket)=1;
		carbon_ketone(sum_len_group+2,iket)=1;
		carbon_taken(sum_len_group+1)++;
		iket++;

		carbon_ketone(sum_len_group+2,iket)=1;
		carbon_ketone(sum_len_group+3,iket)=1;
		carbon_taken(sum_len_group+2)++;
		iket++;*/
		
		len_group=3;
		previous="C";
	      }
	    else if (surrogate[i].smile.substr(0,2)=="KK")
	      {
		/*
		if (debug_mode==1) cout << "merde" << endl;
		carbon_ketone(sum_len_group+1,iket)=1;
		carbon_ketone(sum_len_group+2,iket)=1;
		carbon_taken(sum_len_group+1)++;
		iket++;*/
		len_group=2;
		previous="C";
	      }
	    else if (surrogate[i].smile.substr(0,2)=="EK" or surrogate[i].smile.substr(0,2)=="EL" or surrogate[i].smile.substr(0,2)=="EA" or surrogate[i].smile.substr(0,2)=="EP" or surrogate[i].smile.substr(0,2)=="EB")
	      {
		if (debug_mode==1) cout << "EK2 " << endl;
		if (carbon_anhydre(sum_len_group)<0)
		  {
		    //carbon_taken(sum_len_group+1)=1;
		    carbon_taken(sum_len_group)=1;
		    will_be_ester=1;
		  }
		/*
		carbon_ketone(sum_len_group+1,iket)=1;
		carbon_ketone(sum_len_group+2,iket)=1;
		carbon_taken(sum_len_group+1)=1;
		iket++;*/
		len_group=-1;
		previous="C";
	      }
	    else if (surrogate[i].smile.substr(0,4)=="EC=O")
	      {
		if (debug_mode==1) cout << "EKC=O " << endl;
		if (carbon_anhydre(sum_len_group)<0)
		  {
		    //carbon_taken(sum_len_group+1)=1;
		    carbon_taken(sum_len_group)=1;
		    will_be_ester=1;
		  }
		/*
		carbon_ketone(sum_len_group+1,iket)=1;
		carbon_ketone(sum_len_group+2,iket)=1;
		carbon_taken(sum_len_group+1)=1;
		iket++;*/
		len_group=-1;
		previous="C";
	      }
	    /*
	    else if (surrogate[i].smile.substr(0,3)=="GKO")
	      {
		if (debug_mode==1) cout << "GKO " << endl;
		len_group=-2;
		previous="C";
		}	  */  
	    else if (surrogate[i].smile.substr(0,2)=="GK")
	      {
		len_group=2;
		previous="C";
	      }
	    else if (surrogate[i].smile.substr(0,1)=="K")
	      {
		len_group=1;
		previous="C";
	      }
	    /*
	    else if (surrogate[i].smile.substr(0,4)=="EC=O")
	      {
		if (debug_mode==1) cout << "is aldehyde ester" << endl;
		surrogate[i].groups[31]+=1;
		len_group=4;
		carbon_ester(sum_len_group+1,0)=1;
	      }*/
	    else if (surrogate[i].smile.substr(0,1)=="E")
	      {
		if (carbon_anhydre(sum_len_group)<0)
		  {
		    will_be_ester=1;
		  }
		len_group=1;
		previous="C";
	      }
	    //if (debug_mode==1) cout << "len2: " << len_group << endl;
	    if (surrogate[i].smile.substr(0,3)=="LOO")
	      {
		/*
		if (debug_mode==1) cout << "LOO " << endl;
		len_group=2;
		carbon_ester(sum_len_group)=1;
		int icycle;
		for (icycle=0;icycle<5;icycle++)
		  if (carbon_cycle(sum_len_group,icycle)>0)
		    for (igr=0;igr<ngr;igr++)
		      if (igr!=sum_len_group and carbon_cycle(igr,icycle)>0)			 
			carbon_ester(igr)=1;
		//if (debug_mode==1) cout << "carbon ester " << carbon_ester << endl;
		//iether++;*/
		carbon_nh(sum_len_group)=0;
		carbon_taken(sum_len_group+1)++;
		carbon_taken(sum_len_group+2)++;
		len_group=-1;
		previous="reset";
	      }
	    else if (surrogate[i].smile.substr(0,2)=="LO")
	      {
		if (debug_mode==1) cout << "LO " << endl;
		len_group=2;
		carbon_ester(sum_len_group)=1;
		
		for (icycle=0;icycle<5;icycle++)
		  if (carbon_cycle(sum_len_group,icycle)>0)
		    for (igr=0;igr<ngr;igr++)
		      if (igr!=sum_len_group and carbon_cycle(igr,icycle)>0)			 
			carbon_ester(igr)=1;
		//if (debug_mode==1) cout << "carbon ester " << carbon_ester << endl;
		//iether++;
		previous="C";
	      }
	    /*
	    else if (surrogate[i].smile.substr(0,2)=="LO")
	      {
		will_be_ester=1;
		len_group=2;
		previous="C";
		}	 */
	    else if (surrogate[i].smile.substr(0,1)=="GCOO")
	      {	       
		len_group=1;
		previous="C";		
	      }
	    else if (surrogate[i].smile.substr(0,1)=="G" and surrogate[i].smile.substr(0,2)!="GL" and surrogate[i].smile.substr(0,2)!="GK" and surrogate[i].smile.substr(0,4)!="GC=O")
	      {
		if (debug_mode==1) cout << "G here" << endl; 
		len_group=1;
		previous="C";		
	      }
	    else if (surrogate[i].smile.substr(0,1)=="H")
	      {
		//if (debug_mode==1) cout << "F here" << endl; 
		len_group=1;
		previous="C";		
	      }
	    else if (surrogate[i].smile.substr(0,2)=="OX")
	      {		
		len_group=2;
		previous="C";
	      }
	    else if (surrogate[i].smile.substr(0,2)=="OA")
	      {		
		len_group=1;
		previous="O";
	      }
	    else if (surrogate[i].smile.substr(0,1)=="X")
	      {		
		len_group=1;
		previous="C";
	      }
	    else if (surrogate[i].smile.substr(0,13)=="O[N+](=O)[O-]")
	      {
		len_group=13;
		previous="C";
	      }
	    else if (surrogate[i].smile.substr(0,12)=="[N+](=O)[O-]")
	      {
		len_group=12;
		previous="C";
	      }
	    

	    if (debug_mode==1) cout << "len: " << len_group << endl;
	    if (len_group<=-2)
	      {
		string smile3="C"+surrogate[i].smile.substr(abs(len_group),ngr-abs(len_group));
		if (debug_mode==1) cout << "update2: " << smile3 << endl;
		surrogate[i].smile=smile3;
		sum_len_group=sum_len_group+abs(len_group)-1;
	      }
	    else if (len_group==-1)
	      {
		string smile3="C"+surrogate[i].smile.substr(1,ngr-1);
		if (debug_mode==1) cout << "update: " << smile3 << endl;
		surrogate[i].smile=smile3;
	      }
	    else
	      {
		sum_len_group+=len_group;
		previous=surrogate[i].smile.substr(0,len_group);
		surrogate[i].smile=surrogate[i].smile.substr(len_group,ngr-len_group);
		//if (debug_mode==1) cout << len_group << " " << surrogate[i].smile << endl;

		if (len_group==0)
		  {		
		    cout << surrogate[i].smile.substr(0,5) << " " << sum_len_group << endl;
		    cout << "group not found" << endl;
		    cout << "problem in smiles.cxx. Try manual group decomposition" << endl;
		    if (outFlux==2)
		      {
			fileFlux << surrogate[i].smile.substr(0,5) << endl;
			fileFlux << "group not found" << endl;
			fileFlux << "problem in smiles.cxx. Try manual group decomposition" << endl;
		      }
		    exit(0);
		  }
	      }
	  }

	surrogate[i].smile=smile2;
	ngr=surrogate[i].smile.length();
	if (debug_mode==1) cout << smile2 << " ngr: " << ngr << endl;
	int igr3;
	//cout << ngr << " ether " << endl;
	//cout << carbon_ether << endl;	

	//cout << "molecules has " << sum(carbon_alcool) << " alcool group" << endl;
	surrogate[i].groups[26]+=sum(carbon_alcool);
	
	for (igr=0;igr<ngr;igr++)
	  if (carbon_arom(igr)>0)
	    carbon_taken(igr)++;

	
	if (debug_mode==1)
	  {
	    cout << "Taken carbons: " << endl;
	    for (igr=0;igr<ngr;igr++)
	      cout <<  surrogate[i].smile.substr(igr,1) << " " << carbon_taken(igr) << endl;
	    cout << " " << endl;
	  }
	
	//if (debug_mode==1) cout << "ket " << endl;
	if (debug_mode==1) cout << carbon_ketone << endl;
	//if (debug_mode==1) cout << carbon_taken << endl;
	//if (debug_mode==1) cout << "ester " << endl;
	if (debug_mode==1) cout << carbon_ester << endl;
	//if (debug_mode==1) cout << "ketone " << endl;

	//if (debug_mode==1) cout << carbon_taken << endl;
	/*if (debug_mode==1) cout << carbon_ketone << endl;
	
	for (igr2=0;igr2<ngr_max;igr2++)
	  for (igr=0;igr<ngr;igr++)
	    if (carbon_ketone(igr,igr2)>=0)
	    if (debug_mode==1) cout << igr2 << " " << igr << " " << carbon_ketone(igr,igr2) << endl; //" " << surrogate[i].smile.substr(igr-1,3) << endl; */
	
	//Treatment of peroxide molecules
	if (debug_mode==1) cout << carbon_peroxide << endl;
	for (igr2=0;igr2<ngr_max;igr2++)
	  {
	    int i1=-1;
	    int i2=-1;	    
	    for (igr=0;igr<ngr;igr++)
	      if (carbon_peroxide(igr,igr2)>=0 and i1<0)
		i1=igr;
	      else if (carbon_peroxide(igr,igr2)>=0)
		i2=igr;

	    if (debug_mode==1) cout << i1 << " " << i2 << endl;
	    if (i1>=0 and i2>=0)
	      {
		if (debug_mode==1) cout << "peroxide " << i1 << " " << i2 << endl;
		//if (debug_mode==1) cout << "taken " << carbon_taken(i1) << " " << carbon_taken(i2) << endl;
		if (carbon_taken(i1)>0 or carbon_taken(i2)>0)
		  {
		    if (outFlux==1) 
		      cout << "WARNING: already occupied carbon" << endl;
		    else if (outFlux==2)
		      fileFlux << "WARNING: already occupied carbon" << endl;
		  }
		carbon_taken(i1)++;	      
		carbon_taken(i2)++;

		int nc1=max(max(carbon_nh(i1),carbon_nh(i2)),0);
		int nc2=max(min(carbon_nh(i1),carbon_nh(i2)),0);

		//cout << i1 << " " << nc1 << endl;
		//cout << i2 << " " << nc2 << endl;

		carbon_hydrogen(i1)=-1;
		carbon_hydrogen(i2)=-1;

		if (nc1==3 and nc2==2)		  
		  surrogate[i].groups[45]++;	       
		else if (nc1==3 and nc2==1)		  
		  surrogate[i].groups[46]++;
		else if (nc1==3 and nc2==0)		  
		  surrogate[i].groups[47]++;
		else if (nc1==2 and nc2==2)		  
		  surrogate[i].groups[48]++;
		else if (nc1==2 and nc2==1)		  
		  surrogate[i].groups[49]++;
		else if (nc1==2 and nc2==0)		  
		  surrogate[i].groups[50]++;
		else if (nc1==1 and nc2==1)		  
		  surrogate[i].groups[51]++;
		else if (nc1==1 and nc2==0)		  
		  surrogate[i].groups[52]++;
		else if (nc1==0 and nc2==0)		  
		  surrogate[i].groups[53]++;
		else
		  surrogate[i].groups[53]++;
	      }
	  }

	if (debug_mode==1) cout << "before ketone treatment: " << endl;
	for (igr=0;igr<ngr;igr++)
	  {
	    int iket2,iket3;
	    iket3=-1;
	    for (iket2=0;iket2<ngr_max;iket2++)
	      if (carbon_ketone(igr,iket2)>0)
		iket3=iket2;
	    if (debug_mode==1 and iket3>=0) cout << smile_save.substr(igr,1) << " " << carbon_nh(igr) << " et " << carbon_ketone(igr,iket3) << " " << carbon_taken(igr) << " " << iket3 << endl;
	  }
	if (debug_mode==1) cout << " " << endl;


	//Treatment of ketone molecules
	for (igr2=0;igr2<ngr_max;igr2++)
	  {
	    int i1=-1;
	    int i2=-1;
	    int occupied1=0;
	    int occupied2=0;
	    int has_another_ketone1=0;
	    int has_another_ketone2=0;
	    int ketone1=0;
	    int ketone2=0;
	    for (igr=0;igr<ngr;igr++)
	      if (carbon_ketone(igr,igr2)>=0 and i1<0)
		i1=igr;
	      else if (carbon_ketone(igr,igr2)>=0)
		i2=igr;

	    //if (debug_mode==1) cout << i1 << " " << i2 << endl;
	    if (i1>=0 and i2>=0)
	      {
		//if (debug_mode==1) cout << "ketone " << i1 << " " << i2 << endl;
		//if (debug_mode==1) cout << "taken " << carbon_taken(i1) << " " << carbon_taken(i2) << endl;
		if (carbon_taken(i1)>0)
		  occupied1=1;
		if (carbon_taken(i2)>0)
		  occupied2=1;
	    
		for (igr3=0;igr<ngr_max;igr++)
		  {
		    if (carbon_ketone(i1,igr3)>=0 and igr3!=igr2)
		      has_another_ketone1++;

		    if (carbon_ketone(i1,igr3)>=0)
		      has_another_ketone1++;
		  }
	    
		for (igr3=0;igr<ngr_max;igr++)
		  {
		    if (carbon_ketone(i1,igr3)>=0 and igr3!=igr2)		    
		      has_another_ketone2++;
		    
		    if (carbon_ketone(i1,igr3)>=0)
		      has_another_ketone1++;
		  }
	    

		if (occupied1==1 and occupied2==0)
		  ketone2=1;
		else if (occupied2==1 and occupied1==0)
		  ketone1=1;
		else if (occupied1==0 and occupied2==0)
		  {
		    if (has_another_ketone1>0 and has_another_ketone2==0)
		      ketone2=1;
		    else if (has_another_ketone1==0 and has_another_ketone2>0)
		      ketone1=1;
		    else if (has_another_ketone1==0 and has_another_ketone2==0)
		      {
			if (carbon_nh(i1)>=carbon_nh(i2))
			  ketone1=1;
			else
			  ketone2=1;
		      }
		    else if (has_another_ketone1>0 and has_another_ketone2>0)
		      {
			if (outFlux==1)
			  cout << "Warning: subsequent ketone" << endl;
			else if (outFlux==2)
			  fileFlux << "Warning: subsequent ketone" << endl;
			if (carbon_nh(i1)>=carbon_nh(i2))
			  ketone1=1;
			else
			  ketone2=1;
		      }
		  }
		else
		  {
		    if (outFlux==1)
		      cout << "Warning: already occupied ketone" << endl;
		    else if (outFlux==2)
		      fileFlux << "Warning: already occupied ketone" << endl;
		    if (carbon_nh(i1)>=carbon_nh(i2))
		      ketone1=1;
		    else
		      ketone2=1;
		  }

		int iketone=0;
		if (ketone1>0)	      
		  iketone=i1;
	    
		if (ketone2>0)	      
		  iketone=i2;	      

		carbon_taken(iketone)++;
		carbon_hydrogen(iketone)=-1;
		if (carbon_nh(iketone)==3)
		  {
		    if (debug_mode==1) cout << "ketone with 3 hydrogens" << endl;
		    surrogate[i].groups[29]++;
		  }
		else 
		  {
		    if (debug_mode==1) cout << "ketone with 2 or less hydrogens" << endl;
		    surrogate[i].groups[30]++;
		  }
	      }
	  }

	//if (debug_mode==1) cout << carbon_taken << endl;
     	
	//Treatment of ether molecules
	for (igr2=0;igr2<ngr_max;igr2++)
	  {
	    int i1=-1;
	    int i2=-1;
	    int occupied1=0;
	    int occupied2=0;
	    int has_another_ether1=0;
	    int has_another_ether2=0;
	    int ether1=0;
	    int ether2=0;
	    for (igr=0;igr<ngr;igr++)
	      if (carbon_ether(igr,igr2)>=0 and i1<0)
		i1=igr;
	      else if (carbon_ether(igr,igr2)>=0)
		i2=igr;

	    //if (debug_mode==1) cout << i1 << " " << i2 << endl;
	    if (i1>=0 and i2>=0)
	      {
		//if (debug_mode==1) cout << "ether " << i1 << " " << i2 << endl;
		if (carbon_taken(i1)>0)
		  occupied1=1;
		if (carbon_taken(i2)>0)
		  occupied2=1;
	    
		for (igr3=0;igr<ngr_max;igr++)
		  {
		    if (carbon_ether(i1,igr3)>=0 and igr3!=igr2)
		      has_another_ether1++;

		    if (carbon_ketone(i1,igr3)>=0)
		      has_another_ether1++;
		  }
	    
		for (igr3=0;igr<ngr_max;igr++)
		  {
		    if (carbon_ether(i1,igr3)>=0 and igr3!=igr2)		    
		      has_another_ether2++;
		    
		    if (carbon_ketone(i1,igr3)>=0)
		      has_another_ether1++;
		  }
	    

		if (occupied1==1 and occupied2==0)
		  ether2=1;
		else if (occupied2==1 and occupied1==0)
		  ether1=1;
		else if (occupied1==0 and occupied2==0)
		  {
		    if (has_another_ether1>0 and has_another_ether2==0)
		      ether2=1;
		    else if (has_another_ether1==0 and has_another_ether2>0)
		      ether1=1;
		    else if (has_another_ether1==0 and has_another_ether2==0)
		      {
			if (carbon_nh(i1)>=carbon_nh(i2))
			  ether1=1;
			else
			  ether2=1;
		      }
		    else if (has_another_ether1>0 and has_another_ether2>0)
		      {
			//if (debug_mode==1) cout << "Warning: subsequent ether" << endl;
			if (carbon_nh(i1)>=carbon_nh(i2))
			  ether1=1;
			else
			  ether2=1;
		      }
		  }
		else
		  {
		    //if (debug_mode==1) cout << "Warning: already occupied ether" << endl;
		    if (carbon_nh(i1)>=carbon_nh(i2))
		      ether1=1;
		    else
		      ether2=1;
		  }

		iether=0;
		if (ether1>0)	      
		  iether=i1;
	    
		if (ether2>0)	      
		  iether=i2;	      

		carbon_taken(iether)++;
		carbon_hydrogen(iether)=-1;
		if (carbon_nh(iether)==3)
		  {
		    //if (debug_mode==1) cout << "ether with 3 hydrogens" << endl;
		    surrogate[i].groups[34]++;
		  }
		else if (carbon_nh(iether)==2)
		  {
		    //if (debug_mode==1) cout << "ether with 2 hydrogens" << endl;
		    surrogate[i].groups[35]++;
		  }
		else
		  {	   
		    //if (debug_mode==1) cout << "ether with 0 or 1 hydrogen" << endl;
		    surrogate[i].groups[36]++;	      
		  }
	      }
	  }


	//if (debug_mode==1) cout << linked_to_arom << endl;
	//if (debug_mode==1) cout << carbon_hydrogen << endl;
	//if (debug_mode==1) cout << carbon_arom << endl;
	//Treatment of arom-CHx groups
	for (igr=0;igr<ngr;igr++)
	  if (carbon_hydrogen(igr)>=0 and linked_to_arom(igr)==1 and carbon_taken(igr)==0)	    
	    {
	      if (carbon_nh(igr)==3)
		{
		  //if (debug_mode==1) cout << "arom-CH3" << endl;
		  surrogate[i].groups[23]++;
		  surrogate[i].groups[22]--;
		}
	      else if (carbon_nh(igr)==2)
		{
		  //if (debug_mode==1) cout << "arom-CH2" << endl;
		  surrogate[i].groups[24]++;
		  surrogate[i].groups[22]--;
		}
	      else
		{	   
		  //if (debug_mode==1) cout << "arom-CH" << endl;
		  surrogate[i].groups[25]++;
		  surrogate[i].groups[22]--;
		}
	      carbon_hydrogen(igr)=-1;
	      carbon_taken(igr)=1;
	    }

	//cout << carbon_nh << " " << carbon_arom << endl;
	//exit(0);
	for (igr=0;igr<ngr;igr++)
	  if (carbon_arom(igr)==1)
	    {
	      //if (debug_mode==1) cout << "ici" << endl;
	      if (carbon_nh(igr)==1)
		{
		  if (debug_mode==1) cout << "aromatic carbon with one hydrogen" << " " << igr << endl;		 
		  surrogate[i].groups[21]++;
		  //exit(0);
		}
	      else 
		{
		  if (debug_mode==1) cout << "aromatic carbon with no hydrogen" << " " << igr << endl;	     
		  surrogate[i].groups[22]++;
		}
	    }

	//for (igr3=0;igr3<ngr;igr3++)
	for (igr=0;igr<ngr;igr++)
	  {
	    if (carbon_hydrogen(igr)>=0)
	      if (surrogate[i].smile.substr(igr,2)=="C" or surrogate[i].smile.substr(igr,2)=="C)" or igr==1)
		for (icycle=0;icycle<5;icycle++)
		  if (carbon_cycle(igr,icycle)==0)
		    carbon_tail(igr)=1;	    
	  }

	
		
	for (igr3=0;igr3<ngr;igr3++)
	  {
	    for (igr=0;igr<ngr;igr++)
	      {
		if (carbon_hydrogen(igr)>=0)
		  {
		    //if (debug_mode==1) cout << surrogate[i].smile.substr(igr,2) << endl;
		    if (surrogate[i].smile.substr(igr+1,1)=="C" and carbon_hydrogen(igr+1)>=0)
		      {
			if (carbon_alcool(igr+1)>0 or carbon_near_alcool_tmp(igr+1)>0)		      
			  if (carbon_near_alcool_save(igr)==0) carbon_near_alcool_save(igr)=igr3;
			if (carbon_alcool(igr)>0 or carbon_near_alcool_tmp(igr)>0)
			  if (carbon_near_alcool_save(igr+1)==0) carbon_near_alcool_save(igr+1)=igr3;
		      }

		    for (icycle=0;icycle<5;icycle++)
		      if (carbon_cycle(igr,icycle)>0)
			{
			  int found=0;
			  for (igr2=0;igr2<ngr;igr2++)
			    if (carbon_cycle(igr2,icycle)>0 and igr2!=igr)
			      if (found==0)
				{
				  found=1;
				  if (carbon_alcool(igr2)>0 or carbon_near_alcool_tmp(igr2)>0)
				    if (carbon_near_alcool_save(igr)==0) carbon_near_alcool_save(igr)=igr3;
				  if (carbon_alcool(igr)>0 or carbon_near_alcool_tmp(igr)>0)
				    if (carbon_near_alcool_save(igr2)==0) carbon_near_alcool_save(igr2)=igr3;
				}
			}

		    if (carbon_hydrogen(igr+1)>0)
		      carbon_near_alcool_bef(igr+1)=carbon_near_alcool_tmp(igr);

		    if (surrogate[i].smile.substr(igr+1,1)=="(")
		      {
			if (carbon_hydrogen(igr+2)>=0)
			  {
			    /*
			      if (carbon_near_alcool_tmp(igr+2)>0 and carbon_near_alcool(igr)>0)
			      {
			      if (debug_mode==1) cout << "the two branch has an alcool" << endl;
			      carbon_near_alcool(igr)=1;
			      carbon_near_alcool(igr)=1;
			      }*/
			
			    if (carbon_alcool(igr+2)>0 or carbon_near_alcool_tmp(igr+2)>0)
			      if (carbon_near_alcool_save(igr)==0) carbon_near_alcool_save(igr)=igr3;
			    if (carbon_alcool(igr)>0 or carbon_near_alcool_tmp(igr)>0)
			      if (carbon_near_alcool_save(igr+2)==0) carbon_near_alcool_save(igr+2)=igr3;
			  }

			carbon_near_alcool_bef(igr+2)=carbon_near_alcool_tmp(igr);
		    
			int ipar=0;
			for (igr2=igr+2;igr2<ngr;igr2++)
			  if (ipar==0 and surrogate[i].smile.substr(igr2,1)==")")			
			    ipar=igr2;
			  else if (ipar<0 and surrogate[i].smile.substr(igr2,1)==")")			
			    ipar=ipar+1;
			  else if (ipar<=0 and surrogate[i].smile.substr(igr2,1)=="(" and igr2>igr+2)
			    {
			      //if (debug_mode==1) cout << "une ramification s'est ouverte" << endl;
			      ipar=ipar-1;
			    }
			

			if (carbon_hydrogen(ipar+1)>=0)
			  {
			    if (carbon_alcool(ipar+1)>0 or carbon_near_alcool_tmp(ipar+1)>0)
			      if (carbon_near_alcool_save(igr)==0) carbon_near_alcool_save(igr)=igr3;
			    if (carbon_alcool(igr)>0 or carbon_near_alcool_tmp(igr)>0)
			      if (carbon_near_alcool_save(ipar+1)==0) carbon_near_alcool_save(ipar+1)=igr3;

			    carbon_near_alcool_bef(ipar+1)=carbon_near_alcool_tmp(igr);
			  }
		    
			if (surrogate[i].smile.substr(ipar+1,1)=="(")
			  {
			    if (carbon_hydrogen(ipar+2)>=0)
			      {
				if (carbon_alcool(ipar+2)>0 or carbon_near_alcool_tmp(ipar+2)>0)
				  if (carbon_near_alcool_save(igr)==0) carbon_near_alcool_save(igr)=igr3;
				if (carbon_alcool(igr)>0 or carbon_near_alcool_tmp(igr)>0)
				  if (carbon_near_alcool_save(ipar+2)==0)carbon_near_alcool_save(ipar+2)=igr3;

				carbon_near_alcool_bef(ipar+2)=carbon_near_alcool_tmp(igr);
			      }
					     
			    int ipar2=0;
			    for (igr2=ipar+1;igr2<ngr;igr2++)
			      if (ipar2==0 and surrogate[i].smile.substr(igr2,1)==")")			
				ipar2=igr2;
			      else if (ipar2<0 and surrogate[i].smile.substr(igr2,1)==")")			
				ipar2=ipar2+1;
			      else if (ipar2<=0 and surrogate[i].smile.substr(igr2,1)=="(" and igr2>ipar+1)
				{
				  //if (debug_mode==1) cout << "une ramification s'est ouverte" << endl;
				  ipar2=ipar2-1;
				}
			
			    if (carbon_hydrogen(ipar2+1)>=0)
			      {
				if (carbon_alcool(ipar2+1)>0 or carbon_near_alcool_tmp(ipar2+1)>0)
				  if (carbon_near_alcool_save(igr)==0) carbon_near_alcool_save(igr)=igr3;
				if (carbon_alcool(igr)>0 or carbon_near_alcool_tmp(igr)>0)
				  if (carbon_near_alcool_save(ipar2+1)==0) carbon_near_alcool_save(ipar2+1)=igr3;

				carbon_near_alcool_bef(ipar2+1)=carbon_near_alcool_tmp(igr);
			      }
			       		    
			  }
		    
		      }
	    

		
		  }
	      }
	    
	    for (igr=0;igr<ngr;igr++)
	      carbon_near_alcool_tmp(igr)=carbon_near_alcool_save(igr);
	  }

	for (igr=0;igr<ngr;igr++)
	  if (carbon_alcool(igr)>0)
	    carbon_near_alcool_tmp(igr)=0;
	  else if (carbon_near_alcool_tmp(igr)==0.0)
	    carbon_near_alcool_tmp(igr)=-1;
	 
	

	//if (debug_mode==1) cout << carbon_near_alcool_tmp << endl;

	int near_alcool_save=0;
	int alcool_save=0;
	for (igr=0;igr<ngr;igr++)
	  {
	    near_alcool_save=carbon_near_alcool_bef(igr);
	    //if (debug_mode==1) cout << igr << " " << near_alcool_save << " " << carbon_near_alcool_tmp(igr) << endl;
	    if (carbon_near_alcool_tmp(igr)>=0)
	      {
		//if (debug_mode==1) cout << surrogate[i].smile.substr(igr,2) << endl;
		if (carbon_near_alcool_tmp(igr+1)==carbon_near_alcool_tmp(igr))
		  {
		    carbon_near_alcool(igr)=1;
		    carbon_near_alcool(igr+1)=1;
		  }
		else if (near_alcool_save>=0 and carbon_near_alcool_tmp(igr+1)>=0)		  
		  if (near_alcool_save==carbon_near_alcool_tmp(igr)-1 and carbon_near_alcool_tmp(igr+1)==carbon_near_alcool_tmp(igr)-1)
		    {
		      carbon_near_alcool(igr)=1;		     
		    }

		for (icycle=0;icycle<5;icycle++)
		  if (carbon_cycle(igr,icycle)>0)
		    {
		      int found=0;
		      for (igr2=0;igr2<ngr;igr2++)
			if (carbon_cycle(igr2,icycle)>0 and igr2!=igr)
			  if (found==0)
			    {
			      found=1;
			      if (carbon_near_alcool_tmp(igr2)==carbon_near_alcool_tmp(igr))
				{
				  carbon_near_alcool(igr)=1;
				  carbon_near_alcool(igr2)=1;
				}
			      else if (near_alcool_save>=0 and carbon_near_alcool_tmp(igr2)>=0)		  
				if (near_alcool_save==carbon_near_alcool_tmp(igr)-1 and carbon_near_alcool_tmp(igr2)==carbon_near_alcool_tmp(igr)-1)
				  {
				    carbon_near_alcool(igr)=1;		     
				  }
			    }
		    }
		  	       
		if (surrogate[i].smile.substr(igr+1,1)=="(")
		  {		    		    
		    int ipar=0;
		    for (igr2=igr+2;igr2<ngr;igr2++)
		      if (ipar==0 and surrogate[i].smile.substr(igr2,1)==")")			
			ipar=igr2;
		      else if (ipar<0 and surrogate[i].smile.substr(igr2,1)==")")			
			ipar=ipar+1;
		      else if (ipar<=0 and surrogate[i].smile.substr(igr2,1)=="(" and igr2>igr+2)
			{
			  //if (debug_mode==1) cout << "une ramification s'est ouverte" << endl;
			  ipar=ipar-1;
			}

		    if (carbon_near_alcool_tmp(igr+2)==carbon_near_alcool_tmp(igr))
		      {
			carbon_near_alcool(igr)=1;
			carbon_near_alcool(igr+2)=1;
		      }
		    else if (near_alcool_save>=0 and carbon_near_alcool_tmp(igr+2)>=0)		  
		      if (near_alcool_save==carbon_near_alcool_tmp(igr)-1 and carbon_near_alcool_tmp(igr+2)==carbon_near_alcool_tmp(igr)-1)			
			carbon_near_alcool(igr)=1;
		    
		    if (carbon_near_alcool_tmp(ipar+1)==carbon_near_alcool_tmp(igr))
		      {
			carbon_near_alcool(igr)=1;
			carbon_near_alcool(ipar+1)=1;
		      }
		    else if (near_alcool_save>=0 and carbon_near_alcool_tmp(ipar+1)>=0)		  
		      if (near_alcool_save==carbon_near_alcool_tmp(igr)-1 and carbon_near_alcool_tmp(ipar+1)==carbon_near_alcool_tmp(igr)-1)			
			carbon_near_alcool(igr)=1;
		    
		    if (surrogate[i].smile.substr(ipar+1,1)=="(")
		      {
			if (carbon_near_alcool_tmp(ipar+2)==carbon_near_alcool_tmp(igr))
			  {
			    carbon_near_alcool(igr)=1;
			    carbon_near_alcool(ipar+2)=1;
			  }
			else if (near_alcool_save>=0 and carbon_near_alcool_tmp(ipar+2)>=0)		  
			  if (near_alcool_save==carbon_near_alcool_tmp(igr)-1 and carbon_near_alcool_tmp(ipar+2)==carbon_near_alcool_tmp(igr)-1)			
			    carbon_near_alcool(igr)=1;
										     
			int ipar2=0;
			for (igr2=ipar+1;igr2<ngr;igr2++)
			  if (ipar2==0 and surrogate[i].smile.substr(igr2,1)==")")			
			    ipar2=igr2;
			  else if (ipar2<0 and surrogate[i].smile.substr(igr2,1)==")")			
			    ipar2=ipar2+1;
			  else if (ipar2<=0 and surrogate[i].smile.substr(igr2,1)=="(" and igr2>ipar+1)
			    {
			      //if (debug_mode==1) cout << "une ramification s'est ouverte" << endl;
			      ipar2=ipar2-1;
			    }
			if (carbon_near_alcool_tmp(ipar2+1)==carbon_near_alcool_tmp(igr))
			  {
			    carbon_near_alcool(igr)=1;
			    carbon_near_alcool(ipar2+1)=1;
			  }
			else if (near_alcool_save>=0 and carbon_near_alcool_tmp(ipar2+1)>=0)		  
			  if (near_alcool_save==carbon_near_alcool_tmp(igr)-1 and carbon_near_alcool_tmp(ipar2+1)==carbon_near_alcool_tmp(igr)-1)			
			    carbon_near_alcool(igr)=1;					       
			
		      }
		    
		  }
	    
	    
	      }	 
	  }

	//if (debug_mode==1) cout << "after near:" << endl;
	//if (debug_mode==1) cout << carbon_near_alcool << endl;
	//if (debug_mode==1) cout << "la :" << carbon_near_alcool_tmp << endl;
	//if (debug_mode==1) cout << carbon_near_alcool << endl;
	if (debug_mode==1) cout << carbon_taken << endl;

	for (igr3=0;igr3<ngr;igr3++)
	for (igr=0;igr<ngr;igr++)
	  {
	    if (carbon_near_alcool_tmp(igr)>0)
	      {
		//if (debug_mode==1) cout << surrogate[i].smile.substr(igr,2) << endl;
		if (carbon_near_alcool(igr+1)==1 and carbon_near_alcool_tmp(igr)<carbon_near_alcool_tmp(igr+1))
		  carbon_near_alcool(igr)=1;
		else if (carbon_near_alcool(igr)==1 and carbon_near_alcool_tmp(igr+1)<carbon_near_alcool_tmp(igr))
		  carbon_near_alcool(igr+1)=1;
			        
		for (icycle=0;icycle<5;icycle++)
		  if (carbon_cycle(igr,icycle)>0)
		    {
		      int found=0;
		      for (igr2=0;igr2<ngr;igr2++)
			if (carbon_cycle(igr2,icycle)>0 and igr2!=igr)
			  if (found==0)
			    {
			      found=1;
			      if (carbon_near_alcool(igr2)==1 and carbon_near_alcool_tmp(igr)<carbon_near_alcool_tmp(igr2))
				carbon_near_alcool(igr)=1;
			      else if (carbon_near_alcool(igr)==1 and carbon_near_alcool_tmp(igr2)<carbon_near_alcool_tmp(igr))
				carbon_near_alcool(igr2)=1;
			    }
		    }
		  	       
		if (surrogate[i].smile.substr(igr+1,1)=="(")
		  {		    		    
		    int ipar=0;
		    for (igr2=igr+2;igr2<ngr;igr2++)
		      if (ipar==0 and surrogate[i].smile.substr(igr2,1)==")")			
			ipar=igr2;
		      else if (ipar<0 and surrogate[i].smile.substr(igr2,1)==")")			
			ipar=ipar+1;
		      else if (ipar<=0 and surrogate[i].smile.substr(igr2,1)=="(" and igr2>igr+2)
			{
			  //if (debug_mode==1) cout << "une ramification s'est ouverte" << endl;
			  ipar=ipar-1;
			}

		     if (carbon_near_alcool(igr+2)==1 and carbon_near_alcool_tmp(igr)<carbon_near_alcool_tmp(igr+2))
		      carbon_near_alcool(igr)=1;
		    else if (carbon_near_alcool(igr)==1 and carbon_near_alcool_tmp(igr+2)<carbon_near_alcool_tmp(igr))
		      carbon_near_alcool(igr+2)=1;

		     if (carbon_near_alcool(ipar+1)==1 and carbon_near_alcool_tmp(igr)<carbon_near_alcool_tmp(ipar+1))
		      carbon_near_alcool(igr)=1;
		    else if (carbon_near_alcool(igr)==1 and carbon_near_alcool_tmp(ipar+1)<carbon_near_alcool_tmp(igr))
		      carbon_near_alcool(ipar+1)=1;
		     
		    if (surrogate[i].smile.substr(ipar+1,1)=="(")
		      {
			if (carbon_near_alcool(ipar+2)==1 and carbon_near_alcool_tmp(igr)<carbon_near_alcool_tmp(ipar+2))
			  carbon_near_alcool(igr)=1;
			else if (carbon_near_alcool(igr)==1 and carbon_near_alcool_tmp(ipar+2)<carbon_near_alcool_tmp(igr))
			  carbon_near_alcool(ipar+2)=1;					        
										     
			int ipar2=0;
			for (igr2=ipar+1;igr2<ngr;igr2++)
			  if (ipar2==0 and surrogate[i].smile.substr(igr2,1)==")")			
			    ipar2=igr2;
			  else if (ipar2<0 and surrogate[i].smile.substr(igr2,1)==")")			
			    ipar2=ipar2+1;
			  else if (ipar2<=0 and surrogate[i].smile.substr(igr2,1)=="(" and igr2>ipar+1)
			    {
			      //if (debug_mode==1) cout << "une ramification s'est ouverte" << endl;
			      ipar2=ipar2-1;
			    }

			if (carbon_near_alcool(ipar2+1)==1 and carbon_near_alcool_tmp(igr)<carbon_near_alcool_tmp(ipar2+1))
			  carbon_near_alcool(igr)=1;
			else if (carbon_near_alcool(igr)==1 and carbon_near_alcool_tmp(ipar2+1)<carbon_near_alcool_tmp(igr))
			  carbon_near_alcool(ipar2+1)=1;					       			
		      }
		    
		  }
	    
	    
	      }
	  }
	  
	//if (debug_mode==1) cout << carbon_near_alcool << endl;


	//if (debug_mode==1) cout << smile2 << endl;
	//if (debug_mode==1) cout << carbon_near_alcool << endl;

	int excess_carbon=0;
	for (igr=0;igr<ngr;igr++)
	  excess_carbon+=max(carbon_taken(igr)-1,0);

	if (debug_mode==1) cout << " " << ngr << endl;
	if (debug_mode==1) cout << smile_save << " " << smile_save.length() << endl;
	if (debug_mode==1) cout << smile2 << " " << smile2.length() << endl;
	for (igr=0;igr<ngr;igr++)
	  {
	    int iket2,iket3;
	    iket3=-1;
	    for (iket2=0;iket2<10;iket2++)
	      if (carbon_ketone(igr,iket2)>0)
		iket3=iket2;
	    if (debug_mode==1 and iket3>=0) cout << smile_save.substr(igr,1) << " " << carbon_hydrogen(igr) << " et " << carbon_ketone(igr,iket3) << " " << carbon_taken(igr) << endl;
	  }

	if (excess_carbon>0)
	  {
	    if (outFlux == 1)
	      {
		cout << "WARNING: " << excess_carbon << " carbon to be removed" << endl;
		cout << "WARNING: exact structure not possible, use of an approximate one" << endl;
	      }
	    else if (outFlux == 2)
	      {
		fileFlux << "WARNING: " << excess_carbon << " carbon to be removed" << endl;
		fileFlux << "WARNING: exact structure not possible, use of an approximate one" << endl;
	      }

	    for (igr=0;igr<ngr;igr++)
	      if (carbon_hydrogen(igr)==2 and carbon_near_alcool(igr)==0 and carbon_alcool(igr)==0 and excess_carbon>0)
		{	      
		  excess_carbon--;
		  carbon_hydrogen(igr)=-1;
		}

	    for (igr=0;igr<ngr;igr++)
	      if (carbon_hydrogen(igr)==1 and carbon_near_alcool(igr)==0 and carbon_alcool(igr)==0 and excess_carbon>0)
		{	      
		  excess_carbon--;
		  carbon_hydrogen(igr)=-1;
		}

	    for (igr=0;igr<ngr;igr++)
	      if (carbon_hydrogen(igr)==3 and carbon_near_alcool(igr)==0 and carbon_alcool(igr)==0 and excess_carbon>0)
		{	      
		  excess_carbon--;
		  carbon_hydrogen(igr)=-1;
		}

	    for (igr=0;igr<ngr;igr++)
	      if (carbon_hydrogen(igr)==0 and carbon_near_alcool(igr)==0 and carbon_alcool(igr)==0 and excess_carbon>0)
		{	      
		  excess_carbon--;
		  carbon_hydrogen(igr)=-1;
		}
	    
	    for (igr=0;igr<ngr;igr++)
	      if (carbon_hydrogen(igr)==2 and carbon_alcool(igr)==0 and excess_carbon>0)
		{	      
		  excess_carbon--;
		  carbon_hydrogen(igr)=-1;
		}

	    for (igr=0;igr<ngr;igr++)
	      if (carbon_hydrogen(igr)==1 and carbon_alcool(igr)==0 and excess_carbon>0)
		{	      
		  excess_carbon--;
		  carbon_hydrogen(igr)=-1;
		}

	    for (igr=0;igr<ngr;igr++)
	      if (carbon_hydrogen(igr)==3 and carbon_alcool(igr)==0 and excess_carbon>0)
		{	      
		  excess_carbon--;
		  carbon_hydrogen(igr)=-1;
		}

	    for (igr=0;igr<ngr;igr++)
	      if (carbon_hydrogen(igr)==0 and carbon_alcool(igr)==0 and excess_carbon>0)
		{	      
		  excess_carbon--;
		  carbon_hydrogen(igr)=-1;
		}

	    for (igr=0;igr<ngr;igr++)
	      if (carbon_hydrogen(igr)==2 and excess_carbon>0)
		{	      
		  excess_carbon--;
		  carbon_hydrogen(igr)=-1;
		}

	    for (igr=0;igr<ngr;igr++)
	      if (carbon_hydrogen(igr)==1 and excess_carbon>0)
		{	      
		  excess_carbon--;
		  carbon_hydrogen(igr)=-1;
		}

	    for (igr=0;igr<ngr;igr++)
	      if (carbon_hydrogen(igr)==3 and excess_carbon>0)
		{	      
		  excess_carbon--;
		  carbon_hydrogen(igr)=-1;
		}

	    for (igr=0;igr<ngr;igr++)
	      if (carbon_hydrogen(igr)==0 and excess_carbon>0)
		{	      
		  excess_carbon--;
		  carbon_hydrogen(igr)=-1;
		}

	    if (excess_carbon>0)
	      {
		if (outFlux == 1)
		  cout << "WARNING: " << excess_carbon << " excess carbon cannot be removed" << endl;
		else if (outFlux == 2)
		  fileFlux << "WARNING: " << excess_carbon << " excess carbon cannot be removed" << endl;
	      }

	    
	  }

	if (debug_mode==1) cout << " " << endl;
	for (igr=0;igr<ngr;igr++)
	  {
	    int iket2,iket3;
	    iket3=-1;
	    for (iket2=0;iket2<10;iket2++)
	      if (carbon_ketone(igr,iket2)>0)
		iket3=iket2;
	    if (debug_mode==1 and iket3>=0) cout << smile_save.substr(igr,1) << " " << carbon_hydrogen(igr) << " et " << carbon_ketone(igr,iket3) << " " << carbon_taken(igr) << endl;
	  }
	
	for (igr=0;igr<ngr;igr++)
	  if (carbon_hydrogen(igr)>=0)
	    {
	      if (carbon_arom(sum_len_group)>0 and carbon_alcool(igr)>0)
		{
		  cout << "CRITICAL: alcohol on aromatic. Should be phenol "  << endl;
		  cout << "problem in smiles.cxx. Try manual group decomposition" << endl;
		  if (outFlux == 2)
		    {
		      cout << "CRITICAL: alcohol on aromatic. Should be phenol "  << endl;
		      cout << "problem in smiles.cxx. Try manual group decomposition" << endl;
		    }
		  exit(0);
		}
	      
	      if (carbon_alcool(igr)>0)
		{
		  if (carbon_hydrogen(igr)==3) surrogate[i].groups[4]+=1;
		  else if (carbon_hydrogen(igr)==2) surrogate[i].groups[5]+=1;
		  else if (carbon_hydrogen(igr)==1) surrogate[i].groups[6]+=1;
		  else if (carbon_hydrogen(igr)==0) surrogate[i].groups[7]+=1;
		}
	      else if (carbon_near_alcool(igr)>0)
		{
		  if (carbon_hydrogen(igr)==3) surrogate[i].groups[8]+=1;
		  else if (carbon_hydrogen(igr)==2) surrogate[i].groups[9]+=1;
		  else if (carbon_hydrogen(igr)==1) surrogate[i].groups[10]+=1;
		  else if (carbon_hydrogen(igr)==0) surrogate[i].groups[11]+=1;
		}
	      else if (sum(carbon_alcool)>0) //carbon_tail(igr)>0)
		{
		  if (carbon_hydrogen(igr)==3) surrogate[i].groups[12]+=1;
		  else if (carbon_hydrogen(igr)==2) surrogate[i].groups[13]+=1;
		  else if (carbon_hydrogen(igr)==1) surrogate[i].groups[14]+=1;
		  else if (carbon_hydrogen(igr)==0) surrogate[i].groups[15]+=1;
		}
	      else
		{
		  if (carbon_hydrogen(igr)==3) surrogate[i].groups[0]+=1;
		  else if (carbon_hydrogen(igr)==2) surrogate[i].groups[1]+=1;
		  else if (carbon_hydrogen(igr)==1) surrogate[i].groups[2]+=1;
		  else if (carbon_hydrogen(igr)==0) surrogate[i].groups[3]+=1;	      
		}
	    }

	if (debug_mode==1)
	  {
	    cout << "Taken carbons (end): " << endl;
	    for (igr=0;igr<ngr;igr++)
	      cout <<  surrogate[i].smile.substr(igr,1) << " " << carbon_taken(igr) << endl;
	    cout << " " << endl;
	  }
	
        int out_c=0;
	int out_o=0;
	for (j=0;j<60;j++)
	  if (surrogate[i].groups[j]>0)
	    {
	      if (outFlux == 1)
		cout <<j<<" "<<name_group(j) << " " << surrogate[i].groups[j] << " " << nc_group(j) << " " << no_group(j) << endl;
	      else if (outFlux == 2)
		fileFlux <<j<<" "<<name_group(j) << " " << surrogate[i].groups[j] << " " << nc_group(j) << " " << no_group(j) << endl;
	      out_c=out_c+int(surrogate[i].groups[j]*nc_group(j));
	      out_o=out_o+int(surrogate[i].groups[j]*no_group(j));
	    }

	if (out_c!=total_c)
	  {
	    if (outFlux == 1)
	      cout << "WARNING: not the good number of carbons " << out_c << " instead of " << total_c << endl;
	    else if (outFlux == 2)
	      fileFlux << "WARNING: not the good number of carbons " << out_c << " instead of " << total_c << endl;
	  }

	if (out_c-excess_carbon!=total_c)
	  {	   
	    cout << "CRITICAL: not the good number of carbons " << out_c-excess_carbon << " (without excess) instead of " << total_c << endl;
	    cout << "problem in smiles.cxx. Try manual group decomposition" << endl;
	    if (outFlux == 2)
	      {
		fileFlux << "CRITICAL: not the good number of carbons " << out_c-excess_carbon << " (without excess) instead of " << total_c << endl;
		cout << "problem in smiles.cxx. Try manual group decomposition" << endl;
	      }
	    exit(0);
	  }

	if (out_o!=total_o)
	  {	    
	    cout << "CRITICAL: not the good number of oxygens " << out_o << " instead of " << total_o << endl;
	    cout << "problem in smiles.cxx. Try manual group decomposition" << endl;
	    if (outFlux == 2)
	      {
		fileFlux << "CRITICAL: not the good number of oxygens " << out_o << " instead of " << total_o << endl;
		fileFlux << "problem in smiles.cxx. Try manual group decomposition" << endl;
	      }
	    exit(0);
	  }

        if (surrogate[i].is_generic)
	  if (surrogate[i].groups[37]>=2)
	    {
	      surrogate[i].aq_type="diacid";
	      surrogate[i].Kacidity1=3.95e-4;    // First acidity constant
	      surrogate[i].Kacidity2=7.70e-6;    // Second acidity constant
	    }
	  else if (surrogate[i].groups[37]>=1)
	    {
	      surrogate[i].aq_type="monoacid";	      
	      surrogate[i].Kacidity1=6.52e-4;    // First acidity constant
	    }
	
      }

  for (i=0;i<n;i++)
    if (surrogate[i].is_monomer and surrogate[i].is_organic)
      {
	int j=surrogate[i].ioligo;	
	for (int igr=0;igr<60;igr++)
	  surrogate[j].groups[igr]=surrogate[i].moligo*surrogate[i].groups[igr];
      }

  if (debug_mode==1 or config.SOAPlog==3)
    {
      cout << "SUCCESS" << endl;
      exit(0);
    }
}


void get_vectors(model_config &config, vector<species>& surrogate)
{
  int i,j;
  int n=surrogate.size();
  double tmp_group;
  // for output checking
  Array<string,1> name_group;
  name_group.resize(60);
  name_group(0)="CH3: ";
  name_group(1)="CH2: ";
  name_group(2)="CH: ";
  name_group(3)="C: ";
  name_group(4)="CH3 linked to alcohol: ";
  name_group(5)="CH2 linked to alcohol: ";
  name_group(6)="CH linked to alcohol: ";  
  name_group(7)="C linked to alcohol: ";
  name_group(8)="CH3 between several alcohol groups: ";
  name_group(9)="CH2 between several alcohol groups: ";
  name_group(10)="CH between several alcohol groups: ";
  name_group(11)="C between several alcohol groups: ";
  name_group(12)="CH3 in tails of alcohol molecule: ";
  name_group(13)="CH2 in tails of alcohol molecule: ";
  name_group(14)="CH in tails of alcohol molecule: ";
  name_group(15)="C in tails of alcohol molecule: ";
  name_group(16)="CH2=CH double bound: ";
  name_group(17)="CH=CH double bound: ";
  name_group(18)="CH2=C double bound: ";
  name_group(19)="CH=C double bound: ";
  name_group(20)="C=C double bound: ";
  name_group(21)="aromatic carbon with 1 hydrogen: ";
  name_group(22)="aromatic carbon with 0 hydrogen: ";
  name_group(23)="aromatic carbon linked to CH3: ";
  name_group(24)="aromatic carbon linked to CH2: ";
  name_group(25)="aromatic carbon linked to CH: ";
  name_group(26)="alcohol group: ";
  name_group(27)="water: ";
  name_group(28)="phenol: ";
  name_group(29)="CH3CO ketone: ";
  name_group(30)="CH2CO ketone: ";
  name_group(31)="aldehyde: ";
  name_group(32)="CH3COO ester: "; 
  name_group(33)="CH2COO ester: ";
  name_group(34)="CH3O ether: ";
  name_group(35)="CH2O ether: ";
  name_group(36)="CHO ether: ";
  name_group(37)="acid: ";
  name_group(38)="Nitro aromatic: ";
  name_group(39)="CH2ONO2 nitrate: ";
  name_group(40)="CHONO2 nitrate: ";
  name_group(41)="CONO2 nitrate: ";
  name_group(42)="CH2OOH hydroxyperoxide: ";
  name_group(43)="CHOOH hydroxyperoxide: ";
  name_group(44)="COOH hydroxyperoxide: ";
  name_group(45)="CH3OOCH2 peroxide: ";
  name_group(46)="CH3OOCH peroxide: ";
  name_group(47)="CH3OOC peroxide: ";
  name_group(48)="CH2OOCH2 peroxide: ";
  name_group(49)="CH2OOCH peroxide: ";
  name_group(50)="CH2OOC peroxide: ";
  name_group(51)="CHOOCH peroxide: ";
  name_group(52)="CHOOC peroxide: ";
  name_group(53)="COOC peroxide: ";
  name_group(54)="PAN: ";
  name_group(55)="Peroxyacetyl acid: ";
  name_group(56)="O=COC=O group: ";
  name_group(57)="CH3NO2 group: ";
  name_group(58)="CH2NO2 group: ";
  name_group(59)="CHNO2 group: ";
  /*
  name_group={
        "CH3","CH2","CH1","CH0", // group C
        "CH3-OH","CH2-OH","CH1-OH","CH0-OH", //group C[OH]
        "OHCH3OH","OHCH2OH","OHCHOH","OHCOH", //group Calcohol between several alcohol groups
        "OHCH3","OHCH2","OHCH","OHC", //group Calcohol-tail in tails of alcohol molecule
        "CH2=CH","CH=CH","CH2=C","CH=C","C=C", //group C=C
        "AC-H","AC-0H", //group aromatic carbon (AC)
        "AC-CH3","AC-CH2","AC-CH", // group AC-C
        "OH",  //group OH [CX4,CX3;!$([CX3]=[//7,//8])][OX2;H1] CH3OH
        "H2O", // group H2O
        "ACOH", // group ACOH
        "CH3CO (ketone)","CH2CO (ketone)", //group ketone
        "CHO",   //group aldehyde  
        "CH3COO (ester)","CH2COO (ester)", //group ester
        "CH3O (ether)","CH2O (ether)","CHO (ether)", //group ether  
        "COOH",  //group acid
        "ACNO2",   //group ACNO2
        "CH2ONO2","CHONO2","CONO2", //group NO3 //
        "CH2OOH","CHOOH","COOH", //group CO-OH hydroxyperoxide
        "CH3OOCH2","CH3OOCH","CH3OOC","CH2OOCH2","CH2OOCH","CH2OOC","CHOOCH","CHOOC","COOC", //group CO-OC peroxide //
        "PAN",  //group PAN
        "COOOH",//Peroxyacetyl acid
        "O=COC=O",
        "CH3NO2","CH2NO2","CHNO2"}
  */

  for (i=0;i<n;i++)
    if (surrogate[i].smile!="" and surrogate[i].smile[0]=='&')
      {
        //cout <<"===="<< i <<"===="<<endl;
        //cout << surrogate[i].name << " is constructed from vectors: " << surrogate[i].smile << endl;
        //cout << surrogate[i].smile.substr(2, 5) << endl;
        for (j=0; j<60; j++)
          {
            surrogate[i].groups[j]=0.;
            tmp_group = atof(surrogate[i].smile.substr(j*9+1,8).c_str());
            if (tmp_group > 0.)
              {
                //cout<<j<<" "<<name_group(j)<<" read from "<<surrogate[i].smile.substr(j*9+1,8)<<" as "<<tmp_group<<endl;
                surrogate[i].groups[j]=tmp_group;
              }
            else if (tmp_group < 0.)
              {
                cout<<"Warning!!! aerosol vector input < 0.!!! set to 0. Please check aerosol species list."<<endl;
                cout<<j<<" "<<name_group(j)<<" read from "<<surrogate[i].smile.substr(j*9+1,8)<<" as "<<endl;
                surrogate[i].groups[j]=0.;
              }
          }

        if (surrogate[i].is_generic)
	    if (surrogate[i].groups[37]>=2.)
	      {
	        surrogate[i].aq_type="diacid";
	        surrogate[i].Kacidity1=3.95e-4;    // First acidity constant
	        surrogate[i].Kacidity2=7.70e-6;    // Second acidity constant
	      }
	    else if (surrogate[i].groups[37]>=1.)
	      {
	        surrogate[i].aq_type="monoacid";	      
	        surrogate[i].Kacidity1=6.52e-4;    // First acidity constant
	      }

        // same as in smiles function
	/*
        if (surrogate[i].is_monomer and surrogate[i].is_organic)
	    {
	      int j=surrogate[i].ioligo;	
	      for (int igr=0;igr<60;igr++)
	        surrogate[j].groups[igr]=surrogate[i].moligo*surrogate[i].groups[igr];
	    }*/
      }
    else if (surrogate[i].smile!="" and surrogate[i].smile[0]=='$')
      {
	int ngr=surrogate[i].smile.length();
	int igr;
	string smile2=surrogate[i].smile.substr(1,ngr-1);
	string sgroup="";
	string sigroup="";
        for (j=0; j<60; j++)
	  surrogate[i].groups[j]=0.;

	j=0;
	for (igr=0;igr<ngr;igr++)
	  {
	    if (igr==ngr-1 or smile2.substr(igr,1)=="|")
	      {
		//cout << sigroup << " " << sgroup << endl;
		surrogate[i].groups[atoi(sigroup.c_str())-1]=atof(sgroup.c_str());
		sgroup="";
		sigroup="";
		j=0;
	      }
	    else if (smile2.substr(igr,1)=="=")
	      j=1;
	    else
	      if (j==0)
		sigroup+=smile2.substr(igr,1);
	      else
		sgroup+=smile2.substr(igr,1);
	  }

	if (config.SOAPlog == 1)
	  {
	    cout <<"===="<< i <<"===="<<endl;
	    cout << surrogate[i].name << " is constructed from vectors: " << surrogate[i].smile << endl;
	    for (j=0; j<60; j++)
	      if (surrogate[i].groups[j]>0.)
		cout<<j<<" "<<name_group(j)<<" equal to "<< surrogate[i].groups[j] << endl;
	  }
      }

  /*
  for (i=0;i<n;i++)
    {
      cout << surrogate[i].name;
      for (j=0; j<60; j++)
	cout<< ";"<< surrogate[i].groups[j];
      cout << endl;
    }*/
}

