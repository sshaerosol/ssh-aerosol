//!!-----------------------------------------------------------------------
//!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
//!!     SSH-aerosol is distributed under the GNU General Public License v3
//!!-----------------------------------------------------------------------

#include <fstream>
using namespace ssh_soap;

void get_smiles(model_config &config, vector<species>& surrogate)
{
  int i,j;
  int n=surrogate.size();
  
  // output flux: 1=screen, 2=file
  int outFlux;
  outFlux=2;
  
  //if (outFlux==2)
  //  {
      //Déclaration d'un flux permettant d'écrire dans un fichier.
  string const nomFichier("smile2UNIFAC.decomp");
  ofstream fileFlux(nomFichier.c_str());
  //  }	    
  
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
    //if (surrogate[i].smile!="C=C")
    //  {
    //    cout <<"===="<< i <<"===="<<endl;
    //	cout << surrogate[i].name << " is constructed from smiles: " << surrogate[i].smile << endl;
//	cout << "WARNING: CH2=CH2 not available" << endl;
//	cout << "Group CH2=CH used instead" << endl;
//	surrogate[i].groups[16]+=1	
//      }
    if (surrogate[i].smile!="" and surrogate[i].smile[0]!='&')
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
	//cout << surrogate[i].smile.substr(2, 5) << endl;
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
	Array<int, 1> linked_to_arom,carbon_near_alcool_tmp,carbon_near_alcool_save,carbon_near_alcool_bef;
	Array<int, 2> carbon_cycle,carbon_ether,carbon_ketone,carbon_peroxide;	
	int igr,igr2;
	int ngr=surrogate[i].smile.length();
        int total_c=0;
        int total_o=0;
        for (igr=0;igr<ngr;igr++)
           {
	     if (surrogate[i].smile.substr(igr,1)=="O")
	       total_o++;
	     if (surrogate[i].smile.substr(igr,1)=="c")
	       total_c++;
	     if (surrogate[i].smile.substr(igr,1)=="C")
	       total_c++;
	   }

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
	carbon_ketone.resize(ngr+3,ngr_max);
	carbon_peroxide.resize(ngr+3,ngr_max);
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
	linked_to_arom=0;
	carbon_near_alcool_bef=0;
	carbon_peroxide=-1;
	int iether=0;
	int iket=0;
	int ipero=0;
	
	string smile2="";   
	int last_pos=0;      
	igr=0;       
	while (igr<ngr)
	  {
	    //cout << "ipos: " << igr << " " << surrogate[i].smile.substr(igr,1) << endl;
	    if (surrogate[i].smile.substr(igr,4)=="O=C1")
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
		//cout << "found" << endl;
		//carbon_alcool(1)=1;
		smile2+="C(O)";
		last_pos=igr+3;
		carbon_arom(0)=1;
		ipos=ipos+3;
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
	    
	    if (surrogate[i].smile.substr(igr,19)=="C(=O)OO[N+](=O)[O-]")
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
	    else if (surrogate[i].smile.substr(igr,7)=="OC(=O)C" and igr!=0)
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"F";	      
		last_pos=igr+6;		
	      }
	    else if (surrogate[i].smile.substr(igr,15)=="OO[N+](=O)[O-])")
	      {
		if (outFlux==1)
	          {
	            cout << "WARNING: group OO[N+](=O)[O-]) not available" << endl;
		    cout << "Group O[N+](=O)[O-] used instead" << endl;
	          }
	        else if (outFlux==2)
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
		if (outFlux==1)
	          {
	            cout << "WARNING: group OON(=O)=O not available" << endl;
		    cout << "Group ON(=O)=O used instead" << endl;
	          }
	        else if (outFlux==2)
	          {
	            fileFlux << "WARNING: group OO[N+](=O)[O-]) not available" << endl;
	            fileFlux << "Group ON(=O)=O used instead" << endl;
	          }
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"ON(=O)=O";	
		last_pos=igr+9;
		ipos=ipos+8;
		total_o--;
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
	    else if (surrogate[i].smile.substr(igr,2)=="OC" and igr==0)
	      {
		//carbon_alcool(1)=1;
		smile2+="C(O)";
		last_pos=igr+2;
		ipos=ipos+3;
	      }
	    else if (surrogate[i].smile.substr(igr,7)=="C(=O)OC")
	      {
		smile2+=surrogate[i].smile.substr(last_pos,igr-last_pos)+"E";	      
		last_pos=igr+6;
	      }
	    else if (surrogate[i].smile.substr(igr,7)=="C(=O)OO")
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
	    else if (surrogate[i].smile.substr(igr,6)=="C(=O)O" and surrogate[i].smile.substr(igr,7)!="C(=O)OC")
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


	//cout << surrogate[i].smile << endl;
	//cout << "apres " << smile2 << endl;
	//cout << carbon_cycle << endl;
	
	ngr=smile2.length();	
	//if (smile2.substr(0,2)=="OC
	//cout << smile2 << endl;
	/*
	if (smile2.substr(0,6)=="OC(=O)")
	  {
	    //cout << "found" << endl;
	    //carbon_alcool(1)=1;
	    smile2="A"+smile2.substr(6,ngr-6);
	    ngr=smile2.length();	      
	  }
	else if (smile2.substr(0,2)=="OC")
	  {
	    //cout << "found" << endl;
	    //carbon_alcool(1)=1;
	    smile2="C(O)"+smile2.substr(2,ngr-2);
	    ngr=smile2.length();	      
	  }*/
	/*
	  for (igr=0;igr<ngr-4;igr++)
	  if (smile2.substr(igr,3)=="(O)")
	  {
	  cout << "found" << endl;
	  carbon_alcool(igr)+=1;
	  }*/
	//cout << carbon_alcool << endl;
	if (ngr>=2)
	  if (smile2.substr(ngr-2,2)=="CO" and carbon_arom(ngr-2)==0)
	    {
	      carbon_alcool(ngr-2)=1;	      
	    }
	 
		
	/*
	  if (sum(carbon_alcool)>0)	  
	  carbon_tail=1;			*/

	int will_be_hydroperoxide=0;
	int will_be_nitrate=0;
	int will_be_nitrite=0;
	int sum_len_group=0;
	surrogate[i].smile=smile2;
	
	if (outFlux==1)
	  {
	    cout << smile2 << endl;
	  }
	else if (outFlux==2)
	  {
	    fileFlux << smile2 << endl;
	  }

	//cout << carbon_arom << endl;
	while (surrogate[i].smile!="" and surrogate[i].smile[0]!='&')
	  {	  
	    int len_group=0;
	    //int ngr=surrogate[i].smile.length();
	    //cout << surrogate[i].smile << endl;
	    if (surrogate[i].smile.substr(0,1)=="A")
	      {		
		//cout << "is acid" << endl;
		surrogate[i].groups[37]+=1;
		len_group=1;
	      }
	    else if (surrogate[i].smile.substr(0,1)=="B")
	      {		
		//cout << "is acid" << endl;
		surrogate[i].groups[55]+=1;
		len_group=1;
	      }	 
	    else if (surrogate[i].smile.substr(0,1)=="P")
	      {		
		//cout << "is acid" << endl;
		surrogate[i].groups[54]+=1;
		len_group=1;
	      }
	    else if (surrogate[i].smile.substr(0,1)=="L")
	      {
		int icycle;
		len_group=1;
		//cout << carbon_cycle << endl;
		
		carbon_ketone(sum_len_group+1,iket)=1;
		for (icycle=0;icycle<5;icycle++)
		  if (carbon_cycle(sum_len_group,icycle)>0)
		    for (igr=0;igr<ngr;igr++)
		      if (igr!=sum_len_group and carbon_cycle(igr,icycle)>0)			 
			carbon_ketone(igr,iket)=1;			      
		iket++;
	      }	    
	    else if (surrogate[i].smile.substr(0,15)=="C([N+](=O)[O-])" and carbon_arom(sum_len_group)==1)
	      {		
		//cout << "is nitro aromatique" << endl;
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
		//cout << "is aldehyde" << endl;
		surrogate[i].groups[31]+=1;
		len_group=3;
		if (surrogate[i].smile.substr(0,4)=="O=CK")
		  {
		    carbon_ketone(sum_len_group+2,iket)=1;
		    carbon_ketone(sum_len_group+4,iket)=1;
		    carbon_taken(sum_len_group+2)=1;
		    iket++;
		  }
		if (surrogate[i].smile.substr(0,3)=="C=O")
		  carbon_taken(sum_len_group)++;
	      }
	    else if (surrogate[i].smile.substr(0,3)=="OOC" and sum_len_group==0)	      
	      {
		//cout << "is hydroxyperoxide" << endl;
		len_group=2;
		will_be_hydroperoxide=1;
		/*
		if (surrogate[i].smile.substr(3,1)=="" or previous=="" or surrogate[i].smile.substr(0,4)=="OOC)")
		  {
		    //cout << "is terminal" <<endl;
		    surrogate[i].groups[42]+=1;
		  }
		else if  (surrogate[i].smile.substr(0,4)=="OOCC")
		  {
		    //cout << "has one oxygen" <<endl;
		    surrogate[i].groups[43]+=1;
		  }
		else
		  {
		    //cout << "has no oxygen" <<endl;
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
	    else if (surrogate[i].smile.substr(0,9)=="O=N(=O)C" and sum_len_group==0)	      
	      {
		len_group=7;
		will_be_nitrite=1;      		
	      }

	    /*
	      else if (surrogate[i].smile.substr(0,5)=="C(OO)")	      
	      {
	      cout << "is hydroperoxide" << endl;
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
	      cout << "hydroperoxyde with tho hydrogens" <<endl;
	      surrogate[i].groups[42]+=1;
	      }
	      else if  (nC==1)
	      {
	      cout << "hydroperoxyde with one hydrogen" <<endl;
	      surrogate[i].groups[43]+=1;
	      }
	      else
	      {
	      cout << "hydroperoxyde with no hydrogen" <<endl;
	      surrogate[i].groups[44]+=1;
	      }	      		
	      }*/
	    else if (surrogate[i].smile.substr(0,1)=="C" or surrogate[i].smile.substr(0,1)=="c")
	      {
		int hydroxyperoxide=0;
		int nitrate=0;
		int nitrite=0;
		int ester=0;
		int arom=0;		
		int icycle;
		if (will_be_ester==1)
		  {
		    will_be_ester=0;
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

		if (carbon_arom(sum_len_group)==1)
		  arom=1;
		
		//cout << "is carbon" << endl;
		int nC;
		if (arom==0)
		  {
		    nC=4;	      
		    nC=nC-(carbon_cycle(sum_len_group,0)+carbon_cycle(sum_len_group,1)+carbon_cycle(sum_len_group,2)+carbon_cycle(sum_len_group,3)+carbon_cycle(sum_len_group,4));				    
		    if (previous!="")
		      nC=nC-1;		    
		  }
		else
		  {
		    nC=1;
		    if (carbon_arom(sum_len_group+1)==1)
		      nC=nC+1;
		    if (surrogate[i].smile.substr(0,3)=="C(C" and carbon_arom(sum_len_group+2)==1)
		      nC=nC+1;
		    if (linked_to_arom(sum_len_group)==1)
		      {
			//cout << "found" << endl;
			nC=0;
		      }
		      
		    //nC=nC-(carbon_cycle(sum_len_group,0)+carbon_cycle(sum_len_group,1)+carbon_cycle(sum_len_group,2)+carbon_cycle(sum_len_group,3)+carbon_cycle(sum_len_group,4));
		    /*
		      if (previous!="")
		      if (carbon_arom(sum_len_group-1)<1)
		      nC=nC-1;*/
		  }
		int ipar;
		int ipar2;
		int igr;
		int is_double=0;
		//if (surrogate[i].smile.substr(0,1)=="=")
		//  is_double=1;

		//cout << "nC is " << nC << endl;
		if (surrogate[i].smile.substr(1,1)!="(" and surrogate[i].smile.substr(1,1)!=")")		  
		  if (arom==1 and carbon_arom(sum_len_group+1)==0)
		    {
		      linked_to_arom(sum_len_group)=1;
		      linked_to_arom(sum_len_group+1)=1;
		    }
		  else if (arom==0 and carbon_arom(sum_len_group+1)==1)
		    {
		      linked_to_arom(sum_len_group)=1;
		      linked_to_arom(sum_len_group+1)=1;
		    }

		if (surrogate[i].smile.substr(1,2)=="O)")
		  {		      
		    carbon_alcool(sum_len_group)=1;		 
		    iether++;
		  }
		
		if (surrogate[i].smile.substr(1,2)=="OC")
		  {		      
		    carbon_ether(sum_len_group,iether)=1;
		    carbon_ether(sum_len_group+2,iether)=1;
		    iether++;
		  }
		if (surrogate[i].smile.substr(1,1)=="X")
		  {		      
		    carbon_ether(sum_len_group,iether)=1;
		    for (icycle=0;icycle<5;icycle++)
		      if (carbon_cycle(sum_len_group+1,icycle)>0)
			for (igr=0;igr<ngr;igr++)
			  if (igr!=sum_len_group+1 and carbon_cycle(igr,icycle)>0)			 
			    carbon_ether(igr,iether)=1;
		    iether++;
		  }
		if (surrogate[i].smile.substr(1,1)=="K")
		  {		      
		    carbon_ketone(sum_len_group,iket)=1;
		    carbon_ketone(sum_len_group+2,iket)=1;
		    iket++;
		  }
		if (surrogate[i].smile.substr(1,3)=="OOC")
		  {		      
		    carbon_peroxide(sum_len_group,ipero)=1;
		    carbon_peroxide(sum_len_group+3,ipero)=1;		    
		    ipero++;
		  }
		if (surrogate[i].smile.substr(1,2)=="OX")
		  {		      
		    carbon_peroxide(sum_len_group,ipero)=1;
		    for (icycle=0;icycle<5;icycle++)
		      if (carbon_cycle(sum_len_group+2,icycle)>0)
			for (igr=0;igr<ngr;igr++)
			  if (igr!=sum_len_group+2 and carbon_cycle(igr,icycle)>0)			 
			    carbon_peroxide(igr,ipero)=1;
		    ipero++;
		  }
		if (surrogate[i].smile.substr(1,1)=="F")		  		      			   
		  ester++;
		  
		if (surrogate[i].smile.substr(1,surrogate[i].smile.length()-1)=="OO")
		  hydroxyperoxide++;
		if (surrogate[i].smile.substr(1,8)=="ON(=O)=O" or surrogate[i].smile.substr(1,13)=="O[N+](=O)[O-]")
		  nitrate++;
		if (surrogate[i].smile.substr(1,7)=="N(=O)=O" or surrogate[i].smile.substr(1,12)=="[N+](=O)[O-]")
		  nitrite++;		  		
		//cout << "nitrate: " << nitrate << endl;
		
		if (surrogate[i].smile.substr(0,2)=="C(")
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
		    
		    //cout << "branche " << endl;
		    if (surrogate[i].smile.substr(1,3)=="(O)")
		      carbon_alcool(sum_len_group)++;
		    if (surrogate[i].smile.substr(1,4)=="(OO)")
		      hydroxyperoxide++;
		    if (surrogate[i].smile.substr(1,10)=="(ON(=O)=O)" or surrogate[i].smile.substr(1,14)=="(O[N+](=O)[O-]")
		      nitrate++;
		    if (surrogate[i].smile.substr(1,9)=="(N(=O)=O)" or surrogate[i].smile.substr(1,13)=="([N+](=O)[O-]")
		      nitrite++;
		    if (surrogate[i].smile.substr(2,2)=="OC")
		      {		      
			carbon_ether(sum_len_group,iether)=1;
			carbon_ether(sum_len_group+3,iether)=1;
			iether++;
		      }
		    if (surrogate[i].smile.substr(2,1)=="X")
		      {		      
			carbon_ether(sum_len_group,iether)=1;
			for (icycle=0;icycle<5;icycle++)
			  if (carbon_cycle(sum_len_group+2,icycle)>0)
			    for (igr=0;igr<ngr;igr++)
			      if (igr!=sum_len_group+2 and carbon_cycle(igr,icycle)>0)			 
				carbon_ether(igr,iether)=1;
			iether++;
		      }
		    if (surrogate[i].smile.substr(2,1)=="K")
		      {		      
			carbon_ketone(sum_len_group,iket)=1;
			carbon_ketone(sum_len_group+3,iket)=1;
			iket++;
		      }
		    if (surrogate[i].smile.substr(2,3)=="OOC")
		      {		      
			carbon_peroxide(sum_len_group,ipero)=1;
			carbon_peroxide(sum_len_group+4,ipero)=1;		    
			ipero++;
		      }
		    if (surrogate[i].smile.substr(2,2)=="OX")
		      {		      
			carbon_peroxide(sum_len_group,ipero)=1;
			for (icycle=0;icycle<5;icycle++)
			  if (carbon_cycle(sum_len_group+3,icycle)>0)
			    for (igr=0;igr<ngr;igr++)
			      if (igr!=sum_len_group+3 and carbon_cycle(igr,icycle)>0)			 
				carbon_peroxide(igr,ipero)=1;
			ipero++;
		      }
		    if (surrogate[i].smile.substr(2,1)=="F")		  		      			   
		      ester++;
		    
		    
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
			    //cout << "une ramification s'est ouverte" << endl;
			    ipar=ipar-1;
			  }
		      }
		 
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
			carbon_alcool(sum_len_group)++;
		      }
		    if (surrogate[i].smile.substr(ipar+1,2)=="OC")
		      {		      
			carbon_ether(sum_len_group,iether)=1;
			carbon_ether(sum_len_group+ipar+2,iether)=1;
			iether++;
		      }
		    if (surrogate[i].smile.substr(ipar+1,1)=="X")
		      {		      
			carbon_ether(sum_len_group,iether)=1;
			for (icycle=0;icycle<5;icycle++)
			  if (carbon_cycle(sum_len_group+ipar+1,icycle)>0)
			    for (igr=0;igr<ngr;igr++)
			      if (igr!=sum_len_group+ipar+1 and carbon_cycle(igr,icycle)>0)			 
				carbon_ether(igr,iether)=1;
			iether++;
		      }
		    if (surrogate[i].smile.substr(ipar+1,1)=="K")
		      {
			//cout << "la " << endl;
			carbon_ketone(sum_len_group,iket)=1;
			carbon_ketone(sum_len_group+ipar+2,iket)=1;
			iket++;
		      }
		    if (surrogate[i].smile.substr(ipar+1,3)=="OOC")
		      {		      
			carbon_peroxide(sum_len_group,ipero)=1;
			carbon_peroxide(sum_len_group+ipar+3,ipero)=1;		    
			ipero++;
		      }
		    if (surrogate[i].smile.substr(ipar+1,2)=="OX")
		      {		      
			carbon_peroxide(sum_len_group,ipero)=1;
			for (icycle=0;icycle<5;icycle++)
			  if (carbon_cycle(sum_len_group+ipar+2,icycle)>0)
			    for (igr=0;igr<ngr;igr++)
			      if (igr!=sum_len_group+ipar+2 and carbon_cycle(igr,icycle)>0)			 
				carbon_peroxide(igr,ipero)=1;
			ipero++;
		      }
		    
		    if (surrogate[i].smile.substr(ipar+1,1)=="F")		  		      			   
		      ester++;
		    if (surrogate[i].smile.substr(ipar+1,3)=="(O)")
		      carbon_alcool(sum_len_group)++;
		    if (surrogate[i].smile.substr(ipar+1,surrogate[i].smile.length()-1-ipar)=="OO")
		      hydroxyperoxide++;
		    //cout << "la " << surrogate[i].smile.substr(ipar+1,surrogate[i].smile.length()-1-ipar) << endl;
		    if (surrogate[i].smile.substr(ipar+1,8)=="ON(=O)=O" or surrogate[i].smile.substr(ipar+1,13)=="O[N+](=O)[O-]")
		      nitrate++;
		    if (surrogate[i].smile.substr(ipar+1,7)=="N(=O)=O" or surrogate[i].smile.substr(ipar+1,12)=="[N+](=O)[O-]")
		      nitrite++;
		    //cout << ipar << endl;
		    //cout << "nCa is " << nC << endl;
		    if (surrogate[i].smile.substr(ipar+1,1)=="(")
		      {

			if (surrogate[i].smile.substr(ipar+1,4)=="(OO)")
			  hydroxyperoxide++;			
			if (surrogate[i].smile.substr(ipar+1,10)=="(ON(=O)=O)" or surrogate[i].smile.substr(ipar+1,14)=="(O[N+](=O)[O-]")
			  nitrate++;
			if (surrogate[i].smile.substr(ipar+1,9)=="(N(=O)=O)" or surrogate[i].smile.substr(ipar+1,13)=="([N+](=O)[O-]")
			  nitrite++;
			  
			if (surrogate[i].smile.substr(ipar+1,2)=="(F")		  		      			   
			  ester++;
			if (surrogate[i].smile.substr(ipar+2,2)=="OC")
			  {		      
			    carbon_ether(sum_len_group,iether)=1;
			    carbon_ether(sum_len_group+ipar+3,iether)=1;
			    iether++;
			  }
			if (surrogate[i].smile.substr(ipar+2,1)=="X")
			  {		      
			    carbon_ether(sum_len_group,iether)=1;
			    for (icycle=0;icycle<5;icycle++)
			      if (carbon_cycle(sum_len_group+ipar+2,icycle)>0)
				for (igr=0;igr<ngr;igr++)
				  if (igr!=sum_len_group+ipar+2 and carbon_cycle(igr,icycle)>0)			 
				    carbon_ether(igr,iether)=1;
			    iether++;
			  }
			if (surrogate[i].smile.substr(ipar+2,1)=="K")
			  {		      
			    carbon_ketone(sum_len_group,iket)=1;
			    carbon_ketone(sum_len_group+ipar+3,iket)=1;
			    iket++;
			  }
			if (surrogate[i].smile.substr(ipar+2,3)=="OOC")
			  {		      
			    carbon_peroxide(sum_len_group,ipero)=1;
			    carbon_peroxide(sum_len_group+ipar+4,ipero)=1;		    
			    ipero++;
			  }
			
			if (surrogate[i].smile.substr(ipar+2,2)=="OX")
			  {		      
			    carbon_peroxide(sum_len_group,ipero)=1;
			    for (icycle=0;icycle<5;icycle++)
			      if (carbon_cycle(sum_len_group+ipar+3,icycle)>0)
				for (igr=0;igr<ngr;igr++)
				  if (igr!=sum_len_group+ipar+3 and carbon_cycle(igr,icycle)>0)			 
				    carbon_peroxide(igr,ipero)=1;
			    ipero++;
			  }
			
			nC=nC-1;
			ipar2=0;
			for (igr=ipar+1;igr<int(surrogate[i].smile.length());igr++)
			  if (ipar2==0 and surrogate[i].smile.substr(igr,1)==")")			
			    ipar2=igr;
			  else if (ipar2<0 and surrogate[i].smile.substr(igr,1)==")")			
			    ipar2=ipar2+1;
			  else if (ipar2<=0 and surrogate[i].smile.substr(igr,1)=="(" and igr>ipar+1)
			    {
			      //cout << "une ramification s'est ouverte" << endl;
			      ipar2=ipar2-1;
			    }

			//cout << "I: " << ipar2 << " " << surrogate[i].smile.substr(ipar2+1,1) << " " << carbon_cycle(sum_len_group+ipar2+1,icycle) << endl;

			
			if (surrogate[i].smile.substr(ipar2+1,1)=="O" and ipar2+1==int(surrogate[i].smile.length())-1)
			  {		      			
			    carbon_alcool(sum_len_group)++;
			  }
			if (surrogate[i].smile.substr(ipar2+1,2)=="OC")
			  {		      
			    carbon_ether(sum_len_group,iether)=1;
			    carbon_ether(sum_len_group+ipar2+2,iether)=1;
			    iether++;
			  }
			if (surrogate[i].smile.substr(ipar2+1,1)=="X")
			  {
			    //cout << "iciii" << endl;
			   // cout << carbon_cycle << endl;
			    carbon_ether(sum_len_group,iether)=1;
			    for (icycle=0;icycle<5;icycle++)
			      if (carbon_cycle(sum_len_group+ipar2+1,icycle)>0)
				{
				  //cout << icycle << endl;
				  //cout << carbon_cycle << endl;
				  for (igr=0;igr<ngr;igr++)
				    if (igr!=sum_len_group+ipar2+1 and carbon_cycle(igr,icycle)>0)			 
				      carbon_ether(igr,iether)=1;
				}
			    iether++;
			  }
			if (surrogate[i].smile.substr(ipar2+1,1)=="K")
			  {		      
			    carbon_ketone(sum_len_group,iket)=1;
			    carbon_ketone(sum_len_group+ipar2+2,iket)=1;
			    iket++;
			  }
			if (surrogate[i].smile.substr(ipar2+1,3)=="OOC")
			  {		      
			    carbon_peroxide(sum_len_group,ipero)=1;
			    carbon_peroxide(sum_len_group+ipar2+2,ipero)=1;		    
			    ipero++;
			  }
			if (surrogate[i].smile.substr(ipar2+1,2)=="OX")
			  {		      
			    carbon_peroxide(sum_len_group,ipero)=1;
			    for (icycle=0;icycle<5;icycle++)
			      if (carbon_cycle(sum_len_group+ipar2+2,icycle)>0)
				for (igr=0;igr<ngr;igr++)
				  if (igr!=sum_len_group+ipar2+2 and carbon_cycle(igr,icycle)>0)			 
				    carbon_peroxide(igr,ipero)=1;
			    ipero++;
			  }
			if (surrogate[i].smile.substr(ipar2+1,surrogate[i].smile.length()-1-ipar2)=="OO")
			  hydroxyperoxide++;
			if (surrogate[i].smile.substr(ipar2+1,8)=="ON(=O)=O" or surrogate[i].smile.substr(ipar2+1,13)=="O[N+](=O)[O-]")
			  nitrate++;
			if (surrogate[i].smile.substr(ipar2+1,7)=="N(=O)=O" or surrogate[i].smile.substr(ipar2+1,13)=="[N+](=O)[O-]")
			  nitrite++;
			if (surrogate[i].smile.substr(ipar2+1,1)=="F")		  		      			   
			  ester++;
			if (surrogate[i].smile.substr(ipar2+1,3)=="(O)")
			  carbon_alcool(sum_len_group)++;
		    

			if (surrogate[i].smile.substr(ipar2+1,1)=="(")
			  {
		            if (outFlux==1)
	                      {
	                        cout << "encore un groupe" << endl;
	                      }
	                    else if (outFlux==2)
	                      {
	                        fileFlux << "encore un groupe" << endl;
	                      }
			    break;
			    nC=nC-1;
			  }

			if (surrogate[i].smile.substr(ipar2+1,1)=="=")
			  {
			    //cout << "la" << endl;
			    is_double=1;
			  }
			if (surrogate[i].smile.substr(ipar2+1,1)!="=" and surrogate[i].smile.substr(ipar2+1,1)!="" and surrogate[i].smile.substr(ipar2+1,1)!=")")
			  {
			    //cout << "B" << endl;
			    nC=nC-1;
			  }
		      
		      }
		    else
		      {
			if (surrogate[i].smile.substr(ipar+1,1)=="=")
			  {
			    //cout << "la2" << endl;
			    is_double=1;
			  }
			if (surrogate[i].smile.substr(ipar+1,1)!="=" and surrogate[i].smile.substr(ipar+1,1)!="" and surrogate[i].smile.substr(ipar+1,1)!=")")
			  {
			    //cout << "A" << endl;
			    nC=nC-1;
			  }
		      }
		       
		  }
		else
		  {
		    if (surrogate[i].smile.substr(1,1)=="=")
		      {
			//cout << "la" << endl;
			is_double=1;
		      }
		    if (surrogate[i].smile.substr(1,1)!="=" and surrogate[i].smile.substr(1,1)!="" and surrogate[i].smile.substr(1,1)!=")")
		      {
			//cout << "C " << endl;
			nC=nC-1;
		      }

		  }
		//cout << "nCb is "<< nC << endl;		
		carbon_nh(sum_len_group)=nC;
		if (ester>0)
		  {
		    if (nC==3)
		      {
			//cout << "ester with three hydrogens" <<endl;
			surrogate[i].groups[32]+=ester;
		      }
		    else 
		      {			
			//cout << "ester with two or less hydrogens" <<endl;
			surrogate[i].groups[32]+=ester;
		      }
		    carbon_taken(sum_len_group)+=ester;
		  }		
		if (hydroxyperoxide>0)
		  {
		    if (nC>=2)
		      {
			//cout << "hydroperoxyde with tho hydrogens" <<endl;
			surrogate[i].groups[42]+=hydroxyperoxide;
		      }
		    else if  (nC==1)
		      {
			//cout << "hydroperoxyde with one hydrogen" <<endl;
			surrogate[i].groups[43]+=hydroxyperoxide;
		      }
		    else
		      {
			//cout << "hydroperoxyde with no hydrogen" <<endl;
			surrogate[i].groups[44]+=hydroxyperoxide;
		      }
		    
		    if (carbon_taken(sum_len_group)>0 or hydroxyperoxide>1)
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
		    
		    carbon_taken(sum_len_group)+=hydroxyperoxide;
		  }
		if (nitrate>0)
		  {
		    if (nC>=2)
		      {
			//cout << "nitrate with tho hydrogens" <<endl;
			surrogate[i].groups[39]+=nitrate;
		      }
		    else if  (nC==1)
		      {
			//cout << "nitrate with one hydrogen" <<endl;
			surrogate[i].groups[40]+=nitrate;
		      }
		    else
		      {
			//cout << "nitrate with no hydrogen" <<endl;
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
		    if (nC==3)
		      {
		        if (outFlux==1)
	                  {
			    cout << "nitrite with three hydrogens" <<endl;
			  }
	                else if (outFlux==2)
	                  {
	                    fileFlux << "nitrite with three hydrogens" << endl;
	                  }
			surrogate[i].groups[57]+=nitrite;
		      }
		    else if  (nC==2)
		      {
		        if (outFlux==1)
	                  {
			    cout << "nitrite with two hydrogen" <<endl;
			  }
	                else if (outFlux==2)
	                  {
	                    fileFlux << "nitrite with two hydrogen" << endl;
	                  }
			surrogate[i].groups[58]+=nitrite;
		      }
		    else
		      {
		        if (outFlux==1)
	                  {
			    cout << "nitrite with one or zero hydrogen" <<endl;
			  }
	                else if (outFlux==2)
	                  {
	                    fileFlux << "nitrite with one or zero hydrogen" << endl;
	                  }
			surrogate[i].groups[59]+=nitrite;
		      }

		    if (carbon_taken(sum_len_group)>0 or nitrite>1)
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
		    
		    carbon_taken(sum_len_group)+=nitrite;
		    
		  }


		if (is_double==1 and nc_double_save>=0)
		  {
		    if (outFlux==1)
	              {
		        cout << "another double bound is currently treated" << endl;
		      }
	            else if (outFlux==2)
	              {
	                fileFlux << "another double bound is currently treated" << endl;
	              }

		    //break;
		  }
		else if (is_double==1)
		  {		    
		    nc_double_save=nC-2;
		    //cout << "beg double with " << nC-2 << endl;
		    carbon_taken(sum_len_group)++;
		  }
		else if (will_be_double==1)
		  {
		    nC=nC-1;
		    //cout << "double bounds between " << nc_double_save <<  "hydrogen carbon and " << nC << " hydrogen carbon" << endl;

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
	                    fileFlux << "WARNING: CH2=CH2 not available" << endl;
	                    fileFlux << "Group CH2=CH used instead" << endl;
	                  }
                      }
		    ///////////////////  
		    else
		      {
		        if (outFlux==1)
	                  {
	                    cout << "error double bounds type " << nC1 << " " << nC2 << endl;
                          }
	                else if (outFlux==2)
	                  {
                            fileFlux << "error double bounds type " << nC1 << " " << nC2 << endl;
	                  }
			break;
		      }
		    will_be_double=0;
		    nc_double_save=-1;
		    previous="C=C";
		    carbon_taken(sum_len_group)++;
		  }
		else if (arom==1)
		  {
		    //cout << "is aromatic carbon with " << nC << " hydrogens" << endl;
		    previous="c";
		  }
		else if (hydroxyperoxide==0 and nitrate==0 and ester==0 and nitrite==0)
		  {
		    //cout << "is carbon with " << nC << " hydrogens" << endl;
		    previous="C";
		    carbon_hydrogen(sum_len_group)=nC;
		  }		
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
	    else if (surrogate[i].smile.substr(0,2)=="OO")
	      {
		len_group=2;
		previous="OO";
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
	    else if (surrogate[i].smile.substr(0,4)=="KKKK")
	      {
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
		iket++;				
		
		len_group=4;
		previous="C";
	      }
	    else if (surrogate[i].smile.substr(0,3)=="KKK")
	      {
		carbon_ketone(sum_len_group+1,iket)=1;
		carbon_ketone(sum_len_group+2,iket)=1;
		carbon_taken(sum_len_group+1)++;
		iket++;

		carbon_ketone(sum_len_group+2,iket)=1;
		carbon_ketone(sum_len_group+3,iket)=1;
		carbon_taken(sum_len_group+2)++;
		iket++;
		
		len_group=3;
		previous="C";
	      }
	    else if (surrogate[i].smile.substr(0,2)=="KK")
	      {
		carbon_ketone(sum_len_group+1,iket)=1;
		carbon_ketone(sum_len_group+2,iket)=1;
		carbon_taken(sum_len_group+1)++;
		iket++;
		len_group=2;
		previous="C";
	      }
	    else if (surrogate[i].smile.substr(0,2)=="FK")
	      {
		carbon_ketone(sum_len_group+1,iket)=1;
		carbon_ketone(sum_len_group+2,iket)=1;
		carbon_taken(sum_len_group+1)=1;
		iket++;
		len_group=2;
		previous="C";
	      }
	    else if (surrogate[i].smile.substr(0,1)=="K")
	      {
		len_group=1;
		previous="C";
	      }
	    else if (surrogate[i].smile.substr(0,1)=="E")
	      {
		will_be_ester=1;
		len_group=1;
		previous="C";
	      }
	    else if (surrogate[i].smile.substr(0,1)=="F")
	      {		
		len_group=1;
		previous="C";
	      }
	    else if (surrogate[i].smile.substr(0,2)=="OX")
	      {		
		len_group=2;
		previous="C";
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
	    

	    sum_len_group+=len_group;
	    previous=surrogate[i].smile.substr(0,len_group);
	    surrogate[i].smile=surrogate[i].smile.substr(len_group,ngr-len_group);	    
	    if (len_group==0)
	      {
	        if (outFlux==1)
	          {
	            cout << surrogate[i].smile.substr(0,5) << endl;
		    cout << "group not found" << endl;
                  }
	        else if (outFlux==2)
	          {
                    fileFlux << surrogate[i].smile.substr(0,5) << endl;
		    fileFlux << "group not found" << endl;
	          }
		break;
	      }
	  }

	surrogate[i].smile=smile2;
	ngr=surrogate[i].smile.length();
	int igr3;
	//cout << ngr << " ether " << endl;
	//cout << carbon_ether << endl;	

	//cout << "molecules has " << sum(carbon_alcool) << " alcool group" << endl;
	surrogate[i].groups[26]+=sum(carbon_alcool);

	//cout << "ket " << endl;
	//cout << carbon_ketone << endl;
	//cout << carbon_taken << endl;       
	//cout << carbon_ketone << endl;
	//Treatment of peroxide molecules
	for (igr2=0;igr2<ngr_max;igr2++)
	  {
	    int i1=-1;
	    int i2=-1;	    
	    for (igr=0;igr<ngr;igr++)
	      if (carbon_peroxide(igr,igr2)>=0 and i1<0)
		i1=igr;
	      else if (carbon_peroxide(igr,igr2)>=0)
		i2=igr;

	    //cout << i1 << " " << i2 << endl;
	    if (i1>=0 and i2>=0)
	      {
		//cout << "ketone " << i1 << " " << i2 << endl;
		//cout << "taken " << carbon_taken(i1) << " " << carbon_taken(i2) << endl;
		if (carbon_taken(i1)>0 or carbon_taken(i2)>0)
		  if (outFlux==1)
	            {
	              cout << "WARNING: already occupied carbon" << endl;	
                    }
	          else if (outFlux==2)
	            {
                      fileFlux << "WARNING: already occupied carbon" << endl;	
	            }

		carbon_taken(i1)++;	      
		carbon_taken(i2)++;

		int nc1=max(carbon_nh(i1),carbon_nh(i2));
		int nc2=min(carbon_nh(i1),carbon_nh(i2));

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
	      }
	  }

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

	    //cout << i1 << " " << i2 << endl;
	    if (i1>=0 and i2>=0)
	      {
		//cout << "ketone " << i1 << " " << i2 << endl;
		//cout << "taken " << carbon_taken(i1) << " " << carbon_taken(i2) << endl;
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
	                  {
	                    cout << "Warning: subsequent ketone" << endl;
                          }
	                else if (outFlux==2)
	                  {
                            fileFlux << "Warning: subsequent ketone" << endl;	
	                  }

			if (carbon_nh(i1)>=carbon_nh(i2))
			  ketone1=1;
			else
			  ketone2=1;
		      }
		  }
		else
		  {
		    if (outFlux==1)
	              {
	                cout << "Warning: already occupied ketone" << endl;
                      }
	            else if (outFlux==2)
	              {
                        fileFlux << "Warning: already occupied ketone" << endl;	
	              }

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
		    //cout << "ketone with 3 hydrogens" << endl;
		    surrogate[i].groups[29]++;
		  }
		else 
		  {
		    //cout << "ketone with 2 or less hydrogens" << endl;
		    surrogate[i].groups[30]++;
		  }
	      }
	  }

	//cout << carbon_taken << endl;
     	
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

	    //cout << i1 << " " << i2 << endl;
	    if (i1>=0 and i2>=0)
	      {
		//cout << "ether " << i1 << " " << i2 << endl;
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
			//cout << "Warning: subsequent ether" << endl;
			if (carbon_nh(i1)>=carbon_nh(i2))
			  ether1=1;
			else
			  ether2=1;
		      }
		  }
		else
		  {
		    //cout << "Warning: already occupied ether" << endl;
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
		    //cout << "ether with 3 hydrogens" << endl;
		    surrogate[i].groups[34]++;
		  }
		else if (carbon_nh(iether)==2)
		  {
		    //cout << "ether with 2 hydrogens" << endl;
		    surrogate[i].groups[35]++;
		  }
		else
		  {	   
		    //cout << "ether with 0 or 1 hydrogen" << endl;
		    surrogate[i].groups[36]++;	      
		  }
	      }
	  }


	//cout << linked_to_arom << endl;
	//cout << carbon_hydrogen << endl;
	//cout << carbon_arom << endl;
	//Treatment of arom-CHx groups
	for (igr=0;igr<ngr;igr++)
	  if (carbon_hydrogen(igr)>=0 and linked_to_arom(igr)==1 and carbon_taken(igr)==0)	    
	    {
	      if (carbon_nh(igr)==3)
		{
		  //cout << "arom-CH3" << endl;
		  surrogate[i].groups[23]++;
		  surrogate[i].groups[22]--;
		}
	      else if (carbon_nh(igr)==2)
		{
		  //cout << "arom-CH2" << endl;
		  surrogate[i].groups[24]++;
		  surrogate[i].groups[22]--;
		}
	      else
		{	   
		  //cout << "arom-CH" << endl;
		  surrogate[i].groups[25]++;
		  surrogate[i].groups[22]--;
		}
	      carbon_hydrogen(igr)=-1;
	      carbon_taken(igr)=1;
	    }
	
	for (igr=0;igr<ngr;igr++)
	  if (carbon_arom(igr)==1)
	    {
	      //cout << "ici" << endl;
	      if (carbon_nh(igr)==1)
		{
		  //cout << "aromatic carbon with one hydrogen" << endl;		 
		  surrogate[i].groups[21]++;
		}
	      else 
		{
		  //cout << "aromatic carbon with no hydrogen" << endl;	     
		  surrogate[i].groups[22]++;
		}
	    }


	int icycle;
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
		    //cout << surrogate[i].smile.substr(igr,2) << endl;
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
			      cout << "the two branch has an alcool" << endl;
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
			      //cout << "une ramification s'est ouverte" << endl;
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
				  //cout << "une ramification s'est ouverte" << endl;
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
	 
	

	//cout << carbon_near_alcool_tmp << endl;

	int near_alcool_save=0;
	int alcool_save=0;
	for (igr=0;igr<ngr;igr++)
	  {
	    near_alcool_save=carbon_near_alcool_bef(igr);
	    //cout << igr << " " << near_alcool_save << " " << carbon_near_alcool_tmp(igr) << endl;
	    if (carbon_near_alcool_tmp(igr)>=0)
	      {
		//cout << surrogate[i].smile.substr(igr,2) << endl;
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
			  //cout << "une ramification s'est ouverte" << endl;
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
			      //cout << "une ramification s'est ouverte" << endl;
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

	//cout << "after near:" << endl;
	//cout << carbon_near_alcool << endl;
	//cout << "la :" << carbon_near_alcool_tmp << endl;
	//cout << carbon_near_alcool << endl;

	for (igr3=0;igr3<ngr;igr3++)
	for (igr=0;igr<ngr;igr++)
	  {
	    if (carbon_near_alcool_tmp(igr)>0)
	      {
		//cout << surrogate[i].smile.substr(igr,2) << endl;
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
			  //cout << "une ramification s'est ouverte" << endl;
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
			      //cout << "une ramification s'est ouverte" << endl;
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
	  
	//cout << carbon_near_alcool << endl;


	//cout << smile2 << endl;
	//cout << carbon_near_alcool << endl;

	int excess_carbon=0;
	for (igr=0;igr<ngr;igr++)
	  excess_carbon+=max(carbon_taken(igr)-1,0);

	if (excess_carbon>0)
	  {
	    if (outFlux==1)
	      {
	        cout << "WARNING: " << excess_carbon << " carbon to be removed" << endl;
	        cout << "WARNING: exact structure not possible, use of an approximate one" << endl;
              }
	    else if (outFlux==2)
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
	      if (outFlux==1)
	        {
	          cout << "WARNING: " << excess_carbon << " excess carbon cannot be removed" << endl;
                }
	      else if (outFlux==2)
	        {
                  fileFlux << "WARNING: " << excess_carbon << " excess carbon cannot be removed" << endl;
	        }


	    
	  }
	
	for (igr=0;igr<ngr;igr++)
	  if (carbon_hydrogen(igr)>=0)
	    {	      
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

	
        int out_c=0;
	int out_o=0;
	for (j=0;j<60;j++)
	  if (surrogate[i].groups[j]>0)
	    {
	      if (outFlux==1)
	        {
	          cout <<j<<" "<<name_group(j) << " " << surrogate[i].groups[j] << " " << nc_group(j) << " " << no_group(j) << endl;
                }
	      else if (outFlux==2)
	        {
                  fileFlux <<j<<" "<<name_group(j) << " " << surrogate[i].groups[j] << " " << nc_group(j) << " " << no_group(j) << endl;
	        }
	      out_c=out_c+int(surrogate[i].groups[j]*nc_group(j));
	      out_o=out_o+int(surrogate[i].groups[j]*no_group(j));
	    }

	if (out_c!=total_c)
	  if (outFlux==1)
	    {
	      cout << "WARNING: not the good number of carbons " << out_c << " instead of " << total_c << endl;
            }
	  else if (outFlux==2)
	    {
              fileFlux << "WARNING: not the good number of carbons " << out_c << " instead of " << total_c << endl;
	    }

	if (out_o!=total_o)
	  if (outFlux==1)
	    {
	      cout << "WARNING: not the good number of oxygens " << out_o << " instead of " << total_o << endl;
            }
	  else if (outFlux==2)
	    {
              fileFlux << "WARNING: not the good number of oxygens " << out_o << " instead of " << total_o << endl;
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
}


void get_vectors(model_config &config, vector<species>& surrogate)
{
  int i,j;
  
  // output flux: 1=screen, 2=file
  //int outFlux;
  //outFlux=2;
  
  //if (outFlux==2)
  //  {
      //Déclaration d'un flux permettant d'écrire dans un fichier.
  //string const nomFichier("smile2UNIFAC.vect");
  //ofstream fileFlux(nomFichier.c_str());
  //  }	
  
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
        cout <<"===="<< i <<"===="<<endl;
        cout << surrogate[i].name << " is constructed from vectors: " << surrogate[i].smile << endl;
 /*       
        if (outFlux==1)
	  {
	    cout <<"===="<< i <<"===="<<endl;
            cout << surrogate[i].name << " is constructed from vectors: " << surrogate[i].smile << endl;
          }
	else if (outFlux==2)
	  {
            fileFlux <<"===="<< i <<"===="<<endl;
            fileFlux << surrogate[i].name << " is constructed from vectors: " << surrogate[i].smile << endl;
	  } 
*/
        //cout << surrogate[i].smile.substr(2, 5) << endl;
        for (j=0; j<60; j++)
          {
            surrogate[i].groups[j]=0.;
            tmp_group = atof(surrogate[i].smile.substr(j*9+1,8).c_str());
            if (tmp_group > 0.)
              {
                cout <<j<<" "<<name_group(j)<<" read from "<<surrogate[i].smile.substr(j*9+1,8)<<" as "<<tmp_group<<endl;
/*                if (outFlux==1)
	          {
	            cout <<j<<" "<<name_group(j)<<" read from "<<surrogate[i].smile.substr(j*9+1,8)<<" as "<<tmp_group<<endl;
                  }
	        else if (outFlux==2)
	          {
                    fileFlux <<j<<" "<<name_group(j)<<" read from "<<surrogate[i].smile.substr(j*9+1,8)<<" as "<<tmp_group<<endl;
	          }
*/	          
                surrogate[i].groups[j]=tmp_group;
              }
            else if (tmp_group < 0.)
              {
                cout <<"Warning!!! aerosol vector input < 0.!!! set to 0. Please check aerosol species list."<<endl;
                cout <<j<<" "<<name_group(j)<<" read from "<<surrogate[i].smile.substr(j*9+1,8)<<" as "<<endl;
/*                if (outFlux==1)
	          {
                    cout <<"Warning!!! aerosol vector input < 0.!!! set to 0. Please check aerosol species list."<<endl;
                    cout <<j<<" "<<name_group(j)<<" read from "<<surrogate[i].smile.substr(j*9+1,8)<<" as "<<endl;
                  }
	        else if (outFlux==2)
	          {
                    fileFlux <<"Warning!!! aerosol vector input < 0.!!! set to 0. Please check aerosol species list."<<endl;
                    fileFlux <<j<<" "<<name_group(j)<<" read from "<<surrogate[i].smile.substr(j*9+1,8)<<" as "<<endl;
	          }
*/
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
        if (surrogate[i].is_monomer and surrogate[i].is_organic)
	    {
	      int j=surrogate[i].ioligo;	
	      for (int igr=0;igr<60;igr++)
	        surrogate[j].groups[igr]=surrogate[i].moligo*surrogate[i].groups[igr];
	    }
      }
}

