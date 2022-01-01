/*
"Small protein complex prediction algorithm based on protein-protein interaction network segmentation" IEEE/ACM Transactions on Computational Biology and Bioinformatics (TCBB) The project has been submitted and is under review.
Developer:
JiaqingLv Dalian University of Technology jiaqinglv@foxmail.com
ZhenYao Dalian University of Technology 22151303@zju.edu.cn
BingLiang Dalian University of Technology liangbing@dlut.edu.cn
YijiaZhang Dalian University of Technology zhyj@dmu.edu.cn
*/ 
#include <cstdio>
#include <vector>
#include <algorithm>
#include <map>
#include <iostream>
#include <set>
#include <fstream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <queue>
#include <cstring>
using namespace std;
const int MAXN = 20010,MAXP = 20;
struct Interaction
{
    int proteina,proteinb;
    double interaction,balanced_interaction;
    Interaction (int _proteina = 0,int _proteinb = 0,double _interaction = 0.0,double _balanced_interaction = 0.0):proteina(_proteina),proteinb(_proteinb),interaction(_interaction),balanced_interaction(_balanced_interaction){}
    friend bool operator < (Interaction a,Interaction b)
    {
    	return a.balanced_interaction < b.balanced_interaction;
    }
};

struct PPI
{
    vector <int> protein;
    vector <Interaction> interaction;
    bool mp[MAXP][MAXP]; 
};

struct Result
{
    double cohesion;
    vector <int> protein;
    friend bool operator < (Result a,Result b)
    {
    	return a.cohesion > b.cohesion;
    }
};

map <string,int> Protein2id;
string Id2protein[MAXN];
int Protein_count;

void read_proteins(PPI &Current_ppi,string Ppidata_file)
{//This function reads PPI from Ppidata_file
    ifstream fin(Ppidata_file);
    string Protein_a,Protein_b;
    double Protein_interaction;
    set <int> Protein_set;
    while (fin >> Protein_a >> Protein_b >>  Protein_interaction)
    {
	    if (Protein2id.count(Protein_a) == 0)
        {
            Protein2id[Protein_a] = ++Protein_count;
            Id2protein[Protein_count] = Protein_a;
        }
        if (Protein2id.count(Protein_b) == 0)
        {
            Protein2id[Protein_b] = ++Protein_count;
            Id2protein[Protein_count] = Protein_b;
        }
        if (Protein_set.find(Protein2id[Protein_a]) == Protein_set.end())
        {
            Protein_set.insert(Protein2id[Protein_a]);
            Current_ppi.protein.push_back(Protein2id[Protein_a]);
        }
        if (Protein_set.find(Protein2id[Protein_b]) == Protein_set.end())
        {
            Protein_set.insert(Protein2id[Protein_b]);
            Current_ppi.protein.push_back(Protein2id[Protein_b]);
        }
        Current_ppi.interaction.push_back(Interaction(Protein2id[Protein_a],Protein2id[Protein_b],Protein_interaction));
    }
    fin.close();
    return;
}

void get_balanced_interaction(PPI &Current_ppi,double Balanced_index)
{//This function calculates the balanced weight of interaction
    map <int,double> Sum;
	for (int i = 0;i < Current_ppi.interaction.size();i++)
	{
	    Sum[Current_ppi.interaction[i].proteina] += Current_ppi.interaction[i].interaction;
	    Sum[Current_ppi.interaction[i].proteinb] += Current_ppi.interaction[i].interaction;
	}
	for (int i = 0;i < Current_ppi.interaction.size();i++)
	{
	    Current_ppi.interaction[i].balanced_interaction = pow(Current_ppi.interaction[i].interaction,Balanced_index) / pow(Sum[Current_ppi.interaction[i].proteina],Balanced_index - 1) + pow(Current_ppi.interaction[i].interaction,Balanced_index) / pow(Sum[Current_ppi.interaction[i].proteinb],Balanced_index - 1);
	}
    return;
}

int getfa(int x,int fa[])
{//This is a data structure, called Disjoint Set Union to check if two proteins are in the same set
    int a = x;
    while (x != fa[x]) 
	{
        x = fa[x];
    }
    while (a != fa[a]) 
	{
        int z = a;
        a = fa[a];
        fa[z] = x;
    }
    return x;
}

void split_ppi(queue <PPI> &Ppi_queue,vector <PPI> &Splitted_ppi)
{//This function is to split huge ppi into small ppi
    PPI Current_ppi = Ppi_queue.front();
	Ppi_queue.pop();
	
    int Count = 0,location = -1;
    if (Current_ppi.protein.size() <= MAXP)
    {
        Splitted_ppi.push_back(Current_ppi);
        return;
    }
    
    int fa[MAXN + 1];
    //Initializes Disjoint Set Union
    for (int i = 1;i <= MAXN;i++)
    {
        fa[i] = i;
    }
    
    sort(Current_ppi.interaction.begin(),Current_ppi.interaction.end());
    for (int i = Current_ppi.interaction.size() - 1;i >= 0;i--)
    {
        if (getfa(Current_ppi.interaction[i].proteina,fa) == getfa(Current_ppi.interaction[i].proteinb,fa))
        {
            continue;
        }
        fa[getfa(Current_ppi.interaction[i].proteina,fa)] = getfa(Current_ppi.interaction[i].proteinb,fa);
        Count++;
        if (Count == Current_ppi.protein.size() - 2)
        {
            break;
        }
    }
    
    while (true)
    {
        bool success = false;
        PPI New_ppi;
        set <int> Protein_set;
        for (int i = location + 1;i < Current_ppi.protein.size();i++)
        {
            if (getfa(Current_ppi.protein[i],fa) == Current_ppi.protein[i])
            {
                location = i;
                success = true;
                break;
            }
        }
        if (success == false)
            break;
        for (int i = 0;i < Current_ppi.protein.size();i++)
        {
            if (getfa(Current_ppi.protein[i],fa) == Current_ppi.protein[location])
            { 
                New_ppi.protein.push_back(Current_ppi.protein[i]);
                Protein_set.insert(Current_ppi.protein[i]);
            }
        }
        for (int i = 0;i < Current_ppi.interaction.size();i++)
        {
            if (Protein_set.find(Current_ppi.interaction[i].proteina) != Protein_set.end() && Protein_set.find(Current_ppi.interaction[i].proteinb) != Protein_set.end())
            {
                New_ppi.interaction.push_back(Current_ppi.interaction[i]);
            }
        }
        Ppi_queue.push(New_ppi);
    }
    return;
}

double calculate_similarity(vector <int> Complexa,vector <int> Complexb)
{//This function caculate the similarity of two complexs
    set <int> Complexb_set;
    int count = 0;
    
    for (int i = 0;i < Complexb.size();i++)
    {
        Complexb_set.insert(Complexb[i]);
    }
    for (int i = 0;i < Complexa.size();i++)
    {
        if (Complexb_set.find(Complexa[i]) != Complexb_set.end())
        {
		    count += 1;
		}
    }
    return (double)count / max(Complexa.size(),Complexb.size());
}

double calculate_complex_cohesion(PPI &Current_ppi,vector <int> &Complex)
{//This function caculate the cohesion of one complex
	set <int> Complex_set;
    map <int,double> Sum;
    map <int,int> Count;
    double Cohesion = 0.0;
    
    for (int i = 0;i < Complex.size();i++)
    {
        Complex_set.insert(Complex[i]);
    }
    for (int i = 0;i < Current_ppi.interaction.size();i++)
    {
        if (Complex_set.find(Current_ppi.interaction[i].proteina) == Complex_set.end() || Complex_set.find(Current_ppi.interaction[i].proteinb) == Complex_set.end())
        {
            continue;
    	}
        Count[Current_ppi.interaction[i].proteina] += 1;
        Count[Current_ppi.interaction[i].proteinb] += 1;
        Sum[Current_ppi.interaction[i].proteina] += Current_ppi.interaction[i].balanced_interaction;
        Sum[Current_ppi.interaction[i].proteinb] += Current_ppi.interaction[i].balanced_interaction;
    }
	for (int i = 0;i < Complex.size();i++)
	{
	    Cohesion += Sum[Complex[i]] * (Count[Complex[i]] + 1) / Complex.size();
	}
	Cohesion /= Complex.size();
    return Cohesion;
}



void update_result(Result &Complex,vector <Result> &result,double Similarity_threshold)
{//This function determines whether the connected subset is a protein complex and removes protein complexes that are too similar.
    for (int i = 0;i < result.size();i++)
    {
        if (calculate_similarity(result[i].protein,Complex.protein) >= Similarity_threshold)
        {
            return;
        }
    }
    result.push_back(Complex);
    
    return;
}

int get_Current_ppi_protein(PPI &Current_ppi,int x)
{//This function gets the corresponding serial number of the protein in the current PPI
	return lower_bound(Current_ppi.protein.begin(),Current_ppi.protein.end(),x) - Current_ppi.protein.begin();
}

void get_complexs(PPI &Current_ppi,vector <Result> &result,double Similarity_threshold)
{//This function enumerates the connected subset of each small PPIN
	bool Connectivity[MAXP][MAXP];
	map <int,bool> Complex_record;
	queue <int> Complex_queue;
	memset(Connectivity,false,sizeof(Connectivity));
	
	sort(Current_ppi.protein.begin(),Current_ppi.protein.end());
	for (int i = 0;i < Current_ppi.interaction.size();i++)
	{
		int proteina_id = get_Current_ppi_protein(Current_ppi,Current_ppi.interaction[i].proteina);
		int proteinb_id = get_Current_ppi_protein(Current_ppi,Current_ppi.interaction[i].proteinb);
		Connectivity[proteina_id][proteinb_id] = Connectivity[proteinb_id][proteina_id] = true;
	}
	for (int i = 0;i < Current_ppi.protein.size();i++)
	{
		Complex_queue.push(1 << i);
		Complex_record[1 << i] = true;
	} 
	while (Complex_queue.size() != 0)
	{
		int Complex_state = Complex_queue.front();
		Complex_queue.pop();
		for (int i = 0;i < Current_ppi.protein.size();i++)
		{
			if ((Complex_state & (1 << i)) != 0)
			{
				for (int j = 0;j < Current_ppi.protein.size();j++)
				{
					if (Connectivity[i][j] == false)
					{
						continue;
					}
					if ((Complex_state & (1 << j)) != 0)
					{
						continue;
					}
					if (Complex_record.count(Complex_state | (1 << j)) != 0)
					{
						continue;
					}
					Complex_record[Complex_state | (1 << j)] = true;
					Complex_queue.push(Complex_state | (1 << j));
				}
			}
		}
	}
	
	vector <int> Complex;
	map <int,bool> :: iterator it;
	vector <Result> Complex_result;
	it = Complex_record.begin();
	while (it != Complex_record.end())
	{
		Complex.clear();
		for (int i = 0;i < Current_ppi.protein.size();i++)
		{
			if ((it->first & (1 << i)) != 0)
			{ 
				Complex.push_back(Current_ppi.protein[i]);
			} 
		}
        if (Complex.size() >= 3) 
		{
			Result tmp;
			tmp.protein = Complex;
			tmp.cohesion = calculate_complex_cohesion(Current_ppi,Complex);
			Complex_result.push_back(tmp);
		}
        it++;
	}
	
	if (Complex_result.size() == 0)
		return;
	sort(Complex_result.begin(),Complex_result.end());
	for (int i = 0;i < Complex_result.size();i++)
	{
		update_result(Complex_result[i],result,Similarity_threshold);
	}
	return;
}

vector <Result> get_result(PPI &Current_ppi,double Similarity_threshold)
{//This function divides the PPIN into smaller PPIN and calls the function to solve each smaller PPIN
    queue <PPI> Ppi_queue;
    vector <PPI> Splitting_ppi;
    
    Ppi_queue.push(Current_ppi);
    while (!Ppi_queue.empty())
    {
		split_ppi(Ppi_queue,Splitting_ppi);     
    }
    
    vector <Result> result;
    for (int i = 0;i < Splitting_ppi.size();i++)
    {
        vector <Result> Temporary_res;
        get_complexs(Splitting_ppi[i],Temporary_res,Similarity_threshold); 
      
		for (int j = 0;j < Temporary_res.size();j++)
        {
            result.push_back(Temporary_res[j]);
        }
    }
    return result;
}
void write_proteins(vector <Result> Result_complex,string Result_file)
{////This function writes predicted protein complexes in result
  	sort(Result_complex.begin(),Result_complex.end()); 
	ofstream fout(Result_file);
	int Tmpk = 0;
	string Tmps;
	    
	int Threshold_complex = min(int(Result_complex.size() * 0.05),int(Result_complex.size() - 1));
	double Cohesion_threshold = Result_complex[Threshold_complex].cohesion * 0.5;
		
    for (int i = 0;i < Result_complex.size();i++)
    {
        if (Result_complex[i].protein.size() <= 1)
        {
            continue;
        }
        if (Result_complex[i].cohesion < Cohesion_threshold)
        {
        	break;
        }
        Tmpk += 1;
        Tmps = "";
        for (int j = 0;j < Result_complex[i].protein.size();j++)
        {
            Tmps = Tmps + " " + Id2protein[Result_complex[i].protein[j]];
        }
        fout << Tmps << endl;
    }
    fout.close();
    return;
}



void print_information(string Ppidata_file,string Result_file,double Balanced_index)
{
	printf("_________________________________________\n");
	cout << "The PPI_file is "<< Ppidata_file << endl;
	cout << "The result_file is "<< Result_file << endl;
	printf("The balanced index is %.3lf\n",Balanced_index);
	printf("_________________________________________\n");
	printf("It will takes tens of minutes.\n");
}

int main(int argc, char *argv[])
{
	double Balanced_index,Similarity_threshold = 0.45;
	string Ppidata_file,Result_file;
	PPI Original_ppi;
	
	Ppidata_file = argv[1];
	Result_file = argv[2];
	read_proteins(Original_ppi,Ppidata_file);
	if (argc >= 4)
	{
		sscanf(argv[3],"%lf",&Balanced_index);
	}else
	{
		Balanced_index = 1.5;
	}
	
	print_information(Ppidata_file,Result_file,Balanced_index);
    get_balanced_interaction(Original_ppi,Balanced_index);
    write_proteins(get_result(Original_ppi,Similarity_threshold),Result_file);
    return 0;
}

