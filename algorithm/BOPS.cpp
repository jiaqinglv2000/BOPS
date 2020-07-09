/*
"Small protein complex prediction algorithm based on protein-protein interaction network segmentation" IEEE/ACM Transactions on Computational Biology and Bioinformatics (TCBB) The project has been submitted and is under review.
Developer:
JiaqingLv Dalian University of Technology jiaqinglv@foxmail.com
ZhenYao Dalian University of Technology yaozhen@mail.dlut.edu.cn
BingLiang Dalian University of Technology liangbing@dlut.edu.cn
YijiaZhang Dalian University of Technology zhyj@dlut.edu.cn
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

void update_result(PPI &Current_ppi,vector <int> &Complex,vector <Result> &result,double Similarity_threshold,double Cohesion_threshold)
{//This function determines whether the connected subset is a protein complex and removes protein complexes that are too similar.
    double Cohesion = calculate_complex_cohesion(Current_ppi,Complex);
	int Count = 0,location = -1;
    if (Cohesion <= Cohesion_threshold)
    { 
        return;
    } 
    for (int i = 0;i < result.size();i++)
    {
        if (calculate_similarity(result[i].protein,Complex) >= Similarity_threshold)
        {
            Count += 1;
            location = i;
        }
    }
    if (Count >= 2)
    { 
        return;
    } 
    if (Count == 1)
    {
        if (Cohesion < result[location].cohesion)
        { 
            return;
        } 
        result[location].cohesion = Cohesion;
		result[location].protein = Complex;
        return;
    }
    if (Count == 0)
    {
        result.push_back(Result());
		result[result.size() - 1].cohesion = Cohesion;
		result[result.size() - 1].protein = Complex;
        return;
    }
    
    return;
}

int get_Current_ppi_protein(PPI &Current_ppi,int x)
{//This function gets the corresponding serial number of the protein in the current PPI
	return lower_bound(Current_ppi.protein.begin(),Current_ppi.protein.end(),x) - Current_ppi.protein.begin();
}

void get_complexs(PPI &Current_ppi,vector <Result> &result,double Similarity_threshold,double Cohesion_threshold)
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
			update_result(Current_ppi,Complex,result,Similarity_threshold,Cohesion_threshold);
		}
        it++;
	}
	return;
}

vector <Result> get_result(PPI &Current_ppi,double Similarity_threshold,double Cohesion_threshold)
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
        get_complexs(Splitting_ppi[i],Temporary_res,Similarity_threshold,Cohesion_threshold); 
		for (int j = 0;j < Temporary_res.size();j++)
        {
            result.push_back(Temporary_res[j]);
        }
    }
    return result;
}
void write_proteins(vector <Result> Result_complex,string Result_file)
{////This function writes predicted protein complexes in result
    ofstream fout(Result_file);
    int Tmpk = 0;
    string Tmps;
    for (int i = 0;i < Result_complex.size();i++)
    {
        if (Result_complex[i].protein.size() <= 1)
        {
            continue;
        }
        Tmpk += 1;
        Tmps = "C" + to_string(Tmpk) + ":";
        for (int j = 0;j < Result_complex[i].protein.size();j++)
        {
            Tmps = Tmps + " " + Id2protein[Result_complex[i].protein[j]];
        }
        fout << Tmps << endl;
    }
    fout.close();
    return;
}

void calculate_cohesion_threshold(PPI &Current_ppi,double &Cohesion_threshold,double ta = 0.03,double tb = -0.6)
{//This function calculates the similarity threshold based on the ppi density
	double Network_density = 2.0 * Current_ppi.interaction.size() / Current_ppi.protein.size() / (Current_ppi.protein.size() - 1);
	Cohesion_threshold = ta * pow(Network_density,tb);
	return;
}

void calculate_sim_threshold(PPI &Current_ppi,double &sim)
{//This function calculates the similarity threshold based on the ppi density
	double Network_density = 2.0 * Current_ppi.interaction.size() / Current_ppi.protein.size() / (Current_ppi.protein.size() - 1);
	if (Network_density >= 0.005)
		sim = 0.25;
	else
		sim = 0.45;
	return;
}

void print_information(string Ppidata_file,string Result_file,double Balanced_index,double Cohesion_threshold)
{
	printf("_________________________________________\n");
	cout << "The PPI_file is "<< Ppidata_file << endl;
	cout << "The result_file is "<< Result_file << endl;
	printf("The balanced index is %.3lf\n",Balanced_index);
	printf("The cohesion threshold is %.3lf\n",Cohesion_threshold);
	printf("_________________________________________\n");
	printf("It will takes tens of minutes.\n");
}

int main(int argc, char *argv[])
{
	double Balanced_index,Similarity_threshold,Cohesion_threshold,ta,tb;
	string Ppidata_file,Result_file;
	PPI Original_ppi;
	
	Ppidata_file = argv[1];
	Result_file = argv[2];
	read_proteins(Original_ppi,Ppidata_file);
	if (argc >= 4)
	{
		if (argv[3][0] == 'd')
			Balanced_index = 1.6;
		else 
			sscanf(argv[3],"%lf",&Balanced_index);;
		
	}
	if (argc >= 5)
	{
		if (argv[4][0] == 'd')
			calculate_cohesion_threshold(Original_ppi,Cohesion_threshold);
		else
			sscanf(argv[4],"%lf",&Cohesion_threshold);
	}
    
    if (argc < 4)
	{
		Balanced_index = 1.6;
	}
    if (argc < 5)
    {
		calculate_cohesion_threshold(Original_ppi,Cohesion_threshold);
	}
	calculate_sim_threshold(Original_ppi,Similarity_threshold);
	print_information(Ppidata_file,Result_file,Balanced_index,Cohesion_threshold);
    get_balanced_interaction(Original_ppi,Balanced_index);
    write_proteins(get_result(Original_ppi,Similarity_threshold,Cohesion_threshold),Result_file);
    return 0;
}

