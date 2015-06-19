#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <ctime>

using namespace std;

const int SEED = 5;
const int MAX = 2501; // N'=(2m+x)N
const int m = 2;
const int x = 1; 
const float y = 0; // add ymk extra links
const float e = 0;
struct network
{
    int k;
    int edge[MAX];
    int geodis[MAX];
    float c;
};
network node[MAX];
int M[MAX][MAX];
int M0[MAX][MAX];
int N=0;

float myran();
void printnet(network *);
void init();
void grow(int n);// n times growth

int search(int, int, network *);
int neighbor(int, int, network *);
float pathl2(network *); //method 2; improved.
float cluster(network *);
float pathl1(network *);  //method 1; better!

int main()
{
    srand(time(NULL));
    //cout << rand()%N << endl;

    ofstream outFile;
    outFile.open("mode2.txt");
    outFile << "SEED= " << SEED << endl;
    outFile << "m= " << m << endl;
    outFile << "x= " << x << endl;
    outFile << "y= " << y << endl;
    outFile << "e= " << e << endl;

    init();//initialize node[MAX] & M[MAX][MAX]
    //printnet(node);
    cout << "initialization complete!" << endl;

    grow(4);
    //grow(1);
    outFile << "N= " << N << endl;


    //printnet(node);
    //outFile << "*vertices " << N << endl;
    //outFile << "*edges " << endl;

    for (int i=0;i<N;i++)
        for (int j=i+1;j<N;j++)
            if (M[i][j]==1)
            outFile << i << "  " << j << endl;

    //cout << "C = " << cluster(node) << endl;
    //cout << "L = " << pathl1(node) << endl;

    outFile.close();
    return 0;
}

void grow(int n)
{
    int q=0;
    for (int i=0;i<MAX;i++)
        if (node[i].k)
            q++;
    N=q;
    cout <<  "# of nodes: " << N << endl;

    int y1,y2;
    for (int u=0;u<n;u++)
    {
        for (int i=0; i<N;i++)
            for (int j=0; j<N; j++) 
                M0[i][j]=M[i][j];

        for (int i=0;i<N;i++)
            if (node[i].k)
            {
                node[i].k*=m;
                for (int j=0;j<node[i].k;j++)
                {
                    node[i].edge[j]=q;
                    M[i][q]=M[q][i]=1;
                    q++;
                }
                for (int t=0;t<int(y*node[i].k);t++)
                {
                    y1=rand()%node[i].k;
                    y2=rand()%node[i].k;
                    if (y1==y2)
                        t--;
                    else if (M[node[i].edge[y1]][node[i].edge[y2]]==1)
                        t--;
                    else
                    {
                        M[node[i].edge[y1]][node[i].edge[y2]]=1;
                        M[node[i].edge[y2]][node[i].edge[y1]]=1;
                    }
                }
                node[i].k=0;
            }

        for (int i=0; i<N;i++)
                for (int j=i+1; j<N; j++)
                    if (M0[i][j]==1)
                        M[i][j]=M[j][i]=-1;

        for (int i=0; i<q;i++){
            for (int j=i+1; j<q; j++) 
                if (M[i][j]==1){
                    node[i].edge[node[i].k]=j;
                    node[j].edge[node[j].k]=i;
                    node[i].k++;
                    node[j].k++;
                }
        }

        int f,g;
        for (int i=0; i<N;i++)
            {
                for (int j=i+1; j<N; j++)
                    if (M0[i][j]==1)
                    {
                        for (int t=0;t<x;t++)
                        {
                            f=rand()%(node[i].k);
                            g=rand()%(node[j].k);
                            if (M[node[i].edge[f]][node[j].edge[g]]==1)
                                t--;
                            else
                            {
                                M[node[i].edge[f]][node[j].edge[g]]=1;
                                M[node[j].edge[g]][node[i].edge[f]]=1;
                            }
                        }
                    }
            }
        
        for (int i=0; i<q;i++)
            node[i].k=0;

        for (int i=0; i<q;i++)
            for (int j=i+1; j<q; j++) 
                if (M[i][j]==1){
                    node[i].edge[node[i].k]=j;
                    node[j].edge[node[j].k]=i;
                    node[i].k++;
                    node[j].k++;
                }
        N=q;
        cout << "After " << u+1 << "-time growth, # of nodes : " << N << endl;

        int sum_k=0;
        for (int i=0;i<N;i++)
            sum_k+=node[i].k;
        sum_k/=2;
        cout << "# of links : " << sum_k << endl;

       }
    
}

float myran()
{
    return float(rand())/RAND_MAX;
}

// --- pathl1()

float pathl1(network * node)
{
    for (int i=0; i<N;i++)
        for (int j=0; j<N; j++) {
            M[i][j] = -1;
            if( i == j) M[i][j] = 0;
        }
    for (int i=0; i<N; i++) 
        for( int t=0; t<node[i].k; t++) {
            int j =  node[i].edge[t]; 
            M[i][j] = 1;
        }
    int n = 2 ;
    bool flag=false;
    while (!flag) {
        flag=true;
        for (int i=0;i<N;i++)
            for (int j=0;j<N;j++)
                if(M[i][j] < 0)
                {
                    flag=false;
                    for(int k=0; k<N; k++)
                        if(M[i][k]+M[k][j]==n)
                        {
                            M[i][j]=n;
                            break;
                        }
                }
        if (!flag){ 
            //cout << n << endl;
            n++;
        } else n--;
    }
    cout << "The largest shortest path is " << n << "." << endl;
    int sum_dis=0;
    for (int i=0; i<N;i++)
        for (int j=0; j<N; j++) {
            node[i].geodis[j] = M[i][j];
            sum_dis+=M[i][j];
        }
    float L=float(sum_dis)/(N*(N-1));
    return L;
}


// --- search(i,j)

int search(int i, int j, network * node)
{
    int dis=1;
    int found=0;
    bool flag;
    int candidate[N];
    for (int c=0;c<N;c++)
    {
        candidate[c]=-1;
    }
    candidate[0]=i;
    //for (int c=0;c<N;c++)
    //{
      //  cout << candidate[c] << " ";
    //}
    found=neighbor(i,j,node);
    while (found==0)
    {
        int candidate0[N];
        for (int c=0;c<N;c++)
        {
            candidate0[c]=candidate[c];
            candidate[c]=-1;
        }
        int q=0;
        int m=0;
        while (candidate0[q]!=-1) 
        {
            int p=0;
            while (node[candidate0[q]].edge[p]!=-1) 
            {   
                flag=true;
                for (int u=0;u<m;u++)
                    if (candidate[u]==node[candidate0[q]].edge[p])
                    {
                        flag=false;
                        break;
                    }

                if (flag==true) 
                {
                    candidate[m]=node[candidate0[q]].edge[p];
                    m++;
                }
                p++;
            }
            q++;
        }
        for (int d=0;d<m;d++){
            found+=neighbor(candidate[d],j,node);
            if (node[i].geodis[candidate[d]]==0)
            node[i].geodis[candidate[d]]=node[candidate[d]].geodis[i]=dis;
        }
        dis++;
    }
    return dis;
}

// --- neighbor(i,j)

int neighbor(int i, int j, network * node)
{
    bool flag=false;
    for (int p=0;p<node[i].k;p++)
    {
        if (node[i].edge[p]==-1) break;
        if (node[i].edge[p]==j)
            flag=true;
    }
    if (flag==true)
        return 1;
    else return 0; 
}

// --- calculate L(0) (=averaged geodesic distance)

float pathl2(network * node)
{
    int sum_dis=0;
    for (int i=0;i<N;i++)
    {     
        for (int j=i+1;j<N;j++)
        {   
            if (neighbor(i,j,node))
                node[i].geodis[j]=node[j].geodis[i]=1;
            else if (node[i].geodis[j]==0){
            node[i].geodis[j]=node[j].geodis[i]=search(i,j,node);
            }
            sum_dis+=2*node[i].geodis[j];
        }
    }
    float L=float(sum_dis)/(N*(N-1));

    //display dis[i][j] & L(0)

    /*
    for (int i=0;i<N;i++)
    {
        for (int j=0;j<N;j++)
            cout << node[i].geodis[j] << " ";
        cout << endl;
    }*/
    return L;
}

// --- calculate C(0) (=averaged clustering coefficient)

float cluster(network * node)
{
    float sum_c=0;
    for (int i=0;i<N;i++)
        if (node[i].k>1)
        {
            float e=0;
            for (int j=0;j<node[i].k;j++)
                for (int t=j+1;t<node[i].k;t++)
                    e+=neighbor(node[i].edge[t],node[i].edge[j],node);
            node[i].c=e/(node[i].k*(node[i].k-1)/2);
            sum_c+=node[i].c;
        }
    float C=sum_c/N;

    //display C(0)

    
    /*
    for (int i=0;i<N;i++){
        cout << "Node " << i << "'s clustering coefficient is ";
        cout << node[i].c << "." << endl;}
        */
    
    return C;
}



// --- initialization --- node[MAX],M[MAX][MAX],M0[MAX][MAX]

void init()
{ 
    for (int i=0;i<MAX;i++)
    {
        node[i].k = 0;
        node[i].c = 0;
        for (int j=0;j<MAX;j++)
        {
            node[i].edge[j] = -1;
            node[i].geodis[j] = 0;
        }  
    }
    
    for (int i=1;i<SEED;i++){
        node[0].edge[i-1]=i;
        node[i].edge[0]=0;
        node[i].k=1;
    }
    node[0].k=SEED-1;

    for (int i=0; i<MAX;i++)
        for (int j=0; j<MAX; j++) {
            M[i][j] = -1;
            if( i == j) M[i][j] = 0;
        }
    for (int i=0; i<SEED; i++) 
        for( int t=0; t<node[i].k; t++) {
            int j = node[i].edge[t]; 
            M[i][j] = 1;
        }
}

void printnet(network * node)    
{
    for (int i=0;i<N;i++)
    {   
        cout << "Node " << setw(2) << i << " has " << node[i].k << " edges:"; 
        for (int j=0;j<node[i].k;j++)
            cout << setw(3) << node[i].edge[j] << "  ";
        cout << endl;
    }      
}
