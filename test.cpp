#include <bits/stdc++.h>
using namespace std;
#define m 3
#define M 5
#define D 10
struct node
{
    int type; // 0 for point , 1 for leaf node, 2 for normal rectangle node
    int n;
    double *L;
    double *R;
    node **children;
};
struct node_list
{
    node *N;
    node_list *next;
};
struct duo
{
    node *N1;
    node *N2;
};
double max(double a, double b)
{
    if (a > b)
    {
        return a;
    }
    return b;
}
double min(double a, double b)
{
    if (a < b)
    {
        return a;
    }
    return b;
}
node *create_node()
{
    node *N = new node[1];
    N->n = 0;
    N->type = 0;
    N->L = new double[D];
    N->R = new double[D];
    N->children = new node *[M + 1];
    return N;
}
bool overlap(node *N1, node *N2)
{
    for (int i = 0; i < D; i++)
    {
        if (((N1->R)[i] < (N2->L)[i]) || (N1->L)[i] > (N2->R)[i])
        {
            return false;
        }
    }
    return true;
}
node_list *insert(node_list *L, node *N)
{
    L->N = N;
    L->next = new node_list[1];
    L->next->next = NULL;
    return L->next;
}
void search(node *N, node *Q, node_list **L)
{
    if (N->type == 1)
    {
        for (int i = 0; i < N->n; i++)
        {
            if (overlap((N->children)[i], Q))
            {
                (*L) = insert((*L), (N->children)[i]);
            }
        }
        return;
    }
    for (int i = 0; i < N->n; i++)
    {
        if (overlap((N->children)[i], Q))
        {
            search(N, Q, L);
        }
    }
    return;
}
double area_increase(node *N, node *Q)
{
    double area = 1;
    for (int i = 0; i < D; i++)
    {
        area *= ((N->R)[i] - (N->L)[i]);
    }
    double area2 = 1;
    for (int i = 0; i < D; i++)
    {
        area2 *= ((max((N->R)[i], (Q->R)[i])) - (min((Q->L)[i], (N->L)[i])));
    }
    return area2 - area;
}
double area(node *N)
{
    double area = 1;
    for (int i = 0; i < D; i++)
    {
        area *= ((N->R)[i] - (N->L)[i]);
    }
    return area;
}
duo quad_split(node *N)
{
    int N1 = 0;
    int N2 = 1;
    double A_max = area_increase((N->children)[N1], (N->children)[N2]) - area((N->children)[N2]);
    for (int i = 0; i < N->n; i++)
    {
        for (int j = i + 1; j < N->n; j++)
        {
            double A = area_increase((N->children)[i], (N->children)[j]) - area((N->children)[j]);
            if (A > A_max)
            {
                N1 = i;
                N2 = j;
                A = A_max;
            }
        }
    }
    int size1 = 1;
    int size2 = 1;
    node *New1 = create_node();
    node *New2 = create_node();
    New1->type = N->type;
    New2->type = N->type;
    (New1->children)[0] = (N->children)[N1];
    (New2->children)[0] = (N->children)[N2];
    for (int i = 0; i < D; i++)
    {
        (New1->L)[i] = ((New1->children)[0]->L)[i];
        (New1->R)[i] = ((New1->children)[0]->R)[i];
        (New2->L)[i] = ((New2->children)[0]->L)[i];
        (New2->R)[i] = ((New2->children)[0]->R)[i];
    }
    for (int i = 0; i < N->n; i++)
    {
        if (i != N1 && i != N2)
        {
            if ((N->n) - size2 == m)
            {
                (New1->children)[size1++] = (N->children)[i];
                for (int j = 0; j < D; j++)
                {
                    (New1->L)[j] = min(((New1->children)[size1 - 1]->L)[j], (New1->L)[j]);
                    (New1->R)[j] = max(((New1->children)[size1 - 1]->R)[j], (New1->R)[j]);
                }
            }
            else if ((N->n) - size1 == m)
            {
                (New2->children)[size2++] = (N->children)[i];
                for (int j = 0; j < D; j++)
                {
                    (New2->L)[j] = min(((New1->children)[size1 - 1]->L)[j], (New1->L)[j]);
                    (New2->R)[j] = max(((New1->children)[size1 - 1]->R)[j], (New1->R)[j]);
                }
            }
            else
            {
                double delta_area_1 = area_increase((N->children)[i], (New1->children)[0]);
                double delta_area_2 = area_increase((N->children)[i], (New2->children)[0]);
                if (delta_area_1 < delta_area_2)
                {
                    (New1->children)[size1++] = (N->children)[i];
                    for (int j = 0; j < D; j++)
                    {
                        (New1->L)[j] = min(((New1->children)[size1 - 1]->L)[j], (New1->L)[j]);
                        (New1->R)[j] = max(((New1->children)[size1 - 1]->R)[j], (New1->R)[j]);
                    }
                }
                else
                {
                    (New2->children)[size2++] = (N->children)[i];
                    for (int j = 0; j < D; j++)
                    {
                        (New2->L)[j] = min(((New2->children)[size2 - 1]->L)[j], (New2->L)[j]);
                        (New2->R)[j] = max(((New2->children)[size2 - 1]->R)[j], (New2->R)[j]);
                    }
                }
            }
        }
    }
    (New1->n) = size1;
    (New2->n) = size2;
    duo Two;
    Two.N1 = New1;
    Two.N2 = New2;
    delete N;
    return Two;
}
duo insert(node *N, node *Q)
{
    if (N->type == 1)
    {
        (N->children)[N->n] = Q;
        (N->n)++;
        for (int i = 0; i < D; i++)
        {
            (N->L)[i] = min((Q->L)[i], (N->L)[i]);
            (N->R)[i] = max((Q->R)[i], (N->R)[i]);
        }
        duo Two;
        if (N->n > M)
        {
            Two = quad_split(N);
        }
        else
        {
            Two.N1 = N;
            Two.N2 = NULL;
        }
        return Two;
    }
    if (N->type == 0)
    {
        cout << "LOL" << endl;
    }
    double delta_area_min = area_increase((N->children)[0], Q);
    int N1 = 0;
    for (int i = 1; i < N->n; i++)
    {
        double delta_area = area_increase((N->children)[i], Q);
        if (delta_area < delta_area_min)
        {
            delta_area_min = delta_area;
            N1 = i;
        }
    }
    duo Two = insert((N->children)[N1], Q);
    (N->children)[N1] = Two.N1;
    if (Two.N2 != NULL)
    {
        (N->children)[(N->n)++] = Two.N2;
    }
    for (int i = 0; i < D; i++)
    {
        double Left = ((N->children)[0]->L)[i];
        double Right = ((N->children)[0]->R)[i];
        for (int j = 1; j < N->n; j++)
        {
            Left = min(Left, ((N->children)[j]->L)[i]);
            Right = max(Right, ((N->children)[j]->R)[i]);
        }
        (N->L)[i] = Left;
        (N->R)[i] = Right;
    }
    if (N->n > M)
    {
        return Two = quad_split(N);
    }
    else
    {
        Two.N1 = N;
        Two.N2 = NULL;
    }
    return Two;
}
node *insert_main(node *N, node *temp)
{
    duo Two = insert(N, temp);
    if (Two.N2 != NULL)
    {
        N = create_node();
        N->type = 2;
        N->n = 2;
        (N->children)[0] = Two.N1;
        (N->children)[1] = Two.N2;
    }
    else
    {
        N = Two.N1;
    }
    return N;
}
double dist(node *Q1, node *Q2, int type)
{
    if (type == 1)
    {
        double dist = 0;
        for (int i = 0; i < D; i++)
        {
            dist += abs((Q1->L)[i] - (Q2->L)[i]);
        }
        return dist;
    }
    if (type == 2)
    {
        double dist = 0;
        for (int i = 0; i < D; i++)
        {
            dist += pow(abs((Q1->L)[i] - (Q2->L)[i]), 2);
        }
        return sqrt(dist);
    }
    if(type == 0)
    {
        double dist = abs((Q1->L) - (Q2->L));
        for (int i = 1; i < D; i++)
        {
            dist = max(dist, abs((Q1->L)[i] - (Q2->L)[i]));
        }
        return dist;
    }
    return 0;
}

double mindist(node *N, node *Q, int type)
{
    node* Q1 = create_node();
    Q1 -> type = 0;
    for(int i = 0; i < D; i++)
    {
        if((Q -> L)[i] < (N -> L)[i])
        {
            (Q1 -> L)[i] = (N -> L)[i];
        }
        else if((Q -> L)[i] > (N -> R)[i])
        {
            (Q1 -> L)[i] = (N -> R)[i];
        }
        else 
        {
            (Q1 -> L)[i] = (Q -> L)[i];
        }
    }
    double Dist = dist(Q1,Q,type);
    delete Q1;   
    return Dist;
}
double minmaxdist(node *N, node *Q, int type)
{
    node* Q1 = create_node();
    Q1 -> type = 0;
    node* Q2 = create_node();
    Q2 -> type = 0;
    for(int i = 0; i < D; i++)
    {
        if(2*(Q -> L)[i] <= (N -> L)[i] + (N -> R)[i])
        {
            (Q1 -> L)[i] = (N -> L)[i];
            (Q1 -> R)[i] = (N -> R)[i]; 
        }
        else
        {
            (Q1 -> L)[i] = (N -> R)[i];
            (Q1 -> R)[i] = (N -> L)[i];
        }
    }
    double Dist;
    for(int j = 0; j < D; j++)
    {
        if(j != 0)
        {
            (Q2 -> L)[j] = (Q1 -> R)[j];
        }
        else
        {
            (Q2 -> L)[j] = (Q1 -> L)[j];
        }
         
    }
    Dist = dist(Q,Q2,type);
    for(int i = 0; i < D; i++)
    {
        for(int j = 0; j < D; j++)
        {
            if(j != i)
            {
                (Q2 -> L)[j] = (Q1 -> R)[j];
            }
            else
            {
                (Q2 -> L)[j] = (Q1 -> L)[j];
            }
        }
        Dist = min(Dist,dist(Q,Q2,type));
    }
    delete Q1;
    delete Q2;
    return Dist;
}
void Nearest_Neighbour(node* N,node* Q,double* Dist,node** NN,int type)
{
    if(N -> type == 1)
    {
        double Minmax = minmaxdist(N,Q,type);
        double Min = mindist(N,Q,type);
        if(Min >= (*Dist))
        {
            return;
        }
        for(int i = 0; i < N -> n; i++)
        {
            double temp = dist((N -> children)[i],Q,type);
            if(temp < (*Dist))
            {
                *NN = (N -> children)[i];
                (*Dist) = temp;
            }
        }
        return;
    }
    double low_Minmax = minmaxdist((N -> children)[0],Q,type);
    int N1 = 0;
    double Minmax[N -> n];
    double Min[N -> n];
    for(int i = 0; i < N -> n; i++)
    {
        Minmax[i] = minmaxdist((N->children)[i],Q,type);
        if(Minmax[i] < low_Minmax)
        {
            low_Minmax = Minmax[i];
            N1 = i;
        }
        Min[i] = mindist((N -> children)[i],Q,type);
    }
    for(int i = 0; i < N -> n; i++)
    {
        if(Min[i] > low_Minmax && i != N1)
        {
            continue;
        }
        else if(Min[i] > (*Dist))
        {
            continue;
        }
        else
        {
            Nearest_Neighbour((N->children)[i],Q,Dist,NN,type);
        }
    }
}
void dfs(node *N)
{
    if(N == NULL)
    {
        return;
    }
    if (N->type == 0)
    {
        for (int i = 0; i < D; i++)
        {
            cout << (N->L)[i] << " ";
        }
        cout << endl;
        return;
    }
    for (int i = 0; i < N->n; i++)
    {
        dfs(N->children[i]);
    }
    return;
}
int main()
{
    int n = 200;
    node* P[n];
    for(int i = 0; i < n; i++)
    {
        P[i] = create_node();
        P[i] -> type = 0;
        for(int j = 0; j < D; j++)
        {
            cin >> (P[i] -> R)[j];
            (P[i] -> L)[j] = (P[i] -> R)[j];
        } 
    }
    int q = 10;
    node* Q = create_node();
    Q -> type = 0;
    int t = 1;
    for(int i = 0; i < q; i++)
    {
        for(int k = 0; k < D; k++)
        {
            cin >> (Q -> L)[k];
            (Q -> R)[k] = (Q -> L)[k];
        }
        int N = 0;
        double Dist = dist(P[0],Q,t);
        for(int l = 1; l < n; l++)
        {
            double temp = dist(P[l],Q,t);
            if(temp < Dist)
            {
                N = l;
                Dist = temp;
            }
        }
        dfs(Q);
        cout << Dist << ": ";
        dfs(P[N]);
    }
    return 0;
}
