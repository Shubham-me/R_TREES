#include <bits/stdc++.h>
#include <random>
using namespace std;
// Every node (other than root node) has between number of children between [m,M]
#define m 3
#define M 5

// D is the dimension of the data points
#define D 6

int metric; // 2 for L_2(Euclidian), 1 for L_1(Manhattan), 0 for L_inf(max Norm)
class node {
public:
    int type;                         // type of node
    int numberChildren;               // number of children
    vector<vector<double>> Rectangle; // size(2 x D): bounds of the rectangle the node represents [[l_1,l_2,..l_D],[u_1,u_2,..u_D]]
    vector<node *> children;          // children array

    node() : numberChildren(0), type(0) {
        Rectangle.resize(2, vector<double>(D));
    }
    node(int nC, int type) : numberChildren(nC), type(type) {
        children.resize(M + 1);
        Rectangle.resize(2, vector<double>(D));
    }

    // checks if 2 MBR overlap
    bool overlap(node &N1, node &N2) {
        for (int i = 0; i < D; i++) {
            if (N1.Rectangle[1][i] < N2.Rectangle[0][i]) {
                return false;
            }
            if (N2.Rectangle[1][i] < N1.Rectangle[0][i]) {
                return false;
            }
        }
        return true;
    }

    // searchs all the points in Tree inside the Query MBR Q
    void search(node &Q, vector<node *> &nodeList) {
        if (type == 1) {
            for (int i = 0; i < numberChildren; i++) {
                if (overlap(*children[i], Q)) {
                    nodeList.push_back(children[i]);
                }
            }
            return;
        }
        for (int i = 0; i < numberChildren; i++) {
            if (overlap(*children[i], Q)) {
                children[i]->search(Q, nodeList);
            }
        }
        return;
    }

    // calculates the area of the MBR
    double area() {
        double area = 1;
        for (int i = 0; i < D; i++) {
            area *= (Rectangle[1][i] - Rectangle[0][i]);
        }
        return area;
    }

    // calculates the area Increase if we add Q to current MBR
    double areaIncrease(node &Q) {
        double area = 1;
        for (int i = 0; i < D; i++) {
            area *= (Rectangle[1][i] - Rectangle[0][i]);
        }
        double areaNew = 1;
        for (int i = 0; i < D; i++) {
            areaNew *= (max(Rectangle[1][i], Q.Rectangle[1][i]) - min(Rectangle[0][i], Q.Rectangle[0][i]));
        }
        return areaNew - area;
    }

    // previous function overloaded, now input is a point instead of MBR
    double areaIncrease(vector<double> &P) {
        double area = 1;
        for (int i = 0; i < D; i++) {
            area *= (Rectangle[1][i] - Rectangle[0][i]);
        }
        double areaNew = 1;
        for (int i = 0; i < D; i++) {
            areaNew *= (max(Rectangle[1][i], P[i]) - min(Rectangle[0][i], P[i]));
        }
        return areaNew - area;
    }

    // splits the MBR into two sub-MBRs.
    node *split() {
        // For testing
        // Shift half of them to a new node
        int nC1 = numberChildren / 2;
        int nC2 = numberChildren - nC1;
        node *Sibling = new node(nC2, type);
        for (int i = 0; i < nC2; i++) {
            (Sibling->children)[i] = children[nC1 + i];
        }
        numberChildren = nC1;
        // redefine boundries for Sibling
        for (int i = 0; i < D; i++) {
            Sibling->Rectangle[0][i] = Sibling->children[0]->Rectangle[0][i];
            Sibling->Rectangle[1][i] = Sibling->children[0]->Rectangle[1][i];
            for (int j = 1; j < nC2; j++) {
                Sibling->Rectangle[0][i] = min(Sibling->Rectangle[0][i], Sibling->children[j]->Rectangle[0][i]);
                Sibling->Rectangle[1][i] = max(Sibling->Rectangle[1][i], Sibling->children[j]->Rectangle[1][i]);
            }
        }

        // redefine boundries for current
        for (int i = 0; i < D; i++) {
            Rectangle[0][i] = children[0]->Rectangle[0][i];
            Rectangle[1][i] = children[0]->Rectangle[1][1];
            for (int j = 1; j < nC1; j++) {
                Rectangle[0][i] = min(Rectangle[0][i], children[j]->Rectangle[0][i]);
                Rectangle[1][i] = max(Rectangle[1][i], children[j]->Rectangle[1][i]);
            }
        }

        return Sibling;
    }

    // inserts a point in MBR (return pointer to splitted node if split happens)
    node *insert(vector<double> &P) {
        if (type == 1) {
            // we are in leaf
            node *newNode = new node();
            for (int i = 0; i < D; i++) {
                newNode->Rectangle[0][i] = newNode->Rectangle[1][i] = P[i];
            }
            children[numberChildren] = newNode;
            numberChildren++;

            // redefine boundries
            if (numberChildren <= M) {
                if (numberChildren == 1) {
                    for (int i = 0; i < D; i++) {
                        Rectangle[0][i] = P[i];
                        Rectangle[1][i] = P[i];
                    }
                } else {
                    for (int i = 0; i < D; i++) {
                        Rectangle[0][i] = min(Rectangle[0][i], P[i]);
                        Rectangle[1][i] = max(Rectangle[1][i], P[i]);
                    }
                }
            }

            node *retval = NULL;
            if (numberChildren > M) {
                retval = split();
            }

            return retval;
        }

        // we are in non-Leaf node
        double minareaIncrease = (children[0]->areaIncrease)(P);
        int indexChildren = 0;
        for (int i = 1; i < numberChildren; i++) {
            double areaIncrease_ = (children[i]->areaIncrease)(P);
            if (areaIncrease_ < minareaIncrease) {
                minareaIncrease = areaIncrease_;
                indexChildren = i;
            }
        }
        node *Extra = children[indexChildren]->insert(P);
        // redefine boundries
        for (int i = 0; i < D; i++) {
            Rectangle[0][i] = min(Rectangle[0][i], children[indexChildren]->Rectangle[0][i]);
            Rectangle[1][i] = max(Rectangle[1][i], children[indexChildren]->Rectangle[1][i]);
        }
        if (Extra == NULL)
            return NULL;
        children[numberChildren] = Extra;
        numberChildren++;

        if (numberChildren <= M) {
            // redefine boundries
            for (int i = 0; i < D; i++) {
                Rectangle[0][i] = min(Rectangle[0][i], children[numberChildren - 1]->Rectangle[0][i]);
                Rectangle[1][i] = max(Rectangle[1][i], children[numberChildren - 1]->Rectangle[1][i]);
            }
            return NULL;
        }

        return split();
    }

    // main insert function takes care of new root node emergence condition while inserting
    node *insertMain(node *Tree, vector<double> &P) {
        node *Sibling = Tree->insert(P);
        if (Sibling == NULL) {
            return Tree;
        }
        node *newRoot = new node(2, 2);
        // making the MBR for new root
        for (int i = 0; i < D; i++) {
            newRoot->Rectangle[0][i] = min(Tree->Rectangle[0][i], Sibling->Rectangle[0][i]);
            newRoot->Rectangle[1][i] = max(Tree->Rectangle[1][i], Sibling->Rectangle[1][i]);
        }
        newRoot->children[0] = Tree;
        newRoot->children[1] = Sibling;

        return newRoot;
    }

    // returns distance between 2 points according to metric
    double distanceBetween(vector<double> &P1, vector<double> &P2) {
        if (metric == 2) {
            // Euclidian
            double Dist = 0;
            for (int i = 0; i < D; i++) {
                Dist += (P1[i] - P2[i]) * (P1[i] - P2[i]);
            }
            return sqrt(Dist);
        }
        if (metric == 1) {
            // Manhattan
            double Dist = 0;
            for (int i = 0; i < D; i++) {
                Dist += abs(P1[i] - P2[i]);
            }
            return Dist;
        }
        if (metric == 0) {
            // max norm
            double Dist = abs(P1[0] - P2[0]);
            for (int i = 1; i < D; i++) {
                Dist = max(Dist, abs(P1[i] - P2[i]));
            }
            return Dist;
        }
        return -1;
    }

    // returns min Distance (defined in Documentation) of MBR from point Q
    double minDist(vector<double> &Q) {
        // return min distance between Query point Q and the Rectangle of the node
        vector<double> Q_close(D); // the point in Rectangle closest to Q
        for (int i = 0; i < D; i++) {
            if (Q[i] < Rectangle[0][i]) {
                Q_close[i] = Rectangle[0][i];
            } else if (Q[i] > Rectangle[1][i]) {
                Q_close[i] = Rectangle[1][i];
            } else {
                Q_close[i] = Q[i];
            }
        }
        return distanceBetween(Q, Q_close);
    }

    // returns minmaxDistance (defined in documentation) of MBR from point Q
    double minmaxDist(vector<double> &Q) {
        // returns minmax Distance between the MBR and the point Q
        // minmax distance is the minimum distance from Q such that definitely a point of the MBR lies in that distance
        vector<double> Dist_min(D);
        vector<double> Dist_max(D);
        for (int i = 0; i < D; i++) {
            if (2 * Q[i] <= Rectangle[0][i]) {
                Dist_min[i] = abs(Q[i] - Rectangle[0][i]);
                Dist_max[i] = abs(Q[i] - Rectangle[1][i]);
            } else {
                Dist_min[i] = abs(Q[i] - Rectangle[1][i]);
                Dist_max[i] = abs(Q[i] - Rectangle[0][i]);
            }
        }
        if (metric == 0) {
            // max Norm
            // Basically we have to return 2nd maximum (from Dist_min and Dist_max combined.. proof in documentation)
            double max1 = Dist_max[0];
            double max2 = Dist_min[0];
            for (int i = 1; i < D; i++) {
                if (Dist_max[i] > max1) {
                    max2 = max(max1, Dist_min[i]);
                    max1 = Dist_max[i];
                } else if (Dist_max[i] > max2) {
                    max2 = Dist_max[i];
                }
            }
            return max2;
        }
        if (metric == 1) {
            // manhattan
            // We require minimum of (Dist_max[0] + Dist_max[1] ....Dist_max[D - 1]) - (Dist_max[i] - Dist_min[i])
            // So we require sum(Dist_max) - maximum(Dist_max[i] - Dist_min[i])

            double Delta = Dist_max[0] - Dist_min[0];
            double sum = Dist_max[0];
            for (int i = 1; i < D; i++) {
                sum += Dist_max[i];
                Delta = max(Delta, Dist_max[i] - Dist_min[i]);
            }
            return sum - Delta;
        }
        if (metric == 2) {
            // euclidian
            // We require minimum of (Dist_max[0]^2 + Dist_max[1]^2 ....Dist_max[D - 1]^2) - (Dist_max[i]^2 - Dist_min[i]^2)
            // So we require sum(Dist_max^2) - maximum(Dist_max[i]^2 - Dist_min[i]^2)
            double Delta = Dist_max[0] * Dist_max[0] - Dist_min[0] * Dist_min[0];
            double sum = Dist_max[0] * Dist_max[0];
            for (int i = 1; i < D; i++) {
                sum += Dist_max[i] * Dist_max[i];
                Delta = max(Delta, Dist_max[i] * Dist_max[i] - Dist_min[i] * Dist_min[i]);
            }
            return sqrt(sum - Delta);
        }
        return -1;
    }
    void Nearest_Neighbour(vector<double> &Q, double *Dist, vector<double> &Nearest) {
        if (type == 1) {
            // we are in leaf
            double minmaxDist_ = minmaxDist(Q);
            double minDist_ = minDist(Q);
            if (minDist_ > (*Dist)) {
                return;
            }
            for (int i = 0; i < numberChildren; i++) {
                double dist = distanceBetween(children[i]->Rectangle[0], Q);
                if (dist < (*Dist)) {
                    for (int j = 0; j < D; j++) {
                        Nearest[j] = (children[i]->Rectangle[0])[j];
                    }
                    (*Dist) = dist;
                }
            }
            return;
        }

        // we are in non - Leaf node
        // Algorithm (Find the children with minimum minmaxDist(call it cutoffDist) and prune the search for children with minDist > cutoffDist)
        double cutoffDist = (children[0]->minmaxDist)(Q);
        for (int i = 1; i < numberChildren; i++) {
            cutoffDist = min(cutoffDist, (children[i]->minmaxDist)(Q));
        }
        for (int i = 0; i < numberChildren; i++) {
            if ((children[i]->minDist)(Q) <= min(cutoffDist, *Dist)) {
                children[i]->Nearest_Neighbour(Q, Dist, Nearest);
            }
        }
        return;
    }
};
node *insertInTree(vector<vector<double>> &Points, node *Tree) {
    int n = (int)Points.size();
    node G; // to access the insertMain
    for (int i = 0; i < n; i++) {
        Tree = G.insertMain(Tree, Points[i]);
    }
    return Tree;
}
vector<double> generatePoint() {
    vector<double> Point(D);
    for (int i = 0; i < D; i++) {
        Point[i] = (double)rand() / RAND_MAX;
        Point[i] = -10 + Point[i] * 20;
    }
    return Point;
}
vector<vector<double>> generatePoints(int n) {
    vector<vector<double>> Points(n);
    for (int i = 0; i < n; i++) {
        Points[i] = generatePoint();
    }
    return Points;
}
int main() {
    int n = 50; // number of points to insert in tree
    metric = 2; // Setting metric to 2 for L_2 (Euclidian)
    vector<vector<double>> Points;
    Points = generatePoints(n);
    node *Tree = new node(0, 1);
    Tree = insertInTree(Points, Tree);

    int q = 10; // number of query points
    vector<vector<double>> queryPoints;
    queryPoints = generatePoints(q);
    vector<vector<double>> actualNearest(q);   // Actual nearest neighbour
    vector<double> actualDistance(q, DBL_MAX); // Distance of actual nearest neighbour

    node G; // for using the class function
    for (int i = 0; i < q; i++) {
        int nearestPointIndex = 0;
        actualDistance[i] = G.distanceBetween(queryPoints[i], Points[0]);
        for (int j = 1; j < n; j++) {
            double dist = G.distanceBetween(queryPoints[i], Points[j]);
            if (dist < actualDistance[i]) {
                actualDistance[i] = dist;
                nearestPointIndex = j;
            }
        }
        actualNearest[i] = Points[nearestPointIndex];
    }

    vector<vector<double>> NearestPoint(q, vector<double>(D, 0)); // to store nearest neighbour (obtained from Algo)
    vector<double> Distance(q, DBL_MAX);                          // to store distance of nearest neighbour (obtained from Algo)
    for (int i = 0; i < q; i++) {
        Tree->Nearest_Neighbour(queryPoints[i], &Distance[i], NearestPoint[i]);
    }

    cout << "POINTS : "
         << "\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < D; j++) {
            cout << Points[i][j] << " ";
        }
        cout << '\n';
    }

    cout << '\n';
    cout << "QUERY POINTS : "
         << "\n";
    for (int i = 0; i < q; i++) {
        for (int j = 0; j < D; j++) {
            cout << queryPoints[i][j] << " ";
        }
        cout << '\n';
    }

    cout << '\n';
    cout << "NEAREST NEIGHBOUR (As per Algo)" << '\n';
    for (int i = 0; i < q; i++) {
        cout << Distance[i] << ": ";
        for (int j = 0; j < D; j++) {
            cout << NearestPoint[i][j] << " ";
        }
        cout << '\n';
    }

    cout << '\n';
    cout << "ACTUAL NEIGHBOUR" << '\n';
    for (int i = 0; i < q; i++) {
        cout << actualDistance[i] << ": ";
        for (int j = 0; j < D; j++) {
            cout << actualNearest[i][j] << " ";
        }
        cout << '\n';
    }

    return 0;
}
