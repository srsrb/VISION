// Imagine++ project
// Project:  Fundamental
// Author:   Pascal Monasse

#include "./Imagine/Features.h"
#include <Imagine/Graphics.h>
#include <Imagine/LinAlg.h>
#include <vector>
#include <cstdlib>
#include <ctime>
using namespace Imagine;
using namespace std;

static const float BETA = 0.01f; // Probability of failure

struct Match {
    float x1, y1, x2, y2;
};

// Display SIFT points and fill vector of point correspondences
void algoSIFT(Image<Color,2> I1, Image<Color,2> I2,
              vector<Match>& matches) {
    // Find interest points
    SIFTDetector D;
    D.setFirstOctave(-1);
    Array<SIFTDetector::Feature> feats1 = D.run(I1);
    drawFeatures(feats1, Coords<2>(0,0));
    cout << "Im1: " << feats1.size() << flush;
    Array<SIFTDetector::Feature> feats2 = D.run(I2);
    drawFeatures(feats2, Coords<2>(I1.width(),0));
    cout << " Im2: " << feats2.size() << flush;

    const double MAX_DISTANCE = 100.0*100.0;
    for(size_t i=0; i < feats1.size(); i++) {
        SIFTDetector::Feature f1=feats1[i];
        for(size_t j=0; j < feats2.size(); j++) {
            double d = squaredDist(f1.desc, feats2[j].desc);
            if(d < MAX_DISTANCE) {
                Match m;
                m.x1 = f1.pos.x();
                m.y1 = f1.pos.y();
                m.x2 = feats2[j].pos.x();
                m.y2 = feats2[j].pos.y();
                matches.push_back(m);
            }
        }
    }
}

// RANSAC algorithm to compute F from point matches (8-point algorithm)
// Parameter matches is filtered to keep only inliers as output.
FMatrix<float,3,3> computeF(vector<Match>& matches) {
    const float distMax = 1.5f; // Pixel error for inlier/outlier discrimination
    int Niter=100000; // Adjusted dynamically
    FMatrix<float,3,3> bestF;
    vector<Match> bestInliers;
    // --------------- TODO ------------
    // DO NOT FORGET NORMALIZATION OF POINTS

    FMatrix<float,3,3> N(0.0);  // Normalization mat
    N(0,0) = 0.001;
    N(1,1) = 0.001;
    N(2,2) = 1;

    for(int iter = 0; iter < Niter; iter++) {
        vector<Match> currentInliers;
        FMatrix<float,9,9> A(0.0);
        FMatrix<float,3,3> E(0.0);
        int randIndex = rand(); // to select random matches

        // Create random A mat
        for(int i = 0; i < 8; i++) {
            FVector<float,3> pl(0), pr(0);
            randIndex = (randIndex+1) % matches.size();
            pl[0] = matches.at(randIndex).x1;
            pl[1] = matches.at(randIndex).y1;
            pl[2] = 1;
            pr[0] = matches.at(randIndex).x2;
            pr[1] = matches.at(randIndex).y2;
            pr[2] = 1;
            pl = N * pl;
            pr = N * pr;
            A(i,0) = pl[0] * pr[0];
            A(i,1) = pl[0] * pr[1];
            A(i,2) = pl[0];
            A(i,3) = pl[1] * pr[0];
            A(i,4) = pl[1] * pr[1];
            A(i,5) = pl[1];
            A(i,6) = pr[0];
            A(i,7) = pr[1];
            A(i,8) = 1;
        }
        for(int i = 0; i < 9; i++) {
            A(8,i) = 0;
        }

        // Fundamental mat estimation using 8-point algorithm
        FMatrix<float,3,3> F(0.0);
        FMatrix<float,9,9> U(0.0), Vt(0.0);
        FVector<float,9> S(0.0);
        FMatrix<float,3,3> F2(0.0);
        FMatrix<float,3,3> U2(0.0), Vt2(0.0);
        FVector<float,3> S2(0.0);
        A = transpose(A) * A;
        svd(A, U, S, Vt);
        for(int i = 0; i < 9; i++) {
            F2(i/3,i%3) = Vt(8,i);
        }
        svd(F2, U2, S2, Vt2);
        S2[2] = 0;
        F = U2 * (Diagonal(S2) * Vt2);
        F = transpose(N) * (F * N);
        E = F / norm(F);

        // Get the Inliers
        for(int i = 0; i < matches.size(); i++) {
            FVector<float,3> pl(0), pr(0);
            pl[0] = matches.at(i).x1;
            pl[1] = matches.at(i).y1;
            pl[2] = 1;
            pr[0] = matches.at(i).x2;
            pr[1] = matches.at(i).y2;
            pr[2] = 1;
            FVector<float,3> line = E * pr;
            line = line / (sqrt(line[0] * line[0] + line[1] * line[1]));
            float distance = pl * line;
            distance *= distance;
            if(distance < distMax) {
                currentInliers.push_back(matches.at(i));
            }
        }

        // more Inliers then before ? -> update
        int currentInliersSize = currentInliers.size();
        int InliersSize = bestInliers.size();
        if(currentInliersSize > InliersSize) {
            bestInliers.clear();
            for(int i = 0; i < currentInliersSize; i++) {
                bestInliers.push_back(currentInliers.back());
                currentInliers.pop_back();
            }

            int Niter_tmp;

            float l= log(1- pow( (float)currentInliersSize/(float)matches.size(), 8));

            if(l < 0){
                Niter_tmp = ceil(log(BETA)/l);

                if (Niter_tmp < 1 || Niter_tmp > 1000000) {
                    Niter_tmp = 1000;
                }
            }
            else{
                Niter_tmp = -1;
            }
            if (Niter_tmp != -1){
                Niter = Niter_tmp;
            }
        }
    }

    // Updating matches with inliers only
    matches.clear();
    for(size_t i=0; i<bestInliers.size(); i++)
        matches.push_back(bestInliers[i]);

    // Refine the fundamental mat
    int maxNbInliers = bestInliers.size();
    Matrix<float> rA(maxNbInliers,9);
    Matrix<float> Ura(maxNbInliers,maxNbInliers);
    Vector<float> Sra(maxNbInliers);
    Matrix<float> Vtra(9,9);
    Vector<float> fVector(9);
    FMatrix<float,3,3> reformedF;
    FMatrix<float,3,3> U3(0.0), Vt3(0.0);
    FVector<float,3> S3(0.0);
    for(int i = 0; i < maxNbInliers; i++) {
        FVector<float,3> pl(0), pr(0);
        Match m = bestInliers.back();
        bestInliers.pop_back();
        pl[0] = m.x1;
        pl[1] = m.y1;
        pl[2] = 1;
        pr[0] = m.x2;
        pr[1] = m.y2;
        pr[2] = 1;
        pl = N * pl;
        pr = N * pr;
        rA(i,0) = pl[0] * pr[0];
        rA(i,1) = pl[0] * pr[1];
        rA(i,2) = pl[0];
        rA(i,3) = pl[1] * pr[0];
        rA(i,4) = pl[1] * pr[1];
        rA(i,5) = pl[1];
        rA(i,6) = pr[0];
        rA(i,7) = pr[1];
        rA(i,8) = 1;
    }
    svd(rA, Ura, Sra, Vtra);
    fVector = Vtra.getRow(8);
    for(int i = 0; i < 9; i++) {
        reformedF(i/3, i%3) = fVector[i];
    }
    svd(reformedF, U3, S3, Vt3);
    S3[2] = 0;
    bestF = U3 * (Diagonal(S3) * Vt3);
    bestF = transpose(N) * bestF * N;
    bestF = bestF / norm(bestF);

    return bestF;
}

// Expects clicks in one image and show corresponding line in other image.
// Stop at right-click.
void displayEpipolar(Image<Color> I1, Image<Color> I2, const FMatrix<float,3,3>& F) {
    while(true) {
        int x, y;
        if(getMouse(x, y) == 3) break;
        // --------------- TODO ------------

        fillCircle(x, y, 2, RED, false);
        if(x < I1.width()) {
            FVector<float,3> pr(0.0), pl(0.0);
            FMatrix<float,3,3> Ft = transpose(F);
            pl[0] = x;
            pl[1] = y;
            pl[2] = 1;
            pr = Ft * pl;
            FVector<float,3> epl1(0), epl2(0);
            epl1[0] = 0;
            epl1[1] = (-pr[2] - pr[0] * epl1[0]) / pr[1];
            epl2[0] = I1.width();
            epl2[1] = (-pr[2] - pr[0] * epl2[0]) / pr[1];
            epl1[0] = epl1[0] + I1.width();
            epl2[0] = epl2[0] + I1.width();
            drawLine(epl1[0], epl1[1], epl2[0], epl2[1], RED, 2);
        }
        else {
            FVector<float,3> pl(0), pr(0);
            pr[0] = x - I1.width();
            pr[1] = y;
            pr[2] = 1;
            pl = F * pr;
            FVector<float,3> epl1(0), epl2(0);
            epl1[0] = 0;
            epl1[1] = (-pl[2] - pl[0] * epl1[0]) / pl[1];
            epl2[0] = I1.width();
            epl2[1] = (-pl[2] - pl[0] * epl2[0]) / pl[1];
            drawLine(epl1[0], epl1[1], epl2[0], epl2[1], RED, 2);
        }
    }
}

int main(int argc, char* argv[])
{
    srand((unsigned int)time(0));

    const char* s1 = argc>1? argv[1]: srcPath("im1.jpg");
    const char* s2 = argc>2? argv[2]: srcPath("im2.jpg");

    // Load and display images
    Image<Color,2> I1, I2;
    if( ! load(I1, s1) ||
        ! load(I2, s2) ) {
        cerr<< "Unable to load images" << endl;
        return 1;
    }
    int w = I1.width();
    openWindow(2*w, I1.height());
    display(I1,0,0);
    display(I2,w,0);

    vector<Match> matches;
    algoSIFT(I1, I2, matches);
    const int n = (int)matches.size();
    cout << " matches: " << n << endl;
    drawString(100,20,std::to_string(n)+ " matches",RED);
    click();

    FMatrix<float,3,3> F = computeF(matches);
    cout << "F="<< endl << F;

    // Redisplay with matches
    display(I1,0,0);
    display(I2,w,0);
    for(size_t i=0; i<matches.size(); i++) {
        Color c(rand()%256,rand()%256,rand()%256);
        fillCircle(matches[i].x1+0, matches[i].y1, 2, c);
        fillCircle(matches[i].x2+w, matches[i].y2, 2, c);
    }
    drawString(100, 20, to_string(matches.size())+"/"+to_string(n)+" inliers", RED);
    click();

    // Redisplay without SIFT points
    display(I1,0,0);
    display(I2,w,0);
    displayEpipolar(I1, I2, F);

    endGraphics();
    return 0;
}
