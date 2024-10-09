// Imagine++ project
// Project:  Panorama
// Author:   Pascal Monasse
// Date:     2013/10/08

#include <Imagine/Graphics.h>
#include <Imagine/Images.h>
#include <Imagine/LinAlg.h>
#include <vector>
#include <sstream>
using namespace Imagine;
using namespace std;

// Record clicks in two images, until right button click
void getClicks(Window w1, Window w2,
               vector<IntPoint2>& pts1, vector<IntPoint2>& pts2) {
    IntPoint2 currentPoint;
    Window currentWindow;
    int subWin = NULL; // we dont need
    int button;
    int pointNumber;

    while(true)
    {
        button = anyGetMouse(currentPoint, currentWindow, subWin);  // return 1=left 2=middle 3=right
        if(button==3) { // right click
            return;
        }
        if(button==2){ // middle click
            cout << "I said : \"left\" click !" << endl;
        }
        else if(currentWindow == w1) { // left click on image 1
            pts1.push_back(currentPoint); // add currentPoint in the pts1 vector
            pointNumber = pts1.size();
            // draw a point
            setActiveWindow(w1);
            drawString(currentPoint, to_string(pointNumber), RED);
            drawCircle(currentPoint,4,Color(0,0,255),3,false);
            cout << "Point number " << pointNumber << " added on first image" << endl;
        } else { // left click on image 2
            pts2.push_back(currentPoint); // add currentPoint in the pts2 vector
            pointNumber = pts2.size();
            // draw a point
            setActiveWindow(w2);
            drawString(currentPoint, to_string(pointNumber), RED);
            drawCircle(currentPoint,4,Color(0,0,255),3,false);
            cout << "Point number " << pointNumber << " added on second image" << endl;
        }
    }
}

// Return homography compatible with point matches
Matrix<float> getHomography(const vector<IntPoint2>& pts1,
                            const vector<IntPoint2>& pts2) {
    size_t n = min(pts1.size(), pts2.size());
    if(n<4) {
        cout << "Not enough correspondences: " << n << endl;
        return Matrix<float>::Identity(3);
    }
    Matrix<float> A(2*n,8);
    Vector<float> B(2*n);


    for(int i=0; i<n; i++){
        IntPoint2 pt1 = pts1[i];
        IntPoint2 pt2 = pts2[i];

        // compute the homography
        A(2*i,0) = pt1.x();
        A(2*i,1) = pt1.y();
        A(2*i,2) = 1;
        A(2*i,3) = 0;
        A(2*i,4) = 0;
        A(2*i,5) = 0;
        A(2*i,6) = - (pt2.x()*pt1.x());
        A(2*i,7) = - (pt2.x()*pt1.y());
        B[2*i] = pt2.x();

        A(2*i+1,0) = 0;
        A(2*i+1,1) = 0;
        A(2*i+1,2) = 0;
        A(2*i+1,3) = pt1.x();
        A(2*i+1,4) = pt1.y();
        A(2*i+1,5) = 1;
        A(2*i+1,6) = - (pt2.y()*pt1.x());
        A(2*i+1,7) = - (pt2.y()*pt1.y());
        B[2*i+1] = pt2.y();
    }

    B = linSolve(A, B);
    Matrix<float> H(3, 3);
    H(0,0)=B[0]; H(0,1)=B[1]; H(0,2)=B[2];
    H(1,0)=B[3]; H(1,1)=B[4]; H(1,2)=B[5];
    H(2,0)=B[6]; H(2,1)=B[7]; H(2,2)=1;

    // Sanity check
    cout << "Sanity check"<< endl;
    for(size_t i=0; i<n; i++) {
        float v1[]={(float)pts1[i].x(), (float)pts1[i].y(), 1.0f};
        float v2[]={(float)pts2[i].x(), (float)pts2[i].y(), 1.0f};
        Vector<float> x1(v1,3);
        Vector<float> x2(v2,3);
        x1 = H*x1;
        cout << x1[1]*x2[2]-x1[2]*x2[1] << ' '
             << x1[2]*x2[0]-x1[0]*x2[2] << ' '
             << x1[0]*x2[1]-x1[1]*x2[0] << endl;
    }
    return H;
}

// Grow rectangle of corners (x0,y0) and (x1,y1) to include (x,y)
void growTo(float& x0, float& y0, float& x1, float& y1, float x, float y) {
    if(x<x0) x0=x;
    if(x>x1) x1=x;
    if(y<y0) y0=y;
    if(y>y1) y1=y;
}

// Panorama construction
void panorama(const Image<Color,2>& I1, const Image<Color,2>& I2,
              Matrix<float> H) {
    Vector<float> v(3);
    // Initializes image at size of second image
    float x0=0, y0=0, x1=I2.width(), y1=I2.height();

    v[0]=0; v[1]=0; v[2]=1;
    v=H*v; v/=v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    v[0]=I1.width(); v[1]=0; v[2]=1;
    v=H*v; v/=v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    v[0]=I1.width(); v[1]=I1.height(); v[2]=1;
    v=H*v; v/=v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    v[0]=0; v[1]=I1.height(); v[2]=1;
    v=H*v; v/=v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    cout << "x0 x1 y0 y1=" << x0 << ' ' << x1 << ' ' << y0 << ' ' << y1<<endl;

    Image<Color> I(int(x1-x0), int(y1-y0));
    setActiveWindow( openWindow(I.width(), I.height()) );
    I.fill(WHITE);

    Matrix<float> invH = inverse(H);
    Vector<float> coord(3), prec_coord;

    for (int x = 0; x < I.width(); x++) {
        for (int y = 0; y < I.height(); y++) {
            coord[0] = x+x0; coord[1] = y+y0; coord[2] = 1; // upper left corner
            if (coord[0] > 0 && coord[1] > 0 && coord[0] < I2.width() && coord[1] < I2.height()) {
                I(x,y) = I2(coord[0], coord[1]); // if the coord are inside I2, the color of the pixel is from I2
            }
            prec_coord = coord;
            coord = invH*coord;
            coord /= coord[2]; // normalization
            if (coord[0] > 0 && coord[1] > 0 && coord[0] < I1.width() && coord[1] < I1.height()) {
                if(I(x,y) != WHITE){ // if I(x,y) != WHITE, it means this pixel is in I1 and I2, so we take the average of colors on overlapping area
                    I(x, y)[0] = (I2(prec_coord[0], coord[1])[0] + I1.interpolate(coord[0], coord[1]).r()) / 2;
                    I(x, y)[1] = (I2(prec_coord[0], coord[1])[1] + I1.interpolate(coord[0], coord[1]).g()) / 2;
                    I(x, y)[2] = (I2(prec_coord[0], coord[1])[2] + I1.interpolate(coord[0], coord[1]).b()) / 2;
                }
                else{
                    I(x,y) = I1.interpolate(coord[0], coord[1]); // if the coord are inside I1, the color of the pixel is from I1
                }
            }
        }
    }
    display(I,0,0);
}

// Main function
int main(int argc, char* argv[]) {
    const char* s1 = argc>1? argv[1]: srcPath("image0006.jpg"); //pointer on image 1
    const char* s2 = argc>2? argv[2]: srcPath("image0007.jpg"); //pointer on image 2

    // Load and display images
    Image<Color> I1, I2;
    if( ! load(I1, s1) ||
        ! load(I2, s2) ) {
        cerr<< "Unable to load the images" << endl;
        return 1;
    }
    Window w1 = openWindow(I1.width(), I1.height(), s1);
    display(I1,0,0);
    Window w2 = openWindow(I2.width(), I2.height(), s2);
    setActiveWindow(w2);
    display(I2,0,0);
    cout << "Please select at least 5 matching points on the two images by left clicking on them. " << endl;
    // Get user's clicks in images
    vector<IntPoint2> pts1, pts2;
    getClicks(w1, w2, pts1, pts2);

    vector<IntPoint2>::const_iterator it;
    cout << "pts1="<<endl;
    for(it=pts1.begin(); it != pts1.end(); it++)
        cout << *it << endl;
    cout << "pts2="<<endl;
    for(it=pts2.begin(); it != pts2.end(); it++)
        cout << *it << endl;

    // Compute homography
    Matrix<float> H = getHomography(pts1, pts2);
    cout << "H=" << H/H(2,2);

    // Apply homography
    panorama(I1, I2, H);

    endGraphics();
    return 0;
}
