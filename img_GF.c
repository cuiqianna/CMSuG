#include<mex.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <time.h>
#include <matrix.h>

double RGB2rgb(double val) {
    if (val <= 0.04045) {
        return val / 12.92;
    }
    else {
        return pow((val + 0.055) / 1.055, 2.4);
    }
}

void Rgb2Lab(int* rvec, int* gvec, int* bvec, int length, double* luminosityvec,
    double* channel_Avec, double* channel_Bvec) {
    int index;
    int tempR, tempG, tempB;
    const double eps = 0.008856;
    const double kap = 903.3;
    const double XRef = 0.950456;
    const double YRef = 1.0;
    const double ZRef = 1.088754;
    const double coeX[3] = { 0.4124564,0.3575761,0.1804375 };
    const double coeY[3] = { 0.2126729,0.7151522,0.0721750 };
    const double coeZ[3] = { 0.0193339,0.1191920,0.9503041 };
    double X, Y, Z;
    double Red, Green, Blue;
    double red, green, blue;
    double xRef, yRef, zRef;
    double fx, fy, fz;
    double lValue, aValue, bValue;

    for (index = 0; index < length; ++index) {
        tempR = rvec[index];
        tempG = gvec[index];
        tempB = bvec[index];
        Red = tempR / 255.0;
        Green = tempG / 255.0;
        Blue = tempB / 255.0;

        red = RGB2rgb(Red);
        green = RGB2rgb(Green);
        blue = RGB2rgb(Blue);

        X = red * coeX[0] + green * coeX[1] + blue * coeX[2];
        Y = red * coeY[0] + green * coeY[1] + blue * coeY[2];
        Z = red * coeZ[0] + green * coeZ[1] + blue * coeZ[2];

        xRef = X / XRef;
        yRef = Y / YRef;
        zRef = Z / ZRef;

        if (xRef > eps) { fx = pow(xRef, 1.0 / 3.0); }
        else { fx = (kap * xRef + 16.0) / 116.0; }
        if (yRef > eps) { fy = pow(yRef, 1.0 / 3.0); }
        else { fy = (kap * yRef + 16.0) / 116.0; }
        if (zRef > eps) { fz = pow(zRef, 1.0 / 3.0); }
        else { fz = (kap * zRef + 16.0) / 116.0; }

        lValue = 116.0 * fy - 16.0;
        aValue = 500.0 * (fx - fy);
        bValue = 200.0 * (fy - fz);

        luminosityvec[index] = lValue;
        channel_Avec[index] = aValue;
        channel_Bvec[index] = bValue;


    }

}

void getInitialClusterCentroids(int GridStep, int width, int height,
    int* clusterlabels, int* numclusters) {
    int temp, index;
    int xsteps, ysteps;
    int error1, error2;
    double error1perstep, error2perstep;
    int xthreshold, ythreshold;
    int axisX, axisY;
    int xerror1, yerror2;
    int clusterx, clustery;

    xsteps = (int)(0.5 + (double)(width) / (double)(GridStep));
    ysteps = (int)(0.5 + (double)(height) / (double)(GridStep));

    error1 = width - GridStep * xsteps;
    if (error1 < 0) {
        xsteps--;
        error1 = width - GridStep * xsteps;
    }
    error2 = height - GridStep * ysteps;
    if (error2 < 0) {
        ysteps--;
        error2 = height - GridStep * ysteps;
    }
    error1perstep = (double)(error1) / (double)(xsteps);
    error2perstep = (double)(error2) / (double)(ysteps);

    xthreshold = GridStep / 2;
    ythreshold = GridStep / 2;

    index = 0;
    for (axisY = 0; axisY < ysteps; ++axisY) {
        yerror2 = (int)axisY * error2perstep;
        for (axisX = 0; axisX < xsteps; ++axisX) {
            xerror1 = (int)axisX * error1perstep;
            clusterx = (axisX * GridStep + xthreshold + xerror1);
            clustery = (axisY * GridStep + ythreshold + yerror2);
            temp = clustery * width + clusterx;
            clusterlabels[index] = temp;
            ++index;
        }
    }
    *numclusters = index;

}

void FuzzySLIC(double* luminosityvec, double* channel_Avec, double* channel_Bvec, double* clusterl, double* clustera, double* clusterb, double* clusterx, double* clustery, int width, int height, int numclusters, int* oldlabels, int GridStep)
{
    int x1, y1, x2, y2;
    double l, a, b, mid1, mid2, mid3, mid4, mid5, mid6,a1,a2,a3,a4;
    double dist;
    double distxy;
    int itr;
    int n;
    int x, y;
    int i, i2;
    int k, pos;
    double possi;
    int sz = width * height;
    const int numk = numclusters;
    int offset = GridStep;
    double alpha = 0.9;


    double* Hinv = (double*)mxMalloc(sizeof(double) * numk);
    double* U1 = (double*)mxMalloc(sizeof(double) * numk);
    double* Ux = (double*)mxMalloc(sizeof(double) * numk);
    double* Uy = (double*)mxMalloc(sizeof(double) * numk);
    double* Ul = (double*)mxMalloc(sizeof(double) * numk);
    double* Ua = (double*)mxMalloc(sizeof(double) * numk);
    double* Ub = (double*)mxMalloc(sizeof(double) * numk);

    double(*U)[3] = (double(*)[3])mxMalloc((sizeof(double) * sz) * 3);
    double(*H)[3] = (double(*)[3])mxMalloc(sizeof(double) * sz * 3);
    int(*memb)[3] = (int(*)[3])mxMalloc(sizeof(int) * sz * 3);
    int* count = (int*)mxMalloc(sizeof(int) * sz);
    double(*Lab_mat_c)[3] = (double(*)[3])mxMalloc(sizeof(double) * sz * 3);
    int* x11 = (int(*))mxMalloc(sizeof(int) * sz);
    int* y11 = (int(*))mxMalloc(sizeof(int) * sz);

    //////////////////////
    double invxywt = 1.0 / (GridStep * GridStep);
    double* maxlab = (double*)mxMalloc(sizeof(double) * numk);//variable M, the compactness factor or color distnce normalization factor
    double* distlab = (double*)mxMalloc(sizeof(double) * sz);

    for (i = 0; i < sz; i++)
    {
        distlab[i] = DBL_MAX;
    }
    for (n = 0; n < numk; n++)
    {
        maxlab[n] = 10.0 * 10.0;//initialize with some reasonable compactness value
    }
    ///////////////////////
    for (itr = 0; itr < 10; itr++)
    {
        for (k = 0; k < sz; k++)
        {
            count[k] = 0;
            Lab_mat_c[k][0] = 0;
            Lab_mat_c[k][1] = 0;
            Lab_mat_c[k][2] = 0;
            H[k][0] = 0;
            H[k][1] = 0;
            H[k][2] = 0;
            U[k][0] = 0;
            U[k][1] = 0;
            U[k][2] = 0;
        }

        for (k = 0; k < numk; k++)
        {
            Hinv[k] = 0;
            U1[k] = 0;
            Ux[k] = 0;
            Uy[k] = 0;
            Ul[k] = 0;
            Ua[k] = 0;
            Ub[k] = 0;

        }
        for (n = 0; n < numk; n++)
        {
            x1 = clusterx[n] - offset; if (x1 < 0) x1 = 0;
            y1 = clustery[n] - offset; if (y1 < 0) y1 = 0;
            x2 = clusterx[n] + offset; if (x2 > width)  x2 = width;
            y2 = clustery[n] + offset; if (y2 > height) y2 = height;

            for (y = y1; y < y2; y++)
            {
                for (x = x1; x < x2; x++)
                {
                    i = y * width + x;
                    x11[i] = x + 1;
                    y11[i] = y + 1;
                    l = luminosityvec[i];
                    a = channel_Avec[i];
                    b = channel_Bvec[i];
                    ////////////////
                    distlab[i] = (l - clusterl[n]) * (l - clusterl[n]) +
                        (a - clustera[n]) * (a - clustera[n]) +
                        (b - clusterb[n]) * (b - clusterb[n]);
                    ////////////////
                    distxy = (x - clusterx[n]) * (x - clusterx[n]) + (y - clustery[n]) * (y - clustery[n]);

                    dist = distlab[i] / maxlab[n] + distxy * invxywt;

                    dist = sqrt(dist);

                    if (dist < 0.000000001)
                    {
                        if (count[i] < 3)
                        {
                            if (count[i] == 0) {
                                memb[i][count[i]] = n;
                                Lab_mat_c[i][count[i]] = 0.000000001;
                                count[i] += 1;
                            }
                            if (count[i] == 1) {
                                memb[i][count[i] - 1] = n;
                                Lab_mat_c[i][count[i] - 1] = 0.000000001;
                            }
                            if (count[i] == 2) {

                                memb[i][count[i] - 2] = n;
                                Lab_mat_c[i][count[i] - 2] = 0.000000001;
                                count[i] -= 1;

                            }


                        }
                        else {

                            memb[i][0] = n;
                            Lab_mat_c[i][0] = 0.000000001;
                            memb[i][1] = n;
                            Lab_mat_c[i][1] = 0.000000001;
                            memb[i][2] = n;
                            Lab_mat_c[i][2] = 0.000000001;

                        }



                    }
                    else {
                        if (count[i] < 3)
                        {
                            memb[i][count[i]] = n;
                            Lab_mat_c[i][count[i]] = dist;
                            count[i] += 1;

                        }
                        else {
                            if (Lab_mat_c[i][0] < Lab_mat_c[i][1]) {

                                pos = 1;

                                if (Lab_mat_c[i][1] < Lab_mat_c[i][2]) {

                                    pos = 2;

                                }
                            }
                            else {

                                pos = 0;

                                if (Lab_mat_c[i][0] < Lab_mat_c[i][2]) {

                                    pos = 2;

                                }
                            }

                            if (dist < Lab_mat_c[i][pos]) {

                                memb[i][pos] = n;
                                Lab_mat_c[i][pos] = dist;

                            }
                            else {


                            }



                        }

                    }
                }
            }
        }


        for (i = 0; i < sz; i++) {

            if (count[i] == 1) {

                U[i][0] = 1.0;

            }
            if (count[i] == 2) {

                if (Lab_mat_c[i][0] < Lab_mat_c[i][1]) {
                    a1 = Lab_mat_c[i][0] * Lab_mat_c[i][0] * alpha;
                    mid1 = Lab_mat_c[i][0] * Lab_mat_c[i][0] - a1;
                    mid2 = Lab_mat_c[i][1]* Lab_mat_c[i][1] - a1;
                    U[i][0] = 1.0/(1.0 + mid1 / mid2);
                    U[i][1] = 1.0 / (1.0 + mid2 / mid1);
                }

                else {
                    a2 = Lab_mat_c[i][1] * Lab_mat_c[i][1] * alpha;             
                    mid1 = Lab_mat_c[i][0] * Lab_mat_c[i][0] - a2;
                    mid2 = Lab_mat_c[i][1] * Lab_mat_c[i][1] - a2;
                    U[i][0] = 1.0 / (1.0 + mid1 / mid2);
                    U[i][1] = 1.0 / (1.0 + mid2 / mid1);
                }


            }
            if (count[i] == 3) {

                if (Lab_mat_c[i][0] < Lab_mat_c[i][1]) {
                    if (Lab_mat_c[i][0] < Lab_mat_c[i][2]) {
                        a1 = Lab_mat_c[i][0] * Lab_mat_c[i][0] * alpha;
                        mid1 = Lab_mat_c[i][0] * Lab_mat_c[i][0] - a1;
                        mid2 = Lab_mat_c[i][1] * Lab_mat_c[i][1] - a1;
                        mid3 = Lab_mat_c[i][2] * Lab_mat_c[i][2] - a1;
                        U[i][0] = 1.0 / (1.0 + mid1 / mid2 + mid1/mid3);
                        U[i][1] = 1.0 / (1.0 + mid2 / mid1 + mid2/mid3);
                        U[i][2] = 1.0 / (1.0 + mid3 / mid1 + mid3 / mid2);
                    }
                    else {
                        a2 = Lab_mat_c[i][2] * Lab_mat_c[i][2] * alpha;
                        mid1 = Lab_mat_c[i][0] * Lab_mat_c[i][0] - a2;
                        mid2 = Lab_mat_c[i][1] * Lab_mat_c[i][1] - a2;
                        mid3 = Lab_mat_c[i][2] * Lab_mat_c[i][2] - a2;
                        U[i][0] = 1.0 / (1.0 + mid1 / mid2 + mid1 / mid3);
                        U[i][1] = 1.0 / (1.0 + mid2 / mid1 + mid2 / mid3);
                        U[i][2] = 1.0 / (1.0 + mid3 / mid1 + mid3 / mid2);
                    }
                }
                else {
                    if (Lab_mat_c[i][1] < Lab_mat_c[i][2]) {
                        a3 = Lab_mat_c[i][1] * Lab_mat_c[i][1] * alpha;
                        mid1 = Lab_mat_c[i][0] * Lab_mat_c[i][0] - a3;
                        mid2 = Lab_mat_c[i][1] * Lab_mat_c[i][1] - a3;
                        mid3 = Lab_mat_c[i][2] * Lab_mat_c[i][2] - a3;
                        U[i][0] = 1.0 / (1.0 + mid1 / mid2 + mid1 / mid3);
                        U[i][1] = 1.0 / (1.0 + mid2 / mid1 + mid2 / mid3);
                        U[i][2] = 1.0 / (1.0 + mid3 / mid1 + mid3 / mid2);
                    }
                    else {
                        a4 = Lab_mat_c[i][2] * Lab_mat_c[i][2] * alpha;
                        mid1 = Lab_mat_c[i][0] * Lab_mat_c[i][0] - a4;
                        mid2 = Lab_mat_c[i][1] * Lab_mat_c[i][1] - a4;
                        mid3 = Lab_mat_c[i][2] * Lab_mat_c[i][2] - a4;
                        U[i][0] = 1.0 / (1.0 + mid1 / mid2 + mid1 / mid3);
                        U[i][1] = 1.0 / (1.0 + mid2 / mid1 + mid2 / mid3);
                        U[i][2] = 1.0 / (1.0 + mid3 / mid1 + mid3 / mid2);
                    }
                }
                
            }


        }

        for (i = 0; i < sz; i++) {

            if (count[i] == 1) {

                if (((x11[i] - 1) >= 1) && ((y11[i] - 1) >= 1)) {
                    i2 = (y11[i] - 2) * width + x11[i] - 2;
                    if (count[i2] == 3) {
                        if (memb[i2][0] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][0];

                        }
                        if (memb[i2][1] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][1];

                        }
                        if (memb[i2][2] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][2];

                        }
                    }
                    if (count[i2] == 2) {
                        if (memb[i2][0] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][0];

                        }
                        if (memb[i2][1] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][1];

                        }
                    }
                    if (count[i2] == 1) {
                        if (memb[i2][0] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][0];

                        }
                    }
                }
                if (((x11[i] + 1) <= width) && ((y11[i] + 1) <= height)) {
                    i2 = (y11[i]) * width + x11[i];
                    if (count[i2] == 3) {
                        if (memb[i2][0] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][0];

                        }
                        if (memb[i2][1] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][1];

                        }
                        if (memb[i2][2] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][2];

                        }
                    }
                    if (count[i2] == 2) {
                        if (memb[i2][0] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][0];

                        }
                        if (memb[i2][1] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][1];

                        }
                    }
                    if (count[i2] == 1) {
                        if (memb[i2][0] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][0];

                        }
                    }
                }
                if (((x11[i] + 1) <= width) && ((y11[i] - 1) >= 1)) {
                    i2 = (y11[i] - 2) * width + x11[i];
                    if (count[i2] == 3) {
                        if (memb[i2][0] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][0];

                        }
                        if (memb[i2][1] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][1];

                        }
                        if (memb[i2][2] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][2];

                        }
                    }
                    if (count[i2] == 2) {
                        if (memb[i2][0] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][0];

                        }
                        if (memb[i2][1] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][1];

                        }
                    }
                    if (count[i2] == 1) {
                        if (memb[i2][0] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][0];

                        }
                    }
                }
                if (((x11[i] - 1) >= 1) && ((y11[i] + 1) <= height)) {
                    i2 = (y11[i]) * width + x11[i] - 2;
                    if (count[i2] == 3) {
                        if (memb[i2][0] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][0];

                        }
                        if (memb[i2][1] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][1];

                        }
                        if (memb[i2][2] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][2];

                        }
                    }
                    if (count[i2] == 2) {
                        if (memb[i2][0] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][0];

                        }
                        if (memb[i2][1] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][1];

                        }
                    }
                    if (count[i2] == 1) {
                        if (memb[i2][0] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][0];

                        }
                    }
                }
                if ((x11[i] - 1) >= 1) {
                    i2 = (y11[i] - 1) * width + x11[i] - 2;
                    if (count[i2] == 3) {
                        if (memb[i2][0] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][0];

                        }
                        if (memb[i2][1] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][1];

                        }
                        if (memb[i2][2] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][2];

                        }
                    }
                    if (count[i2] == 2) {
                        if (memb[i2][0] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][0];

                        }
                        if (memb[i2][1] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][1];

                        }
                    }
                    if (count[i2] == 1) {
                        if (memb[i2][0] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][0];

                        }
                    }
                }
                if ((y11[i] - 1) >= 1) {
                    i2 = (y11[i] - 2) * width + x11[i] - 1;
                    if (count[i2] == 3) {
                        if (memb[i2][0] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][0];

                        }
                        if (memb[i2][1] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][1];

                        }
                        if (memb[i2][2] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][2];

                        }
                    }
                    if (count[i2] == 2) {
                        if (memb[i2][0] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][0];

                        }
                        if (memb[i2][1] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][1];

                        }
                    }
                    if (count[i2] == 1) {
                        if (memb[i2][0] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][0];

                        }
                    }
                }
                if ((x11[i] + 1) <= width) {
                    i2 = (y11[i] - 1) * width + x11[i];
                    if (count[i2] == 3) {
                        if (memb[i2][0] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][0];

                        }
                        if (memb[i2][1] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][1];

                        }
                        if (memb[i2][2] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][2];

                        }
                    }
                    if (count[i2] == 2) {
                        if (memb[i2][0] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][0];

                        }
                        if (memb[i2][1] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][1];

                        }
                    }
                    if (count[i2] == 1) {
                        if (memb[i2][0] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][0];

                        }
                    }
                }
                if ((y11[i] + 1) <= height) {
                    i2 = (y11[i]) * width + x11[i] - 1;
                    if (count[i2] == 3) {
                        if (memb[i2][0] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][0];

                        }
                        if (memb[i2][1] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][1];

                        }
                        if (memb[i2][2] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][2];

                        }
                    }
                    if (count[i2] == 2) {
                        if (memb[i2][0] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][0];

                        }
                        if (memb[i2][1] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][1];

                        }
                    }
                    if (count[i2] == 1) {
                        if (memb[i2][0] == memb[i][0]) {
                            H[i][0] = H[i][0] + U[i2][0];

                        }
                    }
                }

                Hinv[memb[i][0]] = Hinv[memb[i][0]] + H[i][0] * H[i][0];
            }
            if (count[i] == 2) {

                for (k = 0; k < 2; k++) {
                    if (((x11[i] - 1) >= 1) && ((y11[i] - 1) >= 1)) {
                        i2 = (y11[i] - 2) * width + x11[i] - 2;
                        if (count[i2] == 3) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                            if (memb[i2][2] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][2];

                            }
                        }
                        if (count[i2] == 2) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                        }
                        if (count[i2] == 1) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                        }
                    }
                    if (((x11[i] + 1) <= width) && ((y11[i] + 1) <= height)) {
                        i2 = (y11[i]) * width + x11[i];
                        if (count[i2] == 3) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                            if (memb[i2][2] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][2];

                            }
                        }
                        if (count[i2] == 2) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                        }
                        if (count[i2] == 1) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                        }
                    }
                    if (((x11[i] + 1) <= width) && ((y11[i] - 1) >= 1)) {
                        i2 = (y11[i] - 2) * width + x11[i];
                        if (count[i2] == 3) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                            if (memb[i2][2] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][2];

                            }
                        }
                        if (count[i2] == 2) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                        }
                        if (count[i2] == 1) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                        }
                    }
                    if (((x11[i] - 1) >= 1) && ((y11[i] + 1) <= height)) {
                        i2 = (y11[i]) * width + x11[i] - 2;
                        if (count[i2] == 3) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                            if (memb[i2][2] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][2];

                            }
                        }
                        if (count[i2] == 2) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                        }
                        if (count[i2] == 1) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                        }
                    }
                    if ((x11[i] - 1) >= 1) {
                        i2 = (y11[i] - 1) * width + x11[i] - 2;
                        if (count[i2] == 3) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                            if (memb[i2][2] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][2];

                            }
                        }
                        if (count[i2] == 2) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                        }
                        if (count[i2] == 1) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                        }
                    }
                    if ((y11[i] - 1) >= 1) {
                        i2 = (y11[i] - 2) * width + x11[i] - 1;
                        if (count[i2] == 3) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                            if (memb[i2][2] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][2];

                            }
                        }
                        if (count[i2] == 2) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                        }
                        if (count[i2] == 1) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                        }
                    }
                    if ((x11[i] + 1) <= width) {
                        i2 = (y11[i] - 1) * width + x11[i];
                        if (count[i2] == 3) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                            if (memb[i2][2] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][2];

                            }
                        }
                        if (count[i2] == 2) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                        }
                        if (count[i2] == 1) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                        }
                    }
                    if ((y11[i] + 1) <= height) {
                        i2 = (y11[i]) * width + x11[i] - 1;
                        if (count[i2] == 3) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                            if (memb[i2][2] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][2];

                            }
                        }
                        if (count[i2] == 2) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                        }
                        if (count[i2] == 1) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                        }
                    }





                }
                Hinv[memb[i][0]] = Hinv[memb[i][0]] + H[i][0] * H[i][0];
                Hinv[memb[i][1]] = Hinv[memb[i][1]] + H[i][1] * H[i][1];
            }
            if (count[i] == 3) {
                for (k = 0; k < 3; k++) {
                    if (((x11[i] - 1) >= 1) && ((y11[i] - 1) >= 1)) {
                        i2 = (y11[i] - 2) * width + x11[i] - 2;
                        if (count[i2] == 3) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                            if (memb[i2][2] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][2];

                            }
                        }
                        if (count[i2] == 2) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                        }
                        if (count[i2] == 1) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                        }
                    }
                    if (((x11[i] + 1) <= width) && ((y11[i] + 1) <= height)) {
                        i2 = (y11[i]) * width + x11[i];
                        if (count[i2] == 3) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                            if (memb[i2][2] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][2];

                            }
                        }
                        if (count[i2] == 2) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                        }
                        if (count[i2] == 1) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                        }
                    }
                    if (((x11[i] + 1) <= width) && ((y11[i] - 1) >= 1)) {
                        i2 = (y11[i] - 2) * width + x11[i];
                        if (count[i2] == 3) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                            if (memb[i2][2] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][2];

                            }
                        }
                        if (count[i2] == 2) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                        }
                        if (count[i2] == 1) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                        }
                    }
                    if (((x11[i] - 1) >= 1) && ((y11[i] + 1) <= height)) {
                        i2 = (y11[i]) * width + x11[i] - 2;
                        if (count[i2] == 3) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                            if (memb[i2][2] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][2];

                            }
                        }
                        if (count[i2] == 2) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                        }
                        if (count[i2] == 1) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                        }
                    }
                    if ((x11[i] - 1) >= 1) {
                        i2 = (y11[i] - 1) * width + x11[i] - 2;
                        if (count[i2] == 3) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                            if (memb[i2][2] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][2];

                            }
                        }
                        if (count[i2] == 2) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                        }
                        if (count[i2] == 1) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                        }
                    }
                    if ((y11[i] - 1) >= 1) {
                        i2 = (y11[i] - 2) * width + x11[i] - 1;
                        if (count[i2] == 3) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                            if (memb[i2][2] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][2];

                            }
                        }
                        if (count[i2] == 2) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                        }
                        if (count[i2] == 1) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                        }
                    }
                    if ((x11[i] + 1) <= width) {
                        i2 = (y11[i] - 1) * width + x11[i];
                        if (count[i2] == 3) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                            if (memb[i2][2] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][2];

                            }
                        }
                        if (count[i2] == 2) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                        }
                        if (count[i2] == 1) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                        }
                    }
                    if ((y11[i] + 1) <= height) {
                        i2 = (y11[i]) * width + x11[i] - 1;
                        if (count[i2] == 3) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                            if (memb[i2][2] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][2];

                            }
                        }
                        if (count[i2] == 2) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                            if (memb[i2][1] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][1];

                            }
                        }
                        if (count[i2] == 1) {
                            if (memb[i2][0] == memb[i][k]) {
                                H[i][k] = H[i][k] + U[i2][0];

                            }
                        }
                    }

                }

                Hinv[memb[i][0]] = Hinv[memb[i][0]] + H[i][0] * H[i][0];
                Hinv[memb[i][1]] = Hinv[memb[i][1]] + H[i][1] * H[i][1];
                Hinv[memb[i][2]] = Hinv[memb[i][2]] + H[i][2] * H[i][2];
            }


        }

        for (i = 0; i < sz; i++) {

            if (count[i] == 1) {
                U[i][0] = H[i][0] * H[i][0] / Hinv[memb[i][0]];
                U1[memb[i][0]] = U1[memb[i][0]] + U[i][0] * U[i][0];
                Ux[memb[i][0]] = Ux[memb[i][0]] + U[i][0] * U[i][0] * (x11[i] - 1);
                Uy[memb[i][0]] = Uy[memb[i][0]] + U[i][0] * U[i][0] * (y11[i] - 1);
                Ul[memb[i][0]] = Ul[memb[i][0]] + U[i][0] * U[i][0] * luminosityvec[i];
                Ua[memb[i][0]] = Ua[memb[i][0]] + U[i][0] * U[i][0] * channel_Avec[i];
                Ub[memb[i][0]] = Ub[memb[i][0]] + U[i][0] * U[i][0] * channel_Bvec[i];
                oldlabels[i] = memb[i][0];

            }
            if (count[i] == 2) {
                U[i][0] = H[i][0] * H[i][0] / Hinv[memb[i][0]];
                U[i][1] = H[i][1] * H[i][1] / Hinv[memb[i][1]];

                U1[memb[i][0]] = U1[memb[i][0]] + U[i][0] * U[i][0];
                Ux[memb[i][0]] = Ux[memb[i][0]] + U[i][0] * U[i][0] * (x11[i] - 1);
                Uy[memb[i][0]] = Uy[memb[i][0]] + U[i][0] * U[i][0] * (y11[i] - 1);
                Ul[memb[i][0]] = Ul[memb[i][0]] + U[i][0] * U[i][0] * luminosityvec[i];
                Ua[memb[i][0]] = Ua[memb[i][0]] + U[i][0] * U[i][0] * channel_Avec[i];
                Ub[memb[i][0]] = Ub[memb[i][0]] + U[i][0] * U[i][0] * channel_Bvec[i];

                U1[memb[i][1]] = U1[memb[i][1]] + U[i][1] * U[i][1];
                Ux[memb[i][1]] = Ux[memb[i][1]] + U[i][1] * U[i][1] * (x11[i] - 1);
                Uy[memb[i][1]] = Uy[memb[i][1]] + U[i][1] * U[i][1] * (y11[i] - 1);
                Ul[memb[i][1]] = Ul[memb[i][1]] + U[i][1] * U[i][1] * luminosityvec[i];
                Ua[memb[i][1]] = Ua[memb[i][1]] + U[i][1] * U[i][1] * channel_Avec[i];
                Ub[memb[i][1]] = Ub[memb[i][1]] + U[i][1] * U[i][1] * channel_Bvec[i];



                if (U[i][0] > U[i][1]) {

                    oldlabels[i] = memb[i][0];

                }
                else {

                    oldlabels[i] = memb[i][1];

                }


            }
            if (count[i] == 3) {

                U[i][0] = H[i][0] * H[i][0] / Hinv[memb[i][0]];
                U[i][1] = H[i][1] * H[i][1] / Hinv[memb[i][1]];
                U[i][2] = H[i][2] * H[i][2] / Hinv[memb[i][2]];

                U1[memb[i][0]] = U1[memb[i][0]] + U[i][0] * U[i][0];
                Ux[memb[i][0]] = Ux[memb[i][0]] + U[i][0] * U[i][0] * (x11[i] - 1);
                Uy[memb[i][0]] = Uy[memb[i][0]] + U[i][0] * U[i][0] * (y11[i] - 1);
                Ul[memb[i][0]] = Ul[memb[i][0]] + U[i][0] * U[i][0] * luminosityvec[i];
                Ua[memb[i][0]] = Ua[memb[i][0]] + U[i][0] * U[i][0] * channel_Avec[i];
                Ub[memb[i][0]] = Ub[memb[i][0]] + U[i][0] * U[i][0] * channel_Bvec[i];

                U1[memb[i][1]] = U1[memb[i][1]] + U[i][1] * U[i][1];
                Ux[memb[i][1]] = Ux[memb[i][1]] + U[i][1] * U[i][1] * (x11[i] - 1);
                Uy[memb[i][1]] = Uy[memb[i][1]] + U[i][1] * U[i][1] * (y11[i] - 1);
                Ul[memb[i][1]] = Ul[memb[i][1]] + U[i][1] * U[i][1] * luminosityvec[i];
                Ua[memb[i][1]] = Ua[memb[i][1]] + U[i][1] * U[i][1] * channel_Avec[i];
                Ub[memb[i][1]] = Ub[memb[i][1]] + U[i][1] * U[i][1] * channel_Bvec[i];

                U1[memb[i][2]] = U1[memb[i][2]] + U[i][2] * U[i][2];
                Ux[memb[i][2]] = Ux[memb[i][2]] + U[i][2] * U[i][2] * (x11[i] - 1);
                Uy[memb[i][2]] = Uy[memb[i][2]] + U[i][2] * U[i][2] * (y11[i] - 1);
                Ul[memb[i][2]] = Ul[memb[i][2]] + U[i][2] * U[i][2] * luminosityvec[i];
                Ua[memb[i][2]] = Ua[memb[i][2]] + U[i][2] * U[i][2] * channel_Avec[i];
                Ub[memb[i][2]] = Ub[memb[i][2]] + U[i][2] * U[i][2] * channel_Bvec[i];

                if (U[i][0] > U[i][1]) {

                    oldlabels[i] = memb[i][0];

                    if (U[i][2] > U[i][0]) {
                        oldlabels[i] = memb[i][2];
                    }
                }
                else {

                    oldlabels[i] = memb[i][1];
                    if (U[i][2] > U[i][1]) {
                        oldlabels[i] = memb[i][2];
                    }
                }


            }


        }
        /////////////////////
        if (0 == itr)
        {
            for (n = 0; n < numk; n++) maxlab[n] = 1.0;
        }
        for (i = 0; i < sz; i++)
        {
            if (maxlab[oldlabels[i]] < distlab[i]) maxlab[oldlabels[i]] = distlab[i];
        }
        //////////////////////

        for (k = 0; k < numk; k++)
        {

            clusterl[k] = Ul[k] / U1[k];
            clustera[k] = Ua[k] / U1[k];
            clusterb[k] = Ub[k] / U1[k];
            clusterx[k] = Ux[k] / U1[k];
            clustery[k] = Uy[k] / U1[k];

        }


    }

    mxFree(Hinv);
    mxFree(H);
    mxFree(U);
    mxFree(U1);
    mxFree(Ul);
    mxFree(Ua);
    mxFree(Ub);
    mxFree(Ux);
    mxFree(Uy);
    mxFree(Lab_mat_c);
    mxFree(count);
    mxFree(y11);
    mxFree(x11);
    mxFree(memb);
    mxFree(maxlab);
    mxFree(distlab);

}

void ConnectivityPostProcessing(int* preLabels, int width, int height,
    int numSuperpixels, int* postLabels, int* numOfPostLabels) {
    int i, j;
    int inner, outer, total;
    int x, y;
    int index;
    int oldIndex, newIndex, adjacentLabel;
    int numOfLabels;
    const int neighborx4[4] = { -1,0,1,0 };
    const int neighbory4[4] = { 0,-1,0,1 };
    const int size = width * height;
    const int superpixelSize = size / numSuperpixels;
    int* xvec = (int*)mxMalloc(sizeof(int) * superpixelSize * 10);
    int* yvec = (int*)mxMalloc(sizeof(int) * superpixelSize * 10);

    for (i = 0; i < size; ++i) {
        postLabels[i] = -1;
    }
    oldIndex = 0;
    adjacentLabel = 0;
    numOfLabels = 0;
    for (i = 0; i < height; ++i) {
        for (j = 0; j < width; ++j) {
            if (postLabels[oldIndex] < 0) {
                postLabels[oldIndex] = numOfLabels;
                xvec[0] = j;
                yvec[0] = i;

                for (inner = 0; inner < 4; inner++) {
                    x = xvec[0] + neighborx4[inner];
                    y = yvec[0] + neighbory4[inner];
                    if ((x >= 0 && x < width) && (y >= 0 && y < height)) {
                        newIndex = y * width + x;
                        if (postLabels[newIndex] >= 0) {
                            adjacentLabel = postLabels[newIndex];
                        }
                    }
                }
                total = 1;
                for (outer = 0; outer < total; ++outer) {
                    for (inner = 0; inner < 4; ++inner) {
                        x = xvec[outer] + neighborx4[inner];
                        y = yvec[outer] + neighbory4[inner];
                        if ((x >= 0 && x < width) && (y >= 0 && y < height)) {
                            newIndex = y * width + x;
                            if (postLabels[newIndex] < 0 && preLabels[oldIndex] == preLabels[newIndex]) {
                                xvec[total] = x;
                                yvec[total] = y;
                                postLabels[newIndex] = numOfLabels;
                                total++;
                            }
                        }
                    }
                }
                if (total <= superpixelSize >> 2) {
                    for (outer = 0; outer < total; ++outer) {
                        index = yvec[outer] * width + xvec[outer];
                        postLabels[index] = adjacentLabel;
                    }
                    numOfLabels--;
                }
                numOfLabels++;
            }
            oldIndex++;
        }
    }
    *numOfPostLabels = numOfLabels;

    mxFree(xvec);
    mxFree(yvec);
}

void mexFunction(int nlhs, mxArray* plhs[],
    int nrhs, const mxArray* prhs[]) {
    int width;
    int height;
    int size;
    int indexOfVec, indexOfImg, i;
    int x, y;
    int* rvec, * gvec, * bvec;
    double* luminosityvec, * channel_Avec, * channel_Bvec;
    double* clusterx, * clustery;
    double* clusterl, * clustera, * clusterb;
    const mwSize* dimensions;
    int* clusterlabels;
    int numclusters;
    int* ouputLabels;
    int numOfPostLabels;
    unsigned char* imgbytes;
    int numParameters;
    int GridStep;
    int numOfSuperpixels;
    int* oldlabels;
    int* newlabels;
    mwSize numOfDims;

    
    imgbytes = (unsigned char*)mxGetPr(prhs[0]);
    numParameters = (int)mxGetNumberOfElements(prhs[0]);
    numOfDims = mxGetNumberOfDimensions(prhs[0]);
    dimensions = mxGetDimensions(prhs[0]);
    numOfSuperpixels = (int)mxGetScalar(prhs[1]);
    
    height = dimensions[0];
    width = dimensions[1];
    size = width * height;

 

    rvec = (int*)mxMalloc(sizeof(int) * size);
    gvec = (int*)mxMalloc(sizeof(int) * size);
    bvec = (int*)mxMalloc(sizeof(int) * size);

    luminosityvec = (double*)mxMalloc(sizeof(double) * size);
    channel_Avec = (double*)mxMalloc(sizeof(double) * size);
    channel_Bvec = (double*)mxMalloc(sizeof(double) * size);
    oldlabels = (int*)mxMalloc(sizeof(int) * size);
    newlabels = (int*)mxMalloc(sizeof(int) * size);
    clusterlabels = (int*)mxMalloc(sizeof(int) * size);
 

    for (x = 0, indexOfImg = 0; x < width; ++x) {
        for (y = 0; y < height; ++y) {
                indexOfVec = y * width + x;
                rvec[indexOfVec] = imgbytes[indexOfImg];
                gvec[indexOfVec] = imgbytes[indexOfImg + size];
                bvec[indexOfVec] = imgbytes[indexOfImg + size + size];
                ++indexOfImg;
        }
    }
    Rgb2Lab(rvec, gvec, bvec, size, luminosityvec, channel_Avec, channel_Bvec);
    
    GridStep = sqrt((double)(size) / (double)(numOfSuperpixels)) + 0.5;
    getInitialClusterCentroids(GridStep, width, height, clusterlabels, &numclusters);

    clusterx = (double*)mxMalloc(sizeof(double) * numclusters);
    clustery = (double*)mxMalloc(sizeof(double) * numclusters);
    clusterl = (double*)mxMalloc(sizeof(double) * numclusters);
    clustera = (double*)mxMalloc(sizeof(double) * numclusters);
    clusterb = (double*)mxMalloc(sizeof(double) * numclusters);

    for (i = 0; i < numclusters; ++i) {
        clusterx[i] = clusterlabels[i] % width;
        clustery[i] = clusterlabels[i] / width;
        clusterl[i] = luminosityvec[clusterlabels[i]];
        clustera[i] = channel_Avec[clusterlabels[i]];
        clusterb[i] = channel_Bvec[clusterlabels[i]];
    }

    FuzzySLIC(luminosityvec, channel_Avec, channel_Bvec, clusterl, clustera, clusterb, clusterx,
        clustery, width, height, numclusters, oldlabels, GridStep);
    ConnectivityPostProcessing(oldlabels, width, height, numOfSuperpixels,
        newlabels, &numOfPostLabels);


    plhs[0] = mxCreateNumericMatrix(height, width, mxINT32_CLASS, mxREAL);
    ouputLabels = (int*)mxGetPr(plhs[0]);
    for (x = 0, indexOfImg = 0; x < width; ++x) {
        for (y = 0; y < height; ++y) {
            indexOfVec = y * width + x;
            ouputLabels[indexOfImg] = newlabels[indexOfVec];
            indexOfImg++;
        }
    }



    mxFree(rvec);
    mxFree(gvec);
    mxFree(bvec);
    mxFree(oldlabels);
    mxFree(newlabels);
    mxFree(luminosityvec);
    mxFree(channel_Avec);
    mxFree(channel_Bvec);
    mxFree(clusterx);
    mxFree(clustery);
    mxFree(clusterl);
    mxFree(clustera);
    mxFree(clusterb);
}
