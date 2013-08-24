/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package autolight;

/**
 *
 * @author Andrey
 */
public class Retinex {
static void HS(int[] rgb, int width, int height, int[] hs1, int[] hs2, int[] hs3, int[] r, int[] s){
    int R = 0;
    int G = 0;
    int B = 0;
    int R1 = 255; int R2 = 0;int G1 = 255; int G2 = 0; int B1 = 255; int B2 = 0;
    int S_1 = 0; int S_2 = 0; int S_3 = 0;
    int size = width*height;
    int c = 0;
    for (int j = 1; j < height - 1; j++){
        for (int i = 1; i < width - 1; i++){
            c = i + j*width; 
            R = (rgb[c] >> 16)&0xff;
            G = (rgb[c] >> 8)&0xff;
            B = rgb[c]&0xff;
            if (R1 > R) R1 = R;
            if (R2 < R) R2 = R;
            if (G1 > G) G1 = G;
            if (G2 < G) G2 = G;
            if (B1 > B) B1 = B;
            if (B2 < B) B2 = B;
            S_1 += R;
            S_2 += G;
            S_3 += B;
            hs1[R] ++;
            hs2[G] ++;
            hs3[B] ++;
                }
    }
    s[0] = S_1/size;
    s[1] = S_2/size;
    s[2] = S_3/size;
    r[0] = R1;
    r[1] = R2;
    r[2] = G1;
    r[3] = G2;
    r[4] = B1;
    r[5] = B2;
}

    
    static void correctionl(int[] rgb, int width, int height, int[] hs1, int[] hs2, int[] hs3, int r1, int r2, int g1, int g2, int b1, int b2){
    double b = 0;

    double b_1 = 0;
    double b_2 = 0;
    double b_3 = 0;

    double a_1 = 0;
    double a_2 = 0;
    double a_3 = 0;


    if ((r2 - r1) !=0 ) b_1 = 255 /(r2 - r1);
    if ((g2 - g1) !=0 ) b_2 = 255 /(g2 - g1);
    if ((b2 - b1) !=0 ) b_3 = 255 /(b2 - b1);
    a_1 = -b_1*r1;
    a_2 = -b_2*g1;
    a_3 = -b_3*b1;

    int A;int R;int G;int B;
    double L = 0;
    int c;
    for (int j = 1; j < height; j++){
        for (int i = 1; i < width; i++){
            c = i + j*width;
            R = (rgb[c] >> 16)&0xff;
            G = (rgb[c] >> 8)&0xff;
            B = rgb[c]&0xff;
            if ((r2 - r1) != 0) L = a_1 + b_1*R;
            R = ((int) L);
            if ((g2 - g1) != 0) L = a_2 + b_2*G;;
            G = ((int) L);
            if ((b2 - b1) != 0) L = a_3 + b_3*B;
            B = ((int) L);
            A = (int) (255 & 0xff);
            if (R > 255) R = 255;
            if (R < 0) R = 0;
            if (G > 255) G = 255;
            if (G < 0) G = 0;
            if (B > 255) B = 255;
            if (B < 0) B = 0;
            rgb[c] = 0xff000000 | (A << 24) | (R << 16) | (G << 8) | B;
          }
    }
    }
    
    static void normalization(int[] rgb, int width, int height, int[] hs1, int[] hs2, int[] hs3, int r1, int r2, int g1, int g2, int b1, int b2){
	int R; int G; int B; int A;
	double[]  WS1 = new double[256];
        double[]  WS2 = new double[256];
        double[]  WS3 = new double[256];
        double S1 = 0.0;
	double R_1 = 0, R_2 = 0, R_3 = 0;
	int size = width*height;
	for (int i = 0; i < 255; i++){
            R_1 += (double)hs1[i]/(double)size;
            R_2 += (double)hs2[i]/(double)size;
            R_3 += (double)hs3[i]/(double)size;
            WS1[i] = R_1;
            WS2[i] = R_2;
            WS3[i] = R_3;
	}
	int c = 0;
	double L;
	for (int j = 1; j < height; j++){
            for (int i = 1; i < width; i++){
                c = i + j*width;
                R = (rgb[c] >> 16)&0xff;
                G = (rgb[c] >> 8)&0xff;
                B = rgb[c]&0xff;
                L = (r2 - r1)*WS1[R] + r1;
                R = ((int) L);
                L = (g2 - g1)*WS2[G] + g1;
                G = ((int) L);
                L = (b2 - b1)*WS3[B] + b1;
                B = ((int) L);
                A = (int) (255 & 0xff);
                if (R > 255) R = 255;
                if (R < 0) R = 0;
                if (G > 255) G = 255;
                if (G < 0) G = 0;
                if (B > 255) B = 255;
                if (B < 0) B = 0;
                rgb[c] = 0xff000000 | (A << 24) | (R << 16) | (G << 8) | B;
              }
	}
    
}
    
    static int[] runProcess(int[] rgb, int width, int height){
	int[] HSR = new int[256];
	int[] HSG = new int[256];
	int[] HSB = new int[256];
	int[] r = new int[6];
	int[] s = new int[3];
	int R1 = 0;int R2 = 0;int G1 = 0;int G2 = 0;int B1 = 0;int B2 = 0;
	int s1 = 0, s2 = 0, s3 = 0;
	HS(rgb, width, height, HSR, HSG, HSB, r, s);
	R1 = r[0];R2 = r[1];G1 = r[2];G2 = r[3];B1 = r[4];B2 = r[5];
	correctionl(rgb, width, height, HSR, HSG, HSB, R1, R2, G1, G2, B1, B2);
	HSR = new int[256];
	HSG = new int[256];
	HSB = new int[256];
	r = new int[6];
	s = new int[3];
	HS(rgb, width, height, HSR, HSG, HSB, r, s);
	R1 = r[0];R2 = r[1];G1 = r[2];G2 = r[3];B1 = r[4];B2 = r[5];
	normalization(rgb, width, height, HSR, HSG, HSB, R1, R2, G1, G2, B1, B2);
	return rgb;
}
    
    static double[] getMatrix(double sigma){
	int N = (int) (2*sigma + 1);
	double mat[] = new double[N*N];
	int k = (N - 1)/2;
	double NN = N*N;
	double S = (double) (1/(2*Math.PI*sigma*sigma));
	int c = 0;
	for (int i = -k; i <= k; i++){
		for (int j = -k; j <=k; j++){
			mat[c] = (Math.exp(-(i*i + j*j)/(2*sigma*sigma)))/(NN);
			c++;
		}
	}
	return mat;
}
    
static double[] getMatrix1(double sigma, int size){
    int N = size;
    double mat[] = new double[N];
    int k = (N - 1)/2;
    int c = 0;
    for (int i = -k; i <= k; i++){
        mat[c] = Math.exp(-i*i/(2*sigma*sigma))/(double)N;
        c++;
    }
    return mat;
}
    
static int[] arr_t = null;
static int[] arr_t1 = null;

static int[] gauss5(int[] rgb, int width, int height, double sigma, int N){
    int[] arr = rgb.clone();
    //System.arraycopy(rgb, 0, arr, 0, rgb.length);
    //double sigma = 31;
    //int N = (int) (2*sigma + 1);
    //int N = 31;
    int k1 = (N - 1)/2;
    if (arr_t == null) arr_t = new int[(width + 2*k1)*(height + 2*k1)];
    if (arr_t1 == null) arr_t1 = new int[(width + 2*k1)*(height + 2*k1)];

    double[] mat = getMatrix1(sigma, N);
    
    int R = 0;
    int G = 0;
    int B = 0;
    int A = (255 & 0xff);
    int c = 0;
    int k = 0;
    int k_1 = 0;
    int k_2 = 0;
    int a1 = 0;
    int a2 = 0;
    int c1;
    int c2;
    int kk = 0;
    for (int j = 0; j < height; j++){
            for (int i = 0; i < width; i++){
                    c = i + j*width;
                    k = (i + k1) + (j + k1)*width;
                    arr_t[k] = rgb[c];
            }
    }
    for (int j = 0; j < k1; j++){
            for (int i = 0; i < k1; i++){
                    c1 = i + j*width;                            
                    c2 = (width - i - 1) + (height - j - 1)*width;
                    k_1 = i + j*(width + 2*k1);
                    k_2 = width + 2*k1 - i - 1 + (height + 2*k1 - j - 1)*(width + 2*k1);
                    arr_t[k_1] = rgb[c1];
                    arr_t[k_2] = rgb[c2];
            }
    }
    System.arraycopy(arr_t, 0, arr_t1, 0, (width + 2*k1)*(height + 2*k1));
    int width1 = width + 2*k1;
    int height1 = height + 2*k1;

    for (int j = k1; j < height1 - k1; j++){
        for (int i = k1; i < width1 - k1; i++){
            c = i + j*width1;
            double S1 = 0;
            double S2 = 0;
            double S3 = 0;
            for (int i1 = -k1; i1 <= k1; i1++){
                c1 = (i1 + k1);
                c = (i + i1) + j*width;
                R = (arr_t[c] >> 16) & 0xff;
                G = (arr_t[c] >> 8) & 0xff;
                B = arr_t[c] & 0xff;
                S1 += (double)R*mat[c1];
                S2 += (double)G*mat[c1];
                S3 += (double)B*mat[c1];
            }

            c = i + j*width1;
            R = (int)S1;
            G = (int)S2;
            B = (int)S3;
            arr_t1[c] = 0xff000000 | (A << 24) | (R << 16) | (G << 8) | B;
        }
    }
    for (int j = k1; j < height1 - k1; j++){
        for (int i = k1; i < width1 - k1; i++){
            c = i + j*width1;
            double S1 = 0;
            double S2 = 0;
            double S3 = 0;
            for (int i1 = -k1; i1 <= k1; i1++){
                c1 = (i1 + k1);
                c = i + (j + i1)*width1;
                R = (arr_t1[c] >> 16) & 0xff;
                G = (arr_t1[c] >> 8) & 0xff;
                B = arr_t1[c] & 0xff;
                S1 += (double)R*mat[c1];
                S2 += (double)G*mat[c1];
                S3 += (double)B*mat[c1];
            }
            c = (i - k1) + (j - k1)*width;
            R = (int)S1;
            G = (int)S2;
            B = (int)S3;
            rgb[c] = 0xff000000 | (A << 24) | (R << 16) | (G << 8) | B;
        }
}
return rgb;
	
}


static int[] retinex1(int[] rgb, int width, int height){
    int[] arr = null;
    //int[] arr1 = new int[rgb.length];
    System.gc();
    arr = rgb.clone();
    int[] rgb1 = gauss5(rgb,  width, height, 1000, 51);
    int R, G, B, A;
    int RG, GG, BG, AG;
    A = (255 & 0xff);
    int c;
    for (int j = 0; j < height; j++){
        for (int i = 0; i < width; i++){
            c = i + j*width;
            R = (arr[c] >> 16) & 0xff;
            G = (arr[c] >> 8) & 0xff;
            B = arr[c] & 0xff;
            RG = (rgb1[c] >> 16) & 0xff;
            GG = (rgb1[c] >> 8) & 0xff;
            BG = rgb1[c] & 0xff;
            R = (int) ((255.0*(Math.log10(R) - Math.log10(RG))) + 127.0);
            G = (int) ((255.0*(Math.log10(G) - Math.log10(GG))) + 127.0);
            B = (int) ((255.0*(Math.log10(B) - Math.log10(BG))) + 127.0);
            if (R > 255) R = 255;
            if (R < 0) R = 0;
            if (G > 255) G = 255;
            if (G < 0) G = 0;
            if (B > 255) B = 255;
            if (B < 0) B = 0;
            rgb[c] = 0xff000000 | (A << 24) | (R << 16) | (G << 8) | B;
        }
    }
    return rgb;
}

private static int [] arr = null;  
static int[] gauss(int[] rgb, int width, int height, double sigma, int N){
	int k1 = (N - 1)/2;
	double[] mat = getMatrix1(sigma, N);
	int R = 0;
	int G = 0;
	int B = 0;
	int A = (255 & 0xff);
	int c = 0;
	int c1;
	int kk = 0;
	if  (arr == null)  arr = rgb.clone();
	else if (arr.length != rgb.length){
		 arr = rgb.clone();
	}
	else {
		for (int j = 0; j < height; j++){
			for (int i = 0; i < width; i++){
				c = i + j*width;
				arr[c] = rgb[c];
			}
		}
	}
	for (int j = 0; j < height; j++){
            for (int i = 0; i < width; i++){
                double S1 = 0;
                double S2 = 0;
                double S3 = 0;
                for (int i1 = -k1; i1 <= k1; i1++){
                    c1 = (i1 + k1);
                    if ((i + i1) < 0) kk = k1 + i1;
                    else if ((i + i1) >= width) kk = i1 - k1;
                    else kk = i1;		
                    c = (i + kk) + j*width;
                    R = (rgb[c] >> 16) & 0xff;
                    G = (rgb[c] >> 8) & 0xff;
                    B = rgb[c] & 0xff;
                    S1 += (double)R*mat[c1];
                    S2 += (double)G*mat[c1];
                    S3 += (double)B*mat[c1];
                }
                c = i + j*width;
                R = (int)S1;
                G = (int)S2;
                B = (int)S3;
                arr[c] = 0xff000000 | (A << 24) | (R << 16) | (G << 8) | B;
            }
	}
	for (int j = 0; j < height; j++){
		for (int i = 0; i < width; i++){
			double S1 = 0;
			double S2 = 0;
			double S3 = 0;
			for (int i1 = -k1; i1 <= k1; i1++){
					c1 = (i1 + k1);
					if ((j + i1) < 0) kk = k1 + i1;
					else if ((j + i1) >= height) kk = i1 - k1;
					else kk = i1;
					c = i + (j + kk)*width;
					R = (arr[c] >> 16) & 0xff;
					G = (arr[c] >> 8) & 0xff;
					B = arr[c] & 0xff;
					S1 += (double)R*mat[c1];
					S2 += (double)G*mat[c1];
					S3 += (double)B*mat[c1];
			}
			c = i + j*width;
			R = (int)S1;
			G = (int)S2;
			B = (int)S3;
			rgb[c] = 0xff000000 | (A << 24) | (R << 16) | (G << 8) | B;
		}
	}
	return rgb;
	
}

static int[] retinex_ms(int[] rgb, int width, int height){
	int[] arr = null;
	double[] ls = {100, 500, 5000};
	double[] ws = {0.1, 0.3, 0.6};
	float arrR[] = new float[rgb.length];
	float arrG[] = new float[rgb.length];
	float arrB[] = new float[rgb.length];
	int R, G, B;
	float R1, G1, B1;
	int RG, GG, BG, AG;
	double a = 0;
	int A = (255 & 0xff);
	int c;
	int N = ls.length;
	if (arr == null) arr = rgb.clone();
	for (int j = 0; j < height; j++){
		for (int i = 0; i < width; i++){
			c = i + j*width;
			arrR[c] = 0;
			arrG[c] = 0;
			arrB[c] = 0;
		}
	}
	for (int i1 = 0; i1 < N; i1++){
		rgb = gauss(rgb, width, height, ls[i1], 25);
		for (int j = 0; j < height; j++){
			for (int i = 0; i < width; i++){
				c = i + j*width;
				R = (arr[c] >> 16) & 0xff;
				G = (arr[c] >> 8) & 0xff;
				B = arr[c] & 0xff;
				RG = (rgb[c] >> 16) & 0xff;
				GG = (rgb[c] >> 8) & 0xff;
				BG = rgb[c] & 0xff;
				//a = Math.log10((double)3*R/(R + G + B));
				R1 = (float) (ws[i1]*(Math.log10(R) - Math.log10(RG)));
				//a = Math.log10((double)3*G/(R + G + B));
				G1 = (float) (ws[i1]*(Math.log10(G) - Math.log10(GG)));
				//a = Math.log10((double)3*B/(R + G + B));
				B1 = (float) (ws[i1]*(Math.log10(B) - Math.log10(BG)));
				arrR[c] += R1;
				arrG[c] += R1;
				arrB[c] += R1;
			}
		}
		for (int j = 0; j < height; j++){
			for (int i = 0; i < width; i++){
				c = i + j*width;
				rgb[c] = arr[c];
			}
		}
	}
	for (int j = 0; j < height; j++){
		for (int i = 0; i < width; i++){
			c = i + j*width;
			R = (int) (255.0*arrR[c] + 127);
			G = (int) (255.0*arrG[c] + 127);
			B = (int) (255.0*arrB[c] + 127);
			RG = (arr[c] >> 16) & 0xff;
			GG = (arr[c] >> 8) & 0xff;
			BG = arr[c] & 0xff;
			//a = 46*Math.log10((double)125*RG/(RG + GG + BG));
			//R *= a;
			//a = 46*Math.log10((double)125*GG/(RG + GG + BG));
			//G *= a;
			//a = 46*Math.log10((double)125*BG/(RG + GG + BG));
			//B *= a;
			if (R > 255) R = 255;
			if (R < 0) R = 0;
			if (G > 255) G = 255;
			if (G < 0) G = 0;
			if (B > 255) B = 255;
			if (B < 0) B = 0;
			rgb[c] = 0xff000000 | (A << 24) | (R << 16) | (G << 8) | B;
		}
	}
	return rgb;
}

static int[] retinex_ms1(int[] rgb, int width, int height){
    final int hx = 4;
    final int hy = 4;
    if ((width > 100) && (height > 100)){
        int w = width/hx;
        int h = height/hy;
        int a = 0;
        int b = 0;
        int w1 = w;
        int h1 = h;
        int[] rgb1 = new int[w*h + 5]; 
        for (int j = 0; j < hy; j++){
                a = 0;
                w1 = w;
                for (int i = 0; i < hx; i++){
                        int k = 0;
                        for (int j1 = b; j1 < h1; j1++){
                                for (int i1 = a; i1 < w1; i1++){
                                        int c = i1 + j1*width;
                                        rgb1[k] = rgb[c];
                                        k++;
                                }
                        }
                        k = 0;
                        for (int j1 = b; j1 < h1; j1++){
                                for (int i1 = a; i1 < w1; i1++){
                                        int c = i1 + j1*width;
                                        rgb[c] = rgb1[k];
                                        k++;
                                }
                        }

                        a += w;
                        w1 += w;

                        if (i == hx - 2){
                                w1 = width;
                                rgb1 = new int[(width - a)*h + 1];
                        }

                    }

                    b += h;
                    h1 += h;

                    if (j == hy - 2){
                            h1 = height;
                            rgb1 = new int [w*(height - b) + 1];
                            h = height - b;
                    }


            }
            return rgb;
	}
	rgb = retinex_ms(rgb, width, height);
	return rgb;
}


}
