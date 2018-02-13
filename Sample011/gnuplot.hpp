#pragma once
#include <iostream>
#include <cmath>
#include <array>
#include "parameter.hpp"
#include "extendedArray.hpp"

using extendedArray = sksat::array_wrapper<N + 1, sksat::array_wrapper<N + 1>>;
using Array = sksat::array_wrapper<N + 1, double>;

namespace mino2357{
    class gnuplot{
    private:
        FILE* gp;
    public:
        gnuplot(){
            gp = popen("gnuplot -persist", "w");
            //fprintf(gp, "set pm3d\n");
            //fprintf(gp, "set pm3d map\n");
            //fprintf(gp, "set contour\n");
            //fprintf(gp, "set cbtics 0.05\n");
            fprintf(gp, "set xr [%f:%f]\n", xstart, xstart + Lx);
            fprintf(gp, "set yr [-1:2]\n");
            //fprintf(gp, "set yr [%f:%f]\n", ystart, ystart + Ly);
          	//fprintf(gp, "set zr [-0.1:1.1]\n");
            //fprintf(gp, "set size square\n");
            fprintf(gp, "set key outside center top samplen 0\n");
            fprintf(gp, "set tics font 'Times New Roman,12'\n");
            fprintf(gp, "set grid\n");
            fprintf(gp, "set term png size 1280, 720\n");
            //fprintf(gp, "unset key\n");
        }
        ~gnuplot(){ pclose(gp);}

        void print(extendedArray& f){
            fprintf(gp, "splot '-' w l\n");
            for(int i=0; i<=N; i = i + plus){
                for(int j=0; j<=N; j = j + plus){
                    //fprintf(gp, "%f %f %f\n", xstart + i * dx, ystart + j * dy, mino2357::speed(i, j, u1, v1));
                    fprintf(gp, "%f %f %f\n", xstart + i * dx, ystart + j * dy, f[i][j]);
                }
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\n");
            fflush(gp);
        }
        
		void print(Array& f){
            static int i = 0;
            fprintf(gp, "plot '-' w l lw 2\n");
            for(int i=0; i<=N; ++i){
                fprintf(gp, "%f %f\n", xstart + i * dx, f[i]);
            }
            fprintf(gp, "e\n");
            fprintf(gp, "set output '%06d.png'\n", i);
            fflush(gp);
            i++;
        }
        
        void printP(extendedArray& p){
            fprintf(gp, "splot '-' w l\n");
            for(int i=0; i<=N; i = i + plus){
                for(int j=0; j<=N; j = j + plus){
                    fprintf(gp, "%f %f %f\n", xstart + i * dx, ystart + j * dy, std::log10(p[i][j]));
                }
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\n");
            fflush(gp);
        }
        
        void speed(extendedArray& f, extendedArray& g){
            fprintf(gp, "splot '-' w l\n");
            for(int i=0; i<=N; i = i + plus){
                for(int j=0; j<=N; j = j + plus){
                    //fprintf(gp, "%f %f %f\n", xstart + i * dx, ystart + j * dy, mino2357::speed(i, j, u1, v1));
                    fprintf(gp, "%f %f %f\n", xstart+i*dx, ystart+j*dy, std::sqrt(f[i][j] * f[i][j] + g[i][j] * g[i][j]));
                }
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\n");
            fflush(gp);
        }
        
        void plot(double x, double y){
            fprintf(gp, "plot '-' w l\n");
            fprintf(gp, "%f %f\n", x, y);
            fprintf(gp, "\n");
            fprintf(gp, "e\n");
            fflush(gp);
        }
        
        void vector(extendedArray& f, extendedArray& g){
            fprintf(gp, "plot '-' with vector linecolor palette\n");
            for(int i=0; i<=N; i = i + plus){
                for(int j=0; j<=N; j = j + plus){
                    fprintf(gp, "%f %f %f %f %f\n", xstart+i*dx, ystart+j*dy, f[i][j], g[i][j],
                                                    std::sqrt(f[i][j] * f[i][j] + g[i][j] * g[i][j]));
                }
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\n");
            fflush(gp);
        }
        
        void vectorWithP(extendedArray& f, extendedArray& g, extendedArray& p){
            fprintf(gp, "plot '-' with vector linecolor palette\n");
            for(int i=0; i<=N; i = i + plus){
                for(int j=0; j<=N; j = j + plus){
                    double L = vecLen * std::sqrt(f[i][j] * f[i][j] + g[i][j] * g[i][j]);
                    fprintf(gp, "%f %f %f %f %f\n", xstart+i*dx, ystart+j*dy, f[i][j] / L, g[i][j] / L, p[i][j]);
                }
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\n");
            fflush(gp);
        }

        void vectorWithLog10P(extendedArray& f, extendedArray& g, extendedArray& p){
            fprintf(gp, "plot '-' with vector linecolor palette\n");
            for(int i=istart; i<=N; i = i + plus){
                for(int j=jstart; j<=N; j = j + plus){
                    double L = vecLen * std::sqrt(f[i][j] * f[i][j] + g[i][j] * g[i][j]);
                    fprintf(gp, "%f %f %f %f %f\n", xstart+i*dx, ystart+j*dy, f[i][j] / L, g[i][j] / L, std::log10(std::abs(p[i][j])));
                }
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\n");
            fflush(gp);
        }
        
        void vectorWithSpeed(extendedArray& f, extendedArray& g){
            fprintf(gp, "plot '-' with vector linecolor palette\n");
            for(int i=0; i<=N; i = i + plus){
                for(int j=0; j<=N; j = j + plus){
                    double L = vecLen * std::sqrt(f[i][j] * f[i][j] + g[i][j] * g[i][j]);
                    fprintf(gp, "%f %f %f %f %f\n", xstart+i*dx, ystart+j*dy, f[i][j] / L, g[i][j] / L, L);
                }
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\n");
            fflush(gp);
        }
        
		void vectorWithSpeedLog10(extendedArray& f, extendedArray& g){
            fprintf(gp, "plot '-' with vector linecolor palette\n");
            for(int i=0; i<=N; i = i + plus){
                for(int j=0; j<=N; j = j + plus){
                    double L = vecLen * std::sqrt(f[i][j] * f[i][j] + g[i][j] * g[i][j]);
                    fprintf(gp, "%f %f %f %f %f\n", xstart+i*dx, ystart+j*dy, f[i][j] / L, g[i][j] / L, std::log10(L));
                }
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\n");
            fflush(gp);
        }
		
		void printRho(extendedArray& rho){
            fprintf(gp, "splot '-' w l\n");
            double x, y;
			for(int i=1; i<N; i = i + plus){
                for(int j=1; j<N; j = j + plus){
                    x = xstart + i * dx;
                    y = ystart + j * dy;
					fprintf(gp, "%f %f %f\n", x, y, rho[i][j]);
                }
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\n");
            fflush(gp);
        }

        void multiplot(extendedArray& f, extendedArray& g, extendedArray& p){
            fprintf(gp, "unset key\n");
            fprintf(gp, "set label center at screen 0.5, 0.9 'label 1'\n");
            fprintf(gp, "set multiplot layout 1,2\n");
            fprintf(gp, "set title 'velocity'\n");
            fprintf(gp, "set lmargin screen 0.1\n");
            fprintf(gp, "set rmargin screen 0.45\n");
            fprintf(gp, "set tmargin screen 0.9\n");
            fprintf(gp, "set bmargin screen 0.1\n");
            fprintf(gp, "plot '-' with vector linecolor palette\n");
            for(int i=0; i<=N; i = i + plus){
                for(int j=0; j<=N; j = j + plus){
                    //double M = 0.01 * vecLen;
                    double L = std::sqrt(f[i][j] * f[i][j] + g[i][j] * g[i][j]);
                    fprintf(gp, "%f %f %f %f %f\n", xstart+i*dx, ystart+j*dy, 0.05 * f[i][j] / L, 0.05 * g[i][j] / L, std::log10(L));
                }
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\n");

            fprintf(gp, "set title 'pressure'\n");
            fprintf(gp, "set lmargin screen 0.55\n");
            fprintf(gp, "set rmargin screen 0.9\n");
            fprintf(gp, "set tmargin screen 0.9\n");
            fprintf(gp, "set bmargin screen 0.1\n");
            fprintf(gp, "splot '-' title 'pressure' w l\n");
            for(int i=0; i<=N; i = i + plus){
                for(int j=0; j<=N; j = j + plus){
                    fprintf(gp, "%f %f %f\n", xstart+i*dx, ystart+j*dy, std::log10(p[i][j]));
                }
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\n");
            fflush(gp);
            fprintf(gp, "unset multiplot\n");
        }
        
        void multiplot(extendedArray& f, extendedArray& g, extendedArray& p, double t){
            fprintf(gp, "unset key\n");
            fprintf(gp, "unset label\n");
            fprintf(gp, "set label center at screen 0.5, 0.9 'Re = %f, t = %f'\n", Re, t);
            fprintf(gp, "set multiplot layout 1,2\n");
            fprintf(gp, "set title 'velocity (log_10 scale)'\n");
            fprintf(gp, "set lmargin screen 0.1\n");
            fprintf(gp, "set rmargin screen 0.45\n");
            fprintf(gp, "set tmargin screen 0.9\n");
            fprintf(gp, "set bmargin screen 0.1\n");
            fprintf(gp, "plot '-' with vector linecolor palette\n");
            for(int i=0; i<=N; i = i + plus){
                for(int j=0; j<=N; j = j + plus){
                    //double M = 0.01 * vecLen;
                    double L = std::sqrt(f[i][j] * f[i][j] + g[i][j] * g[i][j]);
                    fprintf(gp, "%f %f %f %f %f\n", xstart+i*dx, ystart+j*dy, 0.05 * f[i][j] / L, 0.05 * g[i][j] / L, std::log10(L));
                }
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\n");

            fprintf(gp, "set title 'pressure (log_10 scale)'\n");
            fprintf(gp, "set lmargin screen 0.55\n");
            fprintf(gp, "set rmargin screen 0.9\n");
            fprintf(gp, "set tmargin screen 0.9\n");
            fprintf(gp, "set bmargin screen 0.1\n");
            fprintf(gp, "splot '-' title 'pressure' w l\n");
            for(int i=0; i<=N; i = i + plus){
                for(int j=0; j<=N; j = j + plus){
                    fprintf(gp, "%f %f %f\n", xstart+i*dx, ystart+j*dy, std::log10(p[i][j]));
                }
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\n");
            fflush(gp);
            fprintf(gp, "unset multiplot\n");
        }
       /* 
        void multiplotMakePNG(extendedArray& f, extendedArray& g, extendedArray& p, double t){
            static int i = 0;
            fprintf(gp, "set term png size 1280, 720\n");
            //fprintf(gp, "cd './test'\n");
            fprintf(gp, "unset key\n");
            fprintf(gp, "unset label\n");
            fprintf(gp, "set label center at screen 0.5, 0.9 'Re = %8.1f, t = %3.4f'\n", Re, t);
            fprintf(gp, "set multiplot layout 1,2\n");
            fprintf(gp, "set title 'velocity (log_{10} scale)'\n");
            fprintf(gp, "set lmargin screen 0.1\n");
            fprintf(gp, "set rmargin screen 0.42\n");
            fprintf(gp, "set tmargin screen 0.9\n");
            fprintf(gp, "set bmargin screen 0.1\n");
            fprintf(gp, "plot '-' with vector linecolor palette\n");
            for(int i=0; i<=N; i = i + plus){
                for(int j=0; j<=N; j = j + plus){
                    //double M = 0.01 * vecLen;
                    double L = std::sqrt(f[i][j] * f[i][j] + g[i][j] * g[i][j]);
                    fprintf(gp, "%f %f %f %f %f\n", xstart+i*dx, ystart+j*dy, 0.05 * f[i][j] / L, 0.05 * g[i][j] / L, std::log10(L));
                }
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\n");

            fprintf(gp, "set title 'pressure'\n");
            fprintf(gp, "set lmargin screen 0.58\n");
            fprintf(gp, "set rmargin screen 0.9\n");
            fprintf(gp, "set tmargin screen 0.9\n");
            fprintf(gp, "set bmargin screen 0.1\n");
            fprintf(gp, "splot '-' title 'pressure' w l\n");
            for(int i=1; i<N; i = i + plus/2){
                for(int j=1; j<N; j = j + plus/2){
                    fprintf(gp, "%f %f %f\n", xstart+i*dx, ystart+j*dy, (p[i][j]));
                }
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\n");
            fflush(gp);
            fprintf(gp, "unset multiplot\n");
            fprintf(gp, "set output '%06d.png'\n", i);
            i++;
        }
        */
		void multiplotMakePNGWithDensity(extendedArray& f, extendedArray& g, extendedArray& d, double t){
            static int i = 0;
            fprintf(gp, "set term png size 1280, 720\n");
            //fprintf(gp, "cd './test'\n");
            fprintf(gp, "unset key\n");
            fprintf(gp, "unset label\n");
            fprintf(gp, "set label center at screen 0.5, 0.9 'Re=%5.0f,     t = %3.4f'\n", Re, t);
            fprintf(gp, "set multiplot layout 1,2\n");
            fprintf(gp, "set title 'velocity (log scale)'\n");
            fprintf(gp, "set lmargin screen 0.1\n");
            fprintf(gp, "set rmargin screen 0.42\n");
            fprintf(gp, "set tmargin screen 0.9\n");
            fprintf(gp, "set bmargin screen 0.1\n");
            fprintf(gp, "plot '-' with vector linecolor palette\n");
            for(int i=0; i<=N; i = i + plus){
                for(int j=0; j<=N; j = j + plus){
                    //double M = 0.01 * vecLen;
                    double L = std::sqrt(f[i][j] * f[i][j] + g[i][j] * g[i][j]);
                    fprintf(gp, "%f %f %f %f %f\n", xstart+i*dx, ystart+j*dy, 0.05 * f[i][j] / L, 0.05 * g[i][j] / L, std::log10(L));
                }
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\n");

            fprintf(gp, "set title 'Density'\n");
            fprintf(gp, "set lmargin screen 0.58\n");
            fprintf(gp, "set rmargin screen 0.9\n");
            fprintf(gp, "set tmargin screen 0.9\n");
            fprintf(gp, "set bmargin screen 0.1\n");
            fprintf(gp, "splot '-' w l\n");
            for(int i=1; i<N; i = i + plus){
                for(int j=1; j<N; j = j + plus){
                    fprintf(gp, "%f %f %f\n", xstart+i*dx, ystart+j*dy, (d[i][j]));
                }
                fprintf(gp, "\n");
            }
            fprintf(gp, "e\n");
            fflush(gp);
            fprintf(gp, "unset multiplot\n");
            fprintf(gp, "set output '%06d.png'\n", i);
            i++;
        }
    };
}
