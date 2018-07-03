/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                 //
//                                   ///////  //     //   //////   //   ////    //////  /////                      //
//                                   //       //// ////   //       //   //  //  //      //                         //
//                                   ////     // //  //   //////   //   //  //  ////    /////                      //
//                                   //       //     //   //       //   //  //  //         //                      //
//                                   //////   //     //   //       //   ////    //////  /////                      //
//                                                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                 //
//                        Numeryczne rozwiazywanie rownan Maxwella metoda FDTD na siatce Yee                       //
//                                 z wykorzystaniem UPML, TF/SF, ADE i innych                                      //
//                                                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// v 1.744
// - dzialajace kasowanie pierwszeg pliku t= w przypadku zapisu tymczasowego (temp)
// - kosmetyka w dacie kompilacji 
// - poprawiony bug z blednym uzyciem funkcji dla struktury 3D w przypadku struktury 2D
// - zapis kroku czasu w ciagu znakow w pliku pdp (przdatne dla tmp i podgladu)
// - poprawki z itoa niedzialajcym pod gcc - itoa zastapione przez sprintf
// v 1.745
// - dalsze ciecie struktury po przesunieciu na UPML
// - poprawka w polaryzacji E z przy strukturze 2D i rezerwacji pamieci
// - poprawki w nazwie - nazwy fizyczne zamiast M mamy Jh
// v 1.746
// - wiele zrodel (w trakcie implementacji) ale skoñczona max. liczba równa max_n_SRCS
// v 1.747
// - poprawki w ustawieniu zrodla (w main() - literowki)
// - warunki brzegowe symetryczne antysymetryczne i lustrzane wprowadzone poprzez dodatkowy parametr do symetrycznych odwrocenie 
//   skladowych transwersalnych i normalnych wzgl. granicy  
// - kompilacja do wersji real albo bez zmian w kodzie, complex (wersja real g++ ..., wersja complex g++ -D COMPLEX_FIELD ...) 
// - zaleznosc czasowa zrodla takze wczytywana z zewnatrz (nieprzetestowane)
// v 1.748
// - liczne poprawki w warunkach symetrycznych, oraz uwzglednienie odwrotnego zanku symetrii dla pola H
// v 1.749
// - poprawiony bug przy konstrukcj symulacji w H i warunkami x y z
// v 1.750
// - a few corrections and clarifications in code for Portland Group Compilers
// v 1.751
// - a few style corrections 
// v 1.752
// - changing in Bloch parameters from the exact value of k to fraction of 2M_PI shift ()
// v 1.753
// - a horrible bug in behaviour of UPML has been removed by putting a line that had been accidentaly removed before
// v 1.754 - 1.755
// - apart from cross sections obsevation points are added for time series recording
// - Poynting vector components name changed from P to more common S
// - a few small changes in code
// - faster calculations for well bounded in space dispersion media
// - correction in average calculation
// v 1.756
// - calculation also on double numbers with option at compilation -D DOUBLE_PRECISION 
// v 1.757 - v 1.758
// - rewrited Bloch boundary conditions for higher numercial accuracy and stability in long simulations 
// - important change from cos to sinus in modulated source
// - a few minor corrections 
// v 1.759
// - an old Bloch form boundary conditions returns because of slowying the program flow by the new one
//   so the choose is made during complilation with option -D ACCURATE_BLOCH (but differences in result are very small...)
// v 1.760
// - complex sources added : complex ramped sinus and complex sinus gauss
// v 1.761
// - small details in code : style and priniting options with which exec was compiled and changed iosilent mode 
// v 1.762
// - further small corrections in style
// - missed earlier corrections of name Pz to Sz etc
// - bug in 2D-structure and dispersion media x-component function eliminated
// v 1.763
// - bug with detectors reserving memory eliminated
// v 1.764
// - corrected time dependency readed from file
// - corrected setting trigonometric parameters (setting ampl could affected dipol sources in some cases)
// v 1.765
// - optimised periodic parameters
// v 1.766
// - further reduction of computational domain in periodic boundary conditions for 2D Hz and Ez polarization
// v 1.767
// - small correction in buffer size  for sprintf in detectors saving/loading
// - correction in ETA calculations and behaviour of iosilent mode
// v 1.768
// - removed bug with error at saving Poynting vector components (changed the point of cenfun array initialisation)
// v 1.769
// - dipol magnetic sources added (Hx - 200, Hy-201, Hz-202)
// - corrected one more time ETA
// v 1.770
// - serious bug in speed-up of dispersive media calculation found and eliminated : now the beg end end indices are set properly
// v 1.771
// - corrected one more time ETA (didn't I write it before ?)
// v 1.772
// - differnt setting of symmetry : explictie factor for each component
// v 1.773
// - debugging the dispersion media - a line 1155 displacement of Ey with Ex was fixed
// v 1.774
// - the bug the unproper setting of lower index on yee grid of disperion medium was fixed
// v 1.775-1780
// - single poles Drude plus Lorentz model implemented (medium number 12 : 1 for Drude 2 for Lorentz)
// - fixed error in dbyEy calculation
// v 1.781-1.783
// - implemented media 13 (drude-debeye) and 23 (lorentz-debeye)
// ---------------------------------------------KNOWN BUGS------------------------------------------------------------------------------
//  - description: wiazka o rozmiarach (szerokosci) wykraczajcych nieco poza symulacje moze zwiekszac blad obliczen; gdy symulacja zawiera jedynie 
//     centralna czasc wiazki bledny rezultat jest nieuniknionynie dotyczy plane wave (zredukowane w duzym stopniu od 1.70 - 1.74)
//     source: any finite source employed in EMFIDES use only a very simple approximation for calculating the propagation of finite beam that 
//     is used for setting values on TF/SF interface; approximiation is valid for wide beams; it's not used also when incident beam don't 
//     intearct with inteface apart from the origin
//   ad-hoc solution: wiazka nie powinna "trzec" o sciany , znaczy nie powinna z dala od zrodla lezec na scianie - nie dotyczy padania na sciane 
// - produces error messages when compiled under pgCC (Portland Group Compilers) - solved in 1.750 with clarification of pow() arguments type
// - description : compiled under intel or portland group compiler ezec file operating on the complex numbers produced (at least at one case)
//  a flat field distribution , it shows some dependencies on level of optimization however the details are unknown - unsolved ; using of 
//  gcc is highly recommended
//----------------------------------------------OLSNIENIA-------------------------------------------------------------------------------
// - materialy anizotropowe nie sa potrzebne , trzeba tylko odpowiednio rozmiescic material w strukturze  (rozny dla roznych skladowych)
// - wiazke polarazyowana ko³owo uzyskuje sie jako dwa nalozone spolaryzowane liniowo zrodla rozsuniete w fazie o pol okresu


#define VERSION 1.783

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstring>

// precyzja
#ifdef DOUBLE_PRECISION
 typedef double typ_prec;
#else
 typedef float typ_prec;
#endif
// typ pola 
#ifdef COMPLEX_FIELD
 #include <complex>
 typedef std::complex<typ_prec> typ_pola;
 // for avoid rewriting the code
 #ifdef DOUBLE_PRECISION
  inline std::complex<double> operator*(float a,std::complex<double> b) {return ((double)a)*b;}
  inline std::complex<double> operator*(std::complex<double> b,float a) {return ((double)a)*b;}
 #else
  inline std::complex<float> operator*(double a,std::complex<float> b) {return ((float)a)*b;}
  inline std::complex<float> operator*(std::complex<float> b,double a) {return ((float)a)*b;}
 #endif
#else
 typedef typ_prec typ_pola;
#endif

// stale fizyczne
const double light_speed=299792458.0;
const double eps_0=8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e-12;
const double mi_0=1.2566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271e-6;

// odleglosc PML od interfejsu TF/SF
const int disPML=0;

// max liczba roznych zrodel
const int max_n_SRCS=6;

// rozmiary symulowanej siatki Yee
int xr,yr,zr;

// liczba komorek UPML (ABC)
int bound;
int boundx,boundy,boundz;
// inne parametry UPML
const typ_prec wykladnik=4.0;
const typ_prec kappa_max=1.0;
typ_prec sigma_max;

// dodatkowa przestrzeñ na TF/SF przy UPML
int bndx,bndy,bndz;

// liczba symulowanych materialow
int nmat;

// stale materialowe (wzgledne przenikalnosci magnetyczne i elektryczne)
typ_prec* eps; typ_prec* mi;
double* sig;double* sih;
typ_prec* depst;typ_prec* dmit;
typ_prec* dsigt;typ_prec* dsiht;

// liczba zrodel i ich przesiecie czasowe (bezwzg. lub wzgl.)
int n_SRCS;
double SRCS_shift[max_n_SRCS];
bool is_SRCS_shift_rel[max_n_SRCS];

// dlugosc fali
double lambda[max_n_SRCS];
double lambda2[max_n_SRCS];

// krok przestrzenny symulacji
double dr;

// krok czasowy symulacji;
long double dt;

// ilosc krokow czasowych
int maxt;

// krok startowy
int ts;

// dlugosc symulowanego impulsu
int pulse_length;

// material na brzegach
int bound_mat;

// rozmiary obliczanej siatki Yee ( obszar symulowany + brzeg (ABC) + interface TF/SF)
int xs,ys,zs;

// pomocnicze zmienne - koniec symulowanego pola / poczatek brzeg
int xk,yk,zk;

// pomocnicze zmienne - polozenie interface TF/SF w kazdym z kierunkow
int i_0,j_0,k_0;
int i_1,j_1,k_1;

// pomocnicze zmienne - "typ" sciany 0 - UPML, 1 - periodic ,2  - symmetric , 3 - Bloch
short int wx,wy,wz;

// skladowe pola EM
typ_pola*** Ex;typ_pola*** Ey;typ_pola*** Ez;
typ_pola*** Hx;typ_pola*** Hy;typ_pola*** Hz;

// dodatkowe skladowe pola na UPML

//bottom
typ_pola*** bot_Dx;typ_pola*** bot_Dy;typ_pola*** bot_Dz_lef;typ_pola*** bot_Dz_rig;
typ_pola*** bot_Bx;typ_pola*** bot_By;typ_pola*** bot_Bz_lef;typ_pola*** bot_Bz_rig;

//top
typ_pola*** top_Dx;typ_pola*** top_Dy;typ_pola*** top_Dz_lef;typ_pola*** top_Dz_rig;
typ_pola*** top_Bx;typ_pola*** top_By;typ_pola*** top_Bz_lef;typ_pola*** top_Bz_rig;

//left
typ_pola*** lef_Dx;typ_pola*** lef_Dz;
typ_pola*** lef_Bx;typ_pola*** lef_Bz;

//right
typ_pola*** rig_Dx;typ_pola*** rig_Dz;
typ_pola*** rig_Bx;typ_pola*** rig_Bz;

//front
typ_pola*** fro_Dx_lef;typ_pola*** fro_Dx_rig;
typ_pola*** fro_Dx_bot;typ_pola*** fro_Dx_top;
typ_pola*** fro_Dy;typ_pola*** fro_Dz;
typ_pola*** fro_Bx_lef;typ_pola*** fro_Bx_rig;
typ_pola*** fro_Bx_bot;typ_pola*** fro_Bx_top;
typ_pola*** fro_By;typ_pola*** fro_Bz;

//back
typ_pola*** bac_Dx_lef;typ_pola*** bac_Dx_rig;
typ_pola*** bac_Dx_bot;typ_pola*** bac_Dx_top;
typ_pola*** bac_Dy;typ_pola*** bac_Dz;
typ_pola*** bac_Bx_lef;typ_pola*** bac_Bx_rig;
typ_pola*** bac_Bx_bot;typ_pola*** bac_Bx_top;
typ_pola*** bac_By;typ_pola*** bac_Bz;


// tablice materialowe - okreslaja typ materialu w kazdym punkcie siatki
// ogolna tablica z ktorej tworzymy inne tablice materialowe dla 2D i 3D
short int** mat;short int*** mat3;

// tablice materialowe dla poszczegolnych skladowych dla 2D i 3D
short int** matEx;short int** matEy;short int** matEz;
short int** matHx;short int** matHy;short int** matHz;

short int*** matEx3;short int*** matEy3;short int*** matEz3;
short int*** matHx3;short int*** matHy3;short int*** matHz3;

// dodatkowe tablice materialowe na brzegu x
// pole elektryczne
double* c1x_E;double* c2x_E;double* c3x_E;
double* c4x_E;double* c5x_E;double* c6x_E;

// pole magnetyczne
double* c1x_H;double* c2x_H;double* c3x_H;
double* c4x_H;double* c5x_H;double* c6x_H;

// dodatkowe tablice materialowe na brzegu y
// pole elektryczne
double* c1y_E;double* c2y_E;double* c3y_E;
double* c4y_E;double* c5y_E;double* c6y_E;

// pole magnetyczne
double* c1y_H;double* c2y_H;double* c3y_H;
double* c4y_H;double* c5y_H;double* c6y_H;

// dodatkowe tablice materialowe na brzegu z
// pole elektryczne
double* c1z_E;double* c2z_E;double* c3z_E;
double* c4z_E;double* c5z_E;double* c6z_E;

// pole magnetyczne
double* c1z_H;double* c2z_H;double* c3z_H;
double* c4z_H;double* c5z_H;double* c6z_H;

// stala skalowania pola magnetycznego
typ_prec H_scale;

// pomocnicza stala dla zrodla (stymulowanego pola)
double freq_cst[max_n_SRCS];
double freq_cst2[max_n_SRCS];

// amplituda fali
typ_pola ampl[max_n_SRCS];

// typ zrodla
short int source_type[max_n_SRCS];

// parametry zrodla ograniczonego
typ_prec (*shape[max_n_SRCS])(typ_prec,typ_prec,typ_prec,int);
typ_prec param1[max_n_SRCS],param2[max_n_SRCS],param3[max_n_SRCS];

// przesuniecie czasowe impulsu (konieczne by impuls startowal od zera)
typ_prec t_0[max_n_SRCS];

// dodatkowe tablice na liczona propagacje
typ_pola* i_E[max_n_SRCS];typ_pola* i_H[max_n_SRCS];

// katy okreslajace kierunek propagacji i charakter fali
double teta[max_n_SRCS],fi[max_n_SRCS],psi[max_n_SRCS];

// stosunek dlugosci wektora falowego dla kierunku z i kata teta,fi
double angle_ratio[max_n_SRCS];
// przelicznik odleglosc - czas
double st_ratio[max_n_SRCS];
// maks. przesuniecie zrodla wzg. najblizszego i najdalszego wierzcholka
double max_t_shift[max_n_SRCS];

// przesuniecie o pol okresu
double half_period[max_n_SRCS];

// wsp. trygonometryczne dla obliczenia pola
double tr_x[max_n_SRCS],tr_y[max_n_SRCS],tr_z[max_n_SRCS];
double tr_ex[max_n_SRCS],tr_ey[max_n_SRCS],tr_ez[max_n_SRCS];
double tr_hx[max_n_SRCS],tr_hy[max_n_SRCS],tr_hz[max_n_SRCS];
double tr_sin_teta[max_n_SRCS],tr_sin_fi[max_n_SRCS],tr_cos_fi[max_n_SRCS];
double tr_xx[max_n_SRCS],tr_yy[max_n_SRCS];

// wskaznik do funkcji zrodla
typ_pola (*source_func[max_n_SRCS])(typ_prec,int);
// wskazniki do funkcji zrodla
typ_pola (*inc_Ex[max_n_SRCS])(typ_prec,typ_prec,typ_prec,int);typ_pola (*inc_Ey[max_n_SRCS])(typ_prec,typ_prec,typ_prec,int);typ_pola (*inc_Ez[max_n_SRCS])(typ_prec,typ_prec,typ_prec,int);
typ_pola (*inc_Hx[max_n_SRCS])(typ_prec,typ_prec,typ_prec,int);typ_pola (*inc_Hy[max_n_SRCS])(typ_prec,typ_prec,typ_prec,int);typ_pola (*inc_Hz[max_n_SRCS])(typ_prec,typ_prec,typ_prec,int);

// polozenie zrodla
typ_prec i_source[max_n_SRCS],j_source[max_n_SRCS],k_source[max_n_SRCS];
int floor_i_source[max_n_SRCS],floor_j_source[max_n_SRCS],floor_k_source[max_n_SRCS];
// pomocnicze zmienne do gaussa i innych zrodel ograniczonych
typ_prec gauss_c[max_n_SRCS],gauss_distance[max_n_SRCS],gauss_z0x[max_n_SRCS],gauss_z0y[max_n_SRCS],gauss_amp[max_n_SRCS];
// zrodlo wczytwyane z zewnatrz - tablice skladowych
double** shape_ex_real[max_n_SRCS];double** shape_ey_real[max_n_SRCS];
double** shape_ex_imag[max_n_SRCS];double** shape_ey_imag[max_n_SRCS];
// wzgledny krok siatki zrodla w stosunku do wlasciwej
double drs[max_n_SRCS];
// zaleznosc czasowa zrodla wczytywana z pliku
double* time_dependency[max_n_SRCS];
int time_dependency_length[max_n_SRCS];

// rzutowanie skladowych x i y w ukladzie zrodla na odpowiednie kierunki
double tr_exx[max_n_SRCS],tr_eyx[max_n_SRCS],tr_ezx[max_n_SRCS];
double tr_hxx[max_n_SRCS],tr_hyx[max_n_SRCS],tr_hzx[max_n_SRCS];
double tr_exy[max_n_SRCS],tr_eyy[max_n_SRCS],tr_ezy[max_n_SRCS];
double tr_hxy[max_n_SRCS],tr_hyy[max_n_SRCS],tr_hzy[max_n_SRCS];
double sh_cx[max_n_SRCS],sh_cy[max_n_SRCS];

// ile roznych outputow
int n_outputs;

// polozenie plaszczyzn ktore wyrzucamy
int* out_min_t;int* out_max_t;int* out_t_step;
int* slice;
int* out_min_slice;int* out_max_slice;
int* di;int* dj;int* dk;
bool** skladowe;
bool* is_fourie;bool* is_averaged;bool* is_temp;
double* ni_fourie;

// ilosc krokow czasowych w danej chwili
int t;

// ilosc krokow czasowych co ktore zrzucamy wszystko
int backup_delta_t;

// last but not least - zmienne pomocnicze
// indeksy
int i,j,k,p,q,r;
// do ABC
typ_prec c1,c2,c3,c4,c5,c6;
// zmienne pomocnicze do kreowania raportu
char struct_name[100];
time_t strtime,endtime,prevtime;

// zmienna pomocnicza
typ_prec temp;typ_prec temp2;

// czy program  uruchomiony w trybie "cichym"
bool silentmode=false;
bool iosilentmode=false;

// czy struktura 3D
bool mode3D=false;

// czy restart
bool ressurect=false;

// czy nie tworzymy raportu
bool noreport=false;

// funkcje z maina - wskazniki
void (*stage_0)();

void (*stage_E1x)(), (*stage_E1y)(), (*stage_E1z)();

void (*stage_E2x_1)(), (*stage_E3x_1)(), (*stage_E4x_1)();
void (*stage_E2y_1)(), (*stage_E3y_1)(), (*stage_E4y_1)();
void (*stage_E2z_1)(), (*stage_E3z_1)(), (*stage_E4z_1)();

void (*stage_E2x_2)(), (*stage_E3x_2)(), (*stage_E4x_2)();
void (*stage_E2y_2)(), (*stage_E3y_2)(), (*stage_E4y_2)();
void (*stage_E2z_2)(), (*stage_E3z_2)(), (*stage_E4z_2)();

void (*stage_E5_1_1[max_n_SRCS])(int), (*stage_E6_1_1[max_n_SRCS])(int), (*stage_E7_1_1[max_n_SRCS])(int);
void (*stage_E5_1_2[max_n_SRCS])(int), (*stage_E6_1_2[max_n_SRCS])(int), (*stage_E7_1_2[max_n_SRCS])(int);
void (*stage_E5_2_1[max_n_SRCS])(int), (*stage_E6_2_1[max_n_SRCS])(int), (*stage_E7_2_1[max_n_SRCS])(int);
void (*stage_E5_2_2[max_n_SRCS])(int), (*stage_E6_2_2[max_n_SRCS])(int), (*stage_E7_2_2[max_n_SRCS])(int);

void (*stage_E8x)(),(*stage_E8y)(), (*stage_E8z)();
void (*stage_E9x)(),(*stage_E9y)(), (*stage_E9z)();
void (*stage_E10x)(),(*stage_E10y)(), (*stage_E10z)();

void (*stage_J8x)(),(*stage_J8y)(), (*stage_J8z)();
void (*stage_J9x)(),(*stage_J9y)(), (*stage_J9z)();
void (*stage_J10x)(),(*stage_J10y)(), (*stage_J10z)();


void (*stage_E11x)(),(*stage_E11y)(), (*stage_E11z)();
void (*stage_E12y)(),(*stage_E12x)(), (*stage_E12z)();
void (*stage_E13x)(),(*stage_E13y)(), (*stage_E13z)();

void (*stage_E14[max_n_SRCS])(int s);


void (*stage_H1x)(), (*stage_H1y)(), (*stage_H1z)();

void (*stage_H2x_1)(), (*stage_H3x_1)(), (*stage_H4x_1)();
void (*stage_H2y_1)(), (*stage_H3y_1)(), (*stage_H4y_1)();
void (*stage_H2z_1)(), (*stage_H3z_1)(), (*stage_H4z_1)();

void (*stage_H2x_2)(), (*stage_H3x_2)(), (*stage_H4x_2)();
void (*stage_H2y_2)(), (*stage_H3y_2)(), (*stage_H4y_2)();
void (*stage_H2z_2)(), (*stage_H3z_2)(), (*stage_H4z_2)();

void (*stage_H5_1_1[max_n_SRCS])(int), (*stage_H6_1_1[max_n_SRCS])(int), (*stage_H7_1_1[max_n_SRCS])(int);
void (*stage_H5_1_2[max_n_SRCS])(int), (*stage_H6_1_2[max_n_SRCS])(int), (*stage_H7_1_2[max_n_SRCS])(int);
void (*stage_H5_2_1[max_n_SRCS])(int), (*stage_H6_2_1[max_n_SRCS])(int), (*stage_H7_2_1[max_n_SRCS])(int);
void (*stage_H5_2_2[max_n_SRCS])(int), (*stage_H6_2_2[max_n_SRCS])(int), (*stage_H7_2_2[max_n_SRCS])(int);

void (*stage_H8x)(),(*stage_H8y)(), (*stage_H8z)();
void (*stage_H9x)(),(*stage_H9y)(), (*stage_H9z)();
void (*stage_H10x)(),(*stage_H10y)(), (*stage_H10z)();

void (*stage_Jh8x)(),(*stage_Jh8y)(), (*stage_Jh8z)();
void (*stage_Jh9x)(),(*stage_Jh9y)(), (*stage_Jh9z)();
void (*stage_Jh10x)(),(*stage_Jh10y)(), (*stage_Jh10z)();


void (*stage_H11x)(),(*stage_H11y)(), (*stage_H11z)();
void (*stage_H12x)(),(*stage_H12y)(), (*stage_H12z)();
void (*stage_H13x)(),(*stage_H13y)(), (*stage_H13z)();

void (*stage_H14[max_n_SRCS])(int s);

// osrodki dyspersyjne
// drude (elektryczny)
bool* is_drd_E;
int n_drd_Ex,n_drd_Ey,n_drd_Ez;
typ_prec* drd_J_E;
long double* drd_omp_E;long double* drd_gam_E;
typ_prec* drd_be_E;typ_prec* drd_ka_E;
typ_pola* drd_Jex;typ_pola* drd_Jey;typ_pola* drd_Jez;
typ_pola* drd_Ex_n_1;typ_pola* drd_Ey_n_1;typ_pola* drd_Ez_n_1;

// osrodki dyspersyjne - drude (magnetyczny)
bool* is_drd_H;
int n_drd_Hx,n_drd_Hy,n_drd_Hz;
long double* drd_omp_H;long double* drd_gam_H;
typ_prec* drd_J_H;
typ_prec* drd_be_H;typ_prec* drd_ka_H;
typ_pola* drd_Jhx;typ_pola* drd_Jhy;typ_pola* drd_Jhz;
typ_pola* drd_Hx_n_1;typ_pola* drd_Hy_n_1;typ_pola* drd_Hz_n_1;

// lorentz (elektryczny)
bool* is_lor_E;
int n_lor_Ex,n_lor_Ey,n_lor_Ez;

long double* lor_epd_E;long double* lor_omp_E;long double* lor_del_E;
typ_prec* lor_alf_E;typ_prec* lor_ksi_E;typ_prec* lor_gam_E;
typ_prec* lor_a1_E;typ_prec* lor_a2_E;typ_prec* lor_a3_E;

typ_pola* lor_Jex;typ_pola* lor_Jey;typ_pola* lor_Jez;
typ_pola* lor_Jex_n_1;typ_pola* lor_Jey_n_1;typ_pola* lor_Jez_n_1;
typ_pola* lor_Ex_n_1;typ_pola* lor_Ey_n_1;typ_pola* lor_Ez_n_1;
typ_pola* lor_Ex_n_2;typ_pola* lor_Ey_n_2;typ_pola* lor_Ez_n_2;

// lorentz (magnetyczny)
bool* is_lor_H;
int n_lor_Hx,n_lor_Hy,n_lor_Hz;

long double* lor_epd_H;long double* lor_omp_H;long double* lor_del_H;
typ_prec* lor_alf_H;typ_prec* lor_ksi_H;typ_prec* lor_gam_H;
typ_prec* lor_a1_H;typ_prec* lor_a2_H;typ_prec* lor_a3_H;

typ_pola* lor_Jhx;typ_pola* lor_Jhy;typ_pola* lor_Jhz;
typ_pola* lor_Jhx_n_1;typ_pola* lor_Jhy_n_1;typ_pola* lor_Jhz_n_1;
typ_pola* lor_Hx_n_1;typ_pola* lor_Hy_n_1;typ_pola* lor_Hz_n_1;
typ_pola* lor_Hx_n_2;typ_pola* lor_Hy_n_2;typ_pola* lor_Hz_n_2;

// osrodki dyspersyjne
// debeye (elektryczny)
bool* is_dby_E;
int n_dby_Ex,n_dby_Ey,n_dby_Ez;
typ_prec* dby_J_E;
long double* dby_deps_E;long double* dby_tau_E;
typ_prec* dby_be_E;typ_prec* dby_ka_E;
typ_pola* dby_Jex;typ_pola* dby_Jey;typ_pola* dby_Jez;
typ_pola* dby_Ex_n_1;typ_pola* dby_Ey_n_1;typ_pola* dby_Ez_n_1;

// osrodki dyspersyjne - debeye (magnetyczny)
bool* is_dby_H;
int n_dby_Hx,n_dby_Hy,n_dby_Hz;
long double* dby_deps_H;long double* dby_tau_H;
typ_prec* dby_J_H;
typ_prec* dby_be_H;typ_prec* dby_ka_H;
typ_pola* dby_Jhx;typ_pola* dby_Jhy;typ_pola* dby_Jhz;
typ_pola* dby_Hx_n_1;typ_pola* dby_Hy_n_1;typ_pola* dby_Hz_n_1;

// osrodki dyspersyjne kombinowane 
// Drude + Lorentz 
bool* is_drd_lor_E;bool* is_drd_lor_H;
// Drude + Debye
bool* is_drd_dby_E;bool* is_drd_dby_H;
// Drude + Lorentz 
bool* is_lor_dby_E;bool* is_lor_dby_H;

// idealne przewodniki
bool* is_PEC;bool* is_PMC;

// osrodki dyspersyjne - granice;
int drd_Ex_is,drd_Ex_js,drd_Ex_ks;
int drd_Ex_ik,drd_Ex_jk,drd_Ex_kk;
int drd_Ey_is,drd_Ey_js,drd_Ey_ks;
int drd_Ey_ik,drd_Ey_jk,drd_Ey_kk;
int drd_Ez_is,drd_Ez_js,drd_Ez_ks;
int drd_Ez_ik,drd_Ez_jk,drd_Ez_kk;

int lor_Ex_is,lor_Ex_js,lor_Ex_ks;
int lor_Ex_ik,lor_Ex_jk,lor_Ex_kk;
int lor_Ey_is,lor_Ey_js,lor_Ey_ks;
int lor_Ey_ik,lor_Ey_jk,lor_Ey_kk;
int lor_Ez_is,lor_Ez_js,lor_Ez_ks;
int lor_Ez_ik,lor_Ez_jk,lor_Ez_kk;

int dby_Ex_is,dby_Ex_js,dby_Ex_ks;
int dby_Ex_ik,dby_Ex_jk,dby_Ex_kk;
int dby_Ey_is,dby_Ey_js,dby_Ey_ks;
int dby_Ey_ik,dby_Ey_jk,dby_Ey_kk;
int dby_Ez_is,dby_Ez_js,dby_Ez_ks;
int dby_Ez_ik,dby_Ez_jk,dby_Ez_kk;

int drd_Hx_is,drd_Hx_js,drd_Hx_ks;
int drd_Hx_ik,drd_Hx_jk,drd_Hx_kk;
int drd_Hy_is,drd_Hy_js,drd_Hy_ks;
int drd_Hy_ik,drd_Hy_jk,drd_Hy_kk;
int drd_Hz_is,drd_Hz_js,drd_Hz_ks;
int drd_Hz_ik,drd_Hz_jk,drd_Hz_kk;

int lor_Hx_is,lor_Hx_js,lor_Hx_ks;
int lor_Hx_ik,lor_Hx_jk,lor_Hx_kk;
int lor_Hy_is,lor_Hy_js,lor_Hy_ks;
int lor_Hy_ik,lor_Hy_jk,lor_Hy_kk;
int lor_Hz_is,lor_Hz_js,lor_Hz_ks;
int lor_Hz_ik,lor_Hz_jk,lor_Hz_kk;

int dby_Hx_is,dby_Hx_js,dby_Hx_ks;
int dby_Hx_ik,dby_Hx_jk,dby_Hx_kk;
int dby_Hy_is,dby_Hy_js,dby_Hy_ks;
int dby_Hy_ik,dby_Hy_jk,dby_Hy_kk;
int dby_Hz_is,dby_Hz_js,dby_Hz_ks;
int dby_Hz_ik,dby_Hz_jk,dby_Hz_kk;

// warunki brzegowe symmetryczne - okreslaja symetrie
typ_prec wx_Ex,wx_Ey,wx_Ez;
typ_prec wx_Hx,wx_Hy,wx_Hz;
typ_prec wy_Ex,wy_Ey,wy_Ez;
typ_prec wy_Hx,wy_Hy,wy_Hz;
typ_prec wz_Ex,wz_Ey,wz_Ez;
typ_prec wz_Hx,wz_Hy,wz_Hz;

// warunki brzegowe blocha
typ_prec bloch_kx,bloch_ky,bloch_kz;

#ifndef ACCURATE_BLOCH
        typ_pola bloch_expikx,bloch_expiky,bloch_expikz; 
        typ_pola bloch_expmikx,bloch_expmiky,bloch_expmikz;
#endif

// liczenie skladowych (we wnetrzu)
bool is_calc_ex,is_calc_ey,is_calc_ez;
bool is_calc_hx,is_calc_hy,is_calc_hz;

// wska¿niki do funkcji natezenia
typ_pola (*cenSx)(int,int,int);
typ_pola (*cenSy)(int,int,int);
typ_pola (*cenSz)(int,int,int);

// deklaracje funkcji
typ_pola cenEx(int,int,int);
typ_pola cenEy(int,int,int);
typ_pola cenEz(int,int,int);

typ_pola cenHx(int,int,int);
typ_pola cenHy(int,int,int);
typ_pola cenHz(int,int,int);

// detectors positions and records
int detector_n;
int* detector_skladowa;
int *detector_x,*detector_y,*detector_z;
int *detector_t_start,*detector_t_end,*detector_t_step;
typ_pola **recorded_at_detector;

// funkcje zrzutu i nazwy skladowych
typ_pola (*cenfun[9])(int,int,int);
char* skladname[]={"Ex","Ey","Ez","Hx","Hy","Hz","Sx","Sy","Sz","S"}; 

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// pomocnicza funkcja - tworzenie tablicy 3D,2D i 1D
bool init_array3D(typ_pola*** &array,int siz_x ,int siz_y,int siz_z)
{
    array=new(std::nothrow) typ_pola**[siz_x];
    if (!array) return false;

	array[0]=new(std::nothrow) typ_pola*[siz_x*siz_y];
    if (!array[0]) return false;

	array[0][0]=new(std::nothrow) typ_pola[siz_x*siz_y*siz_z];
 	if (!array[0][0]) return false;

    for(j=1;j<siz_y;j++) array[0][j]=array[0][j-1]+siz_z;

    for(i=1;i<siz_x;i++)
    {
			array[i]=array[i-1]+siz_y;
			array[i][0]=array[i-1][0]+siz_y*siz_z;
            for(j=1;j<siz_y;j++) array[i][j]=array[i][j-1]+siz_z;
    }

	for(i=0;i<siz_x;i++) for(j=0;j<siz_y;j++) for(k=0;k<siz_z;k++) array[i][j][k]=0.0;

    return true;
}

bool init_array3D(short int*** &array,int siz_x ,int siz_y,int siz_z)
{
    array=new(std::nothrow) short int**[siz_x];
    if (!array) return false;

	array[0]=new(std::nothrow) short int*[siz_x*siz_y];
    if (!array[0]) return false;

	array[0][0]=new(std::nothrow) short int[siz_x*siz_y*siz_z];
 	if (!array[0][0]) return false;

    for(j=1;j<siz_y;j++) array[0][j]=array[0][j-1]+siz_z;

    for(i=1;i<siz_x;i++)
    {
			array[i]=array[i-1]+siz_y;
			array[i][0]=array[i-1][0]+siz_y*siz_z;
            for(j=1;j<siz_y;j++) array[i][j]=array[i][j-1]+siz_z;
    }

	for(i=0;i<siz_x;i++) for(j=0;j<siz_y;j++) for(k=0;k<siz_z;k++) array[i][j][k]=bound_mat;

    return true;
}

bool init_array2D(short int** &array,int siz_x ,int siz_y)
{
    array=new(std::nothrow) short int*[siz_x];
    if (!array) return false;

	array[0]=new(std::nothrow) short int[siz_x*siz_y];
 	if (!array[0]) return false;

    for(i=1;i<siz_x;i++) array[i]=array[i-1]+siz_y;

	for(i=0;i<siz_x;i++) for(j=0;j<siz_y;j++) array[i][j]=bound_mat;
	return true;
}

bool init_array2D(double** &array,int siz_x ,int siz_y)
{
    array=new(std::nothrow) double*[siz_x];
    if (!array) return false;

	array[0]=new(std::nothrow) double[siz_x*siz_y];
 	if (!array[0]) return false;

    for(i=1;i<siz_x;i++) array[i]=array[i-1]+siz_y;

	for(i=0;i<siz_x;i++) for(j=0;j<siz_y;j++) array[i][j]=0.0;
	return true;
}

bool init_array1D(typ_prec* &array,int siz_x)
{
    array=new(std::nothrow) typ_prec[siz_x];
    if (!array) return false;
    for(i=0;i<siz_x;i++) array[i]=0.0;
    return true;
}

#ifdef COMPLEX_FIELD
bool init_array1D(std::complex<typ_prec>* &array,int siz_x)
{
    array=new(std::nothrow) std::complex<typ_prec>[siz_x];
    if (!array) return false;
    for(i=0;i<siz_x;i++) array[i]=std::complex<typ_prec>(0.0);
    return true;
}
#endif

#ifndef DOUBLE_PRECISION
bool init_array1D(double* &array,int siz_x)
{
    array=new(std::nothrow) double[siz_x];
    if (!array) return false;
    for(i=0;i<siz_x;i++) array[i]=0.0;
    return true;
}
#endif
bool init_array1D(bool* &array,int siz_x)
{
    array=new(std::nothrow) bool[siz_x];
    if (!array) return false;
    return true;
}

// pomocnicze funkcje - usuwanie tablic
void del_array3D(typ_pola*** &array)
{
    if (!array) return;

	delete[] array[0][0];
	delete[] array[0];

    delete[] array;
}

void del_array3D(short int*** &array)
{
    if (!array) return;

	delete[] array[0][0];
	delete[] array[0];

    delete[] array;
}

void del_array2D(short int** &array)
{
    if (!array) return;

    delete[] array[0];
    delete[] array;
}

void del_array2D(double** &array)
{
    if (!array) return;

    delete[] array[0];
    delete[] array;
}

void del_array2D(typ_pola** &array)
{
    if (!array) return;

    delete[] array[0];
    delete[] array;
}

void del_array1D(long double* &array) {if (array) delete[] array;}
#ifndef DOUBLE_PRECISION
void del_array1D(double* &array) {if (array) delete[] array;}
#endif
void del_array1D(typ_prec* &array) {if (array) delete[] array;}
void del_array1D(int* &array) {if (array) delete[] array;}
void del_array1D(bool* &array) {if (array) delete[] array;}
#ifdef COMPLEX_FIELD
       void del_array1D(std::complex<typ_prec>* &array) {if (array) delete[] array;}
#endif

// bardzo uzyteczne funkcje :)
void please_do_nothing(int s) {;}
void please_do_nuthing() {;}
typ_pola return_zero(typ_prec x,typ_prec y ,typ_prec z,int s) {return (typ_pola)0;}
typ_pola return_zer(int x,int y,int z) {return (typ_pola)0;} 
// definicja sprzezenia dla typ_prec
inline typ_prec conj(typ_prec x) {return x;}
// definicja min i max dla typu int
inline int min(int a,int b) {return a>b ? b:a;}
inline int max(int a,int b) {return a>b ? a:b;}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// wczytujemy strukture  i wpisujemy ja do tablic
bool create_structure(char* struct_file_name)
{
	int x,y,z;
	int tmp1;short int tmp2;

	// wczytujemy plik z binarna zapisana struktura (numerem materialu) - por. specyfikacja pdp
    	std::ifstream struct_file(struct_file_name,std::ios::in | std::ios::binary);
	if (!struct_file) return false;
	struct_file.read((char*)&tmp1,sizeof(int));
	struct_file.read(struct_name,tmp1);

	struct_file.read((char*)&tmp1,sizeof(int));x=tmp1;
	struct_file.read((char*)&tmp1,sizeof(int));y=tmp1;
	struct_file.read((char*)&tmp1,sizeof(int));z=tmp1;

	n_drd_Ex=0;n_drd_Ey=0;n_drd_Ez=0;
	n_drd_Hx=0;n_drd_Hy=0;n_drd_Hz=0;

	n_lor_Ex=0;n_lor_Ey=0;n_lor_Ez=0;
	n_lor_Hx=0;n_lor_Hy=0;n_lor_Hz=0;

	n_dby_Ex=0;n_dby_Ey=0;n_dby_Ez=0;
	n_dby_Hx=0;n_dby_Hy=0;n_dby_Hz=0;
	
     drd_Ex_is=0;drd_Ex_js=0;drd_Ex_ks=0;
     drd_Ex_ik=0;drd_Ex_jk=0;drd_Ex_kk=0;
     drd_Ey_is=0;drd_Ey_js=0;drd_Ey_ks=0;
     drd_Ey_ik=0;drd_Ey_jk=0;drd_Ey_kk=0;
     drd_Ez_is=0;drd_Ez_js=0;drd_Ez_ks=0;
     drd_Ez_ik=0;drd_Ez_jk=0;drd_Ez_kk=0;

     lor_Ex_is=0;lor_Ex_js=0;lor_Ex_ks=0;
     lor_Ex_ik=0;lor_Ex_jk=0;lor_Ex_kk=0;
     lor_Ey_is=0;lor_Ey_js=0;lor_Ey_ks=0;
     lor_Ey_ik=0;lor_Ey_jk=0;lor_Ey_kk=0;
     lor_Ez_is=0;lor_Ez_js=0;lor_Ez_ks=0;
     lor_Ez_ik=0;lor_Ez_jk=0;lor_Ez_kk=0;

     dby_Ex_is=0;dby_Ex_js=0;dby_Ex_ks=0;
     dby_Ex_ik=0;dby_Ex_jk=0;dby_Ex_kk=0;
     dby_Ey_is=0;dby_Ey_js=0;dby_Ey_ks=0;
     dby_Ey_ik=0;dby_Ey_jk=0;dby_Ey_kk=0;
     dby_Ez_is=0;dby_Ez_js=0;dby_Ez_ks=0;
     dby_Ez_ik=0;dby_Ez_jk=0;dby_Ez_kk=0;

     drd_Hx_is=0;drd_Hx_js=0;drd_Hx_ks=0;
     drd_Hx_ik=0;drd_Hx_jk=0;drd_Hx_kk=0;
     drd_Hy_is=0;drd_Hy_js=0;drd_Hy_ks=0;
     drd_Hy_ik=0;drd_Hy_jk=0;drd_Hy_kk=0;
     drd_Hz_is=0;drd_Hz_js=0;drd_Hz_ks=0;
     drd_Hz_ik=0;drd_Hz_jk=0;drd_Hz_kk=0;

     lor_Hx_is=0;lor_Hx_js=0;lor_Hx_ks=0;
     lor_Hx_ik=0;lor_Hx_jk=0;lor_Hx_kk=0;
     lor_Hy_is=0;lor_Hy_js=0;lor_Hy_ks=0;
     lor_Hy_ik=0;lor_Hy_jk=0;lor_Hy_kk=0;
     lor_Hz_is=0;lor_Hz_js=0;lor_Hz_ks=0;
     lor_Hz_ik=0;lor_Hz_jk=0;lor_Hz_kk=0;

     dby_Hx_is=0;dby_Hx_js=0;dby_Hx_ks=0;
     dby_Hx_ik=0;dby_Hx_jk=0;dby_Hx_kk=0;
     dby_Hy_is=0;dby_Hy_js=0;dby_Hy_ks=0;
     dby_Hy_ik=0;dby_Hy_jk=0;dby_Hy_kk=0;
     dby_Hz_is=0;dby_Hz_js=0;dby_Hz_ks=0;
     dby_Hz_ik=0;dby_Hz_jk=0;dby_Hz_kk=0;
	
	

	if (mode3D==false)
	{
		// jesli struktura 2d to stosujemy uproszczone tablice materialowe 2D
		if((z!=1)||(x!=2*xr+1)||(y!=2*yr+1)) return false; // sprawdzamy czy rozmiary struktury zgodne z zadeklarowanymi;
		struct_file.read((char*)&tmp1,sizeof(int));if (tmp1!=2) return false; // zwracamy blad jesli pole !=short int

    	for(i=0;i<=2*xr;i++) for(j=0;j<=2*yr;j++) { struct_file.read((char*)&tmp2,sizeof(short int )); mat[i][j]=tmp2; }

    	struct_file.close();

    	//matEx
        if (is_calc_ex)
        {
    	for (i=0;i<xr;i++) for (j=0;j<=yr;j++)
		{
			tmp2=mat[i+i+1][j+j];
			matEx[i+bndx+disPML][j+bndy+disPML]=tmp2;

			if (is_drd_E[tmp2]) 
            {
             n_drd_Ex++;
             drd_Ex_is=min(drd_Ex_is,i); drd_Ex_js=min(drd_Ex_js,j); 
             drd_Ex_ik=max(drd_Ex_ik,i); drd_Ex_jk=max(drd_Ex_jk,j); 
            }
			if (is_lor_E[tmp2])
            { 
             n_lor_Ex++;
             lor_Ex_is=min(lor_Ex_is,i); lor_Ex_js=min(lor_Ex_js,j); 
             lor_Ex_ik=max(lor_Ex_ik,i); lor_Ex_jk=max(lor_Ex_jk,j); 
            }
			if (is_dby_E[tmp2])
            { 
             n_dby_Ex++;
             dby_Ex_is=min(dby_Ex_is,i); lor_Ex_js=min(dby_Ex_js,j);
             dby_Ex_ik=max(dby_Ex_ik,i); dby_Ex_jk=max(dby_Ex_ik,j); 
            } 
		}
		n_drd_Ex=n_drd_Ex*(zr+2*bndz+1+2*disPML);
		n_lor_Ex=n_lor_Ex*(zr+2*bndz+1+2*disPML);
        n_dby_Ex=n_dby_Ex*(zr+2*bndz+1+2*disPML);
        drd_Ex_is+=boundx+bndx+disPML; drd_Ex_js+=boundy+bndy+disPML;
        lor_Ex_is+=boundx+bndx+disPML; lor_Ex_js+=boundy+bndy+disPML;
        dby_Ex_is+=boundx+bndx+disPML; dby_Ex_js+=boundy+bndy+disPML;
        drd_Ex_ik+=boundx+bndx+disPML; drd_Ex_jk+=boundy+bndy+disPML;
        lor_Ex_ik+=boundx+bndx+disPML; lor_Ex_jk+=boundy+bndy+disPML;
        dby_Ex_ik+=boundx+bndx+disPML; dby_Ex_jk+=boundy+bndy+disPML;
        }
        if (is_calc_ey)
        {
		// matEy
    	for (i=0;i<=xr;i++) for (j=0;j<yr;j++)
		{
			tmp2=mat[i+i][j+j+1];
			matEy[i+bndx+disPML][j+bndy+disPML]=tmp2;
			
            if (is_drd_E[tmp2]) 
            {
             n_drd_Ey++;
             drd_Ey_is=min(drd_Ey_is,i); drd_Ey_js=min(drd_Ey_js,j); 
             drd_Ey_ik=max(drd_Ey_ik,i); drd_Ey_jk=max(drd_Ey_jk,j); 
            }
			if (is_lor_E[tmp2])
            { 
             n_lor_Ey++;
             lor_Ey_is=min(lor_Ey_is,i); lor_Ey_js=min(lor_Ey_js,j);
             lor_Ey_ik=max(lor_Ey_ik,i); lor_Ey_jk=max(lor_Ey_jk,j); 
            }
			if (is_dby_E[tmp2])
            { 
             n_dby_Ey++;
             dby_Ey_is=min(dby_Ey_is,i); dby_Ey_js=min(dby_Ey_js,j);
             dby_Ey_ik=max(dby_Ey_ik,i); dby_Ey_jk=max(dby_Ey_jk,j); 
            }
		}
		
		n_drd_Ey=n_drd_Ey*(zr+2*bndz+1+2*disPML);
		n_lor_Ey=n_lor_Ey*(zr+2*bndz+1+2*disPML);
		n_dby_Ey=n_dby_Ey*(zr+2*bndz+1+2*disPML);
        drd_Ey_is+=boundx+bndx+disPML; drd_Ey_js+=boundy+bndy+disPML;
        lor_Ey_is+=boundx+bndx+disPML; lor_Ey_js+=boundy+bndy+disPML;
        dby_Ey_is+=boundx+bndx+disPML; dby_Ey_js+=boundy+bndy+disPML;
        drd_Ey_ik+=boundx+bndx+disPML; drd_Ey_jk+=boundy+bndy+disPML;
        lor_Ey_ik+=boundx+bndx+disPML; lor_Ey_jk+=boundy+bndy+disPML;
        dby_Ey_ik+=boundx+bndx+disPML; dby_Ey_jk+=boundy+bndy+disPML;
     }
     
        if (is_calc_ez)
        {     
		// matEz
    	for (i=0;i<=xr;i++) for (j=0;j<=yr;j++)
		{
			tmp2=mat[i+i][j+j];
			matEz[i+bndx+disPML][j+bndy+disPML]=tmp2;
			
            if (is_drd_E[tmp2]) 
            {
             n_drd_Ez++;
             drd_Ez_is=min(drd_Ez_is,i); drd_Ez_js=min(drd_Ez_js,j); 
             drd_Ez_ik=max(drd_Ez_ik,i); drd_Ez_jk=max(drd_Ez_jk,j); 
            }
			if (is_lor_E[tmp2])
            { 
             n_lor_Ez++;
             lor_Ez_is=min(lor_Ez_is,i); lor_Ez_js=min(lor_Ez_js,j); 
             lor_Ez_ik=max(lor_Ez_ik,i); lor_Ez_jk=max(lor_Ez_jk,j); 
            }
			if (is_dby_E[tmp2])
            { 
             n_dby_Ez++;
             dby_Ez_is=min(dby_Ez_is,i); dby_Ez_js=min(dby_Ez_js,j);
             dby_Ez_ik=max(dby_Ez_ik,i); dby_Ez_jk=max(dby_Ez_jk,j); 
            }
		}
		n_drd_Ez=n_drd_Ez*(zr+2*bndz+2*disPML);
		n_lor_Ez=n_lor_Ez*(zr+2*bndz+2*disPML);
		n_dby_Ez=n_dby_Ez*(zr+2*bndz+2*disPML);
        drd_Ez_is+=boundx+bndx+disPML; drd_Ez_js+=boundy+bndy+disPML;
        lor_Ez_is+=boundx+bndx+disPML; lor_Ez_js+=boundy+bndy+disPML;
        dby_Ez_is+=boundx+bndx+disPML; dby_Ez_js+=boundy+bndy+disPML;
        drd_Ez_ik+=boundx+bndx+disPML; drd_Ez_jk+=boundy+bndy+disPML;
        lor_Ez_ik+=boundx+bndx+disPML; lor_Ez_jk+=boundy+bndy+disPML;
        dby_Ez_ik+=boundx+bndx+disPML; dby_Ez_jk+=boundy+bndy+disPML;
        }
        
        if (is_calc_hx)
        {
    	// matHx
    	for (i=0-1;i<xr-1+1;i++) for (j=0;j<yr;j++)
		{
			tmp2=mat[i+i+2][j+j+1];
			matHx[i+bndx+disPML][j+bndy+disPML]=tmp2;
            if (is_drd_H[tmp2]) 
            {
             n_drd_Hx++;
             drd_Hx_is=min(drd_Hx_is,i); drd_Hx_js=min(drd_Hx_js,j); 
             drd_Hx_ik=max(drd_Hx_ik,i); drd_Hx_jk=max(drd_Hx_jk,j); 
            }
			if (is_lor_H[tmp2])
            { 
             n_lor_Hx++;
             lor_Hx_is=min(lor_Hx_is,i); lor_Hx_js=min(lor_Hx_js,j); 
             lor_Hx_ik=max(lor_Hx_ik,i); lor_Hx_jk=max(lor_Hx_jk,j);  
            }
			if (is_dby_H[tmp2])
            { 
             n_dby_Hx++;
             dby_Hx_is=min(dby_Hx_is,i); dby_Hx_js=min(dby_Hx_js,j); 
             dby_Hx_ik=max(dby_Hx_ik,i); dby_Hx_jk=max(dby_Hx_jk,j); 
            }
		}
		n_drd_Hx=n_drd_Hx*(zr+2*bndz+2*disPML);
		n_lor_Hx=n_lor_Hx*(zr+2*bndz+2*disPML);
		n_dby_Hx=n_dby_Hx*(zr+2*bndz+2*disPML);
        drd_Hx_is+=boundx+bndx+disPML; drd_Hx_js+=boundy+bndy+disPML;
        lor_Hx_is+=boundx+bndx+disPML; lor_Hx_js+=boundy+bndy+disPML;
        dby_Hx_is+=boundx+bndx+disPML; dby_Hx_js+=boundy+bndy+disPML;
        drd_Hx_ik+=boundx+bndx+disPML; drd_Hx_jk+=boundy+bndy+disPML;
        lor_Hx_ik+=boundx+bndx+disPML; lor_Hx_jk+=boundy+bndy+disPML;
        dby_Hx_ik+=boundx+bndx+disPML; dby_Hx_jk+=boundy+bndy+disPML;
        }
        if (is_calc_hy)
        {
		// matHy
    	for (i=0;i<xr;i++) for (j=0-1;j<yr-1+1;j++)
		{
			tmp2=mat[i+i+1][j+j+2];
			matHy[i+bndx+disPML][j+bndy+disPML]=tmp2;
            if (is_drd_H[tmp2]) 
            {
             n_drd_Hy++;
             drd_Hy_is=min(drd_Hy_is,i); drd_Hy_js=min(drd_Hy_js,j); 
             drd_Hy_ik=max(drd_Hy_ik,i); drd_Hy_jk=max(drd_Hy_jk,j); 
            }
			if (is_lor_H[tmp2])
            { 
             n_lor_Hy++;
             lor_Hy_is=min(lor_Hy_is,i); lor_Hy_js=min(lor_Hy_js,j);              
             lor_Hy_ik=max(lor_Hy_ik,i); lor_Hy_jk=max(lor_Hy_jk,j); 
            }
			if (is_dby_H[tmp2])
            { 
             n_dby_Hy++;
             dby_Hy_is=min(dby_Hy_is,i); dby_Hy_js=min(dby_Hy_js,j);
             dby_Hy_ik=max(dby_Hy_ik,i); dby_Hy_jk=max(dby_Hy_jk,j); 
            }
		}
        drd_Hy_is+=boundx+bndx+disPML; drd_Hy_js+=boundy+bndy+disPML;
        lor_Hy_is+=boundx+bndx+disPML; lor_Hy_js+=boundy+bndy+disPML;
        dby_Hy_is+=boundx+bndx+disPML; dby_Hy_js+=boundy+bndy+disPML;
        drd_Hy_ik+=boundx+bndx+disPML; drd_Hy_jk+=boundy+bndy+disPML;
        lor_Hy_ik+=boundx+bndx+disPML; lor_Hy_jk+=boundy+bndy+disPML;
        dby_Hy_ik+=boundx+bndx+disPML; dby_Hy_jk+=boundy+bndy+disPML;
		n_drd_Hy=n_drd_Hy*(zr+2*bndz+2*disPML);
		n_lor_Hy=n_lor_Hy*(zr+2*bndz+2*disPML);
		n_dby_Hy=n_dby_Hy*(zr+2*bndz+2*disPML);
        }
        if (is_calc_hz)
        {
    	// matHz
    	for (i=0;i<xr;i++) for (j=0;j<yr;j++)
		{
			tmp2=mat[i+i+1][j+j+1];
			matHz[i+bndx+disPML][j+bndy+disPML]=tmp2;
            if (is_drd_H[tmp2]) 
            {
             n_drd_Hz++;
             drd_Hz_is=min(drd_Hz_is,i); drd_Hz_js=min(drd_Hz_js,j); 
             drd_Hz_ik=max(drd_Hz_ik,i); drd_Hz_jk=max(drd_Hz_jk,j); 
            }
			if (is_lor_H[tmp2])
            { 
             n_lor_Hz++;
             lor_Hz_is=min(lor_Hz_is,i); lor_Hz_js=min(lor_Hz_js,j); 
             lor_Hz_ik=max(lor_Hz_ik,i); lor_Hz_jk=max(lor_Hz_jk,j);  
            }
			if (is_dby_H[tmp2])
            { 
             n_dby_Hz++;
             dby_Hz_is=min(dby_Hz_is,i); dby_Hz_js=min(dby_Hz_js,j);
             dby_Hz_ik=max(dby_Hz_ik,i); dby_Hz_jk=max(dby_Hz_jk,j); 
            }
		}
		n_drd_Hz=n_drd_Hz*(zr+2*bndz-1+2*disPML);
		n_lor_Hz=n_lor_Hz*(zr+2*bndz-1+2*disPML);
		n_dby_Hz=n_dby_Hz*(zr+2*bndz-1+2*disPML);
        drd_Hz_is+=boundx+bndx+disPML; drd_Hz_js+=boundy+bndy+disPML;
        lor_Hz_is+=boundx+bndx+disPML; lor_Hz_js+=boundy+bndy+disPML;
        dby_Hz_is+=boundx+bndx+disPML; dby_Hz_js+=boundy+bndy+disPML;
        drd_Hz_ik+=boundx+bndx+disPML; drd_Hz_jk+=boundy+bndy+disPML;
        lor_Hz_ik+=boundx+bndx+disPML; lor_Hz_jk+=boundy+bndy+disPML;
        dby_Hz_ik+=boundx+bndx+disPML; dby_Hz_jk+=boundy+bndy+disPML;
        }
    	del_array2D(mat);
	}
	else
	{
		// jesli struktura 3d to i takiez tablice
		if((z!=2*zr+1)||(x!=2*xr+1)||(y!=2*yr+1)) return false; // sprawdzamy czy rozmiary struktury zgodne z zadeklarowanymi
		struct_file.read((char*)&tmp1,sizeof(int));if (tmp1!=2) return false; // zwracamy blad jesli pole !=short int

    	for(i=0;i<=2*xr;i++) for(j=0;j<=2*yr;j++) for(k=0;k<=2*zr;k++) {struct_file.read((char*)&tmp2,sizeof(short int )); mat3[i][j][k]=tmp2;}
    	struct_file.close();

    	//matEx3
    	if (is_calc_ex)
    	{
    	for(i=0;i<xr;i++) for(j=0;j<=yr;j++) for(k=0;k<=zr;k++)
		{
			tmp2=mat3[i+i+1][j+j][k+k];
			matEx3[i+bndx+disPML][j+bndy+disPML][k+bndz+disPML]=tmp2;
			if (is_drd_E[tmp2]) 
            {
             n_drd_Ex++;
             drd_Ex_is=min(drd_Ex_is,i); drd_Ex_js=min(drd_Ex_js,j); drd_Ex_ks=min(drd_Ex_ks,k); 
             drd_Ex_ik=max(drd_Ex_ik,i); drd_Ex_jk=max(drd_Ex_jk,j); drd_Ex_kk=max(drd_Ex_kk,k);  
            }
			if (is_lor_E[tmp2])
            { 
             n_lor_Ex++;
             lor_Ex_is=min(lor_Ex_is,i); lor_Ex_js=min(lor_Ex_js,j); lor_Ex_ks=min(lor_Ex_ks,k);
             lor_Ex_ik=max(lor_Ex_ik,i); lor_Ex_jk=max(lor_Ex_jk,j); lor_Ex_kk=max(lor_Ex_kk,k); 
            }
			if (is_dby_E[tmp2])
            { 
             n_dby_Ex++;
             dby_Ex_is=min(dby_Ex_is,i); dby_Ex_js=min(dby_Ex_js,j); dby_Ex_ks=min(dby_Ex_ks,k);
             dby_Ex_ik=max(dby_Ex_ik,i); dby_Ex_jk=max(dby_Ex_jk,j); dby_Ex_kk=max(dby_Ex_kk,k); 
            }
    	}
        drd_Ex_is+=boundx+bndx+disPML; drd_Ex_js+=boundy+bndy+disPML; drd_Ex_ks+=boundz+bndz+disPML;
        lor_Ex_is+=boundx+bndx+disPML; lor_Ex_js+=boundy+bndy+disPML; lor_Ex_ks+=boundz+bndz+disPML;
        dby_Ex_is+=boundx+bndx+disPML; dby_Ex_js+=boundy+bndy+disPML; dby_Ex_ks+=boundz+bndz+disPML;
        drd_Ex_ik+=boundx+bndx+disPML; drd_Ex_jk+=boundy+bndy+disPML; drd_Ex_kk+=boundz+bndz+disPML;
        lor_Ex_ik+=boundx+bndx+disPML; lor_Ex_jk+=boundy+bndy+disPML; lor_Ex_kk+=boundz+bndz+disPML;
        dby_Ex_ik+=boundx+bndx+disPML; dby_Ex_jk+=boundy+bndy+disPML; dby_Ex_kk+=boundz+bndz+disPML;
        }
        
        if (is_calc_ey)
        {
		// matEy3
    	for(i=0;i<=xr;i++) for(j=0;j<yr;j++) for(k=0;k<=zr;k++)
		{
			tmp2=mat3[i+i][j+j+1][k+k];
			matEy3[i+bndx+disPML][j+bndy+disPML][k+bndz+disPML]=tmp2;
			if (is_drd_E[tmp2]) 
            {
             n_drd_Ey++;
             drd_Ey_is=min(drd_Ey_is,i); drd_Ey_js=min(drd_Ey_js,j); drd_Ey_ks=min(drd_Ey_ks,k); 
             drd_Ey_ik=max(drd_Ey_ik,i); drd_Ey_jk=max(drd_Ey_jk,j); drd_Ey_kk=max(drd_Ey_kk,k);  
            }
			if (is_lor_E[tmp2])
            { 
             n_lor_Ey++;
             lor_Ey_is=min(lor_Ey_is,i); lor_Ey_js=min(lor_Ey_js,j); lor_Ey_ks=min(lor_Ey_ks,k);  
             lor_Ey_ik=max(lor_Ey_ik,i); lor_Ey_jk=max(lor_Ey_jk,j); lor_Ey_kk=max(lor_Ey_kk,k); 
            }
			if (is_dby_E[tmp2])
            {  
             n_dby_Ey++;
             dby_Ey_is=min(dby_Ey_is,i); dby_Ey_js=min(dby_Ey_js,j); dby_Ey_ks=min(dby_Ey_ks,k);
             dby_Ey_ik=max(dby_Ey_ik,i); dby_Ey_jk=max(dby_Ey_jk,j); dby_Ey_kk=max(dby_Ey_kk,k); 
            }
		}
        drd_Ey_is+=boundx+bndx+disPML; drd_Ey_js+=boundy+bndy+disPML; drd_Ey_ks+=boundz+bndz+disPML;
        lor_Ey_is+=boundx+bndx+disPML; lor_Ey_js+=boundy+bndy+disPML; lor_Ey_ks+=boundz+bndz+disPML;
        dby_Ey_is+=boundx+bndx+disPML; dby_Ey_js+=boundy+bndy+disPML; dby_Ey_ks+=boundz+bndz+disPML;
        drd_Ey_ik+=boundx+bndx+disPML; drd_Ey_jk+=boundy+bndy+disPML; drd_Ey_kk+=boundz+bndz+disPML;
        lor_Ey_ik+=boundx+bndx+disPML; lor_Ey_jk+=boundy+bndy+disPML; lor_Ey_kk+=boundz+bndz+disPML;
        dby_Ey_ik+=boundx+bndx+disPML; dby_Ey_jk+=boundy+bndy+disPML; dby_Ey_kk+=boundz+bndz+disPML;
        }
        
        if (is_calc_ez)
        {
    	// matEz3
    	for(i=0;i<=xr;i++) for(j=0;j<=yr;j++) for(k=0;k<zr;k++)
		{
			tmp2=mat3[i+i][j+j][k+k+1];
			matEz3[i+bndx+disPML][j+bndy+disPML][k+bndz+disPML]=tmp2;
			if (is_drd_E[tmp2]) 
            {
             n_drd_Ez++;
             drd_Ez_is=min(drd_Ez_is,i); drd_Ez_js=min(drd_Ez_js,j); drd_Ez_ks=min(drd_Ez_ks,k); 
             drd_Ez_ik=max(drd_Ez_ik,i); drd_Ez_jk=max(drd_Ez_jk,j); drd_Ez_kk=max(drd_Ez_kk,k);  
            }
			if (is_lor_E[tmp2])
            { 
             n_lor_Ez++;
             lor_Ez_is=min(lor_Ez_is,i); lor_Ez_js=min(lor_Ez_js,j); lor_Ez_ks=min(lor_Ez_ks,k);
             lor_Ez_ik=max(lor_Ez_ik,i); lor_Ez_jk=max(lor_Ez_jk,j); lor_Ez_kk=max(lor_Ez_kk,k); 
            }
			if (is_dby_E[tmp2])
            { 
             n_dby_Ez++;
             dby_Ez_is=min(dby_Ez_is,i); dby_Ez_js=min(dby_Ez_js,j); dby_Ez_ks=min(dby_Ez_ks,k);
             dby_Ez_ik=max(dby_Ez_ik,i); dby_Ez_jk=max(dby_Ez_jk,j); dby_Ez_kk=max(dby_Ez_kk,k);  
            }
		}
        drd_Ez_is+=boundx+bndx+disPML; drd_Ez_js+=boundy+bndy+disPML; drd_Ez_ks+=boundz+bndz+disPML;
        lor_Ez_is+=boundx+bndx+disPML; lor_Ez_js+=boundy+bndy+disPML; lor_Ez_ks+=boundz+bndz+disPML;
        dby_Ez_is+=boundx+bndx+disPML; dby_Ez_js+=boundy+bndy+disPML; dby_Ez_ks+=boundz+bndz+disPML;
        drd_Ez_ik+=boundx+bndx+disPML; drd_Ez_jk+=boundy+bndy+disPML; drd_Ez_kk+=boundz+bndz+disPML;
        lor_Ez_ik+=boundx+bndx+disPML; lor_Ez_jk+=boundy+bndy+disPML; lor_Ez_kk+=boundz+bndz+disPML;
        dby_Ez_ik+=boundx+bndx+disPML; dby_Ez_jk+=boundy+bndy+disPML; dby_Ez_kk+=boundz+bndz+disPML;
        }
        if (is_calc_hx)
        {
    	// matHx3
    	for(i=0-1;i<xr-1+1;i++) for(j=0;j<yr;j++) for(k=0;k<zr;k++)
		{
			tmp2=mat3[i+i+2][j+j+1][k+k+1];
			matHx3[i+bndx+disPML][j+bndy+disPML][k+bndz+disPML]=tmp2;
			if (is_drd_H[tmp2]) 
            {
             n_drd_Hx++;
             drd_Hx_is=min(drd_Hx_is,i); drd_Hx_js=min(drd_Hx_js,j); drd_Hx_ks=min(drd_Hx_ks,k); 
             drd_Hx_ik=max(drd_Hx_ik,i); drd_Hx_jk=max(drd_Hx_jk,j); drd_Hx_kk=max(drd_Hx_kk,k);
            }
			if (is_lor_H[tmp2])
            { 
             n_lor_Hx++;
             lor_Hx_is=min(lor_Hx_is,i); lor_Hx_js=min(lor_Hx_js,j); lor_Hx_ks=min(lor_Hx_ks,k);
             lor_Hx_ik=max(lor_Hx_ik,i); lor_Hx_jk=max(lor_Hx_jk,j); lor_Hx_kk=max(lor_Hx_kk,k);
            }
			if (is_dby_H[tmp2])
            {  
             n_dby_Hx++;
             dby_Hx_is=min(dby_Hx_is,i); dby_Hx_js=min(dby_Hx_js,j); dby_Hx_ks=min(dby_Hx_ks,k);
             dby_Hx_ik=max(dby_Hx_ik,i); dby_Hx_jk=max(dby_Hx_jk,j); dby_Hx_kk=max(dby_Hx_kk,k);  
            }
		}
        drd_Hx_is+=boundx+bndx+disPML; drd_Hx_js+=boundy+bndy+disPML; drd_Hx_ks+=boundz+bndz+disPML;
        lor_Hx_is+=boundx+bndx+disPML; lor_Hx_js+=boundy+bndy+disPML; lor_Hx_ks+=boundz+bndz+disPML;
        dby_Hx_is+=boundx+bndx+disPML; dby_Hx_js+=boundy+bndy+disPML; dby_Hx_ks+=boundz+bndz+disPML;
        drd_Hx_ik+=boundx+bndx+disPML; drd_Hx_jk+=boundy+bndy+disPML; drd_Hx_kk+=boundz+bndz+disPML;
        lor_Hx_ik+=boundx+bndx+disPML; lor_Hx_jk+=boundy+bndy+disPML; lor_Hx_kk+=boundz+bndz+disPML;
        dby_Hx_ik+=boundx+bndx+disPML; dby_Hx_jk+=boundy+bndy+disPML; dby_Hx_kk+=boundz+bndz+disPML;
        }
        
        if (is_calc_hy)
        {
    	// matHy3
    	for(i=0;i<xr;i++) for(j=0-1;j<yr-1+1;j++) for(k=0;k<zr;k++)
		{
			tmp2=mat3[i+i+1][j+j+2][k+k+1];
			matHy3[i+bndx+disPML][j+bndy+disPML][k+bndz+disPML]=tmp2;
			if (is_drd_H[tmp2]) 
            {
             n_drd_Hy++;
             drd_Hy_is=min(drd_Hy_is,i); drd_Hy_js=min(drd_Hy_js,j); drd_Hy_ks=min(drd_Hy_ks,k);
             drd_Hy_ik=max(drd_Hy_ik,i); drd_Hy_jk=max(drd_Hy_jk,j); drd_Hy_kk=max(drd_Hy_kk,k);
            }
			if (is_lor_H[tmp2])
            { 
             n_lor_Hy++;
             lor_Hy_is=min(lor_Hy_is,i); lor_Hy_js=min(lor_Hy_js,j); lor_Hy_ks=min(lor_Hy_ks,k);
             lor_Hy_ik=max(lor_Hy_ik,i); lor_Hy_jk=max(lor_Hy_jk,j); lor_Hy_kk=max(lor_Hy_kk,k);
            }
			if (is_dby_H[tmp2])
            { 
             n_dby_Hy++;
             dby_Hy_is=min(dby_Hy_is,i); dby_Hy_js=min(dby_Hy_js,j); dby_Hy_ks=min(dby_Hy_ks,k);
             dby_Hy_ik=max(dby_Hy_ik,i); dby_Hy_jk=max(dby_Hy_jk,j); dby_Hy_kk=max(dby_Hy_kk,k);
            }
		}
        drd_Hy_is+=boundx+bndx+disPML; drd_Hy_js+=boundy+bndy+disPML; drd_Hy_ks+=boundz+bndz+disPML;
        lor_Hy_is+=boundx+bndx+disPML; lor_Hy_js+=boundy+bndy+disPML; lor_Hy_ks+=boundz+bndz+disPML;
        dby_Hy_is+=boundx+bndx+disPML; dby_Hy_js+=boundy+bndy+disPML; dby_Hy_ks+=boundz+bndz+disPML;
        drd_Hy_ik+=boundx+bndx+disPML; drd_Hy_jk+=boundy+bndy+disPML; drd_Hy_kk+=boundz+bndz+disPML;
        lor_Hy_ik+=boundx+bndx+disPML; lor_Hy_jk+=boundy+bndy+disPML; lor_Hy_kk+=boundz+bndz+disPML;
        dby_Hy_ik+=boundx+bndx+disPML; dby_Hy_jk+=boundy+bndy+disPML; dby_Hy_kk+=boundz+bndz+disPML;
    }
    
       if (is_calc_hz)
       {
    	// matHz3
    	for(i=0;i<xr;i++) for(j=0;j<yr;j++) for(k=0-1;k<zr-1+1;k++)
		{
			tmp2=mat3[i+i+1][j+j+1][k+k+2];
			matHz3[i+bndx+disPML][j+bndy+disPML][k+bndz+disPML]=tmp2;
			tmp2=mat3[i+i+1][j+j+2][k+k+1];
			matHz3[i+bndx+disPML][j+bndy+disPML][k+bndz+disPML]=tmp2;
			if (is_drd_H[tmp2]) 
            {
             n_drd_Hz++;
             drd_Hz_is=min(drd_Hz_is,i); drd_Hz_js=min(drd_Hz_js,j); drd_Hz_ks=min(drd_Hz_ks,k);
             drd_Hz_ik=max(drd_Hz_ik,i); drd_Hz_jk=max(drd_Hz_jk,j); drd_Hz_kk=max(drd_Hz_kk,k);
            }
			if (is_lor_H[tmp2])
            { 
             n_lor_Hz++;
             lor_Hz_is=min(lor_Hz_is,i); lor_Hz_js=min(lor_Hz_js,j); lor_Hz_ks=min(lor_Hz_ks,k);
             lor_Hz_ik=max(lor_Hz_ik,i); lor_Hz_jk=max(lor_Hz_jk,j); lor_Hz_kk=max(lor_Hz_kk,k);
            }
			if (is_dby_H[tmp2])
            { 
             n_dby_Hz++;
             dby_Hz_is=min(dby_Hz_is,i); dby_Hz_js=min(dby_Hz_js,j); dby_Hz_ks=min(dby_Hz_ks,k);
             dby_Hz_ik=max(dby_Hz_ik,i); dby_Hz_jk=max(dby_Hz_jk,j); dby_Hz_kk=max(dby_Hz_kk,k);
            }
		}
        drd_Hz_is+=boundx+bndx+disPML; drd_Hz_js+=boundy+bndy+disPML; drd_Hz_ks+=boundz+bndz+disPML;
        lor_Hz_is+=boundx+bndx+disPML; lor_Hz_js+=boundy+bndy+disPML; lor_Hz_ks+=boundz+bndz+disPML;
        dby_Hz_is+=boundx+bndx+disPML; dby_Hz_js+=boundy+bndy+disPML; dby_Hz_ks+=boundz+bndz+disPML;
        drd_Hz_ik+=boundx+bndx+disPML; drd_Hz_jk+=boundy+bndy+disPML; drd_Hz_kk+=boundz+bndz+disPML;
        lor_Hz_ik+=boundx+bndx+disPML; lor_Hz_jk+=boundy+bndy+disPML; lor_Hz_kk+=boundz+bndz+disPML;
        dby_Hz_ik+=boundx+bndx+disPML; dby_Hz_jk+=boundy+bndy+disPML; dby_Hz_kk+=boundz+bndz+disPML;
       }
    	del_array3D(mat3);
	}
	return true;
}

// indeksowanie tablic materialowych
//
//
//                     0     1     2    3     4
// indeksy  --->     ------------------------/-                            np.:
// w tablicy        /|          /|          /|                             Ex(0,0,0)  -  mat(1,0)
// mat             Ez         Ez |         Ez                              Ex(1,0,0)  -  mat(3,0)
//                /  |        /  |        /  |                             Ex(0,1,0)  -  mat(1,2) itd.
//    |       0  ---- Ex -------- Ex --- /-  |
//    |          |   |       | Hx|       | Hx|                             Hx(0,0,0)  -  mat(2,1)
//    -->        |   |       |   |       |   |               z             Hx(1,0,0)  -  mat(4,1)
//            1  Ey  |------ Ey - ------ Ey -|-              /             Hx(0,1,0)  -  mat(2,3) itd.
//               |  Ez   Hy  |  Ez   Hy  |  Ez              /
//               | / |       | / |       | / |             ------>         stad funkcja przejscia z ogolnej
//            2  |--- Ex --- |---- Ex ---|-  |             |      x        tablicy materialowej do tablic dla
//               |   |       |   |       |   |             |               poszczegolnych skladowych jak wyzej
//               |           |           |                 V y

// obliczanie stalych wykorzystywanych przy obliczeniach pola
// roznica w stosunku do zwyklych stalych uwaga na przeskalowanie pola H (parametr H_scale)
void create_cst()
{
	for(i=0;i<nmat;i++)
	{
		depst[i]=((dt/(eps[i]*eps_0*dr))/H_scale)/(1+(sig[i]*dt/(2*eps[i]*eps_0)));
		dsigt[i]=(1-((sig[i]*dt)/(2*eps[i]*eps_0)))/(1+((sig[i]*dt)/(2*eps[i]*eps_0)));

		dmit[i]=((dt/(mi[i]*mi_0*dr))*H_scale)/(1+((sih[i]*dt)/(2*mi[i]*mi_0)));
		dsiht[i]=(1-((sih[i]*dt)/(2*mi[i]*mi_0)))/(1+((sih[i]*dt)/(2*mi[i]*mi_0)));

		if (is_drd_E[i])
		{
			drd_ka_E[i]=(2-drd_gam_E[i]*dt)/(2+drd_gam_E[i]*dt);
			drd_be_E[i]=0.5*(2*drd_omp_E[i]*dt*drd_omp_E[i])*eps[i]*eps_0/(2+drd_gam_E[i]*dt);

			depst[i]=((dt/(eps[i]*eps_0*dr))/H_scale)/(1+((sig[i]*dt+drd_be_E[i]*dt)/(2*eps[i]*eps_0)));
			dsigt[i]=(1-((sig[i]*dt+drd_be_E[i]*dt)/(2*eps[i]*eps_0)))/(1+((sig[i]*dt+drd_be_E[i]*dt)/(2*eps[i]*eps_0)));

			drd_J_E[i]=0.5*dr*depst[i]*(1+drd_ka_E[i])*H_scale*drd_be_E[i];	
		}

		if (is_lor_E[i])
		{
			lor_alf_E[i]=(2-pow(lor_omp_E[i]*dt,2))/(1+lor_del_E[i]*dt);
			lor_ksi_E[i]=(lor_del_E[i]*dt-1)/(lor_del_E[i]*dt+1);
			lor_gam_E[i]=(eps_0*lor_epd_E[i]*pow(lor_omp_E[i]*dt,2)/(lor_del_E[i]*dt+1));

			depst[i]=((dt/(eps[i]*eps_0*dr))/H_scale)/(1+((sig[i]*dt+0.5*lor_gam_E[i])/(2*eps[i]*eps_0)));
			dsigt[i]=(1-(sig[i]*dt/(2*eps[i]*eps_0)))/(1+((sig[i]*dt+0.5*lor_gam_E[i])/(2*eps[i]*eps_0)));

			lor_a1_E[i]=0.5*lor_gam_E[i]/(2*eps_0*eps[i]+0.5*lor_gam_E[i]+sig[i]*dt);
			lor_a2_E[i]=0.5*dr*depst[i]*H_scale*(1+lor_alf_E[i])*(lor_gam_E[i]/(2*dt));
			lor_a3_E[i]=0.5*dr*depst[i]*H_scale*lor_ksi_E[i]*(lor_gam_E[i]/(2*dt));
		}

		if (is_dby_E[i])
		{
			dby_ka_E[i]=(2-dt/dby_tau_E[i])/(2+dt/dby_tau_E[i]);
			dby_be_E[i]=(2*dby_deps_E[i]*dt/dby_tau_E[i])*eps_0/(2+dt/dby_tau_E[i]);

			depst[i]=((dt/(eps[i]*eps_0*dr))/H_scale)/(1+((sig[i]*dt+dby_be_E[i])/(2*eps[i]*eps_0)));
			dsigt[i]=(1-((sig[i]*dt-dby_be_E[i])/(2*eps[i]*eps_0)))/(1+((sig[i]*dt+dby_be_E[i])/(2*eps[i]*eps_0)));

			dby_J_E[i]=0.5*dr*depst[i]*(1+dby_ka_E[i])*H_scale*(dby_be_E[i]/dt);
		}
		
		if (is_drd_lor_E[i])
		{	
            depst[i]=((dt/(eps[i]*eps_0*dr))/H_scale)/(1+((sig[i]*dt+drd_be_E[i]*dt+0.5*lor_gam_E[i])/(2*eps[i]*eps_0)));
			dsigt[i]=(1-((sig[i]*dt+drd_be_E[i]*dt)/(2*eps[i]*eps_0)))/(1+((sig[i]*dt+drd_be_E[i]*dt+0.5*lor_gam_E[i])/(2*eps[i]*eps_0)));     
        }
        
 		if (is_drd_dby_E[i])
		{	
            depst[i]=((dt/(eps[i]*eps_0*dr))/H_scale)/(1+((sig[i]*dt+drd_be_E[i]*dt+dby_be_E[i])/(2*eps[i]*eps_0)));
			dsigt[i]=(1-((sig[i]*dt+drd_be_E[i]*dt-dby_be_E[i])/(2*eps[i]*eps_0)))/(1+((sig[i]*dt+drd_be_E[i]*dt+dby_be_E[i])/(2*eps[i]*eps_0)));     
        }

		if (is_lor_dby_E[i])
		{	
            depst[i]=((dt/(eps[i]*eps_0*dr))/H_scale)/(1+((sig[i]*dt+0.5*lor_gam_E[i]+dby_be_E[i])/(2*eps[i]*eps_0)));
			dsigt[i]=(1-((sig[i]*dt-dby_be_E[i])/(2*eps[i]*eps_0)))/(1+((sig[i]*dt+0.5*lor_gam_E[i]+dby_be_E[i])/(2*eps[i]*eps_0)));     
        }

        

		if (is_PEC[i])
		{
			depst[i]=0.0;
			dsigt[i]=-1.0;
		}

		if (is_drd_H[i])
		{
			drd_ka_H[i]=(2-drd_gam_H[i]*dt)/(2+drd_gam_H[i]*dt);
			drd_be_H[i]=0.5*(2*drd_omp_H[i]*dt*drd_omp_H[i])*mi[i]*mi_0/(2+drd_gam_H[i]*dt);

			dmit[i]=((dt/(mi[i]*mi_0*dr))*H_scale)/(1+((sih[i]*dt+drd_be_H[i]*dt)/(2*mi[i]*mi_0)));
			dsiht[i]=(1-((sih[i]*dt+drd_be_H[i]*dt)/(2*mi[i]*mi_0)))/(1+((sih[i]*dt+drd_be_H[i]*dt)/(2*mi[i]*mi_0)));

			drd_J_H[i]=0.5*dr*dmit[i]*(1+drd_ka_H[i])*drd_be_H[i]/H_scale;
		}
		if (is_lor_H[i])
		{
			lor_alf_H[i]=(2-pow(lor_omp_H[i]*dt,2))/(1+lor_del_H[i]*dt);
			lor_ksi_H[i]=(lor_del_H[i]*dt-1)/(lor_del_H[i]*dt+1);
			lor_gam_H[i]=(mi_0*lor_epd_H[i]*pow(lor_omp_H[i]*dt,2)/(lor_del_H[i]*dt+1));

			dmit[i]=((dt/(mi[i]*mi_0*dr))*H_scale)/(1+((sih[i]*dt+0.5*lor_gam_H[i])/(2*mi[i]*mi_0)));
			dsiht[i]=(1-(sih[i]*dt/(2*mi[i]*mi_0)))/(1+((sih[i]*dt+0.5*lor_gam_H[i])/(2*mi[i]*mi_0)));

			lor_a1_H[i]=0.5*lor_gam_H[i]/(2*mi_0*mi[i]+0.5*lor_gam_H[i]+sih[i]*dt);
			lor_a2_H[i]=0.5*dr*dmit[i]/H_scale*(1+lor_alf_H[i])*lor_gam_H[i]*(1/(2*dt));
			lor_a3_H[i]=0.5*dr*dmit[i]/H_scale*lor_ksi_H[i]*lor_gam_H[i]*(1/(2*dt));
		}

		if (is_dby_H[i])
		{
			dby_ka_H[i]=(2-dt/dby_tau_H[i])/(2+dt/dby_tau_H[i]);
			dby_be_H[i]=(2*dby_deps_H[i]*dt/dby_tau_H[i])*mi_0/(2+dt/dby_tau_H[i]);

			dmit[i]=((dt/(mi[i]*mi_0*dr))*H_scale)/(1+((sih[i]*dt+dby_be_H[i])/(2*mi[i]*mi_0)));
			dsiht[i]=(1-((sih[i]*dt-dby_be_H[i])/(2*mi[i]*mi_0)))/(1+((sih[i]*dt+dby_be_H[i])/(2*mi[i]*mi_0)));

			dby_J_H[i]=(0.5*dr*dmit[i]*(1+dby_ka_H[i])/H_scale)*dby_be_H[i]/dt;
		}
		
        if (is_drd_lor_H[i])
		{
			dmit[i]=((dt/(mi[i]*mi_0*dr))*H_scale)/(1+((sih[i]*dt+drd_be_H[i]*dt+0.5*lor_gam_H[i])/(2*mi[i]*mi_0)));
			dsiht[i]=(1-((sih[i]*dt+drd_be_H[i]*dt)/(2*mi[i]*mi_0)))/(1+((sih[i]*dt+drd_be_H[i]*dt+0.5*lor_gam_H[i])/(2*mi[i]*mi_0)));
        }
        
        if (is_drd_dby_H[i])
		{
			dmit[i]=((dt/(mi[i]*mi_0*dr))*H_scale)/(1+((sih[i]*dt+drd_be_H[i]*dt+dby_be_H[i])/(2*mi[i]*mi_0)));
			dsiht[i]=(1-((sih[i]*dt+drd_be_H[i]*dt-dby_be_H[i])/(2*mi[i]*mi_0)))/(1+((sih[i]*dt+drd_be_H[i]*dt+dby_be_H[i])/(2*mi[i]*mi_0)));                    
        }

        if (is_lor_dby_H[i])
		{
			dmit[i]=((dt/(mi[i]*mi_0*dr))*H_scale)/(1+((sih[i]*dt+0.5*lor_gam_H[i]+dby_be_H[i])/(2*mi[i]*mi_0)));
			dsiht[i]=(1-((sih[i]*dt-dby_be_H[i])/(2*mi[i]*mi_0)))/(1+((sih[i]*dt+0.5*lor_gam_H[i]+dby_be_H[i])/(2*mi[i]*mi_0)));
        }  
  
		if (is_PMC[i])
		{
			dmit[i]=0.0;
			dsiht[i]=-1.0;
		}
	}
}

// pomocnicza funkcja dla obliczenia pola na UPML - zwraca wartosci c1,c2 itd.
double cstx(double* cstnt,int index,int size)
{
    if (index+boundx>=size) return cstnt[size-index-1];
    if (index>boundx) return cstnt[boundx];
    return cstnt[index];
}
double csty(double* cstnt,int index,int size)
{
    if (index+boundy>=size) return cstnt[size-index-1];
    if (index>boundy) return cstnt[boundy];
    return cstnt[index];
}
double cstz(double* cstnt,int index,int size)
{
    if (index+boundz>=size) return cstnt[size-index-1];
    if (index>boundz) return cstnt[boundz];
    return cstnt[index];
}
// funkcje rozkladu sigma i kappa w warunkcach brzegowych - rozklad wielomianowy
// location to pozycja waznej dla nas skladowej i,j lub k razy 2
// sigma powinno dla max.location wynosic 0 dla zera zas wartosc maksymalna

typ_prec sigma(typ_prec location,typ_prec boundr) { return sigma_max*pow((typ_prec)(1.0-location/boundr),wykladnik); }
typ_prec kappa(typ_prec location,typ_prec boundr) { return 1+((kappa_max-1)*pow((typ_prec)(1.0-location/boundr),wykladnik)); }

// funkcja wypelniajace dodatkowe tablice materialowe na brzegach - czyli stale UPML
// por. Taflove, str.
void fill_ABC()
{
    typ_prec temp=eps_0*eps[bound_mat];

    for (i=0;i<boundx;i++)
    {
        // pole elektryczne
        c1x_E[i]=(2*temp*kappa(i,boundx)-sigma(i,boundx)*dt)/(2*temp*kappa(i,boundx)+sigma(i,boundx)*dt);
        c2x_E[i]=(2*temp*dt/(dr*(2*temp*kappa(i,boundx)+sigma(i,boundx)*dt)))/H_scale;

        c3x_E[i]=(2*temp*kappa(i,boundx)-sigma(i,boundx)*dt)/(2*temp*kappa(i,boundx)+sigma(i,boundx)*dt);
        c4x_E[i]=1/(temp*(2*temp*kappa(i,boundx)+sigma(i,boundx)*dt));

        c5x_E[i]=2*temp*kappa(i+0.5,boundx)+sigma(i+0.5,boundx)*dt;
        c6x_E[i]=2*temp*kappa(i+0.5,boundx)-sigma(i+0.5,boundx)*dt;

        // pole magnetyczne
        c1x_H[i]=(2*temp*kappa(i+0.5,boundx)-sigma(i+0.5,boundx)*dt)/(2*temp*kappa(i+0.5,boundx)+sigma(i+0.5,boundx)*dt);
        c2x_H[i]=(2*temp*dt/(dr*(2*temp*kappa(i+0.5,boundx)+sigma(i+0.5,boundx)*dt)))*H_scale;

        c3x_H[i]=(2*temp*kappa(i+0.5,boundx)-sigma(i+0.5,boundx)*dt)/(2*temp*kappa(i+0.5,boundx)+sigma(i+0.5,boundx)*dt);
        c4x_H[i]=1/(mi_0*mi[bound_mat]*(2*temp*kappa(i+0.5,boundx)+sigma(i+0.5,boundx)*dt));

        c5x_H[i]=2*temp*kappa(i+1,boundx)+sigma(i+1,boundx)*dt;
        c6x_H[i]=2*temp*kappa(i+1,boundx)-sigma(i+1,boundx)*dt;
    }

    // pole elektryczne
    c1x_E[boundx]=1;
    c2x_E[boundx]=(dt/dr)/H_scale;

    c3x_E[boundx]=1;
    c4x_E[boundx]=1/(2*temp*temp);

    c5x_E[boundx]=2*temp;
    c6x_E[boundx]=2*temp;

    // pole magnetyczne
    c1x_H[boundx]=1;
    c2x_H[boundx]=(dt/dr)*H_scale;

    c3x_H[boundx]=1;
    c4x_H[boundx]=1/(2*mi_0*mi[bound_mat]*temp);

    c5x_H[boundx]=2*temp;
    c6x_H[boundx]=2*temp;

    for (i=0;i<boundy;i++)
    {
        // pole elektryczne
        c1y_E[i]=(2*temp*kappa(i,boundy)-sigma(i,boundy)*dt)/(2*temp*kappa(i,boundy)+sigma(i,boundy)*dt);
        c2y_E[i]=(2*temp*dt/(dr*(2*temp*kappa(i,boundy)+sigma(i,boundy)*dt)))/H_scale;

        c3y_E[i]=(2*temp*kappa(i,boundy)-sigma(i,boundy)*dt)/(2*temp*kappa(i,boundy)+sigma(i,boundy)*dt);
        c4y_E[i]=1/(temp*(2*temp*kappa(i,boundy)+sigma(i,boundy)*dt));

        c5y_E[i]=2*temp*kappa(i+0.5,boundy)+sigma(i+0.5,boundy)*dt;
        c6y_E[i]=2*temp*kappa(i+0.5,boundy)-sigma(i+0.5,boundy)*dt;

        // pole magnetyczne
        c1y_H[i]=(2*temp*kappa(i+0.5,boundy)-sigma(i+0.5,boundy)*dt)/(2*temp*kappa(i+0.5,boundy)+sigma(i+0.5,boundy)*dt);
        c2y_H[i]=(2*temp*dt/(dr*(2*temp*kappa(i+0.5,boundy)+sigma(i+0.5,boundy)*dt)))*H_scale;

        c3y_H[i]=(2*temp*kappa(i+0.5,boundy)-sigma(i+0.5,boundy)*dt)/(2*temp*kappa(i+0.5,boundy)+sigma(i+0.5,boundy)*dt);
        c4y_H[i]=1/(mi_0*mi[bound_mat]*(2*temp*kappa(i+0.5,boundy)+sigma(i+0.5,boundy)*dt));

        c5y_H[i]=2*temp*kappa(i+1,boundy)+sigma(i+1,boundy)*dt;
        c6y_H[i]=2*temp*kappa(i+1,boundy)-sigma(i+1,boundy)*dt;
    }

    // pole elektryczne
    c1y_E[boundy]=1;
    c2y_E[boundy]=(dt/dr)/H_scale;

    c3y_E[boundy]=1;
    c4y_E[boundy]=1/(2*temp*temp);

    c5y_E[boundy]=2*temp;
    c6y_E[boundy]=2*temp;

    // pole magnetyczne
    c1y_H[boundy]=1;
    c2y_H[boundy]=(dt/dr)*H_scale;

    c3y_H[boundy]=1;
    c4y_H[boundy]=1/(2*mi_0*mi[bound_mat]*temp);

    c5y_H[boundy]=2*temp;
    c6y_H[boundy]=2*temp;

    for (i=0;i<boundz;i++)
    {
        // pole elektryczne
        c1z_E[i]=(2*temp*kappa(i,boundz)-sigma(i,boundz)*dt)/(2*temp*kappa(i,boundz)+sigma(i,boundz)*dt);
        c2z_E[i]=(2*temp*dt/(dr*(2*temp*kappa(i,boundz)+sigma(i,boundz)*dt)))/H_scale;

        c3z_E[i]=(2*temp*kappa(i,boundz)-sigma(i,boundz)*dt)/(2*temp*kappa(i,boundz)+sigma(i,boundz)*dt);
        c4z_E[i]=1/(temp*(2*temp*kappa(i,boundz)+sigma(i,boundz)*dt));

        c5z_E[i]=2*temp*kappa(i+0.5,boundz)+sigma(i+0.5,boundz)*dt;
        c6z_E[i]=2*temp*kappa(i+0.5,boundz)-sigma(i+0.5,boundz)*dt;

        // pole magnetyczne
        c1z_H[i]=(2*temp*kappa(i+0.5,boundz)-sigma(i+0.5,boundz)*dt)/(2*temp*kappa(i+0.5,boundz)+sigma(i+0.5,boundz)*dt);
        c2z_H[i]=(2*temp*dt/(dr*(2*temp*kappa(i+0.5,boundz)+sigma(i+0.5,boundz)*dt)))*H_scale;

        c3z_H[i]=(2*temp*kappa(i+0.5,boundz)-sigma(i+0.5,boundz)*dt)/(2*temp*kappa(i+0.5,boundz)+sigma(i+0.5,boundz)*dt);
        c4z_H[i]=1/(mi_0*mi[bound_mat]*(2*temp*kappa(i+0.5,boundz)+sigma(i+0.5,boundz)*dt));

        c5z_H[i]=2*temp*kappa(i+1,boundz)+sigma(i+1,boundz)*dt;
        c6z_H[i]=2*temp*kappa(i+1,boundz)-sigma(i+1,boundz)*dt;
    }

    // pole elektryczne
    c1z_E[boundz]=1;
    c2z_E[boundz]=(dt/dr)/H_scale;

    c3z_E[boundz]=1;
    c4z_E[boundz]=1/(2*temp*temp);

    c5z_E[boundz]=2*temp;
    c6z_E[boundz]=2*temp;

    // pole magnetyczne
    c1z_H[boundz]=1;
    c2z_H[boundz]=(dt/dr)*H_scale;

    c3z_H[boundz]=1;
    c4z_H[boundz]=1/(2*mi_0*mi[bound_mat]*temp);

    c5z_H[boundz]=2*temp;
    c6z_H[boundz]=2*temp;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                   //
//                      podstawowe funkcje FDTD, uaktualnienie pol - pojedynczy krok czasowy symulacji               //
//                                                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////
//  wnetrze obszaru symulowanego - "zwykle" FDTD   //
/////////////////////////////////////////////////////

// Funkcje dla struktury 2D
// funkcja dla pola magnetycznego we wnetrzu obszaru symulowanego
void update_inter_Hx()
{
    ////////
    // Hx //
    ////////

    for(i=boundx;i<xk-1;i++)
    {
        for(j=boundy;j<yk;j++)
        {
            temp=dmit[matHx[i-boundx][j-boundy]];
	        temp2=dsiht[matHx[i-boundx][j-boundy]];
            for(k=boundz;k<boundz+bndz+disPML;k++) Hx[i][j][k]+=dmit[bound_mat]*(Ey[i+1][j][k+1]-Ey[i+1][j][k]-Ez[i+1][j+1][k]+Ez[i+1][j][k]);
            for(k=boundz+bndz+disPML;k<zk-bndz-disPML;k++) Hx[i][j][k]=temp2*Hx[i][j][k]+temp*(Ey[i+1][j][k+1]-Ey[i+1][j][k]-Ez[i+1][j+1][k]+Ez[i+1][j][k]);
            for(k=zk-bndz-disPML;k<zk;k++) Hx[i][j][k]+=dmit[bound_mat]*(Ey[i+1][j][k+1]-Ey[i+1][j][k]-Ez[i+1][j+1][k]+Ez[i+1][j][k]);
        }
    }
}
void update_inter_Hy()
{
    ////////
    // Hy //
    ////////

    for(i=boundx;i<xk;i++)
    {
        for(j=boundy;j<yk-1;j++)
        {
            temp=dmit[matHy[i-boundx][j-boundy]];
	    temp2=dsiht[matHy[i-boundx][j-boundy]];
            for(k=boundz;k<boundz+bndz+disPML;k++) Hy[i][j][k]+=dmit[bound_mat]*(Ez[i+1][j+1][k]-Ez[i][j+1][k]-Ex[i][j+1][k+1]+Ex[i][j+1][k]);
            for(k=boundz+bndz+disPML;k<zk-bndz-disPML;k++) Hy[i][j][k]=temp2*Hy[i][j][k]+temp*(Ez[i+1][j+1][k]-Ez[i][j+1][k]-Ex[i][j+1][k+1]+Ex[i][j+1][k]);
            for(k=zk-bndz-disPML;k<zk;k++) Hy[i][j][k]+=dmit[bound_mat]*(Ez[i+1][j+1][k]-Ez[i][j+1][k]-Ex[i][j+1][k+1]+Ex[i][j+1][k]);
        }
    }

}
void update_inter_Hz()
{

    ////////
    // Hz //
    ////////

    for(i=boundx;i<xk;i++)
    {
        for(j=boundy;j<yk;j++)
        {
            temp=dmit[matHz[i-boundx][j-boundy]];
	    temp2=dsiht[matHz[i-boundx][j-boundy]];
            for(k=boundz;k<boundz+bndz+disPML;k++) Hz[i][j][k]+=dmit[bound_mat]*(Ex[i][j+1][k+1]-Ex[i][j][k+1]-Ey[i+1][j][k+1]+Ey[i][j][k+1]);
            for(k=boundz+bndz+disPML;k<zk-1-bndz-disPML;k++) Hz[i][j][k]=temp2*Hz[i][j][k]+temp*(Ex[i][j+1][k+1]-Ex[i][j][k+1]-Ey[i+1][j][k+1]+Ey[i][j][k+1]);
            for(k=zk-1-bndz-disPML;k<zk-1;k++) Hz[i][j][k]+=dmit[bound_mat]*(Ex[i][j+1][k+1]-Ex[i][j][k+1]-Ey[i+1][j][k+1]+Ey[i][j][k+1]);
        }
    }
}

// funkcja dla pola elektrycznego we wnetrzu obszaru symulowanego
void update_inter_Ex()
{

    ////////
    // Ex //
    ////////

    for(i=boundx;i<xk;i++)
    {
        for(j=boundy;j<=yk;j++)
        {
            temp=depst[matEx[i-boundx][j-boundy]];
	    temp2=dsigt[matEx[i-boundx][j-boundy]];
            for(k=boundz;k<boundz+bndz+disPML;k++) Ex[i][j][k]+=depst[bound_mat]*(Hz[i][j][k-1]-Hz[i][j-1][k-1]-Hy[i][j-1][k]+Hy[i][j-1][k-1]);
            for(k=boundz+bndz+disPML;k<=zk-bndz-disPML;k++) Ex[i][j][k]=temp2*Ex[i][j][k]+temp*(Hz[i][j][k-1]-Hz[i][j-1][k-1]-Hy[i][j-1][k]+Hy[i][j-1][k-1]);
            for(k=zk-(bndz-1)-disPML;k<=zk;k++) Ex[i][j][k]+=depst[bound_mat]*(Hz[i][j][k-1]-Hz[i][j-1][k-1]-Hy[i][j-1][k]+Hy[i][j-1][k-1]);
        }
    }

}

// funkcja dla pola elektrycznego we wnetrzu obszaru symulowanego
void update_inter_Ey()
{
    ////////
    // Ey //
    ////////

    for(i=boundx;i<=xk;i++)
    {
        for(j=boundy;j<yk;j++)
        {
            temp=depst[matEy[i-boundx][j-boundy]];
	    temp2=dsigt[matEy[i-boundx][j-boundy]];
            for(k=boundz;k<boundz+bndz+disPML;k++) Ey[i][j][k]+=depst[bound_mat]*(Hx[i-1][j][k]-Hx[i-1][j][k-1]-Hz[i][j][k-1]+Hz[i-1][j][k-1]);
            for(k=boundz+bndz+disPML;k<=zk-bndz-disPML;k++) Ey[i][j][k]=temp2*Ey[i][j][k]+temp*(Hx[i-1][j][k]-Hx[i-1][j][k-1]-Hz[i][j][k-1]+Hz[i-1][j][k-1]);
            for(k=zk-(bndz-1)-disPML;k<=zk;k++) Ey[i][j][k]+=depst[bound_mat]*(Hx[i-1][j][k]-Hx[i-1][j][k-1]-Hz[i][j][k-1]+Hz[i-1][j][k-1]);
        }
    }

}

// funkcja dla pola elektrycznego we wnetrzu obszaru symulowanego
void update_inter_Ez()
{
    ////////
    // Ez //
    ////////

    for(i=boundx;i<=xk;i++)
    {
        for(j=boundy;j<=yk;j++)
        {
            temp=depst[matEz[i-boundx][j-boundy]];
	    temp2=dsigt[matEz[i-boundx][j-boundy]];
            for(k=boundz;k<boundz+bndz+disPML;k++) Ez[i][j][k]+=depst[bound_mat]*(Hy[i][j-1][k]-Hy[i-1][j-1][k]-Hx[i-1][j][k]+Hx[i-1][j-1][k]);
            for(k=boundz+bndz+disPML;k<zk-bndz-disPML;k++) Ez[i][j][k]=temp2*Ez[i][j][k]+temp*(Hy[i][j-1][k]-Hy[i-1][j-1][k]-Hx[i-1][j][k]+Hx[i-1][j-1][k]);
            for(k=zk-bndz-disPML;k<zk;k++) Ez[i][j][k]+=depst[bound_mat]*(Hy[i][j-1][k]-Hy[i-1][j-1][k]-Hx[i-1][j][k]+Hx[i-1][j-1][k]);

        }
    }
}

//////////////////////////////////////////////////////////////////////
// funkcje dla struktury 3D
// funkcja dla pola magnetycznego we wnetrzu obszaru symulowanego
void update_inter_H3x()
{
    ////////
    // Hx //
    ////////

    for(i=boundx;i<xk-1;i++)
    {
        for(j=boundy;j<yk;j++)
        {
            for(k=boundz;k<zk;k++) Hx[i][j][k]=dsiht[matHx3[i-boundx][j-boundy][k-boundz]]*Hx[i][j][k]+dmit[matHx3[i-boundx][j-boundy][k-boundz]]*(Ey[i+1][j][k+1]-Ey[i+1][j][k]-Ez[i+1][j+1][k]+Ez[i+1][j][k]);
        }
    }
}

void update_inter_H3y()
{
    ////////
    // Hy //
    ////////

    for(i=boundx;i<xk;i++)
    {
        for(j=boundy;j<yk-1;j++)
        {
            for(k=boundz;k<zk;k++) Hy[i][j][k]=dsiht[matHy3[i-boundx][j-boundy][k-boundz]]*Hy[i][j][k]+dmit[matHy3[i-boundx][j-boundy][k-boundz]]*(Ez[i+1][j+1][k]-Ez[i][j+1][k]-Ex[i][j+1][k+1]+Ex[i][j+1][k]);
        }
    }
}
void update_inter_H3z()
{

    ////////
    // Hz //
    ////////

    for(i=boundx;i<xk;i++)
    {
        for(j=boundy;j<yk;j++)
        {
            for(k=boundz;k<zk-1;k++) Hz[i][j][k]=dsiht[matHz3[i-boundx][j-boundy][k-boundz]]*Hz[i][j][k]+dmit[matHz3[i-boundx][j-boundy][k-boundz]]*(Ex[i][j+1][k+1]-Ex[i][j][k+1]-Ey[i+1][j][k+1]+Ey[i][j][k+1]);
        }
    }
}

// funkcja dla pola elektrycznego we wnetrzu obszaru symulowanego
void update_inter_E3x()
{

    ////////
    // Ex //
    ////////

    for(i=boundx;i<xk;i++)
    {
        for(j=boundy;j<=yk;j++)
        {
            for(k=boundz;k<=zk;k++) Ex[i][j][k]=dsigt[matEx3[i-boundx][j-boundy][k-boundz]]*Ex[i][j][k]+depst[matEx3[i-boundx][j-boundy][k-boundz]]*(Hz[i][j][k-1]-Hz[i][j-1][k-1]-Hy[i][j-1][k]+Hy[i][j-1][k-1]);
        }
    }
}
void update_inter_E3y()
{
    ////////
    // Ey //
    ////////

    for(i=boundx;i<=xk;i++)
    {
        for(j=boundy;j<yk;j++)
        {
            for(k=boundz;k<=zk;k++) Ey[i][j][k]=dsigt[matEy3[i-boundx][j-boundy][k-boundz]]*Ey[i][j][k]+depst[matEy3[i-boundx][j-boundy][k-boundz]]*(Hx[i-1][j][k]-Hx[i-1][j][k-1]-Hz[i][j][k-1]+Hz[i-1][j][k-1]);
        }
    }
}
void update_inter_E3z()
{

    ////////
    // Ez //
    ////////

    for(i=boundx;i<=xk;i++)
    {
        for(j=boundy;j<=yk;j++)
        {
            for(k=boundz;k<zk;k++) Ez[i][j][k]=dsigt[matEz3[i-boundx][j-boundy][k-boundz]]*Ez[i][j][k]+depst[matEz3[i-boundx][j-boundy][k-boundz]]*(Hy[i][j-1][k]-Hy[i-1][j-1][k]-Hx[i-1][j][k]+Hx[i-1][j-1][k]);
        }
    }

}


/////////////////////////////////////////////////////////////////////
// Funkcje zubo¿one dla obliczeñ 2D 

// Funkcje dla struktury 2D
// funkcja dla pola magnetycznego we wnetrzu obszaru symulowanego
void update_inter_Hx_null_Ey()
{
    ////////
    // Hx //
    ////////

    for(i=boundx;i<xk-1;i++)
    {
        for(j=boundy;j<yk;j++)
        {
            temp=dmit[matHx[i-boundx][j-boundy]];
	    temp2=dsiht[matHx[i-boundx][j-boundy]];
            for(k=boundz;k<boundz+bndz+disPML;k++) Hx[i][j][k]+=dmit[bound_mat]*(-Ez[i+1][j+1][k]+Ez[i+1][j][k]);
            for(k=boundz+bndz+disPML;k<zk-bndz-disPML;k++) Hx[i][j][k]=temp2*Hx[i][j][k]+temp*(-Ez[i+1][j+1][k]+Ez[i+1][j][k]);
            for(k=zk-bndz-disPML;k<zk;k++) Hx[i][j][k]+=dmit[bound_mat]*(-Ez[i+1][j+1][k]+Ez[i+1][j][k]);
        }
    }
}

void update_inter_Hx_null_Ez()
{
    ////////
    // Hx //
    ////////

    for(i=boundx;i<xk-1;i++)
    {
        for(j=boundy;j<yk;j++)
        {
            temp=dmit[matHx[i-boundx][j-boundy]];
	    temp2=dsiht[matHx[i-boundx][j-boundy]];
            for(k=boundz;k<boundz+bndz+disPML;k++) Hx[i][j][k]+=dmit[bound_mat]*(Ey[i+1][j][k+1]-Ey[i+1][j][k]);
            for(k=boundz+bndz+disPML;k<zk-bndz-disPML;k++) Hx[i][j][k]=temp2*Hx[i][j][k]+temp*(Ey[i+1][j][k+1]-Ey[i+1][j][k]);
            for(k=zk-bndz-disPML;k<zk;k++) Hx[i][j][k]+=dmit[bound_mat]*(Ey[i+1][j][k+1]-Ey[i+1][j][k]);
        }
    }
}

void update_inter_Hy_null_Ex()
{
    ////////
    // Hy //
    ////////

    for(i=boundx;i<xk;i++)
    {
        for(j=boundy;j<yk-1;j++)
        {
            temp=dmit[matHy[i-boundx][j-boundy]];
	    temp2=dsiht[matHy[i-boundx][j-boundy]];
            for(k=boundz;k<boundz+bndz+disPML;k++) Hy[i][j][k]+=dmit[bound_mat]*(Ez[i+1][j+1][k]-Ez[i][j+1][k]);
            for(k=boundz+bndz+disPML;k<zk-bndz-disPML;k++) Hy[i][j][k]=temp2*Hy[i][j][k]+temp*(Ez[i+1][j+1][k]-Ez[i][j+1][k]);
            for(k=zk-bndz-disPML;k<zk;k++) Hy[i][j][k]+=dmit[bound_mat]*(Ez[i+1][j+1][k]-Ez[i][j+1][k]);
        }
    }

}

void update_inter_Hy_null_Ez()
{
    ////////
    // Hy //
    ////////

    for(i=boundx;i<xk;i++)
    {
        for(j=boundy;j<yk-1;j++)
        {
            temp=dmit[matHy[i-boundx][j-boundy]];
	    temp2=dsiht[matHy[i-boundx][j-boundy]];
            for(k=boundz;k<boundz+bndz+disPML;k++) Hy[i][j][k]+=dmit[bound_mat]*(-Ex[i][j+1][k+1]+Ex[i][j+1][k]);
            for(k=boundz+bndz+disPML;k<zk-bndz-disPML;k++) Hy[i][j][k]=temp2*Hy[i][j][k]+temp*(-Ex[i][j+1][k+1]+Ex[i][j+1][k]);
            for(k=zk-bndz-disPML;k<zk;k++) Hy[i][j][k]+=dmit[bound_mat]*(-Ex[i][j+1][k+1]+Ex[i][j+1][k]);
        }
    }

}

void update_inter_Hz_null_Ey()
{

    ////////
    // Hz //
    ////////

    for(i=boundx;i<xk;i++)
    {
        for(j=boundy;j<yk;j++)
        {
            temp=dmit[matHz[i-boundx][j-boundy]];
	    temp2=dsiht[matHz[i-boundx][j-boundy]];
            for(k=boundz;k<boundz+bndz+disPML;k++) Hz[i][j][k]+=dmit[bound_mat]*(Ex[i][j+1][k+1]-Ex[i][j][k+1]);
            for(k=boundz+bndz+disPML;k<zk-1-bndz-disPML;k++) Hz[i][j][k]=temp2*Hz[i][j][k]+temp*(Ex[i][j+1][k+1]-Ex[i][j][k+1]);
            for(k=zk-1-bndz-disPML;k<zk-1;k++) Hz[i][j][k]+=dmit[bound_mat]*(Ex[i][j+1][k+1]-Ex[i][j][k+1]);
        }
    }
}

void update_inter_Hz_null_Ex()
{

    ////////
    // Hz //
    ////////

    for(i=boundx;i<xk;i++)
    {
        for(j=boundy;j<yk;j++)
        {
            temp=dmit[matHz[i-boundx][j-boundy]];
	    temp2=dsiht[matHz[i-boundx][j-boundy]];
            for(k=boundz;k<boundz+bndz+disPML;k++) Hz[i][j][k]+=dmit[bound_mat]*(-Ey[i+1][j][k+1]+Ey[i][j][k+1]);
            for(k=boundz+bndz+disPML;k<zk-1-bndz-disPML;k++) Hz[i][j][k]=temp2*Hz[i][j][k]+temp*(-Ey[i+1][j][k+1]+Ey[i][j][k+1]);
            for(k=zk-1-bndz-disPML;k<zk-1;k++) Hz[i][j][k]+=dmit[bound_mat]*(-Ey[i+1][j][k+1]+Ey[i][j][k+1]);
        }
    }
}
// funkcja dla pola elektrycznego we wnetrzu obszaru symulowanego
void update_inter_Ex_null_Hy()
{

    ////////
    // Ex //
    ////////

    for(i=boundx;i<xk;i++)
    {
        for(j=boundy;j<=yk;j++)
        {
            temp=depst[matEx[i-boundx][j-boundy]];
	    temp2=dsigt[matEx[i-boundx][j-boundy]];
            for(k=boundz;k<boundz+bndz+disPML;k++) Ex[i][j][k]+=depst[bound_mat]*(Hz[i][j][k-1]-Hz[i][j-1][k-1]);
            for(k=boundz+bndz+disPML;k<=zk-bndz-disPML;k++) Ex[i][j][k]=temp2*Ex[i][j][k]+temp*(Hz[i][j][k-1]-Hz[i][j-1][k-1]);
            for(k=zk-(bndz-1)-disPML;k<=zk;k++) Ex[i][j][k]+=depst[bound_mat]*(Hz[i][j][k-1]-Hz[i][j-1][k-1]);
        }
    }

}
void update_inter_Ex_null_Hz()
{

    ////////
    // Ex //
    ////////

    for(i=boundx;i<xk;i++)
    {
        for(j=boundy;j<=yk;j++)
        {
            temp=depst[matEx[i-boundx][j-boundy]];
	    temp2=dsigt[matEx[i-boundx][j-boundy]];
            for(k=boundz;k<boundz+bndz+disPML;k++) Ex[i][j][k]+=depst[bound_mat]*(-Hy[i][j-1][k]+Hy[i][j-1][k-1]);
            for(k=boundz+bndz+disPML;k<=zk-bndz-disPML;k++) Ex[i][j][k]=temp2*Ex[i][j][k]+temp*(-Hy[i][j-1][k]+Hy[i][j-1][k-1]);
            for(k=zk-(bndz-1)-disPML;k<=zk;k++) Ex[i][j][k]+=depst[bound_mat]*(-Hy[i][j-1][k]+Hy[i][j-1][k-1]);
        }
    }

}
// funkcja dla pola elektrycznego we wnetrzu obszaru symulowanego
void update_inter_Ey_null_Hz()
{
    ////////
    // Ey //
    ////////

    for(i=boundx;i<=xk;i++)
    {
        for(j=boundy;j<yk;j++)
        {
            temp=depst[matEy[i-boundx][j-boundy]];
	    temp2=dsigt[matEy[i-boundx][j-boundy]];
            for(k=boundz;k<boundz+bndz+disPML;k++) Ey[i][j][k]+=depst[bound_mat]*(Hx[i-1][j][k]-Hx[i-1][j][k-1]);
            for(k=boundz+bndz+disPML;k<=zk-bndz-disPML;k++) Ey[i][j][k]=temp2*Ey[i][j][k]+temp*(Hx[i-1][j][k]-Hx[i-1][j][k-1]);
            for(k=zk-(bndz-1)-disPML;k<=zk;k++) Ey[i][j][k]+=depst[bound_mat]*(Hx[i-1][j][k]-Hx[i-1][j][k-1]);
        }
    }

}
void update_inter_Ey_null_Hx()
{
    ////////
    // Ey //
    ////////

    for(i=boundx;i<=xk;i++)
    {
        for(j=boundy;j<yk;j++)
        {
            temp=depst[matEy[i-boundx][j-boundy]];
	    temp2=dsigt[matEy[i-boundx][j-boundy]];
            for(k=boundz;k<boundz+bndz+disPML;k++) Ey[i][j][k]+=depst[bound_mat]*(-Hz[i][j][k-1]+Hz[i-1][j][k-1]);
            for(k=boundz+bndz+disPML;k<=zk-bndz-disPML;k++) Ey[i][j][k]=temp2*Ey[i][j][k]+temp*(-Hz[i][j][k-1]+Hz[i-1][j][k-1]);
            for(k=zk-(bndz-1)-disPML;k<=zk;k++) Ey[i][j][k]+=depst[bound_mat]*(-Hz[i][j][k-1]+Hz[i-1][j][k-1]);
        }
    }

}
// funkcja dla pola elektrycznego we wnetrzu obszaru symulowanego
void update_inter_Ez_null_Hx()
{
    ////////
    // Ez //
    ////////

    for(i=boundx;i<=xk;i++)
    {
        for(j=boundy;j<=yk;j++)
        {
            temp=depst[matEz[i-boundx][j-boundy]];
	    temp2=dsigt[matEz[i-boundx][j-boundy]];
            for(k=boundz;k<boundz+bndz+disPML;k++) Ez[i][j][k]+=depst[bound_mat]*(Hy[i][j-1][k]-Hy[i-1][j-1][k]);
            for(k=boundz+bndz+disPML;k<zk-bndz-disPML;k++) Ez[i][j][k]=temp2*Ez[i][j][k]+temp*(Hy[i][j-1][k]-Hy[i-1][j-1][k]);
            for(k=zk-bndz-disPML;k<zk;k++) Ez[i][j][k]+=depst[bound_mat]*(Hy[i][j-1][k]-Hy[i-1][j-1][k]);

        }
    }
}

void update_inter_Ez_null_Hy()
{
    ////////
    // Ez //
    ////////

    for(i=boundx;i<=xk;i++)
    {
        for(j=boundy;j<=yk;j++)
        {
            temp=depst[matEz[i-boundx][j-boundy]];
	    temp2=dsigt[matEz[i-boundx][j-boundy]];
            for(k=boundz;k<boundz+bndz+disPML;k++) Ez[i][j][k]+=depst[bound_mat]*(-Hx[i-1][j][k]+Hx[i-1][j-1][k]);
            for(k=boundz+bndz+disPML;k<zk-bndz-disPML;k++) Ez[i][j][k]=temp2*Ez[i][j][k]+temp*(-Hx[i-1][j][k]+Hx[i-1][j-1][k]);
            for(k=zk-bndz-disPML;k<zk;k++) Ez[i][j][k]+=depst[bound_mat]*(-Hx[i-1][j][k]+Hx[i-1][j-1][k]);

        }
    }
}

//////////////////////////////////////////////////////////////////////
// funkcje zubo¿one c.d.
// funkcje dla struktury 3D
// funkcja dla pola magnetycznego we wnetrzu obszaru symulowanego
void update_inter_H3x_null_Ez()
{
    ////////
    // Hx //
    ////////

    for(i=boundx;i<xk-1;i++)
    {
        for(j=boundy;j<yk;j++)
        {
            for(k=boundz;k<zk;k++) Hx[i][j][k]=dsiht[matHx3[i-boundx][j-boundy][k-boundz]]*Hx[i][j][k]+dmit[matHx3[i-boundx][j-boundy][k-boundz]]*(Ey[i+1][j][k+1]-Ey[i+1][j][k]);
        }
    }
}
void update_inter_H3x_null_Ey()
{
    ////////
    // Hx //
    ////////

    for(i=boundx;i<xk-1;i++)
    {
        for(j=boundy;j<yk;j++)
        {
            for(k=boundz;k<zk;k++) Hx[i][j][k]=dsiht[matHx3[i-boundx][j-boundy][k-boundz]]*Hx[i][j][k]+dmit[matHx3[i-boundx][j-boundy][k-boundz]]*(-Ez[i+1][j+1][k]+Ez[i+1][j][k]);
        }
    }
}
void update_inter_H3y_null_Ex()
{
    ////////
    // Hy //
    ////////

    for(i=boundx;i<xk;i++)
    {
        for(j=boundy;j<yk-1;j++)
        {
            for(k=boundz;k<zk;k++) Hy[i][j][k]=dsiht[matHy3[i-boundx][j-boundy][k-boundz]]*Hy[i][j][k]+dmit[matHy3[i-boundx][j-boundy][k-boundz]]*(Ez[i+1][j+1][k]-Ez[i][j+1][k]);
        }
    }
}
void update_inter_H3y_null_Ez()
{
    ////////
    // Hy //
    ////////

    for(i=boundx;i<xk;i++)
    {
        for(j=boundy;j<yk-1;j++)
        {
            for(k=boundz;k<zk;k++) Hy[i][j][k]=dsiht[matHy3[i-boundx][j-boundy][k-boundz]]*Hy[i][j][k]+dmit[matHy3[i-boundx][j-boundy][k-boundz]]*(-Ex[i][j+1][k+1]+Ex[i][j+1][k]);
        }
    }
}

void update_inter_H3z_null_Ey()
{

    ////////
    // Hz //
    ////////

    for(i=boundx;i<xk;i++)
    {
        for(j=boundy;j<yk;j++)
        {
            for(k=boundz;k<zk-1;k++) Hz[i][j][k]=dsiht[matHz3[i-boundx][j-boundy][k-boundz]]*Hz[i][j][k]+dmit[matHz3[i-boundx][j-boundy][k-boundz]]*(Ex[i][j+1][k+1]-Ex[i][j][k+1]);
        }
    }
}

void update_inter_H3z_null_Ex()
{

    ////////
    // Hz //
    ////////

    for(i=boundx;i<xk;i++)
    {
        for(j=boundy;j<yk;j++)
        {
            for(k=boundz;k<zk-1;k++) Hz[i][j][k]=dsiht[matHz3[i-boundx][j-boundy][k-boundz]]*Hz[i][j][k]+dmit[matHz3[i-boundx][j-boundy][k-boundz]]*(-Ey[i+1][j][k+1]+Ey[i][j][k+1]);
        }
    }
}

// funkcja dla pola elektrycznego we wnetrzu obszaru symulowanego
void update_inter_E3x_null_Hy()
{

    ////////
    // Ex //
    ////////

    for(i=boundx;i<xk;i++)
    {
        for(j=boundy;j<=yk;j++)
        {
            for(k=boundz;k<=zk;k++) Ex[i][j][k]=dsigt[matEx3[i-boundx][j-boundy][k-boundz]]*Ex[i][j][k]+depst[matEx3[i-boundx][j-boundy][k-boundz]]*(Hz[i][j][k-1]-Hz[i][j-1][k-1]);
        }
    }
}
void update_inter_E3x_null_Hz()
{

    ////////
    // Ex //
    ////////

    for(i=boundx;i<xk;i++)
    {
        for(j=boundy;j<=yk;j++)
        {
            for(k=boundz;k<=zk;k++) Ex[i][j][k]=dsigt[matEx3[i-boundx][j-boundy][k-boundz]]*Ex[i][j][k]+depst[matEx3[i-boundx][j-boundy][k-boundz]]*(-Hy[i][j-1][k]+Hy[i][j-1][k-1]);
        }
    }
}
void update_inter_E3y_null_Hz()
{
    ////////
    // Ey //
    ////////

    for(i=boundx;i<=xk;i++)
    {
        for(j=boundy;j<yk;j++)
        {
            for(k=boundz;k<=zk;k++) Ey[i][j][k]=dsigt[matEy3[i-boundx][j-boundy][k-boundz]]*Ey[i][j][k]+depst[matEy3[i-boundx][j-boundy][k-boundz]]*(Hx[i-1][j][k]-Hx[i-1][j][k-1]);
        }
    }
}
void update_inter_E3y_null_Hx()
{
    ////////
    // Ey //
    ////////

    for(i=boundx;i<=xk;i++)
    {
        for(j=boundy;j<yk;j++)
        {
            for(k=boundz;k<=zk;k++) Ey[i][j][k]=dsigt[matEy3[i-boundx][j-boundy][k-boundz]]*Ey[i][j][k]+depst[matEy3[i-boundx][j-boundy][k-boundz]]*(-Hz[i][j][k-1]+Hz[i-1][j][k-1]);
        }
    }
}
void update_inter_E3z_null_Hx()
{

    ////////
    // Ez //
    ////////

    for(i=boundx;i<=xk;i++)
    {
        for(j=boundy;j<=yk;j++)
        {
            for(k=boundz;k<zk;k++) Ez[i][j][k]=dsigt[matEz3[i-boundx][j-boundy][k-boundz]]*Ez[i][j][k]+depst[matEz3[i-boundx][j-boundy][k-boundz]]*(Hy[i][j-1][k]-Hy[i-1][j-1][k]);
        }
    }

}

void update_inter_E3z_null_Hy()
{

    ////////
    // Ez //
    ////////

    for(i=boundx;i<=xk;i++)
    {
        for(j=boundy;j<=yk;j++)
        {
            for(k=boundz;k<zk;k++) Ez[i][j][k]=dsigt[matEz3[i-boundx][j-boundy][k-boundz]]*Ez[i][j][k]+depst[matEz3[i-boundx][j-boundy][k-boundz]]*(-Hx[i-1][j][k]+Hx[i-1][j-1][k]);
        }
    }

}


/////////////////////////////////////////////////////////////////
//           brzeg obszaru symulowanego - FDTD+UPML            //
/////////////////////////////////////////////////////////////////
//
//     /z
//    /
//   -----
//   |   x  /|||||||||||||
//   |y    /            /|
//        /    bottom  / |
//       /            /  |
//      |||||||||||||| right
//      |    back    |  /
//      |            | /
//      ||||||||||||||/
//
////////////////////////////////////////////////////////////////


void update_bot_Hx()
{

    typ_pola temp;

    ////////
    // Hx //
    ////////

    c3=c3z_H[boundz];
    c4=c4z_H[boundz];

    for(i=0;i<xs-1;i++)
    {
        c5=cstx(c5x_H,i,xs-1);
        c6=cstx(c6x_H,i,xs-1);

        for(j=0;j<boundy;j++)
        {

			c1=c1y_H[j];
			c2=c2y_H[j];

        	for(k=boundz;k<zk;k++)
        	{
	            	r=k-boundz;

	                temp=bot_Bx[i][j][r];
	                bot_Bx[i][j][r]=c1*bot_Bx[i][j][r]+c2*(Ey[i+1][j][k+1]-Ey[i+1][j][k]-Ez[i+1][j+1][k]+Ez[i+1][j][k]);
	                Hx[i][j][k]=c3*Hx[i][j][k]+c4*(c5*bot_Bx[i][j][r]-c6*temp);

            }
        }
    }
}
void update_bot_Hy()
{
    typ_pola temp;
    ////////
    // Hy //
    ////////

    c1=c1z_H[boundz];
    c2=c2z_H[boundz];

    for(i=0;i<xs;i++)
    {
        c3=cstx(c3x_H,i,xs);
        c4=cstx(c4x_H,i,xs);

	 for(j=0;j<boundy;j++)
         {
         	c5=c5y_H[j];
         	c6=c6y_H[j];

        	for(k=boundz;k<zk;k++)
        	{
	            	r=k-boundz;

	                temp=bot_By[i][j][r];
	                bot_By[i][j][r]=c1*bot_By[i][j][r]+c2*(Ez[i+1][j+1][k]-Ez[i][j+1][k]-Ex[i][j+1][k+1]+Ex[i][j+1][k]);
                	Hy[i][j][k]=c3*Hy[i][j][k]+c4*(c5*bot_By[i][j][r]-c6*temp);

            }
        }
    }

}
void update_bot_Hz()
{
    typ_pola temp;

    ////////
    // Hz //
    ////////

	// w centralnej czesci bottom i top nie musimy zapamietywac Bz i Dz

    c5=c5z_H[boundz];
    c6=c6z_H[boundz];
    c1=c1z_H[boundz];
    c2=c2z_H[boundz];

    for(j=0;j<boundy;j++)
    {
    	c3=c3y_H[j];
     	c4=c4y_H[j];

		for(i=boundx;i<xk;i++)
  		{
    		for(k=boundz;k<zk-1;k++)
      		{
                Hz[i][j][k]=c3*Hz[i][j][k]+c4*c5*c2*(Ex[i][j+1][k+1]-Ex[i][j][k+1]-Ey[i+1][j][k+1]+Ey[i][j][k+1]);
            }
        }
    }

	// po bokach (left i right) musimy pamietac wszystko !!

	c5=c5x_H[boundx];
    c6=c6x_H[boundx];

    for(i=0;i<boundx;i++)
    {
        c1=c1x_H[i];
        c2=c2x_H[i];

		q=xs-i-1;

  		for(j=0;j<boundy;j++)
    	{

            c3=c3y_H[j];
            c4=c4y_H[j];

			for(k=boundz;k<zk-1;k++)
        	{
            	r=k-boundz;

				// update po lewej
                temp=bot_Bz_lef[i][j][r];
                bot_Bz_lef[i][j][r]=c1*bot_Bz_lef[i][j][r]+c2*(Ex[i][j+1][k+1]-Ex[i][j][k+1]-Ey[i+1][j][k+1]+Ey[i][j][k+1]);
                Hz[i][j][k]=c3*Hz[i][j][k]+c4*(c5*bot_Bz_lef[i][j][r]-c6*temp);

				// update po prawej
                temp=bot_Bz_rig[i][j][r];
                bot_Bz_rig[i][j][r]=c1*bot_Bz_rig[i][j][r]+c2*(Ex[q][j+1][k+1]-Ex[q][j][k+1]-Ey[q+1][j][k+1]+Ey[q][j][k+1]);
                Hz[q][j][k]=c3*Hz[q][j][k]+c4*(c5*bot_Bz_rig[i][j][r]-c6*temp);
            }
        }
    }
}

void update_top_Hx()
{

    typ_pola temp;

    ////////
    // Hx //
    ////////

    c3=c3z_H[boundz];
    c4=c4z_H[boundz];

    for(i=0;i<xs-1;i++)
    {
        c5=cstx(c5x_H,i,xs-1);
        c6=cstx(c6x_H,i,xs-1);

        for(j=0;j<boundy;j++)
        {

			c1=c1y_H[j];
			c2=c2y_H[j];

			p=ys-j-1;

        	for(k=boundz;k<zk;k++)
        	{
            	r=k-boundz;

                temp=top_Bx[i][j][r];
                top_Bx[i][j][r]=c1*top_Bx[i][j][r]+c2*(Ey[i+1][p][k+1]-Ey[i+1][p][k]-Ez[i+1][p+1][k]+Ez[i+1][p][k]);
                Hx[i][p][k]=c3*Hx[i][p][k]+c4*(c5*top_Bx[i][j][r]-c6*temp);
            }
        }
    }
}

void update_top_Hy()
{
    typ_pola temp;
    ////////
    // Hy //
    ////////

    c1=c1z_H[boundz];
    c2=c2z_H[boundz];

    for(i=0;i<xs;i++)
    {
        c3=cstx(c3x_H,i,xs);
        c4=cstx(c4x_H,i,xs);

		 for(j=0;j<boundy;j++)
         {
         	c5=c5y_H[j];
         	c6=c6y_H[j];

			p=ys-j-2;

        	for(k=boundz;k<zk;k++)
        	{
            	r=k-boundz;

                temp=top_By[i][j][r];
                top_By[i][j][r]=c1*top_By[i][j][r]+c2*(Ez[i+1][p+1][k]-Ez[i][p+1][k]-Ex[i][p+1][k+1]+Ex[i][p+1][k]);
                Hy[i][p][k]=c3*Hy[i][p][k]+c4*(c5*top_By[i][j][r]-c6*temp);
            }
        }
    }
}

void update_top_Hz()
{
    ////////
    // Hz //
    ////////
    typ_pola temp;
	// w centralnej czesci bottom i top nie musimy zapamietywac Bz i Dz

    c5=c5z_H[boundz];
    c6=c6z_H[boundz];
    c1=c1z_H[boundz];
    c2=c2z_H[boundz];

    for(j=0;j<boundy;j++)
    {
    	c3=c3y_H[j];
     	c4=c4y_H[j];

		p=ys-j-1;

		for(i=boundx;i<xk;i++)
  		{
    		for(k=boundz;k<zk-1;k++)
      		{
                Hz[i][p][k]=c3*Hz[i][p][k]+c4*c5*c2*(Ex[i][p+1][k+1]-Ex[i][p][k+1]-Ey[i+1][p][k+1]+Ey[i][p][k+1]);
            }
        }
    }

	// po bokach (left i right) musimy pamietac wszystko !!

	c5=c5x_H[boundx];
    c6=c6x_H[boundx];

    for(i=0;i<boundx;i++)
    {
        c1=c1x_H[i];
        c2=c2x_H[i];

		q=xs-i-1;

  		for(j=0;j<boundy;j++)
    	{

            c3=c3y_H[j];
            c4=c4y_H[j];

			p=ys-j-1;

			for(k=boundz;k<zk-1;k++)
        	{
            	r=k-boundz;

				// update po lewej
                temp=top_Bz_lef[i][j][r];
                top_Bz_lef[i][j][r]=c1*top_Bz_lef[i][j][r]+c2*(Ex[i][p+1][k+1]-Ex[i][p][k+1]-Ey[i+1][p][k+1]+Ey[i][p][k+1]);
                Hz[i][p][k]=c3*Hz[i][p][k]+c4*(c5*top_Bz_lef[i][j][r]-c6*temp);

				// update po prawej
                temp=top_Bz_rig[i][j][r];
                top_Bz_rig[i][j][r]=c1*top_Bz_rig[i][j][r]+c2*(Ex[q][p+1][k+1]-Ex[q][p][k+1]-Ey[q+1][p][k+1]+Ey[q][p][k+1]);
                Hz[q][p][k]=c3*Hz[q][p][k]+c4*(c5*top_Bz_rig[i][j][r]-c6*temp);
            }
        }
    }

}


void update_bot_Ex()
{

    typ_pola temp;

    ////////
    // Ex //
    ////////

    c3=c3z_E[boundz];
    c4=c4z_E[boundz];

    for(i=0;i<xs;i++)
    {
        c5=cstx(c5x_E,i,xs);
        c6=cstx(c6x_E,i,xs);

        for(j=1;j<boundy;j++)
        {
        	c1=c1y_E[j];
        	c2=c2y_E[j];

        	for(k=boundz;k<zk+1;k++)
        	{
            	r=k-boundz;

                temp=bot_Dx[i][j][r];
                bot_Dx[i][j][r]=c1*bot_Dx[i][j][r]+c2*(Hz[i][j][k-1]-Hz[i][j-1][k-1]-Hy[i][j-1][k]+Hy[i][j-1][k-1]);
                Ex[i][j][k]=c3*Ex[i][j][k]+c4*(c5*bot_Dx[i][j][r]-c6*temp);
            }
        }
    }
}

void update_bot_Ey()
{

    typ_pola temp;

    ////////
    // Ey //
    ////////

    c1=c1z_E[boundz];
    c2=c2z_E[boundz];

    for(i=1;i<xs;i++)
    {
        c3=cstx(c3x_E,i,xs+1);
        c4=cstx(c4x_E,i,xs+1);

        for(j=0;j<boundy;j++)
        {
        	c5=c5y_E[j];
        	c6=c6y_E[j];

        	for(k=boundz;k<zk+1;k++)
        	{
            	r=k-boundz;

                temp=bot_Dy[i][j][r];
                bot_Dy[i][j][r]=c1*bot_Dy[i][j][r]+c2*(Hx[i-1][j][k]-Hx[i-1][j][k-1]-Hz[i][j][k-1]+Hz[i-1][j][k-1]);
                Ey[i][j][k]=c3*Ey[i][j][k]+c4*(c5*bot_Dy[i][j][r]-c6*temp);

            }
        }
    }
}

void update_bot_Ez()
{

    typ_pola temp;
    ////////
    // Ez //
    ////////
	// w centralnej czesci bottom i top nie musimy zapamietywac Bz i Dz

    c5=c5z_E[boundz];
    c6=c6z_E[boundz];

	c1=c1z_E[boundz];
 	c2=c2z_E[boundz];

	for(j=1;j<boundy;j++)
	{
		c3=c3y_E[j];
		c4=c4y_E[j];

		for(i=boundx;i<=xk;i++)
    	{
        	for(k=boundz;k<zk;k++)
        	{
                Ez[i][j][k]=c3*Ez[i][j][k]+c4*c5*c2*(Hy[i][j-1][k]-Hy[i-1][j-1][k]-Hx[i-1][j][k]+Hx[i-1][j-1][k]);
            }
        }
    }

	// na brzegach ( left i right) musimy pamietac wszystko

	//c5=c5_E[bound];
    //c6=c6_E[bound];

	for(i=1;i<boundx;i++)
    {
        c1=c1x_E[i];
        c2=c2x_E[i];

		q=xs-i;

  		for(j=1;j<boundy;j++)
    	{
			c3=c3y_E[j];
			c4=c4y_E[j];

        	for(k=boundz;k<zk;k++)
        	{
            	r=k-boundz;

				// update po lewej
                temp=bot_Dz_lef[i][j][r];
                bot_Dz_lef[i][j][r]=c1*bot_Dz_lef[i][j][r]+c2*(Hy[i][j-1][k]-Hy[i-1][j-1][k]-Hx[i-1][j][k]+Hx[i-1][j-1][k]);
                Ez[i][j][k]=c3*Ez[i][j][k]+c4*(c5*bot_Dz_lef[i][j][r]-c6*temp);

				// update po prawej
                temp=bot_Dz_rig[i][j][r];
                bot_Dz_rig[i][j][r]=c1*bot_Dz_rig[i][j][r]+c2*(Hy[q][j-1][k]-Hy[q-1][j-1][k]-Hx[q-1][j][k]+Hx[q-1][j-1][k]);
                Ez[q][j][k]=c3*Ez[q][j][k]+c4*(c5*bot_Dz_rig[i][j][r]-c6*temp);

			}
        }
    }

}

void update_top_Ex()
{

    typ_pola temp;
 
    ////////
    // Ex //
    ////////

    c3=c3z_E[boundz];
    c4=c4z_E[boundz];

    for(i=0;i<xs;i++)
    {
        c5=cstx(c5x_E,i,xs);
        c6=cstx(c6x_E,i,xs);

        for(j=1;j<boundy;j++)
        {
        	c1=c1y_E[j];
        	c2=c2y_E[j];

			p=ys-j;

        	for(k=boundz;k<zk+1;k++)
        	{
            	r=k-boundz;

                temp=top_Dx[i][j][r];
                top_Dx[i][j][r]=c1*top_Dx[i][j][r]+c2*(Hz[i][p][k-1]-Hz[i][p-1][k-1]-Hy[i][p-1][k]+Hy[i][p-1][k-1]);
                Ex[i][p][k]=c3*Ex[i][p][k]+c4*(c5*top_Dx[i][j][r]-c6*temp);
            }
        }
    }
}


void update_top_Ey()
{

    typ_pola temp;

    ////////
    // Ey //
    ////////

    c1=c1z_E[boundz];
    c2=c2z_E[boundz];

    for(i=1;i<xs;i++)
    {
        c3=cstx(c3x_E,i,xs+1);
        c4=cstx(c4x_E,i,xs+1);

        for(j=0;j<boundy;j++)
        {
        	c5=c5y_E[j];
        	c6=c6y_E[j];

    		p=ys-j-1;

        	for(k=boundz;k<zk+1;k++)
        	{
            	r=k-boundz;

                temp=top_Dy[i][j][r];
                top_Dy[i][j][r]=c1*top_Dy[i][j][r]+c2*(Hx[i-1][p][k]-Hx[i-1][p][k-1]-Hz[i][p][k-1]+Hz[i-1][p][k-1]);
                Ey[i][p][k]=c3*Ey[i][p][k]+c4*(c5*top_Dy[i][j][r]-c6*temp);
            }
        }
    }
}



void update_top_Ez()
{

    typ_pola temp;
    ////////
    // Ez //
    ////////
	// w centralnej czesci bottom i top nie musimy zapamietywac Bz i Dz

    c5=c5z_E[boundz];
    c6=c6z_E[boundz];

	c1=c1z_E[boundz];
 	c2=c2z_E[boundz];

	for(j=1;j<boundy;j++)
	{
		c3=c3y_E[j];
		c4=c4y_E[j];
		p=ys-j;

		for(i=boundx;i<=xk;i++)
    	{
        	for(k=boundz;k<zk;k++)
        	{
                Ez[i][p][k]=c3*Ez[i][p][k]+c4*c5*c2*(Hy[i][p-1][k]-Hy[i-1][p-1][k]-Hx[i-1][p][k]+Hx[i-1][p-1][k]);
            }
        }
    }

	// na brzegach ( left i right) musimy pamietac wszystko

	//c5=c5_E[bound];
    //c6=c6_E[bound];

	for(i=1;i<boundx;i++)
    {
        c1=c1x_E[i];
        c2=c2x_E[i];

		q=xs-i;

  		for(j=1;j<boundy;j++)
    	{
			c3=c3y_E[j];
			c4=c4y_E[j];

			p=ys-j;

        	for(k=boundz;k<zk;k++)
        	{
            	r=k-boundz;

				// update po lewej
                temp=top_Dz_lef[i][j][r];
                top_Dz_lef[i][j][r]=c1*top_Dz_lef[i][j][r]+c2*(Hy[i][p-1][k]-Hy[i-1][p-1][k]-Hx[i-1][p][k]+Hx[i-1][p-1][k]);
                Ez[i][p][k]=c3*Ez[i][p][k]+c4*(c5*top_Dz_lef[i][j][r]-c6*temp);

				// update po prawej
                temp=top_Dz_rig[i][j][r];
                top_Dz_rig[i][j][r]=c1*top_Dz_rig[i][j][r]+c2*(Hy[q][p-1][k]-Hy[q-1][p-1][k]-Hx[q-1][p][k]+Hx[q-1][p-1][k]);
                Ez[q][p][k]=c3*Ez[q][p][k]+c4*(c5*top_Dz_rig[i][j][r]-c6*temp);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void update_lef_Hx()
{
    typ_pola temp;

    ////////
    // Hx //
    ////////

    c1=c1x_H[boundx];
    c2=c2x_H[boundx];

    c3=c3x_H[boundx];
    c4=c4x_H[boundx];

	for(i=0;i<boundx;i++)
	{
		c5=c5x_H[i];
		c6=c6x_H[i];

    	for(j=boundy;j<yk;j++)
    	{
        	q=j-boundy;

        	for(k=boundz;k<zk;k++)
        	{
            	r=k-boundz;

                temp=lef_Bx[i][q][r];
                lef_Bx[i][q][r]=c1*lef_Bx[i][q][r]+c2*(Ey[i+1][j][k+1]-Ey[i+1][j][k]-Ez[i+1][j+1][k]+Ez[i+1][j][k]);
                Hx[i][j][k]=c3*Hx[i][j][k]+c4*(c5*lef_Bx[i][q][r]-c6*temp);

            }
        }
    }
}
void update_lef_Hy()
{
    ////////
    // Hy //
    ////////

	// zauwazmy ,ze nie trzeba przechowywac By,Dy bo c1=1, a c5=c6
	// co oszczedza nam pamiec i upraszcza wzory

    // ustawione przy Hx
    c1=c1x_H[boundx];
    c2=c2x_H[boundx];

	c5=c5x_H[boundx];
    c6=c6x_H[boundx];

	for(i=0;i<boundx;i++)
	{
		c3=c3x_H[i];
		c4=c4x_H[i];

    	for(k=boundz;k<zk;k++)
    	{
        	for(j=boundy;j<yk-1;j++)
        	{
                Hy[i][j][k]=c3*Hy[i][j][k]+c4*c5*c2*(Ez[i+1][j+1][k]-Ez[i][j+1][k]-Ex[i][j+1][k+1]+Ex[i][j+1][k]);
            }
        }
    }
}
void update_lef_Hz()
{
    typ_pola temp;
    ////////
    // Hz //
    ////////

    c3=c3x_H[boundx];
    c4=c4x_H[boundx];

	// ustawione przy Hy
 	c5=c5x_H[boundx];
    c6=c6x_H[boundx];

	for(i=0;i<boundx;i++)
	{
		c1=c1x_H[i];
		c2=c2x_H[i];

		for(j=boundy;j<yk;j++)
    	{
        	q=j-boundy;

        	for(k=boundz;k<zk-1;k++)
        	{
            	r=k-boundz;

                temp=lef_Bz[i][q][r];
                lef_Bz[i][q][r]=c1*lef_Bz[i][q][r]+c2*(Ex[i][j+1][k+1]-Ex[i][j][k+1]-Ey[i+1][j][k+1]+Ey[i][j][k+1]);
                Hz[i][j][k]=c3*Hz[i][j][k]+c4*(c5*lef_Bz[i][q][r]-c6*temp);

            }
        }
    }

}

void update_rig_Hx()
{

    typ_pola temp;

    ////////
    // Hx //
    ////////

    c1=c1x_H[boundx];
    c2=c2x_H[boundx];

    c3=c3x_H[boundx];
    c4=c4x_H[boundx];

	for(i=0;i<boundx;i++)
	{
		c5=c5x_H[i];
		c6=c6x_H[i];

		p=xs-i-2;

    	for(j=boundy;j<yk;j++)
    	{
        	q=j-boundy;

        	for(k=boundz;k<zk;k++)
        	{
            	r=k-boundz;

                temp=rig_Bx[i][q][r];
                rig_Bx[i][q][r]=c1*rig_Bx[i][q][r]+c2*(Ey[p+1][j][k+1]-Ey[p+1][j][k]-Ez[p+1][j+1][k]+Ez[p+1][j][k]);
                Hx[p][j][k]=c3*Hx[p][j][k]+c4*(c5*rig_Bx[i][q][r]-c6*temp);
            }
        }
    }
}

void update_rig_Hy()
{

    
    ////////
    // Hy //
    ////////

	// zauwazmy ,ze nie trzeba przechowywac By,Dy bo c1=1, a c5=c6
	// co oszczedza nam pamiec i upraszcza wzory

    // ustawione przy Hx
    c1=c1x_H[boundx];
    c2=c2x_H[boundx];

	c5=c5x_H[boundx];
    c6=c6x_H[boundx];

	for(i=0;i<boundx;i++)
	{
		c3=c3x_H[i];
		c4=c4x_H[i];

		p=xs-i-1;

    	for(k=boundz;k<zk;k++)
    	{
        	for(j=boundy;j<yk-1;j++)
        	{
                Hy[p][j][k]=c3*Hy[p][j][k]+c4*c5*c2*(Ez[p+1][j+1][k]-Ez[p][j+1][k]-Ex[p][j+1][k+1]+Ex[p][j+1][k]);
            }
        }
    }

}

void update_rig_Hz()
{

    typ_pola temp;
    ////////
    // Hz //
    ////////

    c3=c3x_H[boundx];
    c4=c4x_H[boundx];

	// ustawione przy Hy
 	c5=c5x_H[boundx];
    c6=c6x_H[boundx];

	for(i=0;i<boundx;i++)
	{
		c1=c1x_H[i];
		c2=c2x_H[i];

		p=xs-i-1;

		for(j=boundy;j<yk;j++)
    	{
        	q=j-boundy;

        	for(k=boundz;k<zk-1;k++)
        	{
            	r=k-boundz;

                temp=rig_Bz[i][q][r];
                rig_Bz[i][q][r]=c1*rig_Bz[i][q][r]+c2*(Ex[p][j+1][k+1]-Ex[p][j][k+1]-Ey[p+1][j][k+1]+Ey[p][j][k+1]);
                Hz[p][j][k]=c3*Hz[p][j][k]+c4*(c5*rig_Bz[i][q][r]-c6*temp);
            }
        }
    }

}


void update_lef_Ex()
{

    typ_pola temp;

    ////////
    // Ex //
    ////////

    c1=c1x_E[boundx];
    c2=c2x_E[boundx];

    c3=c3x_E[boundx];
    c4=c4x_E[boundx];

	for(i=0;i<boundx;i++)
	{
		c5=c5x_E[i];
    	c6=c6x_E[i];

    	for(j=boundy;j<yk+1;j++)
    	{
        	q=j-boundy;

        	for(k=boundz;k<zk+1;k++)
        	{
            	r=k-boundz;

                temp=lef_Dx[i][q][r];
                lef_Dx[i][q][r]=c1*lef_Dx[i][q][r]+c2*(Hz[i][j][k-1]-Hz[i][j-1][k-1]-Hy[i][j-1][k]+Hy[i][j-1][k-1]);
                Ex[i][j][k]=c3*Ex[i][j][k]+c4*(c5*lef_Dx[i][q][r]-c6*temp);

            }
        }
    }
}


void update_lef_Ey()
{

    ////////
    // Ey //
    ////////

	// zauwazmy ,ze nie trzeba przechowywac By,Dy bo c1=1, a c5=c6
	// co oszczedza nam pamiec i upraszcza wzory

    //ustawione przy Ex
    c1=c1x_E[boundx];
    c2=c2x_E[boundx];

    c5=c5x_E[boundx];
    c6=c6x_E[boundx];

	for(i=1;i<boundx;i++)
	{
		c3=c3x_E[i];
		c4=c4x_E[i];

    	for(k=boundz;k<zk+1;k++)
    	{
        	for(j=boundy;j<yk;j++)
        	{
                Ey[i][j][k]=c3*Ey[i][j][k]+c4*c5*c2*(Hx[i-1][j][k]-Hx[i-1][j][k-1]-Hz[i][j][k-1]+Hz[i-1][j][k-1]);
            }
        }
    }
}


void update_lef_Ez()
{

    typ_pola temp;
    ////////
    // Ez //
    ////////

    // ustawione przy Ey
    c5=c5x_E[boundx];
    c6=c6x_E[boundx];

    c3=c3x_E[boundx];
    c4=c4x_E[boundx];

	for(i=1;i<boundx;i++)
	{
		c1=c1x_E[i];
		c2=c2x_E[i];

    	for(j=boundy;j<yk+1;j++)
    	{
        	q=j-boundy;

        	for(k=boundz;k<zk;k++)
        	{
            	r=k-boundz;

                temp=lef_Dz[i][q][r];
                lef_Dz[i][q][r]=c1*lef_Dz[i][q][r]+c2*(Hy[i][j-1][k]-Hy[i-1][j-1][k]-Hx[i-1][j][k]+Hx[i-1][j-1][k]);
                Ez[i][j][k]=c3*Ez[i][j][k]+c4*(c5*lef_Dz[i][q][r]-c6*temp);
            }
        }
    }
}

void update_rig_Ex()
{

    typ_pola temp;

    ////////
    // Ex //
    ////////

    c1=c1x_E[boundx];
    c2=c2x_E[boundx];

    c3=c3x_E[boundx];
    c4=c4x_E[boundx];

	for(i=0;i<boundx;i++)
	{
		c5=c5x_E[i];
    	c6=c6x_E[i];

		p=xs-i-1;

    	for(j=boundy;j<yk+1;j++)
    	{
        	q=j-boundy;

        	for(k=boundz;k<zk+1;k++)
        	{
            	r=k-boundz;

                temp=rig_Dx[i][q][r];
                rig_Dx[i][q][r]=c1*rig_Dx[i][q][r]+c2*(Hz[p][j][k-1]-Hz[p][j-1][k-1]-Hy[p][j-1][k]+Hy[p][j-1][k-1]);
                Ex[p][j][k]=c3*Ex[p][j][k]+c4*(c5*rig_Dx[i][q][r]-c6*temp);
            }
        }
    }
}

void update_rig_Ey()
{

    ////////
    // Ey //
    ////////

	// zauwazmy ,ze nie trzeba przechowywac By,Dy bo c1=1, a c5=c6
	// co oszczedza nam pamiec i upraszcza wzory

    //ustawione przy Ex
    c1=c1x_E[boundx];
    c2=c2x_E[boundx];

    c5=c5x_E[boundx];
    c6=c6x_E[boundx];

	for(i=1;i<boundx;i++)
	{
		c3=c3x_E[i];
		c4=c4x_E[i];
		p=xs-i;

    	for(k=boundz;k<zk+1;k++)
    	{
        	for(j=boundy;j<yk;j++)
        	{
                Ey[p][j][k]=c3*Ey[p][j][k]+c4*c5*c2*(Hx[p-1][j][k]-Hx[p-1][j][k-1]-Hz[p][j][k-1]+Hz[p-1][j][k-1]);
            }
        }
    }
}

void update_rig_Ez()
{

    typ_pola temp;
    ////////
    // Ez //
    ////////

    // ustawione przy Ey
    c5=c5x_E[boundx];
    c6=c6x_E[boundx];

    c3=c3x_E[boundx];
    c4=c4x_E[boundx];

	for(i=1;i<boundx;i++)
	{
		c1=c1x_E[i];
		c2=c2x_E[i];

  		p=xs-i;

    	for(j=boundy;j<yk+1;j++)
    	{
        	q=j-boundy;

        	for(k=boundz;k<zk;k++)
        	{
            	r=k-boundz;

                temp=rig_Dz[i][q][r];
                rig_Dz[i][q][r]=c1*rig_Dz[i][q][r]+c2*(Hy[p][j-1][k]-Hy[p-1][j-1][k]-Hx[p-1][j][k]+Hx[p-1][j-1][k]);
                Ez[p][j][k]=c3*Ez[p][j][k]+c4*(c5*rig_Dz[i][q][r]-c6*temp);
            }
        }
    }
}



void update_bac_Hx()
{

    typ_pola temp;

    ////////
    // Hx //
    ////////

	// update na gorze i dole
	for(i=0;i<xs-1;i++)
	{
		c5=cstx(c5x_H,i,xs-1);
		c6=cstx(c6x_H,i,xs-1);

    	for(j=0;j<boundy;j++)
    	{
        	c1=c1y_H[j];
        	c2=c2y_H[j];

			q=ys-j-1;

            for(k=0;k<boundz;k++)
            {
                c3=c3z_H[k];
                c4=c4z_H[k];

				// update na dole
                temp=bac_Bx_bot[i][j][k];
                bac_Bx_bot[i][j][k]=c1*bac_Bx_bot[i][j][k]+c2*(Ey[i+1][j][k+1]-Ey[i+1][j][k]-Ez[i+1][j+1][k]+Ez[i+1][j][k]);
                Hx[i][j][k]=c3*Hx[i][j][k]+c4*(c5*bac_Bx_bot[i][j][k]-c6*temp);

				// update na gorze
                temp=bac_Bx_top[i][j][k];
                bac_Bx_top[i][j][k]=c1*bac_Bx_top[i][j][k]+c2*(Ey[i+1][q][k+1]-Ey[i+1][q][k]-Ez[i+1][q+1][k]+Ez[i+1][q][k]);
                Hx[i][q][k]=c3*Hx[i][q][k]+c4*(c5*bac_Bx_top[i][j][k]-c6*temp);

            }
        }
    }

	// update pomiedzy niebem i pieklem ( gora i dolem)

	c1=c1z_H[boundz];
 	c2=c2z_H[boundz];

	c5=c5z_H[boundz];
	c6=c6z_H[boundz];

	// update centrum
	for(k=0;k<boundz;k++)
	{
		c3=c3z_H[k];
		c4=c4z_H[k];

		for(j=boundy;j<yk;j++)
    	{
        	for(i=boundx;i<xk-1;i++)
        	{
                Hx[i][j][k]=c3*Hx[i][j][k]+c4*c5*c2*(Ey[i+1][j][k+1]-Ey[i+1][j][k]-Ez[i+1][j+1][k]+Ez[i+1][j][k]);
            }
        }
    }


	// update lewicy i prawicy
	//c1=c1_H[bound];
 	//c2=c2_H[bound];

	for(i=0;i<boundx;i++)
	{
		c5=c5x_H[i];
		c6=c6x_H[i];

		q=xs-i-2;

		for(k=0;k<boundz;k++)
		{
			c3=c3z_H[k];
			c4=c4z_H[k];


    		for(j=boundy;j<yk;j++)
    		{
     			r=j-boundy;

				// update left
                temp=bac_Bx_lef[i][r][k];
                bac_Bx_lef[i][r][k]=c1*bac_Bx_lef[i][r][k]+c2*(Ey[i+1][j][k+1]-Ey[i+1][j][k]-Ez[i+1][j+1][k]+Ez[i+1][j][k]);
                Hx[i][j][k]=c3*Hx[i][j][k]+c4*(c5*bac_Bx_lef[i][r][k]-c6*temp);

				// update rig
				temp=bac_Bx_rig[i][r][k];
                bac_Bx_rig[i][r][k]=c1*bac_Bx_rig[i][r][k]+c2*(Ey[q+1][j][k+1]-Ey[q+1][j][k]-Ez[q+1][j+1][k]+Ez[q+1][j][k]);
                Hx[q][j][k]=c3*Hx[q][j][k]+c4*(c5*bac_Bx_rig[i][r][k]-c6*temp);

            }
        }
    }
}

void update_bac_Hy()
{

    typ_pola temp;
    
    ////////
    // Hy //
    ////////

    for(i=0;i<xs;i++)
    {
        c3=cstx(c3x_H,i,xs);
        c4=cstx(c4x_H,i,xs);

        for(j=0;j<ys-1;j++)
        {
            c5=csty(c5y_H,j,ys-1);
            c6=csty(c6y_H,j,ys-1);

            for(k=0;k<boundz;k++)
            {
                c1=c1z_H[k];
                c2=c2z_H[k];

                temp=bac_By[i][j][k];
                bac_By[i][j][k]=c1*bac_By[i][j][k]+c2*(Ez[i+1][j+1][k]-Ez[i][j+1][k]-Ex[i][j+1][k+1]+Ex[i][j+1][k]);
                Hy[i][j][k]=c3*Hy[i][j][k]+c4*(c5*bac_By[i][j][k]-c6*temp);

            }

        }
    }
}


void update_bac_Hz()
{

    typ_pola temp;

    ////////
    // Hz //
    ////////

    for(i=0;i<xs;i++)
    {
        c1=cstx(c1x_H,i,xs);
        c2=cstx(c2x_H,i,xs);

        for(j=0;j<ys;j++)
        {
            c3=csty(c3y_H,j,ys);
            c4=csty(c4y_H,j,ys);

            for(k=0;k<boundz;k++)
            {
                c5=c5z_H[k];
                c6=c6z_H[k];

                temp=bac_Bz[i][j][k];
                bac_Bz[i][j][k]=c1*bac_Bz[i][j][k]+c2*(Ex[i][j+1][k+1]-Ex[i][j][k+1]-Ey[i+1][j][k+1]+Ey[i][j][k+1]);
                Hz[i][j][k]=c3*Hz[i][j][k]+c4*(c5*bac_Bz[i][j][k]-c6*temp);

            }
        }
    }

}


void update_fro_Hx()
{

    typ_pola temp;

    ////////
    // Hx //
    ////////

	// update na gorze i dole
	for(i=0;i<xs-1;i++)
	{
		c5=cstx(c5x_H,i,xs-1);
		c6=cstx(c6x_H,i,xs-1);

    	for(j=0;j<boundy;j++)
    	{
        	c1=c1y_H[j];
        	c2=c2y_H[j];

			q=ys-j-1;

            for(k=0;k<boundz;k++)
            {
                c3=c3z_H[k];
                c4=c4z_H[k];

                p=zs-k-1;

				// update na dole
				temp=fro_Bx_bot[i][j][k];
                fro_Bx_bot[i][j][k]=c1*fro_Bx_bot[i][j][k]+c2*(Ey[i+1][j][p+1]-Ey[i+1][j][p]-Ez[i+1][j+1][p]+Ez[i+1][j][p]);
                Hx[i][j][p]=c3*Hx[i][j][p]+c4*(c5*fro_Bx_bot[i][j][k]-c6*temp);

				// update na gorze
                temp=fro_Bx_top[i][j][k];
                fro_Bx_top[i][j][k]=c1*fro_Bx_top[i][j][k]+c2*(Ey[i+1][q][p+1]-Ey[i+1][q][p]-Ez[i+1][q+1][p]+Ez[i+1][q][p]);
                Hx[i][q][p]=c3*Hx[i][q][p]+c4*(c5*fro_Bx_top[i][j][k]-c6*temp);
            }
        }
    }

	// update pomiedzy niebem i pieklem ( gora i dolem)

	c1=c1z_H[boundz];
 	c2=c2z_H[boundz];

	c5=c5z_H[boundz];
	c6=c6z_H[boundz];

	// update centrum
	for(k=0;k<boundz;k++)
	{
		c3=c3z_H[k];
		c4=c4z_H[k];

		p=zs-k-1;

		for(j=boundy;j<yk;j++)
    	{
        	for(i=boundx;i<xk-1;i++)
        	{
                Hx[i][j][p]=c3*Hx[i][j][p]+c4*c5*c2*(Ey[i+1][j][p+1]-Ey[i+1][j][p]-Ez[i+1][j+1][p]+Ez[i+1][j][p]);
            }
        }
    }


	// update lewicy i prawicy
	//c1=c1_H[bound];
 	//c2=c2_H[bound];

	for(i=0;i<boundx;i++)
	{
		c5=c5x_H[i];
		c6=c6x_H[i];

		q=xs-i-2;

		for(k=0;k<boundz;k++)
		{
			c3=c3z_H[k];
			c4=c4z_H[k];

			p=zs-k-1;

    		for(j=boundy;j<yk;j++)
    		{
     			r=j-boundy;

				// update left
                temp=fro_Bx_lef[i][r][k];
                fro_Bx_lef[i][r][k]=c1*fro_Bx_lef[i][r][k]+c2*(Ey[i+1][j][p+1]-Ey[i+1][j][p]-Ez[i+1][j+1][p]+Ez[i+1][j][p]);
                Hx[i][j][p]=c3*Hx[i][j][p]+c4*(c5*fro_Bx_lef[i][r][k]-c6*temp);

				// update rig
                temp=fro_Bx_rig[i][r][k];
                fro_Bx_rig[i][r][k]=c1*fro_Bx_rig[i][r][k]+c2*(Ey[q+1][j][p+1]-Ey[q+1][j][p]-Ez[q+1][j+1][p]+Ez[q+1][j][p]);
                Hx[q][j][p]=c3*Hx[q][j][p]+c4*(c5*fro_Bx_rig[i][r][k]-c6*temp);

            }
        }
    }
}


void update_fro_Hy()
{

    typ_pola temp;
    ////////
    // Hy //
    ////////

    for(i=0;i<xs;i++)
    {
        c3=cstx(c3x_H,i,xs);
        c4=cstx(c4x_H,i,xs);

        for(j=0;j<ys-1;j++)
        {
            c5=csty(c5y_H,j,ys-1);
            c6=csty(c6y_H,j,ys-1);

            for(k=0;k<boundz;k++)
            {
                c1=c1z_H[k];
                c2=c2z_H[k];

				p=zs-k-1;

                temp=fro_By[i][j][k];
                fro_By[i][j][k]=c1*fro_By[i][j][k]+c2*(Ez[i+1][j+1][p]-Ez[i][j+1][p]-Ex[i][j+1][p+1]+Ex[i][j+1][p]);
                Hy[i][j][p]=c3*Hy[i][j][p]+c4*(c5*fro_By[i][j][k]-c6*temp);
            }

        }
    }
}


void update_fro_Hz()
{

    typ_pola temp;

    ////////
    // Hz //
    ////////

    for(i=0;i<xs;i++)
    {
        c1=cstx(c1x_H,i,xs);
        c2=cstx(c2x_H,i,xs);

        for(j=0;j<ys;j++)
        {
            c3=csty(c3y_H,j,ys);
            c4=csty(c4y_H,j,ys);

            for(k=0;k<boundz;k++)
            {
                c5=c5z_H[k];
                c6=c6z_H[k];

				p=zs-k-2;

                temp=fro_Bz[i][j][k];
                fro_Bz[i][j][k]=c1*fro_Bz[i][j][k]+c2*(Ex[i][j+1][p+1]-Ex[i][j][p+1]-Ey[i+1][j][p+1]+Ey[i][j][p+1]);
                Hz[i][j][p]=c3*Hz[i][j][p]+c4*(c5*fro_Bz[i][j][k]-c6*temp);
            }
        }
    }

}

void update_bac_Ex()
{

    typ_pola temp;

    ////////
    // Ex //
    ////////

	for(i=0;i<xs;i++)
	{
		c5=cstx(c5x_E,i,xs);
		c6=cstx(c6x_E,i,xs);

    	for(j=1;j<boundy;j++)
    	{
        	c1=c1y_E[j];
        	c2=c2y_E[j];

			q=ys-j;

            for(k=1;k<boundz;k++)
            {
                c3=c3z_E[k];
                c4=c4z_E[k];

				p=zs-k;

				temp=bac_Dx_bot[i][j][k];
                bac_Dx_bot[i][j][k]=c1*bac_Dx_bot[i][j][k]+c2*(Hz[i][j][k-1]-Hz[i][j-1][k-1]-Hy[i][j-1][k]+Hy[i][j-1][k-1]);
                Ex[i][j][k]=c3*Ex[i][j][k]+c4*(c5*bac_Dx_bot[i][j][k]-c6*temp);

                temp=bac_Dx_top[i][j][k];
                bac_Dx_top[i][j][k]=c1*bac_Dx_top[i][j][k]+c2*(Hz[i][q][k-1]-Hz[i][q-1][k-1]-Hy[i][q-1][k]+Hy[i][q-1][k-1]);
                Ex[i][q][k]=c3*Ex[i][q][k]+c4*(c5*bac_Dx_top[i][j][k]-c6*temp);

			}
        }
    }

	c1=c1z_E[boundz];
 	c2=c2z_E[boundz];

	c5=c5z_E[boundz];
	c6=c6z_E[boundz];

	for(j=boundy;j<=yk;j++)
	{
        for(i=boundx;i<xk;i++)
        {
            for(k=1;k<boundz;k++)
            {
                c3=c3z_E[k];
                c4=c4z_E[k];

                Ex[i][j][k]=c3*Ex[i][j][k]+c4*c5*c2*(Hz[i][j][k-1]-Hz[i][j-1][k-1]-Hy[i][j-1][k]+Hy[i][j-1][k-1]);
            }
        }
    }

	//c1=c1_E[bound];
	//c2=c2_E[bound];

	for(i=0;i<boundx;i++)
	{
		c5=c5x_E[i];
		c6=c6x_E[i];

		q=xs-i-1;

  		for(k=1;k<boundz;k++)
    	{
			c3=c3z_E[k];
			c4=c4z_E[k];

			p=zs-k;

			for(j=boundy;j<=yk;j++)
    		{
    			r=j-boundy;

				temp=bac_Dx_lef[i][r][k];
                bac_Dx_lef[i][r][k]=c1*bac_Dx_lef[i][r][k]+c2*(Hz[i][j][k-1]-Hz[i][j-1][k-1]-Hy[i][j-1][k]+Hy[i][j-1][k-1]);
                Ex[i][j][k]=c3*Ex[i][j][k]+c4*(c5*bac_Dx_lef[i][r][k]-c6*temp);

				temp=bac_Dx_rig[i][r][k];
                bac_Dx_rig[i][r][k]=c1*bac_Dx_rig[i][r][k]+c2*(Hz[q][j][k-1]-Hz[q][j-1][k-1]-Hy[q][j-1][k]+Hy[q][j-1][k-1]);
                Ex[q][j][k]=c3*Ex[q][j][k]+c4*(c5*bac_Dx_rig[i][r][k]-c6*temp);

            }
        }
    }
}

void update_bac_Ey()
{

    typ_pola temp;

    ////////
    // Ey //
    ////////

    for(i=1;i<xs;i++)
    {
        c3=cstx(c3x_E,i,xs+1);
        c4=cstx(c4x_E,i,xs+1);

        for(j=0;j<ys;j++)
        {
            c5=csty(c5y_E,j,ys);
            c6=csty(c6y_E,j,ys);

            for(k=1;k<boundz;k++)
            {
                c1=c1z_E[k];
                c2=c2z_E[k];

		p=zs-k;

                temp=bac_Dy[i][j][k];
                bac_Dy[i][j][k]=c1*bac_Dy[i][j][k]+c2*(Hx[i-1][j][k]-Hx[i-1][j][k-1]-Hz[i][j][k-1]+Hz[i-1][j][k-1]);
                Ey[i][j][k]=c3*Ey[i][j][k]+c4*(c5*bac_Dy[i][j][k]-c6*temp);

            }
        }
    }
}

void update_bac_Ez()
{

    typ_pola temp;

    ////////
    // Ez //
    ////////

    for(i=1;i<xs;i++)
    {
        c1=cstx(c1x_E,i,xs+1);
        c2=cstx(c2x_E,i,xs+1);

        for(j=1;j<ys;j++)
        {
            c3=csty(c3y_E,j,ys+1);
            c4=csty(c4y_E,j,ys+1);

            for(k=0;k<boundz;k++)
            {
                c5=c5z_E[k];
                c6=c6z_E[k];

				p=zs-k-1;

                temp=bac_Dz[i][j][k];
                bac_Dz[i][j][k]=c1*bac_Dz[i][j][k]+c2*(Hy[i][j-1][k]-Hy[i-1][j-1][k]-Hx[i-1][j][k]+Hx[i-1][j-1][k]);
                Ez[i][j][k]=c3*Ez[i][j][k]+c4*(c5*bac_Dz[i][j][k]-c6*temp);

            }
        }
    }

}

void update_fro_Ex()
{

    typ_pola temp;

    ////////
    // Ex //
    ////////

	for(i=0;i<xs;i++)
	{
		c5=cstx(c5x_E,i,xs);
		c6=cstx(c6x_E,i,xs);

    	for(j=1;j<boundy;j++)
    	{
        	c1=c1y_E[j];
        	c2=c2y_E[j];

			q=ys-j;

            for(k=1;k<boundz;k++)
            {
                c3=c3z_E[k];
                c4=c4z_E[k];

				p=zs-k;

                temp=fro_Dx_bot[i][j][k];
                fro_Dx_bot[i][j][k]=c1*fro_Dx_bot[i][j][k]+c2*(Hz[i][j][p-1]-Hz[i][j-1][p-1]-Hy[i][j-1][p]+Hy[i][j-1][p-1]);
                Ex[i][j][p]=c3*Ex[i][j][p]+c4*(c5*fro_Dx_bot[i][j][k]-c6*temp);

                temp=fro_Dx_top[i][j][k];
                fro_Dx_top[i][j][k]=c1*fro_Dx_top[i][j][k]+c2*(Hz[i][q][p-1]-Hz[i][q-1][p-1]-Hy[i][q-1][p]+Hy[i][q-1][p-1]);
                Ex[i][q][p]=c3*Ex[i][q][p]+c4*(c5*fro_Dx_top[i][j][k]-c6*temp);

			}
        }
    }

	c1=c1z_E[boundz];
 	c2=c2z_E[boundz];

	c5=c5z_E[boundz];
	c6=c6z_E[boundz];

	for(j=boundy;j<=yk;j++)
	{
        for(i=boundx;i<xk;i++)
        {
            for(k=1;k<boundz;k++)
            {
                c3=c3z_E[k];
                c4=c4z_E[k];

                p=zs-k;

                Ex[i][j][p]=c3*Ex[i][j][p]+c4*c5*c2*(Hz[i][j][p-1]-Hz[i][j-1][p-1]-Hy[i][j-1][p]+Hy[i][j-1][p-1]);
            }
        }
    }

	//c1=c1_E[bound];
 	//c2=c2_E[bound];

	for(i=0;i<boundx;i++)
	{
		c5=c5x_E[i];
		c6=c6x_E[i];

		q=xs-i-1;

  		for(k=1;k<boundz;k++)
    	{
			c3=c3z_E[k];
			c4=c4z_E[k];

			p=zs-k;

			for(j=boundy;j<=yk;j++)
    		{
    			r=j-boundy;

                temp=fro_Dx_lef[i][r][k];
                fro_Dx_lef[i][r][k]=c1*fro_Dx_lef[i][r][k]+c2*(Hz[i][j][p-1]-Hz[i][j-1][p-1]-Hy[i][j-1][p]+Hy[i][j-1][p-1]);
                Ex[i][j][p]=c3*Ex[i][j][p]+c4*(c5*fro_Dx_lef[i][r][k]-c6*temp);

				temp=fro_Dx_rig[i][r][k];
                fro_Dx_rig[i][r][k]=c1*fro_Dx_rig[i][r][k]+c2*(Hz[q][j][p-1]-Hz[q][j-1][p-1]-Hy[q][j-1][p]+Hy[q][j-1][p-1]);
                Ex[q][j][p]=c3*Ex[q][j][p]+c4*(c5*fro_Dx_rig[i][r][k]-c6*temp);
            }
        }
    }
}

void update_fro_Ey()
{

    typ_pola temp;

    ////////
    // Ey //
    ////////

    for(i=1;i<xs;i++)
    {
        c3=cstx(c3x_E,i,xs+1);
        c4=cstx(c4x_E,i,xs+1);

        for(j=0;j<ys;j++)
        {
            c5=csty(c5y_E,j,ys);
            c6=csty(c6y_E,j,ys);

            for(k=1;k<boundz;k++)
            {
                c1=c1z_E[k];
                c2=c2z_E[k];

				p=zs-k;

                temp=fro_Dy[i][j][k];
                fro_Dy[i][j][k]=c1*fro_Dy[i][j][k]+c2*(Hx[i-1][j][p]-Hx[i-1][j][p-1]-Hz[i][j][p-1]+Hz[i-1][j][p-1]);
                Ey[i][j][p]=c3*Ey[i][j][p]+c4*(c5*fro_Dy[i][j][k]-c6*temp);
            }
        }
    }
}

void update_fro_Ez()
{

    typ_pola temp;

    ////////
    // Ez //
    ////////

    for(i=1;i<xs;i++)
    {
        c1=cstx(c1x_E,i,xs+1);
        c2=cstx(c2x_E,i,xs+1);

        for(j=1;j<ys;j++)
        {
            c3=csty(c3y_E,j,ys+1);
            c4=csty(c4y_E,j,ys+1);

            for(k=0;k<boundz;k++)
            {
                c5=c5z_E[k];
                c6=c6z_E[k];

				p=zs-k-1;

                temp=fro_Dz[i][j][k];
                fro_Dz[i][j][k]=c1*fro_Dz[i][j][k]+c2*(Hy[i][j-1][p]-Hy[i-1][j-1][p]-Hx[i-1][j][p]+Hx[i-1][j-1][p]);
                Ez[i][j][p]=c3*Ez[i][j][p]+c4*(c5*fro_Dz[i][j][k]-c6*temp);
            }
        }
    }

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// funkcje zubo¿one
//
void update_bot_Hx_null_Ez()
{

    typ_pola temp;

    ////////
    // Hx //
    ////////

    c3=c3z_H[boundz];
    c4=c4z_H[boundz];

    for(i=0;i<xs-1;i++)
    {
        c5=cstx(c5x_H,i,xs-1);
        c6=cstx(c6x_H,i,xs-1);

        for(j=0;j<boundy;j++)
        {

			c1=c1y_H[j];
			c2=c2y_H[j];

        	for(k=boundz;k<zk;k++)
        	{
	            	r=k-boundz;

	                temp=bot_Bx[i][j][r];
	                bot_Bx[i][j][r]=c1*bot_Bx[i][j][r]+c2*(Ey[i+1][j][k+1]-Ey[i+1][j][k]);
	                Hx[i][j][k]=c3*Hx[i][j][k]+c4*(c5*bot_Bx[i][j][r]-c6*temp);

            }
        }
    }
}

void update_bot_Hx_null_Ey()
{

    typ_pola temp;

    ////////
    // Hx //
    ////////

    c3=c3z_H[boundz];
    c4=c4z_H[boundz];

    for(i=0;i<xs-1;i++)
    {
        c5=cstx(c5x_H,i,xs-1);
        c6=cstx(c6x_H,i,xs-1);

        for(j=0;j<boundy;j++)
        {

			c1=c1y_H[j];
			c2=c2y_H[j];

        	for(k=boundz;k<zk;k++)
        	{
	            	r=k-boundz;

	                temp=bot_Bx[i][j][r];
	                bot_Bx[i][j][r]=c1*bot_Bx[i][j][r]+c2*(-Ez[i+1][j+1][k]+Ez[i+1][j][k]);
	                Hx[i][j][k]=c3*Hx[i][j][k]+c4*(c5*bot_Bx[i][j][r]-c6*temp);

            }
        }
    }
}

void update_bot_Hy_null_Ex()
{
    typ_pola temp;
    ////////
    // Hy //
    ////////

    c1=c1z_H[boundz];
    c2=c2z_H[boundz];

    for(i=0;i<xs;i++)
    {
        c3=cstx(c3x_H,i,xs);
        c4=cstx(c4x_H,i,xs);

	 for(j=0;j<boundy;j++)
         {
         	c5=c5y_H[j];
         	c6=c6y_H[j];

        	for(k=boundz;k<zk;k++)
        	{
	            	r=k-boundz;

	                temp=bot_By[i][j][r];
	                bot_By[i][j][r]=c1*bot_By[i][j][r]+c2*(Ez[i+1][j+1][k]-Ez[i][j+1][k]);
                	Hy[i][j][k]=c3*Hy[i][j][k]+c4*(c5*bot_By[i][j][r]-c6*temp);

            }
        }
    }
}
void update_bot_Hy_null_Ez()
{
    typ_pola temp;
    ////////
    // Hy //
    ////////

    c1=c1z_H[boundz];
    c2=c2z_H[boundz];

    for(i=0;i<xs;i++)
    {
        c3=cstx(c3x_H,i,xs);
        c4=cstx(c4x_H,i,xs);

	 for(j=0;j<boundy;j++)
         {
         	c5=c5y_H[j];
         	c6=c6y_H[j];

        	for(k=boundz;k<zk;k++)
        	{
	            	r=k-boundz;

	                temp=bot_By[i][j][r];
	                bot_By[i][j][r]=c1*bot_By[i][j][r]+c2*(-Ex[i][j+1][k+1]+Ex[i][j+1][k]);
                	Hy[i][j][k]=c3*Hy[i][j][k]+c4*(c5*bot_By[i][j][r]-c6*temp);

            }
        }
    }
}
void update_bot_Hz_null_Ex()
{
    typ_pola temp;

    ////////
    // Hz //
    ////////

	// w centralnej czesci bottom i top nie musimy zapamietywac Bz i Dz

    c5=c5z_H[boundz];
    c6=c6z_H[boundz];
    c1=c1z_H[boundz];
    c2=c2z_H[boundz];

    for(j=0;j<boundy;j++)
    {
    	c3=c3y_H[j];
     	c4=c4y_H[j];

		for(i=boundx;i<xk;i++)
  		{
    		for(k=boundz;k<zk-1;k++)
      		{
                Hz[i][j][k]=c3*Hz[i][j][k]+c4*c5*c2*(-Ey[i+1][j][k+1]+Ey[i][j][k+1]);
            }
        }
    }

	// po bokach (left i right) musimy pamietac wszystko !!

	c5=c5x_H[boundx];
    c6=c6x_H[boundx];

    for(i=0;i<boundx;i++)
    {
        c1=c1x_H[i];
        c2=c2x_H[i];

		q=xs-i-1;

  		for(j=0;j<boundy;j++)
    	{

            c3=c3y_H[j];
            c4=c4y_H[j];

			for(k=boundz;k<zk-1;k++)
        	{
            	r=k-boundz;

				// update po lewej
                temp=bot_Bz_lef[i][j][r];
                bot_Bz_lef[i][j][r]=c1*bot_Bz_lef[i][j][r]+c2*(-Ey[i+1][j][k+1]+Ey[i][j][k+1]);
                Hz[i][j][k]=c3*Hz[i][j][k]+c4*(c5*bot_Bz_lef[i][j][r]-c6*temp);

				// update po prawej
                temp=bot_Bz_rig[i][j][r];
                bot_Bz_rig[i][j][r]=c1*bot_Bz_rig[i][j][r]+c2*(-Ey[q+1][j][k+1]+Ey[q][j][k+1]);
                Hz[q][j][k]=c3*Hz[q][j][k]+c4*(c5*bot_Bz_rig[i][j][r]-c6*temp);
            }
        }
    }
}
void update_bot_Hz_null_Ey()
{
    typ_pola temp;

    ////////
    // Hz //
    ////////

	// w centralnej czesci bottom i top nie musimy zapamietywac Bz i Dz

    c5=c5z_H[boundz];
    c6=c6z_H[boundz];
    c1=c1z_H[boundz];
    c2=c2z_H[boundz];

    for(j=0;j<boundy;j++)
    {
    	c3=c3y_H[j];
     	c4=c4y_H[j];

		for(i=boundx;i<xk;i++)
  		{
    		for(k=boundz;k<zk-1;k++)
      		{
                Hz[i][j][k]=c3*Hz[i][j][k]+c4*c5*c2*(Ex[i][j+1][k+1]-Ex[i][j][k+1]);
            }
        }
    }

	// po bokach (left i right) musimy pamietac wszystko !!

	c5=c5x_H[boundx];
    c6=c6x_H[boundx];

    for(i=0;i<boundx;i++)
    {
        c1=c1x_H[i];
        c2=c2x_H[i];

		q=xs-i-1;

  		for(j=0;j<boundy;j++)
    	{

            c3=c3y_H[j];
            c4=c4y_H[j];

			for(k=boundz;k<zk-1;k++)
        	{
            	r=k-boundz;

				// update po lewej
                temp=bot_Bz_lef[i][j][r];
                bot_Bz_lef[i][j][r]=c1*bot_Bz_lef[i][j][r]+c2*(Ex[i][j+1][k+1]-Ex[i][j][k+1]);
                Hz[i][j][k]=c3*Hz[i][j][k]+c4*(c5*bot_Bz_lef[i][j][r]-c6*temp);

				// update po prawej
                temp=bot_Bz_rig[i][j][r];
                bot_Bz_rig[i][j][r]=c1*bot_Bz_rig[i][j][r]+c2*(Ex[q][j+1][k+1]-Ex[q][j][k+1]);
                Hz[q][j][k]=c3*Hz[q][j][k]+c4*(c5*bot_Bz_rig[i][j][r]-c6*temp);
            }
        }
    }
}

void update_top_Hx_null_Ez()
{

    typ_pola temp;

    ////////
    // Hx //
    ////////

    c3=c3z_H[bound];
    c4=c4z_H[bound];

    for(i=0;i<xs-1;i++)
    {
        c5=cstx(c5x_H,i,xs-1);
        c6=cstx(c6x_H,i,xs-1);

        for(j=0;j<boundy;j++)
        {

			c1=c1y_H[j];
			c2=c2y_H[j];

			p=ys-j-1;

        	for(k=boundz;k<zk;k++)
        	{
            	r=k-boundz;

                temp=top_Bx[i][j][r];
                top_Bx[i][j][r]=c1*top_Bx[i][j][r]+c2*(Ey[i+1][p][k+1]-Ey[i+1][p][k]);
                Hx[i][p][k]=c3*Hx[i][p][k]+c4*(c5*top_Bx[i][j][r]-c6*temp);
            }
        }
    }
}
void update_top_Hx_null_Ey()
{

    typ_pola temp;

    ////////
    // Hx //
    ////////

    c3=c3z_H[bound];
    c4=c4z_H[bound];

    for(i=0;i<xs-1;i++)
    {
        c5=cstx(c5x_H,i,xs-1);
        c6=cstx(c6x_H,i,xs-1);

        for(j=0;j<boundy;j++)
        {

			c1=c1y_H[j];
			c2=c2y_H[j];

			p=ys-j-1;

        	for(k=boundz;k<zk;k++)
        	{
            	r=k-boundz;

                temp=top_Bx[i][j][r];
                top_Bx[i][j][r]=c1*top_Bx[i][j][r]+c2*(-Ez[i+1][p+1][k]+Ez[i+1][p][k]);
                Hx[i][p][k]=c3*Hx[i][p][k]+c4*(c5*top_Bx[i][j][r]-c6*temp);
            }
        }
    }
}

void update_top_Hy_null_Ex()
{
    typ_pola temp;
    ////////
    // Hy //
    ////////

    c1=c1z_H[boundz];
    c2=c2z_H[boundz];

    for(i=0;i<xs;i++)
    {
        c3=cstx(c3x_H,i,xs);
        c4=cstx(c4x_H,i,xs);

		 for(j=0;j<boundy;j++)
         {
         	c5=c5y_H[j];
         	c6=c6y_H[j];

			p=ys-j-2;

        	for(k=boundz;k<zk;k++)
        	{
            	r=k-boundz;

                temp=top_By[i][j][r];
                top_By[i][j][r]=c1*top_By[i][j][r]+c2*(Ez[i+1][p+1][k]-Ez[i][p+1][k]);
                Hy[i][p][k]=c3*Hy[i][p][k]+c4*(c5*top_By[i][j][r]-c6*temp);
            }
        }
    }
}
void update_top_Hy_null_Ez()
{
    typ_pola temp;
    ////////
    // Hy //
    ////////

    c1=c1z_H[boundz];
    c2=c2z_H[boundz];

    for(i=0;i<xs;i++)
    {
        c3=cstx(c3x_H,i,xs);
        c4=cstx(c4x_H,i,xs);

		 for(j=0;j<boundy;j++)
         {
         	c5=c5y_H[j];
         	c6=c6y_H[j];

			p=ys-j-2;

        	for(k=boundz;k<zk;k++)
        	{
            	r=k-boundz;

                temp=top_By[i][j][r];
                top_By[i][j][r]=c1*top_By[i][j][r]+c2*(-Ex[i][p+1][k+1]+Ex[i][p+1][k]);
                Hy[i][p][k]=c3*Hy[i][p][k]+c4*(c5*top_By[i][j][r]-c6*temp);
            }
        }
    }
}
void update_top_Hz_null_Ey()
{
    ////////
    // Hz //
    ////////
    typ_pola temp;
	// w centralnej czesci bottom i top nie musimy zapamietywac Bz i Dz

    c5=c5z_H[boundz];
    c6=c6z_H[boundz];
    c1=c1z_H[boundz];
    c2=c2z_H[boundz];

    for(j=0;j<boundy;j++)
    {
    	c3=c3y_H[j];
     	c4=c4y_H[j];

		p=ys-j-1;

		for(i=boundx;i<xk;i++)
  		{
    		for(k=boundz;k<zk-1;k++)
      		{
                Hz[i][p][k]=c3*Hz[i][p][k]+c4*c5*c2*(Ex[i][p+1][k+1]-Ex[i][p][k+1]);
            }
        }
    }

	// po bokach (left i right) musimy pamietac wszystko !!

	//c5=c5_H[bound];
    //c6=c6_H[bound];

    for(i=0;i<boundx;i++)
    {
        c1=c1x_H[i];
        c2=c2x_H[i];

		q=xs-i-1;

  		for(j=0;j<boundy;j++)
    	{

            c3=c3y_H[j];
            c4=c4y_H[j];

			p=ys-j-1;

			for(k=boundz;k<zk-1;k++)
        	{
            	r=k-boundz;

				// update po lewej
                temp=top_Bz_lef[i][j][r];
                top_Bz_lef[i][j][r]=c1*top_Bz_lef[i][j][r]+c2*(Ex[i][p+1][k+1]-Ex[i][p][k+1]);
                Hz[i][p][k]=c3*Hz[i][p][k]+c4*(c5*top_Bz_lef[i][j][r]-c6*temp);

				// update po prawej
                temp=top_Bz_rig[i][j][r];
                top_Bz_rig[i][j][r]=c1*top_Bz_rig[i][j][r]+c2*(Ex[q][p+1][k+1]-Ex[q][p][k+1]);
                Hz[q][p][k]=c3*Hz[q][p][k]+c4*(c5*top_Bz_rig[i][j][r]-c6*temp);
            }
        }
    }

}
void update_top_Hz_null_Ex()
{
    ////////
    // Hz //
    ////////
    typ_pola temp;
	// w centralnej czesci bottom i top nie musimy zapamietywac Bz i Dz

    c5=c5z_H[boundz];
    c6=c6z_H[boundz];
    c1=c1z_H[boundz];
    c2=c2z_H[boundz];

    for(j=0;j<boundy;j++)
    {
    	c3=c3y_H[j];
     	c4=c4y_H[j];

		p=ys-j-1;

		for(i=boundx;i<xk;i++)
  		{
    		for(k=boundz;k<zk-1;k++)
      		{
                Hz[i][p][k]=c3*Hz[i][p][k]+c4*c5*c2*(-Ey[i+1][p][k+1]+Ey[i][p][k+1]);
            }
        }
    }

	// po bokach (left i right) musimy pamietac wszystko !!

	//c5=c5_H[bound];
    //c6=c6_H[bound];

    for(i=0;i<boundx;i++)
    {
        c1=c1x_H[i];
        c2=c2x_H[i];

		q=xs-i-1;

  		for(j=0;j<boundy;j++)
    	{

            c3=c3y_H[j];
            c4=c4y_H[j];

			p=ys-j-1;

			for(k=boundz;k<zk-1;k++)
        	{
            	r=k-boundz;

				// update po lewej
                temp=top_Bz_lef[i][j][r];
                top_Bz_lef[i][j][r]=c1*top_Bz_lef[i][j][r]+c2*(-Ey[i+1][p][k+1]+Ey[i][p][k+1]);
                Hz[i][p][k]=c3*Hz[i][p][k]+c4*(c5*top_Bz_lef[i][j][r]-c6*temp);

				// update po prawej
                temp=top_Bz_rig[i][j][r];
                top_Bz_rig[i][j][r]=c1*top_Bz_rig[i][j][r]+c2*(-Ey[q+1][p][k+1]+Ey[q][p][k+1]);
                Hz[q][p][k]=c3*Hz[q][p][k]+c4*(c5*top_Bz_rig[i][j][r]-c6*temp);
            }
        }
    }

}


void update_bot_Ex_null_Hy()
{

    typ_pola temp;

    ////////
    // Ex //
    ////////

    c3=c3z_E[boundz];
    c4=c4z_E[boundz];

    for(i=0;i<xs;i++)
    {
        c5=cstx(c5x_E,i,xs);
        c6=cstx(c6x_E,i,xs);

        for(j=1;j<boundy;j++)
        {
        	c1=c1y_E[j];
        	c2=c2y_E[j];

        	for(k=boundz;k<zk+1;k++)
        	{
            	r=k-boundz;

                temp=bot_Dx[i][j][r];
                bot_Dx[i][j][r]=c1*bot_Dx[i][j][r]+c2*(Hz[i][j][k-1]-Hz[i][j-1][k-1]);
                Ex[i][j][k]=c3*Ex[i][j][k]+c4*(c5*bot_Dx[i][j][r]-c6*temp);
            }
        }
    }
}
void update_bot_Ex_null_Hz()
{

    typ_pola temp;

    ////////
    // Ex //
    ////////

    c3=c3z_E[boundz];
    c4=c4z_E[boundz];

    for(i=0;i<xs;i++)
    {
        c5=cstx(c5x_E,i,xs);
        c6=cstx(c6x_E,i,xs);

        for(j=1;j<boundy;j++)
        {
        	c1=c1y_E[j];
        	c2=c2y_E[j];

        	for(k=boundz;k<zk+1;k++)
        	{
            	r=k-boundz;

                temp=bot_Dx[i][j][r];
                bot_Dx[i][j][r]=c1*bot_Dx[i][j][r]+c2*(-Hy[i][j-1][k]+Hy[i][j-1][k-1]);
                Ex[i][j][k]=c3*Ex[i][j][k]+c4*(c5*bot_Dx[i][j][r]-c6*temp);
            }
        }
    }
}
void update_bot_Ey_null_Hz()
{

    typ_pola temp;

    ////////
    // Ey //
    ////////

    c1=c1z_E[boundz];
    c2=c2z_E[boundz];

    for(i=1;i<xs;i++)
    {
        c3=cstx(c3x_E,i,xs+1);
        c4=cstx(c4x_E,i,xs+1);

        for(j=0;j<boundy;j++)
        {
        	c5=c5y_E[j];
        	c6=c6y_E[j];

        	for(k=boundz;k<zk+1;k++)
        	{
            	r=k-boundz;

                temp=bot_Dy[i][j][r];
                bot_Dy[i][j][r]=c1*bot_Dy[i][j][r]+c2*(Hx[i-1][j][k]-Hx[i-1][j][k-1]);
                Ey[i][j][k]=c3*Ey[i][j][k]+c4*(c5*bot_Dy[i][j][r]-c6*temp);

            }
        }
    }
}
void update_bot_Ey_null_Hx()
{

    typ_pola temp;

    ////////
    // Ey //
    ////////

    c1=c1z_E[boundz];
    c2=c2z_E[boundz];

    for(i=1;i<xs;i++)
    {
        c3=cstx(c3x_E,i,xs+1);
        c4=cstx(c4x_E,i,xs+1);

        for(j=0;j<boundy;j++)
        {
        	c5=c5y_E[j];
        	c6=c6y_E[j];

        	for(k=boundz;k<zk+1;k++)
        	{
            	r=k-boundz;

                temp=bot_Dy[i][j][r];
                bot_Dy[i][j][r]=c1*bot_Dy[i][j][r]+c2*(-Hz[i][j][k-1]+Hz[i-1][j][k-1]);
                Ey[i][j][k]=c3*Ey[i][j][k]+c4*(c5*bot_Dy[i][j][r]-c6*temp);

            }
        }
    }
}
void update_bot_Ez_null_Hx()
{

    typ_pola temp;
    ////////
    // Ez //
    ////////
	// w centralnej czesci bottom i top nie musimy zapamietywac Bz i Dz

    c5=c5z_E[boundz];
    c6=c6z_E[boundz];

	c1=c1z_E[boundz];
 	c2=c2z_E[boundz];

	for(j=1;j<boundy;j++)
	{
		c3=c3y_E[j];
		c4=c4y_E[j];

		for(i=boundx;i<=xk;i++)
    	{
        	for(k=boundz;k<zk;k++)
        	{
                Ez[i][j][k]=c3*Ez[i][j][k]+c4*c5*c2*(Hy[i][j-1][k]-Hy[i-1][j-1][k]);
            }
        }
    }

	// na brzegach ( left i right) musimy pamietac wszystko

	//c5=c5_E[bound];
    //c6=c6_E[bound];

	for(i=1;i<boundx;i++)
    {
        c1=c1x_E[i];
        c2=c2x_E[i];

		q=xs-i;

  		for(j=1;j<boundy;j++)
    	{
			c3=c3y_E[j];
			c4=c4y_E[j];

        	for(k=boundz;k<zk;k++)
        	{
            	r=k-boundz;

				// update po lewej
                temp=bot_Dz_lef[i][j][r];
                bot_Dz_lef[i][j][r]=c1*bot_Dz_lef[i][j][r]+c2*(Hy[i][j-1][k]-Hy[i-1][j-1][k]);
                Ez[i][j][k]=c3*Ez[i][j][k]+c4*(c5*bot_Dz_lef[i][j][r]-c6*temp);

				// update po prawej
                temp=bot_Dz_rig[i][j][r];
                bot_Dz_rig[i][j][r]=c1*bot_Dz_rig[i][j][r]+c2*(Hy[q][j-1][k]-Hy[q-1][j-1][k]);
                Ez[q][j][k]=c3*Ez[q][j][k]+c4*(c5*bot_Dz_rig[i][j][r]-c6*temp);

			}
        }
    }

}
void update_bot_Ez_null_Hy()
{

    typ_pola temp;
    ////////
    // Ez //
    ////////
	// w centralnej czesci bottom i top nie musimy zapamietywac Bz i Dz

    c5=c5z_E[boundz];
    c6=c6z_E[boundz];

	c1=c1z_E[boundz];
 	c2=c2z_E[boundz];

	for(j=1;j<boundy;j++)
	{
		c3=c3y_E[j];
		c4=c4y_E[j];

		for(i=boundx;i<=xk;i++)
    	{
        	for(k=boundz;k<zk;k++)
        	{
                Ez[i][j][k]=c3*Ez[i][j][k]+c4*c5*c2*(-Hx[i-1][j][k]+Hx[i-1][j-1][k]);
            }
        }
    }

	// na brzegach ( left i right) musimy pamietac wszystko

	//c5=c5_E[bound];
    //c6=c6_E[bound];

	for(i=1;i<boundx;i++)
    {
        c1=c1x_E[i];
        c2=c2x_E[i];

		q=xs-i;

  		for(j=1;j<boundy;j++)
    	{
			c3=c3y_E[j];
			c4=c4y_E[j];

        	for(k=boundz;k<zk;k++)
        	{
            	r=k-boundz;

				// update po lewej
                temp=bot_Dz_lef[i][j][r];
                bot_Dz_lef[i][j][r]=c1*bot_Dz_lef[i][j][r]+c2*(-Hx[i-1][j][k]+Hx[i-1][j-1][k]);
                Ez[i][j][k]=c3*Ez[i][j][k]+c4*(c5*bot_Dz_lef[i][j][r]-c6*temp);

				// update po prawej
                temp=bot_Dz_rig[i][j][r];
                bot_Dz_rig[i][j][r]=c1*bot_Dz_rig[i][j][r]+c2*(-Hx[q-1][j][k]+Hx[q-1][j-1][k]);
                Ez[q][j][k]=c3*Ez[q][j][k]+c4*(c5*bot_Dz_rig[i][j][r]-c6*temp);

			}
        }
    }

}

void update_top_Ex_null_Hy()
{

    typ_pola temp;
 
    ////////
    // Ex //
    ////////

    c3=c3z_E[boundz];
    c4=c4z_E[boundz];

    for(i=0;i<xs;i++)
    {
        c5=cstx(c5x_E,i,xs);
        c6=cstx(c6x_E,i,xs);

        for(j=1;j<boundy;j++)
        {
        	c1=c1y_E[j];
        	c2=c2y_E[j];

			p=ys-j;

        	for(k=boundz;k<zk+1;k++)
        	{
            	r=k-boundz;

                temp=top_Dx[i][j][r];
                top_Dx[i][j][r]=c1*top_Dx[i][j][r]+c2*(Hz[i][p][k-1]-Hz[i][p-1][k-1]);
                Ex[i][p][k]=c3*Ex[i][p][k]+c4*(c5*top_Dx[i][j][r]-c6*temp);
            }
        }
    }
}
void update_top_Ex_null_Hz()
{

    typ_pola temp;
 
    ////////
    // Ex //
    ////////

    c3=c3z_E[boundz];
    c4=c4z_E[boundz];

    for(i=0;i<xs;i++)
    {
        c5=cstx(c5x_E,i,xs);
        c6=cstx(c6x_E,i,xs);

        for(j=1;j<boundy;j++)
        {
        	c1=c1y_E[j];
        	c2=c2y_E[j];

			p=ys-j;

        	for(k=boundz;k<zk+1;k++)
        	{
            	r=k-boundz;

                temp=top_Dx[i][j][r];
                top_Dx[i][j][r]=c1*top_Dx[i][j][r]+c2*(-Hy[i][p-1][k]+Hy[i][p-1][k-1]);
                Ex[i][p][k]=c3*Ex[i][p][k]+c4*(c5*top_Dx[i][j][r]-c6*temp);
            }
        }
    }
}

void update_top_Ey_null_Hz()
{

    typ_pola temp;

    ////////
    // Ey //
    ////////

    c1=c1z_E[boundz];
    c2=c2z_E[boundz];

    for(i=1;i<xs;i++)
    {
        c3=cstx(c3x_E,i,xs+1);
        c4=cstx(c4x_E,i,xs+1);

        for(j=0;j<boundy;j++)
        {
        	c5=c5y_E[j];
        	c6=c6y_E[j];

    		p=ys-j-1;

        	for(k=boundz;k<zk+1;k++)
        	{
            	r=k-boundz;

                temp=top_Dy[i][j][r];
                top_Dy[i][j][r]=c1*top_Dy[i][j][r]+c2*(Hx[i-1][p][k]-Hx[i-1][p][k-1]);
                Ey[i][p][k]=c3*Ey[i][p][k]+c4*(c5*top_Dy[i][j][r]-c6*temp);
            }
        }
    }
}
void update_top_Ey_null_Hx()
{

    typ_pola temp;

    ////////
    // Ey //
    ////////

    c1=c1z_E[boundz];
    c2=c2z_E[boundz];

    for(i=1;i<xs;i++)
    {
        c3=cstx(c3x_E,i,xs+1);
        c4=cstx(c4x_E,i,xs+1);

        for(j=0;j<boundy;j++)
        {
        	c5=c5y_E[j];
        	c6=c6y_E[j];

    		p=ys-j-1;

        	for(k=boundz;k<zk+1;k++)
        	{
            	r=k-boundz;

                temp=top_Dy[i][j][r];
                top_Dy[i][j][r]=c1*top_Dy[i][j][r]+c2*(-Hz[i][p][k-1]+Hz[i-1][p][k-1]);
                Ey[i][p][k]=c3*Ey[i][p][k]+c4*(c5*top_Dy[i][j][r]-c6*temp);
            }
        }
    }
}

void update_top_Ez_null_Hx()
{

    typ_pola temp;
    ////////
    // Ez //
    ////////
	// w centralnej czesci bottom i top nie musimy zapamietywac Bz i Dz

    c5=c5z_E[boundz];
    c6=c6z_E[boundz];

	c1=c1z_E[boundz];
 	c2=c2z_E[boundz];

	for(j=1;j<boundy;j++)
	{
		c3=c3y_E[j];
		c4=c4y_E[j];
		p=ys-j;

		for(i=boundx;i<=xk;i++)
    	{
        	for(k=boundz;k<zk;k++)
        	{
                Ez[i][p][k]=c3*Ez[i][p][k]+c4*c5*c2*(Hy[i][p-1][k]-Hy[i-1][p-1][k]);
            }
        }
    }

	// na brzegach ( left i right) musimy pamietac wszystko

	//c5=c5_E[bound];
    //c6=c6_E[bound];

	for(i=1;i<boundx;i++)
    {
        c1=c1x_E[i];
        c2=c2x_E[i];

		q=xs-i;

  		for(j=1;j<boundy;j++)
    	{
			c3=c3y_E[j];
			c4=c4y_E[j];

			p=ys-j;

        	for(k=boundz;k<zk;k++)
        	{
            	r=k-boundz;

				// update po lewej
                temp=top_Dz_lef[i][j][r];
                top_Dz_lef[i][j][r]=c1*top_Dz_lef[i][j][r]+c2*(Hy[i][p-1][k]-Hy[i-1][p-1][k]);
                Ez[i][p][k]=c3*Ez[i][p][k]+c4*(c5*top_Dz_lef[i][j][r]-c6*temp);

				// update po prawej
                temp=top_Dz_rig[i][j][r];
                top_Dz_rig[i][j][r]=c1*top_Dz_rig[i][j][r]+c2*(Hy[q][p-1][k]-Hy[q-1][p-1][k]);
                Ez[q][p][k]=c3*Ez[q][p][k]+c4*(c5*top_Dz_rig[i][j][r]-c6*temp);
            }
        }
    }
}
void update_top_Ez_null_Hy()
{

    typ_pola temp;
    ////////
    // Ez //
    ////////
	// w centralnej czesci bottom i top nie musimy zapamietywac Bz i Dz

    c5=c5z_E[boundz];
    c6=c6z_E[boundz];

	c1=c1z_E[boundz];
 	c2=c2z_E[boundz];

	for(j=1;j<boundy;j++)
	{
		c3=c3y_E[j];
		c4=c4y_E[j];
		p=ys-j;

		for(i=boundx;i<=xk;i++)
    	{
        	for(k=boundz;k<zk;k++)
        	{
                Ez[i][p][k]=c3*Ez[i][p][k]+c4*c5*c2*(-Hx[i-1][p][k]+Hx[i-1][p-1][k]);
            }
        }
    }

	// na brzegach ( left i right) musimy pamietac wszystko

	//c5=c5_E[bound];
    //c6=c6_E[bound];

	for(i=1;i<boundx;i++)
    {
        c1=c1x_E[i];
        c2=c2x_E[i];

		q=xs-i;

  		for(j=1;j<boundy;j++)
    	{
			c3=c3y_E[j];
			c4=c4y_E[j];

			p=ys-j;

        	for(k=boundz;k<zk;k++)
        	{
            	r=k-boundz;

				// update po lewej
                temp=top_Dz_lef[i][j][r];
                top_Dz_lef[i][j][r]=c1*top_Dz_lef[i][j][r]+c2*(-Hx[i-1][p][k]+Hx[i-1][p-1][k]);
                Ez[i][p][k]=c3*Ez[i][p][k]+c4*(c5*top_Dz_lef[i][j][r]-c6*temp);

				// update po prawej
                temp=top_Dz_rig[i][j][r];
                top_Dz_rig[i][j][r]=c1*top_Dz_rig[i][j][r]+c2*(-Hx[q-1][p][k]+Hx[q-1][p-1][k]);
                Ez[q][p][k]=c3*Ez[q][p][k]+c4*(c5*top_Dz_rig[i][j][r]-c6*temp);
            }
        }
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void update_lef_Hx_null_Ez()
{
    typ_pola temp;

    ////////
    // Hx //
    ////////

    c1=c1x_H[boundx];
    c2=c2x_H[boundx];

    c3=c3x_H[boundx];
    c4=c4x_H[boundx];

	for(i=0;i<boundx;i++)
	{
		c5=c5x_H[i];
		c6=c6x_H[i];

    	for(j=boundy;j<yk;j++)
    	{
        	q=j-boundy;

        	for(k=boundz;k<zk;k++)
        	{
            	r=k-boundz;

                temp=lef_Bx[i][q][r];
                lef_Bx[i][q][r]=c1*lef_Bx[i][q][r]+c2*(Ey[i+1][j][k+1]-Ey[i+1][j][k]);
                Hx[i][j][k]=c3*Hx[i][j][k]+c4*(c5*lef_Bx[i][q][r]-c6*temp);

            }
        }
    }
}

void update_lef_Hx_null_Ey()
{
    typ_pola temp;

    ////////
    // Hx //
    ////////

    c1=c1x_H[boundx];
    c2=c2x_H[boundx];

    c3=c3x_H[boundx];
    c4=c4x_H[boundx];

	for(i=0;i<boundx;i++)
	{
		c5=c5x_H[i];
		c6=c6x_H[i];

    	for(j=boundy;j<yk;j++)
    	{
        	q=j-boundy;

        	for(k=boundz;k<zk;k++)
        	{
            	r=k-boundz;

                temp=lef_Bx[i][q][r];
                lef_Bx[i][q][r]=c1*lef_Bx[i][q][r]+c2*(-Ez[i+1][j+1][k]+Ez[i+1][j][k]);
                Hx[i][j][k]=c3*Hx[i][j][k]+c4*(c5*lef_Bx[i][q][r]-c6*temp);

            }
        }
    }
}

void update_lef_Hy_null_Ex()
{
    ////////
    // Hy //
    ////////

	// zauwazmy ,ze nie trzeba przechowywac By,Dy bo c1=1, a c5=c6
	// co oszczedza nam pamiec i upraszcza wzory

    // ustawione przy Hx
    c1=c1x_H[boundx];
    c2=c2x_H[boundx];

	c5=c5x_H[boundx];
    c6=c6x_H[boundx];

	for(i=0;i<boundx;i++)
	{
		c3=c3x_H[i];
		c4=c4x_H[i];

    	for(k=boundz;k<zk;k++)
    	{
        	for(j=boundy;j<yk-1;j++)
        	{
                Hy[i][j][k]=c3*Hy[i][j][k]+c4*c5*c2*(Ez[i+1][j+1][k]-Ez[i][j+1][k]);
            }
        }
    }
}
void update_lef_Hy_null_Ez()
{
    ////////
    // Hy //
    ////////

	// zauwazmy ,ze nie trzeba przechowywac By,Dy bo c1=1, a c5=c6
	// co oszczedza nam pamiec i upraszcza wzory

    // ustawione przy Hx
    c1=c1x_H[boundx];
    c2=c2x_H[boundx];

	c5=c5x_H[boundx];
    c6=c6x_H[boundx];

	for(i=0;i<boundx;i++)
	{
		c3=c3x_H[i];
		c4=c4x_H[i];

    	for(k=boundz;k<zk;k++)
    	{
        	for(j=boundy;j<yk-1;j++)
        	{
                Hy[i][j][k]=c3*Hy[i][j][k]+c4*c5*c2*(-Ex[i][j+1][k+1]+Ex[i][j+1][k]);
            }
        }
    }
}
void update_lef_Hz_null_Ey()
{
    typ_pola temp;
    ////////
    // Hz //
    ////////

    c3=c3x_H[boundx];
    c4=c4x_H[boundx];

	// ustawione przy Hy
 	c5=c5x_H[boundx];
    c6=c6x_H[boundx];

	for(i=0;i<boundx;i++)
	{
		c1=c1x_H[i];
		c2=c2x_H[i];

		for(j=boundy;j<yk;j++)
    	{
        	q=j-boundy;

        	for(k=boundz;k<zk-1;k++)
        	{
            	r=k-boundz;

                temp=lef_Bz[i][q][r];
                lef_Bz[i][q][r]=c1*lef_Bz[i][q][r]+c2*(Ex[i][j+1][k+1]-Ex[i][j][k+1]);
                Hz[i][j][k]=c3*Hz[i][j][k]+c4*(c5*lef_Bz[i][q][r]-c6*temp);

            }
        }
    }

}
void update_lef_Hz_null_Ex()
{
    typ_pola temp;
    ////////
    // Hz //
    ////////

    c3=c3x_H[boundx];
    c4=c4x_H[boundx];

	// ustawione przy Hy
 	c5=c5x_H[boundx];
    c6=c6x_H[boundx];

	for(i=0;i<boundx;i++)
	{
		c1=c1x_H[i];
		c2=c2x_H[i];

		for(j=boundy;j<yk;j++)
    	{
        	q=j-boundy;

        	for(k=boundz;k<zk-1;k++)
        	{
            	r=k-boundz;

                temp=lef_Bz[i][q][r];
                lef_Bz[i][q][r]=c1*lef_Bz[i][q][r]+c2*(-Ey[i+1][j][k+1]+Ey[i][j][k+1]);
                Hz[i][j][k]=c3*Hz[i][j][k]+c4*(c5*lef_Bz[i][q][r]-c6*temp);

            }
        }
    }

}

void update_rig_Hx_null_Ez()
{

    typ_pola temp;

    ////////
    // Hx //
    ////////

    c1=c1x_H[boundx];
    c2=c2x_H[boundx];

    c3=c3x_H[boundx];
    c4=c4x_H[boundx];

	for(i=0;i<boundx;i++)
	{
		c5=c5x_H[i];
		c6=c6x_H[i];

		p=xs-i-2;

    	for(j=boundy;j<yk;j++)
    	{
        	q=j-boundy;

        	for(k=boundz;k<zk;k++)
        	{
            	r=k-boundz;

                temp=rig_Bx[i][q][r];
                rig_Bx[i][q][r]=c1*rig_Bx[i][q][r]+c2*(Ey[p+1][j][k+1]-Ey[p+1][j][k]);
                Hx[p][j][k]=c3*Hx[p][j][k]+c4*(c5*rig_Bx[i][q][r]-c6*temp);
            }
        }
    }
}

void update_rig_Hx_null_Ey()
{

    typ_pola temp;

    ////////
    // Hx //
    ////////

    c1=c1x_H[boundx];
    c2=c2x_H[boundx];

    c3=c3x_H[boundx];
    c4=c4x_H[boundx];

	for(i=0;i<boundx;i++)
	{
		c5=c5x_H[i];
		c6=c6x_H[i];

		p=xs-i-2;

    	for(j=boundy;j<yk;j++)
    	{
        	q=j-boundy;

        	for(k=boundz;k<zk;k++)
        	{
            	r=k-boundz;

                temp=rig_Bx[i][q][r];
                rig_Bx[i][q][r]=c1*rig_Bx[i][q][r]+c2*(-Ez[p+1][j+1][k]+Ez[p+1][j][k]);
                Hx[p][j][k]=c3*Hx[p][j][k]+c4*(c5*rig_Bx[i][q][r]-c6*temp);
            }
        }
    }
}

void update_rig_Hy_null_Ex()
{
    ////////
    // Hy //
    ////////

	// zauwazmy ,ze nie trzeba przechowywac By,Dy bo c1=1, a c5=c6
	// co oszczedza nam pamiec i upraszcza wzory

    // ustawione przy Hx
    c1=c1x_H[boundx];
    c2=c2x_H[boundx];

	c5=c5x_H[boundx];
    c6=c6x_H[boundx];

	for(i=0;i<boundx;i++)
	{
		c3=c3x_H[i];
		c4=c4x_H[i];

		p=xs-i-1;

    	for(k=boundz;k<zk;k++)
    	{
        	for(j=boundy;j<yk-1;j++)
        	{
                Hy[p][j][k]=c3*Hy[p][j][k]+c4*c5*c2*(Ez[p+1][j+1][k]-Ez[p][j+1][k]);
            }
        }
    }

}
void update_rig_Hy_null_Ez()
{

    ////////
    // Hy //
    ////////

	// zauwazmy ,ze nie trzeba przechowywac By,Dy bo c1=1, a c5=c6
	// co oszczedza nam pamiec i upraszcza wzory

    // ustawione przy Hx
    c1=c1x_H[boundx];
    c2=c2x_H[boundx];

	c5=c5x_H[boundx];
    c6=c6x_H[boundx];

	for(i=0;i<boundx;i++)
	{
		c3=c3x_H[i];
		c4=c4x_H[i];

		p=xs-i-1;

    	for(k=boundz;k<zk;k++)
    	{
        	for(j=boundy;j<yk-1;j++)
        	{
                Hy[p][j][k]=c3*Hy[p][j][k]+c4*c5*c2*(-Ex[p][j+1][k+1]+Ex[p][j+1][k]);
            }
        }
    }

}

void update_rig_Hz_null_Ey()
{

    typ_pola temp;
    ////////
    // Hz //
    ////////

    c3=c3x_H[boundx];
    c4=c4x_H[boundx];

	// ustawione przy Hy
 	c5=c5x_H[boundx];
    c6=c6x_H[boundx];

	for(i=0;i<boundx;i++)
	{
		c1=c1x_H[i];
		c2=c2x_H[i];

		p=xs-i-1;

		for(j=boundy;j<yk;j++)
    	{
        	q=j-boundy;

        	for(k=boundz;k<zk-1;k++)
        	{
            	r=k-boundz;

                temp=rig_Bz[i][q][r];
                rig_Bz[i][q][r]=c1*rig_Bz[i][q][r]+c2*(Ex[p][j+1][k+1]-Ex[p][j][k+1]);
                Hz[p][j][k]=c3*Hz[p][j][k]+c4*(c5*rig_Bz[i][q][r]-c6*temp);
            }
        }
    }

}
void update_rig_Hz_null_Ex()
{

    typ_pola temp;
    ////////
    // Hz //
    ////////

    c3=c3x_H[boundx];
    c4=c4x_H[boundx];

	// ustawione przy Hy
 	c5=c5x_H[boundx];
    c6=c6x_H[boundx];

	for(i=0;i<boundx;i++)
	{
		c1=c1x_H[i];
		c2=c2x_H[i];

		p=xs-i-1;

		for(j=boundy;j<yk;j++)
    	{
        	q=j-boundy;

        	for(k=boundz;k<zk-1;k++)
        	{
            	r=k-boundz;

                temp=rig_Bz[i][q][r];
                rig_Bz[i][q][r]=c1*rig_Bz[i][q][r]+c2*(-Ey[p+1][j][k+1]+Ey[p][j][k+1]);
                Hz[p][j][k]=c3*Hz[p][j][k]+c4*(c5*rig_Bz[i][q][r]-c6*temp);
            }
        }
    }

}

void update_lef_Ex_null_Hy()
{

    typ_pola temp;

    ////////
    // Ex //
    ////////

    c1=c1x_E[boundx];
    c2=c2x_E[boundx];

    c3=c3x_E[boundx];
    c4=c4x_E[boundx];

	for(i=0;i<boundx;i++)
	{
		c5=c5x_E[i];
    	c6=c6x_E[i];

    	for(j=boundy;j<yk+1;j++)
    	{
        	q=j-boundy;

        	for(k=boundz;k<zk+1;k++)
        	{
            	r=k-boundz;

                temp=lef_Dx[i][q][r];
                lef_Dx[i][q][r]=c1*lef_Dx[i][q][r]+c2*(Hz[i][j][k-1]-Hz[i][j-1][k-1]);
                Ex[i][j][k]=c3*Ex[i][j][k]+c4*(c5*lef_Dx[i][q][r]-c6*temp);

            }
        }
    }
}
void update_lef_Ex_null_Hz()
{

    typ_pola temp;

    ////////
    // Ex //
    ////////

    c1=c1x_E[boundx];
    c2=c2x_E[boundx];

    c3=c3x_E[boundx];
    c4=c4x_E[boundx];

	for(i=0;i<boundx;i++)
	{
		c5=c5x_E[i];
    	c6=c6x_E[i];

    	for(j=boundy;j<yk+1;j++)
    	{
        	q=j-boundy;

        	for(k=boundz;k<zk+1;k++)
        	{
            	r=k-boundz;

                temp=lef_Dx[i][q][r];
                lef_Dx[i][q][r]=c1*lef_Dx[i][q][r]+c2*(-Hy[i][j-1][k]+Hy[i][j-1][k-1]);
                Ex[i][j][k]=c3*Ex[i][j][k]+c4*(c5*lef_Dx[i][q][r]-c6*temp);

            }
        }
    }
}

void update_lef_Ey_null_Hz()
{

    ////////
    // Ey //
    ////////

	// zauwazmy ,ze nie trzeba przechowywac By,Dy bo c1=1, a c5=c6
	// co oszczedza nam pamiec i upraszcza wzory

    //ustawione przy Ex
    c1=c1x_E[boundx];
    c2=c2x_E[boundx];

    c5=c5x_E[boundx];
    c6=c6x_E[boundx];

	for(i=1;i<boundx;i++)
	{
		c3=c3x_E[i];
		c4=c4x_E[i];

    	for(k=boundz;k<zk+1;k++)
    	{
        	for(j=boundy;j<yk;j++)
        	{
                Ey[i][j][k]=c3*Ey[i][j][k]+c4*c5*c2*(Hx[i-1][j][k]-Hx[i-1][j][k-1]);
            }
        }
    }
}
void update_lef_Ey_null_Hx()
{

    ////////
    // Ey //
    ////////

	// zauwazmy ,ze nie trzeba przechowywac By,Dy bo c1=1, a c5=c6
	// co oszczedza nam pamiec i upraszcza wzory

    //ustawione przy Ex
    c1=c1x_E[boundx];
    c2=c2x_E[boundx];

    c5=c5x_E[boundx];
    c6=c6x_E[boundx];

	for(i=1;i<boundx;i++)
	{
		c3=c3x_E[i];
		c4=c4x_E[i];

    	for(k=boundz;k<zk+1;k++)
    	{
        	for(j=boundy;j<yk;j++)
        	{
                Ey[i][j][k]=c3*Ey[i][j][k]+c4*c5*c2*(-Hz[i][j][k-1]+Hz[i-1][j][k-1]);
            }
        }
    }
}

void update_lef_Ez_null_Hx()
{

    typ_pola temp;
    ////////
    // Ez //
    ////////

    // ustawione przy Ey
    c5=c5x_E[boundx];
    c6=c6x_E[boundx];

    c3=c3x_E[boundx];
    c4=c4x_E[boundx];

	for(i=1;i<boundx;i++)
	{
		c1=c1x_E[i];
		c2=c2x_E[i];

    	for(j=boundy;j<yk+1;j++)
    	{
        	q=j-boundy;

        	for(k=boundz;k<zk;k++)
        	{
            	r=k-boundz;

                temp=lef_Dz[i][q][r];
                lef_Dz[i][q][r]=c1*lef_Dz[i][q][r]+c2*(Hy[i][j-1][k]-Hy[i-1][j-1][k]);
                Ez[i][j][k]=c3*Ez[i][j][k]+c4*(c5*lef_Dz[i][q][r]-c6*temp);
            }
        }
    }
}

void update_lef_Ez_null_Hy()
{

    typ_pola temp;
    ////////
    // Ez //
    ////////

    // ustawione przy Ey
    c5=c5x_E[boundx];
    c6=c6x_E[boundx];

    c3=c3x_E[boundx];
    c4=c4x_E[boundx];

	for(i=1;i<boundx;i++)
	{
		c1=c1x_E[i];
		c2=c2x_E[i];

    	for(j=boundy;j<yk+1;j++)
    	{
        	q=j-boundy;

        	for(k=boundz;k<zk;k++)
        	{
            	r=k-boundz;

                temp=lef_Dz[i][q][r];
                lef_Dz[i][q][r]=c1*lef_Dz[i][q][r]+c2*(-Hx[i-1][j][k]+Hx[i-1][j-1][k]);
                Ez[i][j][k]=c3*Ez[i][j][k]+c4*(c5*lef_Dz[i][q][r]-c6*temp);
            }
        }
    }
}

void update_rig_Ex_null_Hy()
{

    typ_pola temp;

    ////////
    // Ex //
    ////////

    c1=c1x_E[boundx];
    c2=c2x_E[boundx];

    c3=c3x_E[boundx];
    c4=c4x_E[boundx];

	for(i=0;i<boundx;i++)
	{
		c5=c5x_E[i];
    	c6=c6x_E[i];

		p=xs-i-1;

    	for(j=boundy;j<yk+1;j++)
    	{
        	q=j-boundy;

        	for(k=boundz;k<zk+1;k++)
        	{
            	r=k-boundz;

                temp=rig_Dx[i][q][r];
                rig_Dx[i][q][r]=c1*rig_Dx[i][q][r]+c2*(Hz[p][j][k-1]-Hz[p][j-1][k-1]);
                Ex[p][j][k]=c3*Ex[p][j][k]+c4*(c5*rig_Dx[i][q][r]-c6*temp);
            }
        }
    }
}

void update_rig_Ex_null_Hz()
{

    typ_pola temp;

    ////////
    // Ex //
    ////////

    c1=c1x_E[boundx];
    c2=c2x_E[boundx];

    c3=c3x_E[boundx];
    c4=c4x_E[boundx];

	for(i=0;i<boundx;i++)
	{
		c5=c5x_E[i];
    	c6=c6x_E[i];

		p=xs-i-1;

    	for(j=boundy;j<yk+1;j++)
    	{
        	q=j-boundy;

        	for(k=boundz;k<zk+1;k++)
        	{
            	r=k-boundz;

                temp=rig_Dx[i][q][r];
                rig_Dx[i][q][r]=c1*rig_Dx[i][q][r]+c2*(-Hy[p][j-1][k]+Hy[p][j-1][k-1]);
                Ex[p][j][k]=c3*Ex[p][j][k]+c4*(c5*rig_Dx[i][q][r]-c6*temp);
            }
        }
    }
}

void update_rig_Ey_null_Hz()
{

    ////////
    // Ey //
    ////////

	// zauwazmy ,ze nie trzeba przechowywac By,Dy bo c1=1, a c5=c6
	// co oszczedza nam pamiec i upraszcza wzory

    //ustawione przy Ex
    c1=c1x_E[boundx];
    c2=c2x_E[boundx];

    c5=c5x_E[boundx];
    c6=c6x_E[boundx];

	for(i=1;i<boundx;i++)
	{
		c3=c3x_E[i];
		c4=c4x_E[i];
		p=xs-i;

    	for(k=boundz;k<zk+1;k++)
    	{
        	for(j=boundy;j<yk;j++)
        	{
                Ey[p][j][k]=c3*Ey[p][j][k]+c4*c5*c2*(Hx[p-1][j][k]-Hx[p-1][j][k-1]);
            }
        }
    }
}

void update_rig_Ey_null_Hx()
{

    ////////
    // Ey //
    ////////

	// zauwazmy ,ze nie trzeba przechowywac By,Dy bo c1=1, a c5=c6
	// co oszczedza nam pamiec i upraszcza wzory

    //ustawione przy Ex
    c1=c1x_E[boundx];
    c2=c2x_E[boundx];

    c5=c5x_E[boundx];
    c6=c6x_E[boundx];

	for(i=1;i<boundx;i++)
	{
		c3=c3x_E[i];
		c4=c4x_E[i];
		p=xs-i;

    	for(k=boundz;k<zk+1;k++)
    	{
        	for(j=boundy;j<yk;j++)
        	{
                Ey[p][j][k]=c3*Ey[p][j][k]+c4*c5*c2*(-Hz[p][j][k-1]+Hz[p-1][j][k-1]);
            }
        }
    }
}

void update_rig_Ez_null_Hx()
{

    typ_pola temp;
    ////////
    // Ez //
    ////////

    // ustawione przy Ey
    c5=c5x_E[boundx];
    c6=c6x_E[boundx];

    c3=c3x_E[boundx];
    c4=c4x_E[boundx];

	for(i=1;i<boundx;i++)
	{
		c1=c1x_E[i];
		c2=c2x_E[i];

  		p=xs-i;

    	for(j=boundy;j<yk+1;j++)
    	{
        	q=j-boundy;

        	for(k=boundz;k<zk;k++)
        	{
            	r=k-boundz;

                temp=rig_Dz[i][q][r];
                rig_Dz[i][q][r]=c1*rig_Dz[i][q][r]+c2*(Hy[p][j-1][k]-Hy[p-1][j-1][k]);
                Ez[p][j][k]=c3*Ez[p][j][k]+c4*(c5*rig_Dz[i][q][r]-c6*temp);
            }
        }
    }
}

void update_rig_Ez_null_Hy()
{

    typ_pola temp;
    ////////
    // Ez //
    ////////

    // ustawione przy Ey
    c5=c5x_E[boundx];
    c6=c6x_E[boundx];

    c3=c3x_E[boundx];
    c4=c4x_E[boundx];

	for(i=1;i<boundx;i++)
	{
		c1=c1x_E[i];
		c2=c2x_E[i];

  		p=xs-i;

    	for(j=boundy;j<yk+1;j++)
    	{
        	q=j-boundy;

        	for(k=boundz;k<zk;k++)
        	{
            	r=k-boundz;

                temp=rig_Dz[i][q][r];
                rig_Dz[i][q][r]=c1*rig_Dz[i][q][r]+c2*(-Hx[p-1][j][k]+Hx[p-1][j-1][k]);
                Ez[p][j][k]=c3*Ez[p][j][k]+c4*(c5*rig_Dz[i][q][r]-c6*temp);
            }
        }
    }
}

///////////////

void update_bac_Hx_null_Ez()
{

    typ_pola temp;

    ////////
    // Hx //
    ////////

	// update na gorze i dole
	for(i=0;i<xs-1;i++)
	{
		c5=cstx(c5x_H,i,xs-1);
		c6=cstx(c6x_H,i,xs-1);

    	for(j=0;j<boundy;j++)
    	{
        	c1=c1y_H[j];
        	c2=c2y_H[j];

			q=ys-j-1;

            for(k=0;k<boundz;k++)
            {
                c3=c3z_H[k];
                c4=c4z_H[k];

				// update na dole
                temp=bac_Bx_bot[i][j][k];
                bac_Bx_bot[i][j][k]=c1*bac_Bx_bot[i][j][k]+c2*(Ey[i+1][j][k+1]-Ey[i+1][j][k]);
                Hx[i][j][k]=c3*Hx[i][j][k]+c4*(c5*bac_Bx_bot[i][j][k]-c6*temp);

				// update na gorze
                temp=bac_Bx_top[i][j][k];
                bac_Bx_top[i][j][k]=c1*bac_Bx_top[i][j][k]+c2*(Ey[i+1][q][k+1]-Ey[i+1][q][k]);
                Hx[i][q][k]=c3*Hx[i][q][k]+c4*(c5*bac_Bx_top[i][j][k]-c6*temp);

            }
        }
    }

	// update pomiedzy niebem i pieklem ( gora i dolem)

	c1=c1z_H[boundz];
 	c2=c2z_H[boundz];

	c5=c5z_H[boundz];
	c6=c6z_H[boundz];

	// update centrum
	for(k=0;k<boundz;k++)
	{
		c3=c3z_H[k];
		c4=c4z_H[k];

		for(j=boundy;j<yk;j++)
    	{
        	for(i=boundx;i<xk-1;i++)
        	{
                Hx[i][j][k]=c3*Hx[i][j][k]+c4*c5*c2*(Ey[i+1][j][k+1]-Ey[i+1][j][k]);
            }
        }
    }


	// update lewicy i prawicy
	//c1=c1_H[bound];
 	//c2=c2_H[bound];

	for(i=0;i<boundx;i++)
	{
		c5=c5x_H[i];
		c6=c6x_H[i];

		q=xs-i-2;

		for(k=0;k<boundz;k++)
		{
			c3=c3z_H[k];
			c4=c4z_H[k];


    		for(j=boundy;j<yk;j++)
    		{
     			r=j-boundy;

				// update left
                temp=bac_Bx_lef[i][r][k];
                bac_Bx_lef[i][r][k]=c1*bac_Bx_lef[i][r][k]+c2*(Ey[i+1][j][k+1]-Ey[i+1][j][k]);
                Hx[i][j][k]=c3*Hx[i][j][k]+c4*(c5*bac_Bx_lef[i][r][k]-c6*temp);

				// update rig
				temp=bac_Bx_rig[i][r][k];
                bac_Bx_rig[i][r][k]=c1*bac_Bx_rig[i][r][k]+c2*(Ey[q+1][j][k+1]-Ey[q+1][j][k]);
                Hx[q][j][k]=c3*Hx[q][j][k]+c4*(c5*bac_Bx_rig[i][r][k]-c6*temp);

            }
        }
    }
}

void update_bac_Hx_null_Ey()
{

    typ_pola temp;

    ////////
    // Hx //
    ////////

	// update na gorze i dole
	for(i=0;i<xs-1;i++)
	{
		c5=cstx(c5x_H,i,xs-1);
		c6=cstx(c6x_H,i,xs-1);

    	for(j=0;j<boundy;j++)
    	{
        	c1=c1y_H[j];
        	c2=c2y_H[j];

			q=ys-j-1;

            for(k=0;k<boundz;k++)
            {
                c3=c3z_H[k];
                c4=c4z_H[k];

				// update na dole
                temp=bac_Bx_bot[i][j][k];
                bac_Bx_bot[i][j][k]=c1*bac_Bx_bot[i][j][k]+c2*(-Ez[i+1][j+1][k]+Ez[i+1][j][k]);
                Hx[i][j][k]=c3*Hx[i][j][k]+c4*(c5*bac_Bx_bot[i][j][k]-c6*temp);

				// update na gorze
                temp=bac_Bx_top[i][j][k];
                bac_Bx_top[i][j][k]=c1*bac_Bx_top[i][j][k]+c2*(-Ez[i+1][q+1][k]+Ez[i+1][q][k]);
                Hx[i][q][k]=c3*Hx[i][q][k]+c4*(c5*bac_Bx_top[i][j][k]-c6*temp);

            }
        }
    }

	// update pomiedzy niebem i pieklem ( gora i dolem)

	c1=c1z_H[boundz];
 	c2=c2z_H[boundz];

	c5=c5z_H[boundz];
	c6=c6z_H[boundz];

	// update centrum
	for(k=0;k<boundz;k++)
	{
		c3=c3z_H[k];
		c4=c4z_H[k];

		for(j=boundy;j<yk;j++)
    	{
        	for(i=boundx;i<xk-1;i++)
        	{
                Hx[i][j][k]=c3*Hx[i][j][k]+c4*c5*c2*(-Ez[i+1][j+1][k]+Ez[i+1][j][k]);
            }
        }
    }


	// update lewicy i prawicy
	//c1=c1_H[bound];
 	//c2=c2_H[bound];

	for(i=0;i<boundx;i++)
	{
		c5=c5x_H[i];
		c6=c6x_H[i];

		q=xs-i-2;

		for(k=0;k<boundz;k++)
		{
			c3=c3z_H[k];
			c4=c4z_H[k];


    		for(j=boundy;j<yk;j++)
    		{
     			r=j-boundy;

				// update left
                temp=bac_Bx_lef[i][r][k];
                bac_Bx_lef[i][r][k]=c1*bac_Bx_lef[i][r][k]+c2*(-Ez[i+1][j+1][k]+Ez[i+1][j][k]);
                Hx[i][j][k]=c3*Hx[i][j][k]+c4*(c5*bac_Bx_lef[i][r][k]-c6*temp);

				// update rig
				temp=bac_Bx_rig[i][r][k];
                bac_Bx_rig[i][r][k]=c1*bac_Bx_rig[i][r][k]+c2*(-Ez[q+1][j+1][k]+Ez[q+1][j][k]);
                Hx[q][j][k]=c3*Hx[q][j][k]+c4*(c5*bac_Bx_rig[i][r][k]-c6*temp);

            }
        }
    }
}

void update_bac_Hy_null_Ex()
{

    typ_pola temp;
    
    ////////
    // Hy //
    ////////

    for(i=0;i<xs;i++)
    {
        c3=cstx(c3x_H,i,xs);
        c4=cstx(c4x_H,i,xs);

        for(j=0;j<ys-1;j++)
        {
            c5=csty(c5y_H,j,ys-1);
            c6=csty(c6y_H,j,ys-1);

            for(k=0;k<boundz;k++)
            {
                c1=c1z_H[k];
                c2=c2z_H[k];

                temp=bac_By[i][j][k];
                bac_By[i][j][k]=c1*bac_By[i][j][k]+c2*(Ez[i+1][j+1][k]-Ez[i][j+1][k]);
                Hy[i][j][k]=c3*Hy[i][j][k]+c4*(c5*bac_By[i][j][k]-c6*temp);

            }

        }
    }
}
void update_bac_Hy_null_Ez()
{

    typ_pola temp;
    
    ////////
    // Hy //
    ////////

    for(i=0;i<xs;i++)
    {
        c3=cstx(c3x_H,i,xs);
        c4=cstx(c4x_H,i,xs);

        for(j=0;j<ys-1;j++)
        {
            c5=csty(c5y_H,j,ys-1);
            c6=csty(c6y_H,j,ys-1);

            for(k=0;k<boundz;k++)
            {
                c1=c1z_H[k];
                c2=c2z_H[k];

                temp=bac_By[i][j][k];
                bac_By[i][j][k]=c1*bac_By[i][j][k]+c2*(-Ex[i][j+1][k+1]+Ex[i][j+1][k]);
                Hy[i][j][k]=c3*Hy[i][j][k]+c4*(c5*bac_By[i][j][k]-c6*temp);

            }

        }
    }
}


void update_bac_Hz_null_Ey()
{

    typ_pola temp;

    ////////
    // Hz //
    ////////

    for(i=0;i<xs;i++)
    {
        c1=cstx(c1x_H,i,xs);
        c2=cstx(c2x_H,i,xs);

        for(j=0;j<ys;j++)
        {
            c3=csty(c3y_H,j,ys);
            c4=csty(c4y_H,j,ys);

            for(k=0;k<boundz;k++)
            {
                c5=c5z_H[k];
                c6=c6z_H[k];

                temp=bac_Bz[i][j][k];
                bac_Bz[i][j][k]=c1*bac_Bz[i][j][k]+c2*(Ex[i][j+1][k+1]-Ex[i][j][k+1]);
                Hz[i][j][k]=c3*Hz[i][j][k]+c4*(c5*bac_Bz[i][j][k]-c6*temp);

            }
        }
    }

}
void update_bac_Hz_null_Ex()
{

    typ_pola temp;

    ////////
    // Hz //
    ////////

    for(i=0;i<xs;i++)
    {
        c1=cstx(c1x_H,i,xs);
        c2=cstx(c2x_H,i,xs);

        for(j=0;j<ys;j++)
        {
            c3=csty(c3y_H,j,ys);
            c4=csty(c4y_H,j,ys);

            for(k=0;k<boundz;k++)
            {
                c5=c5z_H[k];
                c6=c6z_H[k];

                temp=bac_Bz[i][j][k];
                bac_Bz[i][j][k]=c1*bac_Bz[i][j][k]+c2*(-Ey[i+1][j][k+1]+Ey[i][j][k+1]);
                Hz[i][j][k]=c3*Hz[i][j][k]+c4*(c5*bac_Bz[i][j][k]-c6*temp);

            }
        }
    }

}


void update_fro_Hx_null_Ez()
{

    typ_pola temp;

    ////////
    // Hx //
    ////////

	// update na gorze i dole
	for(i=0;i<xs-1;i++)
	{
		c5=cstx(c5x_H,i,xs-1);
		c6=cstx(c6x_H,i,xs-1);

    	for(j=0;j<boundy;j++)
    	{
        	c1=c1y_H[j];
        	c2=c2y_H[j];

			q=ys-j-1;

            for(k=0;k<boundz;k++)
            {
                c3=c3z_H[k];
                c4=c4z_H[k];

                p=zs-k-1;

				// update na dole
				temp=fro_Bx_bot[i][j][k];
                fro_Bx_bot[i][j][k]=c1*fro_Bx_bot[i][j][k]+c2*(Ey[i+1][j][p+1]-Ey[i+1][j][p]);
                Hx[i][j][p]=c3*Hx[i][j][p]+c4*(c5*fro_Bx_bot[i][j][k]-c6*temp);

				// update na gorze
                temp=fro_Bx_top[i][j][k];
                fro_Bx_top[i][j][k]=c1*fro_Bx_top[i][j][k]+c2*(Ey[i+1][q][p+1]-Ey[i+1][q][p]);
                Hx[i][q][p]=c3*Hx[i][q][p]+c4*(c5*fro_Bx_top[i][j][k]-c6*temp);
            }
        }
    }

	// update pomiedzy niebem i pieklem ( gora i dolem)

	c1=c1z_H[boundz];
 	c2=c2z_H[boundz];

	c5=c5z_H[boundz];
	c6=c6z_H[boundz];

	// update centrum
	for(k=0;k<boundz;k++)
	{
		c3=c3z_H[k];
		c4=c4z_H[k];

		p=zs-k-1;

		for(j=boundy;j<yk;j++)
    	{
        	for(i=boundx;i<xk-1;i++)
        	{
                Hx[i][j][p]=c3*Hx[i][j][p]+c4*c5*c2*(Ey[i+1][j][p+1]-Ey[i+1][j][p]);
            }
        }
    }


	// update lewicy i prawicy
	//c1=c1_H[bound];
 	//c2=c2_H[bound];

	for(i=0;i<boundx;i++)
	{
		c5=c5x_H[i];
		c6=c6x_H[i];

		q=xs-i-2;

		for(k=0;k<boundz;k++)
		{
			c3=c3z_H[k];
			c4=c4z_H[k];

			p=zs-k-1;

    		for(j=boundy;j<yk;j++)
    		{
     			r=j-boundy;

				// update left
                temp=fro_Bx_lef[i][r][k];
                fro_Bx_lef[i][r][k]=c1*fro_Bx_lef[i][r][k]+c2*(Ey[i+1][j][p+1]-Ey[i+1][j][p]);
                Hx[i][j][p]=c3*Hx[i][j][p]+c4*(c5*fro_Bx_lef[i][r][k]-c6*temp);

				// update rig
                temp=fro_Bx_rig[i][r][k];
                fro_Bx_rig[i][r][k]=c1*fro_Bx_rig[i][r][k]+c2*(Ey[q+1][j][p+1]-Ey[q+1][j][p]);
                Hx[q][j][p]=c3*Hx[q][j][p]+c4*(c5*fro_Bx_rig[i][r][k]-c6*temp);

            }
        }
    }
}
void update_fro_Hx_null_Ey()
{

    typ_pola temp;

    ////////
    // Hx //
    ////////

	// update na gorze i dole
	for(i=0;i<xs-1;i++)
	{
		c5=cstx(c5x_H,i,xs-1);
		c6=cstx(c6x_H,i,xs-1);

    	for(j=0;j<boundy;j++)
    	{
        	c1=c1y_H[j];
        	c2=c2y_H[j];

			q=ys-j-1;

            for(k=0;k<boundz;k++)
            {
                c3=c3z_H[k];
                c4=c4z_H[k];

                p=zs-k-1;

				// update na dole
				temp=fro_Bx_bot[i][j][k];
                fro_Bx_bot[i][j][k]=c1*fro_Bx_bot[i][j][k]+c2*(-Ez[i+1][j+1][p]+Ez[i+1][j][p]);
                Hx[i][j][p]=c3*Hx[i][j][p]+c4*(c5*fro_Bx_bot[i][j][k]-c6*temp);

				// update na gorze
                temp=fro_Bx_top[i][j][k];
                fro_Bx_top[i][j][k]=c1*fro_Bx_top[i][j][k]+c2*(-Ez[i+1][q+1][p]+Ez[i+1][q][p]);
                Hx[i][q][p]=c3*Hx[i][q][p]+c4*(c5*fro_Bx_top[i][j][k]-c6*temp);
            }
        }
    }

	// update pomiedzy niebem i pieklem ( gora i dolem)

	c1=c1z_H[boundz];
 	c2=c2z_H[boundz];

	c5=c5z_H[boundz];
	c6=c6z_H[boundz];

	// update centrum
	for(k=0;k<boundz;k++)
	{
		c3=c3z_H[k];
		c4=c4z_H[k];

		p=zs-k-1;

		for(j=boundy;j<yk;j++)
    	{
        	for(i=boundx;i<xk-1;i++)
        	{
                Hx[i][j][p]=c3*Hx[i][j][p]+c4*c5*c2*(-Ez[i+1][j+1][p]+Ez[i+1][j][p]);
            }
        }
    }


	// update lewicy i prawicy
	//c1=c1_H[bound];
 	//c2=c2_H[bound];

	for(i=0;i<boundx;i++)
	{
		c5=c5x_H[i];
		c6=c6x_H[i];

		q=xs-i-2;

		for(k=0;k<boundz;k++)
		{
			c3=c3z_H[k];
			c4=c4z_H[k];

			p=zs-k-1;

    		for(j=boundy;j<yk;j++)
    		{
     			r=j-boundy;

				// update left
                temp=fro_Bx_lef[i][r][k];
                fro_Bx_lef[i][r][k]=c1*fro_Bx_lef[i][r][k]+c2*(-Ez[i+1][j+1][p]+Ez[i+1][j][p]);
                Hx[i][j][p]=c3*Hx[i][j][p]+c4*(c5*fro_Bx_lef[i][r][k]-c6*temp);

				// update rig
                temp=fro_Bx_rig[i][r][k];
                fro_Bx_rig[i][r][k]=c1*fro_Bx_rig[i][r][k]+c2*(-Ez[q+1][j+1][p]+Ez[q+1][j][p]);
                Hx[q][j][p]=c3*Hx[q][j][p]+c4*(c5*fro_Bx_rig[i][r][k]-c6*temp);

            }
        }
    }
}

void update_fro_Hy_null_Ex()
{

    typ_pola temp;
    ////////
    // Hy //
    ////////

    for(i=0;i<xs;i++)
    {
        c3=cstx(c3x_H,i,xs);
        c4=cstx(c4x_H,i,xs);

        for(j=0;j<ys-1;j++)
        {
            c5=csty(c5y_H,j,ys-1);
            c6=csty(c6y_H,j,ys-1);

            for(k=0;k<boundz;k++)
            {
                c1=c1z_H[k];
                c2=c2z_H[k];

				p=zs-k-1;

                temp=fro_By[i][j][k];
                fro_By[i][j][k]=c1*fro_By[i][j][k]+c2*(Ez[i+1][j+1][p]-Ez[i][j+1][p]);
                Hy[i][j][p]=c3*Hy[i][j][p]+c4*(c5*fro_By[i][j][k]-c6*temp);
            }

        }
    }
}
void update_fro_Hy_null_Ez()
{

    typ_pola temp;
    ////////
    // Hy //
    ////////

    for(i=0;i<xs;i++)
    {
        c3=cstx(c3x_H,i,xs);
        c4=cstx(c4x_H,i,xs);

        for(j=0;j<ys-1;j++)
        {
            c5=csty(c5y_H,j,ys-1);
            c6=csty(c6y_H,j,ys-1);

            for(k=0;k<boundz;k++)
            {
                c1=c1z_H[k];
                c2=c2z_H[k];

				p=zs-k-1;

                temp=fro_By[i][j][k];
                fro_By[i][j][k]=c1*fro_By[i][j][k]+c2*(-Ex[i][j+1][p+1]+Ex[i][j+1][p]);
                Hy[i][j][p]=c3*Hy[i][j][p]+c4*(c5*fro_By[i][j][k]-c6*temp);
            }

        }
    }
}

void update_fro_Hz_null_Ey()
{

    typ_pola temp;

    ////////
    // Hz //
    ////////

    for(i=0;i<xs;i++)
    {
        c1=cstx(c1x_H,i,xs);
        c2=cstx(c2x_H,i,xs);

        for(j=0;j<ys;j++)
        {
            c3=csty(c3y_H,j,ys);
            c4=csty(c4y_H,j,ys);

            for(k=0;k<boundz;k++)
            {
                c5=c5z_H[k];
                c6=c6z_H[k];

				p=zs-k-2;

                temp=fro_Bz[i][j][k];
                fro_Bz[i][j][k]=c1*fro_Bz[i][j][k]+c2*(Ex[i][j+1][p+1]-Ex[i][j][p+1]);
                Hz[i][j][p]=c3*Hz[i][j][p]+c4*(c5*fro_Bz[i][j][k]-c6*temp);
            }
        }
    }

}
void update_fro_Hz_null_Ex()
{

    typ_pola temp;

    ////////
    // Hz //
    ////////

    for(i=0;i<xs;i++)
    {
        c1=cstx(c1x_H,i,xs);
        c2=cstx(c2x_H,i,xs);

        for(j=0;j<ys;j++)
        {
            c3=csty(c3y_H,j,ys);
            c4=csty(c4y_H,j,ys);

            for(k=0;k<boundz;k++)
            {
                c5=c5z_H[k];
                c6=c6z_H[k];

				p=zs-k-2;

                temp=fro_Bz[i][j][k];
                fro_Bz[i][j][k]=c1*fro_Bz[i][j][k]+c2*(-Ey[i+1][j][p+1]+Ey[i][j][p+1]);
                Hz[i][j][p]=c3*Hz[i][j][p]+c4*(c5*fro_Bz[i][j][k]-c6*temp);
            }
        }
    }

}

void update_bac_Ex_null_Hz()
{

    typ_pola temp;

    ////////
    // Ex //
    ////////

	for(i=0;i<xs;i++)
	{
		c5=cstx(c5x_E,i,xs);
		c6=cstx(c6x_E,i,xs);

    	for(j=1;j<boundy;j++)
    	{
        	c1=c1y_E[j];
        	c2=c2y_E[j];

			q=ys-j;

            for(k=1;k<boundz;k++)
            {
                c3=c3z_E[k];
                c4=c4z_E[k];

				p=zs-k;

				temp=bac_Dx_bot[i][j][k];
                bac_Dx_bot[i][j][k]=c1*bac_Dx_bot[i][j][k]+c2*(-Hy[i][j-1][k]+Hy[i][j-1][k-1]);
                Ex[i][j][k]=c3*Ex[i][j][k]+c4*(c5*bac_Dx_bot[i][j][k]-c6*temp);

                temp=bac_Dx_top[i][j][k];
                bac_Dx_top[i][j][k]=c1*bac_Dx_top[i][j][k]+c2*(-Hy[i][q-1][k]+Hy[i][q-1][k-1]);
                Ex[i][q][k]=c3*Ex[i][q][k]+c4*(c5*bac_Dx_top[i][j][k]-c6*temp);

			}
        }
    }

	c1=c1z_E[boundz];
 	c2=c2z_E[boundz];

	c5=c5z_E[boundz];
	c6=c6z_E[boundz];

	for(j=boundy;j<=yk;j++)
	{
        for(i=boundx;i<xk;i++)
        {
            for(k=1;k<boundz;k++)
            {
                c3=c3z_E[k];
                c4=c4z_E[k];

                Ex[i][j][k]=c3*Ex[i][j][k]+c4*c5*c2*(-Hy[i][j-1][k]+Hy[i][j-1][k-1]);
            }
        }
    }

	//c1=c1_E[bound];
	//c2=c2_E[bound];

	for(i=0;i<boundx;i++)
	{
		c5=c5x_E[i];
		c6=c6x_E[i];

		q=xs-i-1;

  		for(k=1;k<boundz;k++)
    	{
			c3=c3z_E[k];
			c4=c4z_E[k];

			p=zs-k;

			for(j=boundy;j<=yk;j++)
    		{
    			r=j-boundy;

				temp=bac_Dx_lef[i][r][k];
                bac_Dx_lef[i][r][k]=c1*bac_Dx_lef[i][r][k]+c2*(-Hy[i][j-1][k]+Hy[i][j-1][k-1]);
                Ex[i][j][k]=c3*Ex[i][j][k]+c4*(c5*bac_Dx_lef[i][r][k]-c6*temp);

				temp=bac_Dx_rig[i][r][k];
                bac_Dx_rig[i][r][k]=c1*bac_Dx_rig[i][r][k]+c2*(-Hy[q][j-1][k]+Hy[q][j-1][k-1]);
                Ex[q][j][k]=c3*Ex[q][j][k]+c4*(c5*bac_Dx_rig[i][r][k]-c6*temp);

            }
        }
    }
}
void update_bac_Ex_null_Hy()
{

    typ_pola temp;

    ////////
    // Ex //
    ////////

	for(i=0;i<xs;i++)
	{
		c5=cstx(c5x_E,i,xs);
		c6=cstx(c6x_E,i,xs);

    	for(j=1;j<boundy;j++)
    	{
        	c1=c1y_E[j];
        	c2=c2y_E[j];

			q=ys-j;

            for(k=1;k<boundz;k++)
            {
                c3=c3z_E[k];
                c4=c4z_E[k];

				p=zs-k;

				temp=bac_Dx_bot[i][j][k];
                bac_Dx_bot[i][j][k]=c1*bac_Dx_bot[i][j][k]+c2*(Hz[i][j][k-1]-Hz[i][j-1][k-1]);
                Ex[i][j][k]=c3*Ex[i][j][k]+c4*(c5*bac_Dx_bot[i][j][k]-c6*temp);

                temp=bac_Dx_top[i][j][k];
                bac_Dx_top[i][j][k]=c1*bac_Dx_top[i][j][k]+c2*(Hz[i][q][k-1]-Hz[i][q-1][k-1]);
                Ex[i][q][k]=c3*Ex[i][q][k]+c4*(c5*bac_Dx_top[i][j][k]-c6*temp);

			}
        }
    }

	c1=c1z_E[boundz];
 	c2=c2z_E[boundz];

	c5=c5z_E[boundz];
	c6=c6z_E[boundz];

	for(j=boundy;j<=yk;j++)
	{
        for(i=boundx;i<xk;i++)
        {
            for(k=1;k<boundz;k++)
            {
                c3=c3z_E[k];
                c4=c4z_E[k];

                Ex[i][j][k]=c3*Ex[i][j][k]+c4*c5*c2*(Hz[i][j][k-1]-Hz[i][j-1][k-1]);
            }
        }
    }

	//c1=c1_E[bound];
	//c2=c2_E[bound];

	for(i=0;i<boundx;i++)
	{
		c5=c5x_E[i];
		c6=c6x_E[i];

		q=xs-i-1;

  		for(k=1;k<boundz;k++)
    	{
			c3=c3z_E[k];
			c4=c4z_E[k];

			p=zs-k;

			for(j=boundy;j<=yk;j++)
    		{
    			r=j-boundy;

				temp=bac_Dx_lef[i][r][k];
                bac_Dx_lef[i][r][k]=c1*bac_Dx_lef[i][r][k]+c2*(Hz[i][j][k-1]-Hz[i][j-1][k-1]);
                Ex[i][j][k]=c3*Ex[i][j][k]+c4*(c5*bac_Dx_lef[i][r][k]-c6*temp);

				temp=bac_Dx_rig[i][r][k];
                bac_Dx_rig[i][r][k]=c1*bac_Dx_rig[i][r][k]+c2*(Hz[q][j][k-1]-Hz[q][j-1][k-1]);
                Ex[q][j][k]=c3*Ex[q][j][k]+c4*(c5*bac_Dx_rig[i][r][k]-c6*temp);

            }
        }
    }
}

void update_bac_Ey_null_Hz()
{

    typ_pola temp;

    ////////
    // Ey //
    ////////

    for(i=1;i<xs;i++)
    {
        c3=cstx(c3x_E,i,xs+1);
        c4=cstx(c4x_E,i,xs+1);

        for(j=0;j<ys;j++)
        {
            c5=csty(c5y_E,j,ys);
            c6=csty(c6y_E,j,ys);

            for(k=1;k<boundz;k++)
            {
                c1=c1z_E[k];
                c2=c2z_E[k];

		p=zs-k;

                temp=bac_Dy[i][j][k];
                bac_Dy[i][j][k]=c1*bac_Dy[i][j][k]+c2*(Hx[i-1][j][k]-Hx[i-1][j][k-1]);
                Ey[i][j][k]=c3*Ey[i][j][k]+c4*(c5*bac_Dy[i][j][k]-c6*temp);

            }
        }
    }
}
void update_bac_Ey_null_Hx()
{

    typ_pola temp;

    ////////
    // Ey //
    ////////

    for(i=1;i<xs;i++)
    {
        c3=cstx(c3x_E,i,xs+1);
        c4=cstx(c4x_E,i,xs+1);

        for(j=0;j<ys;j++)
        {
            c5=csty(c5y_E,j,ys);
            c6=csty(c6y_E,j,ys);

            for(k=1;k<boundz;k++)
            {
                c1=c1z_E[k];
                c2=c2z_E[k];

		        p=zs-k;

                temp=bac_Dy[i][j][k];
                bac_Dy[i][j][k]=c1*bac_Dy[i][j][k]+c2*(-Hz[i][j][k-1]+Hz[i-1][j][k-1]);
                Ey[i][j][k]=c3*Ey[i][j][k]+c4*(c5*bac_Dy[i][j][k]-c6*temp);

            }
        }
    }
}
void update_bac_Ez_null_Hx()
{

    typ_pola temp;

    ////////
    // Ez //
    ////////

    for(i=1;i<xs;i++)
    {
        c1=cstx(c1x_E,i,xs+1);
        c2=cstx(c2x_E,i,xs+1);

        for(j=1;j<ys;j++)
        {
            c3=csty(c3y_E,j,ys+1);
            c4=csty(c4y_E,j,ys+1);

            for(k=0;k<boundz;k++)
            {
                c5=c5z_E[k];
                c6=c6z_E[k];

				p=zs-k-1;

                temp=bac_Dz[i][j][k];
                bac_Dz[i][j][k]=c1*bac_Dz[i][j][k]+c2*(Hy[i][j-1][k]-Hy[i-1][j-1][k]);
                Ez[i][j][k]=c3*Ez[i][j][k]+c4*(c5*bac_Dz[i][j][k]-c6*temp);

            }
        }
    }

}
void update_bac_Ez_null_Hy()
{

    typ_pola temp;

    ////////
    // Ez //
    ////////

    for(i=1;i<xs;i++)
    {
        c1=cstx(c1x_E,i,xs+1);
        c2=cstx(c2x_E,i,xs+1);

        for(j=1;j<ys;j++)
        {
            c3=csty(c3y_E,j,ys+1);
            c4=csty(c4y_E,j,ys+1);

            for(k=0;k<boundz;k++)
            {
                c5=c5z_E[k];
                c6=c6z_E[k];

				p=zs-k-1;

                temp=bac_Dz[i][j][k];
                bac_Dz[i][j][k]=c1*bac_Dz[i][j][k]+c2*(-Hx[i-1][j][k]+Hx[i-1][j-1][k]);
                Ez[i][j][k]=c3*Ez[i][j][k]+c4*(c5*bac_Dz[i][j][k]-c6*temp);

            }
        }
    }

}

void update_fro_Ex_null_Hy()
{

    typ_pola temp;

    ////////
    // Ex //
    ////////

	for(i=0;i<xs;i++)
	{
		c5=cstx(c5x_E,i,xs);
		c6=cstx(c6x_E,i,xs);

    	for(j=1;j<boundy;j++)
    	{
        	c1=c1y_E[j];
        	c2=c2y_E[j];

			q=ys-j;

            for(k=1;k<boundz;k++)
            {
                c3=c3z_E[k];
                c4=c4z_E[k];

				p=zs-k;

                temp=fro_Dx_bot[i][j][k];
                fro_Dx_bot[i][j][k]=c1*fro_Dx_bot[i][j][k]+c2*(Hz[i][j][p-1]-Hz[i][j-1][p-1]);
                Ex[i][j][p]=c3*Ex[i][j][p]+c4*(c5*fro_Dx_bot[i][j][k]-c6*temp);

                temp=fro_Dx_top[i][j][k];
                fro_Dx_top[i][j][k]=c1*fro_Dx_top[i][j][k]+c2*(Hz[i][q][p-1]-Hz[i][q-1][p-1]);
                Ex[i][q][p]=c3*Ex[i][q][p]+c4*(c5*fro_Dx_top[i][j][k]-c6*temp);

			}
        }
    }

	c1=c1z_E[boundz];
 	c2=c2z_E[boundz];

	c5=c5z_E[boundz];
	c6=c6z_E[boundz];

	for(j=boundy;j<=yk;j++)
	{
        for(i=boundx;i<xk;i++)
        {
            for(k=1;k<boundz;k++)
            {
                c3=c3z_E[k];
                c4=c4z_E[k];

                p=zs-k;

                Ex[i][j][p]=c3*Ex[i][j][p]+c4*c5*c2*(Hz[i][j][p-1]-Hz[i][j-1][p-1]);
            }
        }
    }

	//c1=c1_E[bound];
 	//c2=c2_E[bound];

	for(i=0;i<boundx;i++)
	{
		c5=c5x_E[i];
		c6=c6x_E[i];

		q=xs-i-1;

  		for(k=1;k<boundz;k++)
    	{
			c3=c3z_E[k];
			c4=c4z_E[k];

			p=zs-k;

			for(j=boundy;j<=yk;j++)
    		{
    			r=j-boundy;

                temp=fro_Dx_lef[i][r][k];
                fro_Dx_lef[i][r][k]=c1*fro_Dx_lef[i][r][k]+c2*(Hz[i][j][p-1]-Hz[i][j-1][p-1]);
                Ex[i][j][p]=c3*Ex[i][j][p]+c4*(c5*fro_Dx_lef[i][r][k]-c6*temp);

				temp=fro_Dx_rig[i][r][k];
                fro_Dx_rig[i][r][k]=c1*fro_Dx_rig[i][r][k]+c2*(Hz[q][j][p-1]-Hz[q][j-1][p-1]);
                Ex[q][j][p]=c3*Ex[q][j][p]+c4*(c5*fro_Dx_rig[i][r][k]-c6*temp);
            }
        }
    }
}
void update_fro_Ex_null_Hz()
{

    typ_pola temp;

    ////////
    // Ex //
    ////////

	for(i=0;i<xs;i++)
	{
		c5=cstx(c5x_E,i,xs);
		c6=cstx(c6x_E,i,xs);

    	for(j=1;j<boundy;j++)
    	{
        	c1=c1y_E[j];
        	c2=c2y_E[j];

			q=ys-j;

            for(k=1;k<boundz;k++)
            {
                c3=c3z_E[k];
                c4=c4z_E[k];

				p=zs-k;

                temp=fro_Dx_bot[i][j][k];
                fro_Dx_bot[i][j][k]=c1*fro_Dx_bot[i][j][k]+c2*(-Hy[i][j-1][p]+Hy[i][j-1][p-1]);
                Ex[i][j][p]=c3*Ex[i][j][p]+c4*(c5*fro_Dx_bot[i][j][k]-c6*temp);

                temp=fro_Dx_top[i][j][k];
                fro_Dx_top[i][j][k]=c1*fro_Dx_top[i][j][k]+c2*(-Hy[i][q-1][p]+Hy[i][q-1][p-1]);
                Ex[i][q][p]=c3*Ex[i][q][p]+c4*(c5*fro_Dx_top[i][j][k]-c6*temp);

			}
        }
    }

	c1=c1z_E[boundz];
 	c2=c2z_E[boundz];

	c5=c5z_E[boundz];
	c6=c6z_E[boundz];

	for(j=boundy;j<=yk;j++)
	{
        for(i=boundx;i<xk;i++)
        {
            for(k=1;k<boundz;k++)
            {
                c3=c3z_E[k];
                c4=c4z_E[k];

                p=zs-k;

                Ex[i][j][p]=c3*Ex[i][j][p]+c4*c5*c2*(-Hy[i][j-1][p]+Hy[i][j-1][p-1]);
            }
        }
    }

	//c1=c1_E[bound];
 	//c2=c2_E[bound];

	for(i=0;i<boundx;i++)
	{
		c5=c5x_E[i];
		c6=c6x_E[i];

		q=xs-i-1;

  		for(k=1;k<boundz;k++)
    	{
			c3=c3z_E[k];
			c4=c4z_E[k];

			p=zs-k;

			for(j=boundy;j<=yk;j++)
    		{
    			r=j-boundy;

                temp=fro_Dx_lef[i][r][k];
                fro_Dx_lef[i][r][k]=c1*fro_Dx_lef[i][r][k]+c2*(-Hy[i][j-1][p]+Hy[i][j-1][p-1]);
                Ex[i][j][p]=c3*Ex[i][j][p]+c4*(c5*fro_Dx_lef[i][r][k]-c6*temp);

				temp=fro_Dx_rig[i][r][k];
                fro_Dx_rig[i][r][k]=c1*fro_Dx_rig[i][r][k]+c2*(-Hy[q][j-1][p]+Hy[q][j-1][p-1]);
                Ex[q][j][p]=c3*Ex[q][j][p]+c4*(c5*fro_Dx_rig[i][r][k]-c6*temp);
            }
        }
    }
}

void update_fro_Ey_null_Hz()
{

    typ_pola temp;

    ////////
    // Ey //
    ////////

    for(i=1;i<xs;i++)
    {
        c3=cstx(c3x_E,i,xs+1);
        c4=cstx(c4x_E,i,xs+1);

        for(j=0;j<ys;j++)
        {
            c5=csty(c5y_E,j,ys);
            c6=csty(c6y_E,j,ys);

            for(k=1;k<boundz;k++)
            {
                c1=c1z_E[k];
                c2=c2z_E[k];

				p=zs-k;

                temp=fro_Dy[i][j][k];
                fro_Dy[i][j][k]=c1*fro_Dy[i][j][k]+c2*(Hx[i-1][j][p]-Hx[i-1][j][p-1]);
                Ey[i][j][p]=c3*Ey[i][j][p]+c4*(c5*fro_Dy[i][j][k]-c6*temp);
            }
        }
    }
}

void update_fro_Ey_null_Hx()
{

    typ_pola temp;

    ////////
    // Ey //
    ////////

    for(i=1;i<xs;i++)
    {
        c3=cstx(c3x_E,i,xs+1);
        c4=cstx(c4x_E,i,xs+1);

        for(j=0;j<ys;j++)
        {
            c5=csty(c5y_E,j,ys);
            c6=csty(c6y_E,j,ys);

            for(k=1;k<boundz;k++)
            {
                c1=c1z_E[k];
                c2=c2z_E[k];

				p=zs-k;

                temp=fro_Dy[i][j][k];
                fro_Dy[i][j][k]=c1*fro_Dy[i][j][k]+c2*(-Hz[i][j][p-1]+Hz[i-1][j][p-1]);
                Ey[i][j][p]=c3*Ey[i][j][p]+c4*(c5*fro_Dy[i][j][k]-c6*temp);
            }
        }
    }
}

void update_fro_Ez_null_Hx()
{

    typ_pola temp;

    ////////
    // Ez //
    ////////

    for(i=1;i<xs;i++)
    {
        c1=cstx(c1x_E,i,xs+1);
        c2=cstx(c2x_E,i,xs+1);

        for(j=1;j<ys;j++)
        {
            c3=csty(c3y_E,j,ys+1);
            c4=csty(c4y_E,j,ys+1);

            for(k=0;k<boundz;k++)
            {
                c5=c5z_E[k];
                c6=c6z_E[k];

				p=zs-k-1;

                temp=fro_Dz[i][j][k];
                fro_Dz[i][j][k]=c1*fro_Dz[i][j][k]+c2*(Hy[i][j-1][p]-Hy[i-1][j-1][p]);
                Ez[i][j][p]=c3*Ez[i][j][p]+c4*(c5*fro_Dz[i][j][k]-c6*temp);
            }
        }
    }

}

void update_fro_Ez_null_Hy()
{

    typ_pola temp;

    ////////
    // Ez //
    ////////

    for(i=1;i<xs;i++)
    {
        c1=cstx(c1x_E,i,xs+1);
        c2=cstx(c2x_E,i,xs+1);

        for(j=1;j<ys;j++)
        {
            c3=csty(c3y_E,j,ys+1);
            c4=csty(c4y_E,j,ys+1);

            for(k=0;k<boundz;k++)
            {
                c5=c5z_E[k];
                c6=c6z_E[k];

				p=zs-k-1;

                temp=fro_Dz[i][j][k];
                fro_Dz[i][j][k]=c1*fro_Dz[i][j][k]+c2*(-Hx[i-1][j][p]+Hx[i-1][j-1][p]);
                Ez[i][j][p]=c3*Ez[i][j][p]+c4*(c5*fro_Dz[i][j][k]-c6*temp);
            }
        }
    }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Osrodki dyspersyjne
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void update_drdE3x()
{
	short int tempsi;

	// indeks mat. drudego
	p=0;
    for(i=drd_Ex_is;i<=drd_Ex_ik;i++)
    {                           
        for(j=drd_Ex_js;j<=drd_Ex_jk;j++)
        {
            for(k=drd_Ex_ks;k<=drd_Ex_kk;k++)
			{
                                             
				tempsi=matEx3[i-boundx][j-boundy][k-boundz];
				if (is_drd_E[tempsi])
				{
					Ex[i][j][k]+=-drd_J_E[tempsi]*drd_Jex[p];
					p++;
				}
			}
		}
    }
}

void update_drdJ3x()
{
	short int tempsi;

	// indeks mat. drudego
	p=0;
    for(i=drd_Ex_is;i<=drd_Ex_ik;i++)
    {                           
        for(j=drd_Ex_js;j<=drd_Ex_jk;j++)
        {
            for(k=drd_Ex_ks;k<=drd_Ex_kk;k++)
			{
                                             
				tempsi=matEx3[i-boundx][j-boundy][k-boundz];
				if (is_drd_E[tempsi])
				{
					drd_Jex[p]=drd_ka_E[tempsi]*drd_Jex[p]+(Ex[i][j][k]+drd_Ex_n_1[p]);
					drd_Ex_n_1[p]=Ex[i][j][k];
					p++;
				}
			}
		}
    }
}


void update_drdE3y()
{
	short int tempsi;

	p=0;
    for(i=drd_Ey_is;i<=drd_Ey_ik;i++)
    {
        for(j=drd_Ey_js;j<=drd_Ey_jk;j++)
        {
            for(k=drd_Ey_ks;k<=drd_Ey_kk;k++)
			{
				tempsi=matEy3[i-boundx][j-boundy][k-boundz];
				if (is_drd_E[tempsi])
				{
					Ey[i][j][k]+=-drd_J_E[tempsi]*drd_Jey[p];
					p++;
				}
			}
		}
    }
}

void update_drdJ3y()
{
	short int tempsi;

	p=0;
    for(i=drd_Ey_is;i<=drd_Ey_ik;i++)
    {
        for(j=drd_Ey_js;j<=drd_Ey_jk;j++)
        {
            for(k=drd_Ey_ks;k<=drd_Ey_kk;k++)
			{
				tempsi=matEy3[i-boundx][j-boundy][k-boundz];
				if (is_drd_E[tempsi])
				{
					drd_Jey[p]=drd_ka_E[tempsi]*drd_Jey[p]+(Ey[i][j][k]+drd_Ey_n_1[p]);
					drd_Ey_n_1[p]=Ey[i][j][k];
					p++;
				}
			}
		}
    }
}


void update_drdE3z()
{
	short int tempsi;
	p=0;
    for(i=drd_Ez_is;i<=drd_Ez_ik;i++)
    {
        for(j=drd_Ez_js;j<=drd_Ez_jk;j++)
        {
            for(k=drd_Ez_ks;k<=drd_Ez_kk;k++)
			{
				tempsi=matEz3[i-boundx][j-boundy][k-boundz];
				if (is_drd_E[tempsi])
				{
					Ez[i][j][k]+=-drd_J_E[tempsi]*drd_Jez[p];
					p++;
				}
			}
		}
    }

}

void update_drdJ3z()
{
	short int tempsi;
	p=0;
    for(i=drd_Ez_is;i<=drd_Ez_ik;i++)
    {
        for(j=drd_Ez_js;j<=drd_Ez_jk;j++)
        {
            for(k=drd_Ez_ks;k<=drd_Ez_kk;k++)
			{
				tempsi=matEz3[i-boundx][j-boundy][k-boundz];
				if (is_drd_E[tempsi])
				{
					drd_Jez[p]=drd_ka_E[tempsi]*drd_Jez[p]+(Ez[i][j][k]+drd_Ez_n_1[p]);
					drd_Ez_n_1[p]=Ez[i][j][k];
					p++;
				}
			}
		}
    }

}


void update_drdEx()
{
	short int tempsi;
	typ_prec tmpv;

	// indeks mat. drudego
    p=0;
    for(i=drd_Ex_is;i<=drd_Ex_ik;i++)
    {
        for(j=drd_Ex_js;j<=drd_Ex_jk;j++)
        {
			tempsi=matEx[i-boundx][j-boundy];
			if (is_drd_E[tempsi])
			{
				tmpv=drd_J_E[tempsi];
				for(k=boundz;k<=zk;k++)
				{
					Ex[i][j][k]+=-tmpv*drd_Jex[p];
					p++;
				}
			}
			if (p>=n_drd_Ex) break;
		}
    }

}

void update_drdJx()
{
	short int tempsi;

	// indeks mat. drudego
    p=0;
    for(i=drd_Ex_is;i<=drd_Ex_ik;i++)
    {
        for(j=drd_Ex_js;j<=drd_Ex_jk;j++)
        {
			tempsi=matEx[i-boundx][j-boundy];
			if (is_drd_E[tempsi])
			{
				for(k=boundz;k<=zk;k++)
				{
					drd_Jex[p]=drd_ka_E[tempsi]*drd_Jex[p]+(Ex[i][j][k]+drd_Ex_n_1[p]);
					drd_Ex_n_1[p]=Ex[i][j][k];
					p++;
				}
			}
			if (p>=n_drd_Ex) break;
		}
    }

}


void update_drdEy()
{
	short int tempsi;
	typ_prec tmpv;

	p=0;
    for(i=drd_Ey_is;i<=drd_Ey_ik;i++)
    {
        for(j=drd_Ey_js;j<=drd_Ey_jk;j++)
        {
			tempsi=matEy[i-boundx][j-boundy];
			if (is_drd_E[tempsi])
			{
				tmpv=drd_J_E[tempsi];
				for(k=boundz;k<=zk;k++)
				{
					Ey[i][j][k]+=-tmpv*drd_Jey[p];
					p++;
				}
			}
			if (p>=n_drd_Ey) break;
		}
    }

}

void update_drdJy()
{
	short int tempsi;

	p=0;
    for(i=drd_Ey_is;i<=drd_Ey_ik;i++)
    {
        for(j=drd_Ey_js;j<=drd_Ey_jk;j++)
        {
			tempsi=matEy[i-boundx][j-boundy];
			if (is_drd_E[tempsi])
			{
				for(k=boundz;k<=zk;k++)
				{
					drd_Jey[p]=drd_ka_E[tempsi]*drd_Jey[p]+(Ey[i][j][k]+drd_Ey_n_1[p]);
					drd_Ey_n_1[p]=Ey[i][j][k];
					p++;
				}
			}
			if (p>=n_drd_Ey) break;
		}
    }

}


void update_drdEz()
{
	short int tempsi;
	typ_prec tmpv;

	p=0;
    for(i=drd_Ez_is;i<=drd_Ez_ik;i++)
    {
        for(j=drd_Ez_js;j<=drd_Ez_jk;j++)
        {
			tempsi=matEz[i-boundx][j-boundy];
			if (is_drd_E[tempsi])
			{
				tmpv=drd_J_E[tempsi];;
            	for(k=boundz;k<zk;k++)
				{
					Ez[i][j][k]+=-tmpv*drd_Jez[p];
					p++;
				}
			}
			if (p>=n_drd_Ez) break;
		}
    }

}

void update_drdJz()
{
	short int tempsi;

	p=0;
    for(i=drd_Ez_is;i<=drd_Ez_ik;i++)
    {
        for(j=drd_Ez_js;j<=drd_Ez_jk;j++)
        {
			tempsi=matEz[i-boundx][j-boundy];
			if (is_drd_E[tempsi])
			{
            	for(k=boundz;k<zk;k++)
				{
					drd_Jey[p]=drd_ka_E[tempsi]*drd_Jez[p]+Ez[i][j][k]+drd_Ez_n_1[p];
					drd_Ez_n_1[p]=Ez[i][j][k];
					p++;
				}
			}
			if (p>=n_drd_Ez) break;
		}
    }

}


void update_drdH3x()
{
	short int tempsi;
	// indeks mat. drudego
	p=0;
    for(i=drd_Hx_is;i<=drd_Hx_ik;i++)
    {
        for(j=drd_Hx_js;j<=drd_Hx_jk;j++)
        {
            for(k=drd_Hx_ks;k<=drd_Hx_kk;k++)
			{
				tempsi=matHx3[i-boundx][j-boundy][k-boundz];
				if (is_drd_H[tempsi])
				{
					Hx[i][j][k]+=-drd_J_H[tempsi]*drd_Jhx[p];
					p++;
				}
			}
		}
    }

}

void update_drdJh3x()
{
	short int tempsi;
	// indeks mat. drudego
	p=0;
    for(i=drd_Hx_is;i<=drd_Hx_ik;i++)
    {
        for(j=drd_Hx_js;j<=drd_Hx_jk;j++)
        {
            for(k=drd_Hx_ks;k<=drd_Hx_kk;k++)
			{
				tempsi=matHx3[i-boundx][j-boundy][k-boundz];
				if (is_drd_H[tempsi])
				{
					drd_Jhx[p]=drd_ka_H[tempsi]*drd_Jhx[p]+Hx[i][j][k]+drd_Hx_n_1[p];
					drd_Hx_n_1[p]=Hx[i][j][k];
					p++;
				}
			}
		}
    }

}


void update_drdH3y()
{
	short int tempsi;

	p=0;
    for(i=drd_Hy_is;i<=drd_Hy_ik;i++)
    {
        for(j=drd_Hy_js;j<=drd_Hy_jk;j++)
        {
            for(k=drd_Hy_ks;k<=drd_Hy_kk;k++)
			{
				tempsi=matHy3[i-boundx][j-boundy][k-boundz];

				if (is_drd_H[tempsi])
				{
					Hy[i][j][k]+=-drd_J_H[tempsi]*drd_Jhy[p];
					p++;
				}
			}
		}
    }

}

void update_drdJh3y()
{
	short int tempsi;

	p=0;
    for(i=drd_Hy_is;i<=drd_Hy_ik;i++)
    {
        for(j=drd_Hy_js;j<=drd_Hy_jk;j++)
        {
            for(k=drd_Hy_ks;k<=drd_Hy_kk;k++)
			{
				tempsi=matHy3[i-boundx][j-boundy][k-boundz];

				if (is_drd_H[tempsi])
				{
					drd_Jhy[p]=drd_ka_H[tempsi]*drd_Jhy[p]+Hy[i][j][k]+drd_Hy_n_1[p];
					drd_Hy_n_1[p]=Hy[i][j][k];
					p++;
				}
			}
		}
    }

}


void update_drdH3z()
{
	short int tempsi;

	p=0;
    for(i=drd_Hz_is;i<=drd_Hz_ik;i++)
    {
        for(j=drd_Hz_js;j<=drd_Hz_jk;j++)
        {
            for(k=drd_Hz_ks;k<=drd_Hz_kk;k++)
			{
				tempsi=matHz3[i-boundx][j-boundy][k-boundz];
				if (is_drd_H[tempsi])
				{
					Hz[i][j][k]+=-drd_J_H[tempsi]*drd_Jhz[p];
					p++;
				}
			}
		}
    }
}

void update_drdJh3z()
{
	short int tempsi;

	p=0;
    for(i=drd_Hz_is;i<=drd_Hz_ik;i++)
    {
        for(j=drd_Hz_js;j<=drd_Hz_jk;j++)
        {
            for(k=drd_Hz_ks;k<=drd_Hz_kk;k++)
			{
				tempsi=matHz3[i-boundx][j-boundy][k-boundz];
				if (is_drd_H[tempsi])
				{
					drd_Jhz[p]=drd_ka_H[tempsi]*drd_Jhz[p]+(Hz[i][j][k]+drd_Hz_n_1[p]);
					drd_Hz_n_1[p]=Hz[i][j][k];
					p++;
				}
			}
		}
    }
}



void update_drdHx()
{
	short int tempsi;
	typ_prec tmpv;
	// indeks mat. drudego
    p=0;

    for(i=drd_Hx_is;i<=drd_Hx_ik;i++)
    {
        for(j=drd_Hx_js;j<=drd_Hx_jk;j++)
        {
			tempsi=matHx[i-boundx][j-boundy];
			if (is_drd_H[tempsi])
			{
				tmpv=drd_J_H[tempsi];
				for(k=boundz;k<zk;k++)
				{
					Hx[i][j][k]+=-tmpv*drd_Jhx[p];
					p++;
				}
			}
			if (p>=n_drd_Hx) break;
		}
    }
}

void update_drdJhx()
{
	short int tempsi;
	// indeks mat. drudego
    p=0;

    for(i=drd_Hx_is;i<=drd_Hx_ik;i++)
    {
        for(j=drd_Hx_js;j<=drd_Hx_jk;j++)
        {
			tempsi=matHx[i-boundx][j-boundy];
			if (is_drd_H[tempsi])
			{
				for(k=boundz;k<zk;k++)
				{
					drd_Jhx[p]=drd_ka_H[tempsi]*drd_Jhx[p]+(Hx[i][j][k]+drd_Hx_n_1[p]);
					drd_Hx_n_1[p]=Hx[i][j][k];
					p++;
				}
			}
			if (p>=n_drd_Hx) break;
		}
    }
}



void update_drdHy()
{
	short int tempsi;
	typ_prec tmpv;

	p=0;
    for(i=drd_Hy_is;i<=drd_Hy_ik;i++)
    {
        for(j=drd_Hy_js;j<=drd_Hy_jk;j++)
        {
			tempsi=matHy[i-boundx][j-boundy];
			if (is_drd_H[tempsi])
			{
				tmpv=drd_J_H[tempsi];
				for(k=boundz;k<zk;k++)
				{
					Hy[i][j][k]+=-tmpv*drd_Jhy[p];
					p++;
				}
			}
			if (p>=n_drd_Hy) break;
		}
    }
}

void update_drdJhy()
{
	short int tempsi;

	p=0;
    for(i=drd_Hy_is;i<=drd_Hy_ik;i++)
    {
        for(j=drd_Hy_js;j<=drd_Hy_jk;j++)
        {
			tempsi=matHy[i-boundx][j-boundy];
			if (is_drd_H[tempsi])
			{
				for(k=boundz;k<zk;k++)
				{
					drd_Jhy[p]=drd_ka_H[tempsi]*drd_Jhy[p]+(Hy[i][j][k]+drd_Hy_n_1[p]);
					drd_Hy_n_1[p]=Hy[i][j][k];
					p++;
				}
			}
			if (p>=n_drd_Hy) break;
		}
    }
}


void update_drdHz()
{
	short int tempsi;
	typ_prec tmpv;

	p=0;
    for(i=drd_Hz_is;i<=drd_Hz_ik;i++)
    {
        for(j=drd_Hz_js;j<=drd_Hz_jk;j++)
        {
			tempsi=matHz[i-boundx][j-boundy];
			if (is_drd_H[tempsi])
			{
				tmpv=drd_J_H[tempsi];
            	for(k=boundz;k<zk-1;k++)
				{
					Hz[i][j][k]+=-tmpv*drd_Jhz[p];
					p++;
				}
			}
			if (p>=n_drd_Hz) break;
		}
    }

}

void update_drdJhz()
{
	short int tempsi;

	p=0;
    for(i=drd_Hz_is;i<=drd_Hz_ik;i++)
    {
        for(j=drd_Hz_js;j<=drd_Hz_jk;j++)
        {
			tempsi=matHz[i-boundx][j-boundy];
			if (is_drd_H[tempsi])
			{
            	for(k=boundz;k<zk-1;k++)
				{
					drd_Jhz[p]=drd_ka_H[tempsi]*drd_Jhz[p]+(Hz[i][j][k]+drd_Hz_n_1[p]);
					drd_Hz_n_1[p]=Hz[i][j][k];
					p++;
				}
			}
			if (p>=n_drd_Hz) break;
		}
    }

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void update_lorE3x()
{
	short int tempsi;

    p=0;
    for(i=lor_Ex_is;i<=lor_Ex_ik;i++)
    {
        for(j=lor_Ex_js;j<=lor_Ex_jk;j++)
        {
            for(k=lor_Ex_ks;k<=lor_Ex_kk;k++)
			{
				tempsi=matEx3[i-boundx][j-boundy][k-boundz];

				if (is_lor_E[tempsi])
				{
					Ex[i][j][k]+=lor_a1_E[tempsi]*lor_Ex_n_2[p]-lor_a2_E[tempsi]*lor_Jex[p]-lor_a3_E[tempsi]*lor_Jex_n_1[p];

					p++;
				}
			}
		}
    }
}

void update_lorJ3x()
{
	short int tempsi;
	typ_pola temp;

    p=0;
    for(i=lor_Ex_is;i<=lor_Ex_ik;i++)
    {
        for(j=lor_Ex_js;j<=lor_Ex_jk;j++)
        {
            for(k=lor_Ex_ks;k<=lor_Ex_kk;k++)
			{
				tempsi=matEx3[i-boundx][j-boundy][k-boundz];

				if (is_lor_E[tempsi])
				{
					temp=lor_Jex[p];

					lor_Jex[p]=lor_alf_E[tempsi]*lor_Jex[p]+lor_ksi_E[tempsi]*lor_Jex_n_1[p]+(Ex[i][j][k]-lor_Ex_n_2[p]);

					lor_Jex_n_1[p]=temp;     // Jx n-1

					lor_Ex_n_2[p]=lor_Ex_n_1[p]; // ex n-1
					lor_Ex_n_1[p]=Ex[i][j][k];   // ex n

					p++;
				}
			}
		}
    }
}


void update_lorE3y()
{
	short int tempsi;

	p=0;
    for(i=lor_Ey_is;i<=lor_Ey_ik;i++)
    {
        for(j=lor_Ey_js;j<=lor_Ey_jk;j++)
        {
            for(k=lor_Ey_ks;k<=lor_Ey_kk;k++)
			{
				tempsi=matEy3[i-boundx][j-boundy][k-boundz];
				if (is_lor_E[tempsi])
				{

					Ey[i][j][k]+=lor_a1_E[tempsi]*lor_Ey_n_2[p]-lor_a2_E[tempsi]*lor_Jey[p]-lor_a3_E[tempsi]*lor_Jey_n_1[p];

					p++;
				}
			}
		}
    }
}

void update_lorJ3y()
{
	short int tempsi;
	typ_pola temp;

	p=0;
    for(i=lor_Ey_is;i<=lor_Ey_ik;i++)
    {
        for(j=lor_Ey_js;j<=lor_Ey_jk;j++)
        {
            for(k=lor_Ey_ks;k<=lor_Ey_kk;k++)
			{
				tempsi=matEy3[i-boundx][j-boundy][k-boundz];
				if (is_lor_E[tempsi])
				{
					temp=lor_Jey[p];
					lor_Jey[p]=lor_alf_E[tempsi]*lor_Jey[p]+lor_ksi_E[tempsi]*lor_Jey_n_1[p]+(Ey[i][j][k]-lor_Ey_n_2[p]);

					lor_Jey_n_1[p]=temp;     // Jy n-1

					lor_Ey_n_2[p]=lor_Ey_n_1[p]; // ey n-1
					lor_Ey_n_1[p]=Ey[i][j][k];   // ey n

					p++;
				}
			}
		}
    }
}


void update_lorE3z()
{
	short int tempsi;

	p=0;
    for(i=lor_Ez_is;i<=lor_Ez_ik;i++)
    {
        for(j=lor_Ez_js;j<=lor_Ez_jk;j++)
        {
            for(k=lor_Ez_ks;k<=lor_Ez_kk;k++)
			{
				tempsi=matEz3[i-boundx][j-boundy][k-boundz];
				if (is_lor_E[tempsi])
				{

					Ez[i][j][k]+=lor_a1_E[tempsi]*lor_Ez_n_2[p]-lor_a2_E[tempsi]*lor_Jez[p]-lor_a3_E[tempsi]*lor_Jez_n_1[p];

					p++;
				}
			}
		}
    }
}

void update_lorJ3z()
{
	short int tempsi;
	typ_pola temp;

	p=0;
    for(i=lor_Ez_is;i<=lor_Ez_ik;i++)
    {
        for(j=lor_Ez_js;j<=lor_Ez_jk;j++)
        {
            for(k=lor_Ez_ks;k<=lor_Ez_kk;k++)
			{
				tempsi=matEz3[i-boundx][j-boundy][k-boundz];
				if (is_lor_E[tempsi])
				{
					temp=lor_Jez[p];
					lor_Jez[p]=lor_alf_E[tempsi]*lor_Jez[p]+lor_ksi_E[tempsi]*lor_Jez_n_1[p]+(Ez[i][j][k]-lor_Ez_n_2[p]);

					lor_Jez_n_1[p]=temp;

					lor_Ez_n_2[p]=lor_Ez_n_1[p]; // ez n-1
					lor_Ez_n_1[p]=Ez[i][j][k];   // ez n

					p++;
				}
			}
		}
    }
}


void update_lorEx()
{
	short int tempsi;
	typ_prec lor_a1_Et,lor_a2_Et,lor_a3_Et;

    p=0;
    for(i=lor_Ex_is;i<=lor_Ex_ik;i++)
    {
        for(j=lor_Ex_js;j<=lor_Ex_jk;j++)
        {
			tempsi=matEx[i-boundx][j-boundy];
			if (is_lor_E[tempsi])
			{
				lor_a1_Et=lor_a1_E[tempsi];
				lor_a2_Et=lor_a2_E[tempsi];
				lor_a3_Et=lor_a3_E[tempsi];

            	for(k=boundz;k<=zk;k++)
				{

					Ex[i][j][k]+=lor_a1_Et*lor_Ex_n_2[p]-lor_a2_Et*lor_Jex[p]-lor_a3_Et*lor_Jex_n_1[p];
					p++;
				}
			}
		}
    }

}

void update_lorJx()
{
	short int tempsi;
	typ_prec lor_alf_Et,lor_ksi_Et;
	typ_pola temp;

    p=0;
    for(i=lor_Ex_is;i<=lor_Ex_ik;i++)
    {
        for(j=lor_Ex_js;j<=lor_Ex_jk;j++)
        {
			tempsi=matEx[i-boundx][j-boundy];
			if (is_lor_E[tempsi])
			{
				lor_alf_Et=lor_alf_E[tempsi];
				lor_ksi_Et=lor_ksi_E[tempsi];

            	for(k=boundz;k<=zk;k++)
				{
					temp=lor_Jex[p];
					lor_Jex[p]=lor_alf_Et*lor_Jex[p]+lor_ksi_Et*lor_Jex_n_1[p]+(Ex[i][j][k]-lor_Ex_n_2[p]);

					lor_Jex_n_1[p]=temp;     // Jx n-1

					lor_Ex_n_2[p]=lor_Ex_n_1[p]; // ex n-1
					lor_Ex_n_1[p]=Ex[i][j][k];   // ex n

					p++;
				}
			}
		}
    }

}


void update_lorEy()
{
	short int tempsi;
	typ_prec lor_a1_Et,lor_a2_Et,lor_a3_Et;


    p=0;
    for(i=lor_Ey_is;i<=lor_Ey_ik;i++)
    {
        for(j=lor_Ey_js;j<=lor_Ey_jk;j++)
        {
			tempsi=matEy[i-boundx][j-boundy];
			if (is_lor_E[tempsi])
			{
				lor_a1_Et=lor_a1_E[tempsi];
				lor_a2_Et=lor_a2_E[tempsi];
				lor_a3_Et=lor_a3_E[tempsi];


            	for(k=boundz;k<=zk;k++)
				{
					Ey[i][j][k]+=lor_a1_Et*lor_Ey_n_2[p]-lor_a2_Et*lor_Jey[p]-lor_a3_Et*lor_Jey_n_1[p];

					p++;
				}
			}
		}
    }

}

void update_lorJy()
{
	short int tempsi;
	typ_prec lor_alf_Et,lor_ksi_Et;
	typ_pola temp;

    p=0;
    for(i=lor_Ey_is;i<=lor_Ey_ik;i++)
    {
        for(j=lor_Ey_js;j<=lor_Ey_jk;j++)
        {
			tempsi=matEy[i-boundx][j-boundy];
			if (is_lor_E[tempsi])
			{

				lor_alf_Et=lor_alf_E[tempsi];
				lor_ksi_Et=lor_ksi_E[tempsi];

            	for(k=boundz;k<=zk;k++)
				{

					temp=lor_Jey[p];
					lor_Jey[p]=lor_alf_Et*lor_Jey[p]+lor_ksi_Et*lor_Jey_n_1[p]+(Ey[i][j][k]-lor_Ey_n_2[p]);

					lor_Jey_n_1[p]=temp;     // Jy n-1

					lor_Ey_n_2[p]=lor_Ey_n_1[p]; // ey n-1
					lor_Ey_n_1[p]=Ey[i][j][k];   // ey n

					p++;
				}
			}
		}
    }

}


void update_lorEz()
{
	short int tempsi;
	typ_prec lor_a1_Et,lor_a2_Et,lor_a3_Et;


    p=0;
    for(i=lor_Ez_is;i<=lor_Ez_ik;i++)
    {
        for(j=lor_Ez_js;j<=lor_Ez_jk;j++)
        {
			tempsi=matEz[i-boundx][j-boundy];
			if (is_lor_E[tempsi])
			{
				lor_a1_Et=lor_a1_E[tempsi];
				lor_a2_Et=lor_a2_E[tempsi];
				lor_a3_Et=lor_a3_E[tempsi];

            	for(k=boundz;k<zk;k++)
				{
					Ez[i][j][k]+=lor_a1_Et*lor_Ez_n_2[p]-lor_a2_Et*lor_Jez[p]-lor_a3_Et*lor_Jez_n_1[p];

					p++;
				}
			}
		}
    }

}

void update_lorJz()
{
	short int tempsi;

	typ_prec lor_alf_Et,lor_ksi_Et;
	typ_pola temp;

    p=0;
    for(i=lor_Ez_is;i<=lor_Ez_ik;i++)
    {
        for(j=lor_Ez_js;j<=lor_Ez_jk;j++)
        {
			tempsi=matEz[i-boundx][j-boundy];
			if (is_lor_E[tempsi])
			{

				lor_alf_Et=lor_alf_E[tempsi];
				lor_ksi_Et=lor_ksi_E[tempsi];

            	for(k=boundz;k<zk;k++)
				{

					temp=lor_Jez[p];
					lor_Jez[p]=lor_alf_Et*lor_Jez[p]+lor_ksi_Et*lor_Jez_n_1[p]+(Ez[i][j][k]-lor_Ez_n_2[p]);

					lor_Jez_n_1[p]=temp;

					lor_Ez_n_2[p]=lor_Ez_n_1[p]; // ez n-1
					lor_Ez_n_1[p]=Ez[i][j][k];   // ez n

					p++;
				}
			}
		}
    }

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void update_lorH3x()
{
	short int tempsi;
	// indeks mat. debeye
    p=0;

    for(i=lor_Hx_is;i<=lor_Hx_ik;i++)
    {
        for(j=lor_Hx_js;j<=lor_Hx_jk;j++)
        {
            for(k=lor_Hx_ks;k<=lor_Hx_kk;k++)
			{
				tempsi=matHx3[i-boundx][j-boundy][k-boundz];

				if (is_lor_H[tempsi])
				{
					Hx[i][j][k]+=lor_a1_H[tempsi]*lor_Hx_n_2[p]-lor_a2_H[tempsi]*lor_Jhx[p]-lor_a3_H[tempsi]*lor_Jhx_n_1[p];
					p++;
				}
			}
		}
    }
}

void update_lorJh3x()
{
	short int tempsi;
	typ_pola temp;
	// indeks mat. debeye
    p=0;

    for(i=lor_Hx_is;i<=lor_Hx_ik;i++)
    {
        for(j=lor_Hx_js;j<=lor_Hx_jk;j++)
        {
            for(k=lor_Hx_ks;k<=lor_Hx_kk;k++)
			{
				tempsi=matHx3[i-boundx][j-boundy][k-boundz];

				if (is_lor_H[tempsi])
				{
					temp=lor_Jhx[p];
					lor_Jhx[p]=lor_alf_H[tempsi]*lor_Jhx[p]+lor_ksi_H[tempsi]*lor_Jhx_n_1[p]+(Hx[i][j][k]-lor_Hx_n_2[p]);

					lor_Jhx_n_1[p]=temp;

					lor_Hx_n_2[p]=lor_Hx_n_1[p];
					lor_Hx_n_1[p]=Hx[i][j][k];
					p++;
				}
			}
		}
    }
}

void update_lorH3y()
{
	short int tempsi;

	p=0;
    for(i=lor_Hy_is;i<=lor_Hy_ik;i++)
    {
        for(j=lor_Hy_js;j<=lor_Hy_jk;j++)
        {
            for(k=lor_Hy_ks;k<=lor_Hy_kk;k++)
			{
				tempsi=matHy3[i-boundx][j-boundy][k-boundz];
				if (is_lor_H[tempsi])
				{
					Hy[i][j][k]+=lor_a1_H[tempsi]*lor_Hy_n_2[p]-lor_a2_H[tempsi]*lor_Jhy[p]-lor_a3_H[tempsi]*lor_Jhy_n_1[p];

					p++;
				}
			}
		}
    }

}

void update_lorJh3y()
{
	short int tempsi;
	typ_pola temp;

	p=0;
    for(i=lor_Hy_is;i<=lor_Hy_ik;i++)
    {
        for(j=lor_Hy_js;j<=lor_Hy_jk;j++)
        {
            for(k=lor_Hy_ks;k<=lor_Hy_kk;k++)
			{
				tempsi=matHy3[i-boundx][j-boundy][k-boundz];
				if (is_lor_H[tempsi])
				{
					temp=lor_Jhy[p];
					lor_Jhy[p]=lor_alf_H[tempsi]*lor_Jhy[p]+lor_ksi_H[tempsi]*lor_Jhy_n_1[p]+(Hy[i][j][k]-lor_Hy_n_2[p]);

					lor_Jhy_n_1[p]=temp;

					lor_Hy_n_2[p]=lor_Hy_n_1[p];
					lor_Hy_n_1[p]=Hy[i][j][k];

					p++;
				}
			}
		}
    }

}


void update_lorH3z()
{
	short int tempsi;

	p=0;
    for(i=lor_Hz_is;i<=lor_Hz_ik;i++)
    {
        for(j=lor_Hz_js;j<=lor_Hz_jk;j++)
        {
            for(k=lor_Hz_ks;k<=lor_Hz_kk;k++)
			{
				tempsi=matHz3[i-boundx][j-boundy][k-boundz];
				if (is_lor_H[tempsi])
				{
					Hz[i][j][k]+=lor_a1_H[tempsi]*lor_Hz_n_2[p]-lor_a2_H[tempsi]*lor_Jhz[p]-lor_a3_H[tempsi]*lor_Jhz_n_1[p];

					p++;
				}
			}
		}
    }
}

void update_lorJh3z()
{
	short int tempsi;
	typ_pola temp;

	p=0;
    for(i=lor_Hz_is;i<=lor_Hz_ik;i++)
    {
        for(j=lor_Hz_js;j<=lor_Hz_jk;j++)
        {
            for(k=lor_Hz_ks;k<=lor_Hz_kk;k++)
			{
				tempsi=matHz3[i-boundx][j-boundy][k-boundz];
				if (is_lor_H[tempsi])
				{

					temp=lor_Jhz[p];
					lor_Jhz[p]=lor_alf_H[tempsi]*lor_Jhz[p]+lor_ksi_H[tempsi]*lor_Jhz_n_1[p]+(Hz[i][j][k]-lor_Hz_n_2[p]);

					lor_Jhz_n_1[p]=temp;

					lor_Hz_n_2[p]=lor_Hz_n_1[p];
					lor_Hz_n_1[p]=Hz[i][j][k];

					p++;
				}
			}
		}
    }
}



void update_lorHx()
{
	short int tempsi;
	typ_prec lor_a1_Ht,lor_a2_Ht,lor_a3_Ht;


	p=0;
    for(i=lor_Hx_is;i<=lor_Hx_ik;i++)
    {
        for(j=lor_Hx_js;j<=lor_Hx_jk;j++)
        {
			tempsi=matHx[i-boundx][j-boundy];
			if (is_lor_H[tempsi])
			{
				lor_a1_Ht=lor_a1_H[tempsi];
				lor_a2_Ht=lor_a2_H[tempsi];
				lor_a3_Ht=lor_a3_H[tempsi];


            	for(k=boundz;k<zk;k++)
				{
					Hx[i][j][k]+=lor_a1_Ht*lor_Hx_n_2[p]-lor_a2_Ht*lor_Jhx[p]-lor_a3_Ht*lor_Jhx_n_1[p];

					p++;

				}
			}
		}
    }
}

void update_lorJhx()
{
	short int tempsi;

	typ_prec lor_alf_Ht,lor_ksi_Ht;
	typ_pola temp;

	p=0;
    for(i=lor_Hx_is;i<=lor_Hx_ik;i++)
    {
        for(j=lor_Hx_js;j<=lor_Hx_jk;j++)
        {
			tempsi=matHx[i-boundx][j-boundy];
			if (is_lor_H[tempsi])
			{

				lor_alf_Ht=lor_alf_H[tempsi];
				lor_ksi_Ht=lor_ksi_H[tempsi];

            	for(k=boundz;k<zk;k++)
				{
					temp=lor_Jhx[p];
					lor_Jhx[p]=lor_alf_Ht*lor_Jhx[p]+lor_ksi_Ht*lor_Jhx_n_1[p]+(Hx[i][j][k]-lor_Hx_n_2[p]);

					lor_Jhx_n_1[p]=temp;

					lor_Hx_n_2[p]=lor_Hx_n_1[p];
					lor_Hx_n_1[p]=Hx[i][j][k];

					p++;

				}
			}
		}
    }
}


void update_lorHy()
{
	short int tempsi;
	typ_prec lor_a1_Ht,lor_a2_Ht,lor_a3_Ht;

    p=0;
    for(i=lor_Hy_is;i<=lor_Hy_ik;i++)
    {
        for(j=lor_Hy_js;j<=lor_Hy_jk;j++)
        {

			tempsi=matHy[i-boundx][j-boundy];
			if (is_lor_H[tempsi])
			{
				lor_a1_Ht=lor_a1_H[tempsi];
				lor_a2_Ht=lor_a2_H[tempsi];
				lor_a3_Ht=lor_a3_H[tempsi];


            	for(k=boundz;k<zk;k++)
				{
					Hy[i][j][k]+=lor_a1_Ht*lor_Hy_n_2[p]-lor_a2_Ht*lor_Jhy[p]-lor_a3_Ht*lor_Jhy_n_1[p];

					p++;
				}
			}
		}
    }

}

void update_lorJhy()
{
	short int tempsi;

	typ_prec lor_alf_Ht,lor_ksi_Ht;
	typ_pola temp;
    p=0;
    for(i=lor_Hy_is;i<=lor_Hy_ik;i++)
    {
        for(j=lor_Hy_js;j<=lor_Hy_jk;j++)
        {

			tempsi=matHy[i-boundx][j-boundy];
			if (is_lor_H[tempsi])
			{

				lor_alf_Ht=lor_alf_H[tempsi];
				lor_ksi_Ht=lor_ksi_H[tempsi];

            	for(k=boundz;k<zk;k++)
				{
					temp=lor_Jhy[p];
					lor_Jhy[p]=lor_alf_Ht*lor_Jhy[p]+lor_ksi_Ht*lor_Jhy_n_1[p]+(Hy[i][j][k]-lor_Hy_n_2[p]);

					lor_Jhy_n_1[p]=temp;

					lor_Hy_n_2[p]=lor_Hy_n_1[p];
					lor_Hy_n_1[p]=Hy[i][j][k];

					p++;
				}
			}
		}
    }

}


void update_lorHz()
{
	short int tempsi;
	typ_prec lor_a1_Ht,lor_a2_Ht,lor_a3_Ht;


    p=0;
    for(i=lor_Hz_is;i<=lor_Hz_ik;i++)
    {
        for(j=lor_Hz_js;j<=lor_Hz_jk;j++)
        {
			tempsi=matHz[i-boundx][j-boundy];
			if (is_lor_H[tempsi])
			{
				lor_a1_Ht=lor_a1_H[tempsi];
				lor_a2_Ht=lor_a2_H[tempsi];
				lor_a3_Ht=lor_a3_H[tempsi];


            	for(k=boundz;k<zk-1;k++)
				{
					Hz[i][j][k]+=lor_a1_Ht*lor_Hz_n_2[p]-lor_a2_Ht*lor_Jhz[p]-lor_a3_Ht*lor_Jhz_n_1[p];

					p++;
				}
			}
		}
    }

}

void update_lorJhz()
{
	short int tempsi;

	typ_prec lor_alf_Ht,lor_ksi_Ht;
	typ_pola temp;

    p=0;
    for(i=lor_Hz_is;i<=lor_Hz_ik;i++)
    {
        for(j=lor_Hz_js;j<=lor_Hz_jk;j++)
        {
			tempsi=matHz[i-boundx][j-boundy];
			if (is_lor_H[tempsi])
			{

				lor_alf_Ht=lor_alf_H[tempsi];
				lor_ksi_Ht=lor_ksi_H[tempsi];

            	for(k=boundz;k<zk-1;k++)
				{
					temp=lor_Jhz[p];
					lor_Jhz[p]=lor_alf_Ht*lor_Jhz[p]+lor_ksi_Ht*lor_Jhz_n_1[p]+(Hz[i][j][k]-lor_Hz_n_2[p]);

					lor_Jhz_n_1[p]=temp;

					lor_Hz_n_2[p]=lor_Hz_n_1[p];
					lor_Hz_n_1[p]=Hz[i][j][k];

					p++;
				}
			}
		}
    }

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void update_dbyE3x()
{
	short int tempsi;

	// indeks mat. debeya
	p=0;
    for(i=dby_Ex_is;i<=dby_Ex_ik;i++)
    {
        for(j=dby_Ex_js;j<=dby_Ex_jk;j++)
        {
            for(k=dby_Ex_ks;k<=dby_Ex_kk;k++)
			{
				tempsi=matEx3[i-boundx][j-boundy][k-boundz];
				if (is_dby_E[tempsi])
				{
					Ex[i][j][k]+=-dby_J_E[tempsi]*dby_Jex[p];
					p++;
				}
			}
		}
    }

}

void update_dbyJ3x()
{
	short int tempsi;

	// indeks mat. debeya
	p=0;
    for(i=dby_Ex_is;i<=dby_Ex_ik;i++)
    {
        for(j=dby_Ex_js;j<=dby_Ex_jk;j++)
        {
            for(k=dby_Ex_ks;k<=dby_Ex_kk;k++)
			{
				tempsi=matEx3[i-boundx][j-boundy][k-boundz];
				if (is_dby_E[tempsi])
				{
					dby_Jex[p]=dby_ka_E[tempsi]*dby_Jex[p]+(Ex[i][j][k]-dby_Ex_n_1[p]);
					dby_Ex_n_1[p]=Ex[i][j][k];
					p++;
				}
			}
		}
    }

}


void update_dbyE3y()
{
	short int tempsi;

	p=0;
    for(i=dby_Ey_is;i<=dby_Ey_ik;i++)
    {
        for(j=dby_Ey_js;j<=dby_Ey_jk;j++)
        {
            for(k=dby_Ey_ks;k<=dby_Ey_kk;k++)
			{
				tempsi=matEy3[i-boundx][j-boundy][k-boundz];
				if (is_dby_E[tempsi])
				{
					Ey[i][j][k]+=-dby_J_E[tempsi]*dby_Jey[p];
					p++;
				}
			}
		}
    }
}

void update_dbyJ3y()
{
	short int tempsi;

	p=0;
    for(i=dby_Ey_is;i<=dby_Ey_ik;i++)
    {
        for(j=dby_Ey_js;j<=dby_Ey_jk;j++)
        {
            for(k=dby_Ey_ks;k<=dby_Ey_kk;k++)
			{
				tempsi=matEy3[i-boundx][j-boundy][k-boundz];
				if (is_dby_E[tempsi])
				{
					dby_Jey[p]=dby_ka_E[tempsi]*dby_Jey[p]+(Ey[i][j][k]-dby_Ey_n_1[p]);
					dby_Ey_n_1[p]=Ey[i][j][k];
					p++;
				}
			}
		}
    }
}


void update_dbyE3z()
{
	short int tempsi;

	p=0;
    for(i=dby_Ez_is;i<=dby_Ez_ik;i++)
    {
        for(j=dby_Ez_js;j<=dby_Ez_jk;j++)
        {
            for(k=dby_Ez_ks;k<=dby_Ez_kk;k++)
			{
				tempsi=matEz3[i-boundx][j-boundy][k-boundz];
				if (is_dby_E[tempsi])
				{
					Ez[i][j][k]+=-dby_J_E[tempsi]*dby_Jez[p];
					p++;
				}
			}
		}
    }

}

void update_dbyJ3z()
{
	short int tempsi;

	p=0;
    for(i=dby_Ez_is;i<=dby_Ez_ik;i++)
    {
        for(j=dby_Ez_js;j<=dby_Ez_jk;j++)
        {
            for(k=dby_Ez_ks;k<=dby_Ez_kk;k++)
			{
				tempsi=matEz3[i-boundx][j-boundy][k-boundz];
				if (is_dby_E[tempsi])
				{
					dby_Jez[p]=dby_ka_E[tempsi]*dby_Jez[p]+(Ez[i][j][k]-dby_Ez_n_1[p]);
					dby_Ez_n_1[p]=Ez[i][j][k];
					p++;
				}
			}
		}
    }

}


void update_dbyEx()
{
	short int tempsi;
	typ_prec tmpv;

	// indeks mat. debeye
    p=0;
    for(i=dby_Ex_is;i<=dby_Ex_ik;i++)
    {
        for(j=dby_Ex_js;j<=dby_Ex_jk;j++)
        {
			tempsi=matEx[i-boundx][j-boundy];
			if (is_dby_E[tempsi])
			{
				tmpv=dby_J_E[tempsi];
				for(k=boundz;k<=zk;k++)
				{
					Ex[i][j][k]+=-tmpv*dby_Jex[p];
					p++;
				}
			}
			if (p>=n_dby_Ex) break;
		}
    }
}

void update_dbyJx()
{
	short int tempsi;


	// indeks mat. debeye
    p=0;
    for(i=dby_Ex_is;i<=dby_Ex_ik;i++)
    {
        for(j=dby_Ex_js;j<=dby_Ex_jk;j++)
        {
			tempsi=matEx[i-boundx][j-boundy];
			if (is_dby_E[tempsi])
			{

				for(k=boundz;k<=zk;k++)
				{
					dby_Jex[p]=dby_ka_E[tempsi]*dby_Jex[p]+(Ex[i][j][k]-dby_Ex_n_1[p]);
					dby_Ex_n_1[p]=Ex[i][j][k];
					p++;
				}
			}
			if (p>=n_dby_Ex) break;
		}
    }
}


void update_dbyEy()
{
	short int tempsi;
	typ_prec tmpv;
	p=0;
    for(i=dby_Ey_is;i<=dby_Ey_ik;i++)
    {
        for(j=dby_Ey_js;j<=dby_Ey_jk;j++)
        {
			tempsi=matEy[i-boundx][j-boundy];
			if (is_dby_E[tempsi])
			{
				tmpv=dby_J_E[tempsi];
				for(k=boundz;k<=zk;k++)
				{
					Ey[i][j][k]+=-tmpv*dby_Jey[p];
					p++;
				}
			}
			if (p>=n_dby_Ey) break;
		}
    }
}

void update_dbyJy()
{
	short int tempsi;

	p=0;
    for(i=dby_Ey_is;i<=dby_Ey_ik;i++)
    {
        for(j=dby_Ey_js;j<=dby_Ey_jk;j++)
        {
			tempsi=matEy[i-boundx][j-boundy];
			if (is_dby_E[tempsi])
			{

				for(k=boundz;k<=zk;k++)
				{
					dby_Jey[p]=dby_ka_E[tempsi]*dby_Jey[p]+(Ey[i][j][k]-dby_Ey_n_1[p]);
					dby_Ey_n_1[p]=Ey[i][j][k];
					p++;
				}
			}
			if (p>=n_dby_Ey) break;
		}
    }
}


void update_dbyEz()
{
	short int tempsi;
	typ_prec tmpv;

	p=0;
    for(i=dby_Ez_is;i<=dby_Ez_ik;i++)
    {
        for(j=dby_Ez_js;j<=dby_Ez_jk;j++)
        {
			tempsi=matEz[i-boundx][j-boundy];
			if (is_dby_E[tempsi])
			{
				tmpv=dby_J_E[tempsi];;
            	for(k=boundz;k<zk;k++)
				{
					Ez[i][j][k]+=-tmpv*dby_Jez[p];
					p++;
				}
			}
			if (p>=n_dby_Ez) break;
		}
    }

}

void update_dbyJz()
{
	short int tempsi;


	p=0;
    for(i=dby_Ez_is;i<=dby_Ez_ik;i++)
    {
        for(j=dby_Ez_js;j<=dby_Ez_jk;j++)
        {
			tempsi=matEz[i-boundx][j-boundy];
			if (is_dby_E[tempsi])
			{

            	for(k=boundz;k<zk;k++)
				{
					dby_Jez[p]=dby_ka_E[tempsi]*dby_Jez[p]+Ez[i][j][k]-dby_Ez_n_1[p];
					dby_Ez_n_1[p]=Ez[i][j][k];
					p++;
				}
			}
			if (p>=n_dby_Ez) break;
		}
    }

}


void update_dbyH3x()
{
	short int tempsi;
	// indeks mat. drudego
	p=0;
    for(i=dby_Hx_is;i<=dby_Hx_ik;i++)
    {
        for(j=dby_Hx_js;j<=dby_Hx_jk;j++)
        {
            for(k=dby_Hx_ks;k<=dby_Hx_kk;k++)
			{
				tempsi=matHx3[i-boundx][j-boundy][k-boundz];
				if (is_dby_H[tempsi])
				{
					Hx[i][j][k]+=-dby_J_H[tempsi]*dby_Jhx[p];
					p++;
				}
			}
		}
    }

}

void update_dbyJh3x()
{
	short int tempsi;
	// indeks mat. drudego
	p=0;
    for(i=dby_Hx_is;i<=dby_Hx_ik;i++)
    {
        for(j=dby_Hx_js;j<=dby_Hx_jk;j++)
        {
            for(k=dby_Hx_ks;k<=dby_Hx_kk;k++)
			{
				tempsi=matHx3[i-boundx][j-boundy][k-boundz];
				if (is_dby_H[tempsi])
				{
					dby_Jhx[p]=dby_ka_H[tempsi]*dby_Jhx[p]+Hx[i][j][k]-dby_Hx_n_1[p];
					dby_Hx_n_1[p]=Hx[i][j][k];
					p++;
				}
			}
		}
    }

}


void update_dbyH3y()
{
	short int tempsi;

	p=0;
    for(i=dby_Hy_is;i<=dby_Hy_ik;i++)
    {
        for(j=dby_Hy_js;j<=dby_Hy_jk;j++)
        {
            for(k=dby_Hy_ks;k<=dby_Hy_kk;k++)
			{
				tempsi=matHy3[i-boundx][j-boundy][k-boundz];

				if (is_dby_H[tempsi])
				{
					Hy[i][j][k]+=-dby_J_H[tempsi]*dby_Jhy[p];
					p++;
				}
			}
		}
    }

}

void update_dbyJh3y()
{
	short int tempsi;

	p=0;
    for(i=dby_Hy_is;i<=dby_Hy_ik;i++)
    {
        for(j=dby_Hy_js;j<=dby_Hy_jk;j++)
        {
            for(k=dby_Hy_ks;k<=dby_Hy_kk;k++)
			{
				tempsi=matHy3[i-boundx][j-boundy][k-boundz];

				if (is_dby_H[tempsi])
				{
					dby_Jhy[p]=dby_ka_H[tempsi]*dby_Jhy[p]+Hy[i][j][k]-dby_Hy_n_1[p];
					dby_Hy_n_1[p]=Hy[i][j][k];
					p++;
				}
			}
		}
    }

}


void update_dbyH3z()
{
	short int tempsi;

	p=0;
    for(i=dby_Hz_is;i<=dby_Hz_ik;i++)
    {
        for(j=dby_Hz_js;j<=dby_Hz_jk;j++)
        {
            for(k=dby_Hz_ks;k<=dby_Hz_kk;k++)
			{
				tempsi=matHz3[i-boundx][j-boundy][k-boundz];
				if (is_dby_H[tempsi])
				{
					Hz[i][j][k]+=-dby_J_H[tempsi]*dby_Jhz[p];
					p++;
				}
			}
		}
    }
}

void update_dbyJh3z()
{
	short int tempsi;

	p=0;
    for(i=dby_Hz_is;i<=dby_Hz_ik;i++)
    {
        for(j=dby_Hz_js;j<=dby_Hz_jk;j++)
        {
            for(k=dby_Hz_ks;k<=dby_Hz_kk;k++)
			{
				tempsi=matHz3[i-boundx][j-boundy][k-boundz];
				if (is_dby_H[tempsi])
				{
					dby_Jhz[p]=dby_ka_H[tempsi]*dby_Jhz[p]+(Hz[i][j][k]-dby_Hz_n_1[p]);
					dby_Hz_n_1[p]=Hz[i][j][k];
					p++;
				}
			}
		}
    }
}



void update_dbyHx()
{
	short int tempsi;
	typ_prec tmpv;
	// indeks mat. debeyea
    p=0;

    for(i=dby_Hx_is;i<=dby_Hx_ik;i++)
    {
        for(j=dby_Hx_js;j<=dby_Hx_jk;j++)
        {
			tempsi=matHx[i-boundx][j-boundy];
			if (is_dby_H[tempsi])
			{
				tmpv=dby_J_H[tempsi];
				for(k=boundz;k<zk;k++)
				{
					Hx[i][j][k]+=-tmpv*dby_Jhx[p];
					p++;
				}
			}
			if (p>=n_dby_Hx) break;
		}
    }
}

void update_dbyJhx()
{
	short int tempsi;

	// indeks mat. debeyea
    p=0;

    for(i=dby_Hx_is;i<=dby_Hx_ik;i++)
    {
        for(j=dby_Hx_js;j<=dby_Hx_jk;j++)
        {
			tempsi=matHx[i-boundx][j-boundy];
			if (is_dby_H[tempsi])
			{

				for(k=boundz;k<zk;k++)
				{
					dby_Jhx[p]=dby_ka_H[tempsi]*dby_Jhx[p]+(Hx[i][j][k]-dby_Hx_n_1[p]);
					dby_Hx_n_1[p]=Hx[i][j][k];
					p++;
				}
			}
			if (p>=n_dby_Hx) break;
		}
    }
}



void update_dbyHy()
{
	short int tempsi;
	typ_prec tmpv;

	p=0;
    for(i=dby_Hy_is;i<=dby_Hy_ik;i++)
    {
        for(j=dby_Hy_js;j<=dby_Hy_jk;j++)
        {
			tempsi=matHy[i-boundx][j-boundy];
			if (is_dby_H[tempsi])
			{
				tmpv=dby_J_H[tempsi];
				for(k=boundz;k<zk;k++)
				{
					Hy[i][j][k]+=-tmpv*dby_Jhy[p];
					p++;
				}
			}
			if (p>=n_dby_Hy) break;
		}
    }

}

void update_dbyJhy()
{
	short int tempsi;


	p=0;
    for(i=dby_Hy_is;i<=dby_Hy_ik;i++)
    {
        for(j=dby_Hy_js;j<=dby_Hy_jk;j++)
        {
			tempsi=matHy[i-boundx][j-boundy];
			if (is_dby_H[tempsi])
			{

				for(k=boundz;k<zk;k++)
				{
					dby_Jhy[p]=dby_ka_H[tempsi]*dby_Jhy[p]+(Hy[i][j][k]-dby_Hy_n_1[p]);
					dby_Hy_n_1[p]=Hy[i][j][k];
					p++;
				}
			}
			if (p>=n_dby_Hy) break;
		}
    }

}


void update_dbyHz()
{
	short int tempsi;
	typ_prec tmpv;

	p=0;
    for(i=dby_Hz_is;i<=dby_Hz_ik;i++)
    {
        for(j=dby_Hz_js;j<=dby_Hz_jk;j++)
        {
			tempsi=matHz[i-boundx][j-boundy];
			if (is_dby_H[tempsi])
			{
				tmpv=dby_J_H[tempsi];
            	for(k=boundz;k<zk-1;k++)
				{
					Hz[i][j][k]+=-tmpv*dby_Jhz[p];
					p++;
				}
			}
			if (p>=n_dby_Hz) break;
		}
    }

}

void update_dbyJhz()
{
	short int tempsi;

	p=0;
    for(i=dby_Hz_is;i<=dby_Hz_ik;i++)
    {
        for(j=dby_Hz_js;j<=dby_Hz_jk;j++)
        {
			tempsi=matHz[i-boundx][j-boundy];
			if (is_dby_H[tempsi])
			{
            	for(k=boundz;k<zk-1;k++)
				{
					dby_Jhz[p]=dby_ka_H[tempsi]*dby_Jhz[p]+(Hz[i][j][k]-dby_Hz_n_1[p]);
					dby_Hz_n_1[p]=Hz[i][j][k];
					p++;
				}
			}
			if (p>=n_dby_Hz) break;
		}
    }

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// dodatkowe funkcje do zrodla
// zrodlo sinusoidalne (ramped sinus - max po trzech okresach)
typ_pola sinus(typ_prec shift=0.0,int n=0) {return ((t-shift)>max_t_shift[n] ? ((freq_cst[n]*(t-max_t_shift[n]))<6*M_PI ? ampl[n]*sin(freq_cst[n]*(t-shift-max_t_shift[n]))*sin(freq_cst[n]*(t-shift-max_t_shift[n])/12) : ampl[n]*sin(freq_cst[n]*(t-shift-max_t_shift[n]))) : 0);}

// impuls gaussowski
typ_pola gauss(typ_prec shift=0.0,int n=0) {return ampl[n]*exp(-freq_cst[n]*pow(t-shift-max_t_shift[n],2));}

// impuls singauss
typ_pola singauss(typ_prec shift=0.0,int n=0) {return ampl[n]*sin(freq_cst[n]*(t-shift-max_t_shift[n]))*exp(-freq_cst2[n]*pow(t-shift-max_t_shift[n],2));}

// impuls diffgauss
typ_pola diffgauss(typ_prec shift=0.0,int n=0) {return ampl[n]*(freq_cst2[n]*(t-shift-max_t_shift[n]))*exp(-(freq_cst[n]*pow(t-shift-max_t_shift[n],2)-0.5));}

// sinus (stary ,bez lagodnego startu)
typ_pola sinus_old(typ_prec shift=0.0,int n=0) {return ((t-shift)>max_t_shift[n] ? ampl[n]*sin(freq_cst[n]*(t-shift-max_t_shift[n])) : 0);}

#ifdef COMPLEX_FIELD
// zrodlo expikx (ramped sinus - max po trzech okresach)
typ_pola sinus_complex(typ_prec shift=0.0,int n=0) {return ((t-shift)>max_t_shift[n] ? ((freq_cst[n]*(t-max_t_shift[n]))<6*M_PI ? ampl[n]*std::polar((typ_prec)1,(typ_prec)(freq_cst[n]*(t-shift-max_t_shift[n])))*sin(freq_cst[n]*(t-shift-max_t_shift[n])/12) : ampl[n]*std::polar((typ_prec)1,(typ_prec)(freq_cst[n]*(t-shift-max_t_shift[n])))) : 0);}
// impuls singauss
typ_pola singauss_complex(typ_prec shift=0.0,int n=0) {return ampl[n]*std::polar((typ_prec)1,(typ_prec)(freq_cst[n]*(t-shift-max_t_shift[n])))*exp(-freq_cst2[n]*pow(t-shift-max_t_shift[n],2));}
#else
// == rzeczywiste odpowiedniki
typ_pola sinus_complex(typ_prec shift=0.0,int n=0) {return ((t-shift)>max_t_shift[n] ? ((freq_cst[n]*(t-max_t_shift[n]))<6*M_PI ? ampl[n]*sin(freq_cst[n]*(t-shift-max_t_shift[n]))*sin(freq_cst[n]*(t-shift-max_t_shift[n])/12) : ampl[n]*sin(freq_cst[n]*(t-shift-max_t_shift[n]))) : 0);}
typ_pola singauss_complex(typ_prec shift=0.0,int n=0) {return ampl[n]*sin(freq_cst[n]*(t-shift-max_t_shift[n]))*exp(-freq_cst2[n]*pow(t-shift-max_t_shift[n],2));}

#endif

//double 
typ_pola time_dependency_from_file(typ_prec shift=0.0,int n=0) 
{
       if ((t-shift)< max_t_shift[n]) return 0;
       if (t>=time_dependency_length[n]+shift+max_t_shift[n]) return 0;
       int tmp=(int)floor(t-shift-max_t_shift[n]);
       typ_prec d=t-shift-max_t_shift[n]-tmp;
       return ampl[n]*((1.0-d)*time_dependency[n][tmp]+d*time_dependency[n][tmp+1]);
}

// 1D FDTD pomocnicze dla zadawanego zrodla
// propagacja fali plaskiej na siatce dyskretnej pod pewnym katem
// od ktorego zalezy angle_ratio
// update 1D FDTD - pole elektryczne
void update_i_E(int s)
{
     // tutak check if needed
     i_E[s][0]=source_func[s](0.0,s);
     temp=depst[bound_mat]/angle_ratio[s];
     for(i=1;i<pulse_length+1;i++) i_E[s][i]+=temp*(i_H[s][i-1]-i_H[s][i]);
}
// update 1D FDTD -  pole magnetyczne
void update_i_H(int s)
{
     // tutaj check if needed
     temp=dmit[bound_mat]/angle_ratio[s];
     for(i=0;i<pulse_length+1;i++) i_H[s][i]+=temp*(i_E[s][i]-i_E[s][i+1]);
}

// funkcja pozwalajace obliczyc stosunek vi(fi,teta)/vi(0,0)
// stala k dla propagacji pod katami fi,teta
typ_prec calc_k_angle(int s)
{

    const double A=dr*tr_x[s]/2;
    const double B=dr*tr_y[s]/2;
    const double C=dr*tr_z[s]/2;
    const double D=pow((dr/(dt*light_speed))*sin(freq_cst[s]/2),2);

    typ_prec tmp=2*M_PI/lambda[s];
    for(i=0;i<100;i++) tmp-=(pow(sin(A*tmp),2)+pow(sin(B*tmp),2)+pow(sin(C*tmp),2)-D)/(A*sin(2*A*tmp)+B*sin(2*B*tmp)+C*sin(2*C*tmp));
    return tmp;
}

// stala k dla propagacji pod katem 0,0
typ_prec calc_k(int s)
{
    const double D=pow((dr/(dt*light_speed))*sin(freq_cst[s]/2),2);

    typ_prec tmp=2*M_PI/lambda[s];
    for(i=0;i<100;i++) tmp-=2*(pow(sin(dr*tmp/2),2)-D)/(dr*sin(dr*tmp));
    return tmp;
}

// przejscie z ukladu symulacji do ukladu zrodla
inline typ_prec euler_2x(typ_prec x,typ_prec y,typ_prec z,int s) {return (x-i_source[s])*tr_xx[s]+(y-j_source[s])*tr_yy[s]-(z-k_source[s])*tr_sin_teta[s];}
inline typ_prec euler_2y(typ_prec x,typ_prec y,typ_prec z,int s) {return (y-j_source[s])*tr_cos_fi[s]-(x-i_source[s])*tr_sin_fi[s];}
inline typ_prec euler_2z(typ_prec x,typ_prec y,typ_prec z,int s) {return (x-i_source[s])*tr_x[s]+(y-j_source[s])*tr_y[s]+(z-k_source[s])*tr_z[s];}

// obliczanie stalych trygonometrycznych
void calc_trig()
{
    for(int s=0;s<n_SRCS;s++)
    {
      if (source_type[s]<3)
      {
      tr_x[s]=cos(fi[s])*sin(teta[s]);
      tr_y[s]=sin(fi[s])*sin(teta[s]);
      tr_z[s]=cos(teta[s]);

	  tr_sin_teta[s]=sin(teta[s]);
	  tr_sin_fi[s]=sin(fi[s]);
	  tr_cos_fi[s]=cos(fi[s]);

      tr_xx[s]=cos(fi[s])*cos(teta[s]);
      tr_yy[s]=sin(fi[s])*cos(teta[s]);
      tr_hx[s]=sin(psi[s])*sin(fi[s])+cos(psi[s])*cos(teta[s])*cos(fi[s]);
      tr_hy[s]=-sin(psi[s])*cos(fi[s])+cos(psi[s])*cos(teta[s])*sin(fi[s]);
      tr_hz[s]=-cos(psi[s])*sin(teta[s]);

      tr_ex[s]=cos(psi[s])*sin(fi[s])-sin(psi[s])*cos(teta[s])*cos(fi[s]);
      tr_ey[s]=-cos(psi[s])*cos(fi[s])-sin(psi[s])*cos(teta[s])*sin(fi[s]);
      tr_ez[s]=sin(psi[s])*sin(teta[s]);

    	temp=calc_k(s);
    	angle_ratio[s]=calc_k_angle(s)/temp;
        ampl[s]*=(angle_ratio[s]*angle_ratio[s]);
    }
	if (source_type[s]==0)
    {
		// ustala odpowiednie parametry do liczenia d
    	i_source[s]=(((fi[s]<M_PI/2)||(fi[s]>1.5*M_PI)) ? (i_0-1):(i_1+1));
    	j_source[s]=((fi[s]<M_PI) ? (j_0-1):(j_1+1));
		k_source[s]=((teta[s]<M_PI/2) ? (k_0-1):(k_1+1));

	}

	if (source_type[s]==1 || source_type[s]==2)
	{

		st_ratio[s]=sqrt(mi[bound_mat]*eps[bound_mat])*dr/(light_speed*dt)/angle_ratio[s];

		max_t_shift[s]=SRCS_shift[s]+t_0[s]-st_ratio[s]*euler_2z((((fi[s]<M_PI/2)||(fi[s]>1.5*M_PI)) ? (i_0-1):(i_1+1)),((fi[s]<M_PI) ? (j_0-1):(j_1+1)),((teta[s]<M_PI/2) ? (k_0-1):(k_1+1)),s);

        // parametry gaussa
        // szybkie stworznie tablicy wiercho³ków ijej posortowanie wzgl. odleg³oœci od Ÿród³a
        {
         typ_prec tab_corner[]={euler_2z(i_0-1,j_0-1,k_0-1,s),euler_2z(i_0-1,j_0-1,k_1+1,s),euler_2z(i_0-1,j_1+1,k_1+1,s),euler_2z(i_0-1,j_1+1,k_0-1,s),euler_2z(i_1+1,j_0-1,k_0-1,s),euler_2z(i_1+1,j_0-1,k_1+1,s),euler_2z(i_1+1,j_1+1,k_1+1,s),euler_2z(i_1+1,j_1+1,k_0-1,s)};
         for(i=0;i<8;i++) for (j=0;j<7;j++){k=j; while ((k<7)&&(tab_corner[k]>tab_corner[k+1])) {typ_prec temp=tab_corner[k];tab_corner[k]=tab_corner[k+1];tab_corner[k+1]=temp;k++;}} 
         gauss_c[s]=tab_corner[3]+1;
         }
		if (source_type[s]==1)
		{      
         gauss_z0x[s]=pow(M_PI*param1[s]/(lambda[s]/dr),2);
         gauss_z0y[s]=pow(M_PI*param2[s]/(lambda[s]/dr),2);
         gauss_amp[s]=sqrt(param1[s]+param2[s]);
         
        }
        
		if (source_type[s]==2)
		{
			// rzutowanie E wzdluz x na pozostale osie, wiec psi = -pi/2
    		tr_hxx[s]=-sin(fi[s]);
    		tr_hyx[s]=cos(fi[s]);
    		tr_hzx[s]=0;

    		tr_exx[s]=cos(teta[s])*cos(fi[s]);
    		tr_eyx[s]=cos(teta[s])*sin(fi[s]);
    		tr_ezx[s]=-sin(teta[s]);

			// rzutowanie E wzdluz y na pozostale osie , wiec psi = 0
    		tr_hxy[s]=cos(teta[s])*cos(fi[s]);
    		tr_hyy[s]=cos(teta[s])*sin(fi[s]);
    		tr_hzy[s]=-sin(teta[s]);

    		tr_exy[s]=sin(fi[s]);
    		tr_eyy[s]=-cos(fi[s]);
    		tr_ezy[s]=0;
		}
	}
	else
	{
		max_t_shift[s]=SRCS_shift[s]+t_0[s];
	}
	
    if (is_SRCS_shift_rel[s]==false) max_t_shift[s]=SRCS_shift[s];

}
	
}

// oblicznie odleglosci od zrodla == tozsame z euler_2z
inline typ_prec calc_distance(typ_prec a, typ_prec b,typ_prec c,int s) {return ((a-i_source[s])*tr_x[s]+(b-j_source[s])*tr_y[s]+(c-k_source[s])*tr_z[s]);}

// funkcje wyliczajace z 1D FDTD pola E,H dla dowolnej odleglosci d przez interpolacje
typ_pola interp_E_field(typ_prec distance,int s)
{
    int tmp=(int)floor(distance);
    typ_prec d=distance-tmp;

    return ((1.0-d)*i_E[s][tmp]+d*i_E[s][tmp+1]);
}

typ_pola interp_H_field(typ_prec distance,int s)
{
    int tmp=(int)floor(distance+0.5);
    typ_prec d=distance+0.5-tmp;

    return ((1.0-d)*i_H[s][tmp-1]+d*i_H[s][tmp]);
}

// funkcje zwracajace wartosc padajacej fali  w danym pkcie przestrzenii
inline typ_pola plane_inc_Ex(typ_prec a,typ_prec b,typ_prec c,int s) {return interp_E_field(calc_distance(a,b,c,s),s)*tr_ex[s];}
inline typ_pola plane_inc_Ey(typ_prec a,typ_prec b,typ_prec c,int s) {return interp_E_field(calc_distance(a,b,c,s),s)*tr_ey[s];}
inline typ_pola plane_inc_Ez(typ_prec a,typ_prec b,typ_prec c,int s) {return interp_E_field(calc_distance(a,b,c,s),s)*tr_ez[s];}

inline typ_pola plane_inc_Hx(typ_prec a,typ_prec b,typ_prec c, int s) {return interp_H_field(calc_distance(a,b,c,s),s)*tr_hx[s];}
inline typ_pola plane_inc_Hy(typ_prec a,typ_prec b,typ_prec c, int s) {return interp_H_field(calc_distance(a,b,c,s),s)*tr_hy[s];}
inline typ_pola plane_inc_Hz(typ_prec a,typ_prec b,typ_prec c, int s) {return interp_H_field(calc_distance(a,b,c,s),s)*tr_hz[s];}

// zrodlo gausssowskie i hermitowsko-gaussowskie
typ_prec gaussshape(typ_prec a,typ_prec b,typ_prec c,int s) 
{
	if (c>gauss_c[s]) return 0;
	
	typ_prec wzx2=1+c*c/gauss_z0x[s];
	typ_prec wzy2=1+c*c/gauss_z0y[s];

	return (gauss_amp[s]/sqrt(param1[s]*wzx2+param2[s]*wzy2))*exp(-((a*a/(param1[s]*wzx2))+(b*b/(param2[s]*wzy2))));
}

typ_prec hermite_gauss01(typ_prec a,typ_prec b,typ_prec c,int s) {return (2*b/param3[s])*gaussshape(a,b,c,s);}
typ_prec hermite_gauss10(typ_prec a,typ_prec b,typ_prec c,int s) {return (2*a/param3[s])*gaussshape(a,b,c,s);}
typ_prec hermite_gauss11(typ_prec a,typ_prec b,typ_prec c,int s) {return (4*a*b/param3[s])*gaussshape(a,b,c,s);}
typ_prec hermite_gauss02(typ_prec a,typ_prec b,typ_prec c,int s) {return ((4*b*b-2)/param3[s])*gaussshape(a,b,c,s);}
typ_prec hermite_gauss20(typ_prec a,typ_prec b,typ_prec c,int s) {return ((4*a*a-2)/param3[s])*gaussshape(a,b,c,s);}
typ_prec hermite_gauss12(typ_prec a,typ_prec b,typ_prec c,int s) {return ((2*a)*(4*b*b-2)/param3[s])*gaussshape(a,b,c,s);}
typ_prec hermite_gauss21(typ_prec a,typ_prec b,typ_prec c,int s) {return ((4*a*a-2)*(2*b)/param3[s])*gaussshape(a,b,c,s);}
typ_prec hermite_gauss22(typ_prec a,typ_prec b,typ_prec c,int s) {return ((4*a*a-2)*(4*b*b-2)/param3[s])*gaussshape(a,b,c,s);}

typ_prec bessel0(typ_prec a,typ_prec b,typ_prec c,int s) {typ_prec r=sqrt(a*a+b*b)/param1[s];if ((r<2.4048) && (c<gauss_c[s])) return y0(r);else return 0.0;}


// dowolne zrodlo (wczytywane z zewntarz)

inline typ_pola beam_inc_Ex(typ_prec a, typ_prec b,typ_prec c,int s) {return shape[s](euler_2x(a,b,c,s),euler_2y(a,b,c,s),euler_2z(a,b,c,s),s)*source_func[s](euler_2z(a,b,c,s)*st_ratio[s]-0.5,s)*tr_ex[s];}
inline typ_pola beam_inc_Ey(typ_prec a, typ_prec b,typ_prec c,int s) {return shape[s](euler_2x(a,b,c,s),euler_2y(a,b,c,s),euler_2z(a,b,c,s),s)*source_func[s](euler_2z(a,b,c,s)*st_ratio[s]-0.5,s)*tr_ey[s];}
inline typ_pola beam_inc_Ez(typ_prec a, typ_prec b,typ_prec c,int s) {return shape[s](euler_2x(a,b,c,s),euler_2y(a,b,c,s),euler_2z(a,b,c,s),s)*source_func[s](euler_2z(a,b,c,s)*st_ratio[s]-0.5,s)*tr_ez[s];}

inline typ_pola beam_inc_Hx(typ_prec a, typ_prec b,typ_prec c,int s) {return shape[s](euler_2x(a,b,c,s),euler_2y(a,b,c,s),euler_2z(a,b,c,s),s)*source_func[s](euler_2z(a,b,c,s)*st_ratio[s],s)*tr_hx[s];}
inline typ_pola beam_inc_Hy(typ_prec a, typ_prec b,typ_prec c,int s) {return shape[s](euler_2x(a,b,c,s),euler_2y(a,b,c,s),euler_2z(a,b,c,s),s)*source_func[s](euler_2z(a,b,c,s)*st_ratio[s],s)*tr_hy[s];}
inline typ_pola beam_inc_Hz(typ_prec a, typ_prec b,typ_prec c,int s) {return shape[s](euler_2x(a,b,c,s),euler_2y(a,b,c,s),euler_2z(a,b,c,s),s)*source_func[s](euler_2z(a,b,c,s)*st_ratio[s],s)*tr_hz[s];}

inline double shape_interp(double*** shape_tab,typ_prec a,typ_prec b,int s)
{
	if (fabs(a)<sh_cx[s] && fabs(b)<sh_cy[s])
	{
		a=(a+sh_cx[s])/drs[s];b=(b+sh_cy[s])/drs[s];int a1=(int)ceil(a);int a2=(int)floor(a);int b1=(int)ceil(b);int b2=(int)floor(b);
		return 0.25*(shape_tab[s][a1][b1]+shape_tab[s][a1][b2]+shape_tab[s][a2][b1]+shape_tab[s][a2][b2]);
	}
	return 0;
}

inline typ_prec shape_x_real(typ_prec a,typ_prec b,int s) {return shape_interp(shape_ex_real,a,b,s);}
inline typ_prec shape_x_imag(typ_prec a,typ_prec b,int s) {return shape_interp(shape_ex_imag,a,b,s);}
inline typ_prec shape_y_real(typ_prec a,typ_prec b,int s) {return shape_interp(shape_ey_real,a,b,s);}
inline typ_prec shape_y_imag(typ_prec a,typ_prec b,int s) {return shape_interp(shape_ey_imag,a,b,s);}

inline typ_pola obeam_inc_Ex(typ_prec a, typ_prec b,typ_prec c,int s)
{
	typ_prec xe=euler_2x(a,b,c,s);typ_prec ye=euler_2y(a,b,c,s);typ_prec zet=euler_2z(a,b,c,s)*st_ratio[s];
	return ((zet<gauss_c[s]) ? (shape_x_real(xe,ye,s)*tr_exx[s]+shape_y_real(xe,ye,s)*tr_exy[s])*source_func[s](zet-0.5,s)+(shape_x_imag(xe,ye,s)*tr_exx[s]+shape_y_imag(xe,ye,s)*tr_exy[s])*source_func[s](zet+half_period[s]-0.5,s) : 0);
}
inline typ_pola obeam_inc_Ey(typ_prec a, typ_prec b,typ_prec c,int s)
{
	typ_prec xe=euler_2x(a,b,c,s);typ_prec ye=euler_2y(a,b,c,s);typ_prec zet=euler_2z(a,b,c,s)*st_ratio[s];
	return ((zet<gauss_c[s]) ? (shape_x_real(xe,ye,s)*tr_eyx[s]+shape_y_real(xe,ye,s)*tr_eyy[s])*source_func[s](zet-0.5,s)+(shape_x_imag(xe,ye,s)*tr_eyx[s]+shape_y_imag(xe,ye,s)*tr_eyy[s])*source_func[s](zet+half_period[s]-0.5,s) : 0);
}
inline typ_pola obeam_inc_Ez(typ_prec a, typ_prec b,typ_prec c,int s)
{
	typ_prec xe=euler_2x(a,b,c,s);typ_prec ye=euler_2y(a,b,c,s);typ_prec zet=euler_2z(a,b,c,s)*st_ratio[s];
	return ((zet<gauss_c[s]) ? (shape_x_real(xe,ye,s)*tr_ezx[s]+shape_y_real(xe,ye,s)*tr_ezy[s])*source_func[s](zet-0.5,s)+(shape_x_imag(xe,ye,s)*tr_ezx[s]+shape_y_imag(xe,ye,s)*tr_ezy[s])*source_func[s](zet+half_period[s]-0.5,s) : 0);
}

inline typ_pola obeam_inc_Hx(typ_prec a, typ_prec b,typ_prec c,int s)
{
	typ_prec xe=euler_2x(a,b,c,s);typ_prec ye=euler_2y(a,b,c,s);typ_prec zet=euler_2z(a,b,c,s)*st_ratio[s];
	return ((zet<gauss_c[s]) ? (shape_x_real(xe,ye,s)*tr_hxx[s]+shape_y_real(xe,ye,s)*tr_hxy[s])*source_func[s](zet,s)+(shape_x_imag(xe,ye,s)*tr_hxx[s]+shape_y_imag(xe,ye,s)*tr_hxy[s])*source_func[s](zet+half_period[s],s) : 0);
}
inline typ_pola obeam_inc_Hy(typ_prec a, typ_prec b,typ_prec c,int s)
{
	typ_prec xe=euler_2x(a,b,c,s);typ_prec ye=euler_2y(a,b,c,s);typ_prec zet=euler_2z(a,b,c,s)*st_ratio[s];
	return ((zet<gauss_c[s]) ? (shape_x_real(xe,ye,s)*tr_hyx[s]+shape_y_real(xe,ye,s)*tr_hyy[s])*source_func[s](zet,s)+(shape_x_imag(xe,ye,s)*tr_hyx[s]+shape_y_imag(xe,ye,s)*tr_hyy[s])*source_func[s](zet+half_period[s],s) : 0);
}
inline typ_pola obeam_inc_Hz(typ_prec a, typ_prec b,typ_prec c,int s)
{
	typ_prec xe=euler_2x(a,b,c,s);typ_prec ye=euler_2y(a,b,c,s);typ_prec zet=euler_2z(a,b,c,s)*st_ratio[s];
	return ((zet<gauss_c[s]) ? (shape_x_real(xe,ye,s)*tr_hzx[s]+shape_y_real(xe,ye,s)*tr_hzy[s])*source_func[s](zet,s)+(shape_x_imag(xe,ye,s)*tr_hzx[s]+shape_y_imag(xe,ye,s)*tr_hzy[s])*source_func[s](zet+half_period[s],s) : 0);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
void electric_dipolx(int s) {Ex[floor_i_source[s]][floor_j_source[s]][floor_k_source[s]]=source_func[s](0.0,s);}
void electric_dipoly(int s) {Ey[floor_i_source[s]][floor_j_source[s]][floor_k_source[s]]=source_func[s](0.0,s);}
void electric_dipolz(int s) {Ez[floor_i_source[s]][floor_j_source[s]][floor_k_source[s]]=source_func[s](0.0,s);}

void magnetic_dipolx(int s) {Hx[floor_i_source[s]][floor_j_source[s]][floor_k_source[s]]=source_func[s](0.0,s);}
void magnetic_dipoly(int s) {Hy[floor_i_source[s]][floor_j_source[s]][floor_k_source[s]]=source_func[s](0.0,s);}
void magnetic_dipolz(int s) {Hz[floor_i_source[s]][floor_j_source[s]][floor_k_source[s]]=source_func[s](0.0,s);}


// zadawanie pola zgodnie z TF/SF
void TF_SF_bac_H_Hx_Ey(int s)
{
    temp=dmit[bound_mat];
    // check
    for(i=i_0-1;i<i_1;i++) for(j=j_0;j<j_1;j++) Hx[i][j][k_0-1]+=-temp*inc_Ey[s](i+1,j+0.5,k_0,s);
}
void TF_SF_bac_H_Hy_Ex(int s)
{
    temp=dmit[bound_mat];
   
    for(i=i_0;i<i_1;i++) for(j=j_0-1;j<j_1;j++) Hy[i][j][k_0-1]+=temp*inc_Ex[s](i+0.5,j+1,k_0,s);
}
void TF_SF_fro_H_Hx_Ey(int s)
{
    temp=dmit[bound_mat];
    for(i=i_0-1;i<i_1;i++) for(j=j_0;j<j_1;j++) Hx[i][j][k_1]+=temp*inc_Ey[s](i+1,j+0.5,k_1,s);
}

void TF_SF_fro_H_Hy_Ex(int s)
{
    temp=dmit[bound_mat];
	 for(i=i_0;i<i_1;i++) for(j=j_0-1;j<j_1;j++) Hy[i][j][k_1]+=-temp*inc_Ex[s](i+0.5,j+1,k_1,s);
}

void TF_SF_bot_H_Hx_Ez(int s)
{
    temp=dmit[bound_mat];
    for(i=i_0-1;i<i_1;i++) for(k=k_0;k<k_1;k++) Hx[i][j_0-1][k]+=temp*inc_Ez[s](i+1,j_0,k+0.5,s);
}

void TF_SF_bot_H_Hz_Ex(int s)
{
    temp=dmit[bound_mat];
    for(i=i_0;i<i_1;i++) for(k=k_0-1;k<k_1;k++) Hz[i][j_0-1][k]+=-temp*inc_Ex[s](i+0.5,j_0,k+1,s);
}

void TF_SF_top_H_Hx_Ez(int s)
{
    temp=dmit[bound_mat];
    for(i=i_0-1;i<i_1;i++) for(k=k_0;k<k_1;k++) Hx[i][j_1][k]+=-temp*inc_Ez[s](i+1,j_1,k+0.5,s);
}

void TF_SF_top_H_Hz_Ex(int s)
{
    temp=dmit[bound_mat];
    for(i=i_0;i<i_1;i++) for(k=k_0-1;k<k_1;k++) Hz[i][j_1][k]+=temp*inc_Ex[s](i+0.5,j_1,k+1,s);
}

void TF_SF_lef_H_Hy_Ez(int s)
{
	temp=dmit[bound_mat];
    for(j=j_0-1;j<j_1;j++) for(k=k_0;k<k_1;k++) Hy[i_0-1][j][k]+=-temp*inc_Ez[s](i_0,j+1,k+0.5,s);
}

void TF_SF_lef_H_Hz_Ey(int s)
{
	temp=dmit[bound_mat];
    for(j=j_0;j<j_1;j++) for(k=k_0-1;k<k_1;k++) Hz[i_0-1][j][k]+=temp*inc_Ey[s](i_0,j+0.5,k+1,s);
}
void TF_SF_rig_H_Hy_Ez(int s)
{
	temp=dmit[bound_mat];
    for(j=j_0-1;j<j_1;j++) for(k=k_0;k<k_1;k++) Hy[i_1][j][k]+=temp*inc_Ez[s](i_1,j+1,k+0.5,s);
}
void TF_SF_rig_H_Hz_Ey(int s)
{
	temp=dmit[bound_mat];
    for(j=j_0;j<j_1;j++) for(k=k_0-1;k<k_1;k++) Hz[i_1][j][k]+=-temp*inc_Ey[s](i_1,j+0.5,k+1,s);
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void TF_SF_bac_E_Ex_Hy(int s)
{
    temp=depst[bound_mat];
    for(i=i_0;i<i_1;i++) for(j=j_0;j<=j_1;j++) Ex[i][j][k_0]+=temp*inc_Hy[s](i+0.5,j,k_0-0.5,s);
}

void TF_SF_bac_E_Ey_Hx(int s)
{
    temp=depst[bound_mat];
    for(i=i_0;i<=i_1;i++) for(j=j_0;j<j_1;j++) Ey[i][j][k_0]+=-temp*inc_Hx[s](i,j+0.5,k_0-0.5,s);
}

void TF_SF_fro_E_Ex_Hy(int s)
{
    temp=depst[bound_mat];
    for(i=i_0;i<i_1;i++) for(j=j_0;j<=j_1;j++) Ex[i][j][k_1]+=-temp*inc_Hy[s](i+0.5,j,k_1+0.5,s);
}

void TF_SF_fro_E_Ey_Hx(int s)
{
    temp=depst[bound_mat];
    for(i=i_0;i<=i_1;i++) for(j=j_0;j<j_1;j++) Ey[i][j][k_1]+=temp*inc_Hx[s](i,j+0.5,k_1+0.5,s);
}

void TF_SF_lef_E_Ey_Hz(int s)
{
    temp=depst[bound_mat];
    for(j=j_0;j<j_1;j++) for(k=k_0;k<=k_1;k++) Ey[i_0][j][k]+=temp*inc_Hz[s](i_0-0.5,j+0.5,k,s);
}

void TF_SF_lef_E_Ez_Hy(int s)
{
    temp=depst[bound_mat];
    for(j=j_0;j<=j_1;j++) for(k=k_0;k<k_1;k++) Ez[i_0][j][k]+=-temp*inc_Hy[s](i_0-0.5,j,k+0.5,s);
}

void TF_SF_rig_E_Ey_Hz(int s)
{
    temp=depst[bound_mat];
    for(j=j_0;j<j_1;j++) for(k=k_0;k<=k_1;k++) Ey[i_1][j][k]+=-temp*inc_Hz[s](i_1+0.5,j+0.5,k,s);
}

void TF_SF_rig_E_Ez_Hy(int s)
{
    temp=depst[bound_mat];
    for(j=j_0;j<=j_1;j++) for(k=k_0;k<k_1;k++) Ez[i_1][j][k]+=temp*inc_Hy[s](i_1+0.5,j,k+0.5,s);
}

void TF_SF_bot_E_Ex_Hz(int s)
{
	temp=depst[bound_mat];
    for(i=i_0;i<i_1;i++) for(k=k_0;k<=k_1;k++) Ex[i][j_0][k]+=-temp*inc_Hz[s](i+0.5,j_0-0.5,k,s);
}

void TF_SF_bot_E_Ez_Hx(int s)
{
	temp=depst[bound_mat];
    for(i=i_0;i<=i_1;i++) for(k=k_0;k<k_1;k++) Ez[i][j_0][k]+=temp*inc_Hx[s](i,j_0-0.5,k+0.5,s);
}

void TF_SF_top_E_Ez_Hx(int s)
{
	temp=depst[bound_mat];
    for(i=i_0;i<=i_1;i++) for(k=k_0;k<k_1;k++) Ez[i][j_1][k]+=-temp*inc_Hx[s](i,j_1+0.5,k+0.5,s);
}

void TF_SF_top_E_Ex_Hz(int s)
{
	temp=depst[bound_mat];
    for(i=i_0;i<i_1;i++) for(k=k_0;k<=k_1;k++) Ex[i][j_1][k]+=temp*inc_Hz[s](i+0.5,j_1+0.5,k,s);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// periodyczne warunki brzegowe
void periodic_bacfro_Ex()
{
    for(i=0;i<xs;i++) for(j=0;j<=ys;j++) for(k=0;k<k_0;k++)  {Ex[i][j][k_1+k+1]=Ex[i][j][k_0+k+1];Ex[i][j][k_0-k-1]=Ex[i][j][k_1-k-1];}
}
void periodic_bacfro_Ey()
{
        for(i=0;i<=xs;i++) for(j=0;j<ys;j++) for(k=0;k<k_0;k++)  {Ey[i][j][k_1+k+1]=Ey[i][j][k_0+k+1];Ey[i][j][k_0-k-1]=Ey[i][j][k_1-k-1];}
}
void periodic_bacfro_Ez()
{
	for(i=0;i<=xs;i++) for(j=0;j<=ys;j++) for(k=0;k<k_0;k++) {Ez[i][j][k_1+k]=Ez[i][j][k_0+k];Ez[i][j][k_0-k-1]=Ez[i][j][k_1-k-1];}
}

void periodic_bacfro_Hx()
{
	for(i=0;i<xs-1;i++) for(j=0;j<ys;j++) for(k=0;k<k_0;k++) {Hx[i][j][k_1+k]=Hx[i][j][k_0+k];Hx[i][j][k_0-k-1]=Hx[i][j][k_1-k-1];}
}
void periodic_bacfro_Hy()
{
	for(i=0;i<xs;i++) for(j=0;j<ys-1;j++) for(k=0;k<k_0;k++) {Hy[i][j][k_1+k]=Hy[i][j][k_0+k];Hy[i][j][k_0-k-1]=Hy[i][j][k_1-k-1];}
}
void periodic_bacfro_Hz()
{
	for(i=0;i<xs;i++) for(j=0;j<ys;j++) for(k=0;k<k_0-1;k++) {Hz[i][j][k_1+k]=Hz[i][j][k_0+k];Hz[i][j][k_0-k-2]=Hz[i][j][k_1-k-2];}
}

void periodic_lefrig_Ex()
{
    for(i=0;i<i_0;i++) for(j=0;j<=ys;j++) for(k=0;k<=zs;k++) {Ex[i_1+i][j][k]=Ex[i_0+i][j][k];Ex[i_0-i-1][j][k]=Ex[i_1-i-1][j][k];}
}	
void periodic_lefrig_Ey()
{
    for(i=0;i<i_0;i++) for(j=0;j<ys;j++)  for(k=0;k<=zs;k++) {Ey[i_1+i+1][j][k]=Ey[i_0+i+1][j][k];Ey[i_0-i-1][j][k]=Ey[i_1-i-1][j][k];}
}
void periodic_lefrig_Ez()
{
    for(i=0;i<i_0;i++) for(j=0;j<=ys;j++) for(k=0;k<zs;k++)  {Ez[i_1+i+1][j][k]=Ez[i_0+i+1][j][k];Ez[i_0-i-1][j][k]=Ez[i_1-i-1][j][k];}
}

void periodic_lefrig_Hx()
{
    for(i=0;i<i_0-1;i++) for(j=0;j<ys;j++) for(k=0;k<zs;k++) {Hx[i_1+i][j][k]=Hx[i_0+i][j][k];Hx[i_0-i-2][j][k]=Hx[i_1-i-2][j][k];}
}
void periodic_lefrig_Hy()
{
	for(i=0;i<i_0;i++) for(j=0;j<ys-1;j++) for(k=0;k<zs;k++) {Hy[i_1+i][j][k]=Hy[i_0+i][j][k];Hy[i_0-i-1][j][k]=Hy[i_1-i-1][j][k];}
}
void periodic_lefrig_Hz()
{
	for(i=0;i<i_0;i++) for(j=0;j<ys;j++) for(k=0;k<zs-1;k++) {Hz[i_1+i][j][k]=Hz[i_0+i][j][k];Hz[i_0-i-1][j][k]=Hz[i_1-i-1][j][k];}
}
void periodic_bottop_Ex()
{
    for(i=0;i<xs;i++) for(j=0;j<j_0;j++) for(k=0;k<=zs;k++)  {Ex[i][j_1+j+1][k]=Ex[i][j_0+j+1][k];Ex[i][j_0-j-1][k]=Ex[i][j_1-j-1][k];}
}
void periodic_bottop_Ey()
{
    for(i=0;i<=xs;i++) for(j=0;j<j_0;j++) for(k=0;k<=zs;k++) {Ey[i][j_1+j][k]=Ey[i][j_0+j][k];Ey[i][j_0-j-1][k]=Ey[i][j_1-j-1][k];}	
}
void periodic_bottop_Ez()
{
	for(i=0;i<=xs;i++) for(j=0;j<j_0;j++) for(k=0;k<zs;k++)  {Ez[i][j_1+j+1][k]=Ez[i][j_0+j+1][k];Ez[i][j_0-j-1][k]=Ez[i][j_1-j-1][k];}
}

void periodic_bottop_Hx()
{
	for(i=0;i<xs-1;i++) for(j=0;j<j_0;j++) for(k=0;k<zs;k++) {Hx[i][j_1+j][k]=Hx[i][j_0+j][k];Hx[i][j_0-j-1][k]=Hx[i][j_1-j-1][k];}
}
void periodic_bottop_Hy()
{
	for(i=0;i<xs;i++) for(j=0;j<j_0-1;j++) for(k=0;k<zs;k++) {Hy[i][j_1+j][k]=Hy[i][j_0+j][k];Hy[i][j_0-j-2][k]=Hy[i][j_1-j-2][k];}
}
void periodic_bottop_Hz()
{
	for(i=0;i<xs;i++) for(j=0;j<j_0;j++) for(k=0;k<zs-1;k++) {Hz[i][j_1+j][k]=Hz[i][j_0+j][k];Hz[i][j_0-j-1][k]=Hz[i][j_1-j-1][k];}
}

inline void set_periodic_lefrig() 
{
       i_0=boundx+bndx+disPML;i_1=i_0+xr;
       if (is_calc_ex) periodic_lefrig_Ex();if (is_calc_ey) periodic_lefrig_Ey();if (is_calc_ez) periodic_lefrig_Ez();
       if (is_calc_hx) periodic_lefrig_Hx();if (is_calc_hy) periodic_lefrig_Hy();if (is_calc_hz) periodic_lefrig_Hz();
}
inline void set_periodic_bottop() 
{
       j_0=boundy+bndy+disPML;j_1=j_0+yr;
       if (is_calc_ex) periodic_bottop_Ex();if (is_calc_ey) periodic_bottop_Ey();if (is_calc_ez) periodic_bottop_Ez();
       if (is_calc_hx) periodic_bottop_Hx();if (is_calc_hy) periodic_bottop_Hy();if (is_calc_hz) periodic_bottop_Hz();
}
inline void set_periodic_bacfro() 
{
       k_0=boundz+bndz+disPML;k_1=k_0+zr;
       if (is_calc_ex) periodic_bacfro_Ex();if (is_calc_ey) periodic_bacfro_Ey();if (is_calc_ez) periodic_bacfro_Ez();
       if (is_calc_hx) periodic_bacfro_Hx();if (is_calc_hy) periodic_bacfro_Hy();if (is_calc_hz) periodic_bacfro_Hz();
}

// symetryczne warunki brzegowe
void symmetric_bacfro_Ex()
{
	for(i=0;i<xs;i++) for(j=0;j<=ys;j++) for(k=0;k<=k_0-1;k++)  Ex[i][j][k_1+k+1]=wz_Ex*Ex[i][j][k_1-k-1];
	if (wz_Ex==-1) for(i=0;i<xs;i++) for(j=0;j<=ys;j++) Ex[i][j][k_1]=0;
}
void symmetric_bacfro_Ey()
{
	for(i=0;i<=xs;i++) for(j=0;j<ys;j++) for(k=0;k<=k_0-1;k++)  Ey[i][j][k_1+k+1]=wz_Ey*Ey[i][j][k_1-k-1];
	if(wz_Ey==-1) for(i=0;i<=xs;i++) for(j=0;j<ys;j++) Ey[i][j][k_1]=0;
}
void symmetric_bacfro_Ez()
{
	for(i=0;i<=xs;i++) for(j=0;j<=ys;j++) for(k=0;k<=k_0-1-1;k++) Ez[i][j][k_1+k]=wz_Ez*Ez[i][j][k_1-k-1];
}

void symmetric_bacfro_Hx()
{
	for(i=0;i<xs-1;i++) for(j=0;j<ys;j++) for(k=0;k<k_0;k++) Hx[i][j][k_1+k]=wz_Hx*Hx[i][j][k_1-k-1];
}
void symmetric_bacfro_Hy()
{
	for(i=0;i<xs;i++) for(j=0;j<ys-1;j++) for(k=0;k<k_0;k++) Hy[i][j][k_1+k]=wz_Hy*Hy[i][j][k_1-k-1];

}
void symmetric_bacfro_Hz()
{
	for(i=0;i<xs;i++) for(j=0;j<ys;j++) for(k=0;k<k_0-1;k++) Hz[i][j][k_1+k]=wz_Hz*Hz[i][j][k_1-k-2];
	if(wz_Hz==-1) for(i=0;i<xs;i++) for(j=0;j<ys;j++) Hz[i][j][k_1-1]=0;
}

void symmetric_lefrig_Ex()
{
    	for(i=0;i<i_0;i++) for(j=0;j<=ys;j++) for(k=0;k<=zs;k++) Ex[i_1+i][j][k]=wx_Ex*Ex[i_1-i-1][j][k];
}
void symmetric_lefrig_Ey()
{
	for(i=0;i<i_0;i++) for(j=0;j<ys;j++)  for(k=0;k<=zs;k++) Ey[i_1+i+1][j][k]=wx_Ey*Ey[i_1-i-1][j][k];
	if(wx_Ey==-1) for(j=0;j<ys;j++)  for(k=0;k<=zs;k++) Ey[i_1][j][k]=0;
}
void symmetric_lefrig_Ez()
{
	for(i=0;i<i_0;i++) for(j=0;j<=ys;j++) for(k=0;k<zs;k++)  Ez[i_1+i+1][j][k]=wx_Ez*Ez[i_1-i-1][j][k];
}

void symmetric_lefrig_Hx()
{
	for(i=0;i<i_0-1;i++) for(j=0;j<ys;j++) for(k=0;k<zs;k++) Hx[i_1+i][j][k]=wx_Hx*Hx[i_1-i-2][j][k];
	if(wx_Ex==-1) for(j=0;j<ys;j++) for(k=0;k<zs;k++) Hx[i_1-1][j][k]=0;
}
void symmetric_lefrig_Hy()
{
	for(i=0;i<i_0;i++) for(j=0;j<ys-1;j++) for(k=0;k<zs;k++) Hy[i_1+i][j][k]=wx_Hy*Hy[i_1-i-1][j][k];
}
void symmetric_lefrig_Hz()
{
	for(i=0;i<i_0;i++) for(j=0;j<ys;j++) for(k=0;k<zs-1;k++) Hz[i_1+i][j][k]=wx_Hz*Hz[i_1-i-1][j][k];
}

void symmetric_bottop_Ex()
{
	for(i=0;i<xs;i++) for(j=0;j<j_0;j++) for(k=0;k<=zs;k++)  Ex[i][j_1+j+1][k]=wy_Ex*Ex[i][j_1-j-1][k];
	if(wy_Ex==-1) for(i=0;i<xs;i++) for(k=0;k<=zs;k++) Ex[i][j_1][k]=0;
}
void symmetric_bottop_Ey()
{
	for(i=0;i<=xs;i++) for(j=0;j<j_0;j++) for(k=0;k<=zs;k++) Ey[i][j_1+j][k]=wy_Ey*Ey[i][j_1-j-1][k];
}
void symmetric_bottop_Ez()
{
	for(i=0;i<=xs;i++) for(j=0;j<j_0;j++) for(k=0;k<zs;k++)  Ez[i][j_1+j+1][k]=wy_Ez*Ez[i][j_1-j-1][k];
	if(wy_Ez==-1) for(i=0;i<=xs;i++) for(k=0;k<zs;k++)  Ez[i][j_1][k]=0;
}

void symmetric_bottop_Hx()
{
	for(i=0;i<xs-1;i++) for(j=0;j<j_0;j++) for(k=0;k<zs;k++) Hx[i][j_1+j][k]=wy_Hx*Hx[i][j_1-j-1][k];
}
void symmetric_bottop_Hy()
{
	for(i=0;i<xs;i++) for(j=0;j<j_0-1;j++) for(k=0;k<zs;k++) Hy[i][j_1+j][k]=wy_Hy*Hy[i][j_1-j-2][k];
	if(wy_Hy==-1) for(i=0;i<xs;i++) for(k=0;k<zs;k++) Hy[i][j_1-1][k]=0;
}
void symmetric_bottop_Hz()
{
	for(i=0;i<xs;i++) for(j=0;j<j_0;j++) for(k=0;k<zs-1;k++) Hz[i][j_1+j][k]=wy_Hz*Hz[i][j_1-j-1][k];
}

inline void set_symmetric_lefrig() 
{
       i_1=boundx+bndx+disPML+xr;
       if (is_calc_ex) symmetric_lefrig_Ex();if (is_calc_ey) symmetric_lefrig_Ey();if (is_calc_ez) symmetric_lefrig_Ez();
       if (is_calc_hx) symmetric_lefrig_Hx();if (is_calc_hy) symmetric_lefrig_Hy();if (is_calc_hz) symmetric_lefrig_Hz();
}
inline void set_symmetric_bottop() 
{
       j_1=boundy+bndy+disPML+yr;
       if (is_calc_ex) symmetric_bottop_Ex();if (is_calc_ey) symmetric_bottop_Ey();if (is_calc_ez) symmetric_bottop_Ez();
       if (is_calc_hx) symmetric_bottop_Hx();if (is_calc_hy) symmetric_bottop_Hy();if (is_calc_hz) symmetric_bottop_Hz();
}
inline void set_symmetric_bacfro() 
{
       k_1=boundz+bndz+disPML+zr;
       if (is_calc_ex) symmetric_bacfro_Ex();if (is_calc_ey) symmetric_bacfro_Ey();if (is_calc_ez) symmetric_bacfro_Ez();
       if (is_calc_hx) symmetric_bacfro_Hx();if (is_calc_hy) symmetric_bacfro_Hy();if (is_calc_hz) symmetric_bacfro_Hz();
}

//////////////////////////////////
// periodyczne warunki brzegowe
#ifdef ACCURATE_BLOCH
#ifdef COMPLEX_FIELD
void bloch_bacfro_Ex()
{
     for(i=0;i<xs;i++) for(j=0;j<=ys;j++) for(k=0;k<k_0;k++)  
	{Ex[i][j][k_1+k+1]=std::polar(abs(Ex[i][j][k_0+k+1]),(typ_prec)(arg(Ex[i][j][k_0+k+1])+bloch_kz));
     Ex[i][j][k_0-k-1]=std::polar(abs(Ex[i][j][k_1-k-1]),(typ_prec)(arg(Ex[i][j][k_1-k-1])-bloch_kz));}
}
void bloch_bacfro_Ey()
{
	for(i=0;i<=xs;i++) for(j=0;j<ys;j++) for(k=0;k<k_0;k++)  
	{Ey[i][j][k_1+k+1]=std::polar(abs(Ey[i][j][k_0+k+1]),(typ_prec)(arg(Ey[i][j][k_0+k+1])+bloch_kz));
    Ey[i][j][k_0-k-1]=std::polar(abs(Ey[i][j][k_1-k-1]),(typ_prec)(arg(Ey[i][j][k_1-k-1])-bloch_kz));}
}	
void bloch_bacfro_Ez()
{   
    for(i=0;i<=xs;i++) for(j=0;j<=ys;j++) for(k=0;k<k_0;k++) 
	{Ez[i][j][k_1+k]=std::polar(abs(Ez[i][j][k_0+k]),(typ_prec)(arg(Ez[i][j][k_0+k])+bloch_kz));
    Ez[i][j][k_0-k-1]=std::polar(abs(Ez[i][j][k_1-k-1]),(typ_prec)(arg(Ez[i][j][k_1-k-1])-bloch_kz));}
}

void bloch_bacfro_Hx()
{
	for(i=0;i<xs-1;i++) for(j=0;j<ys;j++) for(k=0;k<k_0;k++) 
	{Hx[i][j][k_1+k]=std::polar(abs(Hx[i][j][k_0+k]),(typ_prec)(arg(Hx[i][j][k_0+k])+bloch_kz));
    Hx[i][j][k_0-k-1]=std::polar(abs(Hx[i][j][k_1-k-1]),(typ_prec)(arg(Hx[i][j][k_1-k-1])-bloch_kz));}
}
void bloch_bacfro_Hy()
{
	for(i=0;i<xs;i++) for(j=0;j<ys-1;j++) for(k=0;k<k_0;k++) 
	{Hy[i][j][k_1+k]=std::polar(abs(Hy[i][j][k_0+k]),(typ_prec)(arg(Hy[i][j][k_0+k])+bloch_kz));    
    Hy[i][j][k_0-k-1]=std::polar(abs(Hy[i][j][k_1-k-1]),(typ_prec)(arg(Hy[i][j][k_1-k-1])-bloch_kz));}
}
void bloch_bacfro_Hz()
{
	for(i=0;i<xs;i++) for(j=0;j<ys;j++) for(k=0;k<k_0-1;k++) 
	{Hz[i][j][k_1+k]=std::polar(abs(Hz[i][j][k_0+k]),(typ_prec)(arg(Hz[i][j][k_0+k])+bloch_kz));
    Hz[i][j][k_0-k-2]=std::polar(abs(Hz[i][j][k_1-k-2]),(typ_prec)(arg(Hz[i][j][k_1-k-2])-bloch_kz));}
}

void bloch_lefrig_Ex()
{
    for(i=0;i<i_0;i++) for(j=0;j<=ys;j++) for(k=0;k<=zs;k++) 
	{Ex[i_1+i][j][k]=std::polar(abs(Ex[i_0+i][j][k]),(typ_prec)(arg(Ex[i_0+i][j][k])+bloch_kx));
    Ex[i_0-i-1][j][k]=std::polar(abs(Ex[i_1-i-1][j][k]),(typ_prec)(arg(Ex[i_1-i-1][j][k])-bloch_kx));}
}
void bloch_lefrig_Ey()
{
    for(i=0;i<i_0;i++) for(j=0;j<ys;j++)  for(k=0;k<=zs;k++) 
	{Ey[i_1+i+1][j][k]=std::polar(abs(Ey[i_0+i+1][j][k]),(typ_prec)(arg(Ey[i_0+i+1][j][k])+bloch_kx));
    Ey[i_0-i-1][j][k]=std::polar(abs(Ey[i_1-i-1][j][k]),(typ_prec)(arg(Ey[i_1-i-1][j][k])-bloch_kx));}
}
void bloch_lefrig_Ez()
{
    for(i=0;i<i_0;i++) for(j=0;j<=ys;j++) for(k=0;k<zs;k++)
	{Ez[i_1+i+1][j][k]=std::polar(abs(Ez[i_0+i+1][j][k]),(typ_prec)(arg(Ez[i_0+i+1][j][k])+bloch_kx));
    Ez[i_0-i-1][j][k]=std::polar(abs(Ez[i_1-i-1][j][k]),(typ_prec)(arg(Ez[i_1-i-1][j][k])-bloch_kx));}
}

void bloch_lefrig_Hx()
{
    for(i=0;i<i_0-1;i++) for(j=0;j<ys;j++) for(k=0;k<zs;k++) 
	{Hx[i_1+i][j][k]=std::polar(abs(Hx[i_0+i-1][j][k]),(typ_prec)(arg(Hx[i_0+i][j][k])+bloch_kx));
    Hx[i_0-i-2][j][k]=std::polar(abs(Hx[i_1-i-2][j][k]),(typ_prec)(arg(Hx[i_1-i-2][j][k])-bloch_kx));}
}
void bloch_lefrig_Hy()
{
	for(i=0;i<i_0;i++) for(j=0;j<ys-1;j++) for(k=0;k<zs;k++)
	{Hy[i_1+i][j][k]=std::polar(abs(Hy[i_0+i][j][k]),(typ_prec)(arg(Hy[i_0+i][j][k])+bloch_kx));
    Hy[i_0-i-1][j][k]=std::polar(abs(Hy[i_1-i-1][j][k]),(typ_prec)(arg(Hy[i_1-i-1][j][k])-bloch_kx));}
}
void bloch_lefrig_Hz()
{
	for(i=0;i<i_0;i++) for(j=0;j<ys;j++) for(k=0;k<zs-1;k++) 
	{Hz[i_1+i][j][k]=std::polar(abs(Hz[i_0+i][j][k]),(typ_prec)(arg(Hz[i_0+i][j][k])+bloch_kx));
    Hz[i_0-i-1][j][k]=std::polar(abs(Hz[i_1-i-1][j][k]),(typ_prec)(arg(Hz[i_1-i-1][j][k])-bloch_kx));}
}

void bloch_bottop_Ex()
{
    for(i=0;i<xs;i++) for(j=0;j<j_0;j++) for(k=0;k<=zs;k++)  
	{Ex[i][j_1+j+1][k]=std::polar(abs(Ex[i][j_0+j+1][k]),(typ_prec)(arg(Ex[i][j_0+j+1][k])+bloch_ky));
    Ex[i][j_0-j-1][k]=std::polar(abs(Ex[i][j_1-j-1][k]),(typ_prec)(arg(Ex[i][j_1-j-1][k])-bloch_ky));}
}
void bloch_bottop_Ey()
{   
    for(i=0;i<=xs;i++) for(j=0;j<j_0;j++) for(k=0;k<=zs;k++) 
	{Ey[i][j_1+j][k]=std::polar(abs(Ey[i][j_0+j][k]),(typ_prec)(arg(Ey[i][j_0+j][k])+bloch_ky));
    Ey[i][j_0-j-1][k]=std::polar(abs(Ey[i][j_1-j-1][k]),(typ_prec)(arg(Ey[i][j_1-j-1][k])-bloch_ky));}
}
void bloch_bottop_Ez()
{
	for(i=0;i<=xs;i++) for(j=0;j<j_0;j++) for(k=0;k<zs;k++)
	{Ez[i][j_1+j+1][k]=std::polar(abs(Ez[i][j_0+j+1][k]),(typ_prec)(arg(Ez[i][j_0+j+1][k])+bloch_ky));
    Ez[i][j_0-j-1][k]=std::polar(abs(Ez[i][j_1-j-1][k]),(typ_prec)(arg(Ez[i][j_1-j-1][k])-bloch_ky));}
}

void bloch_bottop_Hx()
{
	for(i=0;i<xs-1;i++) for(j=0;j<j_0;j++) for(k=0;k<zs;k++) 
	{Hx[i][j_1+j][k]=std::polar(abs(Hx[i][j_0+j][k]),(typ_prec)(arg(Hx[i][j_0+j][k])+bloch_ky));
    Hx[i][j_0-j-1][k]=std::polar(abs(Hx[i][j_1-j-1][k]),(typ_prec)(arg(Hx[i][j_1-j-1][k])-bloch_ky));}
}
void bloch_bottop_Hy()
{
	for(i=0;i<xs;i++) for(j=0;j<j_0-1;j++) for(k=0;k<zs;k++) 
	{Hy[i][j_1+j][k]=std::polar(abs(Hy[i][j_0+j][k]),(typ_prec)(arg(Hy[i][j_0+j][k])+bloch_ky));
    Hy[i][j_0-j-2][k]=std::polar(abs(Hy[i][j_1-j-2][k]),(typ_prec)(arg(Hy[i][j_1-j-2][k])-bloch_ky));}
}
void bloch_bottop_Hz()
{
	for(i=0;i<xs;i++) for(j=0;j<j_0;j++) for(k=0;k<zs-1;k++) 
	{Hz[i][j_1+j][k]=std::polar(abs(Hz[i][j_0+j][k]),(typ_prec)(arg(Hz[i][j_0+j][k])+bloch_ky));
    Hz[i][j_0-j-1][k]=std::polar(abs(Hz[i][j_1-j-1][k]),(typ_prec)(arg(Hz[i][j_1-j-1][k])-bloch_ky));}
}
#endif
#else
void bloch_bacfro_Ex()
{

    for(i=0;i<xs;i++) for(j=0;j<=ys;j++) for(k=0;k<k_0;k++)  
	{Ex[i][j][k_1+k+1]=bloch_expikz*Ex[i][j][k_0+k+1];Ex[i][j][k_0-k-1]=(bloch_expmikz)*Ex[i][j][k_1-k-1];}
}
void bloch_bacfro_Ey()
{
	for(i=0;i<=xs;i++) for(j=0;j<ys;j++) for(k=0;k<k_0;k++)  
	{Ey[i][j][k_1+k+1]=bloch_expikz*Ey[i][j][k_0+k+1];Ey[i][j][k_0-k-1]=(bloch_expmikz)*Ey[i][j][k_1-k-1];}
}	
void bloch_bacfro_Ez()
{   
    for(i=0;i<=xs;i++) for(j=0;j<=ys;j++) for(k=0;k<k_0;k++) 
	{Ez[i][j][k_1+k]=bloch_expikz*Ez[i][j][k_0+k];Ez[i][j][k_0-k-1]=(bloch_expmikz)*Ez[i][j][k_1-k-1];}
}

void bloch_bacfro_Hx()
{
	for(i=0;i<xs-1;i++) for(j=0;j<ys;j++) for(k=0;k<k_0;k++) 
	{Hx[i][j][k_1+k]=bloch_expikz*Hx[i][j][k_0+k];Hx[i][j][k_0-k-1]=(bloch_expmikz)*Hx[i][j][k_1-k-1];}
}
void bloch_bacfro_Hy()
{
	for(i=0;i<xs;i++) for(j=0;j<ys-1;j++) for(k=0;k<k_0;k++) 
	{Hy[i][j][k_1+k]=bloch_expikz*Hy[i][j][k_0+k];Hy[i][j][k_0-k-1]=(bloch_expmikz)*Hy[i][j][k_1-k-1];}
}
void bloch_bacfro_Hz()
{
	for(i=0;i<xs;i++) for(j=0;j<ys;j++) for(k=0;k<k_0-1;k++) 
	{Hz[i][j][k_1+k]=bloch_expikz*Hz[i][j][k_0+k];Hz[i][j][k_0-k-2]=(bloch_expmikz)*Hz[i][j][k_1-k-2];}
}

void bloch_lefrig_Ex()
{
    for(i=0;i<i_0;i++) for(j=0;j<=ys;j++) for(k=0;k<=zs;k++) 
	{Ex[i_1+i][j][k]=bloch_expikx*Ex[i_0+i][j][k];Ex[i_0-i-1][j][k]=(bloch_expmikx)*Ex[i_1-i-1][j][k];}
}
void bloch_lefrig_Ey()
{
    for(i=0;i<i_0;i++) for(j=0;j<ys;j++)  for(k=0;k<=zs;k++) 
	{Ey[i_1+i+1][j][k]=bloch_expikx*Ey[i_0+i+1][j][k];Ey[i_0-i-1][j][k]=(bloch_expmikx)*Ey[i_1-i-1][j][k];}
}
void bloch_lefrig_Ez()
{ 
    for(i=0;i<i_0;i++) for(j=0;j<=ys;j++) for(k=0;k<zs;k++)
	{Ez[i_1+i+1][j][k]=bloch_expikx*Ez[i_0+i+1][j][k];Ez[i_0-i-1][j][k]=(bloch_expmikx)*Ez[i_1-i-1][j][k];}
}

void bloch_lefrig_Hx()
{
    for(i=0;i<i_0-1;i++) for(j=0;j<ys;j++) for(k=0;k<zs;k++) 
	{Hx[i_1+i][j][k]=bloch_expikx*Hx[i_0+i][j][k];Hx[i_0-i-2][j][k]=(bloch_expmikx)*Hx[i_1-i-2][j][k];}
}
void bloch_lefrig_Hy()
{
	for(i=0;i<i_0;i++) for(j=0;j<ys-1;j++) for(k=0;k<zs;k++)
	{Hy[i_1+i][j][k]=bloch_expikx*Hy[i_0+i][j][k];Hy[i_0-i-1][j][k]=(bloch_expmikx)*Hy[i_1-i-1][j][k];}
}
void bloch_lefrig_Hz()
{
	for(i=0;i<i_0;i++) for(j=0;j<ys;j++) for(k=0;k<zs-1;k++) 
	{Hz[i_1+i][j][k]=bloch_expikx*Hz[i_0+i][j][k];Hz[i_0-i-1][j][k]=(bloch_expmikx)*Hz[i_1-i-1][j][k];}
}

void bloch_bottop_Ex()
{
    for(i=0;i<xs;i++) for(j=0;j<j_0;j++) for(k=0;k<=zs;k++)  
	{Ex[i][j_1+j+1][k]=bloch_expiky*Ex[i][j_0+j+1][k];Ex[i][j_0-j-1][k]=(bloch_expmiky)*Ex[i][j_1-j-1][k];}
}
void bloch_bottop_Ey()
{   
    for(i=0;i<=xs;i++) for(j=0;j<j_0;j++) for(k=0;k<=zs;k++) 
	{Ey[i][j_1+j][k]=bloch_expiky*Ey[i][j_0+j][k];Ey[i][j_0-j-1][k]=(bloch_expmiky)*Ey[i][j_1-j-1][k];}
}
void bloch_bottop_Ez()
{
	for(i=0;i<=xs;i++) for(j=0;j<j_0;j++) for(k=0;k<zs;k++)
	{Ez[i][j_1+j+1][k]=bloch_expiky*Ez[i][j_0+j+1][k];Ez[i][j_0-j-1][k]=(bloch_expmiky)*Ez[i][j_1-j-1][k];}
}

void bloch_bottop_Hx()
{
	for(i=0;i<xs-1;i++) for(j=0;j<j_0;j++) for(k=0;k<zs;k++) 
	{Hx[i][j_1+j][k]=bloch_expiky*Hx[i][j_0+j][k];Hx[i][j_0-j-1][k]=(bloch_expmiky)*Hx[i][j_1-j-1][k];}
}
void bloch_bottop_Hy()
{
	for(i=0;i<xs;i++) for(j=0;j<j_0-1;j++) for(k=0;k<zs;k++) 
	{Hy[i][j_1+j][k]=bloch_expiky*Hy[i][j_0+j][k];Hy[i][j_0-j-2][k]=(bloch_expmiky)*Hy[i][j_1-j-2][k];}
}
void bloch_bottop_Hz()
{
	for(i=0;i<xs;i++) for(j=0;j<j_0;j++) for(k=0;k<zs-1;k++) 
	{Hz[i][j_1+j][k]=bloch_expiky*Hz[i][j_0+j][k];Hz[i][j_0-j-1][k]=(bloch_expmiky)*Hz[i][j_1-j-1][k];}
}
#endif


inline void set_bloch_lefrig() 
{
       i_0=boundx+bndx+disPML;i_1=i_0+xr;
       if (is_calc_ex) bloch_lefrig_Ex();if (is_calc_ey) bloch_lefrig_Ey();if (is_calc_ez) bloch_lefrig_Ez();
       if (is_calc_hx) bloch_lefrig_Hx();if (is_calc_hy) bloch_lefrig_Hy();if (is_calc_hz) bloch_lefrig_Hz();
}
inline void set_bloch_bottop() 
{
       j_0=boundy+bndy+disPML;j_1=j_0+yr;
       if (is_calc_ex) bloch_bottop_Ex();if (is_calc_ey) bloch_bottop_Ey();if (is_calc_ez) bloch_bottop_Ez();
       if (is_calc_hx) bloch_bottop_Hx();if (is_calc_hy) bloch_bottop_Hy();if (is_calc_hz) bloch_bottop_Hz();
}
inline void set_bloch_bacfro() 
{
       k_0=boundz+bndz+disPML;k_1=k_0+zr;
       if (is_calc_ex) bloch_bacfro_Ex();if (is_calc_ey) bloch_bacfro_Ey();if (is_calc_ez) bloch_bacfro_Ez();
       if (is_calc_hx) bloch_bacfro_Hx();if (is_calc_hy) bloch_bacfro_Hy();if (is_calc_hz) bloch_bacfro_Hz();
}
//////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                              //
// katy teta,fi,psi jak w Taflovie s.215 - rys. 5.14a, s.221-224                                //
//  teta,fi - okreslaja kierunek wektora falowego w uk;adzie                                    //
//  sferycznym                                                                                  //
//  psi - kat miedzy wektorem k x z i wektorem E pola                                           //
//                                                                                              //
// wektor falowy                                                                                //
//  kx=k*cos(fi)*sin(teta);                                                                     //
//  ky=k*sin(fi)*sin(teta);                                                                     //
//  kz=k*cos(teta);                                                                             //
//                                                                                              //
// Skladowe pola  (za Taflove)                                                                  //
//                                                                                              //
//  Hx=H*(sin(psi)*sin(fi)+cos(psi)*cos(teta)*cos(fi));                                         //
//  Hy=H*(-sin(psi)*cos(fi)+cos(psi)*cos(teta)*sin(fi));                                        //
//  Hz=H*(-cos(psi)*sin(teta));                                                                 //
//                                                                                              //
//  Ex=E*(cos(psi)*sin(fi)-sin(psi)*cos(teta)*cos(fi));                                         //
//  Ey=E*(-cos(psi)*cos(fi)-sin(psi)*cos(teta)*sin(fi));                                        //
//  Ez=E*(sin(psi)*sin(teta);                                                                   //
//                                                                                              //
// gdzie pola En i E w odleglosi d od zrodla ( frontu falowego w chwili t=0)                    //
//                                                                                              //
//  d=k*r                                                                                       //
//  r zalezy od teta i fi ( Taflove, s.222)                                                     //
//////////////////////////////////////////////////////////////////////////////////////////////////


// tworzenie tablic wykorzystywanych przez program

bool init_structure_arrays()
{
	if (mode3D)
	{
		if (!init_array3D(mat3,2*xr+1,2*yr+1,2*zr+1)) return false;
    	
        if (is_calc_ex) if (!init_array3D(matEx3,xr+2*bndx+2*disPML,yr+1+2*bndy+2*disPML,zr+2*bndz+1+2*disPML)) return false;
        if (is_calc_ey) if (!init_array3D(matEy3,xr+1+2*bndx+2*disPML,yr+2*bndy+2*disPML,zr+2*bndz+1+2*disPML)) return false;
        if (is_calc_ez) if (!init_array3D(matEz3,xr+1+2*bndx+2*disPML,yr+1+2*bndy+2*disPML,zr+2*bndz+2*disPML)) return false;
    	
        if (is_calc_hx) if (!init_array3D(matHx3,xr+2*bndx-1+2*disPML,yr+2*bndy+2*disPML,zr+2*bndz+2*disPML)) return false;
        if (is_calc_hy) if (!init_array3D(matHy3,xr+2*bndx+2*disPML,yr+2*bndy-1+2*disPML,zr+2*bndz+2*disPML)) return false;
        if (is_calc_hz) if (!init_array3D(matHz3,xr+2*bndx+2*disPML,yr+2*bndy+2*disPML,zr+2*bndz-1+2*disPML)) return false;
	}
	else
	{
    	if (!init_array2D(mat,2*xr+1,2*yr+1)) return false;
	    
        if (is_calc_ex) if (!init_array2D(matEx,xr+2*bndx+2*disPML,yr+1+2*bndy+2*disPML)) return false;
        if (is_calc_ey) if (!init_array2D(matEy,xr+1+2*bndx+2*disPML,yr+2*bndy+2*disPML)) return false;
        if (is_calc_ez) if (!init_array2D(matEz,xr+1+2*bndx+2*disPML,yr+1+2*bndy+2*disPML)) return false;
	    
        if (is_calc_hx) if (!init_array2D(matHx,xr+2*bndx-1+2*disPML,yr+2*bndy+2*disPML)) return false;
        if (is_calc_hy) if (!init_array2D(matHy,xr+2*bndx+2*disPML,yr+2*bndy-1+2*disPML)) return false;
        if (is_calc_hz) if (!init_array2D(matHz,xr+2*bndx+2*disPML,yr+2*bndy+2*disPML)) return false;
	}
	return true;
}

bool init_arrays()
{

    if (is_calc_ex) if (!init_array3D(Ex,xs,ys+1,zs+1)) return false;
    if (is_calc_ey) if (!init_array3D(Ey,xs+1,ys,zs+1)) return false;
    if (is_calc_ez) if (!init_array3D(Ez,xs+1,ys+1,zs)) return false;
    
    if (is_calc_hx) if (!init_array3D(Hx,xs-1,ys,zs)) return false;
    if (is_calc_hy) if (!init_array3D(Hy,xs,ys-1,zs)) return false;
    if (is_calc_hz) if (!init_array3D(Hz,xs,ys,zs-1)) return false;

if (wy!=1)
{
    if (is_calc_ex) if (!init_array3D(bot_Dx,xs,boundy,zr+1+2*bndz+2*disPML)) return false;
    
    if (is_calc_ey) if (!init_array3D(bot_Dy,xs+1,boundy,zr+1+2*bndz+2*disPML)) return false;
	
	if (is_calc_ez)
	{
     if (!init_array3D(bot_Dz_lef,boundx,boundy,zr+2*bndz+2*disPML)) return false;
     if (!init_array3D(bot_Dz_rig,boundx,boundy,zr+2*bndz+2*disPML)) return false;
    }
    
    if (is_calc_hx) if (!init_array3D(bot_Bx,xs-1,boundy,zr+2*bndz+2*disPML)) return false;
    if (is_calc_hy) if (!init_array3D(bot_By,xs,boundy,zr+2*bndz+2*disPML)) return false;
	
	if (is_calc_hz)
	{
     if (!init_array3D(bot_Bz_lef,boundx,boundy,zr-1+2*bndz+2*disPML)) return false;
     if (!init_array3D(bot_Bz_rig,boundx,boundy,zr-1+2*bndz+2*disPML)) return false;
    }
}

if (wy==0)
{
    if (is_calc_ex) if (!init_array3D(top_Dx,xs,boundy,zr+1+2*bndz+2*disPML)) return false;
    if (is_calc_ey) if (!init_array3D(top_Dy,xs+1,boundy,zr+1+2*bndz+2*disPML)) return false;
	
	if (is_calc_ez)
	{
     if (!init_array3D(top_Dz_lef,boundx,boundy,zr+2*bndz+2*disPML)) return false;
     if (!init_array3D(top_Dz_rig,boundx,boundy,zr+2*bndz+2*disPML)) return false;
    }
     
    if (is_calc_hx) if (!init_array3D(top_Bx,xs-1,boundy,zr+2*bndz+2*disPML)) return false;    
    if (is_calc_hy) if (!init_array3D(top_By,xs,boundy,zr+2*bndz+2*disPML)) return false;
    
    if (is_calc_hz)
    {
	 if (!init_array3D(top_Bz_lef,boundx,boundy,zr-1+2*bndz+2*disPML)) return false;
     if (!init_array3D(top_Bz_rig,boundx,boundy,zr-1+2*bndz+2*disPML)) return false;
    }
}

if (wx!=1)
{
    if (is_calc_ex) if (!init_array3D(lef_Dx,boundx,yr+1+2*bndy+2*disPML,zr+1+2*bndz+2*disPML)) return false;
    if (is_calc_ez) if (!init_array3D(lef_Dz,boundx,yr+1+2*bndy+2*disPML,zr+2*bndz+2*disPML)) return false;
    if (is_calc_hx) if (!init_array3D(lef_Bx,boundx,yr+2*bndy+2*disPML,zr+2*bndz+2*disPML)) return false;
    if (is_calc_hz) if (!init_array3D(lef_Bz,boundx,yr+2*bndy+2*disPML,zr-1+2*bndz+2*disPML)) return false;
}

if (wx==0)
{
    if (is_calc_ex) if (!init_array3D(rig_Dx,boundx,yr+1+2*bndy+2*disPML,zr+1+2*bndz+2*disPML)) return false;
    if (is_calc_ez) if (!init_array3D(rig_Dz,boundx,yr+1+2*bndy+2*disPML,zr+2*bndz+2*disPML)) return false;
    if (is_calc_hx) if (!init_array3D(rig_Bx,boundx,yr+2*bndy+2*disPML,zr+2*bndz+2*disPML)) return false; 
    if (is_calc_hz) if (!init_array3D(rig_Bz,boundx,yr+2*bndy+2*disPML,zr-1+2*bndz+2*disPML)) return false;
}

if (wz==0)
{
    if (is_calc_ex)
    {
      if (!init_array3D(fro_Dx_bot,xs,boundy,boundz)) return false;
	  if (!init_array3D(fro_Dx_top,xs,boundy,boundz)) return false;

	  if (!init_array3D(fro_Dx_lef,boundx,yr+1+2*bndy+2*disPML,boundz)) return false;
	  if (!init_array3D(fro_Dx_rig,boundx,yr+1+2*bndy+2*disPML,boundz)) return false;
     }
    if (is_calc_ey) if (!init_array3D(fro_Dy,xs+1,ys,boundz)) return false;   
    if (is_calc_ez) if (!init_array3D(fro_Dz,xs+1,ys+1,boundz)) return false;

    if (is_calc_hx)
    {
	 if (!init_array3D(fro_Bx_bot,xs-1,boundy,boundz)) return false;
	 if (!init_array3D(fro_Bx_top,xs-1,boundy,boundz)) return false;
     if (!init_array3D(fro_Bx_lef,boundx,yr+2*bndy+2*disPML,boundz)) return false;
     if (!init_array3D(fro_Bx_rig,boundx,yr+2*bndy+2*disPML,boundz)) return false;
    }
	if (is_calc_hy)	if (!init_array3D(fro_By,xs,ys-1,boundz)) return false;
	if (is_calc_hz) if (!init_array3D(fro_Bz,xs,ys,boundz)) return false;
}

if (wz!=1)
{
    if (is_calc_ex)
    {
     if (!init_array3D(bac_Dx_bot,xs,boundy,boundz)) return false;
	 if (!init_array3D(bac_Dx_top,xs,boundy,boundz)) return false;
	 if (!init_array3D(bac_Dx_lef,boundx,yr+1+2*bndy+2*disPML,boundz)) return false;
	 if (!init_array3D(bac_Dx_rig,boundx,yr+1+2*bndy+2*disPML,boundz)) return false;
    }
    
    if (is_calc_ey) if (!init_array3D(bac_Dy,xs+1,ys,boundz)) return false;
   	if (is_calc_ez) if (!init_array3D(bac_Dz,xs+1,ys+1,boundz)) return false;

	if (is_calc_hx)
	{
	 if (!init_array3D(bac_Bx_bot,xs-1,boundy,boundz)) return false;
	 if (!init_array3D(bac_Bx_top,xs-1,boundy,boundz)) return false;
	 if (!init_array3D(bac_Bx_lef,boundx,yr+2*bndy+2*disPML,boundz)) return false;
     if (!init_array3D(bac_Bx_rig,boundx,yr+2*bndy+2*disPML,boundz)) return false;
    }
    
	if (is_calc_hy) if (!init_array3D(bac_By,xs,ys-1,boundz)) return false;
	if (is_calc_hz) if (!init_array3D(bac_Bz,xs,ys,boundz)) return false;
}

	if (!init_array1D(c1x_E,bound+1)) return false;if (!init_array1D(c2x_E,bound+1)) return false;if (!init_array1D(c3x_E,bound+1)) return false;
    if (!init_array1D(c4x_E,bound+1)) return false;if (!init_array1D(c5x_E,bound+1)) return false;if (!init_array1D(c6x_E,bound+1)) return false;

    if (!init_array1D(c1x_H,bound+1)) return false;if (!init_array1D(c2x_H,bound+1)) return false;if (!init_array1D(c3x_H,bound+1)) return false;
    if (!init_array1D(c4x_H,bound+1)) return false;if (!init_array1D(c5x_H,bound+1)) return false;if (!init_array1D(c6x_H,bound+1)) return false;

	if (!init_array1D(c1y_E,bound+1)) return false;if (!init_array1D(c2y_E,bound+1)) return false;if (!init_array1D(c3y_E,bound+1)) return false;
    if (!init_array1D(c4y_E,bound+1)) return false;if (!init_array1D(c5y_E,bound+1)) return false;if (!init_array1D(c6y_E,bound+1)) return false;

    if (!init_array1D(c1y_H,bound+1)) return false;if (!init_array1D(c2y_H,bound+1)) return false;if (!init_array1D(c3y_H,bound+1)) return false;
    if (!init_array1D(c4y_H,bound+1)) return false;if (!init_array1D(c5y_H,bound+1)) return false;if (!init_array1D(c6y_H,bound+1)) return false;


	if (!init_array1D(c1z_E,bound+1)) return false;if (!init_array1D(c2z_E,bound+1)) return false;if (!init_array1D(c3z_E,bound+1)) return false;
    if (!init_array1D(c4z_E,bound+1)) return false;if (!init_array1D(c5z_E,bound+1)) return false;if (!init_array1D(c6z_E,bound+1)) return false;

    if (!init_array1D(c1z_H,bound+1)) return false;if (!init_array1D(c2z_H,bound+1)) return false;if (!init_array1D(c3z_H,bound+1)) return false;
    if (!init_array1D(c4z_H,bound+1)) return false;if (!init_array1D(c5z_H,bound+1)) return false;if (!init_array1D(c6z_H,bound+1)) return false;



    if (!init_array1D(dmit,nmat)) return false;if (!init_array1D(depst,nmat)) return false;
	if (!init_array1D(dsigt,nmat)) return false;if (!init_array1D(dsiht,nmat)) return false;

    for (int s=0;s<n_SRCS;s++)
    {
     if (!init_array1D(i_E[s],pulse_length+2)) return false;if (!init_array1D(i_H[s],pulse_length+1)) return false;
    }
    
    if (!init_array1D(drd_ka_E,nmat)) return false;if (!init_array1D(drd_be_E,nmat)) return false;if (!init_array1D(drd_J_E,nmat)) return false;
    if (!init_array1D(drd_ka_H,nmat)) return false;if (!init_array1D(drd_be_H,nmat)) return false;if (!init_array1D(drd_J_H,nmat)) return false;

    if (!init_array1D(lor_alf_E,nmat)) return false;if (!init_array1D(lor_ksi_E,nmat)) return false;if (!init_array1D(lor_gam_E,nmat)) return false;
	if (!init_array1D(lor_a1_E,nmat)) return false;if (!init_array1D(lor_a2_E,nmat)) return false;if (!init_array1D(lor_a3_E,nmat)) return false;

	if (!init_array1D(lor_alf_H,nmat)) return false;if (!init_array1D(lor_ksi_H,nmat)) return false;if (!init_array1D(lor_gam_H,nmat)) return false;
	if (!init_array1D(lor_a1_H,nmat)) return false;if (!init_array1D(lor_a2_H,nmat)) return false;if (!init_array1D(lor_a3_H,nmat)) return false;

    if (!init_array1D(dby_ka_E,nmat)) return false;if (!init_array1D(dby_be_E,nmat)) return false;if (!init_array1D(dby_J_E,nmat)) return false;
    if (!init_array1D(dby_ka_H,nmat)) return false;if (!init_array1D(dby_be_H,nmat)) return false;if (!init_array1D(dby_J_H,nmat)) return false;

	if (n_drd_Ex>0 && is_calc_ex)
	{
		if (!init_array1D(drd_Jex,n_drd_Ex)) return false;if (!init_array1D(drd_Ex_n_1,n_drd_Ex)) return false;
	}
	if (n_drd_Ey>0 && is_calc_ey)
	{
		if (!init_array1D(drd_Jey,n_drd_Ey)) return false;if (!init_array1D(drd_Ey_n_1,n_drd_Ey)) return false;
	}
	if (n_drd_Ez>0 && is_calc_ez)
	{
		if (!init_array1D(drd_Jez,n_drd_Ez)) return false;if (!init_array1D(drd_Ez_n_1,n_drd_Ez)) return false;
	}
	
	if (n_drd_Hx>0 && is_calc_hx)
	{
		if (!init_array1D(drd_Jhx,n_drd_Hx)) return false;if (!init_array1D(drd_Hx_n_1,n_drd_Hx)) return false;
	}
	if (n_drd_Hy>0 && is_calc_hy)
	{
		if (!init_array1D(drd_Jhy,n_drd_Hy)) return false;if (!init_array1D(drd_Hy_n_1,n_drd_Hy)) return false;
	}
	if (n_drd_Hz>0 && is_calc_hz)
	{
		if (!init_array1D(drd_Jhz,n_drd_Hz)) return false;if (!init_array1D(drd_Hz_n_1,n_drd_Hz)) return false;
	}

	if (n_lor_Ex>0 && is_calc_ex)
	{
		if (!init_array1D(lor_Jex,n_lor_Ex)) return false;if (!init_array1D(lor_Jex_n_1,n_lor_Ex)) return false;
		if (!init_array1D(lor_Ex_n_1,n_lor_Ex)) return false;if (!init_array1D(lor_Ex_n_2,n_lor_Ex)) return false;
	}
	
    if (n_lor_Ey>0 && is_calc_ey)
	{
		if (!init_array1D(lor_Jey,n_lor_Ey)) return false;if (!init_array1D(lor_Jey_n_1,n_lor_Ey)) return false;
		if (!init_array1D(lor_Ey_n_1,n_lor_Ey)) return false;if (!init_array1D(lor_Ey_n_2,n_lor_Ey)) return false;
	}
	
    if (n_lor_Ez>0 && is_calc_ez)
	{
		if (!init_array1D(lor_Jez,n_lor_Ez)) return false;if (!init_array1D(lor_Jez_n_1,n_lor_Ez)) return false;
		if (!init_array1D(lor_Ez_n_1,n_lor_Ez)) return false;if (!init_array1D(lor_Ez_n_2,n_lor_Ez)) return false;
	}
	
	if (n_lor_Hx>0 && is_calc_hx)
	{
		if (!init_array1D(lor_Jhx,n_lor_Hx)) return false;if (!init_array1D(lor_Jhx_n_1,n_lor_Hx)) return false;
		if (!init_array1D(lor_Hx_n_1,n_lor_Hx)) return false;if (!init_array1D(lor_Hx_n_2,n_lor_Hx)) return false;
	}
	
    if (n_lor_Hy>0 && is_calc_hy)
	{
		if (!init_array1D(lor_Jhy,n_lor_Hy)) return false;if (!init_array1D(lor_Jhy_n_1,n_lor_Hy)) return false;
		if (!init_array1D(lor_Hy_n_1,n_lor_Hy)) return false;if (!init_array1D(lor_Hy_n_2,n_lor_Hy)) return false;
	}
	
    if (n_lor_Hz>0 && is_calc_hz)
	{
		if (!init_array1D(lor_Jhz,n_lor_Hz)) return false;if (!init_array1D(lor_Jhz_n_1,n_lor_Hz)) return false;
		if (!init_array1D(lor_Hz_n_1,n_lor_Hz)) return false;if (!init_array1D(lor_Hz_n_2,n_lor_Hz)) return false;
	}	

	if (n_dby_Ex>0 && is_calc_ex)
	{
		if (!init_array1D(dby_Jex,n_dby_Ex)) return false;if (!init_array1D(dby_Ex_n_1,n_dby_Ex)) return false;
    }
	
    if (n_dby_Ey>0 && is_calc_ey)
	{
		if (!init_array1D(dby_Jey,n_dby_Ey)) return false;if (!init_array1D(dby_Ey_n_1,n_dby_Ey)) return false;
    }
   	
    if (n_dby_Ez>0 && is_calc_ez)
	{
		if (!init_array1D(dby_Jez,n_dby_Ez)) return false;if (!init_array1D(dby_Ez_n_1,n_dby_Ez)) return false;
    }
  
	if (n_dby_Hx>0 && is_calc_hx)
	{
		if (!init_array1D(dby_Jhx,n_dby_Hx)) return false;if (!init_array1D(dby_Hx_n_1,n_dby_Hx)) return false;
    }
	
    if (n_dby_Hy>0 && is_calc_hy)
	{
		if (!init_array1D(dby_Jhy,n_dby_Hy)) return false;if (!init_array1D(dby_Hy_n_1,n_dby_Hy)) return false;
    }
   	
    if (n_dby_Hz>0 && is_calc_hz)
	{
		if (!init_array1D(dby_Jhz,n_dby_Hz)) return false;if (!init_array1D(dby_Hz_n_1,n_dby_Hz)) return false;
    }  

    return true;

}

// usuwanie tablic tworzonych przez program

void deinit_arrays()
{
    if (is_calc_ex) del_array3D(Ex);if (is_calc_ey) del_array3D(Ey);if (is_calc_ez) del_array3D(Ez);
    if (is_calc_hx) del_array3D(Hx);if (is_calc_hy) del_array3D(Hy);if (is_calc_hz) del_array3D(Hz);

if (wy!=1 && wy!=3)
{
    if (is_calc_ex) del_array3D(bot_Dx);if (is_calc_ey) del_array3D(bot_Dy);
	if (is_calc_ez) {del_array3D(bot_Dz_lef);del_array3D(bot_Dz_rig);}
    if (is_calc_hx) del_array3D(bot_Bx);if (is_calc_hy) del_array3D(bot_By);
	if (is_calc_hz) {del_array3D(bot_Bz_lef);del_array3D(bot_Bz_rig);}
}
if (wy==0)
{
	if (is_calc_ex) del_array3D(top_Dx);if (is_calc_ey) del_array3D(top_Dy);
	if (is_calc_ez) {del_array3D(top_Dz_lef);del_array3D(top_Dz_rig);}
    if (is_calc_hx) del_array3D(top_Bx);if (is_calc_hy) del_array3D(top_By);
	if (is_calc_hz) {del_array3D(top_Bz_lef);del_array3D(top_Bz_rig);}
}

if (wx!=1 && wx!=3)
{
    if (is_calc_ex) del_array3D(lef_Dx);if (is_calc_ez) del_array3D(lef_Dz);
    if (is_calc_hx) del_array3D(lef_Bx);if (is_calc_hz) del_array3D(lef_Bz);
}

if (wx==0)
{
    if (is_calc_ex) del_array3D(rig_Dx);if (is_calc_ez) del_array3D(rig_Dz);
    if (is_calc_hx) del_array3D(rig_Bx);if (is_calc_hz) del_array3D(rig_Bz);
}

if (wz==0)
{
     if (is_calc_ex)
     {      
      del_array3D(fro_Dx_bot);del_array3D(fro_Dx_top);
	  del_array3D(fro_Dx_lef);del_array3D(fro_Dx_rig);
     }
     
	if (is_calc_ey) del_array3D(fro_Dy);if (is_calc_ez) del_array3D(fro_Dz);

    if (is_calc_hx)
    { 
	 del_array3D(fro_Bx_bot);del_array3D(fro_Bx_top);
	 del_array3D(fro_Bx_lef);del_array3D(fro_Bx_rig);
    }
    
	if (is_calc_hy) del_array3D(fro_By);if (is_calc_hz) del_array3D(fro_Bz);
}

if (wz!=1 && wz!=3)
{
    if (is_calc_ex) 
	{
     del_array3D(bac_Dx_bot);del_array3D(bac_Dx_top);
	 del_array3D(bac_Dx_lef);del_array3D(bac_Dx_rig);
    }
	if (is_calc_ey) del_array3D(bac_Dy);if (is_calc_ez) del_array3D(bac_Dz);

    if (is_calc_hx)
    {
	 del_array3D(bac_Bx_bot);del_array3D(bac_Bx_top);
	 del_array3D(bac_Bx_lef);del_array3D(bac_Bx_rig);
    }
	if (is_calc_hy) del_array3D(bac_By);if (is_calc_hz) del_array3D(bac_Bz);
}

    del_array1D(c1x_E);del_array1D(c2x_E);del_array1D(c3x_E);
    del_array1D(c4x_E);del_array1D(c5x_E);del_array1D(c6x_E);

    del_array1D(c1x_H);del_array1D(c2x_H);del_array1D(c3x_H);
    del_array1D(c4x_H);del_array1D(c5x_H);del_array1D(c6x_H);


    del_array1D(c1y_E);del_array1D(c2y_E);del_array1D(c3y_E);
    del_array1D(c4y_E);del_array1D(c5y_E);del_array1D(c6y_E);

    del_array1D(c1y_H);del_array1D(c2y_H);del_array1D(c3y_H);
    del_array1D(c4y_H);del_array1D(c5y_H);del_array1D(c6y_H);


    del_array1D(c1z_E);del_array1D(c2z_E);del_array1D(c3z_E);
    del_array1D(c4z_E);del_array1D(c5z_E);del_array1D(c6z_E);

    del_array1D(c1z_H);del_array1D(c2z_H);del_array1D(c3z_H);
    del_array1D(c4z_H);del_array1D(c5z_H);del_array1D(c6z_H);


    del_array1D(dmit);del_array1D(depst);del_array1D(dsigt);del_array1D(dsiht);
    for (int s=0;s<n_SRCS;s++)
    {
     del_array1D(i_E[s]);del_array1D(i_H[s]);
     if(time_dependency_length[s]>0) del_array1D(time_dependency[s]);
    }
    del_array1D(drd_ka_E);del_array1D(drd_be_E);del_array1D(drd_J_E);
    del_array1D(drd_ka_H);del_array1D(drd_be_H);del_array1D(drd_J_H);

	del_array1D(lor_a1_E);del_array1D(lor_a2_E);del_array1D(lor_a3_E);
    del_array1D(lor_alf_E);del_array1D(lor_ksi_E);del_array1D(lor_gam_E);
	
	del_array1D(lor_a1_H);del_array1D(lor_a2_H);del_array1D(lor_a3_H);
    del_array1D(lor_alf_H);del_array1D(lor_ksi_H);del_array1D(lor_gam_H);

    del_array1D(dby_ka_E);del_array1D(dby_be_E);del_array1D(dby_J_E);
    del_array1D(dby_ka_H);del_array1D(dby_be_H);del_array1D(dby_J_H);

	if (n_drd_Ex>0 && is_calc_ex) {del_array1D(drd_Jex);del_array1D(drd_Ex_n_1);}
	if (n_drd_Ey>0 && is_calc_ey) {del_array1D(drd_Jey);del_array1D(drd_Ey_n_1);}
	if (n_drd_Ez>0 && is_calc_ez) {del_array1D(drd_Jez);del_array1D(drd_Ez_n_1);}

	if (n_drd_Hx>0 && is_calc_hx) {del_array1D(drd_Jhx);del_array1D(drd_Hx_n_1);}
	if (n_drd_Hy>0 && is_calc_hy) {del_array1D(drd_Jhy);del_array1D(drd_Hy_n_1);}
	if (n_drd_Hz>0 && is_calc_hz) {del_array1D(drd_Jhz);del_array1D(drd_Hz_n_1);} 
        	
	if (n_lor_Ex>0 && is_calc_ex) {del_array1D(lor_Jex);del_array1D(lor_Jex_n_1);del_array1D(lor_Ex_n_1);del_array1D(lor_Ex_n_2);}
	if (n_lor_Ey>0 && is_calc_ey) {del_array1D(lor_Jey);del_array1D(lor_Jey_n_1);del_array1D(lor_Ey_n_1);del_array1D(lor_Ey_n_2);}
	if (n_lor_Ez>0 && is_calc_ez) {del_array1D(lor_Jez);del_array1D(lor_Jez_n_1);del_array1D(lor_Ez_n_1);del_array1D(lor_Ez_n_2);}
 
 	if (n_lor_Hx>0 && is_calc_hx) {del_array1D(lor_Jhx);del_array1D(lor_Jhx_n_1);del_array1D(lor_Hx_n_1);del_array1D(lor_Hx_n_2);}
	if (n_lor_Hy>0 && is_calc_hy) {del_array1D(lor_Jhy);del_array1D(lor_Jhy_n_1);del_array1D(lor_Hy_n_1);del_array1D(lor_Hy_n_2);}
	if (n_lor_Hz>0 && is_calc_hz) {del_array1D(lor_Jhz);del_array1D(lor_Jhz_n_1);del_array1D(lor_Hz_n_1);del_array1D(lor_Hz_n_2);}

	if (n_dby_Ex>0 && is_calc_ex) {del_array1D(dby_Jex);del_array1D(dby_Ex_n_1);}
	if (n_dby_Ey>0 && is_calc_ey) {del_array1D(dby_Jey);del_array1D(dby_Ey_n_1);}
	if (n_dby_Ez>0 && is_calc_ez) {del_array1D(dby_Jez);del_array1D(dby_Ez_n_1);}

	if (n_dby_Hx>0 && is_calc_hx) {del_array1D(dby_Jhx);del_array1D(dby_Hx_n_1);}
	if (n_dby_Hy>0 && is_calc_hy) {del_array1D(dby_Jhy);del_array1D(dby_Hy_n_1);}
	if (n_dby_Hz>0 && is_calc_hz) {del_array1D(dby_Jhz);del_array1D(dby_Hz_n_1);}       

}

void deinit_structure_arrays()
{
	if (mode3D)
	{
        if (is_calc_ex) del_array3D(matEx3);if (is_calc_ey) del_array3D(matEy3);if (is_calc_ez) del_array3D(matEz3);
	    if (is_calc_hx) del_array3D(matHx3);if (is_calc_hy) del_array3D(matHy3);if (is_calc_hz) del_array3D(matHz3);
	}
	else
	{
	    if (is_calc_ex) del_array2D(matEx);if (is_calc_ey) del_array2D(matEy);if (is_calc_ez) del_array2D(matEz);
	    if (is_calc_hx) del_array2D(matHx);if (is_calc_hy) del_array2D(matHy);if (is_calc_hz) del_array2D(matHz);
	}
}

//usuwanie tablic tworzonych podczas czytania pliku z parametrami;
void deinit_other_arrays()
{   
	del_array1D(eps);del_array1D(mi);del_array1D(sig);del_array1D(sih);

	del_array1D(drd_omp_E);del_array1D(drd_gam_E);del_array1D(is_drd_E);
	del_array1D(drd_omp_H);del_array1D(drd_gam_H);del_array1D(is_drd_H);
	
	del_array1D(is_drd_lor_E);del_array1D(is_drd_dby_E);del_array1D(is_lor_dby_E);
	del_array1D(is_drd_lor_H);del_array1D(is_drd_dby_H);del_array1D(is_lor_dby_H);
	
	del_array1D(lor_epd_E);del_array1D(lor_omp_E);del_array1D(lor_del_E);del_array1D(is_lor_E);
	del_array1D(lor_epd_H);del_array1D(lor_omp_H);del_array1D(lor_del_H);del_array1D(is_lor_H);

	del_array1D(dby_deps_E);del_array1D(dby_tau_E);del_array1D(is_dby_E);
	del_array1D(dby_deps_H);del_array1D(dby_tau_H);del_array1D(is_dby_H);

	del_array1D(is_PEC);del_array1D(is_PMC);
	
    if (n_outputs>0)
    {
	 del_array1D(out_min_t);del_array1D(out_max_t);
	 del_array1D(slice);
	 del_array1D(out_min_slice);del_array1D(out_max_slice);
	 del_array1D(di);del_array1D(dj);del_array1D(dk);
     if (skladowe) {delete[] skladowe[0];delete[] skladowe;}
	 del_array1D(is_averaged);del_array1D(is_fourie);del_array1D(ni_fourie);del_array1D(is_temp);
    }
    
	for (int s=0;s<n_SRCS;s++)
	{
	 if (source_type[s]==2) {del_array2D(shape_ex_real[s]);del_array2D(shape_ey_real[s]);del_array2D(shape_ex_imag[s]);del_array2D(shape_ey_imag[s]);}
    }

	if (detector_n>0)
	{
     del_array1D(detector_skladowa);
     del_array1D(detector_x);del_array1D(detector_y);del_array1D(detector_z);
	 del_array1D(detector_t_start);del_array1D(detector_t_end);del_array1D(detector_t_step);
     del_array2D(recorded_at_detector);
    } 
}

////////////////////////////////////////////////////////////////////////////////
// Struktura pliku PDP
//
// plik sluzy do przechowywania wartosci pola skalarnego w 3D,2D i 1D
//
// <int> length_of_name		dlugosc std::stringu z nazwa
// <char[length_of_name]> name	nazwa
//
// <int> x_size			rozmiary pola
// <int> y_size
// <int> z_size
//
// <int> format_spec               deklaracja formatu pola:
// 				0 float
//				1 double
//				2 short
//				3 int
//				4 complex<float>
//              5 complex<double>
//
// <format_spec[xr*yr*zr]> field   pole
//
//
// kolejno:
// pole[0,0,0] ... pole[0,0,zr-1],pole[0,1,0] ... pole[0,1,zr-1],..
// az do
// pole[xr-1,yr-1,0] ...pole[xr-1][yr-1][zr-1]
//
//////////////////////////////////////////////////////////////////////////////////////

// znajduje ( przez interpolacje) wartosc pola w srodku komorki Yee
typ_pola cenEx(int ci,int cj,int ck)
{
	return 0.25*(Ex[ci+i_0+1][cj+j_0+1][ck+k_0+1]+Ex[ci+i_0+1][cj+j_0+1+1][ck+k_0+1]+Ex[ci+i_0+1][cj+j_0+1][ck+k_0+1+1]+Ex[ci+i_0+1][cj+j_0+1+1][ck+k_0+1+1]);
}

typ_pola cenEy(int ci,int cj,int ck)
{
	return 0.25*(Ey[ci+i_0+1][cj+j_0+1][ck+k_0+1]+Ey[ci+i_0+1+1][cj+j_0+1][ck+k_0+1]+Ey[ci+i_0+1+1][cj+j_0+1][ck+k_0+1+1]+Ey[ci+i_0+1][cj+j_0+1][ck+k_0+1+1]);
}

typ_pola cenEz(int ci,int cj,int ck)
{
	return 0.25*(Ez[ci+i_0+1][cj+j_0+1][ck+k_0+1]+Ez[ci+i_0+1+1][cj+j_0+1+1][ck+k_0+1]+Ez[ci+i_0+1+1][cj+j_0+1][ck+k_0+1]+Ez[ci+i_0+1][cj+j_0+1+1][ck+k_0+1]);
}

typ_pola cenHx(int ci,int cj,int ck)
{
	return 0.5*(Hx[ci+i_0+1-1][cj+j_0+1][ck+k_0+1]+Hx[ci+i_0+1+1-1][cj+j_0+1][ck+k_0+1])/H_scale;
}

typ_pola cenHy(int ci,int cj,int ck)
{
	return 0.5*(Hy[ci+i_0+1][cj+j_0+1-1][ck+k_0+1]+Hy[ci+i_0+1][cj+j_0+1+1-1][ck+k_0+1])/H_scale;
}

typ_pola cenHz(int ci,int cj,int ck)
{
	return 0.5*(Hz[ci+i_0+1][cj+j_0+1][ck+k_0+1-1]+Hz[ci+i_0+1][cj+j_0+1][ck+k_0+1+1-1])/H_scale;
}

#ifdef COMPLEX_FIELD
       typ_pola cenSx_all(int ci,int cj,int ck) {return (typ_pola)(conj(cenEy(ci,cj,ck))*cenHz(ci,cj,ck)-conj(cenEz(ci,cj,ck))*cenHy(ci,cj,ck));}
       typ_pola cenSy_all(int ci,int cj,int ck) {return (typ_pola)(conj(cenEz(ci,cj,ck))*cenHx(ci,cj,ck)-conj(cenEx(ci,cj,ck))*cenHz(ci,cj,ck));}
       typ_pola cenSz_all(int ci,int cj,int ck) {return (typ_pola)(conj(cenEx(ci,cj,ck))*cenHy(ci,cj,ck)-conj(cenEy(ci,cj,ck))*cenHx(ci,cj,ck));}

       typ_pola cenSx_null_Ez_Hy(int ci,int cj,int ck) {return (typ_pola)(conj(cenEy(ci,cj,ck))*cenHz(ci,cj,ck));}
       typ_pola cenSy_null_Ex_Hz(int ci,int cj,int ck) {return (typ_pola)(conj(cenEz(ci,cj,ck))*cenHx(ci,cj,ck));}
       typ_pola cenSz_null_Ey_Hx(int ci,int cj,int ck) {return (typ_pola)(conj(cenEx(ci,cj,ck))*cenHy(ci,cj,ck));}

       typ_pola cenSx_null_Ey_Hz(int ci,int cj,int ck) {return (typ_pola)(-conj(cenEz(ci,cj,ck))*cenHy(ci,cj,ck));}
       typ_pola cenSy_null_Ez_Hx(int ci,int cj,int ck) {return (typ_pola)(-conj(cenEx(ci,cj,ck))*cenHz(ci,cj,ck));}
       typ_pola cenSz_null_Ex_Hy(int ci,int cj,int ck) {return (typ_pola)(-conj(cenEy(ci,cj,ck))*cenHx(ci,cj,ck));}
#else
     typ_pola cenSx_all(int ci,int cj,int ck) {return cenEy(ci,cj,ck)*cenHz(ci,cj,ck)-cenEz(ci,cj,ck)*cenHy(ci,cj,ck);}
     typ_pola cenSy_all(int ci,int cj,int ck) {return cenEz(ci,cj,ck)*cenHx(ci,cj,ck)-cenEx(ci,cj,ck)*cenHz(ci,cj,ck);}
     typ_pola cenSz_all(int ci,int cj,int ck) {return cenEx(ci,cj,ck)*cenHy(ci,cj,ck)-cenEy(ci,cj,ck)*cenHx(ci,cj,ck);}

     typ_pola cenSx_null_Ez_Hy(int ci,int cj,int ck) {return cenEy(ci,cj,ck)*cenHz(ci,cj,ck);}
     typ_pola cenSy_null_Ex_Hz(int ci,int cj,int ck) {return cenEz(ci,cj,ck)*cenHx(ci,cj,ck);}
     typ_pola cenSz_null_Ey_Hx(int ci,int cj,int ck) {return cenEx(ci,cj,ck)*cenHy(ci,cj,ck);}

     typ_pola cenSx_null_Ey_Hz(int ci,int cj,int ck) {return -cenEz(ci,cj,ck)*cenHy(ci,cj,ck);}
     typ_pola cenSy_null_Ez_Hx(int ci,int cj,int ck) {return -cenEx(ci,cj,ck)*cenHz(ci,cj,ck);}
     typ_pola cenSz_null_Ex_Hy(int ci,int cj,int ck) {return -cenEy(ci,cj,ck)*cenHx(ci,cj,ck);}
#endif
// wyrzucanie wynikow do pliku pdp - por. specyfikacja pliku PDP powyzej
bool put_out_field1(int n_out,char* output_files_prefix,char* name,typ_pola (*cen)(int,int,int),int slice,int slicenumber1,int slicenumber2)
{
	std::ostringstream outss;
	std::string outs;
	std::ofstream output_file;

	int tmp1;
	typ_pola tmp2;
	char buffer[30];
	
	int startx=0,endx=xr,starty=0,endy=yr,startz=0,endz=zr;

	outss<<output_files_prefix<<"_"<<n_out<<"_"<<name<<"_"<< ((slice==0) ? "x" : ( (slice==1) ? "y":"z") )<<"="<<slicenumber1<<"-"<<slicenumber2<<"_t="<<t<<".pdp";
	outs=outss.str();

	output_file.open(outs.c_str(),std::ios::out | std::ios::binary);
	if (!output_file) return false;

    sprintf(buffer,"%s,  t = %d",name,t); 

	tmp1=strlen(buffer);output_file.write((char*)&tmp1,sizeof(int));
 	output_file.write(buffer,tmp1);

	tmp1=slicenumber2-slicenumber1+1;

	if (slice==0) {startx=slicenumber1;endx=slicenumber2+1;} else { if (slice==1) {starty=slicenumber1;endy=slicenumber2+1;} else {startz=slicenumber1;endz=slicenumber2+1;}}

	tmp1=(endx-startx)/di[n_out];output_file.write((char*)&tmp1,sizeof(int));
	tmp1=(endy-starty)/dj[n_out];output_file.write((char*)&tmp1,sizeof(int));
	tmp1=(endz-startz)/dk[n_out];output_file.write((char*)&tmp1,sizeof(int));

	//zapis pola typu typ_prec lub complex<typ_prec>
	#ifdef COMPLEX_FIELD
	       #ifdef DOUBLE_PRECISION
                tmp1=5;
	       #else
	            tmp1=4;
	       #endif
    #else
           #ifdef DOUBLE_PRECISION
                tmp1=1;
	       #else
	            tmp1=0;
	       #endif
    #endif
	output_file.write((char*)&tmp1,sizeof(int));

	for(i=startx;i<endx;i+=di[n_out]) for(j=starty;j<endy;j+=dj[n_out]) for(k=startz;k<endz;k+=dk[n_out])
	{
		tmp2=cen(i,j,k);
		output_file.write((char*)&tmp2,sizeof(typ_pola));
	}

	output_file.close();

	return true;

}

bool put_out_field2(int n_out,char* output_files_prefix,char* name,typ_pola (*cen)(int,int,int),int slice,int slicenumber1,int slicenumber2)
{
	std::ostringstream outss;
	std::string outs;
	std::ifstream input_file;

	std::ostringstream outss2;
	std::string outs2;
	std::ofstream output_file;

	int tmp1;
	typ_pola tmp2;
    char buffer[30];

	int startx=0,endx=xr,starty=0,endy=yr,startz=0,endz=zr;

	outss<<output_files_prefix<<"_"<<n_out<<"_"<<name<<"_"<< ((slice==0) ? "x" : ( (slice==1) ? "y":"z") )<<"="<<slicenumber1<<"-"<<slicenumber2<<"_t="<<t<<".pdp";
	outs=outss.str();

	outss2<<output_files_prefix<<"_"<<n_out<<"_"<<name<<"_"<< ((slice==0) ? "x" : ( (slice==1) ? "y":"z") )<<"="<<slicenumber1<<"-"<<slicenumber2<<"_temp.pdp";
	outs2=outss2.str();

	input_file.open(outs.c_str(),std::ios::in | std::ios::binary);if (!input_file) return false;
	output_file.open(outs2.c_str(),std::ios::out | std::ios::binary);if (!output_file) return false;

    sprintf(buffer,"%s,  t = %d",name,t); 

	tmp1=strlen(buffer);output_file.write((char*)&tmp1,sizeof(int));
 	output_file.write(buffer,tmp1);
 	

	tmp1=slicenumber2-slicenumber1+1;

	if (slice==0) {startx=slicenumber1;endx=slicenumber2+1;} else { if (slice==1) {starty=slicenumber1;endy=slicenumber2+1;} else {startz=slicenumber1;endz=slicenumber2+1;}}

	tmp1=(endx-startx)/di[n_out];output_file.write((char*)&tmp1,sizeof(int));
	tmp1=(endy-starty)/dj[n_out];output_file.write((char*)&tmp1,sizeof(int));
	tmp1=(endz-startz)/dk[n_out];output_file.write((char*)&tmp1,sizeof(int));


	//zapis pola typu typ_prec
	#ifdef COMPLEX_FIELD
	       #ifdef DOUBLE_PRECISION
                tmp1=5;
	       #else
	            tmp1=4;
	       #endif
    #else
           #ifdef DOUBLE_PRECISION
                tmp1=1;
	       #else
	            tmp1=0;
	       #endif
    #endif
	output_file.write((char*)&tmp1,sizeof(int));
	
	input_file.seekg(output_file.tellp());

	for(i=startx;i<endx;i+=di[n_out]) for(j=starty;j<endy;j+=dj[n_out]) for(k=startz;k<endz;k+=dk[n_out])
	{
		input_file.read((char*)&tmp2,sizeof(typ_pola));
		tmp2=(tmp2+cen(i,j,k))/typ_prec(2.0);
		output_file.write((char*)&tmp2,sizeof(typ_pola));
	}

	input_file.close();
	output_file.close();
	
	if (remove(outs.c_str())!=0) return false;
	if (rename(outs2.c_str(),outs.c_str())!=0) return false;

	return true;

}

bool put_out_intensity(int n_out,char* output_files_prefix,int slice,int slicenumber1,int slicenumber2)
{
	std::ostringstream outss;std::ostringstream outss1;std::ostringstream outss2;std::ostringstream outss3;
	std::string outs1;std::string outs2;std::string outs3;
	std::ofstream output_file;
	std::ifstream input1_file;std::ifstream input2_file;std::ifstream input3_file;

	int tmp1;
	typ_pola tmp2;
	char buffer[30];

	int startx=0,endx=xr,starty=0,endy=yr,startz=0,endz=zr;

	outss<<output_files_prefix<<"_"<<n_out<<"_S_"<< ((slice==0) ? "x" : ( (slice==1) ? "y":"z") )<<"="<<slicenumber1<<"-"<<slicenumber2<<"_t="<<t<<".pdp";
	outs1=outss.str();
	output_file.open(outs1.c_str(),std::ios::out | std::ios::binary);if (!output_file) return false;

	outss1<<output_files_prefix<<"_"<<n_out<<"_Sx_"<< ((slice==0) ? "x" : ( (slice==1) ? "y":"z") )<<"="<<slicenumber1<<"-"<<slicenumber2<<"_t="<<t<<".pdp";
	outs1=outss1.str();
	input1_file.open(outs1.c_str(),std::ios::in | std::ios::binary);if (!input1_file) return false;

	outss2<<output_files_prefix<<"_"<<n_out<<"_Sy_"<< ((slice==0) ? "x" : ( (slice==1) ? "y":"z") )<<"="<<slicenumber1<<"-"<<slicenumber2<<"_t="<<t<<".pdp";
	outs2=outss2.str();
	input2_file.open(outs2.c_str(),std::ios::in | std::ios::binary);if (!input2_file) return false;

	outss3<<output_files_prefix<<"_"<<n_out<<"_Sz_"<< ((slice==0) ? "x" : ( (slice==1) ? "y":"z") )<<"="<<slicenumber1<<"-"<<slicenumber2<<"_t="<<t<<".pdp";
	outs3=outss3.str();
	input3_file.open(outs3.c_str(),std::ios::in | std::ios::binary);if (!input3_file) return false;

    sprintf(buffer,"S,  t = %d",t); 

	tmp1=strlen(buffer);output_file.write((char*)&tmp1,sizeof(int));
 	output_file.write(buffer,tmp1);

	tmp1=slicenumber2-slicenumber1+1;

	if (slice==0) {startx=slicenumber1;endx=slicenumber2+1;} else { if (slice==1) {starty=slicenumber1;endy=slicenumber2+1;} else {startz=slicenumber1;endz=slicenumber2+1;}}

	tmp1=(endx-startx)/di[n_out];output_file.write((char*)&tmp1,sizeof(int));
	tmp1=(endy-starty)/dj[n_out];output_file.write((char*)&tmp1,sizeof(int));
	tmp1=(endz-startz)/dk[n_out];output_file.write((char*)&tmp1,sizeof(int));


 	//zapis pola typu typ_prec lub complex<float>
	#ifdef COMPLEX_FIELD
	       #ifdef DOUBLE_PRECISION
                tmp1=5;
	       #else
	            tmp1=4;
	       #endif
    #else
           #ifdef DOUBLE_PRECISION
                tmp1=1;
	       #else
	            tmp1=0;
	       #endif
    #endif
	output_file.write((char*)&tmp1,sizeof(int));

	input1_file.read((char*)&tmp1,sizeof(int));
	input1_file.seekg(tmp1*sizeof(char)+4*sizeof(int),std::ios_base::cur);
	
	input2_file.seekg(input1_file.tellg());
	input3_file.seekg(input1_file.tellg());

	typ_pola tmp3;

	for(i=startx;i<endx;i+=di[n_out]) for(j=starty;j<endy;j+=dj[n_out]) for(k=startz;k<endz;k+=dk[n_out])
	{
		input1_file.read((char*)&tmp2,sizeof(typ_pola));
		tmp3=pow((tmp2+cenSx(i,j,k))/typ_prec(2.0),2);

		input2_file.read((char*)&tmp2,sizeof(typ_pola));
		tmp3=tmp3+pow((tmp2+cenSy(i,j,k))/typ_prec(2.0),2);

		input3_file.read((char*)&tmp2,sizeof(typ_pola));
		tmp3=tmp3+pow((tmp2+cenSz(i,j,k))/typ_prec(2.0),2);

		tmp2=sqrt(tmp3);

		output_file.write((char*)&tmp2,sizeof(typ_pola));
	}
	output_file.close();

	input1_file.close();
	input2_file.close();
	input3_file.close();

	if (!skladowe[n_out][6]) {if (remove(outs1.c_str())!=0) return false;}
	if (!skladowe[n_out][7]) {if (remove(outs2.c_str())!=0) return false;}
	if (!skladowe[n_out][8]) {if (remove(outs3.c_str())!=0) return false;}

	return true;
}


bool put_out_modifed(const char* output_file_name,const char* input_file_name,const char* tmp_file_name,typ_prec czynnik,int tk)
{
	std::fstream input_file;std::fstream output_file;std::ofstream tmp_file;
	int tmp1;typ_pola tmp2;typ_pola tmp2o;int endi;int endj;int endk;

	tmp_file.open(tmp_file_name,std::ios::out | std::ios::binary);if (!tmp_file) return false;

	input_file.open(input_file_name,std::ios::in |std::ios::binary);if (!input_file) return false;
	input_file.read((char*)&tmp1,sizeof(int));
	{
		char tmp1a[tmp1];input_file.read((char*)&tmp1a,tmp1);
		char buffer[5]="mod";tmp1=strlen(buffer);
		tmp_file.write((char*)&tmp1,sizeof(int));
		tmp_file.write(buffer,tmp1);
	}

	input_file.read((char*)&endi,sizeof(int));tmp_file.write((char*)&endi,sizeof(int));
	input_file.read((char*)&endj,sizeof(int));tmp_file.write((char*)&endj,sizeof(int));
	input_file.read((char*)&endk,sizeof(int));tmp_file.write((char*)&endk,sizeof(int));
	input_file.read((char*)&tmp1,sizeof(int));tmp_file.write((char*)&tmp1,sizeof(int));

	if (tk==0)
	{
		for(i=0;i<endi;i++) for(j=0;j<endj;j++) for(k=0;k<endk;k++) 
		{
			input_file.read((char*)&tmp2,sizeof(typ_pola));
			tmp2=czynnik*tmp2;tmp_file.write((char*)&tmp2,sizeof(typ_pola));
		}
		remove(output_file_name); // kasujemy ewentualny poprzedni plik
	}
	else
	{
		output_file.open(output_file_name,std::ios::in |std::ios::out | std::ios::binary);if (!output_file) return false;
		output_file.seekp(tmp_file.tellp());
		for(i=0;i<endi;i++) for(j=0;j<endj;j++) for(k=0;k<endk;k++)
		{

			output_file.read((char*)&tmp2o,sizeof(typ_pola));input_file.read((char*)&tmp2,sizeof(typ_pola));
			tmp2=tmp2o+czynnik*tmp2;tmp_file.write((char*)&tmp2,sizeof(typ_pola));
		}
		output_file.close();
		if (remove(output_file_name)!=0) return false;
	}
	
	tmp_file.close();input_file.close();
	if (rename(tmp_file_name,output_file_name)!=0) return false;
	return true;
}


bool put_out_averaged(int n_out,char* output_files_prefix,char* name,int slice,int slicenumber1,int slicenumber2)
{
	std::ostringstream outss1;std::ostringstream outss3;std::ostringstream outss4;
	std::string outs1;std::string outs3;std::string outs4;

	outss1<<output_files_prefix<<"_"<<n_out<<"_"<<name<<"_"<< ((slice==0) ? "x" : ( (slice==1) ? "y":"z") )<<"="<<slicenumber1<<"-"<<slicenumber2<<"_AVG.pdp";
	outss3<<output_files_prefix<<"_"<<n_out<<"_"<<name<<"_"<< ((slice==0) ? "x" : ( (slice==1) ? "y":"z") )<<"="<<slicenumber1<<"-"<<slicenumber2<<"_t="<<t<<".pdp";
	outss4<<output_files_prefix<<"_"<<n_out<<"_"<<name<<"_"<< ((slice==0) ? "x" : ( (slice==1) ? "y":"z") )<<"="<<slicenumber1<<"-"<<slicenumber2<<"_temp.pdp";
	outs1=outss1.str();outs3=outss3.str();outs4=outss4.str();
	if (!put_out_modifed(outs1.c_str(),outs3.c_str(),outs4.c_str(),1.0/floor(1.0+(out_max_t[n_out]-out_min_t[n_out])/out_t_step[n_out]),(t-out_min_t[n_out]))) return false;
	if (remove(outs3.c_str())!=0) return false;
	return true;
}

bool put_out_fourie(int n_out,char* output_files_prefix,char* name,int slice,int slicenumber1,int slicenumber2)
{
	std::ostringstream outss1;std::ostringstream outss2;std::ostringstream outss3;std::ostringstream outss4;
	std::string outs1;std::string outs2;std::string outs3;std::string outs4;

	outss1<<output_files_prefix<<"_"<<n_out<<"_"<<name<<"_"<< ((slice==0) ? "x" : ( (slice==1) ? "y":"z") )<<"="<<slicenumber1<<"-"<<slicenumber2<<"_Freal.pdp";
	outss2<<output_files_prefix<<"_"<<n_out<<"_"<<name<<"_"<< ((slice==0) ? "x" : ( (slice==1) ? "y":"z") )<<"="<<slicenumber1<<"-"<<slicenumber2<<"_Fimag.pdp";
	outss3<<output_files_prefix<<"_"<<n_out<<"_"<<name<<"_"<< ((slice==0) ? "x" : ( (slice==1) ? "y":"z") )<<"="<<slicenumber1<<"-"<<slicenumber2<<"_t="<<t<<".pdp";
	outss4<<output_files_prefix<<"_"<<n_out<<"_"<<name<<"_"<< ((slice==0) ? "x" : ( (slice==1) ? "y":"z") )<<"="<<slicenumber1<<"-"<<slicenumber2<<"_temp.pdp";
	outs1=outss1.str();outs2=outss2.str();outs3=outss3.str();outs4=outss4.str();

	if (!put_out_modifed(outs1.c_str(),outs3.c_str(),outs4.c_str(),cos(ni_fourie[n_out]*t)/floor((out_max_t[n_out]-out_min_t[n_out]+1.0)/out_t_step[n_out]),(t-out_min_t[n_out]))) return false;
	if (!put_out_modifed(outs2.c_str(),outs3.c_str(),outs4.c_str(),sin(ni_fourie[n_out]*t)/floor((out_max_t[n_out]-out_min_t[n_out]+1.0)/out_t_step[n_out]),(t-out_min_t[n_out]))) return false;
	if (remove(outs3.c_str())!=0) return false;
	return true;
}

bool put_out_temp(int n_out,char* output_files_prefix,char* name,int slice,int slicenumber1,int slicenumber2)
{
	std::ostringstream outss3;
	std::string outs3;

	std::ostringstream outss4;
	std::string outs4;

	outss3<<output_files_prefix<<"_"<<n_out<<"_"<<name<<"_"<<((slice==0) ? "x" : ( (slice==1) ? "y":"z") )<<"="<<slicenumber1<<"-"<<slicenumber2<<"_t="<<t<<".pdp";
	outs3=outss3.str();

	outss4<<output_files_prefix<<"_"<<n_out<<"_"<<name<<"_"<<((slice==0) ? "x" : ( (slice==1) ? "y":"z") )<<"="<<slicenumber1<<"-"<<slicenumber2<<"_TMP.pdp";
	outs4=outss4.str();
    remove(outs4.c_str());
	
	if (rename(outs3.c_str(),outs4.c_str())!=0) return false;
	
	return true;
}

bool put_out_detector(char* output_files_prefix,int detector_numer)
{

	std::ostringstream outss;
	std::string outs;
	std::ofstream output_file;

	int tmp1;
	typ_pola tmp2;
    int size_x=(int)((detector_t_end[detector_numer]-detector_t_start[detector_numer])/detector_t_step[detector_numer]);
    outss<<output_files_prefix<<"_detector_"<<detector_numer<<"_"<<skladname[detector_skladowa[detector_numer]]<<".pdp";
    
    outs=outss.str();
	output_file.open(outs.c_str(),std::ios::out | std::ios::binary);
	if (!output_file) return false;

	tmp1=strlen("detector");output_file.write((char*)&tmp1,sizeof(int));
 	output_file.write("detector",tmp1);
	tmp1=size_x;output_file.write((char*)&tmp1,sizeof(int));
	tmp1=1;output_file.write((char*)&tmp1,sizeof(int));
	tmp1=1;output_file.write((char*)&tmp1,sizeof(int));
 	
	#ifdef COMPLEX_FIELD
	       #ifdef DOUBLE_PRECISION
                tmp1=5;
	       #else
	            tmp1=4;
	       #endif
    #else
           #ifdef DOUBLE_PRECISION
                tmp1=1;
	       #else
	            tmp1=0;
	       #endif
    #endif
    
	
	output_file.write((char*)&tmp1,sizeof(int));

	for(i=0;i<size_x;i++)
	{
		tmp2=recorded_at_detector[detector_numer][i];
		output_file.write((char*)&tmp2,sizeof(typ_pola));
	}
	output_file.close();

	return true;
}



// funkcje do zrzucania calych tablic
bool put_out_bin_whole(typ_pola*** tab,char* output_files_prefix,char* name,int size_x,int size_y,int size_z)
{
	if (!silentmode && !iosilentmode) std::cout<<"  *writing: "<<name<<" ...\n";

	std::ostringstream outss;
	std::string outs;
	std::ofstream output_file;

	int tmp1;
	typ_pola tmp2;

	outss<<output_files_prefix<<"_"<<name<<"_backup.pdp";
	outs=outss.str();
	output_file.open(outs.c_str(),std::ios::out | std::ios::binary);
	if (!output_file) return false;

	tmp1=strlen(name);output_file.write((char*)&tmp1,sizeof(int));
 	output_file.write(name,tmp1);
	tmp1=size_x;output_file.write((char*)&tmp1,sizeof(int));
	tmp1=size_y;output_file.write((char*)&tmp1,sizeof(int));
	tmp1=size_z;output_file.write((char*)&tmp1,sizeof(int));

	#ifdef COMPLEX_FIELD
	       #ifdef DOUBLE_PRECISION
                tmp1=5;
	       #else
	            tmp1=4;
	       #endif
    #else
           #ifdef DOUBLE_PRECISION
                tmp1=1;
	       #else
	            tmp1=0;
	       #endif
    #endif
    
	output_file.write((char*)&tmp1,sizeof(int));

	for(i=0;i<size_x;i++) for(j=0;j<size_y;j++) for(k=0;k<size_z;k++) {tmp2=tab[i][j][k]; output_file.write((char*)&tmp2,sizeof(typ_pola));}
	output_file.close();

	return true;

}

bool put_out_bin_whole(typ_prec* tab,char* output_files_prefix,char* name,int size_x)
{
	if (!silentmode && !iosilentmode) std::cout<<"  *writing: "<<name<<" ...\n";

	std::ostringstream outss;
	std::string outs;
	std::ofstream output_file;

	int tmp1;
	typ_prec tmp2;

	outss<<output_files_prefix<<"_"<<name<<"_backup.pdp";
	outs=outss.str();
	output_file.open(outs.c_str(),std::ios::out | std::ios::binary);
	if (!output_file) return false;

	tmp1=strlen(name);output_file.write((char*)&tmp1,sizeof(int));
 	output_file.write(name,tmp1);
	tmp1=size_x;output_file.write((char*)&tmp1,sizeof(int));
	tmp1=1;output_file.write((char*)&tmp1,sizeof(int));
	tmp1=1;output_file.write((char*)&tmp1,sizeof(int));

	//zapis pola typu typ_prec

           #ifdef DOUBLE_PRECISION
                tmp1=1;
	       #else
	            tmp1=0;
	       #endif

	output_file.write((char*)&tmp1,sizeof(int));

	for(i=0;i<size_x;i++)
	{
		tmp2=tab[i];
		output_file.write((char*)&tmp2,sizeof(typ_prec));
	}
	output_file.close();

	return true;
}

#ifdef COMPLEX_FIELD
bool put_out_bin_whole(std::complex<typ_prec>* tab,char* output_files_prefix,char* name,int size_x)
{
	if (!silentmode && !iosilentmode) std::cout<<"  *writing: "<<name<<" ...\n";

	std::ostringstream outss;
	std::string outs;
	std::ofstream output_file;

	int tmp1;
	std::complex<typ_prec> tmp2;

	outss<<output_files_prefix<<"_"<<name<<"_backup.pdp";
	outs=outss.str();
	output_file.open(outs.c_str(),std::ios::out | std::ios::binary);
	if (!output_file) return false;

	tmp1=strlen(name);output_file.write((char*)&tmp1,sizeof(int));
 	output_file.write(name,tmp1);
	tmp1=size_x;output_file.write((char*)&tmp1,sizeof(int));
	tmp1=1;output_file.write((char*)&tmp1,sizeof(int));
	tmp1=1;output_file.write((char*)&tmp1,sizeof(int));

	//zapis pola 
	       #ifdef DOUBLE_PRECISION
                tmp1=5;
	       #else
	            tmp1=4;
	       #endif
	output_file.write((char*)&tmp1,sizeof(int));

	for(i=0;i<size_x;i++)
	{
		tmp2=tab[i];
		output_file.write((char*)&tmp2,sizeof(std::complex<typ_prec>));
	}
	output_file.close();

	return true;
}
#endif

bool read_bin_whole(typ_pola*** tab,char* prefix,char* filename,int size_x,int size_y,int size_z)
{
	if (!silentmode && !iosilentmode) std::cout<<"  *reading: "<<filename<<" ...\n";

	int x,y,z;
	int tmp1;typ_pola tmp2;

	std::ostringstream outss;
	std::string outs;

	outss<<prefix<<"_"<<filename<<"_backup.pdp";
	outs=outss.str();

	std::ifstream plik(outs.c_str(),std::ios::in | std::ios::binary); if (!plik) return false;

	plik.read((char*)&tmp1,sizeof(int));
	plik.seekg(tmp1*sizeof(char),std::ios_base::cur);


	plik.read((char*)&tmp1,sizeof(int));x=tmp1;
	plik.read((char*)&tmp1,sizeof(int));y=tmp1;
	plik.read((char*)&tmp1,sizeof(int));z=tmp1;
	if((x!=size_x)||(y!=size_y)||(z!=size_z)) return false;
	plik.read((char*)&tmp1,sizeof(int));

	for(i=0;i<size_x;i++) for(j=0;j<size_y;j++) for(k=0;k<size_z;k++) { plik.read((char*)&tmp2,sizeof(typ_pola)); tab[i][j][k]=tmp2; }
	plik.close();

	return true;
}

bool read_bin_whole(typ_prec* tab,char* prefix,char* filename,int size_x)
{
	if (!silentmode && !iosilentmode) std::cout<<"  *reading: "<<filename<<" ...\n";

	int x,y,z;
	int tmp1;typ_prec tmp2;

	std::ostringstream outss;
	std::string outs;

	outss<<prefix<<"_"<<filename<<"_backup.pdp";
	outs=outss.str();

    	std::ifstream plik(outs.c_str(),std::ios::in | std::ios::binary);
    	if (!plik) return false;

	plik.read((char*)&tmp1,sizeof(int));
	plik.seekg(tmp1*sizeof(char),std::ios_base::cur);

	plik.read((char*)&tmp1,sizeof(int));x=tmp1;
	plik.read((char*)&tmp1,sizeof(int));y=tmp1;
	plik.read((char*)&tmp1,sizeof(int));z=tmp1;

	if((x>size_x)||(y!=1)||(z!=1)) return false;
	plik.read((char*)&tmp1,sizeof(int));

	for(i=0;i<x;i++) {plik.read((char*)&tmp2,sizeof(typ_prec)); tab[i]=tmp2; }
	plik.close();

	return true;

}

#ifdef COMPLEX_FIELD
bool read_bin_whole(std::complex<typ_prec>* tab,char* prefix,char* filename,int size_x)
{
	if (!silentmode && !iosilentmode) std::cout<<"  *reading: "<<filename<<" ...\n";

	int x,y,z;
	int tmp1;std::complex<typ_prec> tmp2;

	std::ostringstream outss;
	std::string outs;

	outss<<prefix<<"_"<<filename<<"_backup.pdp";
	outs=outss.str();

    std::ifstream plik(outs.c_str(),std::ios::in | std::ios::binary);
    if (!plik) return false;

	plik.read((char*)&tmp1,sizeof(int));
	plik.seekg(tmp1*sizeof(char),std::ios_base::cur);

	plik.read((char*)&tmp1,sizeof(int));x=tmp1;
	plik.read((char*)&tmp1,sizeof(int));y=tmp1;
	plik.read((char*)&tmp1,sizeof(int));z=tmp1;

	if((x>size_x)||(y!=1)||(z!=1)) return false;
	plik.read((char*)&tmp1,sizeof(int));

	for(i=0;i<x;i++) {plik.read((char*)&tmp2,sizeof(std::complex<typ_prec>)); tab[i]=tmp2; }
	plik.close();

	return true;

}
#endif

bool read_bin(double** tab,char* filename,int size_x,int size_y)
{
	if (!silentmode&&!iosilentmode) std::cout<<"    *reading source shape file: "<<filename<<" ...\n";

	int x,y,z;
	int tmp1;double tmp2;typ_prec tmp2f;

    std::ifstream plik(filename,std::ios::in | std::ios::binary);
    if (!plik) return false;
	plik.read((char*)&tmp1,sizeof(int));
	plik.seekg(tmp1*sizeof(char),std::ios_base::cur);

	plik.read((char*)&tmp1,sizeof(int));x=tmp1;
	plik.read((char*)&tmp1,sizeof(int));y=tmp1;
	plik.read((char*)&tmp1,sizeof(int));z=tmp1;
	if((x!=size_x)||(y!=size_y)||(z!=1)) return false;
	plik.read((char*)&tmp1,sizeof(int));

	if (tmp1==0) {for(i=0;i<size_x;i++) for(j=0;j<size_y;j++) {plik.read((char*)&tmp2f,sizeof(double)); tab[i][j]=tmp2f;}}
	else {for(i=0;i<size_x;i++) for(j=0;j<size_y;j++) {plik.read((char*)&tmp2,sizeof(double)); tab[i][j]=tmp2;}}

	plik.close();

	return true;
}

bool read_bin(double* tab,char* filename,int size_x)
{
	if (!silentmode&&!iosilentmode) std::cout<<"    *reading source time dependency file: "<<filename<<" ...\n";

	int x,y,z;
	int tmp1;double tmp2;typ_prec tmp2f;

    std::ifstream plik(filename,std::ios::in | std::ios::binary);
    if (!plik) return false;
    
	plik.read((char*)&tmp1,sizeof(int));
	plik.seekg(tmp1*sizeof(char),std::ios_base::cur);

	plik.read((char*)&tmp1,sizeof(int));x=tmp1;
	plik.read((char*)&tmp1,sizeof(int));y=tmp1;
	plik.read((char*)&tmp1,sizeof(int));z=tmp1;
	if((x!=size_x)||(y!=1)||(z!=1)) return false;
	plik.read((char*)&tmp1,sizeof(int));

	if (tmp1==0) {for(i=0;i<size_x;i++) {plik.read((char*)&tmp2f,sizeof(float)); tab[i]=tmp2f;}}
	else {for(i=0;i<size_x;i++) {plik.read((char*)&tmp2,sizeof(double)); tab[i]=tmp2;}}

	plik.close();

	return true;
}

// wyrzucanie wynikow do plikow - po kolei poszczegolne skladowe
bool put_out1(char* output_files_prefix,int n_out)
{

	for (int u=3;u<6;u++) if (skladowe[n_out][u]) if (!put_out_field1(n_out,output_files_prefix,skladname[u],cenfun[u],slice[n_out],out_min_slice[n_out],out_max_slice[n_out])) return false;
	for (int u=6;u<9;u++) if (skladowe[n_out][u] || skladowe[n_out][9]) if (!put_out_field1(n_out,output_files_prefix,skladname[u],cenfun[u],slice[n_out],out_min_slice[n_out],out_max_slice[n_out])) return false;

	return true;
}

bool put_out2(char* output_files_prefix,int n_out)
{

	for (int u=0;u<3;u++) if (skladowe[n_out][u]) if (!put_out_field1(n_out,output_files_prefix,skladname[u],cenfun[u],slice[n_out],out_min_slice[n_out],out_max_slice[n_out])) return false;
	for (int u=3;u<6;u++) if (skladowe[n_out][u]) if (!put_out_field2(n_out,output_files_prefix,skladname[u],cenfun[u],slice[n_out],out_min_slice[n_out],out_max_slice[n_out])) return false;

	if (skladowe[n_out][9])	if (!put_out_intensity(n_out,output_files_prefix,slice[n_out],out_min_slice[n_out],out_max_slice[n_out])) return false;

	for (int u=6;u<9;u++) if (skladowe[n_out][u]) if (!put_out_field2(n_out,output_files_prefix,skladname[u],cenfun[u],slice[n_out],out_min_slice[n_out],out_max_slice[n_out])) return false;

	if (is_fourie[n_out]) for (int u=0;u<10;u++) if(skladowe[n_out][u]) if(!put_out_fourie(n_out,output_files_prefix,skladname[u],slice[n_out],out_min_slice[n_out],out_max_slice[n_out])) return false;
	if (is_averaged[n_out]) for (int u=0;u<10;u++) if(skladowe[n_out][u]) if(!put_out_averaged(n_out,output_files_prefix,skladname[u],slice[n_out],out_min_slice[n_out],out_max_slice[n_out])) return false;
	if (is_temp[n_out]) for (int u=0;u<10;u++) if(skladowe[n_out][u]) if(!put_out_temp(n_out,output_files_prefix,skladname[u],slice[n_out],out_min_slice[n_out],out_max_slice[n_out])) return false;

	return true;
}

std::string report(char* prefix,char* struct_filename,char* param_filename,bool ress,bool st=false)
{
	std::ostringstream output;

	for (int u=0;u<99;u++) output<<"-";
	if (ress) output<<"\n| ponowne uruchomienie z kopi w kroku : "<<ts;
	output<<"\n| nazwa pliku z struktura: "<<struct_filename<<" ("<<struct_name<<"), nazwa pliku z parametrami: "<<param_filename<<"\n";
	output<<"| rozmiary struktury: "<<xr<<" x "<<yr<<" x "<<zr<<"\n";
	for (int s=0;s<n_SRCS;s++)
	{
     output<<"| typ zrodla nr "<<s<<": "<<((source_type[s]==0) ? "fala plaska" : ((source_type[s]==1) ? "wiazka" :((source_type[s]==2) ? "wiazka o ksztalcie wczytanym z pliku" : ((source_type[s]>2 && source_type[s]<8) ? "dipol" : "brak" ))))<<"\n";
    }
   	for (int s=0;s<n_SRCS;s++)
	{
     output<<"| dlugosci fali zrodla nr "<<s<<": "<<lambda[s]<<"\n";
    }
    output<<"| krok przestrzenny dr = "<<dr<<", krok czasowy dt = "<<dt<<"\n";
    output<<"| liczba krokow czasowych = "<<maxt<<"\n";
    output<<"| warunki brzegowe :  ";
    output<<"  x : "<<((wx==0) ? "UPML" : ((wx==1) ? "periodyczne" :((wx==2) ? "symetryczne" : "Bloch" )));
    output<<", y : "<<((wy==0) ? "UPML" : ((wy==1) ? "periodyczne" :((wy==2) ? "symetryczne" : "Bloch" )));
    output<<", z : "<<((wz==0) ? "UPML" : ((wz==1) ? "periodyczne" :((wz==2) ? "symetryczne" : "Bloch" )));
    output<<"\n";
	output<<"| pliki z wynikami ("<<prefix<<"):"<<"\n";

	for(i=0;i<n_outputs;i++)
	{
		output<<"|  "<<i<<")"<<" zapis wynikow dla t = "<<out_min_t[i]<<" - "<<out_max_t[i]<<" co " <<out_t_step[i]<<", ";
		output<<((slice[i]==0) ? "x" : ( (slice[i]==1) ? "y":"z") )<<" = "<<out_min_slice[i]<<" -  "<<out_max_slice[i]<<", ";
		output<<" z krokami "<<di[i]<<" x "<<dj[i]<<" x "<<dk[i]<<"\n";
		if (is_averaged[i]) {output<<"|  srednia  :";}
		else if (is_fourie[i]) {output<<"|  DFT dla czestotliwosci [1/dt] = "<<ni_fourie[i]/(2*M_PI)<<" :";}
		else if (is_temp[i]) {output<<"|  zapis tymczasowy z skladowymi :";}
		else{output<<"|  skladowe :";}
		output<<(skladowe[i][0] ? " Ex":"")<<(skladowe[i][1] ? " Ey":"")<<(skladowe[i][2] ? " Ez":"");
		output<<(skladowe[i][3] ? " Hx":"")<<(skladowe[i][4] ? " Hy":"")<<(skladowe[i][5] ? " Hz":"");
		output<<(skladowe[i][6] ? " Sx":"")<<(skladowe[i][7] ? " Sy":"")<<(skladowe[i][8] ? " Sz":"");
		output<<(skladowe[i][9] ? " S":"")<<"\n";
	}
	output<<"| detektory ("<<prefix<<"):"<<"\n";
	for(i=0;i<detector_n;i++)
	{
		output<<"|  "<<i<<")"<<" zapis "<<skladname[detector_skladowa[i]]<<" dla t = "<<detector_t_start[i]<<" - "<<detector_t_end[i]<<" co " <<detector_t_step[i]<<" polozenie: x:"<<detector_x[i]<<" y:"<<detector_y[i]<<" z:"<<detector_z[i]<<"\n";
	}

	output<<"| uruchomienie programu : "<<ctime(&strtime);
	if (!st) output<<"| zakonczenie           : "<<ctime(&endtime);
	for (int u=0;u<99;u++) output<<"-";
	return output.str();
}

// tworzenie pliku z raportem - do rozbudowy
bool create_report(char* prefix,char* struct_filename,char* param_filename,bool ress)
{
	std::ostringstream outss;
	std::string outs;
	std::ofstream output_file;

	outss<<prefix<<"_report";

	outs=outss.str();

	if (!ress) output_file.open(outs.c_str());else output_file.open(outs.c_str(),std::ios::out | std::ios::app);;
	if (!output_file) return false;
	output_file<<report(prefix,struct_filename,param_filename,ress);
	output_file.close();
	return true;
}

bool load_all(char* prefix)
{

    if (is_calc_ex) if (!read_bin_whole(Ex,prefix,"Ex",xs,ys+1,zs+1)) return false;
    if (is_calc_ey) if (!read_bin_whole(Ey,prefix,"Ey",xs+1,ys,zs+1)) return false;
    if (is_calc_ez) if (!read_bin_whole(Ez,prefix,"Ez",xs+1,ys+1,zs)) return false;
    if (is_calc_hx) if (!read_bin_whole(Hx,prefix,"Hx",xs-1,ys,zs)) return false;
    if (is_calc_hy) if (!read_bin_whole(Hy,prefix,"Hy",xs,ys-1,zs)) return false;
    if (is_calc_hz) if (!read_bin_whole(Hz,prefix,"Hz",xs,ys,zs-1)) return false;

if (wy!=1 && wy!=3)
{
    if (is_calc_ex) if (!read_bin_whole(bot_Dx,prefix,"bot_Dx",xs,boundy,zr+1+2*bndz+2*disPML)) return false;
    if (is_calc_ey) if (!read_bin_whole(bot_Dy,prefix,"bot_Dy",xs+1,boundy,zr+1+2*bndz+2*disPML)) return false;
	
	if (is_calc_ez)
    { 
     if (!read_bin_whole(bot_Dz_lef,prefix,"bot_Dz_lef",boundx,boundy,zr+2*bndz+2*disPML)) return false;
     if (!read_bin_whole(bot_Dz_rig,prefix,"bot_Dz_rig",boundx,boundy,zr+2*bndz+2*disPML)) return false;
    }
    
    if (is_calc_hx) if (!read_bin_whole(bot_Bx,prefix,"bot_Bx",xs-1,boundy,zr+2*bndz+2*disPML)) return false;
    if (is_calc_hy) if (!read_bin_whole(bot_By,prefix,"bot_By",xs,boundy,zr+2*bndz+2*disPML)) return false;
	
	if (is_calc_hz)
    { 
     if (!read_bin_whole(bot_Bz_lef,prefix,"bot_Bz_lef",boundx,boundy,zr-1+2*bndz+2*disPML)) return false;
     if (!read_bin_whole(bot_Bz_rig,prefix,"bot_Bz_rig",boundx,boundy,zr-1+2*bndz+2*disPML)) return false;
     }
}

if (wy==0)
{
    if (is_calc_ex) if (!read_bin_whole(top_Dx,prefix,"top_Dx",xs,boundy,zr+1+2*bndz+2*disPML)) return false;
    if (is_calc_ey) if (!read_bin_whole(top_Dy,prefix,"top_Dy",xs+1,boundy,zr+1+2*bndz+2*disPML)) return false;
	
	if (is_calc_ez)
    { 
      if (!read_bin_whole(top_Dz_lef,prefix,"top_Dz_lef",boundx,boundy,zr+2*bndz+2*disPML)) return false;
      if (!read_bin_whole(top_Dz_rig,prefix,"top_Dz_rig",boundx,boundy,zr+2*bndz+2*disPML)) return false;
    }
    
    if (is_calc_hx) if (!read_bin_whole(top_Bx,prefix,"top_Bx",xs-1,boundy,zr+2*bndz+2*disPML)) return false;
    if (is_calc_hy) if (!read_bin_whole(top_By,prefix,"top_By",xs,boundy,zr+2*bndz+2*disPML)) return false;
	
	if (is_calc_hz) 
	{
     if (!read_bin_whole(top_Bz_lef,prefix,"top_Bz_lef",boundx,boundy,zr-1+2*bndz+2*disPML)) return false;
     if (!read_bin_whole(top_Bz_rig,prefix,"top_Bz_rig",boundx,boundy,zr-1+2*bndz+2*disPML)) return false;
    }
}

if (wx!=1 && wx!=3)
{
    if (is_calc_ex) if (!read_bin_whole(lef_Dx,prefix,"lef_Dx",boundx,yr+1+2*bndy+2*disPML,zr+1+2*bndz+2*disPML)) return false;
    if (is_calc_ez) if (!read_bin_whole(lef_Dz,prefix,"lef_Dz",boundx,yr+1+2*bndy+2*disPML,zr+2*bndz+2*disPML)) return false;
    if (is_calc_hx) if (!read_bin_whole(lef_Bx,prefix,"lef_Bx",boundx,yr+2*bndy+2*disPML,zr+2*bndz+2*disPML)) return false;
    if (is_calc_hz) if (!read_bin_whole(lef_Bz,prefix,"lef_Bz",boundx,yr+2*bndy+2*disPML,zr-1+2*bndz+2*disPML)) return false;
}

if (wx==0)
{
    if (is_calc_ex) if (!read_bin_whole(rig_Dx,prefix,"rig_Dx",boundx,yr+1+2*bndy+2*disPML,zr+1+2*bndz+2*disPML)) return false;
    if (is_calc_ez) if (!read_bin_whole(rig_Dz,prefix,"rig_Dz",boundx,yr+1+2*bndy+2*disPML,zr+2*bndz+2*disPML)) return false;
    
    if (is_calc_hx) if (!read_bin_whole(rig_Bx,prefix,"rig_Bx",boundx,yr+2*bndy+2*disPML,zr+2*bndz+2*disPML)) return false;
    if (is_calc_hz) if (!read_bin_whole(rig_Bz,prefix,"rig_Bz",boundx,yr+2*bndy+2*disPML,zr-1+2*bndz+2*disPML)) return false;
}

if (wz==0)
{
    if (is_calc_ex)
    {       
     if (!read_bin_whole(fro_Dx_bot,prefix,"fro_Dx_bot",xs,boundy,boundz)) return false;
	 if (!read_bin_whole(fro_Dx_top,prefix,"fro_Dx_top",xs,boundy,boundz)) return false;
	 if (!read_bin_whole(fro_Dx_lef,prefix,"fro_Dx_lef",boundx,yr+1+2*bndy+2*disPML,boundz)) return false;
	 if (!read_bin_whole(fro_Dx_rig,prefix,"fro_Dx_rig",boundx,yr+1+2*bndy+2*disPML,boundz)) return false;
    }
	if (is_calc_ey) if (!read_bin_whole(fro_Dy,prefix,"fro_Dy",xs+1,ys,boundz)) return false;
    if (is_calc_ez) if (!read_bin_whole(fro_Dz,prefix,"fro_Dz",xs+1,ys+1,boundz)) return false;

    if (is_calc_hx)
    { 
	 if (!read_bin_whole(fro_Bx_bot,prefix,"fro_Bx_bot",xs-1,boundy,boundz)) return false;
	 if (!read_bin_whole(fro_Bx_top,prefix,"fro_Bx_top",xs-1,boundy,boundz)) return false;
	 if (!read_bin_whole(fro_Bx_lef,prefix,"fro_Bx_lef",boundx,yr+2*bndy+2*disPML,boundz)) return false;
	 if (!read_bin_whole(fro_Bx_rig,prefix,"fro_Bx_rig",boundx,yr+2*bndy+2*disPML,boundz)) return false;
    }
	if (is_calc_hy) if (!read_bin_whole(fro_By,prefix,"fro_By",xs,ys-1,boundz)) return false;
    if (is_calc_hz) if (!read_bin_whole(fro_Bz,prefix,"fro_Bz",xs,ys,boundz)) return false;
}

if (wz!=1 && wz!=3)
{
    if (is_calc_ex)
    { 
     if (!read_bin_whole(bac_Dx_bot,prefix,"bac_Dx_bot",xs,boundy,boundz)) return false;
	 if (!read_bin_whole(bac_Dx_top,prefix,"bac_Dx_top",xs,boundy,boundz)) return false;

	 if (!read_bin_whole(bac_Dx_lef,prefix,"bac_Dx_lef",boundx,yr+1+2*bndy+2*disPML,boundz)) return false;
	 if (!read_bin_whole(bac_Dx_rig,prefix,"bac_Dx_rig",boundx,yr+1+2*bndy+2*disPML,boundz)) return false;
    }
    if (is_calc_ey) if (!read_bin_whole(bac_Dy,prefix,"bac_Dy",xs+1,ys,boundz)) return false;
    if (is_calc_ez) if (!read_bin_whole(bac_Dz,prefix,"bac_Dz",xs+1,ys+1,boundz)) return false;

    if (is_calc_hx)
    {
	 if (!read_bin_whole(bac_Bx_bot,prefix,"bac_Bx_bot",xs-1,boundy,boundz)) return false;
	 if (!read_bin_whole(bac_Bx_top,prefix,"bac_Bx_top",xs-1,boundy,boundz)) return false;
     if (!read_bin_whole(bac_Bx_lef,prefix,"bac_Bx_lef",boundx,yr+2*bndy+2*disPML,boundz)) return false;
	 if (!read_bin_whole(bac_Bx_rig,prefix,"bac_Bx_rig",boundx,yr+2*bndy+2*disPML,boundz)) return false;
    }
    
	if (is_calc_hy) if (!read_bin_whole(bac_By,prefix,"bac_By",xs,ys-1,boundz)) return false;
    if (is_calc_hz) if (!read_bin_whole(bac_Bz,prefix,"bac_Bz",xs,ys,boundz)) return false;
}
    for (int s=0;s<n_SRCS;s++)
    {
     if (source_type[s]==0)
     {
                        char buffer[10];	
                        sprintf(buffer,"i_E_%d",s); 
                        if (!read_bin_whole(i_E[s],prefix,buffer,pulse_length+2)) return false;
                        sprintf(buffer,"i_H_%d",s);                         
                        if (!read_bin_whole(i_H[s],prefix,buffer,pulse_length+1)) return false;
     }
    }
     
	if (n_drd_Ex>0 && is_calc_ex)
	{
		if (!read_bin_whole(drd_Jex,prefix,"drd_Jex",n_drd_Ex)) return false;if (!read_bin_whole(drd_Ex_n_1,prefix,"drd_Ex_n_1",n_drd_Ex)) return false;
    } 
	if (n_drd_Ey>0 && is_calc_ey)
	{
		if (!read_bin_whole(drd_Jey,prefix,"drd_Jey",n_drd_Ey)) return false;if (!read_bin_whole(drd_Ey_n_1,prefix,"drd_Ey_n_1",n_drd_Ey)) return false;
    } 
   	if (n_drd_Ez>0 && is_calc_ez)
	{
		if (!read_bin_whole(drd_Jez,prefix,"drd_Jez",n_drd_Ez)) return false;if (!read_bin_whole(drd_Ez_n_1,prefix,"drd_Ex_n_1",n_drd_Ez)) return false;
    } 
    
	if (n_drd_Hx>0 && is_calc_hx)
	{
		if (!read_bin_whole(drd_Jhx,prefix,"drd_Jhx",n_drd_Hx)) return false;if (!read_bin_whole(drd_Hx_n_1,prefix,"drd_Hx_n_1",n_drd_Hx)) return false;
    } 
	if (n_drd_Hy>0 && is_calc_hy)
	{
		if (!read_bin_whole(drd_Jhy,prefix,"drd_Jhy",n_drd_Hy)) return false;if (!read_bin_whole(drd_Hy_n_1,prefix,"drd_Hy_n_1",n_drd_Hy)) return false;
    } 
   	if (n_drd_Hz>0 && is_calc_hz)
	{
		if (!read_bin_whole(drd_Jhz,prefix,"drd_Jhz",n_drd_Hz)) return false;if (!read_bin_whole(drd_Hz_n_1,prefix,"drd_Hx_n_1",n_drd_Hz)) return false;
    }         

	if (n_lor_Ex>0 && is_calc_ex)
	{
		if (!read_bin_whole(lor_Jex,prefix,"lor_Jex",n_lor_Ex)) return false;if (!read_bin_whole(lor_Jex_n_1,prefix,"lor_Jex_n_1",n_lor_Ex)) return false;
        if (!read_bin_whole(lor_Ex_n_1,prefix,"lor_Ex_n_1",n_lor_Ex)) return false;if (!read_bin_whole(lor_Ex_n_2,prefix,"lor_Ex_n_2",n_lor_Ex)) return false;
    }
	if (n_lor_Ey>0 && is_calc_ey)
	{
		if (!read_bin_whole(lor_Jey,prefix,"lor_Jey",n_lor_Ey)) return false;if (!read_bin_whole(lor_Jey_n_1,prefix,"lor_Jey_n_1",n_lor_Ey)) return false;
        if (!read_bin_whole(lor_Ey_n_1,prefix,"lor_Ey_n_1",n_lor_Ey)) return false;if (!read_bin_whole(lor_Ey_n_2,prefix,"lor_Ey_n_2",n_lor_Ey)) return false;
    }
   	if (n_lor_Ez>0 && is_calc_ez)
	{
		if (!read_bin_whole(lor_Jez,prefix,"lor_Jez",n_lor_Ez)) return false;if (!read_bin_whole(lor_Jez_n_1,prefix,"lor_Jez_n_1",n_lor_Ez)) return false;
        if (!read_bin_whole(lor_Ez_n_1,prefix,"lor_Ez_n_1",n_lor_Ez)) return false;if (!read_bin_whole(lor_Ez_n_2,prefix,"lor_Ez_n_2",n_lor_Ez)) return false;
    }

	if (n_lor_Hx>0 && is_calc_hx)
	{
		if (!read_bin_whole(lor_Jhx,prefix,"lor_Jhx",n_lor_Hx)) return false;if (!read_bin_whole(lor_Jhx_n_1,prefix,"lor_Jhx_n_1",n_lor_Hx)) return false;
        if (!read_bin_whole(lor_Hx_n_1,prefix,"lor_Hx_n_1",n_lor_Hx)) return false;if (!read_bin_whole(lor_Hx_n_2,prefix,"lor_Hx_n_2",n_lor_Hx)) return false;
    }
	if (n_lor_Hy>0 && is_calc_hy)
	{
		if (!read_bin_whole(lor_Jhy,prefix,"lor_Jhy",n_lor_Hy)) return false;if (!read_bin_whole(lor_Jhy_n_1,prefix,"lor_Jhy_n_1",n_lor_Hy)) return false;
        if (!read_bin_whole(lor_Hy_n_1,prefix,"lor_Hy_n_1",n_lor_Hy)) return false;if (!read_bin_whole(lor_Hy_n_2,prefix,"lor_Hy_n_2",n_lor_Hy)) return false;
    }
   	if (n_lor_Hz>0 && is_calc_hz)
	{
		if (!read_bin_whole(lor_Jhz,prefix,"lor_Jhz",n_lor_Hz)) return false;if (!read_bin_whole(lor_Jhz_n_1,prefix,"lor_Jhz_n_1",n_lor_Hz)) return false;
        if (!read_bin_whole(lor_Hz_n_1,prefix,"lor_Hz_n_1",n_lor_Hz)) return false;if (!read_bin_whole(lor_Hz_n_2,prefix,"lor_Hz_n_2",n_lor_Hz)) return false;
    }

	if (n_dby_Ex>0 && is_calc_ex)
	{
		if (!read_bin_whole(dby_Jex,prefix,"dby_Jex",n_dby_Ex)) return false;if (!read_bin_whole(dby_Ex_n_1,prefix,"dby_Ex_n_1",n_dby_Ex)) return false;
    } 
	if (n_dby_Ey>0 && is_calc_ey)
	{
		if (!read_bin_whole(dby_Jey,prefix,"dby_Jey",n_dby_Ey)) return false;if (!read_bin_whole(dby_Ey_n_1,prefix,"dby_Ey_n_1",n_dby_Ey)) return false;
    } 
   	if (n_dby_Ez>0 && is_calc_ez)
	{
		if (!read_bin_whole(dby_Jez,prefix,"dby_Jez",n_dby_Ez)) return false;if (!read_bin_whole(dby_Ez_n_1,prefix,"dby_Ex_n_1",n_dby_Ez)) return false;
    } 
    
	if (n_dby_Hx>0 && is_calc_hx)
	{
		if (!read_bin_whole(dby_Jhx,prefix,"dby_Jhx",n_dby_Hx)) return false;if (!read_bin_whole(dby_Hx_n_1,prefix,"dby_Hx_n_1",n_dby_Hx)) return false;
    } 
	if (n_dby_Hy>0 && is_calc_hy)
	{
		if (!read_bin_whole(dby_Jhy,prefix,"dby_Jhy",n_dby_Hy)) return false;if (!read_bin_whole(dby_Hy_n_1,prefix,"dby_Hy_n_1",n_dby_Hy)) return false;
    } 
   	if (n_dby_Hz>0 && is_calc_hz)
	{
		if (!read_bin_whole(dby_Jhz,prefix,"dby_Jhz",n_dby_Hz)) return false;if (!read_bin_whole(dby_Hz_n_1,prefix,"dby_Hx_n_1",n_dby_Hz)) return false;
    }         

    if (detector_n>0)
    {
     for (i=0;i<detector_n;i++)
     {
        char buffer[20];	
        sprintf(buffer,"detector_%d",i); 
        if (!read_bin_whole(recorded_at_detector[i],prefix,buffer,(int)((detector_t_end[i]-detector_t_start[i])/detector_t_step[i]))+1) return false;
      }                
    }


	std::ostringstream outss;
	std::string outs;
	std::ifstream input_file;

	outss<<prefix<<"_backup";
	outs=outss.str();

	input_file.open(outs.c_str());
	input_file>>ts;
	input_file.close();
	return true;

}

bool save_all(char** args)
{
    if (is_calc_ex) if (!put_out_bin_whole(Ex,args[3],"Ex",xs,ys+1,zs+1)) return false;
    if (is_calc_ey) if (!put_out_bin_whole(Ey,args[3],"Ey",xs+1,ys,zs+1)) return false;
    if (is_calc_ez) if (!put_out_bin_whole(Ez,args[3],"Ez",xs+1,ys+1,zs)) return false;
    if (is_calc_hx) if (!put_out_bin_whole(Hx,args[3],"Hx",xs-1,ys,zs)) return false;
    if (is_calc_hy) if (!put_out_bin_whole(Hy,args[3],"Hy",xs,ys-1,zs)) return false;
    if (is_calc_hz) if (!put_out_bin_whole(Hz,args[3],"Hz",xs,ys,zs-1)) return false;

if (wy!=1 && wy!=3)
{
    if (is_calc_ex) if (!put_out_bin_whole(bot_Dx,args[3],"bot_Dx",xs,boundy,zr+1+2*bndz+2*disPML)) return false;
    if (is_calc_ey) if (!put_out_bin_whole(bot_Dy,args[3],"bot_Dy",xs+1,boundy,zr+1+2*bndz+2*disPML)) return false;
	
	if (is_calc_ez)
    { 
     if (!put_out_bin_whole(bot_Dz_lef,args[3],"bot_Dz_lef",boundx,boundy,zr+2*bndz+2*disPML)) return false;
     if (!put_out_bin_whole(bot_Dz_rig,args[3],"bot_Dz_rig",boundx,boundy,zr+2*bndz+2*disPML)) return false;
    }
    
    if (is_calc_hx) if (!put_out_bin_whole(bot_Bx,args[3],"bot_Bx",xs-1,boundy,zr+2*bndz+2*disPML)) return false;
    if (is_calc_hy) if (!put_out_bin_whole(bot_By,args[3],"bot_By",xs,boundy,zr+2*bndz+2*disPML)) return false;
	
	if (is_calc_hz)
    { 
     if (!put_out_bin_whole(bot_Bz_lef,args[3],"bot_Bz_lef",boundx,boundy,zr-1+2*bndz+2*disPML)) return false;
     if (!put_out_bin_whole(bot_Bz_rig,args[3],"bot_Bz_rig",boundx,boundy,zr-1+2*bndz+2*disPML)) return false;
     }
}

if (wy==0)
{
    if (is_calc_ex) if (!put_out_bin_whole(top_Dx,args[3],"top_Dx",xs,boundy,zr+1+2*bndz+2*disPML)) return false;
    if (is_calc_ey) if (!put_out_bin_whole(top_Dy,args[3],"top_Dy",xs+1,boundy,zr+1+2*bndz+2*disPML)) return false;
	
	if (is_calc_ez)
    { 
      if (!put_out_bin_whole(top_Dz_lef,args[3],"top_Dz_lef",boundx,boundy,zr+2*bndz+2*disPML)) return false;
      if (!put_out_bin_whole(top_Dz_rig,args[3],"top_Dz_rig",boundx,boundy,zr+2*bndz+2*disPML)) return false;
    }
    
    if (is_calc_hx) if (!put_out_bin_whole(top_Bx,args[3],"top_Bx",xs-1,boundy,zr+2*bndz+2*disPML)) return false;
    if (is_calc_hy) if (!put_out_bin_whole(top_By,args[3],"top_By",xs,boundy,zr+2*bndz+2*disPML)) return false;
	
	if (is_calc_hz) 
	{
     if (!put_out_bin_whole(top_Bz_lef,args[3],"top_Bz_lef",boundx,boundy,zr-1+2*bndz+2*disPML)) return false;
     if (!put_out_bin_whole(top_Bz_rig,args[3],"top_Bz_rig",boundx,boundy,zr-1+2*bndz+2*disPML)) return false;
    }
}

if (wx!=1 && wx!=3)
{
    if (is_calc_ex) if (!put_out_bin_whole(lef_Dx,args[3],"lef_Dx",boundx,yr+1+2*bndy+2*disPML,zr+1+2*bndz+2*disPML)) return false;
    if (is_calc_ez) if (!put_out_bin_whole(lef_Dz,args[3],"lef_Dz",boundx,yr+1+2*bndy+2*disPML,zr+2*bndz+2*disPML)) return false;
    if (is_calc_hx) if (!put_out_bin_whole(lef_Bx,args[3],"lef_Bx",boundx,yr+2*bndy+2*disPML,zr+2*bndz+2*disPML)) return false;
    if (is_calc_hz) if (!put_out_bin_whole(lef_Bz,args[3],"lef_Bz",boundx,yr+2*bndy+2*disPML,zr-1+2*bndz+2*disPML)) return false;
}

if (wx==0)
{
    if (is_calc_ex) if (!put_out_bin_whole(rig_Dx,args[3],"rig_Dx",boundx,yr+1+2*bndy+2*disPML,zr+1+2*bndz+2*disPML)) return false;
    if (is_calc_ez) if (!put_out_bin_whole(rig_Dz,args[3],"rig_Dz",boundx,yr+1+2*bndy+2*disPML,zr+2*bndz+2*disPML)) return false;
    
    if (is_calc_hx) if (!put_out_bin_whole(rig_Bx,args[3],"rig_Bx",boundx,yr+2*bndy+2*disPML,zr+2*bndz+2*disPML)) return false;
    if (is_calc_hz) if (!put_out_bin_whole(rig_Bz,args[3],"rig_Bz",boundx,yr+2*bndy+2*disPML,zr-1+2*bndz+2*disPML)) return false;
}

if (wz==0)
{
    if (is_calc_ex)
    {       
     if (!put_out_bin_whole(fro_Dx_bot,args[3],"fro_Dx_bot",xs,boundy,boundz)) return false;
	 if (!put_out_bin_whole(fro_Dx_top,args[3],"fro_Dx_top",xs,boundy,boundz)) return false;
	 if (!put_out_bin_whole(fro_Dx_lef,args[3],"fro_Dx_lef",boundx,yr+1+2*bndy+2*disPML,boundz)) return false;
	 if (!put_out_bin_whole(fro_Dx_rig,args[3],"fro_Dx_rig",boundx,yr+1+2*bndy+2*disPML,boundz)) return false;
    }
	if (is_calc_ey) if (!put_out_bin_whole(fro_Dy,args[3],"fro_Dy",xs+1,ys,boundz)) return false;
    if (is_calc_ez) if (!put_out_bin_whole(fro_Dz,args[3],"fro_Dz",xs+1,ys+1,boundz)) return false;

    if (is_calc_hx)
    { 
	 if (!put_out_bin_whole(fro_Bx_bot,args[3],"fro_Bx_bot",xs-1,boundy,boundz)) return false;
	 if (!put_out_bin_whole(fro_Bx_top,args[3],"fro_Bx_top",xs-1,boundy,boundz)) return false;
	 if (!put_out_bin_whole(fro_Bx_lef,args[3],"fro_Bx_lef",boundx,yr+2*bndy+2*disPML,boundz)) return false;
	 if (!put_out_bin_whole(fro_Bx_rig,args[3],"fro_Bx_rig",boundx,yr+2*bndy+2*disPML,boundz)) return false;
    }
	if (is_calc_hy) if (!put_out_bin_whole(fro_By,args[3],"fro_By",xs,ys-1,boundz)) return false;
    if (is_calc_hz) if (!put_out_bin_whole(fro_Bz,args[3],"fro_Bz",xs,ys,boundz)) return false;
}

if (wz!=1 && wz!=3)
{
    if (is_calc_ex)
    { 
     if (!put_out_bin_whole(bac_Dx_bot,args[3],"bac_Dx_bot",xs,boundy,boundz)) return false;
	 if (!put_out_bin_whole(bac_Dx_top,args[3],"bac_Dx_top",xs,boundy,boundz)) return false;

	 if (!put_out_bin_whole(bac_Dx_lef,args[3],"bac_Dx_lef",boundx,yr+1+2*bndy+2*disPML,boundz)) return false;
	 if (!put_out_bin_whole(bac_Dx_rig,args[3],"bac_Dx_rig",boundx,yr+1+2*bndy+2*disPML,boundz)) return false;
    }
    if (is_calc_ey) if (!put_out_bin_whole(bac_Dy,args[3],"bac_Dy",xs+1,ys,boundz)) return false;
    if (is_calc_ez) if (!put_out_bin_whole(bac_Dz,args[3],"bac_Dz",xs+1,ys+1,boundz)) return false;

    if (is_calc_hx)
    {
	 if (!put_out_bin_whole(bac_Bx_bot,args[3],"bac_Bx_bot",xs-1,boundy,boundz)) return false;
	 if (!put_out_bin_whole(bac_Bx_top,args[3],"bac_Bx_top",xs-1,boundy,boundz)) return false;
     if (!put_out_bin_whole(bac_Bx_lef,args[3],"bac_Bx_lef",boundx,yr+2*bndy+2*disPML,boundz)) return false;
	 if (!put_out_bin_whole(bac_Bx_rig,args[3],"bac_Bx_rig",boundx,yr+2*bndy+2*disPML,boundz)) return false;
    }
    
	if (is_calc_hy) if (!put_out_bin_whole(bac_By,args[3],"bac_By",xs,ys-1,boundz)) return false;
    if (is_calc_hz) if (!put_out_bin_whole(bac_Bz,args[3],"bac_Bz",xs,ys,boundz)) return false;
}

    for (int s=0;s<n_SRCS;s++)
    {
     if (source_type[s]==0)
     {
        char buffer[10];	
        sprintf(buffer,"i_E_%d",s); 
        if (!put_out_bin_whole(i_E[s],args[3],buffer,pulse_length+2)) return false;
        sprintf(buffer,"i_H_%d",s);         
        if (!put_out_bin_whole(i_H[s],args[3],buffer,pulse_length+1)) return false;
      }
     }
     
	if (n_drd_Ex>0 && is_calc_ex)
	{
		if (!put_out_bin_whole(drd_Jex,args[3],"drd_Jex",n_drd_Ex)) return false;if (!put_out_bin_whole(drd_Ex_n_1,args[3],"drd_Ex_n_1",n_drd_Ex)) return false;
    } 
	if (n_drd_Ey>0 && is_calc_ey)
	{
		if (!put_out_bin_whole(drd_Jey,args[3],"drd_Jey",n_drd_Ey)) return false;if (!put_out_bin_whole(drd_Ey_n_1,args[3],"drd_Ey_n_1",n_drd_Ey)) return false;
    } 
   	if (n_drd_Ez>0 && is_calc_ez)
	{
		if (!put_out_bin_whole(drd_Jez,args[3],"drd_Jez",n_drd_Ez)) return false;if (!put_out_bin_whole(drd_Ez_n_1,args[3],"drd_Ex_n_1",n_drd_Ez)) return false;
    } 
    
	if (n_drd_Hx>0 && is_calc_hx)
	{
		if (!put_out_bin_whole(drd_Jhx,args[3],"drd_Jhx",n_drd_Hx)) return false;if (!put_out_bin_whole(drd_Hx_n_1,args[3],"drd_Hx_n_1",n_drd_Hx)) return false;
    } 
	if (n_drd_Hy>0 && is_calc_hy)
	{
		if (!put_out_bin_whole(drd_Jhy,args[3],"drd_Jhy",n_drd_Hy)) return false;if (!put_out_bin_whole(drd_Hy_n_1,args[3],"drd_Hy_n_1",n_drd_Hy)) return false;
    } 
   	if (n_drd_Hz>0 && is_calc_hz)
	{
		if (!put_out_bin_whole(drd_Jhz,args[3],"drd_Jhz",n_drd_Hz)) return false;if (!put_out_bin_whole(drd_Hz_n_1,args[3],"drd_Hx_n_1",n_drd_Hz)) return false;
    }         

	if (n_lor_Ex>0 && is_calc_ex)
	{
		if (!put_out_bin_whole(lor_Jex,args[3],"lor_Jex",n_lor_Ex)) return false;if (!put_out_bin_whole(lor_Jex_n_1,args[3],"lor_Jex_n_1",n_lor_Ex)) return false;
        if (!put_out_bin_whole(lor_Ex_n_1,args[3],"lor_Ex_n_1",n_lor_Ex)) return false;if (!put_out_bin_whole(lor_Ex_n_2,args[3],"lor_Ex_n_2",n_lor_Ex)) return false;
    }
	if (n_lor_Ey>0 && is_calc_ey)
	{
		if (!put_out_bin_whole(lor_Jey,args[3],"lor_Jey",n_lor_Ey)) return false;if (!put_out_bin_whole(lor_Jey_n_1,args[3],"lor_Jey_n_1",n_lor_Ey)) return false;
        if (!put_out_bin_whole(lor_Ey_n_1,args[3],"lor_Ey_n_1",n_lor_Ey)) return false;if (!put_out_bin_whole(lor_Ey_n_2,args[3],"lor_Ey_n_2",n_lor_Ey)) return false;
    }
   	if (n_lor_Ez>0 && is_calc_ez)
	{
		if (!put_out_bin_whole(lor_Jez,args[3],"lor_Jez",n_lor_Ez)) return false;if (!put_out_bin_whole(lor_Jez_n_1,args[3],"lor_Jez_n_1",n_lor_Ez)) return false;
        if (!put_out_bin_whole(lor_Ez_n_1,args[3],"lor_Ez_n_1",n_lor_Ez)) return false;if (!put_out_bin_whole(lor_Ez_n_2,args[3],"lor_Ez_n_2",n_lor_Ez)) return false;
    }

	if (n_lor_Hx>0 && is_calc_hx)
	{
		if (!put_out_bin_whole(lor_Jhx,args[3],"lor_Jhx",n_lor_Hx)) return false;if (!put_out_bin_whole(lor_Jhx_n_1,args[3],"lor_Jhx_n_1",n_lor_Hx)) return false;
        if (!put_out_bin_whole(lor_Hx_n_1,args[3],"lor_Hx_n_1",n_lor_Hx)) return false;if (!put_out_bin_whole(lor_Hx_n_2,args[3],"lor_Hx_n_2",n_lor_Hx)) return false;
    }
	if (n_lor_Hy>0 && is_calc_hy)
	{
		if (!put_out_bin_whole(lor_Jhy,args[3],"lor_Jhy",n_lor_Hy)) return false;if (!put_out_bin_whole(lor_Jhy_n_1,args[3],"lor_Jhy_n_1",n_lor_Hy)) return false;
        if (!put_out_bin_whole(lor_Hy_n_1,args[3],"lor_Hy_n_1",n_lor_Hy)) return false;if (!put_out_bin_whole(lor_Hy_n_2,args[3],"lor_Hy_n_2",n_lor_Hy)) return false;
    }
   	if (n_lor_Hz>0 && is_calc_hz)
	{
		if (!put_out_bin_whole(lor_Jhz,args[3],"lor_Jhz",n_lor_Hz)) return false;if (!put_out_bin_whole(lor_Jhz_n_1,args[3],"lor_Jhz_n_1",n_lor_Hz)) return false;
        if (!put_out_bin_whole(lor_Hz_n_1,args[3],"lor_Hz_n_1",n_lor_Hz)) return false;if (!put_out_bin_whole(lor_Hz_n_2,args[3],"lor_Hz_n_2",n_lor_Hz)) return false;
    }

	if (n_dby_Ex>0 && is_calc_ex)
	{
		if (!put_out_bin_whole(dby_Jex,args[3],"dby_Jex",n_dby_Ex)) return false;if (!put_out_bin_whole(dby_Ex_n_1,args[3],"dby_Ex_n_1",n_dby_Ex)) return false;
    } 
	if (n_dby_Ey>0 && is_calc_ey)
	{
		if (!put_out_bin_whole(dby_Jey,args[3],"dby_Jey",n_dby_Ey)) return false;if (!put_out_bin_whole(dby_Ey_n_1,args[3],"dby_Ey_n_1",n_dby_Ey)) return false;
    } 
   	if (n_dby_Ez>0 && is_calc_ez)
	{
		if (!put_out_bin_whole(dby_Jez,args[3],"dby_Jez",n_dby_Ez)) return false;if (!put_out_bin_whole(dby_Ez_n_1,args[3],"dby_Ex_n_1",n_dby_Ez)) return false;
    } 
    
	if (n_dby_Hx>0 && is_calc_hx)
	{
		if (!put_out_bin_whole(dby_Jhx,args[3],"dby_Jhx",n_dby_Hx)) return false;if (!put_out_bin_whole(dby_Hx_n_1,args[3],"dby_Hx_n_1",n_dby_Hx)) return false;
    } 
	if (n_dby_Hy>0 && is_calc_hy)
	{
		if (!put_out_bin_whole(dby_Jhy,args[3],"dby_Jhy",n_dby_Hy)) return false;if (!put_out_bin_whole(dby_Hy_n_1,args[3],"dby_Hy_n_1",n_dby_Hy)) return false;
    } 
   	if (n_dby_Hz>0 && is_calc_hz)
	{
		if (!put_out_bin_whole(dby_Jhz,args[3],"dby_Jhz",n_dby_Hz)) return false;if (!put_out_bin_whole(dby_Hz_n_1,args[3],"dby_Hx_n_1",n_dby_Hz)) return false;
    } 
    
    if (detector_n>0)
    {
     for (i=0;i<detector_n;i++)
     {
        char buffer[20];	
        sprintf(buffer,"detector_%d",i); 
        if (!put_out_bin_whole(recorded_at_detector[i],args[3],buffer,(int)((detector_t_end[i]-detector_t_start[i])/detector_t_step[i]))+1) return false;
      }                
    }
    
	std::ostringstream outss;
	std::string outs;
	std::ofstream output_file;

	outss<<args[3]<<"_backup";
	outs=outss.str();

	output_file.open(outs.c_str());
	if (!output_file) return false;

	output_file<<t+1<<"\n";
	output_file<<args[0]<<" "<<args[1]<<" "<<args[2]<<" "<<args[3]<<"\n";
	output_file.close();

	return true;

}


void print_time()
{
	// wyrzucamy aktualny czas,procent wykonania, krok wykonania symulacji
	// czas jest liczony jako srednia wazona z przewidywan na podstawie 
    // a) czasu wykonania ostatniego kroku 
    // b)czasu wykonania wszystkich poprzednich krokow
	
	int auxd;int auxh;int auxm;int auxs;
	std::cout<<100.0*t/maxt<<"% (" << t <<" step(s) made";	
	std::cout<<" ETA:";
	endtime=time(0);
	double wsp=((double)(maxt-t))/maxt;
	auxs=(int)(((difftime(endtime,strtime)*wsp*(maxt-t+ts-1.0))/(t-ts+1.0)+(1-wsp)*difftime(endtime,prevtime)*(maxt-t+ts-1)));	
	auxd=(int)(auxs/(60*60*24));auxs-=auxd*(60*60*24);
	auxh=(int)(auxs/(60*60));auxs-=auxh*(60*60);
	auxm=(int)(auxs/(60));auxs-=auxm*(60);
	std::cout<<auxd<<"d"<<auxh<<"h"<<auxm<<"m"<<auxs<<"s)\n";
	prevtime=endtime;
}


bool read_params(char* params_file_name)
{
	int s_type;int source_time_dep[max_n_SRCS];bool relativ;

	//otwieramy plik z parametrami
	std::ifstream params_file(params_file_name);
	if (!params_file) return false;

    // liczba zrodel - max. max_n_SRCS
    params_file>>n_SRCS;if (n_SRCS>max_n_SRCS) return false;
    
    for (int s=0;s<n_SRCS;s++)
    {
     params_file>>SRCS_shift[s];params_file>>is_SRCS_shift_rel[s];
	 //typ zrodla
	 params_file>>source_time_dep[s];
	 
     if (source_time_dep[s]!=-1)
	 {
            if (source_time_dep[s]!=127)
            {                   
                 //wczytywanie dlugosc fali
                 params_file>>lambda[s];
                 if (source_time_dep[s]==2 || source_time_dep[s]==22) params_file>>lambda2[s];
            }
            else
            {
                char filename[100];
                params_file>>filename;
                params_file>>time_dependency_length[s];
                time_dependency[s]=new(std::nothrow) double [time_dependency_length[s]];
                if (!time_dependency[s]) {params_file.close();return false;};
                if (!read_bin(time_dependency[s],filename,time_dependency_length[s])) {params_file.close();return false;};                
            }
			//wczytwywanie amplitudy impulsu (fali)
			params_file>>ampl[s];
	
	
		// typ zrodla
		params_file>>s_type;

		if (s_type==0)
		{
			source_type[s]=0; // fala plaska
			//kierunek wektora falowego
			params_file>>teta[s]; params_file>>fi[s]; params_file>>psi[s];
			teta[s]*=M_PI;fi[s]*=M_PI;psi[s]*=M_PI;
		}
		else if (((s_type<103) && (s_type>99)) || ((s_type<203) && (s_type>199))) // dipol
		{
			source_type[s]=s_type%100+3;// dipol i jego polozenie
			params_file>>i_source[s];
			params_file>>j_source[s];
			params_file>>k_source[s];
		}
		else if (s_type==127)
		{
			//kierunek wektora falowego
			params_file>>teta[s]; params_file>>fi[s];
			teta[s]*=M_PI;fi[s]*=M_PI;

			source_type[s]=2;

			params_file>>i_source[s];params_file>>j_source[s];params_file>>k_source[s];

			char filename[100];
			int sx,sy;

			params_file>>drs[s];params_file>>sx;params_file>>sy;

			if (!init_array2D(shape_ex_real[s],sx,sy)) return false;
			params_file>>filename;
			if (!read_bin(shape_ex_real[s],filename,sx,sy)) return false;

			if (!init_array2D(shape_ey_real[s],sx,sy)) return false;
			params_file>>filename;
			if (!read_bin(shape_ey_real[s],filename,sx,sy)) return false;

			if (!init_array2D(shape_ex_imag[s],sx,sy)) return false;
			params_file>>filename;
			if (!read_bin(shape_ex_imag[s],filename,sx,sy)) return false;

			if (!init_array2D(shape_ey_imag[s],sx,sy)) return false;
			params_file>>filename;
			if (!read_bin(shape_ey_imag[s],filename,sx,sy)) return false;

			sh_cx[s]=drs[s]*(sx-1)/2.0;
			sh_cy[s]=drs[s]*(sy-1)/2.0;
		}
		else if(s_type==-1) //brak zrodla
		{
			source_type[s]=-1;
		}
		else // zrodlo od 1 do 9
		{
			//kierunek wektora falowego
			params_file>>teta[s]; params_file>>fi[s]; params_file>>psi[s];
			teta[s]*=M_PI;fi[s]*=M_PI;psi[s]*=M_PI;

			source_type[s]=1;// zrodlo ograniczone wczytywane z programu albo zrodla predefiniowane w programie , ograniczone
	
			typ_prec (*shapes_choice[10])(typ_prec,typ_prec,typ_prec,int) = {gaussshape,hermite_gauss01,hermite_gauss10,hermite_gauss11,hermite_gauss02,hermite_gauss20,hermite_gauss12,hermite_gauss21,hermite_gauss22,bessel0};
			shape[s]=shapes_choice[s_type-1];
			params_file>>param1[s];params_file>>param2[s];params_file>>param3[s];
	  		params_file>>i_source[s];params_file>>j_source[s];params_file>>k_source[s];
	  	}
	
	 }
    }
	// krok przestrzenny i czasowy
	params_file>>dr;
	params_file>>relativ;
	if (relativ)
	{
		if (source_time_dep[0]==2 || source_time_dep[0]==22) // jesli zrodlo to sinus modulowany gausem
		{
			dr*=((1/((1/lambda[0])+(1/lambda2[0])))/20);
		}
		else dr*=(lambda[0]/20);
	}

	params_file>>dt;
	params_file>>relativ;
	if (relativ) dt*=dr/(2*light_speed);

	// liczba krokow czasowych
	params_file>>maxt;

	// rozmiary struktury
	params_file>>xr; params_file>>yr; params_file>>zr;
	
	// war. brzegowe 0 - UPML  1 - periodic 2 - symmetric 3 - Bloch
	bloch_kx=0;bloch_ky=0;bloch_kz=0;
	
    params_file>>wx; 
    if (wx==2) 
    {
     params_file>>wx_Ex; params_file>>wx_Ey;params_file>>wx_Ez;
     params_file>>wx_Hx; params_file>>wx_Hy;params_file>>wx_Hz;
    } 
    if (wx==3) params_file>>bloch_kx; 
	
    params_file>>wy; 
    if (wy==2) 
    {
     params_file>>wy_Ex; params_file>>wy_Ey;params_file>>wy_Ez;
     params_file>>wy_Hx; params_file>>wy_Hy;params_file>>wy_Hz;
    } 
    if (wy==3) params_file>>bloch_ky; 
	
    params_file>>wz;
    if (wz==2) 
    {
     params_file>>wz_Ex; params_file>>wz_Ey;params_file>>wz_Ez;
     params_file>>wz_Hx; params_file>>wz_Hy;params_file>>wz_Hz;
    } 
    if (wz==3) params_file>>bloch_kz;

	// tryb 3D ?
	params_file>>mode3D;

	// liczba symulowanych materialow
	params_file>>nmat;
	// parametry materialow
	eps=new(std::nothrow) typ_prec[nmat]; mi=new(std::nothrow) typ_prec [nmat];sig=new(std::nothrow) double [nmat];sih=new(std::nothrow) double [nmat];
	if ((!eps)||(!mi)||(!sig)||(!sih)) {params_file.close();return false;}

	drd_omp_E=new(std::nothrow) long double[nmat];drd_gam_E=new(std::nothrow) long double[nmat];is_drd_E=new(std::nothrow) bool[nmat];
	if ((!drd_omp_E)||(!drd_gam_E)||(!is_drd_E)) {params_file.close();return false;}

	drd_omp_H=new(std::nothrow) long double[nmat];drd_gam_H=new(std::nothrow) long double[nmat];is_drd_H=new(std::nothrow) bool[nmat];is_drd_lor_H=new(std::nothrow) bool[nmat];
	if ((!drd_omp_H)||(!drd_gam_H)||(!is_drd_H)|| (!is_drd_lor_H)) {params_file.close();return false;}

	lor_epd_E=new(std::nothrow) long double[nmat];lor_omp_E=new(std::nothrow) long double[nmat];lor_del_E=new(std::nothrow) long double[nmat];is_lor_E=new(std::nothrow) bool[nmat];is_PEC=new(std::nothrow) bool[nmat];
	if ((!lor_epd_E)||(!lor_omp_E)||(!lor_del_E)||(!is_lor_E)||(!is_PEC)) {params_file.close();return false;}

	lor_alf_E=new(std::nothrow) typ_prec[nmat];lor_ksi_E=new(std::nothrow) typ_prec[nmat];lor_gam_E=new(std::nothrow) typ_prec[nmat];
	if ((!lor_alf_E)||(!lor_ksi_E)||(!lor_gam_E)) {params_file.close();return false;}

	lor_a1_E=new(std::nothrow) typ_prec[nmat];lor_a2_E=new(std::nothrow) typ_prec[nmat];lor_a3_E=new(std::nothrow) typ_prec[nmat];
	if ((!lor_a1_E)||(!lor_a2_E)||(!lor_a3_E)) {params_file.close();return false;}

	lor_epd_H=new(std::nothrow) long double[nmat];lor_omp_H=new(std::nothrow) long double[nmat];lor_del_H=new(std::nothrow) long double[nmat];is_lor_H=new(std::nothrow) bool[nmat];is_PMC=new(std::nothrow) bool[nmat];
	if ((!lor_epd_H)||(!lor_omp_H)||(!lor_del_H)||(!is_lor_H)||(!is_PMC)) {params_file.close();return false;}

	lor_alf_H=new(std::nothrow) typ_prec[nmat];lor_ksi_H=new(std::nothrow) typ_prec[nmat];lor_gam_H=new(std::nothrow) typ_prec[nmat];
	if ((!lor_alf_H)||(!lor_ksi_H)||(!lor_gam_H)) {params_file.close();return false;}

	lor_a1_H=new(std::nothrow) typ_prec[nmat];lor_a2_H=new(std::nothrow) typ_prec[nmat];lor_a3_H=new(std::nothrow) typ_prec[nmat];
	if ((!lor_a1_H)||(!lor_a2_H)||(!lor_a3_H)) {params_file.close();return false;}

	dby_deps_E=new(std::nothrow) long double[nmat];dby_tau_E=new(std::nothrow) long double[nmat];is_dby_E=new(std::nothrow) bool[nmat];
	if ((!dby_deps_E)||(!dby_tau_E)||(!is_dby_E)) {params_file.close();return false;}

	dby_deps_H=new(std::nothrow) long double[nmat];dby_tau_H=new(std::nothrow) long double[nmat];is_dby_H=new(std::nothrow) bool[nmat];
	if ((!dby_deps_H)||(!dby_tau_H)||(!is_dby_H)) {params_file.close();return false;}

    is_drd_lor_E=new(std::nothrow) bool[nmat];is_drd_dby_E=new(std::nothrow) bool[nmat];is_lor_dby_E=new(std::nothrow) bool[nmat];
    if ((!is_drd_lor_E)||(!is_drd_dby_E)||(!is_lor_dby_E)) {params_file.close();return false;}
    
    is_drd_lor_H=new(std::nothrow) bool[nmat];is_drd_dby_H=new(std::nothrow) bool[nmat];is_lor_dby_H=new(std::nothrow) bool[nmat];
    if ((!is_drd_lor_H)||(!is_drd_dby_H)||(!is_lor_dby_H)) {params_file.close();return false;}
    
    
	for(i=0;i<nmat;i++)
	{
			int disp_type;

        	params_file>>eps[i];
        	params_file>>mi[i];
			params_file>>sig[i];
			params_file>>sih[i];

			// zaleznosc dyspersyjna dla eps
			params_file>>disp_type;

			// drude
			if (disp_type==1 || disp_type==12 || disp_type==13)
			{
				params_file>>drd_omp_E[i];
				params_file>>drd_gam_E[i];
				is_drd_E[i]=true;
			}
			else {is_drd_E[i]=false;}

			// lorentz
			if (disp_type==2 || disp_type==12 || disp_type==23)
			{
				params_file>>lor_epd_E[i];
				params_file>>lor_omp_E[i];
				params_file>>lor_del_E[i];
				is_lor_E[i]=true;
			}
			else {is_lor_E[i]=false;}

			// debeye
			if (disp_type==3 || disp_type==13 || disp_type==23)
			{
				params_file>>dby_deps_E[i];
				params_file>>dby_tau_E[i];
				is_dby_E[i]=true;
			}
			else {is_dby_E[i]=false;}

			if (disp_type==12) is_drd_lor_E[i]=true; else  is_drd_lor_E[i]=false; // drude + lorentz
            if (disp_type==13) is_drd_dby_E[i]=true; else  is_drd_dby_E[i]=false; // drude + debeye
            if (disp_type==23) is_lor_dby_E[i]=true; else  is_lor_dby_E[i]=false; // lorentz + debeye 
			
			if (disp_type==100) is_PEC[i]=true; else is_PEC[i]=false;

			// zaleznosc dyspersyjna dla mi
			params_file>>disp_type;

			// drd
			if (disp_type==1 || disp_type==12 || disp_type==13)
			{
				params_file>>drd_omp_H[i];
				params_file>>drd_gam_H[i];
				is_drd_H[i]=true;
			}
			else {is_drd_H[i]=false;}

			// lorentz
			if (disp_type==2 || disp_type==12 || disp_type==23)
			{
				params_file>>lor_epd_H[i];
				params_file>>lor_omp_H[i];
				params_file>>lor_del_H[i];
				is_lor_H[i]=true;
			}
			else {is_lor_H[i]=false;}

			// dby
			if (disp_type==3 || disp_type==13 || disp_type==23)
			{
				params_file>>dby_deps_H[i];
				params_file>>dby_tau_H[i];
				is_dby_H[i]=true;
			}
			else {is_dby_H[i]=false;}
			
			if (disp_type==12) is_drd_lor_H[i]=true; else  is_drd_lor_H[i]=false; // drude + lorentz
            if (disp_type==13) is_drd_dby_H[i]=true; else  is_drd_dby_H[i]=false; // drude + debeye
            if (disp_type==23) is_lor_dby_H[i]=true; else  is_lor_dby_H[i]=false; // lorentz + debeye
 
			if (disp_type==100) is_PMC[i]=true; else is_PMC[i]=false;

	}

	//ilosc krokow przestrzennych na warunki brzegowe
	params_file>>bound;

	//
	boundx= (((wx!=1)&&(wx!=3)) ? bound:1);
	boundy= (((wy!=1)&&(wy!=3)) ? bound:1);
	boundz= (((wz!=1)&&(wz!=3)) ? bound:0);

	bndx= (((wx!=1)&&(wx!=3)) ? 2:1);
	bndy= (((wy!=1)&&(wy!=3)) ? 2:1);
	bndz= (((wz!=1)&&(wz!=3)) ? 2:1);

	// material wypelniajacy przestrzen
	params_file>>bound_mat;

	// liczba outputow
	params_file>>n_outputs;
    if (n_outputs>0)
    {
	 //outputs
	 out_min_t=new(std::nothrow) int[n_outputs];out_max_t=new(std::nothrow) int[n_outputs];
	 if ((!out_min_t)||(!out_max_t)) {params_file.close();return false;}
  
  	 out_t_step=new(std::nothrow) int[n_outputs];
	 if (!out_t_step) {params_file.close();return false;}

	 out_min_slice=new(std::nothrow) int[n_outputs];out_max_slice=new(std::nothrow) int[n_outputs];
	 if ((!out_min_slice)||(!out_max_slice)) {params_file.close();return false;}

	 slice=new(std::nothrow) int[n_outputs];
	 if (!slice) {params_file.close();return false;}

	 di=new(std::nothrow) int[n_outputs];
	 if (!di) {params_file.close();return false;}

	 dj=new(std::nothrow) int[n_outputs];
	 if (!dj) {params_file.close();return false;}

	 dk=new(std::nothrow) int[n_outputs];
	 if (!dk) {params_file.close();return false;}

	 /////////////////////////////////////////////////////
	 skladowe=new(std::nothrow) bool*[n_outputs];
     if (!skladowe) {params_file.close();return false;}

	 skladowe[0]=new(std::nothrow) bool[10*n_outputs];
	 if (!skladowe[0]) {params_file.close();return false;}

     for(i=1;i<n_outputs;i++) skladowe[i]=skladowe[i-1]+10;
	 /////////////////////////////////////////////////////
	 is_averaged=new(std::nothrow) bool[n_outputs];if(!is_averaged) return false;for(int u=0;u<n_outputs;u++) is_averaged[u]=false;

	 is_temp=new(std::nothrow) bool[n_outputs];if(!is_temp) return false;for(int u=0;u<n_outputs;u++) is_temp[u]=false;

	 is_fourie=new(std::nothrow) bool[n_outputs];if(!is_fourie) return false;for(int u=0;u<n_outputs;u++) is_fourie[u]=false;
	 ni_fourie=new(std::nothrow) double[n_outputs];if(!ni_fourie) return false;

	/////////////////////////////////////////////////////
    }
	for(i=0;i<n_outputs;i++)
	{
		int inp_tmp;

		params_file>>out_min_t[i];params_file>>out_max_t[i];params_file>>out_t_step[i];

		params_file>>slice[i];

		params_file>>out_min_slice[i];params_file>>out_max_slice[i];

		params_file>>di[i];params_file>>dj[i];params_file>>dk[i];

		for(j=0;j<10;j++) params_file>>skladowe[i][j];

		params_file>>inp_tmp;
		
		if (inp_tmp==1) 
		{
			is_averaged[i]=true;
		}
		else if (inp_tmp==2)
		{
			is_temp[i]=true;
		}
		else if (inp_tmp==3)
		{		
			is_fourie[i]=true;		
			params_file>>ni_fourie[i];ni_fourie[i]*=2*M_PI;
		}
		

	}

	// liczba outputow
	params_file>>detector_n;
	if (detector_n>0)
	{
     
     detector_skladowa=new(std::nothrow) int[detector_n];
     if (!detector_skladowa) {params_file.close();return false;}              
                     
     detector_x=new(std::nothrow) int[detector_n];
     detector_y=new(std::nothrow) int[detector_n];
     detector_z=new(std::nothrow) int[detector_n];
	 if ((!detector_x)||(!detector_y)||(!detector_z)) {params_file.close();return false;}
	 
     detector_t_start=new(std::nothrow) int[detector_n];
     detector_t_end=new(std::nothrow) int[detector_n];
     detector_t_step=new(std::nothrow) int[detector_n];
	 if ((!detector_t_start)||(!detector_t_end)||(!detector_t_step)) {params_file.close();return false;}
     recorded_at_detector=new(std::nothrow) typ_pola*[detector_n];
     if (!recorded_at_detector) {params_file.close();return false;}	 
     int sum=0;	 
	 
     for (i=0;i<detector_n;i++)
	 {
         params_file>>detector_skladowa[i];
         params_file>>detector_x[i];params_file>>detector_y[i];params_file>>detector_z[i];
         params_file>>detector_t_start[i];params_file>>detector_t_end[i];params_file>>detector_t_step[i];
         sum+=(int)((detector_t_end[i]-detector_t_start[i])/detector_t_step[i]+1);
     } 

	 recorded_at_detector[0]=new(std::nothrow) typ_pola[sum];
 	 if (!recorded_at_detector[0]) {params_file.close();return false;}
     for(i=1;i<detector_n;i++) recorded_at_detector[i]=recorded_at_detector[i-1]+(int)((detector_t_end[i-1]-detector_t_start[i-1])/detector_t_step[i-1])+1;
         
    }

	// wczytujemy odstep czasowy miedzy backup - ami
	params_file>>backup_delta_t;

	params_file>>is_calc_ex;params_file>>is_calc_ey;params_file>>is_calc_ez;
	params_file>>is_calc_hx;params_file>>is_calc_hy;params_file>>is_calc_hz;	

	// zamykamy plik
	params_file.close();

	// ----------------------------

  for(int s=0;s<n_SRCS;s++)
  {
	if (source_time_dep[s]==0)
	{
 		// ramped sinus
		// dlugosc podawanej fali to pop prostu dlugosc fali zrodla
		source_func[s]=sinus;
		freq_cst[s]=2*M_PI*light_speed*dt/lambda[s];
		t_0[s]=0; 
	}
	else if (source_time_dep[s]==1)
	{
		//gauss
		// dlugosc podawanej fali = 3*sigma gaussa po transformacie Fouriera
		source_func[s]=gauss;
		freq_cst[s]=2*pow(M_PI*light_speed*dt/(3*lambda[s]),2);
		t_0[s]=20*lambda[s]/(2*M_PI*light_speed*dt);
	}
	else if(source_time_dep[s]==2)
	{
		// sinus modulowany gaussem
		source_func[s]=singauss;
		freq_cst[s]=2*M_PI*light_speed*dt/lambda[s];
		freq_cst2[s]=2*pow(M_PI*light_speed*dt/(3*lambda2[s]),2);
		t_0[s]=20*lambda2[s]/(2*M_PI*light_speed*dt);
	}
	else if(source_time_dep[s]==3)
	{
		// diffgauss
		source_func[s]=diffgauss;
		freq_cst[s]=pow(M_PI*light_speed*dt/(3*lambda[s]),2);
		freq_cst2[s]=sqrt(2*freq_cst[s]);
		t_0[s]=20*lambda[s]/(2*M_PI*light_speed*dt);
	}
	else if(source_time_dep[s]==127)
	{
		// zaleznosc czasowa wczytana z pliku
		source_func[s]=time_dependency_from_file;
		t_0[s]=0; 
	}
	else if(source_time_dep[s]==10)
	{
		//sinus -- stary , ukryta funkcja :)
		// dlugosc podawanej fali to po prostu dlugosc fali zrodla
		source_func[s]=sinus_old;
		freq_cst[s]=2*M_PI*light_speed*dt/lambda[s];
		t_0[s]=0; 
	}
		if (source_time_dep[s]==20)
	{
 		// ramped sinus complex
		// dlugosc podawanej fali to pop prostu dlugosc fali zrodla
		source_func[s]=sinus_complex;
		freq_cst[s]=2*M_PI*light_speed*dt/lambda[s];
		t_0[s]=0; 
	}
	else if(source_time_dep[s]==22)
	{
		// sinus modulowany gaussem complex
		source_func[s]=singauss_complex;
		freq_cst[s]=2*M_PI*light_speed*dt/lambda[s];
		freq_cst2[s]=2*pow(M_PI*light_speed*dt/(3*lambda2[s]),2);
		t_0[s]=20*lambda2[s]/(2*M_PI*light_speed*dt);
	}
		// zrodlo c.d.
	i_source[s]+=boundx+disPML+bndx;
	j_source[s]+=boundy+disPML+bndy;
	k_source[s]+=boundz+disPML+bndz;

	floor_i_source[s]=(int)floor(i_source[s]);
	floor_j_source[s]=(int)floor(j_source[s]);
	floor_k_source[s]=(int)floor(k_source[s]);

	half_period[s]=lambda[s]/(2*light_speed*dt);
	
    }
    
	ts=1;

	xs=xr+2*(boundx+disPML+bndx);ys=yr+2*(boundy+disPML+bndy);zs=zr+2*(boundz+disPML+bndz);
	xk=xr+boundx+2*(disPML+bndx);yk=yr+boundy+2*(disPML+bndy);zk=zr+boundz+2*(disPML+bndz);
	i_0=boundx+(bndx-1)+disPML;j_0=boundy+(bndy-1)+disPML;k_0=boundz+(bndz-1)+disPML;
	i_1=xr+boundx+(bndx+1)+disPML;j_1=yr+boundy+(bndy+1)+disPML;k_1=zr+boundz+(bndz+1)+disPML;

	H_scale=sqrt((mi_0*mi[bound_mat])/(eps_0*eps[bound_mat]));

	sigma_max=0.8*(wykladnik+1)/(sqrt(mi_0*mi[bound_mat]/(eps_0*eps[bound_mat]))*dr);

	// obliczmy dlugosc tabel pomocniczych - z sporym zapasem
	pulse_length=(int)ceil(sqrt((long double)(xs*xs+ys*ys+zs*zs)))+(int)ceil(light_speed*dt*maxt/dr);


	bloch_kx=2*(typ_prec)M_PI*bloch_kx;
	bloch_ky=2*(typ_prec)M_PI*bloch_ky;
	bloch_kz=2*(typ_prec)M_PI*bloch_kz;
	
	#ifndef ACCURATE_BLOCH	
	        bloch_expikx=(typ_pola)exp(bloch_kx*sqrt((typ_pola)-1));
	        bloch_expiky=(typ_pola)exp(bloch_ky*sqrt((typ_pola)-1));
	        bloch_expikz=(typ_pola)exp(bloch_kz*sqrt((typ_pola)-1));
	
	        bloch_expmikx=(typ_pola)exp(-bloch_kx*sqrt((typ_pola)-1));
	        bloch_expmiky=(typ_pola)exp(-bloch_ky*sqrt((typ_pola)-1));
	        bloch_expmikz=(typ_pola)exp(-bloch_kz*sqrt((typ_pola)-1));
    #endif

    // simple check of parameters
    { // szybkie sprawdzenie czy nie ma zapisu nie istniejacej skladowej - o jeden crash mniej
      bool check_table[]={is_calc_ex,is_calc_ey,is_calc_ez,is_calc_hx,is_calc_hy,is_calc_hz};
      for(i=0;i<n_outputs;i++) for(j=0;j<6;j++) if (!check_table[j] && skladowe[i][j]) skladowe[i][j]=false;
    }

	return true;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
//                                                Main                                                      //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int n_args,char** args)
{
	// sprawdzamy jak program zostal wywolany

	for (int n=4;n<=7;n++)
	{
		if (n_args<=n) break;
		{
			std::string tmpstr(args[n]);
			if (tmpstr=="ressurect") ressurect=true;
			if (tmpstr=="noreport") noreport=true;
			if (tmpstr=="silent") silentmode=true;
			if (tmpstr=="iosilent") iosilentmode=true;
		}

	}

	if (!silentmode)
    { 
      std::cout<<"\nEMFIDES engine v "<<VERSION<<"         coded by W.M.Saj, compiled: "<<__DATE__<<" \n";
      #ifdef DOUBLE_PRECISION
      std::cout<<"double";
      #else
      std::cout<<"single";
      #endif
      std::cout<<" precision numbers used for ";
      #ifdef COMPLEX_FIELD
      std::cout<<"complex";
      #else
      std::cout<<"real";
      #endif
      std::cout<<" field computations";
      #ifdef ACCURATE_BLOCH	
      std::cout<<" with improved accuracy Bloch condition"; 
      #endif
      std::cout<<"\n\n";
    }
	// jesli zle wywolany to wyrzuc informacje o uzyciu
	if ((n_args!=4) && (silentmode==false) && (ressurect==false) && (noreport==false) && (iosilentmode==false))
	{
        std::cout<<"usage : emfides"<<(int)(1000*VERSION)<<".exe  parameters_file structure_file output_files_prefix [ressurect] [noreport] [silent] [iosilent]\n";
        return 1;
    }

    // wczytywanie pliku z parametrami
	if (!silentmode) std::cout<<" *loading parameters...\n";
	if (!read_params(args[1])) {deinit_other_arrays();std::cerr<<"Error: Parameters file reading failed\n"; return 1;}
	// wczytywanie struktury
	if (!silentmode) std::cout<<" *loading structure...\n";
	if (!init_structure_arrays()) {deinit_other_arrays();std::cerr<<"Error: Memory allocation for structure reading failed\n"; return 1;}
	if (!create_structure(args[2])) {deinit_structure_arrays();deinit_arrays(); deinit_other_arrays();std::cerr<<"Error: Structure file reading failed\n"; return 1;}

	// rezerwacja pamieci pod tablice
	if (!silentmode) std::cout<<" *memory allocation...\n";

    if (!init_arrays()) {deinit_structure_arrays();deinit_arrays();deinit_other_arrays();std::cerr<<"Error: Memory allocation failed\n"; return 1;}
    // odczytujemy czas startu ( do raportu)
    strtime=time(0);prevtime=strtime;
    // ustawiamy wartosc roznych pomocniczych stalych
    fill_ABC();
    create_cst(); 
    calc_trig();

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if (ressurect==true)
	{
		if (!silentmode) std::cout<<" *ressurection...\n";
		if (!load_all(args[3])) {deinit_arrays(); deinit_other_arrays();perror("Ressurection failed:"); return 1;} ;
	}


	/////////////////////////////////////////////////////////////////////
	// konstrukcja symulacji
	
	// ustalenie funkcji zrzutu
	cenSx= ((!is_calc_ey ||!is_calc_hz)&&(!is_calc_ez ||!is_calc_hy)) ? return_zer:((!is_calc_ey ||!is_calc_hz) ? cenSx_null_Ey_Hz : ((!is_calc_ez ||!is_calc_hy) ? cenSx_null_Ez_Hy :cenSx_all));
	cenSy= ((!is_calc_ex ||!is_calc_hz)&&(!is_calc_ez ||!is_calc_hx)) ? return_zer:((!is_calc_ex ||!is_calc_hz) ? cenSy_null_Ex_Hz : ((!is_calc_ez ||!is_calc_hx) ? cenSy_null_Ez_Hx :cenSy_all)); 
	cenSz= ((!is_calc_ey ||!is_calc_hx)&&(!is_calc_ex ||!is_calc_hy)) ? return_zer:((!is_calc_ey ||!is_calc_hx) ? cenSz_null_Ey_Hx : ((!is_calc_ex ||!is_calc_hy) ? cenSz_null_Ex_Hy :cenSz_all));
	
	// filling cenfun 
	cenfun[0]=cenEx;cenfun[1]=cenEy;cenfun[2]=cenEz;
	cenfun[3]=cenHx;cenfun[4]=cenHy;cenfun[5]=cenHz;
	cenfun[6]=cenSx;cenfun[7]=cenSy;cenfun[8]=cenSz;
	
	(wx==1)  ? set_periodic_lefrig() : ((wx==2) ? set_symmetric_lefrig():( (wx==3) ? set_bloch_lefrig():please_do_nuthing()));
	(wy==1)  ? set_periodic_bottop() : ((wy==2) ? set_symmetric_bottop():( (wy==3) ? set_bloch_bottop():please_do_nuthing()));
	(wz==1)  ? set_periodic_bacfro() : ((wz==2) ? set_symmetric_bacfro():( (wz==3) ? set_bloch_bacfro():please_do_nuthing()));

    for(int s=0;s<n_SRCS;s++)
	{
    if (source_type[s]==0)
	{
		inc_Ex[s]= (is_calc_ex ? plane_inc_Ex:return_zero);
		inc_Ey[s]= (is_calc_ey ? plane_inc_Ey:return_zero);
		inc_Ez[s]= (is_calc_ez ? plane_inc_Ez:return_zero);

		inc_Hx[s]= (is_calc_hx ? plane_inc_Hx:return_zero);
		inc_Hy[s]= (is_calc_hy ? plane_inc_Hy:return_zero);
		inc_Hz[s]= (is_calc_hz ? plane_inc_Hz:return_zero);	
	}
	else if (source_type[s]==1)
	{
		inc_Ex[s]= (is_calc_ex ? beam_inc_Ex:return_zero);
		inc_Ey[s]= (is_calc_ey ? beam_inc_Ey:return_zero);
		inc_Ez[s]= (is_calc_ez ? beam_inc_Ez:return_zero);

		inc_Hx[s]= (is_calc_hx ? beam_inc_Hx:return_zero);
		inc_Hy[s]= (is_calc_hy ? beam_inc_Hy:return_zero);
		inc_Hz[s]= (is_calc_hz ? beam_inc_Hz:return_zero);
	}
	else if (source_type[s]==2)
	{
		inc_Ex[s]= (is_calc_ex ? obeam_inc_Ex:return_zero);
		inc_Ey[s]= (is_calc_ey ? obeam_inc_Ey:return_zero);
		inc_Ez[s]= (is_calc_ez ? obeam_inc_Ez:return_zero);

		inc_Hx[s]= (is_calc_hx ? obeam_inc_Hx:return_zero);
		inc_Hy[s]= (is_calc_hy ? obeam_inc_Hy:return_zero);
		inc_Hz[s]= (is_calc_hz ? obeam_inc_Hz:return_zero);
	}
    }
	if (!silentmode) {stage_0=print_time;} else {stage_0=please_do_nuthing;}

	// ustawianie update'u srodka
	if (mode3D) 
	{
		stage_E1x=(!is_calc_ex || (!is_calc_hy && !is_calc_hz)) ? please_do_nuthing:(is_calc_hy ? (is_calc_hz ? update_inter_E3x:update_inter_E3x_null_Hz):update_inter_E3x_null_Hy);
		stage_E1y=(!is_calc_ey || (!is_calc_hx && !is_calc_hz)) ? please_do_nuthing:(is_calc_hx ? (is_calc_hz ? update_inter_E3y:update_inter_E3y_null_Hz):update_inter_E3y_null_Hx);
		stage_E1z=(!is_calc_ez || (!is_calc_hy && !is_calc_hx)) ? please_do_nuthing:(is_calc_hy ? (is_calc_hx ? update_inter_E3z:update_inter_E3z_null_Hx):update_inter_E3z_null_Hy);

		stage_H1x=(!is_calc_hx || (!is_calc_ey && !is_calc_ez)) ? please_do_nuthing:(is_calc_ey ? (is_calc_ez ? update_inter_H3x:update_inter_H3x_null_Ez):update_inter_H3x_null_Ey);
		stage_H1y=(!is_calc_hy || (!is_calc_ex && !is_calc_ez)) ? please_do_nuthing:(is_calc_ex ? (is_calc_ez ? update_inter_H3y:update_inter_H3y_null_Ez):update_inter_H3y_null_Ex);
		stage_H1z=(!is_calc_hz || (!is_calc_ey && !is_calc_ex)) ? please_do_nuthing:(is_calc_ey ? (is_calc_ex ? update_inter_H3z:update_inter_H3z_null_Ex):update_inter_H3z_null_Ey);


	} 
	else 
	{
        stage_E1x=(!is_calc_ex || (!is_calc_hy && !is_calc_hz)) ? please_do_nuthing:(is_calc_hy ? (is_calc_hz ? update_inter_Ex:update_inter_Ex_null_Hz):update_inter_Ex_null_Hy);
		stage_E1y=(!is_calc_ey || (!is_calc_hx && !is_calc_hz)) ? please_do_nuthing:(is_calc_hx ? (is_calc_hz ? update_inter_Ey:update_inter_Ey_null_Hz):update_inter_Ey_null_Hx);
		stage_E1z=(!is_calc_ez || (!is_calc_hy && !is_calc_hx)) ? please_do_nuthing:(is_calc_hy ? (is_calc_hx ? update_inter_Ez:update_inter_Ez_null_Hx):update_inter_Ez_null_Hy);

		stage_H1x=(!is_calc_hx || (!is_calc_ey && !is_calc_ez)) ? please_do_nuthing:(is_calc_ey ? (is_calc_ez ? update_inter_Hx:update_inter_Hx_null_Ez):update_inter_Hx_null_Ey);
		stage_H1y=(!is_calc_hy || (!is_calc_ex && !is_calc_ez)) ? please_do_nuthing:(is_calc_ex ? (is_calc_ez ? update_inter_Hy:update_inter_Hy_null_Ez):update_inter_Hy_null_Ex);
		stage_H1z=(!is_calc_hz || (!is_calc_ey && !is_calc_ex)) ? please_do_nuthing:(is_calc_ey ? (is_calc_ex ? update_inter_Hz:update_inter_Hz_null_Ex):update_inter_Hz_null_Ey);
	}

	// budowanie symulacji
	// UPML ??
	
	stage_E2x_1=((wx==1 || wx==3 || !is_calc_ex || (!is_calc_hz && !is_calc_hy)) ? please_do_nuthing : (is_calc_hz ? (is_calc_hy ? update_lef_Ex:update_lef_Ex_null_Hy):update_lef_Ex_null_Hz));
	stage_E3x_1=((wy==1 || wy==3 || !is_calc_ex || (!is_calc_hz && !is_calc_hy)) ? please_do_nuthing : (is_calc_hz ? (is_calc_hy ? update_bot_Ex:update_bot_Ex_null_Hy):update_bot_Ex_null_Hz));
	stage_E4x_1=((wz==1 || wz==3 || !is_calc_ex || (!is_calc_hz && !is_calc_hy)) ? please_do_nuthing : (is_calc_hz ? (is_calc_hy ? update_bac_Ex:update_bac_Ex_null_Hy):update_bac_Ex_null_Hz));

	stage_E2y_1=((wx==1 || wx==3 || !is_calc_ey || (!is_calc_hz && !is_calc_hx)) ? please_do_nuthing : (is_calc_hz ? (is_calc_hx ? update_lef_Ey:update_lef_Ey_null_Hx):update_lef_Ey_null_Hz));
	stage_E3y_1=((wy==1 || wy==3 || !is_calc_ey || (!is_calc_hz && !is_calc_hx)) ? please_do_nuthing : (is_calc_hz ? (is_calc_hx ? update_bot_Ey:update_bot_Ey_null_Hx):update_bot_Ey_null_Hz));
	stage_E4y_1=((wz==1 || wz==3 || !is_calc_ey || (!is_calc_hz && !is_calc_hx)) ? please_do_nuthing : (is_calc_hz ? (is_calc_hx ? update_bac_Ey:update_bac_Ey_null_Hx):update_bac_Ey_null_Hz));

	stage_E2z_1=((wx==1 || wx==3 || !is_calc_ez || (!is_calc_hy && !is_calc_hx)) ? please_do_nuthing : (is_calc_hx ? (is_calc_hy ? update_lef_Ez:update_lef_Ez_null_Hy):update_lef_Ez_null_Hx));
	stage_E3z_1=((wy==1 || wy==3 || !is_calc_ez || (!is_calc_hy && !is_calc_hx)) ? please_do_nuthing : (is_calc_hx ? (is_calc_hy ? update_bot_Ez:update_bot_Ez_null_Hy):update_bot_Ez_null_Hx));
	stage_E4z_1=((wz==1 || wz==3 || !is_calc_ez || (!is_calc_hy && !is_calc_hx)) ? please_do_nuthing : (is_calc_hx ? (is_calc_hy ? update_bac_Ez:update_bac_Ez_null_Hy):update_bac_Ez_null_Hx));

	stage_E2x_2=((wx==1 || wx==2 || wx==3 || !is_calc_ex || (!is_calc_hz && !is_calc_hy)) ? please_do_nuthing : (is_calc_hz ? (is_calc_hy ? update_rig_Ex:update_rig_Ex_null_Hy):update_rig_Ex_null_Hz));
	stage_E3x_2=((wy==1 || wy==2 || wy==3 || !is_calc_ex || (!is_calc_hz && !is_calc_hy)) ? please_do_nuthing : (is_calc_hz ? (is_calc_hy ? update_top_Ex:update_top_Ex_null_Hy):update_top_Ex_null_Hz));
	stage_E4x_2=((wz==1 || wz==2 || wz==3 || !is_calc_ex || (!is_calc_hz && !is_calc_hy)) ? please_do_nuthing : (is_calc_hz ? (is_calc_hy ? update_fro_Ex:update_fro_Ex_null_Hy):update_fro_Ex_null_Hz));

	stage_E2y_2=((wx==1 || wx==2 || wx==3 || !is_calc_ey || (!is_calc_hz && !is_calc_hx)) ? please_do_nuthing : (is_calc_hz ? (is_calc_hx ? update_rig_Ey:update_rig_Ey_null_Hx):update_rig_Ey_null_Hz));
	stage_E3y_2=((wy==1 || wy==2 || wy==3 || !is_calc_ey || (!is_calc_hz && !is_calc_hx)) ? please_do_nuthing : (is_calc_hz ? (is_calc_hx ? update_top_Ey:update_top_Ey_null_Hx):update_top_Ey_null_Hz));
	stage_E4y_2=((wz==1 || wz==2 || wz==3 || !is_calc_ey || (!is_calc_hz && !is_calc_hx)) ? please_do_nuthing : (is_calc_hz ? (is_calc_hx ? update_fro_Ey:update_fro_Ey_null_Hx):update_fro_Ey_null_Hz));

	stage_E2z_2=((wx==1 || wx==2  || wx==3 || !is_calc_ez || (!is_calc_hy && !is_calc_hx)) ? please_do_nuthing : (is_calc_hx ? (is_calc_hy ? update_rig_Ez:update_rig_Ez_null_Hy):update_rig_Ez_null_Hx));
	stage_E3z_2=((wy==1 || wy==2  || wy==3 || !is_calc_ez || (!is_calc_hy && !is_calc_hx)) ? please_do_nuthing : (is_calc_hx ? (is_calc_hy ? update_top_Ez:update_top_Ez_null_Hy):update_top_Ez_null_Hx));
	stage_E4z_2=((wz==1 || wz==2  || wz==3 || !is_calc_ez || (!is_calc_hy && !is_calc_hx)) ? please_do_nuthing : (is_calc_hx ? (is_calc_hy ? update_fro_Ez:update_fro_Ez_null_Hy):update_fro_Ez_null_Hx));


	stage_H2x_1=((wx==1 || wx==3 || !is_calc_hx || (!is_calc_ez && !is_calc_ey)) ? please_do_nuthing : (is_calc_ez ? (is_calc_ey ? update_lef_Hx:update_lef_Hx_null_Ey):update_lef_Hx_null_Ez));
	stage_H3x_1=((wy==1 || wy==3 || !is_calc_hx || (!is_calc_ez && !is_calc_ey)) ? please_do_nuthing : (is_calc_ez ? (is_calc_ey ? update_bot_Hx:update_bot_Hx_null_Ey):update_bot_Hx_null_Ez));
	stage_H4x_1=((wz==1 || wz==3 || !is_calc_hx || (!is_calc_ez && !is_calc_ey)) ? please_do_nuthing : (is_calc_ez ? (is_calc_ey ? update_bac_Hx:update_bac_Hx_null_Ey):update_bac_Hx_null_Ez));

	stage_H2y_1=((wx==1 || wx==3 || !is_calc_hy || (!is_calc_ez && !is_calc_ex)) ? please_do_nuthing : (is_calc_ez ? (is_calc_ex ? update_lef_Hy:update_lef_Hy_null_Ex):update_lef_Hy_null_Ez));
	stage_H3y_1=((wy==1 || wy==3 || !is_calc_hy || (!is_calc_ez && !is_calc_ex)) ? please_do_nuthing : (is_calc_ez ? (is_calc_ex ? update_bot_Hy:update_bot_Hy_null_Ex):update_bot_Hy_null_Ez));
	stage_H4y_1=((wz==1 || wz==3 || !is_calc_hy || (!is_calc_ez && !is_calc_ex)) ? please_do_nuthing : (is_calc_ez ? (is_calc_ex ? update_bac_Hy:update_bac_Hy_null_Ex):update_bac_Hy_null_Ez));

	stage_H2z_1=((wx==1 || wx==3 || !is_calc_hz || (!is_calc_ey && !is_calc_ex)) ? please_do_nuthing : (is_calc_ex ? (is_calc_ey ? update_lef_Hz:update_lef_Hz_null_Ey):update_lef_Hz_null_Ex));
	stage_H3z_1=((wy==1 || wy==3 || !is_calc_hz || (!is_calc_ey && !is_calc_ex)) ? please_do_nuthing : (is_calc_ex ? (is_calc_ey ? update_bot_Hz:update_bot_Hz_null_Ey):update_bot_Hz_null_Ex));
	stage_H4z_1=((wz==1 || wz==3 || !is_calc_hz || (!is_calc_ey && !is_calc_ex)) ? please_do_nuthing : (is_calc_ex ? (is_calc_ey ? update_bac_Hz:update_bac_Hz_null_Ey):update_bac_Hz_null_Ex));

	stage_H2x_2=((wx==1 || wx==2 || wx==3 || !is_calc_hx || (!is_calc_ez && !is_calc_ey)) ? please_do_nuthing : (is_calc_ez ? (is_calc_ey ? update_rig_Hx:update_rig_Hx_null_Ey):update_rig_Hx_null_Ez));
	stage_H3x_2=((wy==1 || wy==2 || wy==3 || !is_calc_hx || (!is_calc_ez && !is_calc_ey)) ? please_do_nuthing : (is_calc_ez ? (is_calc_ey ? update_top_Hx:update_top_Hx_null_Ey):update_top_Hx_null_Ez));
	stage_H4x_2=((wz==1 || wz==2 || wz==3 || !is_calc_hx || (!is_calc_ez && !is_calc_ey)) ? please_do_nuthing : (is_calc_ez ? (is_calc_ey ? update_fro_Hx:update_fro_Hx_null_Ey):update_fro_Hx_null_Ez));

	stage_H2y_2=((wx==1 || wx==2 || wx==3 || !is_calc_hy || (!is_calc_ez && !is_calc_ex)) ? please_do_nuthing : (is_calc_ez ? (is_calc_ex ? update_rig_Hy:update_rig_Hy_null_Ex):update_rig_Hy_null_Ez));
	stage_H3y_2=((wy==1 || wy==2 || wy==3 || !is_calc_hy || (!is_calc_ez && !is_calc_ex)) ? please_do_nuthing : (is_calc_ez ? (is_calc_ex ? update_top_Hy:update_top_Hy_null_Ex):update_top_Hy_null_Ez));
	stage_H4y_2=((wz==1 || wz==2 || wz==3 || !is_calc_hy || (!is_calc_ez && !is_calc_ex)) ? please_do_nuthing : (is_calc_ez ? (is_calc_ex ? update_fro_Hy:update_fro_Hy_null_Ex):update_fro_Hy_null_Ez));

	stage_H2z_2=((wx==1 || wx==2  || wx==3 || !is_calc_hz || (!is_calc_ey && !is_calc_ex)) ? please_do_nuthing : (is_calc_ex ? (is_calc_ey ? update_rig_Hz:update_rig_Hz_null_Ey):update_rig_Hz_null_Ex));
	stage_H3z_2=((wy==1 || wy==2  || wy==3 || !is_calc_hz || (!is_calc_ey && !is_calc_ex)) ? please_do_nuthing : (is_calc_ex ? (is_calc_ey ? update_top_Hz:update_top_Hz_null_Ey):update_top_Hz_null_Ex));
	stage_H4z_2=((wz==1 || wz==2  || wz==3 || !is_calc_hz || (!is_calc_ey && !is_calc_ex)) ? please_do_nuthing : (is_calc_ex ? (is_calc_ey ? update_fro_Hz:update_fro_Hz_null_Ey):update_fro_Hz_null_Ex));



	//

	// zrodlo ??
	for (int s=0;s<=n_SRCS;s++)
    {
	if (source_type[s]<=2)
	{
 		stage_E5_1_1[s]=(((wx==1 || wx==3) || (!is_calc_ey || !is_calc_hz)) ? please_do_nothing : TF_SF_lef_E_Ey_Hz);
		stage_E6_1_1[s]=(((wy==1 || wy==3) || (!is_calc_ex || !is_calc_hz)) ? please_do_nothing : TF_SF_bot_E_Ex_Hz);
		stage_E7_1_1[s]=(((wz==1 || wz==3) || (!is_calc_ex || !is_calc_hy)) ? please_do_nothing : TF_SF_bac_E_Ex_Hy);
		
	    stage_E5_1_2[s]=(((wx==1 || wx==3) || (!is_calc_ez || !is_calc_hy)) ? please_do_nothing : TF_SF_lef_E_Ez_Hy);
		stage_E6_1_2[s]=(((wy==1 || wy==3) || (!is_calc_ez || !is_calc_hx)) ? please_do_nothing : TF_SF_bot_E_Ez_Hx);
		stage_E7_1_2[s]=(((wz==1 || wz==3) || (!is_calc_ey || !is_calc_hx)) ? please_do_nothing : TF_SF_bac_E_Ey_Hx);

 		stage_E5_2_1[s]=(((wx==1 || wx==2 || wx==3) || (!is_calc_ey || !is_calc_hz)) ? please_do_nothing : TF_SF_rig_E_Ey_Hz);
		stage_E6_2_1[s]=(((wy==1 || wy==2 || wy==3) || (!is_calc_ex || !is_calc_hz)) ? please_do_nothing : TF_SF_top_E_Ex_Hz);
		stage_E7_2_1[s]=(((wz==1 || wz==2 || wz==3) || (!is_calc_ex || !is_calc_hy)) ? please_do_nothing : TF_SF_fro_E_Ex_Hy);
        
        stage_E5_2_2[s]=(((wx==1 || wx==2 || wx==3) || (!is_calc_ez || !is_calc_hy)) ? please_do_nothing : TF_SF_rig_E_Ez_Hy);
		stage_E6_2_2[s]=(((wy==1 || wy==2 || wy==3) || (!is_calc_ez || !is_calc_hx)) ? please_do_nothing : TF_SF_top_E_Ez_Hx);
		stage_E7_2_2[s]=(((wz==1 || wz==2 || wz==3) || (!is_calc_ey || !is_calc_hx)) ? please_do_nothing : TF_SF_fro_E_Ey_Hx);	
	
    }
	else
	{
		stage_E5_1_1[s]=((source_type[s]==3 && is_calc_ex) ? electric_dipolx : ((source_type[s]==4 && is_calc_ey) ? electric_dipoly :((source_type[s]==5 && is_calc_ez) ? electric_dipolz:please_do_nothing)));
		stage_E6_1_1[s]=please_do_nothing;
		stage_E7_1_1[s]=please_do_nothing;
		
		stage_E5_1_2[s]=please_do_nothing;
		stage_E6_1_2[s]=please_do_nothing;
		stage_E7_1_2[s]=please_do_nothing;
		
		stage_E5_2_1[s]=please_do_nothing;
		stage_E6_2_1[s]=please_do_nothing;
		stage_E7_2_1[s]=please_do_nothing;
		
		stage_E5_2_2[s]=please_do_nothing;
		stage_E6_2_2[s]=please_do_nothing;
		stage_E7_2_2[s]=please_do_nothing;		
	}
	
	if (source_type[s]<=2)
	{
 		stage_H5_1_1[s]=(((wx==1 || wx==3) || (!is_calc_hy || !is_calc_ez)) ? please_do_nothing : TF_SF_lef_H_Hy_Ez);
		stage_H6_1_1[s]=(((wy==1 || wy==3) || (!is_calc_hx || !is_calc_ez)) ? please_do_nothing : TF_SF_bot_H_Hx_Ez);
		stage_H7_1_1[s]=(((wz==1 || wz==3) || (!is_calc_hx || !is_calc_ey)) ? please_do_nothing : TF_SF_bac_H_Hx_Ey);
		
	    stage_H5_1_2[s]=(((wx==1 || wx==3) || (!is_calc_hz || !is_calc_ey)) ? please_do_nothing : TF_SF_lef_H_Hz_Ey);
		stage_H6_1_2[s]=(((wy==1 || wy==3) || (!is_calc_hz || !is_calc_ex)) ? please_do_nothing : TF_SF_bot_H_Hz_Ex);
		stage_H7_1_2[s]=(((wz==1 || wz==3) || (!is_calc_hy || !is_calc_ex)) ? please_do_nothing : TF_SF_bac_H_Hy_Ex);

 		stage_H5_2_1[s]=(((wx==1 || wx==2 || wx==3) || (!is_calc_hy || !is_calc_ez)) ? please_do_nothing : TF_SF_rig_H_Hy_Ez);
		stage_H6_2_1[s]=(((wy==1 || wy==2 || wy==3) || (!is_calc_hx || !is_calc_ez)) ? please_do_nothing : TF_SF_top_H_Hx_Ez);
		stage_H7_2_1[s]=(((wz==1 || wz==2 || wz==3) || (!is_calc_hx || !is_calc_ey)) ? please_do_nothing : TF_SF_fro_H_Hx_Ey);
        
        stage_H5_2_2[s]=(((wx==1 || wx==2 || wx==3) || (!is_calc_hz || !is_calc_ey)) ? please_do_nothing : TF_SF_rig_H_Hz_Ey);
		stage_H6_2_2[s]=(((wy==1 || wy==2 || wy==3) || (!is_calc_hz || !is_calc_ex)) ? please_do_nothing : TF_SF_top_H_Hz_Ex);
		stage_H7_2_2[s]=(((wz==1 || wz==2 || wz==3) || (!is_calc_hy || !is_calc_ex)) ? please_do_nothing : TF_SF_fro_H_Hy_Ex);	
                       
	}
	else
	{
 		stage_H5_1_1[s]=((source_type[s]==6 && is_calc_hx) ? magnetic_dipolx : ((source_type[s]==7 && is_calc_hy) ? magnetic_dipoly :((source_type[s]==8 && is_calc_hz) ? magnetic_dipolz:please_do_nothing)));
		stage_H6_1_1[s]=please_do_nothing;
		stage_H7_1_1[s]=please_do_nothing;
		
		stage_H5_1_2[s]=please_do_nothing;
		stage_H6_1_2[s]=please_do_nothing;
		stage_H7_1_2[s]=please_do_nothing;
		
		stage_H5_2_1[s]=please_do_nothing;
		stage_H6_2_1[s]=please_do_nothing;
		stage_H7_2_1[s]=please_do_nothing;
		
		stage_H5_2_2[s]=please_do_nothing;
		stage_H6_2_2[s]=please_do_nothing;
		stage_H7_2_2[s]=please_do_nothing;
	}
}

	// drude,lorentz,debeye
	if (mode3D)
	{
		stage_E8x = (n_drd_Ex>0 && is_calc_ex) ? update_drdE3x : please_do_nuthing;
		stage_E8y = (n_drd_Ey>0 && is_calc_ey) ? update_drdE3y : please_do_nuthing;
		stage_E8z = (n_drd_Ez>0 && is_calc_ez) ? update_drdE3z : please_do_nuthing;

		stage_H8x = (n_drd_Hx>0 && is_calc_hx) ? update_drdH3x : please_do_nuthing;
		stage_H8y = (n_drd_Hy>0 && is_calc_hy) ? update_drdH3y : please_do_nuthing;
		stage_H8z = (n_drd_Hz>0 && is_calc_hz) ? update_drdH3z : please_do_nuthing;

		stage_E9x = (n_lor_Ex>0 && is_calc_ex) ? update_lorE3x : please_do_nuthing;
		stage_E9y = (n_lor_Ey>0 && is_calc_ey) ? update_lorE3y : please_do_nuthing;
		stage_E9z = (n_lor_Ez>0 && is_calc_ez) ? update_lorE3z : please_do_nuthing;

		stage_H9x = (n_lor_Hx>0 && is_calc_hx) ? update_lorH3x : please_do_nuthing;
		stage_H9y = (n_lor_Hy>0 && is_calc_hy) ? update_lorH3y : please_do_nuthing;
		stage_H9z = (n_lor_Hz>0 && is_calc_hz) ? update_lorH3z : please_do_nuthing;

		stage_E10x = (n_dby_Ex>0 && is_calc_ex) ? update_dbyE3x : please_do_nuthing;
		stage_E10y = (n_dby_Ey>0 && is_calc_ey) ? update_dbyE3y : please_do_nuthing;
		stage_E10z = (n_dby_Ey>0 && is_calc_ez) ? update_dbyE3z : please_do_nuthing;

		stage_H10x = (n_dby_Hx>0 && is_calc_hx) ? update_dbyH3x : please_do_nuthing;
		stage_H10y = (n_dby_Hy>0 && is_calc_hy) ? update_dbyH3y : please_do_nuthing;
		stage_H10z = (n_dby_Hz>0 && is_calc_hz) ? update_dbyH3z : please_do_nuthing;
		
		
		stage_J8x = (n_drd_Ex>0 && is_calc_ex) ? update_drdJ3x : please_do_nuthing;
		stage_J8y = (n_drd_Ey>0 && is_calc_ey) ? update_drdJ3y : please_do_nuthing;
		stage_J8z = (n_drd_Ez>0 && is_calc_ez) ? update_drdJ3z : please_do_nuthing;

		stage_Jh8x = (n_drd_Hx>0 && is_calc_hx) ? update_drdJh3x : please_do_nuthing;
		stage_Jh8y = (n_drd_Hy>0 && is_calc_hy) ? update_drdJh3y : please_do_nuthing;
		stage_Jh8z = (n_drd_Hz>0 && is_calc_hz) ? update_drdJh3z : please_do_nuthing;

		stage_J9x = (n_lor_Ex>0 && is_calc_ex) ? update_lorJ3x : please_do_nuthing;
		stage_J9y = (n_lor_Ey>0 && is_calc_ey) ? update_lorJ3y : please_do_nuthing;
		stage_J9z = (n_lor_Ez>0 && is_calc_ez) ? update_lorJ3z : please_do_nuthing;

		stage_Jh9x = (n_lor_Hx>0 && is_calc_hx) ? update_lorJh3x : please_do_nuthing;
		stage_Jh9y = (n_lor_Hy>0 && is_calc_hy) ? update_lorJh3y : please_do_nuthing;
		stage_Jh9z = (n_lor_Hz>0 && is_calc_hz) ? update_lorJh3z : please_do_nuthing;

		stage_J10x = (n_dby_Ex>0 && is_calc_ex) ? update_dbyJ3x : please_do_nuthing;
		stage_J10y = (n_dby_Ey>0 && is_calc_ey) ? update_dbyJ3y : please_do_nuthing;
		stage_J10z = (n_dby_Ey>0 && is_calc_ez) ? update_dbyJ3z : please_do_nuthing;

		stage_Jh10x = (n_dby_Hx>0 && is_calc_hx) ? update_dbyJh3x : please_do_nuthing;
		stage_Jh10y = (n_dby_Hy>0 && is_calc_hy) ? update_dbyJh3y : please_do_nuthing;
		stage_Jh10z = (n_dby_Hz>0 && is_calc_hz) ? update_dbyJh3z : please_do_nuthing;
		
		
	}
	else
	{
		stage_E8x = (n_drd_Ex>0 && is_calc_ex) ? update_drdEx : please_do_nuthing;
		stage_E8y = (n_drd_Ey>0 && is_calc_ey) ? update_drdEy : please_do_nuthing;
		stage_E8z = (n_drd_Ez>0 && is_calc_ez) ? update_drdEz : please_do_nuthing;

		stage_H8x = (n_drd_Hx>0 && is_calc_hx) ? update_drdHx : please_do_nuthing;
		stage_H8y = (n_drd_Hy>0 && is_calc_hy) ? update_drdHy : please_do_nuthing;
		stage_H8z = (n_drd_Hz>0 && is_calc_hz) ? update_drdHz : please_do_nuthing;

		stage_E9x = (n_lor_Ex>0 && is_calc_ex) ? update_lorEx : please_do_nuthing;
		stage_E9y = (n_lor_Ey>0 && is_calc_ey) ? update_lorEy : please_do_nuthing;
		stage_E9z = (n_lor_Ez>0 && is_calc_ez) ? update_lorEz : please_do_nuthing;

		stage_H9x = (n_lor_Hx>0 && is_calc_hx) ? update_lorHx : please_do_nuthing;
		stage_H9y = (n_lor_Hy>0 && is_calc_hy) ? update_lorHy : please_do_nuthing;
		stage_H9z = (n_lor_Hz>0 && is_calc_hz) ? update_lorHz : please_do_nuthing;

		stage_E10x = (n_dby_Ex>0 && is_calc_ex) ? update_dbyEx : please_do_nuthing;
		stage_E10y = (n_dby_Ey>0 && is_calc_ey) ? update_dbyEy : please_do_nuthing;
		stage_E10z = (n_dby_Ey>0 && is_calc_ez) ? update_dbyEz : please_do_nuthing;

		stage_H10x = (n_dby_Hx>0 && is_calc_hx) ? update_dbyHx : please_do_nuthing;
		stage_H10y = (n_dby_Hy>0 && is_calc_hy) ? update_dbyHy : please_do_nuthing;
		stage_H10z = (n_dby_Hz>0 && is_calc_hz) ? update_dbyHz : please_do_nuthing;
		
		stage_J8x = (n_drd_Ex>0 && is_calc_ex) ? update_drdJx : please_do_nuthing;
		stage_J8y = (n_drd_Ey>0 && is_calc_ey) ? update_drdJy : please_do_nuthing;
		stage_J8z = (n_drd_Ez>0 && is_calc_ez) ? update_drdJz : please_do_nuthing;

		stage_Jh8x = (n_drd_Hx>0 && is_calc_hx) ? update_drdJhx : please_do_nuthing;
		stage_Jh8y = (n_drd_Hy>0 && is_calc_hy) ? update_drdJhy : please_do_nuthing;
		stage_Jh8z = (n_drd_Hz>0 && is_calc_hz) ? update_drdJhz : please_do_nuthing;

		stage_J9x = (n_lor_Ex>0 && is_calc_ex) ? update_lorJx : please_do_nuthing;
		stage_J9y = (n_lor_Ey>0 && is_calc_ey) ? update_lorJy : please_do_nuthing;
		stage_J9z = (n_lor_Ez>0 && is_calc_ez) ? update_lorJz : please_do_nuthing;

		stage_Jh9x = (n_lor_Hx>0 && is_calc_hx) ? update_lorJhx : please_do_nuthing;
		stage_Jh9y = (n_lor_Hy>0 && is_calc_hy) ? update_lorJhy : please_do_nuthing;
		stage_Jh9z = (n_lor_Hz>0 && is_calc_hz) ? update_lorJhz : please_do_nuthing;

		stage_J10x = (n_dby_Ex>0 && is_calc_ex) ? update_dbyJx : please_do_nuthing;
		stage_J10y = (n_dby_Ey>0 && is_calc_ey) ? update_dbyJy : please_do_nuthing;
		stage_J10z = (n_dby_Ey>0 && is_calc_ez) ? update_dbyJz : please_do_nuthing;

		stage_Jh10x = (n_dby_Hx>0 && is_calc_hx) ? update_dbyJhx : please_do_nuthing;
		stage_Jh10y = (n_dby_Hy>0 && is_calc_hy) ? update_dbyJhy : please_do_nuthing;
		stage_Jh10z = (n_dby_Hz>0 && is_calc_hz) ? update_dbyJhz : please_do_nuthing;		

	}

	// war brzegowe periodyczne lub symetryczne
	stage_E11x=((wx==0 || !is_calc_ex) ? please_do_nuthing : ((wx==1) ? periodic_lefrig_Ex : ((wx==3) ? bloch_lefrig_Ex : symmetric_lefrig_Ex)));
	stage_E11y=((wx==0 || !is_calc_ey) ? please_do_nuthing : ((wx==1) ? periodic_lefrig_Ey : ((wx==3) ? bloch_lefrig_Ey : symmetric_lefrig_Ey)));
	stage_E11z=((wx==0 || !is_calc_ez) ? please_do_nuthing : ((wx==1) ? periodic_lefrig_Ez : ((wx==3) ? bloch_lefrig_Ez : symmetric_lefrig_Ez)));
	
	stage_E12x=((wy==0 || !is_calc_ex) ? please_do_nuthing : ((wy==1) ? periodic_bottop_Ex : ((wy==3) ? bloch_bottop_Ex : symmetric_bottop_Ex)));
	stage_E12y=((wy==0 || !is_calc_ey) ? please_do_nuthing : ((wy==1) ? periodic_bottop_Ey : ((wy==3) ? bloch_bottop_Ey : symmetric_bottop_Ey)));
	stage_E12z=((wy==0 || !is_calc_ez) ? please_do_nuthing : ((wy==1) ? periodic_bottop_Ez : ((wy==3) ? bloch_bottop_Ez : symmetric_bottop_Ez)));
	
	stage_E13x=((wz==0 || !is_calc_ex) ? please_do_nuthing : ((wz==1) ? periodic_bacfro_Ex : ((wz==3) ? bloch_bacfro_Ex : symmetric_bacfro_Ex)));
	stage_E13y=((wz==0 || !is_calc_ey) ? please_do_nuthing : ((wz==1) ? periodic_bacfro_Ey : ((wz==3) ? bloch_bacfro_Ey : symmetric_bacfro_Ey)));
	stage_E13z=((wz==0 || !is_calc_ez) ? please_do_nuthing : ((wz==1) ? periodic_bacfro_Ez : ((wz==3) ? bloch_bacfro_Ez : symmetric_bacfro_Ez)));
	
	
	stage_H11x=((wx==0 || !is_calc_hx) ? please_do_nuthing : ((wx==1) ? periodic_lefrig_Hx : ((wx==3) ? bloch_lefrig_Hx : symmetric_lefrig_Hx)));
	stage_H11y=((wx==0 || !is_calc_hy) ? please_do_nuthing : ((wx==1) ? periodic_lefrig_Hy : ((wx==3) ? bloch_lefrig_Hy : symmetric_lefrig_Hy)));
	stage_H11z=((wx==0 || !is_calc_hz) ? please_do_nuthing : ((wx==1) ? periodic_lefrig_Hz : ((wx==3) ? bloch_lefrig_Hz : symmetric_lefrig_Hz)));
	
	stage_H12x=((wy==0 || !is_calc_hx) ? please_do_nuthing : ((wy==1) ? periodic_bottop_Hx : ((wy==3) ? bloch_bottop_Hx : symmetric_bottop_Hx)));
	stage_H12y=((wy==0 || !is_calc_hy) ? please_do_nuthing : ((wy==1) ? periodic_bottop_Hy : ((wy==3) ? bloch_bottop_Hy : symmetric_bottop_Hy)));
	stage_H12z=((wy==0 || !is_calc_hz) ? please_do_nuthing : ((wy==1) ? periodic_bottop_Hz : ((wy==3) ? bloch_bottop_Hz : symmetric_bottop_Hz)));
	
	stage_H13x=((wz==0 || !is_calc_hx) ? please_do_nuthing : ((wz==1) ? periodic_bacfro_Hx : ((wz==3) ? bloch_bacfro_Hx : symmetric_bacfro_Hx)));
	stage_H13y=((wz==0 || !is_calc_hy) ? please_do_nuthing : ((wz==1) ? periodic_bacfro_Hy : ((wz==3) ? bloch_bacfro_Hy : symmetric_bacfro_Hy)));
	stage_H13z=((wz==0 || !is_calc_hz) ? please_do_nuthing : ((wz==1) ? periodic_bacfro_Hz : ((wz==3) ? bloch_bacfro_Hz : symmetric_bacfro_Hz)));

	// 1D - FDTD ??
	for(int s=0;s<n_SRCS;s++)
	{
	 stage_E14[s]=(source_type[s]==0) ? update_i_E : please_do_nothing;
	 stage_H14[s]=(source_type[s]==0) ? update_i_H : please_do_nothing;
    }

	/////////////////////////////////////////////////////////////////////

	if (!silentmode)
	{
		if (!noreport) std::cout<<report(args[3],args[2],args[1],ressurect,1)<<"\n";
		std::cout<<" *simulation starts now...\n";
	}

    ///////////////////////////////////////////////////////////////////////
    //                      glowna petla programu                        //
    //          obliczamy pole EM w  kolejnych chwilach czasu            //
    ///////////////////////////////////////////////////////////////////////

    for(t=ts;t<=maxt;t++)
    {
		// update pola el - centrum
        stage_E1x();stage_E1y();stage_E1z();
		// update UPML
        stage_E2x_1();stage_E3x_1();stage_E4x_1();
        stage_E2y_1();stage_E3y_1();stage_E4y_1();
        stage_E2z_1();stage_E3z_1();stage_E4z_1();
     
        stage_E2x_2();stage_E3x_2();stage_E4x_2();
        stage_E2y_2();stage_E3y_2();stage_E4y_2();
        stage_E2z_2();stage_E3z_2();stage_E4z_2();  
        
        // zrodla 1;
        for(int s=0;s<n_SRCS;s++)
        {
           stage_E5_1_1[s](s);stage_E6_1_1[s](s);stage_E7_1_1[s](s);
           stage_E5_1_2[s](s);stage_E6_1_2[s](s);stage_E7_1_2[s](s);
           stage_E5_2_1[s](s);stage_E6_2_1[s](s);stage_E7_2_1[s](s);
           stage_E5_2_2[s](s);stage_E6_2_2[s](s);stage_E7_2_2[s](s);		
        }

        // osrodki dyspersyjne
		stage_E8x();stage_E8y();stage_E8z();
		stage_E9x();stage_E9y();stage_E9z();
		stage_E10x();stage_E10y();stage_E10z();

         // osrodki dyspersyjne 2 : prady
		stage_J8x();stage_J8y();stage_J8z();
		stage_J9x();stage_J9y();stage_J9z();
		stage_J10x();stage_J10y();stage_J10z();

        // war periodyczne i symetryczne
		stage_E11x();stage_E11y();stage_E11z();
        stage_E12x();stage_E12y();stage_E12z();
        stage_E13x();stage_E13y();stage_E13z();

		// dodatkowe 1d FDTD
		for(int s=0;s<n_SRCS;s++) stage_E14[s](s);
		// wyrzucamy wyniki na wyjscie
		for(p=0;p<n_outputs;p++) if ((t>=out_min_t[p])&&(t<=out_max_t[p])&&(!((t-out_min_t[p])%out_t_step[p])))
		{
			if (!silentmode) std::cout<<" *saving results ...\n";
			if (!put_out1(args[3],p)) perror("error(s) during saving results...");
		}
         		
		// aktualizacja detektorow part 1
		for(p=0;p<detector_n;p++) if ((t>=detector_t_start[p])&&(t<=detector_t_end[p])&&(!((t-detector_t_end[p])%detector_t_step[p])))
	    {
         if (detector_skladowa[p]>2)
         {
            recorded_at_detector[p][(int)(t-detector_t_start[p])/detector_t_step[p]]=(typ_pola)0.5*cenfun[detector_skladowa[p]](detector_x[p],detector_y[p],detector_z[p]);
         }
        }
        
		// update pola mag - centrum
        stage_H1x();stage_H1y();stage_H1z();

		// update UPML
        stage_H2x_1();stage_H3x_1();stage_H4x_1();
        stage_H2y_1();stage_H3y_1();stage_H4y_1();
        stage_H2z_1();stage_H3z_1();stage_H4z_1();

        stage_H2x_2();stage_H3x_2();stage_H4x_2();
        stage_H2y_2();stage_H3y_2();stage_H4y_2();
        stage_H2z_2();stage_H3z_2();stage_H4z_2();

		// zrodla 2
        for(int s=0;s<n_SRCS;s++)
        {         
         stage_H5_1_1[s](s);stage_H6_1_1[s](s);stage_H7_1_1[s](s);
         stage_H5_1_2[s](s);stage_H6_1_2[s](s);stage_H7_1_2[s](s);
         stage_H5_2_1[s](s);stage_H6_2_1[s](s);stage_H7_2_1[s](s);
         stage_H5_2_2[s](s);stage_H6_2_2[s](s);stage_H7_2_2[s](s);
        }

        // osrodki dyspersyjne
		stage_H8x();stage_H8y();stage_H8z();
		stage_H9x();stage_H9y();stage_H9z();
		stage_H10x();stage_H10y();stage_H10z();

        // osrodki dyspersyjne 2 : prady
		stage_Jh8x();stage_Jh8y();stage_Jh8z();
		stage_Jh9x();stage_Jh9y();stage_Jh9z();
		stage_Jh10x();stage_Jh10y();stage_Jh10z();

		// war periodyczne i symetryczne
		stage_H11x();stage_H11y();stage_H11z();
        stage_H12x();stage_H12y();stage_H12z();
        stage_H13x();stage_H13y();stage_H13z(); 
		// dodatkowe 1D FDTD
		for(int s=0;s<n_SRCS;s++) stage_H14[s](s);
		// wyrzucamy wyniki na wyjscie
		for(p=0;p<n_outputs;p++) if ((t>=out_min_t[p])&&(t<=out_max_t[p])&&(!((t-out_min_t[p])%out_t_step[p])))
		{
			if (!silentmode) std::cout<<" *saving results ...\n";
			if (!put_out2(args[3],p)) perror("error(s) during saving results...");
		}
		// aktualizacja detektorow part 2
		for(p=0;p<detector_n;p++) if ((t>=detector_t_start[p])&&(t<=detector_t_end[p])&&(!((t-detector_t_end[p])%detector_t_step[p])))
	    {
         if (detector_skladowa[p]>2)
         {
            recorded_at_detector[p][(int)(t-detector_t_start[p])/detector_t_step[p]]+=(typ_pola)0.5*cenfun[detector_skladowa[p]](detector_x[p],detector_y[p],detector_z[p]);
         }
         else recorded_at_detector[p][(int)(t-detector_t_start[p])/detector_t_step[p]]=cenfun[detector_skladowa[p]](detector_x[p],detector_y[p],detector_z[p]);        
        }
		
		// zrzut;;
		if ((!(t%backup_delta_t)) && (t!=maxt) && (t!=0))
		{
			if (!silentmode) std::cout<<" *making backup...\n";
			if (!save_all(args)) perror("error(s) during backup..");
		}
		stage_0();
	}
	////////////////////////////////////////////////////////////////////////////
	//zapis tablic z detektorami
	for(p=0;p<detector_n;p++)
	{       if (t>=detector_t_end[p])
	        {
			 if (!silentmode) std::cout<<" *saving results from detector...\n";
			 if (!put_out_detector(args[3],p)) perror("error(s) during saving results...");
            }
    }
    ///////////////////////////////////////////////////////////////////////
    //kreowanie raportu
	if (!noreport)
	{
		if (!silentmode && !iosilentmode)  std::cout<<" *creating report...\n";
    	endtime=time(0);
    	create_report(args[3],args[2],args[1],ressurect);
	}

    //czyszczenie pamieci
	if (!silentmode) std::cout<<" *freeing memory...\n";
	deinit_structure_arrays();
	deinit_arrays();
	deinit_other_arrays();
    return 0;
}
