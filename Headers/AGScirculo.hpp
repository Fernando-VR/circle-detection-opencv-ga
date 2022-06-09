

#include "AGS.hpp"

//int Genetico_Circulo_Parametros(void)
int Genetico_Circulo_Parametros(float r[], float c[], int nEdges, float &X0, float &Y0, float &R0) ;
void EvaluarPoblacionCirculo(POBLACION *ptr, float r[], float c[], int nEdges);
float FuncionObjetivoCirculo(float *Valor, unsigned int Ne, float r[], float c[], int nEdges);
int syParametrosDelCirculoDadosTresPuntos (float xi, float yi, float xj, float yj, float xk, float yk, float &X0, float &Y0, float &R0);
float DistanciaRadial(float ren, float col, float X0, float Y0);
