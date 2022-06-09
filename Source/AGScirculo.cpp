
#include "../Headers/AGScirculo.hpp"


const float ProbabilidadCruza=0.8;
const float ProbabilidadMuta=0.01;
const unsigned int NumMaxGeneraciones=35;
const unsigned int NUMEROdeINDIVIDUOS=120;
const unsigned int NUMEROdeGENES=3;
unsigned char NUMEROdeBITporGEN[NUMEROdeGENES]={10,10,10};
float LIMITES_SUP[NUMEROdeGENES]={ 10000,  10000, 10000};
float LIMITES_INf[NUMEROdeGENES]={0, 0, 0};



//int main()
//{
//    Genetico_Circulo_Parametros;
//    return 0;
//}


//int mainGenetico( )
//int Genetico_Circulo_Parametros(void)
int Genetico_Circulo_Parametros(float r[], float c[], int nEdges, float &X0, float &Y0, float &R0)
{
    //cambio
  LIMITES_SUP[0] = (float)(nEdges-1);
  LIMITES_SUP[1] = (float)(nEdges-1);
  LIMITES_SUP[2] = (float)(nEdges-1);
  //Fin del cambio


  POBLACION *Old;
  POBLACION *New;
  POBLACION *Aux;
  unsigned int* ListaSeleccionados, /*i,*/ generacion=0;
  Old=CrearPoblacion(NUMEROdeINDIVIDUOS,NUMEROdeGENES,NUMEROdeBITporGEN,LIMITES_SUP,LIMITES_INf);
  New=CrearPoblacion(NUMEROdeINDIVIDUOS,NUMEROdeGENES,NUMEROdeBITporGEN,LIMITES_SUP,LIMITES_INf);
  srand(time(NULL));

  InicializarPoblacion(Old);
  DecodificaReal(Old);

  EvaluarPoblacionCirculo(Old, r, c, nEdges);             // CAMBIO

  Obj2Fit(Old);
//  printf("\npoblacion inicial: ");
//  ImprimePoblacion(Old);

  while(generacion<NumMaxGeneraciones)
       { ListaSeleccionados=Selecion(Old);
         Cruza(ProbabilidadCruza,Old,New,ListaSeleccionados);
         Muta(ProbabilidadMuta,New);
         DecodificaReal(New);
         EvaluarPoblacionCirculo(New, r, c, nEdges);     //CAMBIO
         Elitismo(Old,New);
         EvaluarPoblacionCirculo(New, r, c, nEdges);     //CAMBIO
         Obj2Fit(New);
         generacion++;
         free(ListaSeleccionados);
         printf("\ngeneracion %i",generacion);
//         ImprimePoblacion(New);
         Aux=New;
         New=Old;
         Old=Aux;
       }
  //ImprimePoblacion(New);

  //Decodificacion del resultado
  int i = (int) (New->Pob[New->IdBest].ValorReal[0]);
  int j = (int) (New->Pob[New->IdBest].ValorReal[1]);
  int k = (int) (New->Pob[New->IdBest].ValorReal[2]);
  float xi = r[i];  float yi = c[i];
  float xj = r[j];  float yj = c[j];
  float xk = r[k];  float yk = c[k];

  syParametrosDelCirculoDadosTresPuntos(xi, yi, xj, yj, xk, yk, X0, Y0, R0);

  EliminarPoblacion(Old);
  EliminarPoblacion(New);
//  printf(" \n");
  return 0;
}

void EvaluarPoblacionCirculo(POBLACION *ptr, float r[], float c[], int nEdges)
{ unsigned int i, Id_max=0, Id_min=0;
  for(i=0; i<ptr->NumIndividuos; i++)
     { ptr->Pob[i].ValorObj=FuncionObjetivoCirculo(ptr->Pob[i].ValorReal,ptr->Pob[i].NumGenes, r, c, nEdges);    //CAMBIO
       if(ptr->Pob[Id_max].ValorObj<ptr->Pob[i].ValorObj)
          Id_max=i;
       if(ptr->Pob[Id_min].ValorObj>ptr->Pob[i].ValorObj)
          Id_min=i;
       }
  ptr->IdBest=Id_max;
  ptr->IdMin=Id_min;
}



float FuncionObjetivoCirculo(float *Valor, unsigned int Ne, float r[], float c[], int nEdges)
{
    //Se necesitan 3 puntos con indices i, j, k para con ellos determinar x0, y0 r0
    int i = (int) Valor[0];
    int j = (int) Valor[1];
    int k = (int) Valor[2];

    float xi = r[i];    float yi = c[i];
    float xj = r[j];    float yj = c[j];
    float xk = r[k];    float yk = c[k];

    float X0=0, Y0=0, R0=0;
    syParametrosDelCirculoDadosTresPuntos(xi, yi, xj, yj, xk, yk, X0, Y0, R0);

    if(R0<20)       //Solo quiero circulos de al menos 20 pixeles de radio
        return 0.0;

    float contador = 0.0;
    for(int k=0; k<nEdges; k++)
        if( fabs( DistanciaRadial(r[k], c[k], X0, Y0) - R0 ) <0.5 )
            contador++;

//    return 100.0*contador / (2.0*3.14159216*R0);
    return contador;
}

float DistanciaRadial(float ren, float col, float X0, float Y0)
{
    return sqrt( (ren-X0)*(ren-X0) + (col-Y0)*(col-Y0) );
}

int syParametrosDelCirculoDadosTresPuntos (float xi, float yi, float xj, float yj, float xk, float yk, float &X0, float &Y0, float &R0)
{
    X0=Y0=R0=0.0;

    float xi2 = xi*xi; float yi2 = yi*yi;
    float xj2 = xj*xj; float yj2 = yj*yj;
    float xk2 = xk*xk; float yk2 = yk*yk;

    float denom = 4.0*( (xj-xi)*(yk-yi)-(xk-xi)*(yj-yi) );
    if ( fabs(denom)>0.01 )
    {
        X0 = ( ( ( xj2+yj2-(xi2+yi2))*2.0*(yk-yi) ) - ( (xk2+yk2)-(xi2+yi2)*2.0*(yj-yi) ) )/denom;
        Y0 = ( ( ( xj2+yj2-(xi2+yi2))*2.0*(xj-xi) ) - ( (xj2+yj2)-(xi2+yi2)*2.0*(xk-xi) ) )/denom;
        R0 = sqrt( (xi-X0)*(xi-X0) + (yi-Y0)*(yi-Y0));
    }
    return 0;
}
