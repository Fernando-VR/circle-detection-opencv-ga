
#include "../Headers/AGS.hpp"


//const float ProbabilidadCruza=0.8;
//const float ProbabilidadMuta=0.01;
//const unsigned int NumMaxGeneraciones=100;
//const unsigned int NUMEROdeINDIVIDUOS=20;
//const unsigned int NUMEROdeGENES=2;
//unsigned char NUMEROdeBITporGEN[NUMEROdeGENES]={10,10};
//float LIMITES_SUP[NUMEROdeGENES]={ 10,  10};
//float LIMITES_INf[NUMEROdeGENES]={-10, -10};
//
//
//
//int main()
//{
//    mainGenetico;
//    return 0;
//}


int mainGenetico( )
{ POBLACION *Old;
  POBLACION *New;
  POBLACION *Aux;
  unsigned int* ListaSeleccionados, i, generacion=0;
  Old=CrearPoblacion(NUMEROdeINDIVIDUOS,NUMEROdeGENES,NUMEROdeBITporGEN,LIMITES_SUP,LIMITES_INf);
  New=CrearPoblacion(NUMEROdeINDIVIDUOS,NUMEROdeGENES,NUMEROdeBITporGEN,LIMITES_SUP,LIMITES_INf);
  srand(time(NULL));

  InicializarPoblacion(Old);
  DecodificaReal(Old);
  EvaluarPoblacion(Old);
  Obj2Fit(Old);
  printf("\npoblacion inicial: ");
  ImprimePoblacion(Old);

  while(generacion<NumMaxGeneraciones)
       { ListaSeleccionados=Selecion(Old);
         Cruza(ProbabilidadCruza,Old,New,ListaSeleccionados);
         Muta(ProbabilidadMuta,New);
         DecodificaReal(New);
         EvaluarPoblacion(New);
         Elitismo(Old,New);
         EvaluarPoblacion(New);
         Obj2Fit(New);
         generacion++;
         free(ListaSeleccionados);
         printf("\ngeneracion %i",generacion);
         ImprimePoblacion(New);
         Aux=New;
         New=Old;
         Old=Aux;
       }
  //ImprimePoblacion(New);
  EliminarPoblacion(Old);
  EliminarPoblacion(New);
  printf(" \n");
  return 0;
}

// --------------------------- FUNCIONES ---------------------------------------------------
/*
   Problema de ejemplo:
   Maximizar la funcion:
   f(x,y)= 50 -(x-5)^2 -(y-5)^2;

   Que combinación de Valores de x,y nos proporcina el mayor valor como resultaodo
   del cálculo de la función f(x,y)?
   -10 < x < 10
   -10 < y < 10

   máximo se encuentra para x=5, y=5
 */

/*
   Problema de ejemplo:
   Maximizar la funcion (Rastrigin):
   Z=100 - (( XX.^2 - 10*cos(2*pi*XX) + 10 ) + ( YY.^2 - 10*cos(2*pi*YY) + 10 ));

   Que combinación de Valores de x,y nos proporcina el mayor valor como resultaodo
   del cálculo de la función f(x,y)?
   -5.12 < x < 5.12
   -5.12 < y < 5.12

   máximo se encuentra para x=0, y=0

 */
void Elitismo(POBLACION *OldPop, POBLACION *NewPob)
{ unsigned int k;

  for(k=0; k<NewPob->Pob[NewPob->IdMin].TamCromosoma; k++)
      NewPob->Pob[NewPob->IdMin].Cromosoma[k]=OldPop->Pob[OldPop->IdBest].Cromosoma[k];

  for(k=0; k<NewPob->Pob[NewPob->IdMin].NumGenes; k++)
     {  NewPob->Pob[NewPob->IdMin].ValorEnt[k]=OldPop->Pob[OldPop->IdBest].ValorEnt[k];
        NewPob->Pob[NewPob->IdMin].ValorReal[k]=OldPop->Pob[OldPop->IdBest].ValorReal[k];
      }
  NewPob->Pob[NewPob->IdMin].ValorObj=OldPop->Pob[OldPop->IdBest].ValorObj;
  NewPob->Pob[NewPob->IdMin].Fit=OldPop->Pob[OldPop->IdBest].Fit;
}

void Muta(float Pmuta, POBLACION *NewPop)
{ unsigned int i,k;
  float  ValorAleatorio;
  for(i=0; i<NewPop->NumIndividuos; i++)
  for(k=0; k<NewPop->Pob[i].TamCromosoma; k++)
     { ValorAleatorio=(float)rand()/RAND_MAX;
       if(ValorAleatorio<=Pmuta)
          NewPop->Pob[i].Cromosoma[k]=!NewPop->Pob[i].Cromosoma[k];
     }
}
void Cruza(float Pcruza, POBLACION *OldPop, POBLACION *NewPob, unsigned int* ListaSeleccion)
{  unsigned int i,k, PuntoDeCruza;
   float  ValorAleatorio;

   ValorAleatorio=(float)rand()/RAND_MAX;
   if(ValorAleatorio<=Pcruza)
     { for(i=0; i<OldPop->NumIndividuos; i+=2)
          { PuntoDeCruza=( (float)rand()/RAND_MAX ) * (OldPop->Pob[i].TamCromosoma-2);
            //printf("\nPunto de Cruza: %i",PuntoDeCruza);
            for(k=0; k<=PuntoDeCruza; k++)
               { NewPob->Pob[i].Cromosoma[k]=OldPop->Pob[ ListaSeleccion[i+1] ].Cromosoma[k];
                 NewPob->Pob[i+1].Cromosoma[k]=OldPop->Pob[ ListaSeleccion[i] ] .Cromosoma[k];
               }
            for(k=PuntoDeCruza+1; k<OldPop->Pob[i].TamCromosoma; k++)
               { NewPob->Pob[i].Cromosoma[k]=OldPop->Pob[ ListaSeleccion[i] ].Cromosoma[k];
                 NewPob->Pob[i+1].Cromosoma[k]=OldPop->Pob[ ListaSeleccion[i+1] ].Cromosoma[k];
                }
           }
      }
  else
     { for(i=0; i<OldPop->NumIndividuos; i+=2)
       for(k=0; k<OldPop->Pob[i].TamCromosoma; k++)
         { NewPob->Pob[i].Cromosoma[k]=OldPop->Pob[ ListaSeleccion[i] ].Cromosoma[k];
           NewPob->Pob[i+1].Cromosoma[k]=OldPop->Pob[ ListaSeleccion[i+1] ] .Cromosoma[k];
           }
     }
}


 unsigned int* Selecion(POBLACION *ptr)
{ float SumaTotal=0;
  float Proporciones[ptr->NumIndividuos];
  float suma,ValorAleatorio;
  unsigned int  k,i,Id_Seleccionado;
  unsigned int* IndividuosSeleccionados;

  IndividuosSeleccionados=(unsigned int*) malloc( sizeof(unsigned int)*ptr->NumIndividuos );
  if(IndividuosSeleccionados==NULL)
    { printf("\n Error al reservar la memoria para los individuos seleccionados.");
      exit(0);
    }
  for(i=0; i<ptr->NumIndividuos; i++)
     SumaTotal+=ptr->Pob[i].Fit;
  for(i=0; i<ptr->NumIndividuos; i++)
      Proporciones[i]=ptr->Pob[i].Fit / SumaTotal;

  for(k=0; k<ptr->NumIndividuos; k++)
    { suma=0;
      ValorAleatorio=(float)rand()/RAND_MAX;
      for(i=0; i<ptr->NumIndividuos; i++)
         { suma+=Proporciones[i];
           if(suma>ValorAleatorio)
             { IndividuosSeleccionados[k]=i;
               break;
              }
          }
     }
 return IndividuosSeleccionados;
}

void Obj2Fit(POBLACION *ptr)
{ float Rango;
  unsigned int i;
  Rango= ptr->Pob[ptr->IdBest].ValorObj - ptr->Pob[ptr->IdMin].ValorObj;
  for(i=0; i<ptr->NumIndividuos; i++)
     ptr->Pob[i].Fit=( (ptr->Pob[i].ValorObj - ptr->Pob[ptr->IdMin].ValorObj ) / Rango )*100;
}

float FuncionObjetivo(float *Valor,unsigned int Ne)
{ float R;
  R=50-pow(Valor[0]-5.0,2)-pow(Valor[1]-5.0,2);
  return R;
}

/*
float FuncionObjetivo(float *Valor,unsigned int Ne)
{ unsigned int i;
  float aux,R;
  const float pi=3.14159;

  aux=0;
  for(i=0; i<Ne; i++)
      aux+= pow(Valor[i],2) - 10.0*cos(2.0*pi*Valor[i]) + 10;
  R=100.0 - aux;
  return -1*R;
}
*/
void EvaluarPoblacion(POBLACION *ptr)
{ unsigned int i, Id_max=0, Id_min=0;
  for(i=0; i<ptr->NumIndividuos; i++)
     { ptr->Pob[i].ValorObj=FuncionObjetivo(ptr->Pob[i].ValorReal,ptr->Pob[i].NumGenes);
       if(ptr->Pob[Id_max].ValorObj<ptr->Pob[i].ValorObj)
          Id_max=i;
       if(ptr->Pob[Id_min].ValorObj>ptr->Pob[i].ValorObj)
          Id_min=i;
       }
  ptr->IdBest=Id_max;
  ptr->IdMin=Id_min;
}

void DecodificaReal(POBLACION *ptr)
{  unsigned int i,j,k,Inicio;
   float Rango;

   for(i=0; i<ptr->NumIndividuos; i++)
      { Inicio=0;
        for(k=0; k<ptr->Pob[i].NumGenes; k++)
           {
             ptr->Pob[i].ValorEnt[k]=0;
             for(j=Inicio; j<(ptr->Pob[i].BitPorGen[k]+Inicio); j++)
                 ptr->Pob[i].ValorEnt[k]+=ptr->Pob[i].Cromosoma[j]*pow(2,j-Inicio);
             Inicio+=ptr->Pob[i].BitPorGen[k];
             Rango= ( ptr->Pob[i].LSup[k] - ptr->Pob[i].LInf[k] ) / ( pow(2,ptr->Pob[i].BitPorGen[k]) - 1 );
             ptr->Pob[i].ValorReal[k]= (ptr->Pob[i].ValorEnt[k] * Rango) + ptr->Pob[i].LInf[k];
           }
       }
}

void ImprimeValoresReales(POBLACION *ptr, unsigned int Id)
{ unsigned int j;
  printf("{ ");
  for(j=0; j<ptr->Pob[Id].NumGenes; j++)
      printf("%f ",ptr->Pob[Id].ValorReal[j]);
      printf("}");
}

void ImprimeValoresEntero(POBLACION *ptr, unsigned int Id)
{ unsigned int j;
  printf("{ ");
  for(j=0; j<ptr->Pob[Id].NumGenes; j++)
      printf("%i ",ptr->Pob[Id].ValorEnt[j]);
  printf("}");
}
void DecodificaEntero(POBLACION *ptr)
{  unsigned int i,j,k,Inicio;

   for(i=0; i<ptr->NumIndividuos; i++)
      { Inicio=0;
        for(k=0; k<ptr->Pob[i].NumGenes; k++)
           { ptr->Pob[i].ValorEnt[k]=0;
             for(j=Inicio; j<(ptr->Pob[i].BitPorGen[k]+Inicio); j++)
                 ptr->Pob[i].ValorEnt[k]+=ptr->Pob[i].Cromosoma[j]*pow(2,j-Inicio);
             Inicio+=ptr->Pob[i].BitPorGen[k];
           }
       }
}

void ImprimeIndividuo(POBLACION *ptr, unsigned int Id)
{ unsigned int j,NG,aux;
  NG=0;
  aux=ptr->Pob[Id].BitPorGen[NG];
  printf("\n%i: ",Id);
  for(j=0; j<ptr->Pob[Id].TamCromosoma; j++)
     { printf("%i",ptr->Pob[Id].Cromosoma[j]);
        if(j==(aux-1))
           { printf(":");
             NG++;
             aux+=ptr->Pob[Id].BitPorGen[NG]; }
       }
  printf(" ng:%i {",ptr->Pob[Id].NumGenes);

  for(j=0; j<ptr->Pob[Id].NumGenes; j++)
      printf("%i,",ptr->Pob[Id].BitPorGen[j]);
  printf(" }");
  printf(" Obj:%f Fit: %f",ptr->Pob[Id].ValorObj,ptr->Pob[Id].Fit);
}

void ImprimePoblacion(POBLACION *ptr)
{ int i;

   for(i=0; i<ptr->NumIndividuos; i++)
       ImprimeIndividuo(ptr,i);
  printf("\n Id_Best: %i",ptr->IdBest);
  printf("\n Id_Min: %i",ptr->IdMin);
}

void InicializarPoblacion(POBLACION *ptr)
{  int i,j;
   for(i=0; i<ptr->NumIndividuos; i++)
   for(j=0; j<ptr->Pob[i].TamCromosoma; j++)
   	   ptr->Pob[i].Cromosoma[j]=rand()%2;
}

POBLACION* CrearPoblacion(unsigned int NumInd, unsigned int NumGen, unsigned char* NumBitxGen, float* LimitesSup, float* LimitesInf)
{ POBLACION *ptr;
  int k, TotalBitsDelCromosoma=0;

  ptr=(POBLACION *)malloc(sizeof(POBLACION));
  if(ptr==NULL)
    { printf("\n Error al reservar la memoria para la estrutura poblacion.");
      exit(0);
    }
  ptr->NumIndividuos=NumInd;
  ptr->Pob=(INDIVIDUO*)malloc( NumInd*sizeof(INDIVIDUO) );
  if(ptr->Pob==NULL)
    { printf("\n Error al reservar la memoria para los individuos de la poblacion.");
      exit(0);
    }
  for(k=0; k<NumGen; k++)
     TotalBitsDelCromosoma+=NumBitxGen[k];

  for(k=0; k<NumInd; k++)
  	 { ptr->Pob[k].Cromosoma=(unsigned char*)malloc(TotalBitsDelCromosoma*sizeof(char));
  	   ptr->Pob[k].BitPorGen=NumBitxGen;
       ptr->Pob[k].NumGenes=NumGen;
       ptr->Pob[k].TamCromosoma=TotalBitsDelCromosoma;
       ptr->Pob[k].ValorEnt=(unsigned int*)malloc(NumGen*sizeof(unsigned int));
       ptr->Pob[k].ValorReal=(float *)malloc(NumGen*sizeof(float));
       ptr->Pob[k].LInf=LimitesInf;
       ptr->Pob[k].LSup=LimitesSup;
  	 }
  return ptr;
}

void EliminarPoblacion(POBLACION *ptr)
{ int k;
  if(ptr!=NULL)
    { for(k=0; k<ptr->NumIndividuos; k++)
  	      { free(ptr->Pob[k].Cromosoma);
            free(ptr->Pob[k].ValorEnt);
            free(ptr->Pob[k].ValorReal);
            }
      free(ptr->Pob);
      free(ptr);                     }
  else
  	 printf("\nNo se libero la memoria, poblacion no inicializada.");

}
