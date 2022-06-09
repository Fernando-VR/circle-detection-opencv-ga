#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>



extern const float ProbabilidadCruza;
extern const float ProbabilidadMuta;
extern const unsigned int NumMaxGeneraciones;
extern const unsigned int NUMEROdeINDIVIDUOS;
extern const unsigned int NUMEROdeGENES;
extern unsigned char NUMEROdeBITporGEN[];
extern float LIMITES_SUP[];
extern float LIMITES_INf[];


typedef struct{ unsigned char *Cromosoma;
	              unsigned char *BitPorGen;
	              unsigned int  NumGenes;
	              unsigned int  TamCromosoma;
                unsigned int  *ValorEnt;
                float         *ValorReal;
                float         *LInf;
                float         *LSup;
                float         ValorObj;
                float         Fit;
               }INDIVIDUO;

typedef struct{ INDIVIDUO *Pob;
	              unsigned int NumIndividuos;
                unsigned int IdBest;
                unsigned int IdMin;
               }POBLACION;

int mainGenetico();

POBLACION* CrearPoblacion(unsigned int NumInd, unsigned int NumGen, unsigned char* NumBitxGen, float* LimitesSup, float* LimitesInf);
void EliminarPoblacion(POBLACION *ptr);
void InicializarPoblacion(POBLACION *ptr);
void ImprimePoblacion(POBLACION *ptr);
void ImprimeIndividuo(POBLACION *ptr, unsigned int Id);
void DecodificaEntero(POBLACION *ptr);
void DecodificaReal(POBLACION *ptr);
void ImprimeValoresEntero(POBLACION *ptr, unsigned int Id);
void ImprimeValoresReales(POBLACION *ptr, unsigned int Id);
void EvaluarPoblacion(POBLACION *ptr);
void Obj2Fit(POBLACION *ptr);
float FuncionObjetivo(float *Valor,unsigned int Ne);
unsigned int* Selecion(POBLACION *ptr);
void Cruza(float Pcruza, POBLACION *OldPop, POBLACION *NewPob, unsigned int* ListaSeleccion);
void Muta(float Pmuta, POBLACION *NewPop);
void Elitismo(POBLACION *OldPop, POBLACION *NewPob);
