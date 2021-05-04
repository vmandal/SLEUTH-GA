#define MAIN_MODULE
#define MAIN
#ifdef MPI
#include <mpi.h>
#endif
#define CATCH_SIGNALS
#ifdef CATCH_SIGNALS
#include <signal.h>
#include <errno.h>
#endif
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "coeff_obj.h"
#include "igrid_obj.h"
#include "landclass_obj.h"
#include "globals.h"
#include "output.h"
#include "utilities.h"
#include "random.h"
#include "driver.h"
#include "input.h"
#include "scenario_obj.h"
#include "proc_obj.h"
#include "timer_obj.h"
#include "landclass_obj.h"
#include "pgrid_obj.h"
#include "color_obj.h"
#include "memory_obj.h"
#include "color_obj.h"
#include "stats_obj.h"
#include "transition_obj.h"
#include "ugm_macros.h"

int Generations;
int PopulationSize;
float MutationRate;
int 	MaxEvals;
int NumOffSpring;
int NumReplaced;

FILE * pFile;
FILE * pStatsFile;

/*****************************************************************************\
*******************************************************************************
**                                                                           **
**                        STATIC FUNCTION PROTOTYPES                         **
**                                                                           **
*******************************************************************************
\*****************************************************************************/
static void print_usage (char *binary);
#ifdef CATCH_SIGNALS
void catch (int signo);
#endif

/*****************************************************************************\
*******************************************************************************
**                                                                           **
**                      STATIC MEMORY FOR THIS OBJECT                        **
**                                                                           **
*******************************************************************************
\*****************************************************************************/
static int tracer;

/*****************************************************************************\
*******************************************************************************
**                                                                           **
**                               SCCS ID                                     **
**                                                                           **
*******************************************************************************
\*****************************************************************************/
char main_c_sccs_id[] = "@(#)main.c	1.629	12/4/00";

/******************************************************************************/
/* 
 /*
 /*
 /******************************************************************************/

// Default headers
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

double RandMax1()
{
    return ((double)rand())/((double)(RAND_MAX));
}

unsigned int RandMax(unsigned int Maximum )
{
unsigned int value = -1;
	while ((value < 0) || (value > Maximum)) {value = rand() % Maximum;}
	return value;

}

double randNormal01()
{
	double U1,U2,V1,V2;
	double S=2;
	while ( S>=1)
	{
		U1 = RandMax1();
		U2 = RandMax1();
		V1 = 2.0*U1 - 1.0;
		V2 = 2.0*U2 - 1.0;
		S = pow(V1,2)+pow(V2,2);
	};
	
	return V1*sqrt((-2.0*log(S))/S);
}

unsigned int GenerateMutationStep(unsigned int size )
{
	double random = randNormal01();
	double value = random * size;
	return (unsigned int)value;
}

/******************************************************************************/
/* Structure: 
 /* Members:  
 
 /******************************************************************************/
typedef struct Gene_Type{
    int params[5];
    float fitness;
}Gene_Type;

Gene_Type BEST_GENE;

/******************************************************************************/
/* Structure: 
 /* Members:  
 
 /******************************************************************************/
typedef struct Population_Type{
    Gene_Type     *Genes;
    Gene_Type     *TempGenes;
    unsigned int  PopulationSize;
    unsigned int  NumSubPopulations;
    float         MutationRate;
    float         MaxParamValue;
    Gene_Type     BestGene;
    Gene_Type     EliteGene;
    unsigned int  MaxEvals;
    unsigned int  CurrentEvals;
}Population_Type;

/******************************************************************************/
/* Function: initGene()
 /* Input:  NONE
 /* Output: Gene_Type
 /* Description:
 /*      Initializes new Gene to random weights and fitness to 0
 /******************************************************************************/
Gene_Type initGene()
{
    Gene_Type NewGene;
	
    // set all parameters for a gene to random values
    NewGene.params[0] = (int) RandMax(100);
    NewGene.params[1] = (int) RandMax(100);
    NewGene.params[2] = (int) RandMax(100);
    NewGene.params[3] = (int) RandMax(100);
    NewGene.params[4] = (int) RandMax(100);
    
    // set unevaluated fitness value to 0
    NewGene.fitness = 0.0;

    // SSA printf statements to see what are the assigned values for New Gene
    printf("New Gene Parameters are: [");
    for(int i=0;i<5;++i){
    	printf(" %.2d ",NewGene.params[i]);
    }
    printf(" ] \n");
    // return a copy of the new randomized gene
    return NewGene;
}

/******************************************************************************/
/* Function: 
 /* Input:  
 /* Output: 
 /* Description:
 
 /******************************************************************************/
void EvaluateGene( Gene_Type* Gene ,Population_Type* Pop){
    int i;
    
    assert(Gene);
    if(Gene)
    {
		InitRandom (scen_GetRandomSeed ());
		Pop->CurrentEvals++;
             
        	Gene->fitness = 0.0;
 	        coeff_SetCurrentDiffusion ((int)Gene->params[0]); // SSA double to int modification
                coeff_SetCurrentSpread ((int)Gene->params[1]);
                coeff_SetCurrentBreed ((int)Gene->params[2]);
                coeff_SetCurrentSlopeResist ((int)Gene->params[3]);
                coeff_SetCurrentRoadGravity ((int)Gene->params[4]);

	       drv_driver ();

              proc_IncrementNumRunsExecThisCPU ();
              if (scen_GetLogFlag ())
              {
                if (scen_GetLogTimingsFlag () > 1)
                {
                  scen_Append2Log ();
                  timer_LogIt (scen_GetLogFP ());
                  scen_CloseLog ();
                }
              }



              proc_IncrementCurrentRun ();
              if (proc_GetProcessingType () == TESTING)
              {
                stats_ConcatenateControlFiles ();
                if (scen_GetWriteCoeffFileFlag ())
                {
                  coeff_ConcatenateFiles ();
                }
                if (scen_GetWriteAvgFileFlag ())
                {
                  stats_ConcatenateAvgFiles ();
                }
                if (scen_GetWriteStdDevFileFlag ())
                {
                  stats_ConcatenateStdDevFiles ();
                }

                timer_Stop (TOTAL_TIME);
                if (scen_GetLogFlag ())
                {
                  scen_Append2Log ();
                  if (scen_GetLogTimingsFlag () > 0)
                  {
                    timer_LogIt (scen_GetLogFP ());
                  }
                  mem_LogMinFreeWGrids (scen_GetLogFP ());
                  scen_CloseLog ();
                }
                EXIT (0);
              }   
		 Gene->fitness = getProduct();
		 
		 if( Gene->fitness > BEST_GENE.fitness )
			BEST_GENE = *Gene;
        
    }
}

/******************************************************************************/
/* Function: 
 /* Input:  
 /* Output: 
 /* Description:
 
 /******************************************************************************/
void EvaluatePopulation( Population_Type* Pop ){
	int i;
	
	if( Pop )
	{
        for( i = 0 ; i < Pop->PopulationSize ; i++ )
        {
            EvaluateGene(&(Pop->Genes[i]), Pop);
        }
    }
}


void SortGenes(Gene_Type* Genes, unsigned int PopSize )
{
	Gene_Type TempGene;
	int i,j;
	int done;
	
	for( i = 1 ; i < PopSize ; i ++ )
	{
		TempGene = Genes[i];
		j = i - 1;
		done = 0;
		do{
			if( Genes[j].fitness < TempGene.fitness )
			{
				Genes[j+1] = Genes[j];
				j = j-1 ;
				if( j < 0 )
					done = 1;
				
			}
			else
			{
				done = 1;
			}
		}while( done != 1);
		Genes[j+1] = TempGene;
	}
}


/******************************************************************************/
/* Function: 
 /* Input:  
 /* Output: 
 /* Description:
 
 /******************************************************************************/
Population_Type* initPopulation(
								unsigned int PopSize,
								float         MutRate)
{
    Population_Type*    NewPopulation;
    int             	i;
	
    // allocate memory for population
    NewPopulation = (Population_Type*)malloc(sizeof(Population_Type));
	
    // check for good memory allocation
    if( NewPopulation != NULL )    
    {
        // set population variables by input 
        // Input restrictions:
        //      Population Size must be a multiple of 2
        //      SubPopulations are not currently supported TODO 
        //      Mutation Rate must be between 0 and 1 inclusive
        NewPopulation->PopulationSize    = (PopSize/2)*2 + ((PopSize%2)*2); // round up pop size to multiple of 2
        NewPopulation->MutationRate      = MutRate;
        NewPopulation->MaxParamValue     = 100;
//        NewPopulation->NumSubPopulations = (PopSize/2); //SSA Number of Sub Population Must be initially defined!
		
        NewPopulation->Genes     = (Gene_Type*)malloc( NewPopulation->PopulationSize * sizeof(Gene_Type));      
        NewPopulation->TempGenes = (Gene_Type*)malloc( NewPopulation->PopulationSize * 5  * sizeof(Gene_Type));
		//printf(" Allocations at %llx and %llx\n",(unsigned long long int)NewPopulation->Genes,(unsigned long long int)NewPopulation->TempGenes);
		
		srand(time(NULL));
		
        if( NewPopulation->Genes == NULL || NewPopulation->TempGenes == NULL )    
        {
            assert(0);
            if( NewPopulation->Genes == NULL )
                free(NewPopulation->Genes);
            if( NewPopulation->TempGenes == NULL )
                free(NewPopulation->TempGenes);
            free(NewPopulation);
        }    
        else
        {
		
			for( i = 0 ; i < NewPopulation->PopulationSize/2 ; i++ )
			{
				NewPopulation->Genes[i].params[0] = i*(100/(NewPopulation->PopulationSize/2));
				NewPopulation->Genes[i].params[1] = i*(100/(NewPopulation->PopulationSize/2));
				NewPopulation->Genes[i].params[2] = i*(100/(NewPopulation->PopulationSize/2));
				NewPopulation->Genes[i].params[3] = i*(100/(NewPopulation->PopulationSize/2));
				NewPopulation->Genes[i].params[4] = i*(100/(NewPopulation->PopulationSize/2));
			}
			
			for( ; i < NewPopulation->PopulationSize; i++ )
				NewPopulation->Genes[i] =  initGene();
			

        }        
		
    }
	
    return NewPopulation;    
}

/******************************************************************************/
/* Function: 
 /* Input:  
 /* Output: 
 /* Description:
 
 /******************************************************************************/
void DestroyPopulation( Population_Type* Pop)
{
    printf("Debug");
    if( Pop )
    {
        if( Pop->Genes )
            free( Pop->Genes );
        if( Pop->TempGenes )
            free( Pop->TempGenes );    
        
        free(Pop);
        Pop = NULL;
    }
    else
    {
        assert(0);
    }
}

/******************************************************************************/
/* Function: 
 /* Input:  
 /* Output: 
 /* Description:
 
 /******************************************************************************/
void DebugPrintGene( Gene_Type Gene)
{
    fprintf(pFile,"Params: %d %d %d %d %d Fitness: %f\n", 
            (int) Gene.params[0],
            (int) Gene.params[1],
            (int) Gene.params[2],
            (int) Gene.params[3],
            (int) Gene.params[4],
            Gene.fitness);
}          

void DebugPrintGeneStat( Gene_Type Gene)
{
    fprintf(pStatsFile," %d, %d, %d, %d, %d, %f\n", 
            (int) Gene.params[0],
            (int) Gene.params[1],
            (int) Gene.params[2],
            (int) Gene.params[3],
            (int) Gene.params[4],
            Gene.fitness);
}

void DebugPrintGeneOut( Gene_Type Gene)
{
    printf(" %d, %d, %d, %d, %d, %f\n", 
            (int) Gene.params[0],
            (int) Gene.params[1],
            (int) Gene.params[2],
            (int) Gene.params[3],
            (int) Gene.params[4],
            Gene.fitness);
}


/******************************************************************************/
/* Function: 
/* Input:  
/* Output: 
/* Description:

/******************************************************************************/
void DebugPrintPopulation( Population_Type* Pop)      
{
    int i = 0,j = 0;
    
    for( i = 0 ; i < Pop->PopulationSize ; i++ )
    {
    	printf("Problem Starts Here: \n");
    	printf("Number of Subpopulations: %d\n",Pop->NumSubPopulations);
    	printf("Number of Population Size: %d\n",Pop->PopulationSize);
   		printf("Divison: %d\n",Pop->PopulationSize/Pop->NumSubPopulations);

        if( Pop->NumSubPopulations >= 2 && i % (Pop->PopulationSize/Pop->NumSubPopulations) == 0 )
            fprintf(pFile,"Sub Population:%5d\n",j++);
        fprintf(pFile," Gene: %5d ", i );
        DebugPrintGene( Pop->Genes[i] );
    }
    
}

/******************************************************************************/
/* Function: 
/* Input:  
/* Output: 
/* Description:

/******************************************************************************/
void DebugPrintPopulationTemp( Population_Type* Pop)      
{
    int i = 0,j = 0;
    
    for( i = 0 ; i < Pop->PopulationSize ; i++ )
    {   
        if( Pop->NumSubPopulations >= 2 && i % (Pop->PopulationSize/Pop->NumSubPopulations) == 0 )
            fprintf(pFile,"Sub Population:%5d\n",j++);
        fprintf(pFile," Gene: %5d ", i );
        DebugPrintGene( Pop->TempGenes[i] );
    }
    
}

void DebugPrintPopulationOut( Population_Type* Pop)      
{
    int i = 0,j = 0;
    
    for( i = 0 ; i < Pop->PopulationSize ; i++ )
    {   
        if( Pop->NumSubPopulations >= 2 && i % (Pop->PopulationSize/Pop->NumSubPopulations) == 0 )
            printf("Sub Population:%5d\n",j++);
        printf(" Gene: %5d ", i );
        DebugPrintGeneOut( Pop->Genes[i] );
    }
    
}




/******************************************************************************/
/* Function: 
 /* Input:  
 /* Output: 
 /* Description:
 
 /******************************************************************************/
void MutateGene(Gene_Type* Gene, float MutationRate,Population_Type* Pop )
{
    int CurrentParam = 0;
	int MutationOccured = 0;
	int MutationSize;
    
    if(Gene)
    {
		//  printf("%x , %d ",(unsigned int)Gene,sizeof(Gene_Type));
		//  DebugPrintGene(*Gene);
		
		for( CurrentParam = 0 ; CurrentParam < 5 ; CurrentParam++ )
		{
				
				if( MutationRate > RandMax1() )
				{
					Gene->params[CurrentParam] = RandMax(100);
					MutationOccured = 1;
				}
				/*if( MutationRate > RandMax1() )
				{	
					MutationSize = GenerateMutationStep(25);
					if( rand() % 2 )
						Gene->params[CurrentParam] -= MutationSize;
					else
						Gene->params[CurrentParam] += MutationSize;
					
					if( Gene->params[CurrentParam] > 100 )
						Gene->params[CurrentParam] = 100;
					else if( Gene->params[CurrentParam] < 0 )
						Gene->params[CurrentParam] = 1;
						
					MutationOccured = 1;
				}*/
			
			
			
		}
		
		if( MutationOccured == 1)
		{
			EvaluateGene(Gene,Pop);
		}
		//  DebugPrintGene(*Gene);
    }    
}

void MutatePopulation( Population_Type* Pop )
{
    int i=0;
    
    if(Pop)
    {
		
		//        DebugPrintPopulation(Pop);
        for( i = Pop->PopulationSize/5; i < Pop->PopulationSize ; i ++ )
        {
			
            MutateGene(&(Pop->Genes[i]),Pop->MutationRate,Pop);
        }
    }
}



/******************************************************************************/
/* Function: 
 /* Input:  
 /* Output: 
 /* Description:
 
 /******************************************************************************/
void CrossOverGene( Gene_Type* GeneIn1,
				   Gene_Type* GeneIn2,
				   Gene_Type* GeneOut1,
				   Gene_Type* GeneOut2)
{
    unsigned int randomValue;
    unsigned int CurrentParam = 0;
	
	//    printf("Crossover Mask %x \n ",CrossOverMask );
    if(1)
	{	
		randomValue = rand() %5;
		while( CurrentParam < randomValue )
		{
			GeneOut1->params[CurrentParam] = GeneIn1->params[CurrentParam];
			GeneOut2->params[CurrentParam] = GeneIn2->params[CurrentParam];
			CurrentParam++;
		}
	
		if( 0 ) // mean crossover
		{
			while( CurrentParam < 5 )
			{
				GeneOut1->params[CurrentParam] = (GeneIn2->params[CurrentParam] + GeneIn1->params[CurrentParam])/2.0;
				GeneOut2->params[CurrentParam] = (GeneIn2->params[CurrentParam] + GeneIn1->params[CurrentParam])/2.0;
				CurrentParam++;
			}	
		}
		else // simple airithmetic crossover
		{
			while( CurrentParam < 5 )
			{
				GeneOut1->params[CurrentParam] = GeneIn2->params[CurrentParam];
				GeneOut2->params[CurrentParam] = GeneIn1->params[CurrentParam];
				CurrentParam++;
			}
		}
	}
	else
	{
		randomValue = rand();
			
		if( randomValue & 0x1 )
		{
			GeneOut1->params[0] = GeneIn2->params[0];
			GeneOut2->params[0] = GeneIn1->params[0];
		}
		else if(1)
		{
			GeneOut1->params[0] = (GeneIn2->params[0] + GeneIn1->params[0])/2.0;
		    GeneOut2->params[0] = (GeneIn2->params[0] + GeneIn1->params[0])/2.0;
		}
		else
		{
			GeneOut1->params[0] = GeneIn1->params[0];
			GeneOut2->params[0] = GeneIn2->params[0];
		}
		
		if( randomValue & 0x2 )
		{
			GeneOut1->params[1] = GeneIn2->params[1];
			GeneOut2->params[1] = GeneIn1->params[1];
		}
		else if(1)
		{
			GeneOut1->params[1] = (GeneIn2->params[1] + GeneIn1->params[1])/2.0;
		    GeneOut2->params[1] = (GeneIn2->params[1] + GeneIn1->params[1])/2.0;
		}
		else
		{
			GeneOut1->params[1] = GeneIn1->params[1];
			GeneOut2->params[1] = GeneIn2->params[1];
		}
		
		if( randomValue & 0x4 )
		{
			GeneOut1->params[2] = GeneIn2->params[2];
			GeneOut2->params[2] = GeneIn1->params[2];
		}
		else if(1)
		{
			GeneOut1->params[2] = (GeneIn2->params[2] + GeneIn1->params[2])/2.0;
		    GeneOut2->params[2] = (GeneIn2->params[2] + GeneIn1->params[2])/2.0;
		}
		else
		{
			GeneOut1->params[2] = GeneIn1->params[2];
			GeneOut2->params[2] = GeneIn2->params[2];
		}
		
		if( randomValue & 0x8 )
		{
			GeneOut1->params[3] = GeneIn2->params[3];
			GeneOut2->params[3] = GeneIn1->params[3];
		}
		else if(1)
		{
			GeneOut1->params[3] = (GeneIn2->params[3] + GeneIn1->params[3])/2.0;
		    GeneOut2->params[3] = (GeneIn2->params[3] + GeneIn1->params[3])/2.0;
		}
		else
		{
			GeneOut1->params[3] = GeneIn1->params[3];
			GeneOut2->params[3] = GeneIn2->params[3];
		}
		
		if( randomValue & 0x10 )
		{
			GeneOut1->params[4] = GeneIn2->params[4];
			GeneOut2->params[4] = GeneIn1->params[4];
		}
		else if(1)
		{
			GeneOut1->params[4] = (GeneIn2->params[4] + GeneIn1->params[4])/2.0;
		    GeneOut2->params[4] = (GeneIn2->params[4] + GeneIn1->params[4])/2.0;
		}
		else
		{
			GeneOut1->params[4] = GeneIn1->params[4];
			GeneOut2->params[4] = GeneIn2->params[4];
		}
	}
    //DebugPrintGene(*GeneIn1);
    //DebugPrintGene(*GeneIn2);
    //DebugPrintGene(*GeneOut1);
    //DebugPrintGene(*GeneOut2);
}    

void GenerateOffspring( Population_Type* Pop , unsigned int NumOffSpring ){
    float TotalFitness = 0;
	int	  Parent1, Parent2;
	int   ParentCand[4],tempValue; 
    float CurrentValue;
    float RandomValue;
    //SSA
    float TotalFitness_Test = 0.0;
    double RandomValue_Test;
    int i =0;
    int j;
	
    //currently only supporting fitness proportional selection methods
    for( i = 0 ; i < Pop->PopulationSize ; i++ )
    {
        TotalFitness += Pop->Genes[i].fitness;        
    }
	printf("Total fitness=%f\n", TotalFitness);
      
    for( i = 0 ; i < NumOffSpring/2 ; i++ )
    {
		if(0) //Fitness Proportional Selection
        {
			printf("Running in Fitness Proportional Selection now:\n");
			RandomValue_Test = ((double)rand()/(double)RAND_MAX)*TotalFitness_Test;
			printf("Random Value Double %f\n",RandomValue_Test);

			RandomValue = ((float)rand()/((float)RAND_MAX))*TotalFitness;
			        printf("Random Value Float %f\n",RandomValue);
			CurrentValue = 0.0;
			j = 0;
			while( CurrentValue < RandomValue && j <= Pop->PopulationSize )
			{
				j++;
				CurrentValue += Pop->Genes[j-1].fitness;
			}
			Parent1 = j-1;
	
			RandomValue = ((float)rand()/((float)RAND_MAX))*TotalFitness;
			        printf("Random Value %f\n",RandomValue);
			CurrentValue = 0.0;
			j = 0;
			while( CurrentValue < RandomValue && j <= Pop->PopulationSize )
			{
				j++;
				CurrentValue += Pop->Genes[j-1].fitness;
			}
			Parent2 = j-1;
			if( Parent1 == Parent2 )
			{
				Parent2 = (Parent1 < Pop->PopulationSize - 1 ? Parent1+1 : Parent1-1);
			}
		
			CrossOverGene( &Pop->Genes[Parent1],
						&Pop->Genes[Parent2],
						&Pop->TempGenes[2*i],
						&Pop->TempGenes[(2*i)+1]);
						
			EvaluateGene(&Pop->TempGenes[2*i],Pop);
			EvaluateGene(&Pop->TempGenes[(2*i)+1],Pop);\
		}
		else //Tournament Selection
		{
			printf("Running in Tournament Selection now:\n");
			ParentCand[0] = rand() % Pop->PopulationSize;
			ParentCand[1] = rand() % Pop->PopulationSize;
			ParentCand[2] = rand() % Pop->PopulationSize;
			ParentCand[3] = rand() % Pop->PopulationSize;

			// SSA while loops ensure that all the parent candidate numbers are different!
			while( ParentCand[0] == ParentCand[1] )
				ParentCand[1] = rand() % Pop->PopulationSize;
				
			while( ParentCand[0] == ParentCand[2] || 
				   ParentCand[1] == ParentCand[2] )
				ParentCand[2] = rand() % Pop->PopulationSize;
				
			while( ParentCand[0] == ParentCand[3] || 
				   ParentCand[1] == ParentCand[3] || 
				   ParentCand[2] == ParentCand[3] )
				   ParentCand[3] = rand() % Pop->PopulationSize;
			printf("ParentCandidates: [ %d %d %d %d ] \n",ParentCand[0],ParentCand[1],ParentCand[2],ParentCand[3]);

			if(Pop->Genes[ParentCand[0]].fitness < Pop->Genes[ParentCand[1]].fitness )
			{
				tempValue = ParentCand[0];
				ParentCand[0] = ParentCand[1];
				ParentCand[1] = tempValue;
			}
			
			if(Pop->Genes[ParentCand[1]].fitness < Pop->Genes[ParentCand[2]].fitness )
			{
				tempValue = ParentCand[1];
				ParentCand[1] = ParentCand[2];
				ParentCand[2] = tempValue;
				
				if(Pop->Genes[ParentCand[0]].fitness < Pop->Genes[ParentCand[1]].fitness )
				{
					tempValue = ParentCand[0];
					ParentCand[0] = ParentCand[1];
					ParentCand[1] = tempValue;
				}
			}
			
			if(Pop->Genes[ParentCand[2]].fitness < Pop->Genes[ParentCand[3]].fitness )
			{
				tempValue = ParentCand[2];
				ParentCand[2] = ParentCand[3];
				ParentCand[3] = tempValue;
					
				if(Pop->Genes[ParentCand[1]].fitness < Pop->Genes[ParentCand[2]].fitness )
				{
					tempValue = ParentCand[1];
					ParentCand[1] = ParentCand[2];
					ParentCand[2] = tempValue;
					
					if(Pop->Genes[ParentCand[0]].fitness < Pop->Genes[ParentCand[1]].fitness )
					{
						tempValue = ParentCand[0];
						ParentCand[0] = ParentCand[1];
						ParentCand[1] = tempValue;
					}
				}
			}
			
			
				
			RandomValue = RandMax1();
			
			if( RandomValue > .5 )
				Parent1= ParentCand[0];
			else if( RandomValue > .25 )
				Parent1= ParentCand[1];
			else if( RandomValue > .125 )
				Parent1= ParentCand[2];
			else
				Parent1= ParentCand[2];
			
			Parent2 = 5;
			
			do{
				RandomValue = RandMax1();
				if( RandomValue > .5 )
					Parent2= ParentCand[0];
				else if( RandomValue > .25 )
					Parent2= ParentCand[1];
				else if( RandomValue > .125 )
					Parent2= ParentCand[2];
				else
					Parent2= ParentCand[2];
			}while( Parent2 == Parent1);
					
				
			CrossOverGene( &Pop->Genes[Parent1],
						&Pop->Genes[Parent2],
						&Pop->TempGenes[2*i],
						&Pop->TempGenes[(2*i)+1]);

			Pop->TempGenes[2*i].fitness = 0.0;
			Pop->TempGenes[(2*i)+1].fitness = 0.0;
						
			
			MutateGene(&Pop->TempGenes[2*i],Pop->MutationRate,Pop);
			MutateGene(&Pop->TempGenes[(2*i)+1],Pop->MutationRate,Pop);
			
			if(Pop->TempGenes[2*i].fitness <= 0.00001 )
				EvaluateGene(&Pop->TempGenes[2*i],Pop);
			
			if(Pop->TempGenes[(2*i)+1].fitness <= 0.00001 )
				EvaluateGene(&Pop->TempGenes[(2*i)+1],Pop);
		}
	}   
	
	SortGenes(Pop->TempGenes, NumOffSpring);
	
}

void Replacement( Population_Type *Pop , unsigned int NumReplaced )
{
	unsigned int i = 0;
	unsigned int j = 0;
	
	if(0)
	{
		for( i = 1 ; i < Pop->PopulationSize && j < NumReplaced ; i ++ )
		{
			if( Pop->Genes[i].fitness < Pop->TempGenes[j].fitness ||
				i > (Pop->PopulationSize)*2/3)
			{
				unsigned int k = i;
				Gene_Type TempGene1,TempGene2;
				
				TempGene1 = Pop->Genes[i];
				for( k = i ; k< (Pop->PopulationSize-1) ; k++ )
				{
					TempGene2 = Pop->Genes[k+1];
					Pop->Genes[k+1]=TempGene1;
					TempGene1 =TempGene2;
				}
			
				Pop->Genes[i] = Pop->TempGenes[j];
				j++;
			}
		}
		Pop->Genes[Pop->PopulationSize-1] = Pop->Genes[0];
	}
	else
	{
		
		for( j=0 ;  j < NumReplaced; j ++ )
		{
			if( j > Pop->PopulationSize/2 && Pop->Genes[(Pop->PopulationSize - j) - 1].fitness > Pop->TempGenes[j].fitness )
				break;			
			Pop->Genes[(Pop->PopulationSize - j) - 1] = Pop->TempGenes[j];
		}
		
			//Pop->Genes[(Pop->PopulationSize - j) - 1] = Pop->Genes[0];
			//MutateGene(&Pop->Genes[(Pop->PopulationSize - j) - 1],0.5,Pop);
	}
	
	//printf("%d replaced : ",j);

	
}

void RemoveDuplicates( Population_Type* Pop )
{
	int i , j;
	int NoDuplicates = 0;
	
	
	while( NoDuplicates == 0 )
	{
		NoDuplicates = 1;
		for( i = 0 ; i < Pop->PopulationSize ; i ++)
		{
			for( j = i + 1 ; j < Pop->PopulationSize; j++ )
			{
				if( Pop->Genes[i].params[0] == Pop->Genes[j].params[0] &&
					Pop->Genes[i].params[1] == Pop->Genes[j].params[1] &&
					Pop->Genes[i].params[2] == Pop->Genes[j].params[2] &&
					Pop->Genes[i].params[3] == Pop->Genes[j].params[3] &&
					Pop->Genes[i].params[4] == Pop->Genes[j].params[4] )
				{
					NoDuplicates = 0;
					MutateGene(&Pop->Genes[j], 1,Pop);
				}
			}
		}
	}
}
			
void LowArithmeticDrift( Gene_Type* Gene , Population_Type* Pop)
{
	Gene_Type NewGene = *Gene;
	int MutationSize,CurrentParam;
	for( CurrentParam = 0 ; CurrentParam < 5 ; CurrentParam++ )
	{
		MutationSize = GenerateMutationStep(5);
		if( rand() % 2 )
			NewGene.params[CurrentParam] -= MutationSize;
		else
			NewGene.params[CurrentParam] += MutationSize;
					
		if( NewGene.params[CurrentParam] > 100 )
			NewGene.params[CurrentParam] = 100;
		else if( NewGene.params[CurrentParam] < 0 )
			NewGene.params[CurrentParam] = 1;
	}	
		
	
	EvaluateGene(&NewGene,Pop);	
	
	
	*Gene = NewGene;
		
}

/******************************************************************************/
/* Function: BeginEvolution()
 /* Input:  called within main
 /* Output: uses same output functions
 /* Description: Applies genetic algorithm to calibration of SLEUTH
/* See: http://www.geocomputation.org/2011/papers/clarke-lauer.pdf
 
 /******************************************************************************/
int BeginEvolution()
{
	
    int i ;
    int j,low;
    float fitness;
	float StatSum = 0,OldStatSum = 0;
	char filename[100];
	char filename1[100];
	Gene_Type StoredGene;

	
    Population_Type* Pop;
/* Initialize stored gene */
	for (i = 0; i < 5; i++) StoredGene.params[i] = 0;
	StoredGene.fitness = 0.0;
    
    Pop =  initPopulation(PopulationSize,MutationRate);
    
	Pop->MaxEvals = MaxEvals;
	Pop->CurrentEvals = 0;
	
    EvaluatePopulation(Pop);
	SortGenes(Pop->Genes,Pop->PopulationSize);
	sprintf(filename,"SSA_log_%d_%d_%d_%1.2f",PopulationSize,Generations,NumReplaced,MutationRate);
	pFile = fopen(filename,"w");
	sprintf(filename1,"SSA_stats_%d_%d_%d_%1.2f",PopulationSize,Generations,NumReplaced,MutationRate);
	pStatsFile = fopen(filename1,"w");
	int a = 0;
	printf("Population Number %d\n",a % Pop->PopulationSize/Pop->NumSubPopulations);
	//SSA DebugPrintPopulation(Pop);
	i = 0;
	BEST_GENE = Pop->Genes[0];
	while( Pop->MaxEvals > Pop->CurrentEvals )
	{
        printf("****************GENERATION %d *****************\n",i);    
		printf("Stored Gene:");
		DebugPrintGeneOut(StoredGene);
		StatSum = 0;
       //SSA fprintf(pFile,"****************GENERATION %d *****************\n",i+1);
        //printf("****************GENERATION %d *****************\n",i);    
		GenerateOffspring(Pop,NumOffSpring);
		//printf("****************GENERATION %d *****************\n",i);    
		Replacement( Pop , NumReplaced );
		
		//printf("****************GENERATION %d *****************\n",i);    
        MutatePopulation(Pop);
		//printf("****************GENERATION %d *****************\n",i);    
		SortGenes(Pop->Genes,Pop->PopulationSize);
		//printf("****************GENERATION %d *****************\n",i);    
		Pop->Genes[Pop->PopulationSize-2] = BEST_GENE;
		//printf("****************GENERATION %d *****************\n",i);    
		LowArithmeticDrift(&Pop->Genes[Pop->PopulationSize-2],Pop);
		Pop->Genes[Pop->PopulationSize-1] = BEST_GENE;
		/*for( j = 0 ; j < 1;j++ )
		{
			LowArithmeticDrift(&Pop->Genes[j] ,Pop);
		}*/
		RemoveDuplicates(Pop);
		//printf("****************GENERATION %d *****************\n",i);    
		SortGenes(Pop->Genes,Pop->PopulationSize);	
		//printf("****************GENERATION %d *****************\n",i);    
		
		for( j = 0 ; j < Pop->PopulationSize ;j++ )
		{
			StatSum += Pop->Genes[j].fitness;
		}
		
		{// Standard Deviation 
			float avg = StatSum/(float)Pop->PopulationSize;
			float devSquared=0.0;
			float StandardDev;
			int k;
			
			for( k = 0 ; k < Pop->PopulationSize ; k ++ )
			{
				devSquared+= (Pop->Genes[k].fitness - avg) * (Pop->Genes[k].fitness - avg) ;
			}
			
			StandardDev = sqrt(devSquared/(float)Pop->PopulationSize);
			printf("Standard Dev : %f\n",StandardDev);
		}
		
			
		
		/*if( OldStatSum / StatSum > 1.01 )
		{
			OldStatSum = 0 ;
			
			for( j = 2 ; j < Pop->PopulationSize ;j++ )
			{
				MutateGene( &Pop->Genes[j],0.5, Pop); 
			}
			SortGenes(Pop->Genes,Pop->PopulationSize);
			printf("Population Stagnant: Accelerated Evolution Introduced\n");
		}
		else
		{	
			OldStatSum = StatSum;
		}*/
			
		DebugPrintPopulation(Pop);	
		fprintf(pStatsFile,"%d, %f ,%f, ", (i), StatSum/((float)Pop->PopulationSize), Pop->Genes[0].fitness);
		printf("Generation: %d, Average: %f\n", (i++), StatSum/((float)Pop->PopulationSize));
		DebugPrintGeneStat(Pop->Genes[0]);
		DebugPrintPopulationOut(Pop);
		StoredGene = Pop->Genes[0];
		//DebugPrintPopulation(Pop);
    }

	
   
	//    DestroyPopulation(Pop);
    
}

/******************************************************************************
*******************************************************************************
** FUNCTION NAME: main
** PURPOSE:       
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
int
  main (int argc, char *argv[])
{
  char func[] = "main";
  char fname[MAX_FILENAME_LEN];
  char command[5 * MAX_FILENAME_LEN];
  int restart_run = 0;
  RANDOM_SEED_TYPE random_seed;
  int diffusion_coeff;
  int breed_coeff;
  int spread_coeff;
  int slope_resistance;
  int road_gravity;
  int restart_diffusion;
  int restart_breed;
  int restart_spread;
  int restart_slope_resistance;
  int restart_road_gravity;
  time_t tp;
  int evolution = 0;
  char processing_str[MAX_FILENAME_LEN];
  int i;
#ifdef CATCH_SIGNALS
  struct sigaction act, oact;
  tracer = 1;
  act.sa_handler = catch;
  sigemptyset (&act.sa_mask);
  act.sa_flags = 0;
  sigaction (SIGFPE, &act, &oact);
  sigaction (SIGINT, &act, &oact);
  sigaction (SIGSEGV, &act, &oact);
  sigaction (SIGBUS, &act, &oact);
#endif

#ifdef MPI
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &glb_mype);
  MPI_Comm_size (MPI_COMM_WORLD, &glb_npes);
#else
  glb_mype = 0;
  glb_npes = 1;
#endif

  timer_Init ();
  timer_Start (TOTAL_TIME);

  glb_call_stack_index = -1;
  FUNC_INIT;

  /*
   *
   * PARSE COMMAND LINE
   *
   */
/*  if (argc != 9)
  {
    print_usage (argv[0]);
  }
  if (!((strcmp (argv[1], "predict")) ||
      (strcmp (argv[1], "restart")) ||
      (strcmp (argv[1], "test")) ||
      (strcmp (argv[1], "calibrate")) ||
      (strcmp (argv[1], "evolve"))))
  {
    print_usage (argv[0]);
  }
  if (strcmp (argv[1], "predict") == 0)
  {
    proc_SetProcessingType (PREDICTING);
    strcpy (processing_str, "PREDICTING");
  }
  if (strcmp (argv[1], "restart") == 0)
  {
    proc_SetProcessingType (CALIBRATING);
    proc_SetRestartFlag (TRUE);
    strcpy (processing_str, "restart CALIBRATING");
  }
  if (strcmp (argv[1], "test") == 0)
  {
    proc_SetProcessingType (TESTING);
    strcpy (processing_str, "TESTING");
  }
  if (strcmp (argv[1], "calibrate") == 0)
  {
    proc_SetProcessingType (CALIBRATING);
    strcpy (processing_str, "CALIBRATING");
  }
  if (strcmp (argv[1], "evolve") == 0)
  {
    proc_SetProcessingType (EVOLVING);
    strcpy (processing_str, "EVOLVING");
	evolution = 1;
  }
  scen_init (argv[2]);

  PopulationSize = atoi(argv[3]);
  Generations = atoi(argv[4]);
  MutationRate = atof(argv[5]);
  MaxEvals =atoi(argv[6]);
  NumOffSpring =atoi(argv[7]);
  NumReplaced =atoi(argv[8]);
  */
  strcmp (argv[1], "evolve");
  proc_SetProcessingType (EVOLVING);
  strcpy (processing_str, "EVOLVING");
  evolution = 1;

  scen_init("C:/cygwin64/home/sedas/SLEUTHGA_Test/Scenarios/scenario.demo200_calibrate");
  //scen_init("/home/salap/sleuth_ga/Scenarios/scenario.demo200_calibrate");

  PopulationSize = 80;
  Generations = 100;
  MutationRate = 0.13;
  MaxEvals = 200;
  NumOffSpring = 55;
  NumReplaced = 55;

  
    printf("GA constants\n");
	printf("Population size: %d\n", PopulationSize);
	printf("Generations: %d\n", Generations);
	printf("Mutation Rate: %f\n", MutationRate);
	printf("Number of offspring: %d\n", NumOffSpring);
	printf("Number replaced: %d\n", NumReplaced);
  /*
   *
   * SET SOME VARIABLES
   *
   */
  random_seed = scen_GetRandomSeed ();

/*
 * void landclassSetGrayscale (int index, int val);
 * void landclassSetColor (int index, int val);
 * void landclassSetType (int index, char* string);
 * void landclassSetName (int index, char* string);
 * void landclassSetNumClasses (int val);
 * int scen_GetNumLanduseClasses ();
 * char* scen_GetLanduseClassName (int);
 * char* scen_GetLanduseClassType (int);
 * int scen_GetLanduseClassColor (int);
 * int scen_GetLanduseClassGrayscale (int i);
 * 
 */
  landclassSetNumClasses (scen_GetNumLanduseClasses ());
  for (i = 0; i < scen_GetNumLanduseClasses (); i++)
  {
    landclassSetGrayscale (i, scen_GetLanduseClassGrayscale (i));
    landclassSetName (i, scen_GetLanduseClassName (i));
    landclassSetType (i, scen_GetLanduseClassType (i));
    landclassSetColor (i, scen_GetLanduseClassColor (i));
  }

  /*
   *
   * SET UP COEFFICIENTS
   *
   */
  if (strcmp (argv[1], "restart") == 0)
  {
    if (scen_GetLogFlag ())
    {
      scen_Append2Log ();
      if (scen_GetLogFP ())
      {
        fprintf (scen_GetLogFP (), "%s %u Reading restart file\n",
                 __FILE__, __LINE__);
      }
      scen_CloseLog ();
    }
    inp_read_restart_file (&restart_diffusion,
                           &restart_breed,
                           &restart_spread,
                           &restart_slope_resistance,
                           &restart_road_gravity,
                           &random_seed,
                           &restart_run);
    proc_SetCurrentRun (restart_run);

  }
  else
  {
    proc_SetCurrentRun (0);
  }
  coeff_SetStartDiffusion (scen_GetCoeffDiffusionStart ());
  coeff_SetStartSpread (scen_GetCoeffSpreadStart ());
  coeff_SetStartBreed (scen_GetCoeffBreedStart ());
  coeff_SetStartSlopeResist (scen_GetCoeffSlopeResistStart ());
  coeff_SetStartRoadGravity (scen_GetCoeffRoadGravityStart ());

  coeff_SetStopDiffusion (scen_GetCoeffDiffusionStop ());
  coeff_SetStopSpread (scen_GetCoeffSpreadStop ());
  coeff_SetStopBreed (scen_GetCoeffBreedStop ());
  coeff_SetStopSlopeResist (scen_GetCoeffSlopeResistStop ());
  coeff_SetStopRoadGravity (scen_GetCoeffRoadGravityStop ());

  coeff_SetStepDiffusion (scen_GetCoeffDiffusionStep ());
  coeff_SetStepSpread (scen_GetCoeffSpreadStep ());
  coeff_SetStepBreed (scen_GetCoeffBreedStep ());
  coeff_SetStepSlopeResist (scen_GetCoeffSlopeResistStep ());
  coeff_SetStepRoadGravity (scen_GetCoeffRoadGravityStep ());

  coeff_SetBestFitDiffusion (scen_GetCoeffDiffusionBestFit ());
  coeff_SetBestFitSpread (scen_GetCoeffSpreadBestFit ());
  coeff_SetBestFitBreed (scen_GetCoeffBreedBestFit ());
  coeff_SetBestFitSlopeResist (scen_GetCoeffSlopeResistBestFit ());
  coeff_SetBestFitRoadGravity (scen_GetCoeffRoadGravityBestFit ());

  /*
   *
   * INITIALIZE IGRID
   *
   */
  igrid_init ();

  /*
   *
   * PRINT BANNER
   *
   */
  if (scen_GetEchoFlag ())
  {
    out_banner (stdout);
  }
  if (scen_GetLogFlag ())
  {
    scen_Append2Log ();
    out_banner (scen_GetLogFP ());
    scen_CloseLog ();
  }

  /*
   *
   * LOG SOME STUFF
   *
   */
  if (scen_GetLogFlag ())
  {
    scen_Append2Log ();
    time (&tp);
    fprintf (scen_GetLogFP (), "DATE OF RUN: %s\n",
             asctime (localtime (&tp)));
    fprintf (scen_GetLogFP (), "USER: %s\n", getenv ("USER"));
    fprintf (scen_GetLogFP (), "HOST: %s\n", getenv ("HOST"));
    fprintf (scen_GetLogFP (), "HOSTTYPE: %s\n", getenv ("HOSTTYPE"));
    fprintf (scen_GetLogFP (), "OSTYPE: %s\n", getenv ("OSTYPE"));
    fprintf (scen_GetLogFP (), "Type of architecture: %u bit\n\n",
             BYTES_PER_WORD*8);
    fprintf (scen_GetLogFP (), "Number of CPUs %u \n\n",
             glb_npes);
    fprintf (scen_GetLogFP (), "PWD: %s\n", getenv ("PWD"));
    fprintf (scen_GetLogFP (), "Scenario File: %s\n",
             scen_GetScenarioFilename ());
    fprintf (scen_GetLogFP (), "Type of Processing: %s\n",
             processing_str);
    fprintf (scen_GetLogFP (), "\n\n");

    scen_echo (scen_GetLogFP ());
    coeff_LogStart (scen_GetLogFP ());
    coeff_LogStop (scen_GetLogFP ());
    coeff_LogStep (scen_GetLogFP ());
    coeff_LogBestFit (scen_GetLogFP ());
    scen_CloseLog ();
    tracer = 2;
  }

  /*
   *
   * SET UP FLAT MEMORY
   *
   */
#ifdef MPI
  MPI_Barrier (MPI_COMM_WORLD);
#endif
  mem_Init ();
  if (scen_GetLogFlag ())
  {
    scen_Append2Log ();
    mem_LogPartition (scen_GetLogFP ());
    mem_CheckMemory (scen_GetLogFP (), __FILE__, func, __LINE__);
    scen_CloseLog ();
  }
#ifdef MPI
  MPI_Barrier (MPI_COMM_WORLD);
#endif

  /*
   *
   * INITIALIZE LANDUSE
   *
   */
  if (scen_GetDoingLanduseFlag ())
  {
    landclass_Init ();
    if (scen_GetLogLandclassSummaryFlag ())
    {
      if (scen_GetLogFlag ())
      {
        scen_Append2Log ();
        landclass_LogIt (scen_GetLogFP ());
        scen_CloseLog ();
      }
    }
  }

  /*
   *
   * INITIALIZE COLORTABLES
   *
   */
  color_Init ();

  /*
   *
   * WRITE MEMORY MAPS
   *
   */
  if (scen_GetLogMemoryMapFlag ())
  {
    color_MemoryLog (mem_GetLogFP ());
    coeff_MemoryLog (mem_GetLogFP ());
    timer_MemoryLog (mem_GetLogFP ());
    igrid_MemoryLog (mem_GetLogFP ());
    pgrid_MemoryLog (mem_GetLogFP ());
    stats_MemoryLog (mem_GetLogFP ());
    mem_MemoryLog (mem_GetLogFP ());
    proc_MemoryLog (mem_GetLogFP ());
    landclass_MemoryLog (mem_GetLogFP ());
    scen_MemoryLog (mem_GetLogFP ());
    trans_MemoryLog (mem_GetLogFP ());
    mem_CloseLog ();
  }

  /*
   *
   * READ INPUT DATA FILES
   *
   */
  igrid_ReadFiles ();
  if (scen_GetLogFlag ())
  {
    scen_Append2Log ();
    igrid_ValidateGrids (scen_GetLogFP ());
    scen_CloseLog ();
  }
  else
  {
    igrid_ValidateGrids (NULL);
  }
  igrid_NormalizeRoads ();
  if (scen_GetLogFlag ())
  {
    scen_Append2Log ();
    igrid_LogIt (scen_GetLogFP ());
    igrid_VerifyInputs (scen_GetLogFP ());
    scen_CloseLog ();
  }
  else
  {
    igrid_VerifyInputs (NULL);
  }

  /*
   *
   * INITIALIZE THE PGRID GRIDS
   *
   */
  pgrid_Init ();

  if (scen_GetLogFlag ())
  {
    if (scen_GetLogColortablesFlag ())
    {
      scen_Append2Log ();
      color_LogIt (scen_GetLogFP ());
      scen_CloseLog ();
    }
  }

  /*
   *
   * COUNT THE NUMBER OF RUNS
   *
   */
  proc_SetTotalRuns ();
  if (scen_GetLogFlag ())
  {
    if (proc_GetProcessingType () == CALIBRATING)
    {
      scen_Append2Log ();
      fprintf (scen_GetLogFP (), "%s %u Total Number of Runs = %u\n",
               __FILE__, __LINE__, proc_GetTotalRuns ());
      scen_CloseLog ();
    }
  }

  proc_SetLastMonteCarlo (scen_GetMonteCarloIterations () - 1);
  /*
   *
   * COMPUTE THE TRANSITION MATRIX
   *
   */
  if (scen_GetDoingLanduseFlag ())
  {
    trans_Init ();
    if (scen_GetLogFlag ())
    {
      if (scen_GetLogTransitionMatrixFlag ())
      {
        scen_Append2Log ();
        trans_LogTransition (scen_GetLogFP ());
        scen_CloseLog ();
      }
    }
  }
  /*
   *
   * COMPUTE THE BASE STATISTICS AGAINST WHICH CALIBRATION WILL TAKE PLACE
   *
   */
  stats_Init ();
  if (scen_GetLogFlag ())
  {
    if (scen_GetLogBaseStatsFlag ())
    {
      scen_Append2Log ();
      stats_LogBaseStats (scen_GetLogFP ());
      scen_CloseLog ();
    }
  }
  if (scen_GetLogFlag ())
  {
    if (scen_GetLogDebugFlag ())
    {
      scen_Append2Log ();
      igrid_Debug (scen_GetLogFP (), __FILE__, __LINE__);
      scen_CloseLog ();
    }
  }

  proc_SetNumRunsExecThisCPU (0);
  if (proc_GetCurrentRun () == 0 && glb_mype == 0)
  {
    if (proc_GetProcessingType () != PREDICTING)
    {
      sprintf (fname, "%scontrol_stats.log", scen_GetOutputDir ());
      stats_CreateControlFile (fname);
    }
    if (scen_GetWriteStdDevFileFlag ())
    {
      sprintf (fname, "%sstd_dev.log", scen_GetOutputDir ());
      stats_CreateStatsValFile (fname);
    }
    if (scen_GetWriteAvgFileFlag ())
    {
      sprintf (fname, "%savg.log", scen_GetOutputDir ());
      stats_CreateStatsValFile (fname);
    }
  }

  coeff_CreateCoeffFile ();

#ifdef MPI
  MPI_Barrier (MPI_COMM_WORLD);
#endif
    proc_SetStopYear (igrid_GetUrbanYear (igrid_GetUrbanCount () - 1));


              sprintf (fname, "%s%s%u", scen_GetOutputDir (),
                       RESTART_FILE, glb_mype);
              out_write_restart_data (fname,
                                      diffusion_coeff,
                                      breed_coeff,
                                      spread_coeff,
                                      slope_resistance,
                                      road_gravity,
                                      scen_GetRandomSeed (),
                                      restart_run);

      
    if (evolution) {
    	printf("Evolution begins...\n");
		BeginEvolution();

	}
  

#ifdef MPI
  MPI_Barrier (MPI_COMM_WORLD);
#endif
  if (glb_mype == 0)
  {
    if (scen_GetWriteCoeffFileFlag())
    {
      coeff_ConcatenateFiles();
    }
    if (scen_GetWriteAvgFileFlag())
    {
      stats_ConcatenateAvgFiles();
    }
    if (scen_GetWriteStdDevFileFlag())
    {
      stats_ConcatenateStdDevFiles();
    }
    if (proc_GetProcessingType() != PREDICTING)
    {
      stats_ConcatenateControlFiles();
    }
  }

  if (scen_GetPostprocessingFlag ())
  {
#ifdef MPI
    MPI_Barrier (MPI_COMM_WORLD);
#endif
    if (glb_mype == 0)
    {
      if (strlen (scen_GetWhirlgifBinary()) > 0)
      {
        if (scen_GetViewDeltatronAgingFlag())
        {
          sprintf (command,
            "%s -time 100 -o %sanimated_deltatron.gif %sdeltatron_*.gif",
                   scen_GetWhirlgifBinary(), scen_GetOutputDir (),
                   scen_GetOutputDir());
          system (command);
        }
        if (scen_GetViewGrowthTypesFlag ())
        {
          sprintf (command,
                   "%s -time 100 -o %sanimated_z_growth.gif %sz_growth_types_*.gif",
                   scen_GetWhirlgifBinary (), scen_GetOutputDir (),
                   scen_GetOutputDir ());
          system (command);
        }
        if (proc_GetProcessingType () != CALIBRATING)
        {
          if (scen_GetDoingLanduseFlag ())
          {
            sprintf (command,
                     "%s -time 100 -o %sanimated_land_n_urban.gif %s*_land_n_urban*.gif",
                     scen_GetWhirlgifBinary (), scen_GetOutputDir (),
                     scen_GetOutputDir ());
            system (command);
          }
          else
          {
            sprintf (command,
                  "%s -time 100 -o %sanimated_urban.gif %s*_urban_*.gif",
                     scen_GetWhirlgifBinary (),
                     scen_GetOutputDir (), scen_GetOutputDir ());
            system (command);
          }
        }
      }
    }
  }
#ifdef MPI
  MPI_Finalize ();
#endif
  timer_Stop (TOTAL_TIME);

  if (scen_GetLogFlag ())
  {
    scen_Append2Log ();
    if (scen_GetLogTimingsFlag () > 0)
    {
      timer_LogIt (scen_GetLogFP ());
    }
    mem_LogMinFreeWGrids (scen_GetLogFP ());
    scen_CloseLog ();
  }
  return (0);
}

/******************************************************************************
*******************************************************************************
** FUNCTION NAME: print_usage
** PURPOSE:       help the user
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
static void
  print_usage (char *binary)
{
  printf ("Usage:\n");
  printf ("%s <mode> <scenario file>\n", binary);
  printf ("Allowable modes are:\n");
  printf ("  calibrate\n");
  printf ("  restart\n");
  printf ("  test\n");
  printf ("  predict\n");
  printf ("  evolve\n");
  EXIT (1);
}
#ifdef CATCH_SIGNALS

/******************************************************************************
*******************************************************************************
** FUNCTION NAME: catch
** PURPOSE:       catch signals
** AUTHOR:        Keith Clarke
** PROGRAMMER:    Tommy E. Cathey of NESC (919)541-1500
** CREATION DATE: 11/11/1999
** DESCRIPTION:
**
**
*/
void
  catch (int signo)
{
  int i;
  if (tracer < 2)
  {
    printf ("%s %u pe: %u %s\n", __FILE__, __LINE__, glb_mype,
            "Please make sure the following env variables are defined");
    printf ("%s %u pe: %u %s\n", __FILE__, __LINE__, glb_mype,
            "USER -- set to your username");
    printf ("%s %u pe: %u %s\n", __FILE__, __LINE__, glb_mype,
            "HOST -- set to your machines's name = uname -n");
    printf ("%s %u pe: %u %s\n", __FILE__, __LINE__, glb_mype,
          "HOSTTYPE -- set to your machines's type such as Sparc, Cray");
    printf ("%s %u pe: %u %s\n", __FILE__, __LINE__, glb_mype,
          "OSTYPE -- set to your machines's OS type such as Solaris2.7");
    printf ("%s %u pe: %u %s\n", __FILE__, __LINE__, glb_mype,
            "PWD -- set to your current working directory");
  }
  if (signo == SIGBUS)
  {
    printf ("%s %u caught signo SIGBUS : bus error\n",
            __FILE__, __LINE__);
    printf ("Currently executing function: %s\n",
            glb_call_stack[glb_call_stack_index]);
    for (i = glb_call_stack_index; i >= 0; i--)
    {
      printf ("glb_call_stack[%d]= %s\n", i, glb_call_stack[i]);
    }
    EXIT (1);
  }
  if (signo == SIGSEGV)
  {
    printf ("%s %u pe: %u caught signo SIGSEGV : Invalid storage access\n",
            __FILE__, __LINE__, glb_mype);
    printf ("Currently executing function: %s\n",
            glb_call_stack[glb_call_stack_index]);
    for (i = glb_call_stack_index; i >= 0; i--)
    {
      printf ("glb_call_stack[%d]= %s\n", i, glb_call_stack[i]);
    }
    EXIT (1);
  }
  if (signo == SIGINT)
  {
    printf ("%s %u caught signo SIGINT : Interrupt\n",
            __FILE__, __LINE__);
    printf ("Currently executing function: %s\n",
            glb_call_stack[glb_call_stack_index]);
    for (i = glb_call_stack_index; i >= 0; i--)
    {
      printf ("glb_call_stack[%d]= %s\n", i, glb_call_stack[i]);
    }
    EXIT (1);
  }
  if (signo == SIGFPE)
  {
    printf ("%s %u caught signo SIGFPE : Floating-point exception\n",
            __FILE__, __LINE__);
    printf ("Currently executing function: %s\n",
            glb_call_stack[glb_call_stack_index]);
    for (i = glb_call_stack_index; i >= 0; i--)
    {
      printf ("glb_call_stack[%d]= %s\n", i, glb_call_stack[i]);
    }
    EXIT (1);
  }
  printf ("%s %u caught signo %d\n", __FILE__, __LINE__, signo);
  EXIT (1);
}
#endif
