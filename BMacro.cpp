// #include "Decay.h"

#include <vector>
#include "../git/g-2Simulation/src/Utilities/Vector3D.h"
#include <cmath>
#include <iostream>

void PositionToGrid(const Position3D * position, int & xindex, int & yindex, int & zindex);
void GridToPosition(Position3D * position, const int & xindex, const int & yindex, const int & zindex);

// B Decay::BGen::operator()(Position3D p) const 
//int main()
int BMacro()
{
  bool verboso = true;				//necessario per il debugging, fa capire a che punto si è e evita di bloccarsi nei loop
  bool co = true;				//computational optimization: se questo parametro è settato su true, allora si ignorerà in pratica la componente x
						//dei campi vettoriali, perché non necessaria. Ciò è possibile se non si hanno correnti lungo x ed è un asse di simmetri
						//per il sistema prima della riduzione delle dimensioni
  
  using namespace std;
  
  //Campo magnetico
  vector<Vector3D>   fieldLayer2;		//identifica la z una volta note x e y
  vector<Vector3D*>  fieldLayer1;		//identifica la y una volta nota la x
  vector<Vector3D**> field;			//identificaa la x (asse)
  
  //Potenziale vettore
  vector<Vector3D>   potentialLayer2;
  vector<Vector3D*>  potentialLayer1;
  vector<Vector3D**> potential;

  //Correnti
  vector<Vector3D>   currentsLayer2;
  vector<Vector3D*>  currentsLayer1;
  vector<Vector3D**> currents;
  
  
  //impostazione delle dimensioni del sistema
  double xsyst = 10.;
  double ysyst = 10.;
  double zsyst = 10.;
  
  //impostazione delle dimensioni del solenoide (funziona solo con parallelepipedo non ruotato)
  double xsol = 1.;
  double ysol = 1.;
  double zsol = 1.;
  
  //impostazione della griglia, indica in pratica lo step tra due punti della griglia stessa
  double delta = 0.1;
  
  //impostazione della corrente, indica in pratica la corrente che viaggia sul filo usato (si considera costante)
  double current = 10.;	//mA
  
  //impostazione della lunghezza del filo
//   double lenght = 1000. * 1000; //m, non mm, per quello il fattore 1000
  
  //calcolo dei punti per il processo numerico
  int xsteps = int(xsyst/(delta*2)) + 1; 
  int ysteps = int(ysyst/(delta*2)) + 1; 
  int zsteps = int(zsyst/(delta*2)) + 1;
  
  if(verboso) cout <<"xsteps = " << xsteps << " ysteps = " << ysteps << " zsteps = " << zsteps << endl;
  
  //creazione delle griglie di lavoro tridimentionali e loro inizializzazione: RIDUZIONE SISTEMA A 1/8
  cout << "Inizio creazione griglie..." << endl;
  fieldLayer2.resize(xsteps * ysteps * zsteps, {0.,0.,0.});	//contiene tutto il campo
  fieldLayer1.resize(xsteps * ysteps);				//contiene puntatori a un piano orizzontale
  field.resize(xsteps);						//contiene puntatori alla posizione lungo l'asse del solenoide

  potentialLayer2.resize(xsteps * ysteps * zsteps, {0.,0.,0.});	//contiene tutto il campo
  potentialLayer1.resize(xsteps * ysteps);			//contiene puntatori a un piano orizzontale
  potential.resize(xsteps);					//contiene puntatori alla posizione lungo l'asse del solenoide

  currentsLayer2.resize(xsteps * ysteps * zsteps, {0.,0.,0.});	//contiene tutto il campo
  currentsLayer1.resize(xsteps * ysteps);			//contiene puntatori a un piano orizzontale
  currents.resize(xsteps);					//contiene puntatori alla posizione lungo l'asse del solenoide
  
  for(unsigned int i = 0; i < fieldLayer1.size(); i++)	//sistema i puntatori in modo che puntino al posto giusto
  {
    fieldLayer1.at(i) = & fieldLayer2.at(i * ysteps);
    potentialLayer1.at(i) = & potentialLayer2.at(i * ysteps);
    currentsLayer1.at(i) = & currentsLayer2.at(i * ysteps);
  }
  for(unsigned int i = 0; i < field.size(); i++)
  {
    field.at(i) = & fieldLayer1.at(i * xsteps);
    potential.at(i) = & potentialLayer1.at(i * xsteps);
    currents.at(i) = & currentsLayer1.at(i * xsteps);
  }
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  double perimeter = (2 * ysol) + (2 * zsol);	//occhio a x o y 
  double current_density = current/perimeter;	//sarà la densità di corrente lungo x
  
  //inserimento delle correnti
  if (verboso) cout << "Inizio inserimento delle correnti..." << endl;
  int xindex = 0;
  int yindex = 0;
  int zindex = 0;
  Position3D starting_point = {0., ysol*.5, 0.};
  Position3D walking_point = starting_point;
  Vector3D * modify_inloop;
  
  
  int x_border = 0;			//variabili utili per tornare al posto giusto a fine ciclo
  int y_border = 0;
  int z_border = 0;
  
//   PositionToGrid(& solenoid_starting_point, x_border, y_border, z_border);
  
  unsigned int debug1 = 0;
  unsigned int debug2 = 0;
  unsigned int debug3 = 0;
  
  while(walking_point.x < xsol*.5)		//ciclo sulla lunghezza del solenoide
  {
    walking_point.z = starting_point.z;		//tornare alla posizione di partenza in z e y
    walking_point.y = starting_point.y;
    while(walking_point.z < zsol*.5)		//salita
    {
     PositionToGrid(& walking_point, xindex, yindex, zindex);	//trova il valore nella grid della posizione
//      if(verboso) cout << "xindex = " << xindex << " yindex = " << yindex << " zindex = " << zindex << endl;
     modify_inloop = & currents[xindex][yindex][zindex];	//si fissa nel punto desiderato
     modify_inloop->z = current_density;			//cambia la corrente in tale punto
     zindex++;							//avanza lungo la griglia
     GridToPosition(& walking_point, xindex, yindex, zindex);	//passa da griglia a posizione
     debug1++;
     if( debug1 > 100000 && verboso)
     {
       cout << "ERROR IN CYCLE1" << endl;
       return 0;
     }
    }
//     if(verboso) cout << "FINE PRIMO CICLO CORRENTE" << endl;


    
    debug1 =0;
    while(walking_point.y > 0.)		//tratto orizzontale
    {
      
      PositionToGrid(& walking_point, xindex, yindex, zindex);
      modify_inloop = & currents[xindex][yindex][zindex];
      modify_inloop->y = -current_density;
      
      if(debug1 == 0)				//mette l'angolino con corrente a 45 gradi     
      {
	modify_inloop->y = -current_density*.5;
	modify_inloop->z = current_density*.5;
      }
      
      yindex--;
      GridToPosition(& walking_point, xindex, yindex, zindex);
      debug1++;
      if( debug1 > 100000 && verboso)
      {
	cout << "ERROR IN CYCLE2" << endl;
	return 0;
      }
    }
    debug1 = 0;
    xindex++;					//avanza lungo l'asse
    
    GridToPosition(& walking_point, xindex, yindex, zindex);
    debug2++;
    if( debug2 > 100000 && verboso)
    {
      cout << "ERROR IN CYCLE3" << endl;
      return 0;
    }
//     if(verboso) cout <<"Step done" << endl;
  }
  
  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //si ricorda che si sfrutta la simmetria del problema e si studia solo un quarto del sistema,
  //cioè la parte che va da z = 0 a z = zsyst * .5 e da y = 0 a y = ysyst *.5  e da x = 0 a x = xsyst*.5

  
  //inizializzazione del potenziale vettore
  if (verboso) cout << "Inizio inizializzazione potenziale vettore..." << endl;
  starting_point  = {0.,0.,0.};		//inizia dal centro del sistema
  walking_point   = starting_point;
  double rho = 0;
  double theta = 0;
  
  debug1 = 0;
  debug2 = 0;
  debug3 = 0;
  
//   PositionToGrid(& walking_point, x_border, y_border, z_border);
  
  //set del potenziale vettore all'interno del solenoide
  double mu_0 = 100.;
  while(walking_point.x < xsol*.5)		//ciclo sulla lunghezza del solenoide
  {
    walking_point.z = starting_point.z;
    while(walking_point.z < zsol*.5)		//ciclo sull'altezza del solenoide
    {
      walking_point.y = starting_point.y;
      while(walking_point.y < ysol*.5)		//riempie striscia di larghezza
      {
	PositionToGrid(& walking_point, xindex, yindex, zindex);
	modify_inloop = & potential[xindex][yindex][zindex];
	rho = mu_0 * .5 * sqrt(walking_point.x * walking_point.x + walking_point.y * walking_point.y);
	if(walking_point.z == 0) theta = 0;				//serve a evitare infiniti
	else theta = atan(walking_point.y/walking_point.z);		//ordine giusto? credo di sì, basta mantenere quello sopra qua come quello con seno quando riempio
	* modify_inloop = {0., rho * sin(theta), rho * cos(theta)};
// 	if(verboso) cout << modify_inloop->y << " " << modify_inloop->z << endl;
	yindex++;
	GridToPosition(& walking_point, xindex, yindex, zindex);
	debug1++;
	if( debug1 > 100000 && verboso)
	{
	  cout << "ERROR IN CYCLE4" << endl;
	  return 0;
	}
      }
      zindex++;					//avanza lungo l'altezza
      GridToPosition(& walking_point, xindex, yindex, zindex);
      debug2++;
      if( debug2 > 100000 && verboso)
      {
	cout << "ERROR IN CYCLE5" << endl;
	return 0;
      }
    }
    xindex++;					//avanza lungo l'asse
    GridToPosition(& walking_point, xindex, yindex, zindex);
    debug1++;
    if( debug3 > 100000 && verboso)
    {
      cout << "ERROR IN CYCLE6" << endl;
      return 0;
    }
  }
  
  
  debug1 = 0;
  debug2 = 0;
  debug3 = 0;
  
  //set del potenziale vettore all'esterno del solenoide !!SOPRA!!
  if (verboso) cout << "Inizio inizializzazione esterna al poteziale vettore..." << endl;
  starting_point = {0., 0., zsol * .5};		//inizia da sopra il solenoide
  walking_point  = starting_point;
  
//   PositionToGrid(& walking_point, x_border, y_border, z_border);
  
  
  while(walking_point.x < xsol*.5)		//ciclo sulla lunghezza del solenoide
  {
    walking_point.z = starting_point.z;
    while(walking_point.z < zsyst*.5)		//ciclo sull'altezza del sistema
    {
      walking_point.y = starting_point.y;
      while(walking_point.y < ysol*.5)		//riempie striscia di larghezza del solenoide
      {
	PositionToGrid(& walking_point, xindex, yindex, zindex);
	modify_inloop = & potential[xindex][yindex][zindex];
	rho = mu_0 * .5 * zsol * .5 * zsol * .5 / (sqrt(walking_point.x * walking_point.x + walking_point.y * walking_point.y)+0.0001);
	if(walking_point.z == 0) theta = 0;				//serve a evitare infiniti
	else theta = atan(walking_point.y/walking_point.z);		//ordine giusto? credo di sì, basta mantenere quello sopra qua come quello con seno quando riempio
	* modify_inloop = {0., rho * sin(theta), rho * cos(theta)};
// 	if(verboso) cout << modify_inloop->y << " " << modify_inloop->z << endl;
	yindex++;
	GridToPosition(& walking_point, xindex, yindex, zindex);
	debug1++;
	if( debug1 > 100000 && verboso)
	{
	  cout << "ERROR IN CYCLE7" << endl;
	  return 0;
	}
      }
      zindex++;					//avanza lungo l'altezza
      GridToPosition(& walking_point, xindex, yindex, zindex);
      debug2++;
      if( debug2 > 100000 && verboso)
      {
	cout << "ERROR IN CYCLE8" << endl;
	return 0;
      }
    }
    xindex++;					//avanza lungo l'asse
    GridToPosition(& walking_point, xindex, yindex, zindex);
    debug3++;
    if( debug3 > 100000 && verboso)
    {
      cout << "ERROR IN CYCLE9" << endl;
      return 0;
    }
  }
  
  debug1 = 0;
  debug2 = 0;
  debug3 = 0;
  
  //set del potenziale vettore all'esterno del solenoide !!DI LATO!!
  
  
  starting_point = {0., ysol * .5, 0};		//inizia da di lato al solenoide
  walking_point = starting_point;
//   PositionToGrid(& walking_point, x_border, y_border, z_border);
  
  
  while(walking_point.x < xsol*.5)		//ciclo sulla lunghezza del solenoide
  {
    walking_point.z = starting_point.z;
    while(walking_point.z < zsol*.5)		//ciclo sull'altezza del solenoide
    {
      walking_point.y = starting_point.y;
      while(walking_point.y < ysyst*.5)		//riempie striscia di larghezza del sistema
      {
	PositionToGrid(& walking_point, xindex, yindex, zindex);
	modify_inloop = & potential[xindex][yindex][zindex];
	rho = mu_0 * .5 * ysol * .5 * ysol * .5 / (sqrt(walking_point.x * walking_point.x + walking_point.y * walking_point.y) + 0.0001);
	if(walking_point.z == 0) theta = 0;				//serve a evitare infiniti
	theta = atan(walking_point.y/walking_point.z);		//ordine giusto? credo di sì, basta mantenere quello sopra qua come quello con seno quando riempio
	* modify_inloop = {0., rho * sin(theta), rho * cos(theta)};
// 	if(verboso) cout << modify_inloop->y << " " << modify_inloop->z << endl;
	yindex++;
	GridToPosition(& walking_point, xindex, yindex, zindex);
	debug1++;
	if( debug1 > 100000 && verboso)
	{
	  cout << "ERROR IN CYCLE10"  << endl;
	  cout << "yindex = " << yindex << " ysol*.5 = " << ysol *.5 << " ysyst = " << ysyst << endl;
	  return 0;
	}
      }
      zindex++;					//avanza lungo l'altezza
      GridToPosition(& walking_point, xindex, yindex, zindex);
      debug2++;
      if( debug2 > 100000 && verboso)
      {
	cout << "ERROR IN CYCLE11" << endl;
	return 0;
      }
    }
    xindex++;					//avanza lungo l'asse
    GridToPosition(& walking_point, xindex, yindex, zindex);
    debug3++;
    if( debug3 > 100000 && verboso)
    {
      cout << "ERROR IN CYCLE12" << endl;
      return 0;
    }
  }
  
  debug1 = 0;
  debug2 = 0;
  debug3 = 0;
  
  //set del potenziale vettore all'esterno del solenoide !!ACCANTO!!
  
  
  starting_point = {0., ysol * .5, zsol * .5};	//inizia dal'angolino del solenoide
  walking_point = starting_point;
//   PositionToGrid(& walking_point, x_border, y_border, z_border);
  
  while(walking_point.x < xsol*.5)		//ciclo sulla lunghezza del solenoide
  {
    walking_point.z = starting_point.z;
    while(walking_point.z < zsyst*.5)		//ciclo sull'altezza del sistema
    {
      walking_point.y = starting_point.y;
      while(walking_point.y < ysyst*.5)		//riempie striscia di larghezza del sistema
      {
	PositionToGrid(& walking_point, xindex, yindex, zindex);
	modify_inloop = & potential[xindex][yindex][zindex];
	rho = mu_0 * .5 * zsol * .5 * ysol * .5 / (sqrt(walking_point.x * walking_point.x + walking_point.y * walking_point.y) + 0.0001);	//ibrido tra il sopra e il di lato
	if(walking_point.z == 0) theta = 0;				//serve a evitare infiniti
	else theta = atan(walking_point.y/walking_point.z);		//ordine giusto? credo di sì, basta mantenere quello sopra qua come quello con seno quando riempio
	* modify_inloop = {0., rho * sin(theta), rho * cos(theta)};
// 	if(verboso) cout << modify_inloop->y << " " << modify_inloop->z << endl;
	yindex++;
	GridToPosition(& walking_point, xindex, yindex, zindex);
	debug1++;
	if( debug1 > 100000 && verboso)
	{
	  cout << "ERROR IN CYCLE13" << endl;
	  return 0;
	}
      }
      zindex++;					//avanza lungo l'altezza
      GridToPosition(& walking_point, xindex, yindex, zindex);
      debug2++;
      if( debug2 > 100000 && verboso)
      {
	cout << "ERROR IN CYCLE14" << endl;
	return 0;
      }
    }
    xindex++;					//avanza lungo l'asse
    GridToPosition(& walking_point, xindex, yindex, zindex);
    debug3++;
    if( debug3 > 100000 && verboso)
    {
      cout << "ERROR IN CYCLE15" << endl;
      return 0;
    }
  } 
  //il resto del potenziale vettore è alsciato a zero, verrà poi modificato dall'algoritmo numerico.
  
  
  
  
  
  
  //calcolo del potenziale vettore utilizzando l'algoritmo di Jacobi, usando i due campi field e potential
  //le condizioni al contorno sono miste: quando il punto è alla fine dell'intero sistema allora il potenziale vettore è fissato a zero,
  //quando invece si è nelle "pareti fittizie" (quelle date dalla divisione per simmetria del problema) allora sono riflettenti.
 long ninteractions = 1000;

//  Position3D numeric_start = {0,0,0};
  

//lavoro nel discreto, ignoro il continuo

// PositionToGrid(& numeric_start, xindex, yindex, zindex);
// int x_border = xindex;
// int y_border = yindex;
// int z_border = zindex;
 x_border = 0;
 y_border = 0;
 z_border = 0;

double analyzedpoint;
double beforex, beforey, beforez;
double beyondx, beyondy, beyondz;
double leftx, lefty, leftz;
double rightx, righty, rightz;
double belowx, belowy, belowz;
double abovex, abovey, abovez;

double zero = 0;
double currentj;

if(verboso) cout << "Inizio computazione del potenziale vettore..." << endl;
for(unsigned int i = 0; i < ninteractions; i++)
{
  
  xindex = x_border;
  yindex = y_border;
  zindex = z_border;
  
  while(xindex < xsteps)
  {
    yindex = y_border;
    while(yindex < ysteps)
    {
      zindex = z_border;
      while(zindex < zsteps)
      {
	if(xindex == x_border)                     
	{
	  beforex = potential[xindex+1][yindex][zindex].x;	//riflessione
	  beforey = potential[xindex+1][yindex][zindex].y;	//riflessione
	  beforez = potential[xindex+1][yindex][zindex].z;	//riflessione
	} else
	{
	  beforex = potential[xindex-1][yindex][zindex].x;	//bulk
	  beforey = potential[xindex-1][yindex][zindex].y;	//bulk
	  beforez = potential[xindex-1][yindex][zindex].z;	//bulk
	}
	if(xindex == xsteps - 1)            
	{
	  beyondx = zero;					//nullo
	  beyondy = zero;					//nullo
	  beyondz = zero;					//nullo
	} else 
	{
	  beyondx = potential[xindex+1][yindex][zindex].x;	//bulk
	  beyondy = potential[xindex+1][yindex][zindex].y;	//bulk
	  beyondz = potential[xindex+1][yindex][zindex].z;	//bulk
	}
	if(yindex == y_border)   
	{
	  leftx = potential[xindex][yindex+1][zindex].x;	//riflessione
	  lefty = potential[xindex][yindex+1][zindex].y;	//riflessione
	  leftz = potential[xindex][yindex+1][zindex].z;	//riflessione
	} else
	{
	  leftx = potential[xindex][yindex-1][zindex].x;	//bulk
	  lefty = potential[xindex][yindex-1][zindex].y;	//bulk
	  leftz = potential[xindex][yindex-1][zindex].z;	//bulk
	}
	if(yindex == ysteps - 1)      
	{
	  rightx = zero;					//nullo
	  righty = zero;					//nullo
	  rightz = zero;					//nullo
	} else 
	{
	  rightx = potential[xindex][yindex+1][zindex].x;	//bulk
	  righty = potential[xindex][yindex+1][zindex].y;	//bulk
	  rightz = potential[xindex][yindex+1][zindex].z;	//bulk
	}
	if(zindex == z_border)    
	{
	  belowx = potential[xindex][yindex][zindex+1].x;	//riflessione
	  belowy = potential[xindex][yindex][zindex+1].y;	//riflessione
	  belowz = potential[xindex][yindex][zindex+1].z;	//riflessione
	} else 
	{
	  belowx = potential[xindex][yindex][zindex-1].x;	//bulk
	  belowy = potential[xindex][yindex][zindex-1].y;	//bulk
	  belowz = potential[xindex][yindex][zindex-1].z;	//bulk
	}
	if(zindex == zsteps - 1)      
	{
	  abovex = zero;					//nullo
	  abovey = zero;					//nullo
	  abovez = zero;					//nullo
	} else      
	{
	  abovex = potential[xindex][yindex][zindex+1].x;	//bulk
	  abovey = potential[xindex][yindex][zindex+1].y;	//bulk
	  abovez = potential[xindex][yindex][zindex+1].z;	//bulk
	}
	//X potenziale vettore
	currentj = currents[xindex][yindex][zindex].x;
	analyzedpoint = (abovex + belowx + leftx + rightx + beyondx + beforex - (mu_0*currentj))/6.;
	field[xindex][yindex][zindex].x = analyzedpoint;
	//Y potenziale vettore
	currentj = currents[xindex][yindex][zindex].y;
	analyzedpoint = (abovey + belowy + lefty + righty + beyondy + beforey - (mu_0*currentj))/6.;
	field[xindex][yindex][zindex].y = analyzedpoint;
	//Z potenziale vettore
	currentj = currents[xindex][yindex][zindex].z;
	analyzedpoint = (abovez + belowz + leftz + rightz + beyondz + beforez - (mu_0*currentj))/6.;
	field[xindex][yindex][zindex].z = analyzedpoint;
	
	zindex++;
      }
      yindex++;
    }
    xindex++;
  }
  
  xindex = x_border;
  yindex = y_border;
  zindex = z_border;

  while(xindex < xsteps)
  {
    yindex = y_border;
    while(yindex < ysteps)
    {
      zindex = z_border;
      while(zindex < zsteps)
      {
	if(xindex == x_border)                     
	{
	  beforex = field[xindex+1][yindex][zindex].x;	//riflessione
	  beforey = field[xindex+1][yindex][zindex].y;	//riflessione
	  beforez = field[xindex+1][yindex][zindex].z;	//riflessione
	} else
	{
	  beforex = field[xindex-1][yindex][zindex].x;	//bulk
	  beforey = field[xindex-1][yindex][zindex].y;	//bulk
	  beforez = field[xindex-1][yindex][zindex].z;	//bulk
	}
	if(xindex == xsteps - 1)            
	{
	  beyondx = zero;				//nullo
	  beyondy = zero;				//nullo
	  beyondz = zero;				//nullo
	} else 
	{
	  beyondx = field[xindex+1][yindex][zindex].x;	//bulk
	  beyondy = field[xindex+1][yindex][zindex].y;	//bulk
	  beyondz = field[xindex+1][yindex][zindex].z;	//bulk
	}
	if(yindex == y_border)   
	{
	  leftx = field[xindex][yindex+1][zindex].x;	//riflessione
	  lefty = field[xindex][yindex+1][zindex].y;	//riflessione
	  leftz = field[xindex][yindex+1][zindex].z;	//riflessione
	} else
	{
	  leftx = field[xindex][yindex-1][zindex].x;	//bulk
	  lefty = field[xindex][yindex-1][zindex].y;	//bulk
	  leftz = field[xindex][yindex-1][zindex].z;	//bulk
	}
	if(yindex == ysteps - 1)      
	{
	  rightx = zero;				//nullo
	  righty = zero;				//nullo
	  rightz = zero;				//nullo
	} else 
	{
	  rightx = field[xindex][yindex+1][zindex].x;	//bulk
	  righty = field[xindex][yindex+1][zindex].y;	//bulk
	  rightz = field[xindex][yindex+1][zindex].z;	//bulk
	}
	if(zindex == z_border)    
	{
	  belowx = field[xindex][yindex][zindex+1].x;	//riflessione
	  belowy = field[xindex][yindex][zindex+1].y;	//riflessione
	  belowz = field[xindex][yindex][zindex+1].z;	//riflessione
	} else 
	{
	  belowx = field[xindex][yindex][zindex-1].x;	//bulk
	  belowy = field[xindex][yindex][zindex-1].y;	//bulk
	  belowz = field[xindex][yindex][zindex-1].z;	//bulk
	}
	if(zindex == zsteps - 1)      
	{
	  abovex = zero;				//nullo
	  abovey = zero;				//nullo
	  abovez = zero;				//nullo
	} else      
	{
	  abovex = field[xindex][yindex][zindex+1].x;	//bulk
	  abovey = field[xindex][yindex][zindex+1].y;	//bulk
	  abovez = field[xindex][yindex][zindex+1].z;	//bulk
	}
	//X potenziale vettore
	currentj = currents[xindex][yindex][zindex].x;
	analyzedpoint = (abovex + belowx + leftx + rightx + beyondx + beforex - (mu_0*currentj))/6.;
	potential[xindex][yindex][zindex].x = analyzedpoint;
	//Y potenziale vettore
	currentj = currents[xindex][yindex][zindex].y;
	analyzedpoint = (abovey + belowy + lefty + righty + beyondy + beforey - (mu_0*currentj))/6.;
	potential[xindex][yindex][zindex].y = analyzedpoint;
	//Z potenziale vettore
	currentj = currents[xindex][yindex][zindex].z;
	analyzedpoint = (abovez + belowz + leftz + rightz + beyondz + beforez - (mu_0*currentj))/6.;
	potential[xindex][yindex][zindex].z = analyzedpoint;
	
	zindex++;
      }
      yindex++;
    }
    xindex++;
  }

}
//stima del campo magnetico a partire dal potenziale vettore appena calcolato
if(verboso) cout << "Inizio calcolo numerico del rotore..." << endl;
double dxy, dxz, dyx, dyz, dzx, dzy;



//////////////////////////////////////////////////////////////////////////////////
// solo per visualizzare il campo su root
TH2F * h1 = new TH2F("h1", "MagneticField x", xsteps+2, x_border-1, xsteps+1, ysteps+2, y_border-1, ysteps+1);
TH2F * h2 = new TH2F("h2", "MagneticField y", xsteps+2, x_border-1, xsteps+1, ysteps+2, y_border-1, ysteps+1);
TH2F * h3 = new TH2F("h3", "MagneticField z", xsteps+2, x_border-1, xsteps+1, ysteps+2, y_border-1, ysteps+1);
///////////////////////////////////////////////////////////////////////////////////

double actualx, actualy, actualz;

xindex = x_border;
yindex = y_border;
zindex = z_border;

while(xindex < xsteps)
{
  yindex = y_border;
  while(yindex < ysteps)
  {
    zindex = z_border;
    while(zindex < zsteps)
    {
      if(xindex == x_border)                     
      {
	beforex = potential[xindex+1][yindex][zindex].x;	//riflessione
	beforey = potential[xindex+1][yindex][zindex].y;	//riflessione
	beforez = potential[xindex+1][yindex][zindex].z;	//riflessione
      } else
      {
	beforex = potential[xindex-1][yindex][zindex].x;	//bulk
	beforey = potential[xindex-1][yindex][zindex].y;	//bulk
	beforez = potential[xindex-1][yindex][zindex].z;	//bulk
      }
      if(xindex == xsteps - 1)            
      {
	beyondx = zero;						//nullo
	beyondy = zero;						//nullo
	beyondz = zero;						//nullo
      } else 
      {
	beyondx = potential[xindex+1][yindex][zindex].x;	//bulk
	beyondy = potential[xindex+1][yindex][zindex].y;	//bulk
	beyondz = potential[xindex+1][yindex][zindex].z;	//bulk
      }
      if(yindex == y_border)   
      {
	leftx = potential[xindex][yindex+1][zindex].x;		//riflessione
	lefty = potential[xindex][yindex+1][zindex].y;		//riflessione
	leftz = potential[xindex][yindex+1][zindex].z;		//riflessione
      } else
      {
	leftx = potential[xindex][yindex-1][zindex].x;		//bulk
	lefty = potential[xindex][yindex-1][zindex].y;		//bulk
	leftz = potential[xindex][yindex-1][zindex].z;		//bulk
      }
      if(yindex == ysteps - 1)      
      {
	rightx = zero;						//nullo
	righty = zero;						//nullo
	rightz = zero;						//nullo
      } else 
      {
	rightx = potential[xindex][yindex+1][zindex].x;		//bulk
	righty = potential[xindex][yindex+1][zindex].y;		//bulk
	rightz = potential[xindex][yindex+1][zindex].z;		//bulk
      }
      if(zindex == z_border)    
      {
	belowx = potential[xindex][yindex][zindex+1].x;		//riflessione
	belowy = potential[xindex][yindex][zindex+1].y;		//riflessione
	belowz = potential[xindex][yindex][zindex+1].z;		//riflessione
      } else 
      {
	belowx = potential[xindex][yindex][zindex-1].x;		//bulk
	belowy = potential[xindex][yindex][zindex-1].y;		//bulk
	belowz = potential[xindex][yindex][zindex-1].z;		//bulk
      }
      if(zindex == zsteps - 1)      
      {
	abovex = zero;						//nullo
	abovey = zero;						//nullo
	abovez = zero;						//nullo
      } else      
      {
	abovex = potential[xindex][yindex][zindex+1].x;		//bulk
	abovey = potential[xindex][yindex][zindex+1].y;		//bulk
	abovez = potential[xindex][yindex][zindex+1].z;		//bulk
      }
      
      
      dxy = (rightx - leftx)    / (2 * delta);			//calcolo delle derivate
      dxz = (abovex - belowx)   / (2 * delta);
      dyx = (beyondy - beforey) / (2 * delta);
      dyz = (abovey - belowy)   / (2 * delta);
      dzx = (beyondz - beforez) / (2 * delta);
      dzy = (rightz - leftz)    / (2 * delta);
 
      if(zindex == 0)						//sistema i pèroblemi relativi al calcolo del rotore in presenza di condizioni al contorno riflettenti
      {
	dxz = (abovex - potential[xindex][yindex][zindex].x) /delta;
	dyz = (abovey - potential[xindex][yindex][zindex].y) /delta;
      }
      
      if(yindex == 0)
      {
	dxy = (rightx - potential[xindex][yindex][zindex].x) /delta;
	dzy = (rightz - potential[xindex][yindex][zindex].z) /delta;
      }
      
      if(xindex == 0)
      {
	dyx = (beyondy - potential[xindex][yindex][zindex].y) /delta;
	dzx = (beyondz - potential[xindex][yindex][zindex].z) /delta;
      }
      
      
      actualx = dzy - dyz;					//calcolo del rotore
      actualy = dxz - dzx;
      actualz = dyx - dxy;
      
      
      
      field[xindex][yindex][zindex].x = actualx;		//riempimento effettivo campo
      field[xindex][yindex][zindex].y = actualy;
      field[xindex][yindex][zindex].z = actualz;
      
      //check per il debug
      if(xindex == 0 && zindex == 0 && yindex == 0){
      cout << "Centro: " << actualx << " = " << dzy - dyz << endl;}
      
      if(xindex == 1 && zindex == 0 && yindex == 0){
      cout << "Primo punto: " << actualx << " = " << dzy - dyz << endl;} 
      
      if(potential[xindex][yindex][zindex].x != 0) cout << "POTENZIALE LUNGO X!!" << endl;
      if(verboso && xindex == 0)
      {
// 	cout << actualx << " ";
	h1->SetBinContent(yindex+2, zindex+2, field[xindex][yindex][zindex].x);
// 	h2->SetBinContent(xindex+2, yindex+2, actualx);
// 	h3->SetBinContent(xindex+2, yindex+2, actualx);
      }
      zindex++;
    }
    
//     if(verboso) cout << endl;
    yindex++;
  }
//   if(verboso) cout << endl;
  xindex++;
}



  h1->Draw("colz");
//   h2->Draw("colz");
//   h3->Draw("colz");

  return 0;
//   return B{1,0,0};
};


void PositionToGrid(const Position3D * position, int & xindex, int & yindex, int & zindex)
{
  double delta = 0.1;		//dovrebbe vederla dal programma andrà implementato nella classe
  xindex = int(position->x/delta);
  yindex = int(position->y/delta);
  zindex = int(position->z/delta);
  return;
}

void GridToPosition(Position3D * position, const int & xindex, const int & yindex, const int & zindex)
{
  double delta = 0.1;
  position->x = xindex * delta + (0.5 * delta);
  position->y = yindex * delta + (0.5 * delta);
  position->z = zindex * delta + (0.5 * delta);
  return;
}
