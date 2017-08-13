#include <iostream>
#include <cmath>
#include <sstream>
#include <time.h>

#include "MagneticField.h"

using namespace std;

#include "TFile.h"
#include "TH2F.h"


MagneticField::MagneticField(double xsyst, double ysyst, double zsyst, double xsol, double ysol, double zsol, double delta, double current, long ninteractions):
  xsyst(xsyst),
  ysyst(ysyst),
  zsyst(zsyst),
  xsol(xsol),
  ysol(ysol),
  zsol(zsol),
  delta(delta),
  current(current),
  ninteractions(ninteractions),
  mu_0(1),
  verboso(true),
  co(true),
  zero(0.),
  pot_computed(false),
  field_computed(false)
  {
    xsteps = int(xsyst/(delta*2)) + 1; 
    ysteps = int(ysyst/(delta*2)) + 1; 
    zsteps = int(zsyst/(delta*2)) + 1;
    
    if(verboso) cout <<"xsteps = " << xsteps << " ysteps = " << ysteps << " zsteps = " << zsteps << endl;
    
    
    //creazione delle griglie di lavoro tridimentionali e loro inizializzazione: RIDUZIONE SISTEMA A 1/8
    cout << "Inizio creazione griglie..." << endl;
    fieldLayer2.resize(xsteps * ysteps * zsteps, {0.,0.,0.});		//contiene tutto il campo
    fieldLayer1.resize(xsteps * ysteps);				//contiene puntatori a un piano orizzontale
    field.resize(xsteps);						//contiene puntatori alla posizione lungo l'asse del solenoide
    
    potentialLayer2.resize(xsteps * ysteps * zsteps, {0.,0.,0.});	//contiene tutto il campo
    potentialLayer1.resize(xsteps * ysteps);				//contiene puntatori a un piano orizzontale
    potential.resize(xsteps);						//contiene puntatori alla posizione lungo l'asse del solenoide
    
    /*currentsLayer2.resize(xsteps * ysteps * zsteps, {0.,0.,0.});	//contiene tutto il campo
    currentsLayer1.resize(xsteps * ysteps);				//contiene puntatori a un piano orizzontale
    currents.resize(xsteps);						//contiene puntatori alla posizione lungo l'asse del solenoide */
    
    for(unsigned int i = 0; i < fieldLayer1.size(); i++)		//sistema i puntatori in modo che puntino al posto giusto
    {
      fieldLayer1.at(i) = & fieldLayer2.at(i * ysteps);
      potentialLayer1.at(i) = & potentialLayer2.at(i * ysteps);
      //currentsLayer1.at(i) = & currentsLayer2.at(i * ysteps);
    }
    for(unsigned int i = 0; i < field.size(); i++)
    {
      field.at(i) = & fieldLayer1.at(i * xsteps);
      potential.at(i) = & potentialLayer1.at(i * xsteps);
      //currents.at(i) = & currentsLayer1.at(i * xsteps);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    double perimeter = (2 * ysol) + (2 * zsol); 
    current_density = current/perimeter;		//sarà la densità di corrente
    
    Position3D current_corner = {xsol*.5, ysol*.5, zsol*.5};
    PositionToGrid(&current_corner, xcurrent, ycurrent, zcurrent);
    
    PotentialInitialization();
  }
  
  MagneticField::~MagneticField()
  { }
  
  
  void MagneticField::PotentialInitialization()
  {
    //inizializzazione del potenziale vettore
    if (verboso) cout << "Inizio inizializzazione potenziale vettore..." << endl;
    starting_point  = {0.,0.,0.};		//inizia dal centro del sistema
    walking_point   = starting_point;
    Vector3D * modify_inloop;		//variabile temporaanea
    double rho = 0;
    double theta = 0;
    
    debug1 = 0;
    debug2 = 0;
    debug3 = 0;
    
    //   PositionToGrid(& walking_point, x_border, y_border, z_border);
    
    //set del potenziale vettore all'interno del solenoide
    double mu_0 = 1.;
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
	  rho = mu_0 * current * .5 * sqrt(walking_point.y * walking_point.y + walking_point.z * walking_point.z);
	  if(walking_point.y == 0) theta = M_PI/2;				//serve a evitare infiniti
	  else theta = atan(walking_point.z/walking_point.y);		//ordine giusto? credo di sì, basta mantenere quello sopra qua come quello con seno quando riempio
	  * modify_inloop = {0., -rho * sin(theta), rho * cos(theta)};
	  yindex++;
	  GridToPosition(& walking_point, xindex, yindex, zindex);
	  debug1++;
	  if( debug1 > 100000 && verboso)
	  {
	    cout << "ERROR IN CYCLE4" << endl;
	    return;
	  }
	}
	zindex++;					//avanza lungo l'altezza
	GridToPosition(& walking_point, xindex, yindex, zindex);
	debug2++;
	if( debug2 > 100000 && verboso)
	{
	  cout << "ERROR IN CYCLE5" << endl;
	  return;
	}
      }
      xindex++;					//avanza lungo l'asse
      GridToPosition(& walking_point, xindex, yindex, zindex);
      debug1++;
      if( debug3 > 100000 && verboso)
      {
	cout << "ERROR IN CYCLE6" << endl;
	return;
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
	  rho = mu_0 * current * .5 * zsol * .5 * zsol * .5 / sqrt(walking_point.y * walking_point.y + walking_point.z * walking_point.z);
	  if(walking_point.y == 0) theta = M_PI/2;				//serve a evitare infiniti
	  else theta = atan(walking_point.z/walking_point.y);		//ordine giusto? credo di sì, basta mantenere quello sopra qua come quello con seno quando riempio
	  * modify_inloop = {0., -rho * sin(theta), rho * cos(theta)};	  
	  yindex++;
	  GridToPosition(& walking_point, xindex, yindex, zindex);
	  debug1++;
	  if( debug1 > 100000 && verboso)
	  {
	    cout << "ERROR IN CYCLE7" << endl;
	    return;
	  }
	}
	zindex++;					//avanza lungo l'altezza
	GridToPosition(& walking_point, xindex, yindex, zindex);
	debug2++;
	if( debug2 > 100000 && verboso)
	{
	  cout << "ERROR IN CYCLE8" << endl;
	  return;
	}
      }
      xindex++;					//avanza lungo l'asse
      GridToPosition(& walking_point, xindex, yindex, zindex);
      debug3++;
      if( debug3 > 100000 && verboso)
      {
	cout << "ERROR IN CYCLE9" << endl;
	return;
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
	  rho = mu_0 * current * .5 * ysol * .5 * ysol * .5 / sqrt(walking_point.y * walking_point.y + walking_point.z * walking_point.z);
	  if(walking_point.y == 0) theta = M_PI/2;				//serve a evitare infiniti
	  else theta = atan(walking_point.z/walking_point.y);		//ordine giusto? credo di sì, basta mantenere quello sopra qua come quello con seno quando riempio
	  * modify_inloop = {0., -rho * sin(theta), rho * cos(theta)};
	  if(verboso && modify_inloop->z > 10000 ) cout << modify_inloop->z << " at " << xindex << " " << yindex << " " << zindex << " ; " << rho << " " << theta << endl;
	  yindex++;
	  GridToPosition(& walking_point, xindex, yindex, zindex);
	  debug1++;
	  if( debug1 > 100000 && verboso)
	  {
	    cout << "ERROR IN CYCLE10"  << endl;
	    cout << "yindex = " << yindex << " ysol*.5 = " << ysol *.5 << " ysyst = " << ysyst << endl;
	    return;
	  }
	}
	zindex++;					//avanza lungo l'altezza
	GridToPosition(& walking_point, xindex, yindex, zindex);
	debug2++;
	if( debug2 > 100000 && verboso)
	{
	  cout << "ERROR IN CYCLE11" << endl;
	  return;
	}
      }
      xindex++;					//avanza lungo l'asse
      GridToPosition(& walking_point, xindex, yindex, zindex);
      debug3++;
      if( debug3 > 100000 && verboso)
      {
	cout << "ERROR IN CYCLE12" << endl;
	return;
      }
    }
    
    debug1 = 0;
    debug2 = 0;
    debug3 = 0;
    
    //set del potenziale vettore all'esterno del solenoide !!TRA ALTO E DI LATO, IL QUADRATONE!!
    
    
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
	  rho = mu_0 * current * .5 * zsol * .5 * ysol * .5 / sqrt(walking_point.y * walking_point.y + walking_point.z * walking_point.z);	//ibrido tra il sopra e il di lato
	  if(walking_point.y == 0) theta = M_PI/2;				//serve a evitare infiniti
	  else theta = atan(walking_point.z/walking_point.y);		//ordine giusto? credo di sì, basta mantenere quello sopra qua come quello con seno quando riempio
	  * modify_inloop = {0., -rho * sin(theta), rho * cos(theta)};
	  if(verboso && modify_inloop->z > 10000 ) cout << modify_inloop->z << " at " << xindex << " " << yindex << " " << zindex << " ; " << rho << " " << theta << endl;
	  yindex++;
	  GridToPosition(& walking_point, xindex, yindex, zindex);
	  debug1++;
	  if( debug1 > 100000 && verboso)
	  {
	    cout << "ERROR IN CYCLE13" << endl;
	    return;
	  }
	}
	zindex++;					//avanza lungo l'altezza
	GridToPosition(& walking_point, xindex, yindex, zindex);
	debug2++;
	if( debug2 > 100000 && verboso)
	{
	  cout << "ERROR IN CYCLE14" << endl;
	  return;
	}
      }
      xindex++;					//avanza lungo l'asse
      GridToPosition(& walking_point, xindex, yindex, zindex);
      debug3++;
      if( debug3 > 100000 && verboso)
      {
	cout << "ERROR IN CYCLE15" << endl;
	return;
      }
    } 
    //il resto del potenziale vettore è alsciato a zero, verrà poi modificato dall'algoritmo numerico.
  }
  
  
  void MagneticField::ComputePotential()
  {
    //calcolo del potenziale vettore utilizzando l'algoritmo di Jacobi, usando i due campi field e potential
    //le condizioni al contorno sono miste: quando il punto è alla fine dell'intero sistema allora il potenziale vettore è fissato a zero,
    //quando invece si è nelle "pareti fittizie" (quelle date dalla divisione per simmetria del problema) allora sono riflettenti.

    xborder = 0;
    yborder = 0;
    zborder = 0;
    
    time_t beginning = time(NULL);
    
    Vector3D currentj;
    
    if(verboso) cout << "Inizio computazione del potenziale vettore..." << endl;
    for(unsigned int i = 0; i < ninteractions; i++)
    {
      
      xindex = xborder;
      yindex = yborder;
      zindex = zborder;
      
      while(xindex < xsteps)
      {
	yindex = yborder;
	while(yindex < ysteps)
	{
	  zindex = zborder;
	  while(zindex < zsteps)
	  {
	    if(xindex == xborder)                     
	    {
	      beforey = potential[xindex+1][yindex][zindex].y;	//riflessione
	      beforez = potential[xindex+1][yindex][zindex].z;	//riflessione
	    } else
	    {
	      beforey = potential[xindex-1][yindex][zindex].y;	//bulk
	      beforez = potential[xindex-1][yindex][zindex].z;	//bulk
	    }
	    if(xindex == xsteps - 1)            
	    {
	      beyondy = zero;					//nullo
	      beyondz = zero;					//nullo
	    } else 
	    {
	      beyondy = potential[xindex+1][yindex][zindex].y;	//bulk
	      beyondz = potential[xindex+1][yindex][zindex].z;	//bulk
	    }
	    if(yindex == yborder)   
	    {
	      lefty = potential[xindex][yindex+1][zindex].y;	//riflessione
	      leftz = potential[xindex][yindex+1][zindex].z;	//riflessione
	    } else
	    {
	      lefty = potential[xindex][yindex-1][zindex].y;	//bulk
	      leftz = potential[xindex][yindex-1][zindex].z;	//bulk
	    }
	    if(yindex == ysteps - 1)      
	    {
	      righty = zero;					//nullo
	      rightz = zero;					//nullo
	    } else 
	    {
	      righty = potential[xindex][yindex+1][zindex].y;	//bulk
	      rightz = potential[xindex][yindex+1][zindex].z;	//bulk
	    }
	    if(zindex == zborder)    
	    {
	      belowy = potential[xindex][yindex][zindex+1].y;	//riflessione
	      belowz = potential[xindex][yindex][zindex+1].z;	//riflessione
	    } else 
	    {
	      belowy = potential[xindex][yindex][zindex-1].y;	//bulk
	      belowz = potential[xindex][yindex][zindex-1].z;	//bulk
	    }
	    if(zindex == zsteps - 1)      
	    {
	      abovey = zero;					//nullo
	      abovez = zero;					//nullo
	    } else      
	    {
	      abovey = potential[xindex][yindex][zindex+1].y;	//bulk
	      abovez = potential[xindex][yindex][zindex+1].z;	//bulk
	    }
	    if(co == false) {//X potenziale vettore
	      if(xindex == xborder)    beforex = potential[xindex+1][yindex][zindex].x;		//riflessione
	      else                      beforex = potential[xindex-1][yindex][zindex].x;	//bulk
	      if(xindex == xsteps - 1)  beyondx = zero;						//nullo
	      else                      beyondx = potential[xindex+1][yindex][zindex].x;	//bulk
	      if(yindex == yborder)    leftx = potential[xindex][yindex+1][zindex].x;		//riflessione
	      else                      leftx = potential[xindex][yindex-1][zindex].x;		//bulk
	      if(yindex == ysteps - 1)  rightx = zero;						//nullo
	      else                      rightx = potential[xindex][yindex+1][zindex].x;		//bulk
	      if(zindex == zborder)    belowx = potential[xindex][yindex][zindex+1].x;		//riflessione
	      else                      belowx = potential[xindex][yindex][zindex-1].x;		//bulk
	      if(zindex == zsteps - 1)  abovex = zero;						//nullo
	      else	                    abovex = potential[xindex][yindex][zindex+1].x;	//bulk
	      currentj = CurrentAt(xindex, yindex, zindex);
	      analyzedpoint = (abovex + belowx + leftx + rightx + beyondx + beforex - (mu_0*currentj.x))/6.;
	      field[xindex][yindex][zindex].x = analyzedpoint;  
	    }
	    //Y potenziale vettore
	    currentj = CurrentAt(xindex, yindex, zindex);
	    analyzedpoint = (abovey + belowy + lefty + righty + beyondy + beforey - (mu_0*currentj.y))/6.;
	    field[xindex][yindex][zindex].y = analyzedpoint;
	    //Z potenziale vettore
	    //currentj = currents[xindex][yindex][zindex].z;
	    analyzedpoint = (abovez + belowz + leftz + rightz + beyondz + beforez - (mu_0*currentj.z))/6.;
	    field[xindex][yindex][zindex].z = analyzedpoint;
	    
	    zindex++;
	  }
	  yindex++;
	}
	xindex++;
      }
      
      xindex = xborder;
      yindex = yborder;
      zindex = zborder;
      
      while(xindex < xsteps)
      {
	yindex = yborder;
	while(yindex < ysteps)
	{
	  zindex = zborder;
	  while(zindex < zsteps)
	  {
	    if(xindex == xborder)                     
	    {
	      beforey = field[xindex+1][yindex][zindex].y;	//riflessione
	      beforez = field[xindex+1][yindex][zindex].z;	//riflessione
	    } else
	    {
	      beforey = field[xindex-1][yindex][zindex].y;	//bulk
	      beforez = field[xindex-1][yindex][zindex].z;	//bulk
	    }
	    if(xindex == xsteps - 1)            
	    {
	      beyondy = zero;					//nullo
	      beyondz = zero;					//nullo
	    } else 
	    {
	      beyondy = field[xindex+1][yindex][zindex].y;	//bulk
	      beyondz = field[xindex+1][yindex][zindex].z;	//bulk
	    }
	    if(yindex == yborder)   
	    {
	      lefty = field[xindex][yindex+1][zindex].y;	//riflessione
	      leftz = field[xindex][yindex+1][zindex].z;	//riflessione
	    } else
	    {
	      lefty = field[xindex][yindex-1][zindex].y;	//bulk
	      leftz = field[xindex][yindex-1][zindex].z;	//bulk
	    }
	    if(yindex == ysteps - 1)      
	    {
	      righty = zero;					//nullo
	      rightz = zero;					//nullo
	    } else 
	    {
	      righty = field[xindex][yindex+1][zindex].y;	//bulk
	      rightz = field[xindex][yindex+1][zindex].z;	//bulk
	    }
	    if(zindex == zborder)    
	    {
	      belowy = field[xindex][yindex][zindex+1].y;	//riflessione
	      belowz = field[xindex][yindex][zindex+1].z;	//riflessione
	    } else 
	    {
	      belowy = field[xindex][yindex][zindex-1].y;	//bulk
	      belowz = field[xindex][yindex][zindex-1].z;	//bulk
	    }
	    if(zindex == zsteps - 1)      
	    {
	      abovey = zero;					//nullo
	      abovez = zero;					//nullo
	    } else      
	    {
	      abovey = field[xindex][yindex][zindex+1].y;	//bulk
	      abovez = field[xindex][yindex][zindex+1].z;	//bulk
	    }
	    if(co == false) {//X potenziale vettore
	      if(xindex == xborder)    beforex = field[xindex+1][yindex][zindex].x;	//riflessione
	      else                      beforex = field[xindex-1][yindex][zindex].x;	//bulk
	      if(xindex == xsteps - 1)  beyondx = zero;					//nullo
	      else                      beyondx = field[xindex+1][yindex][zindex].x;	//bulk
	      if(yindex == yborder)    leftx = field[xindex][yindex+1][zindex].x;	//riflessione
	      else                      leftx = field[xindex][yindex-1][zindex].x;	//bulk
	      if(yindex == ysteps - 1)  rightx = zero;					//nullo
	      else                      rightx = field[xindex][yindex+1][zindex].x;	//bulk
	      if(zindex == zborder)    belowx = field[xindex][yindex][zindex+1].x;	//riflessione
	      else                      belowx = field[xindex][yindex][zindex-1].x;	//bulk
	      if(zindex == zsteps - 1)  abovex = zero;					//nullo
	      else	                    abovex = field[xindex][yindex][zindex+1].x;	//bulk
	      currentj = CurrentAt(xindex, yindex, zindex);
	      analyzedpoint = (abovex + belowx + leftx + rightx + beyondx + beforex - (mu_0*currentj.x))/6.;
	      potential[xindex][yindex][zindex].x = analyzedpoint;  
	    }
	    //Y potenziale vettore
	    currentj = CurrentAt(xindex, yindex, zindex);
	    analyzedpoint = (abovey + belowy + lefty + righty + beyondy + beforey - (mu_0*currentj.y))/6.;
	    potential[xindex][yindex][zindex].y = analyzedpoint;
	    //Z potenziale vettore
	    //currentj = CurrentAt(xindex, yindex, zindex);
	    analyzedpoint = (abovez + belowz + leftz + rightz + beyondz + beforez - (mu_0*currentj.z))/6.;
	    potential[xindex][yindex][zindex].z = analyzedpoint;
	    
	    zindex++;
	  }
	  yindex++;
	}
	xindex++;
      }
      
    }
    
    time_t end = time(NULL);
    if(verboso) cout << "Tempo di calcolo = " << difftime(end, beginning) << " s " << endl;
    
    pot_computed = true;
    return;
  }
  
  void MagneticField::ComputeField()
  {
    //stima del campo magnetico a partire dal potenziale vettore appena calcolato
    
    if(pot_computed == false) ComputePotential();
    
    if(verboso) cout << "Inizio calcolo numerico del rotore..." << endl;
    
    double actualx, actualy, actualz;
    
    xborder = 0;
    yborder = 0;
    zborder = 0;
    
    xindex = xborder;
    yindex = yborder;
    zindex = zborder;
    
    while(xindex < xsteps)
    {
      yindex = yborder;
      while(yindex < ysteps)
      {
	zindex = zborder;
	while(zindex < zsteps)
	{
	  if(xindex == xborder)                     
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
	  if(yindex == yborder)   
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
	  if(zindex == zborder)    
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
	  
	  if(zindex == 0)						//sistema i problemi relativi al calcolo del rotore in presenza di condizioni al contorno riflettenti
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
	  
	  zindex++;
	}
	yindex++;
      }
      xindex++;
    }
    
    field_computed = true;
    return;
    
  }
  
  
  Vector3D MagneticField::CurrentAt(Position3D pos)
  {
    unsigned int xtemp, ytemp, ztemp;
    PositionToGrid(& pos, xtemp, ytemp, ztemp);
    return CurrentAt(xtemp, ytemp, ztemp);
  }
  
  
  Vector3D MagneticField::CurrentAt(int x, int y, int z)
  {
    if( x < xcurrent && y == ycurrent && z == zcurrent) return {0., -current_density * .5, current_density * .5};
    if( x < xcurrent && y == ycurrent &&  z < zcurrent) return {0., 0., current_density};
    if( x < xcurrent && y < ycurrent  && z == zcurrent) return {0., -current_density, 0.};
    return {0., 0., 0.,};
  }
  
  
  
  void MagneticField::PositionToGrid(const Position3D * position, unsigned int & xindex, unsigned int & yindex, unsigned int & zindex)
  {
    double delta = 0.1;		//dovrebbe vederla dal programma andrà implementato nella classe
    xindex = int(position->x/delta);
    yindex = int(position->y/delta);
    zindex = int(position->z/delta);
    return;
  }
  
  void MagneticField::GridToPosition(Position3D * position, const unsigned int & xindex, const unsigned int & yindex, const unsigned int & zindex)
  {
    double delta = 0.1;
    position->x = xindex * delta + (0.5 * delta);
    position->y = yindex * delta + (0.5 * delta);
    position->z = zindex * delta + (0.5 * delta);
    return;
  }

  
  
  void MagneticField::CreateRootOutput(string outname)
  {
    
    if (field_computed == false) ComputeField();
    
    
    xborder = 0;
    yborder = 0;
    zborder = 0;
    
    TDirectory* currentDir = gDirectory;
    stringstream sstr; 
    sstr << outname << ".root";
    string str = sstr.str();
    TFile* outfile = new TFile( str.c_str(), "RECREATE" );
    if ( outfile->IsOpen() ) cout << "ROOT File opened successfully" << endl;
    
    TH2F * hxx = new TH2F("hxx", "MagneticField x in x = 0", ysteps+1, yborder, ysteps+1, zsteps+1, zborder, zsteps+1);
    TH2F * hyx = new TH2F("hyx", "MagneticField y in x = 0", ysteps+1, yborder, ysteps+1, zsteps+1, zborder, zsteps+1);
    TH2F * hzx = new TH2F("hzx", "MagneticField z in x = 0", ysteps+1, yborder, ysteps+1, zsteps+1, zborder, zsteps+1);
    
    TH2F * hxy = new TH2F("hxy", "MagneticField x in y = 0", xsteps+1, xborder, xsteps+1, zsteps+1, zborder, zsteps+1);
    TH2F * hyy = new TH2F("hyy", "MagneticField y in y = 0", xsteps+1, xborder, xsteps+1, zsteps+1, zborder, zsteps+1);
    TH2F * hzy = new TH2F("hzy", "MagneticField z in y = 0", xsteps+1, xborder, xsteps+1, zsteps+1, zborder, zsteps+1);
    
    TH2F * hxz = new TH2F("hxz", "MagneticField x in z = 0", xsteps+1, xborder, xsteps+1, ysteps+1, yborder, ysteps+1);
    TH2F * hyz = new TH2F("hyz", "MagneticField y in z = 0", xsteps+1, xborder, xsteps+1, ysteps+1, yborder, ysteps+1);
    TH2F * hzz = new TH2F("hzz", "MagneticField z in z = 0", xsteps+1, xborder, xsteps+1, ysteps+1, yborder, ysteps+1);
    
    
    xindex = xborder;
    yindex = yborder;
    zindex = zborder;
    
    while(xindex < xsteps)
    {
      yindex = yborder;
      while(yindex < ysteps)
      {
	zindex = zborder;
	while(zindex < zsteps)
	{
	  if(xindex == 0)  hxx->SetBinContent(yindex+1, zindex+1, field[xindex][yindex][zindex].x);
	  if(xindex == 0)  hyx->SetBinContent(yindex+1, zindex+1, field[xindex][yindex][zindex].y);
	  if(xindex == 0)  hzx->SetBinContent(yindex+1, zindex+1, field[xindex][yindex][zindex].z);
	  
	  if(yindex == 0)  hxy->SetBinContent(xindex+1, zindex+1, field[xindex][yindex][zindex].x);
	  if(yindex == 0)  hyy->SetBinContent(xindex+1, zindex+1, field[xindex][yindex][zindex].y);
	  if(yindex == 0)  hzy->SetBinContent(xindex+1, zindex+1, field[xindex][yindex][zindex].z);
	  
	  if(zindex == 0)  hxz->SetBinContent(xindex+1, yindex+1, field[xindex][yindex][zindex].x);
	  if(zindex == 0)  hyz->SetBinContent(xindex+1, yindex+1, field[xindex][yindex][zindex].y);
	  if(zindex == 0)  hzz->SetBinContent(xindex+1, yindex+1, field[xindex][yindex][zindex].z);
	  
	  zindex++;
	}
	yindex++;
      }
      xindex++;
    }
    
    hxx->Write();
    hyx->Write();
    hzx->Write();
    
    hxy->Write();
    hyy->Write();
    hzy->Write();
    
    hxz->Write();
    hyz->Write();
    hzz->Write();
    
    outfile->Close();
    delete outfile;
    currentDir->cd();
    
    return;
  }