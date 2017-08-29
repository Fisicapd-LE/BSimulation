#include <string>
#include <vector>
#include "../g-2Simulation/src/Utilities/Vector3D.h"

class MagneticField
{
  
  
public:
  
  MagneticField(double xsyst, double ysyst, double zsyst, double xsol, double ysol, double zsol, double delta, double current, double convergence_rate, double max_time);
  MagneticField(std::string in_name);
  virtual ~MagneticField();
  
  void CreateRootOutput(std::string out_name);
  void CreateTextOutput(std::string out_name);
  
  Vector3D GetB(Position3D pos);	//la funzione per cui tutto il programma è stato scritto
  
private: 
  
  double xsyst;			//system dimensions
  double ysyst;
  double zsyst;
  double xsol;
  double ysol;
  double zsol;
  
  double delta;
  double current;
  
  unsigned int xsteps;
  unsigned int ysteps;
  unsigned int zsteps;
  
  double mu_0;
  double current_density;
  
  bool verboso;				//necessario per il debugging, fa capire a che punto si è e evita di bloccarsi nei loop
  bool co;				//computational optimization: se questo parametro è settato su true, allora si ignorerà in pratica la componente x
  //dei campi vettoriali, perché non necessaria. Ciò è possibile se non si hanno correnti lungo x ed è un asse di simmetri
  //per il sistema prima della riduzione delle dimensioni
  
  
  //Campo magnetico
  std::vector<Vector3D>   fieldLayer2;		//identifica la z una volta note x e y
  std::vector<Vector3D*>  fieldLayer1;		//identifica la y una volta nota la x
  std::vector<Vector3D**> field;			//identificaa la x (asse)
  
  //Potenziale vettore
  std::vector<Vector3D>   potentialLayer2;
  std::vector<Vector3D*>  potentialLayer1;
  std::vector<Vector3D**> potential;
  
  //Correnti
  //vector<Vector3D>   currentsLayer2;
  //vector<Vector3D*>  currentsLayer1;
  //vector<Vector3D**> currents;
  
  //variabili per i cicli
  unsigned int xindex;
  unsigned int yindex;
  unsigned int zindex;
  
  unsigned int xborder;
  unsigned int yborder;
  unsigned int zborder;
  
  Position3D starting_point;
  Position3D walking_point;
  
  //funzioni utili per passaggio da continuo a discreto e viceversa
  void PositionToGrid(const Position3D * position, unsigned int & xindex, unsigned int & yindex, unsigned int & zindex);
  void GridToPosition(Position3D * position, const unsigned int & xindex, const unsigned int & yindex, const unsigned int & zindex);
  
  //variabili per trovare la corrente
  unsigned int xcurrent;
  unsigned int ycurrent;
  unsigned int zcurrent;
  
  //funzioni per avere la corrente
  Vector3D CurrentAt(Position3D pos);
  Vector3D CurrentAt(int x, int y, int z);
  
  //variabili di sicurezza
  unsigned int debug1;
  unsigned int debug2;
  unsigned int debug3;
  
  //funzione di inizializzazione del potenziale vettore
  void PotentialInitialization();
  
  //variabili per l'algoritmo numerico
//   long ninteractions;
  double convergence_rate;
  double max_time;
  
  bool pot_computed;		//se falsa, non è stato calcolato il campo magnetico
  
  double zero;
  
  double analyzedpoint;
  double beforex, beforey, beforez;
  double beyondx, beyondy, beyondz;
  double leftx, lefty, leftz;
  double rightx, righty, rightz;
  double belowx, belowy, belowz;
  double abovex, abovey, abovez;
  
  void ComputePotential();		//stima il potenziale utilizzando l'algoritmo di Jacobi
  
  //variabili per il calcolo del campo magnetico
  bool field_computed;
  double dxy, dxz, dyx, dyz, dzx, dzy;
  
  void ComputeField();
  
};
