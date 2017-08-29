#include "MagneticField.h"

int main()
{
//    MagneticField * obj = new MagneticField(1750., 2200., 80., 1000., 550., 20., 1./1.6, 10., 5 * 1000);	//12399 s
//    MagneticField * obj = new MagneticField(2000., 2000., 80., 1000., 550., 20., 2., 10.,10 * 1000);		//provo ad allungare x, forse quello dà problemi
														//874 s
//    MagneticField * obj = new MagneticField(2000., 2000., 80., 1000., 550., 20., 1., 10., 5 * 1000);		//3171 s PROBLEMATICO
//    MagneticField * obj = new MagneticField(2000., 2000., 500., 1000., 550., 117., 2.5, 10., 0.0, 480);	//08_28
//   MagneticField * obj = new MagneticField(2000., 60., 60., 1000., 10., 20., 1, 10., 8 * 2000);		//funziona abbastanza bene
//   MagneticField * obj = new MagneticField(2000., 2000., 2000., 1000., 550., 20., 20, 10., 1000);			//fisico ma cubico
//   MagneticField * obj = new MagneticField(10., 10., 10., 7., 1., 1., 0.1, 10., 0.1, 10);
//   MagneticField * obj = new MagneticField(1000., 1000., 1000., 50., 50., 50., 10, 10., 1000);
  MagneticField * obj = new MagneticField("field_output");
  obj->CreateRootOutput("field_output");
//   obj->CreateTextOutput("field_output");
  
  return 0;
}
