#include "MagneticField.h"

int main()
{
  MagneticField * obj = new MagneticField(2000., 2000., 200., 1000., 55., 55., 2.5, 10., 2000);
//   MagneticField * obj = new MagneticField(2000., 2000., 2000., 1000., 550., 20., 20, 10., 1000);			//fisico ma cubico
//   MagneticField * obj = new MagneticField(10., 10., 10., 7., 1., 1., 0.01, 10., 1000);
//   MagneticField * obj = new MagneticField(1000., 1000., 1000., 50., 50., 50., 10, 10., 1000);
//   MagneticField * obj = new MagneticField("field_output");
  obj->CreateRootOutput("field_output");
//   obj->CreateTextOutput("field_output");
  
  return 0;
}
