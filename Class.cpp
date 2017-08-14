#include "MagneticField.h"

int main()
{
  MagneticField * obj = new MagneticField(2000., 1100., 50., 1000., 550., 25., .5, 10., 1000);
//   MagneticField * obj = new MagneticField("field_output");
  obj->CreateRootOutput("field_output");
  obj->CreateTextOutput("field_output");
  
  return 0;
}
